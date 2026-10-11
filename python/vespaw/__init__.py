"""Thin Python API for VESPAW.

This wraps the real `vespaw` CLI dispatcher (tools/Scripts/vespaw) and
reads back the HDF5 results container the Fortran binary writes
(src/vespaw_ioroutines.f90, HDF5IO$WRITE) -- it does not reimplement any
physics. This module itself does not generate .cntl/.strc input for you
(see tools/Scripts/vespaw_convert.py for CIF/POSCAR/... -> .strc
import); for a higher-level YAML input frontend that lowers to .cntl/.strc
and hands back a ready-to-run Calculation, see vespaw.yaml_input.prepare_run.

It is the minimum needed to make VESPAW usable from a notebook, a workflow engine, or a
dataset-generation script instead of only from the command line.

Requires the .cntl deck to already contain a

    !HDF5OUT READ=T FILE='results.h5' !END

block under !CONTROL (see src/vespaw_ioroutines.f90, READIN_HDF5OUT) --
this API deliberately does not rewrite input decks on your behalf.

Example
-------
>>> from vespaw import Calculation
>>> calc = Calculation("mgo.cntl", env={"OPENBLAS_NUM_THREADS": "1"})
>>> result = calc.run()
>>> result.energy          # eV
>>> result.forces          # (nat, 3) ndarray, eV/Angstrom
>>> result.provenance["git_revision"]
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

__all__ = ["Calculation", "Result", "VespawError"]


def _sha256_of(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


class VespawError(RuntimeError):
    """The vespaw binary failed, or its HDF5 output is missing/malformed."""


_HDF5OUT_BLOCK_RE = re.compile(r"!HDF5OUT\b", re.IGNORECASE)
_HDF5OUT_FILE_RE = re.compile(
    r"!HDF5OUT\b[^!]*?\bFILE\s*=\s*'([^']+)'", re.IGNORECASE | re.DOTALL
)


@dataclass
class Result:
    """A read-only view of one VESPAW HDF5 results file.

    Units are the user-facing units of VESPAW (since 2026-10-02): eV for
    energies, Angstrom for lengths, eV/Angstrom for forces -- exactly as
    written by HDF5IO$WRITE and exactly as documented on each dataset's
    own 'units' HDF5 attribute (see e.g. `h5dump results.h5`).

    Schema note: this reads the /results/* namespace (results.h5 schema
    revision 2026-09-04) -- a breaking rename of the original
    2026-09-03 MVP's flat /energy, /forces groups. See
    src/vespaw_ioroutines.f90, HDF5IO$WRITE, for the authoritative
    Fortran-side schema.

    eigenvalues/occupations/kpoints/kpoint_weights/homo/lumo/is_metal are
    None whenever the run used dynamical occupations (DYNOCC DYN=T) --
    they are only meaningful for fixed-occupation runs, matching
    WAVES$REPORTEIG's own scope limit on the Fortran side. Not yet
    available at all (deliberately, not silently): stress tensor,
    real-space density, per-atom magnetic moments, SCF iteration
    history -- see the Fortran routine's header comment for why each one
    isn't wired up yet.
    """

    path: Path
    energy: float
    forces: Any            # (nat, 3) ndarray, eV/Angstrom
    positions: Any         # (nat, 3) ndarray, Angstrom, Cartesian; final
                            # (relaxed) position if !RDYN ran, not the
                            # initial one
    lattice_vectors: Any   # (3, 3) ndarray, Angstrom; row i = lattice vector i
    species: List[str]     # ATOM NAME per atom, in atom order
    provenance: Dict[str, str]
    eigenvalues: Optional[Any] = None   # (nspin, nkpt, nb) ndarray, eV
    occupations: Optional[Any] = None   # (nspin, nkpt, nb) ndarray, electrons
    kpoints: Optional[Any] = None       # (nkpt, 3) ndarray, fractional
    kpoint_weights: Optional[Any] = None  # (nkpt,) ndarray
    homo: Optional[float] = None        # eV
    lumo: Optional[float] = None        # eV
    is_metal: Optional[bool] = None

    @property
    def gap(self) -> Optional[float]:
        """LUMO - HOMO in eV, or None if eigenvalue data isn't
        available (dynamical-occupation run) -- never a negative number
        pretending to be a gap: that case is is_metal=True instead."""
        if self.homo is None or self.lumo is None:
            return None
        return self.lumo - self.homo

    @classmethod
    def from_file(cls, path: Union[str, Path]) -> "Result":
        import h5py  # local import: only needed once a result is read back

        path = Path(path)
        if not path.exists():
            raise VespawError(
                f"no HDF5 results file at {path} -- did the .cntl deck "
                "have '!HDF5OUT READ=T' and did the run finish normally?")

        def _str(v):
            return v.decode() if isinstance(v, bytes) else v

        with h5py.File(path, "r") as f:
            r = f["results"]
            energy = float(r["energy"][0])
            forces = r["forces"][()]
            positions = f["structure/positions"][()]
            lattice_vectors = f["structure/lattice_vectors"][()]
            species = _str(f["structure/species"][()]).split(",")
            provenance = {k: _str(f[f"provenance/{k}"][()])
                          for k in f["provenance"]}

            has_eig = "eigenvalues" in r
            eigenvalues = r["eigenvalues"][()] if has_eig else None
            occupations = r["occupations"][()] if has_eig else None
            kpoints = r["kpoints"][()] if has_eig else None
            kpoint_weights = r["kpoint_weights"][()] if has_eig else None
            homo = float(r["homo"][0]) if has_eig else None
            lumo = float(r["lumo"][0]) if has_eig else None
            is_metal = bool(r["is_metal"][0]) if has_eig else None

        return cls(path=path, energy=energy, forces=forces,
                    positions=positions, lattice_vectors=lattice_vectors,
                    species=species, provenance=provenance,
                    eigenvalues=eigenvalues, occupations=occupations,
                    kpoints=kpoints, kpoint_weights=kpoint_weights,
                    homo=homo, lumo=lumo, is_metal=is_metal)

    def write_run_json(self, path: Optional[Union[str, Path]] = None,
                        cntl_path: Optional[Union[str, Path]] = None
                        ) -> Path:
        """Writes the reproducibility record run.json: a flat JSON object of exactly the
        provenance fields HDF5IO$WRITE already captured, plus sha256
        hashes of the .cntl and its companion .strc (same basename, the
        CP-PAW convention -- see src/vespaw_ioroutines.f90) when found.

        No JSON schema/library dependency beyond the stdlib -- this is a
        small, known, flat set of string keys.
        """
        path = Path(path) if path else self.path.with_name("run.json")
        record: Dict[str, Any] = dict(self.provenance)
        record["results_file"] = str(self.path)
        if cntl_path is not None:
            cntl_path = Path(cntl_path)
            record["cntl_file"] = str(cntl_path)
            record["cntl_sha256"] = _sha256_of(cntl_path)
            strc_path = cntl_path.with_suffix(".strc")
            if strc_path.exists():
                record["strc_file"] = str(strc_path)
                record["strc_sha256"] = _sha256_of(strc_path)
        record["written"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        return path

@dataclass
class Calculation:
    """Runs one VESPAW .cntl deck via the `vespaw` CLI and reads back its
    HDF5 results.

    Parameters
    ----------
    cntl:
        Path to the .cntl deck. Must already contain a
        ``!HDF5OUT READ=T FILE='...'`` block -- see the module docstring.
    hdf5:
        Path to the HDF5 file the run will produce. If omitted, parsed
        out of the deck's own ``FILE=`` keyword, defaulting to
        ``RESULTS.H5`` next to the deck (matching HDF5IO_MODULE's own
        Fortran-side default).
    nprocs:
        If > 1, runs under ``mpirun -np NPROCS vespaw ...`` instead of
        plain serial dispatch.
    vespaw_bin:
        Path to the ``vespaw`` dispatcher script. Defaults to looking it
        up on ``PATH``.
    env:
        Extra environment variables merged into the subprocess call
        (e.g. ``{"OPENBLAS_NUM_THREADS": "1"}``, which avoids oversubscription under MPI).
    """

    cntl: Union[str, Path]
    hdf5: Optional[Union[str, Path]] = None
    nprocs: int = 1
    vespaw_bin: Optional[str] = None
    env: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.cntl = Path(self.cntl)
        if not self.cntl.exists():
            raise VespawError(f"no such .cntl file: {self.cntl}")

        text = self.cntl.read_text()
        if not _HDF5OUT_BLOCK_RE.search(text):
            raise VespawError(
                f"{self.cntl} has no !HDF5OUT block -- add "
                "\"!HDF5OUT READ=T FILE='results.h5' !END\" under "
                "!CONTROL (see src/vespaw_ioroutines.f90, "
                "READIN_HDF5OUT) before running it through this API.")

        if self.hdf5 is None:
            m = _HDF5OUT_FILE_RE.search(text)
            name = m.group(1) if m else "RESULTS.H5"
            self.hdf5 = self.cntl.parent / name
        else:
            self.hdf5 = Path(self.hdf5)

    def _resolve_vespaw_bin(self) -> str:
        if self.vespaw_bin:
            return self.vespaw_bin
        found = shutil.which("vespaw")
        if not found:
            raise VespawError(
                "no 'vespaw' dispatcher found on PATH -- pass "
                "vespaw_bin=... explicitly, or add bin/<arch>/ to PATH")
        return found

    def run(self, timeout: Optional[float] = None,
            write_run_json: bool = True, stream: bool = False) -> Result:
        """Runs the deck to completion and returns its Result.

        Raises VespawError (with the tail of stdout/stderr attached) if
        the binary exits non-zero, and again if the promised HDF5 file
        doesn't show up afterwards.

        write_run_json: also writes run.json (provenance + input/structure
        hashes) next to the HDF5 file -- see Result.write_run_json. Set
        False to skip (e.g. many quick calls in a loop where only
        result.energy etc. matter). The human-readable report is the
        <root>.prot file written by VESPAW itself.

        stream: let the binary's stdout/stderr go straight to this
        process's own (e.g. so a terminal watching a `vespaw run` CLI
        invocation sees the live-SCF table as it happens) instead of
        capturing it -- VespawError's failure message then has no
        stdout/stderr tail to attach, since none was captured.
        """
        vespaw_bin = self._resolve_vespaw_bin()
        cmd = ([vespaw_bin] if self.nprocs <= 1
               else ["mpirun", "-np", str(self.nprocs), vespaw_bin])
        cmd.append(self.cntl.name)

        run_env = dict(os.environ)
        run_env.update(self.env)

        proc = subprocess.run(
            cmd, cwd=self.cntl.parent, env=run_env, timeout=timeout,
            **({} if stream else dict(capture_output=True, text=True)),
        )
        if proc.returncode != 0:
            tail = ("" if stream else
                    f"\n--- stdout (tail) ---\n{proc.stdout[-4000:]}\n"
                    f"--- stderr (tail) ---\n{proc.stderr[-4000:]}")
            raise VespawError(
                f"'{' '.join(cmd)}' exited with code {proc.returncode}{tail}")

        result = Result.from_file(self.hdf5)
        if write_run_json:
            result.write_run_json(cntl_path=self.cntl)
        return result
