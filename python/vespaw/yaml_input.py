"""YAML input frontend for VESPAW -- the "simple case" subset of Phase 6's
target input-language design (see the dev notes, "Input language"),
plus a ``raw:`` passthrough that reaches every ``.cntl`` keyword, not just
the ones this module has bespoke schema fields for.

Scope of the structured ("simple case") fields, exactly as confirmed:
``vespaw.schema``, ``system.structure``, ``electronic.xc``/``basis``/
``kpoints``, ``scf``, and ``ions.mode``. Everything the design's full
design also specifies -- ``execution``/``restart``/``output`` policy
sections, the ``vespaw validate``/``inspect-input``/``migrate`` CLI
subcommands, and caret-precise line/column error messages -- is
deliberately NOT implemented here; an input that uses those *structured*
sections fails with a clear "not supported in this pass" message rather
than being silently ignored.

``raw:``, added alongside the structured fields, is the answer to "how do I
reach a `.cntl` keyword this module doesn't have a schema field for" --
rather than hand-mapping CP-PAW's entire (large, long-tailed) `.cntl`
grammar one keyword at a time, ``raw:`` is a direct, literal mirror of the
`.cntl` block tree: each key is a `!BLOCK` name, nested dicts are
sub-blocks, and any non-dict value is a `KEY=VALUE` pair inside that block
(a list value becomes a space-joined multi-value `KEY=V1 V2 V3`, matching
keywords like `DIV=3 3 3`). ``raw:`` values are merged on top of whatever
the structured fields already produced -- raw wins on a conflicting
key/block, and anything raw introduces that the structured lowering never
touched is simply added. Example::

    raw:
      GENERIC:
        RSTRTTYPE: DYNAMIC
      PSIDYN:
        FRIC: 0.01
        AUTO:
          READ: false
      MERMIN:
        READ: false

This module lowers a validated input.yaml straight to a CP-PAW ``.cntl``
deck (string) plus a ``.strc`` file (written by shelling out to the existing
``vespaw_convert.py`` structure importer, reusing its element/species-block
matching against the ``potentials/`` library), then hands back a
:class:`vespaw.Calculation` ready to ``.run()``.

Only PyYAML + hand-rolled validation is used -- no jsonschema/pydantic --
per the confirmed library choice (lightest dependency, full control over
error wording; caret-precise pointing is accepted as a stretch goal not
polished this pass).
"""
from __future__ import annotations

import re
import subprocess
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

from . import Calculation

_REPO_ROOT = Path(__file__).resolve().parents[2]
_CONVERT_SCRIPT = _REPO_ROOT / "tools" / "Scripts" / "vespaw_convert.py"
_VESPAW_POT_DIR = _REPO_ROOT / "potentials"

# Named functional aliases actually wired into `.cntl` (see
# the dev notes, "Named functional aliases", READIN_DFTTYPEFROMNAME in
# src/vespaw_ioroutines.f90). Deliberately only the well-established names,
# not all ~25 numeric DFT_INITIALIZE codes.
_PLAIN_FUNCTIONALS = {
    "NONE", "PZ", "PW92", "PBE", "RPBE", "PW91", "REVPBE", "PBESOL", "SCANL", "SCAN",
    # VASP GGA tags, accepted as synonyms by READIN_DFTTYPEFROMNAME
    "PE", "RP", "91", "RE", "PS",
}
# Onsite hybrids are NOT aliases in the Fortran layer; this frontend expands
# them into (global !DFT TYPE, extra !DFT keys, extra !NTBO keys), with
# !NTBO MODUS='HYBRID' HFWEIGHT=<weight>. Names as of 2026-10-03:
#   PBE0R    : PBE everywhere + onsite exact exchange
#   PBE-SL0R : PBE everywhere; onsite exchange SCAN-L mixed with exact
#              exchange (DCTYPE='SCANL'; this was called SCAN0R before)
#   SCAN0R   : orbital-tau SCAN everywhere (TSCANDELOC=T) + onsite exact
#              exchange, double counting with the same SCAN exchange
_HYBRID_FUNCTIONALS = {
    "PBE0R": ("PBE", {}, {}),
    "PBE-SL0R": ("PBE", {}, {"DCTYPE": "SCANL"}),
    "SCAN0R": ("SCAN", {"TSCANDELOC": True}, {}),
}
# plain orbital-tau SCAN needs the delocalized (grid + one-center) evaluation
_PLAIN_DFT_EXTRA = {"SCAN": {"TSCANDELOC": True}}

_DIAGONALIZER_MAP = {
    "cg": "CG",
    "davidson": "DAV",
    "rmm-diis": "RMM",
    "rmm_diis": "RMM",
    "rmm": "RMM",
}
_MIXER_MAP = {
    "anderson": "ANDERSON",
    "pulay": "PULAY",
    "variational": "VARIATIONAL",
}

_HARTREE_IN_EV = 27.211386245988
_ENERGY_TO_HARTREE = {
    "ha": 1.0,
    "hartree": 1.0,
    "ry": 0.5,
    "rydberg": 0.5,
    "ev": 1.0 / _HARTREE_IN_EV,
}

_QUANTITY_RE = re.compile(
    r"^\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*([A-Za-z]+)?\s*$")


class Block:
    """One `!NAME ... !END` node in the `.cntl` tree -- the internal
    representation both the structured lowering and the `raw:` overlay
    populate before final text rendering. `keys` holds `KEY=VALUE` pairs
    (value: str/int/float/bool, or a list of those for a multi-value
    keyword like `DIV=3 3 3`); `children` holds nested `!SUBBLOCK`s, in
    insertion order (matters only for readability of the generated
    `.cntl`, not for correctness -- LINKEDLIST's tokenizer doesn't care
    about sibling order)."""

    def __init__(self, name: str):
        self.name = name
        self.keys: "OrderedDict[str, Any]" = OrderedDict()
        self.children: "OrderedDict[str, Block]" = OrderedDict()

    def child(self, name: str) -> "Block":
        if name not in self.children:
            self.children[name] = Block(name)
        return self.children[name]


def _merge_raw_into(blocks: "OrderedDict[str, Block]", raw: Dict[str, Any],
                     path: str) -> None:
    """Top-level entry point: `raw` is `cfg['raw']`, `blocks` is the
    structured lowering's top-level block dict (GENERIC/FOURIER/DFT/...)."""
    if not isinstance(raw, dict):
        _fail(path, f"expected a mapping of !BLOCK names, got {raw!r}")
    for name, val in raw.items():
        if not isinstance(val, dict):
            _fail(f"{path}.{name}",
                  "a top-level raw entry must itself be a mapping (it's a "
                  "!BLOCK, so its contents are KEY=VALUE pairs and/or "
                  "nested sub-blocks)")
        if name not in blocks:
            blocks[name] = Block(name)
        _merge_block(blocks[name], val, f"{path}.{name}")


def _merge_block(block: "Block", val: Dict[str, Any], path: str) -> None:
    for k, v in val.items():
        if isinstance(v, dict):
            _merge_block(block.child(k), v, f"{path}.{k}")
        elif isinstance(v, list):
            block.keys[k] = [_check_scalar(x, f"{path}.{k}") for x in v]
        else:
            block.keys[k] = _check_scalar(v, f"{path}.{k}")


def _check_scalar(v: Any, path: str) -> Any:
    if not isinstance(v, (bool, int, float, str)):
        _fail(path, f"expected a plain value (string/number/bool), got {v!r}")
    return v


def _format_scalar(v: Any, path: str) -> str:
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        s = f"{v:.10g}"
        if "." not in s and "e" not in s and "E" not in s:
            s += "."
        return s
    _fail(path, f"unsupported value type: {v!r}")


def _render_kv(key: str, value: Any, path: str) -> str:
    if isinstance(value, list):
        return f"{key}={' '.join(_format_scalar(x, path) for x in value)}"
    return f"{key}={_format_scalar(value, path)}"


def _render_block(block: "Block", indent: str, path: str) -> List[str]:
    kw = " ".join(_render_kv(k, v, f"{path}.{k}") for k, v in block.keys.items())
    header = f"{indent}!{block.name}" + (f" {kw}" if kw else "")
    lines = [header]
    for child in block.children.values():
        lines.extend(_render_block(child, indent + "  ", f"{path}.{child.name}"))
    lines.append(f"{indent}!END")
    return lines


class VespawInputError(Exception):
    """Raised for anything wrong with an input.yaml -- a missing/misspelled
    key, a value of the wrong type or out of range, or a schema section this
    pass doesn't support yet. The message always names the offending
    dotted-path key, e.g. "electronic.xc.functional", so the user knows
    exactly where to look (this substitutes, for now, for the caret-precise
    line/column pointer the design design calls for)."""


def _fail(path: str, message: str) -> None:
    raise VespawInputError(f"{path}: {message}")


def _get(cfg: Dict[str, Any], path: str, required: bool = True,
          default: Any = None) -> Any:
    node: Any = cfg
    parts = path.split(".")
    for i, key in enumerate(parts):
        if not isinstance(node, dict) or key not in node:
            if required:
                _fail(path, "missing required key")
            return default
        node = node[key]
    return node


def _reject_unsupported(cfg: Dict[str, Any]) -> None:
    """Fail clearly on any *structured* section the design design
    specifies but this pass does not implement, instead of silently
    ignoring it. Everything rejected here is still reachable through the
    `raw:` passthrough (see the module docstring) if the underlying
    `.cntl` keyword exists -- these are specifically the not-yet-designed
    Python-side/structured-schema concepts (hardware selection, restart
    policy, output-property selection, magnetism, ...), not gaps in what
    `.cntl` itself can express."""
    for key in ("execution", "restart", "output", "paw", "properties"):
        if key in cfg:
            _fail(key, "not supported as a structured field in this pass "
                        "(see the dev notes, Phase 6 'Input language') "
                        "-- use raw: to reach the underlying .cntl keyword "
                        "directly if there is one")
    system = cfg.get("system", {})
    if isinstance(system, dict):
        for key in ("spin", "charge"):
            if key in system and key == "spin":
                _fail(f"system.{key}", "not supported in this pass")
    electronic = cfg.get("electronic", {})
    if isinstance(electronic, dict) and "occupations" in electronic:
        _fail("electronic.occupations", "not supported in this pass")
    scf = cfg.get("scf", {})
    if isinstance(scf, dict):
        mixing = scf.get("mixing", {})
        if isinstance(mixing, dict):
            for key in ("beta", "history"):
                if key in mixing:
                    _fail(f"scf.mixing.{key}",
                          "not supported as a structured field -- no "
                          "matching .cntl keyword is known to exist "
                          "(only method and kerker.q0 are wired up); use "
                          "raw: PSIDYN if there turns out to be one")
            kerker = mixing.get("kerker", {})
            if isinstance(kerker, dict) and kerker.get("enabled") is False:
                _fail("scf.mixing.kerker.enabled",
                      "not supported as a structured field -- Kerker "
                      "preconditioning cannot be disabled through .cntl")
    ions = cfg.get("ions", {})
    if isinstance(ions, dict):
        extra = set(ions) - {"mode"}
        if extra:
            _fail(f"ions.{sorted(extra)[0]}",
                  "not supported as a structured field in this pass -- "
                  "only ions.mode is implemented (static or relax); use "
                  "raw: RDYN for other ionic-dynamics keywords")


def _parse_quantity(raw: Any, path: str) -> Tuple[float, Optional[str]]:
    """Splits "40 Ry" / "1.0e-8Ha" / 500 into (value, unit-or-None)."""
    if isinstance(raw, (int, float)):
        return float(raw), None
    if not isinstance(raw, str):
        _fail(path, f"expected a number or 'number unit' string, got {raw!r}")
    m = _QUANTITY_RE.match(raw)
    if not m:
        _fail(path, f"could not parse quantity {raw!r} "
                     "(expected e.g. '40 Ry', '500 eV', or a plain number)")
    value = float(m.group(1))
    unit = m.group(2).lower() if m.group(2) else None
    return value, unit


def _energy_to_hartree(raw: Any, path: str, default_unit: str = "ha") -> float:
    value, unit = _parse_quantity(raw, path)
    unit = (unit or default_unit).lower()
    if unit not in _ENERGY_TO_HARTREE:
        _fail(path, f"unknown energy unit '{unit}' (expected Ha, Ry, or eV)")
    return value * _ENERGY_TO_HARTREE[unit]


def _positive_int(raw: Any, path: str) -> int:
    if not isinstance(raw, int) or isinstance(raw, bool) or raw < 1:
        _fail(path, f"expected a positive integer, got {raw!r}")
    return raw


def _int_triple(raw: Any, path: str) -> Tuple[int, int, int]:
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        _fail(path, f"expected a 3-element list, got {raw!r}")
    out = []
    for v in raw:
        if not isinstance(v, int) or isinstance(v, bool):
            _fail(path, f"expected 3 integers, got {raw!r}")
        out.append(v)
    return tuple(out)  # type: ignore[return-value]


@dataclass
class LoweredInput:
    cntl_text: str
    structure_path: Path
    kdiv: Tuple[int, int, int]
    kshift: Tuple[int, int, int]


def load(yaml_path: Path) -> Dict[str, Any]:
    with open(yaml_path) as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise VespawInputError(
            f"{yaml_path}: top-level YAML document must be a mapping "
            f"(got {type(cfg).__name__})")
    return cfg


def validate_and_lower(cfg: Dict[str, Any], yaml_dir: Path) -> LoweredInput:
    """Validates `cfg` against the "simple case" schema subset and lowers it
    to a `.cntl` deck (as text) plus the structure-conversion parameters.
    Raises VespawInputError on the first problem found."""
    schema = _get(cfg, "vespaw.schema")
    if schema != 1:
        _fail("vespaw.schema", f"only schema 1 is supported, got {schema!r}")

    _reject_unsupported(cfg)

    structure = _get(cfg, "system.structure")
    if not isinstance(structure, str):
        _fail("system.structure", "expected a path to a structure file")
    structure_path = (yaml_dir / structure).resolve()
    if not structure_path.exists():
        _fail("system.structure", f"file not found: {structure_path}")

    # -- electronic.xc --------------------------------------------------
    functional_raw = _get(cfg, "electronic.xc.functional")
    if not isinstance(functional_raw, str):
        _fail("electronic.xc.functional", "expected a string")
    functional = functional_raw.upper()
    hybrid_weight = _get(cfg, "electronic.xc.hybrid.weight", required=False)
    if functional in _HYBRID_FUNCTIONALS:
        if hybrid_weight is None:
            _fail("electronic.xc.hybrid.weight",
                  f"required for hybrid functional '{functional_raw}'")
        if not isinstance(hybrid_weight, (int, float)) or not (0.0 <= hybrid_weight <= 1.0):
            _fail("electronic.xc.hybrid.weight",
                  f"expected a number in [0, 1], got {hybrid_weight!r}")
        dft_type, dft_extra, ntbo_extra = _HYBRID_FUNCTIONALS[functional]
    elif functional in _PLAIN_FUNCTIONALS:
        if hybrid_weight is not None:
            _fail("electronic.xc.hybrid.weight",
                  f"given, but '{functional_raw}' is not a hybrid functional "
                  "(use PBE0r, PBE-SL0r or SCAN0r for onsite-hybrid mixing)")
        dft_type = functional
        dft_extra = _PLAIN_DFT_EXTRA.get(functional, {})
        ntbo_extra = None
    else:
        known = ", ".join(sorted(_PLAIN_FUNCTIONALS | set(_HYBRID_FUNCTIONALS)))
        _fail("electronic.xc.functional",
              f"unknown functional '{functional_raw}' (known: {known})")

    cutoff_raw = _get(cfg, "electronic.basis.wavefunction_cutoff")
    cutoff_value, cutoff_unit = _parse_quantity(
        cutoff_raw, "electronic.basis.wavefunction_cutoff")
    cutoff_unit = (cutoff_unit or "ry").lower()
    if cutoff_unit not in _ENERGY_TO_HARTREE:
        _fail("electronic.basis.wavefunction_cutoff",
              f"unknown energy unit '{cutoff_unit}' (expected Ha, Ry, or eV)")
    epwpsi_ry = cutoff_value * _ENERGY_TO_HARTREE[cutoff_unit] / 0.5

    mesh = _int_triple(_get(cfg, "electronic.kpoints.mesh"),
                        "electronic.kpoints.mesh")
    if any(m < 1 for m in mesh):
        _fail("electronic.kpoints.mesh", f"all entries must be >= 1, got {mesh}")
    shift = _int_triple(_get(cfg, "electronic.kpoints.shift", required=False,
                              default=[0, 0, 0]),
                         "electronic.kpoints.shift")

    # -- scf --------------------------------------------------------------
    max_iter = _positive_int(_get(cfg, "scf.max_iterations"),
                              "scf.max_iterations")

    ediff = _get(cfg, "scf.convergence.energy", required=False)
    ediff_ha = (_energy_to_hartree(ediff, "scf.convergence.energy")
                if ediff is not None else None)

    mixtol = _get(cfg, "scf.convergence.density", required=False)
    if mixtol is not None and not isinstance(mixtol, (int, float)):
        _fail("scf.convergence.density", f"expected a number, got {mixtol!r}")

    diag_method_raw = _get(cfg, "scf.diagonalizer.method", required=False,
                            default="rmm-diis")
    diag_method_key = str(diag_method_raw).lower()
    if diag_method_key not in _DIAGONALIZER_MAP:
        _fail("scf.diagonalizer.method",
              f"unknown method '{diag_method_raw}' "
              f"(known: {', '.join(sorted(set(_DIAGONALIZER_MAP)))})")
    optimizer = _DIAGONALIZER_MAP[diag_method_key]

    nrmmx = _get(cfg, "scf.diagonalizer.max_iterations", required=False)
    if nrmmx is not None:
        nrmmx = _positive_int(nrmmx, "scf.diagonalizer.max_iterations")

    rmmdelay = _get(cfg, "scf.diagonalizer.warm_start.iterations",
                     required=False)
    if rmmdelay is not None:
        rmmdelay = _positive_int(rmmdelay,
                                  "scf.diagonalizer.warm_start.iterations")

    mixer_method_raw = _get(cfg, "scf.mixing.method", required=False,
                             default="variational")
    mixer_method_key = str(mixer_method_raw).lower()
    if mixer_method_key not in _MIXER_MAP:
        _fail("scf.mixing.method",
              f"unknown method '{mixer_method_raw}' "
              f"(known: {', '.join(sorted(_MIXER_MAP))})")
    mixer = _MIXER_MAP[mixer_method_key]

    mixg0sq = _get(cfg, "scf.mixing.kerker.q0", required=False)
    if mixg0sq is not None and not isinstance(mixg0sq, (int, float)):
        _fail("scf.mixing.kerker.q0", f"expected a number, got {mixg0sq!r}")

    # -- ions ---------------------------------------------------------
    ions_mode = _get(cfg, "ions.mode", required=False, default="static")
    if ions_mode not in ("static", "relax"):
        _fail("ions.mode", f"expected 'static' or 'relax', got {ions_mode!r}")

    cntl_text = _build_cntl(
        nstep=max_iter, epwpsi_ry=epwpsi_ry, dft_type=dft_type,
        dft_extra=dft_extra, ntbo_extra=ntbo_extra, hybrid_weight=hybrid_weight,
        optimizer=optimizer, mixer=mixer, mixg0sq=mixg0sq,
        mixtol=mixtol, ediff_ha=ediff_ha, nrmmx=nrmmx, rmmdelay=rmmdelay,
        ions_mode=ions_mode, raw=cfg.get("raw"),
    )
    return LoweredInput(cntl_text=cntl_text, structure_path=structure_path,
                         kdiv=mesh, kshift=shift)


def _build_cntl(*, nstep: int, epwpsi_ry: float, dft_type: str,
                 dft_extra: Dict[str, Any],
                 ntbo_extra: Optional[Dict[str, str]],
                 hybrid_weight: Optional[float], optimizer: str, mixer: str,
                 mixg0sq: Optional[float], mixtol: Optional[float],
                 ediff_ha: Optional[float], nrmmx: Optional[int],
                 rmmdelay: Optional[int], ions_mode: str,
                 raw: Optional[Dict[str, Any]]) -> str:
    nwrite = max(1, min(50, nstep))

    control: "OrderedDict[str, Block]" = OrderedDict()

    generic = control["GENERIC"] = Block("GENERIC")
    generic.keys.update(NSTEP=nstep, DT=5.0, NWRITE=nwrite, START=True)

    fourier = control["FOURIER"] = Block("FOURIER")
    fourier.keys.update({"EPWPSI": round(epwpsi_ry, 6), "CDUAL": 2.0})

    dft = control["DFT"] = Block("DFT")
    dft.keys["TYPE"] = dft_type
    dft.keys.update(dft_extra)
    if ntbo_extra is not None:
        ntbo = dft.child("NTBO")
        ntbo.keys["MODUS"] = "HYBRID"
        ntbo.keys["HFWEIGHT"] = round(float(hybrid_weight), 6)
        ntbo.keys.update(ntbo_extra)

    psidyn = control["PSIDYN"] = Block("PSIDYN")
    psidyn.keys.update(OPTIMIZER=optimizer, SAFEORTHO=True, STOP=True,
                        MIXER=mixer)
    if mixg0sq is not None:
        psidyn.keys["MIXG0SQ"] = float(mixg0sq)
    if mixtol is not None:
        psidyn.keys["MIXTOL"] = float(mixtol)
    if ediff_ha is not None:
        psidyn.keys["EDIFF"] = float(ediff_ha)
    if nrmmx is not None:
        psidyn.keys["NRMMX"] = nrmmx
    if rmmdelay is not None:
        psidyn.keys["RMMDELAY"] = rmmdelay

    if ions_mode == "relax":
        rdyn = control["RDYN"] = Block("RDYN")
        rdyn.keys.update(READ=True, STOP=True, FRIC=0.3)

    hdf5out = control["HDF5OUT"] = Block("HDF5OUT")
    hdf5out.keys.update(READ=True, FILE="results.h5")

    if raw is not None:
        _merge_raw_into(control, raw, "raw")

    lines = ["!CONTROL"]
    for block in control.values():
        lines.extend(_render_block(block, "  ", block.name))
    lines.append("!END")
    lines.append("!EOB")
    return "\n".join(lines) + "\n"


def _run_convert(structure_path: Path, strc_path: Path,
                  kdiv: Tuple[int, int, int],
                  kshift: Tuple[int, int, int]) -> None:
    species_files = sorted(_VESPAW_POT_DIR.glob("*.strc-species"))
    if not species_files:
        raise VespawInputError(
            f"no *.strc-species files found in {_VESPAW_POT_DIR} -- the "
            "potentials species library is required to build a .strc file")
    argv = [
        sys.executable, str(_CONVERT_SCRIPT),
        str(structure_path), str(strc_path),
        "--species-from", *[str(p) for p in species_files],
        "--kdiv", str(kdiv[0]), str(kdiv[1]), str(kdiv[2]),
        "--kshift", str(kshift[0]), str(kshift[1]), str(kshift[2]),
    ]
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode != 0:
        raise VespawInputError(
            f"vespaw convert failed while building {strc_path}:\n"
            f"{proc.stdout}\n{proc.stderr}")


def prepare_run(yaml_path: Path, workdir: Optional[Path] = None) -> Calculation:
    """Loads, validates, and lowers `yaml_path` into a `.cntl`+`.strc` pair
    under `workdir` (default: alongside the YAML file), then returns a
    ready-to-run vespaw.Calculation. Does not run anything itself -- call
    .run() on the result."""
    yaml_path = Path(yaml_path)
    cfg = load(yaml_path)
    workdir = Path(workdir) if workdir is not None else yaml_path.parent
    workdir.mkdir(parents=True, exist_ok=True)

    lowered = validate_and_lower(cfg, yaml_path.parent)

    stem = yaml_path.stem
    cntl_path = workdir / f"{stem}.cntl"
    strc_path = workdir / f"{stem}.strc"

    cntl_path.write_text(lowered.cntl_text)
    _run_convert(lowered.structure_path, strc_path, lowered.kdiv, lowered.kshift)

    return Calculation(cntl=cntl_path)
