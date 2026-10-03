# Using VESPAW

This covers running a calculation from the command line and from Python.
For the full `.cntl`/`.strc` input-file reference — every block, every
keyword, the underlying physics — see `doc/manual.pdf` (the upstream
CP-PAW manual: VESPAW inherits nearly all of CP-PAW's input format
unchanged). This page only covers what's specific to VESPAW.

## 1. The `vespaw` dispatcher

Once built (see [INSTALL.md](INSTALL.md)), `bin/<arch>/vespaw` is a
single command that dispatches by file extension, and auto-selects the
parallel build when launched under `mpirun`/`mpiexec`/`srun`:

```
usage: vespaw [--np] FILE.<cntl|dcntl|dpcntl|bcntl|pcntl|gcntl> [args...]
       vespaw run FILE.yaml [args...]   (or just: vespaw FILE.yaml)
```

| Extension | Runs |
|---|---|
| `.cntl` | the SCF/MD driver (the main calculation) |
| `.dcntl` | density-of-states post-processing |
| `.dpcntl` | DOS plotting |
| `.bcntl` | band-structure calculation |
| `.pcntl` | pre-optimizer |
| `.gcntl` | "grab" post-processor |
| `.yaml`/`.yml`, or `vespaw run ...` | the YAML input frontend (§3 below) |
| `vespaw convert ...` | structure-file importer (§2 below) |

```sh
vespaw foo.cntl                          # serial
OPENBLAS_NUM_THREADS=1 mpirun -np 4 vespaw foo.cntl   # parallel — see INSTALL.md
```

Other CP-PAW tools (`strc`, `tostrc`, `toxyz`, `tra`, `wave`,
`murnaghan`, ...) aren't yet wired into this dispatcher — call their
`vespaw_<name>.x` binary directly.

## 2. Building a structure file

`.strc` files can be hand-written (see `doc/manual.pdf` and
`installation/examples/template_full.cntl` for the block syntax) or
generated from any ASE-readable structure file (CIF, POSCAR, ...) plus
one or more species definitions from the [`potentials`](../potentials/README.md)
library:

```sh
vespaw convert POSCAR mgo.strc \
    --species-from potentials/mg_ndlss_v0.strc-species potentials/o_ndlss_v0.strc-species \
    --kdiv 4 4 4 --empty 4
```

## 3. A minimal worked example

The fastest way to see a real, working `.cntl` deck is to copy one of
the checked-in regression cases (`testsuite/pbe_atom/atom.cntl` is a
small single-atom PBE run) rather than starting from a blank page.
`installation/examples/template_full.cntl` is the other reference worth
knowing about — an extensively annotated **template**, not a runnable
deck (no atoms/structure block), meant to be copied from when you need
to see every VESPAW-specific `.cntl` extension in one place.

A hand-written minimal deck, PBE on a pre-built `.strc`:

```
!CONTROL
  !GENERIC NSTEP=60 DT=5.0 NWRITE=10 START=T !END
  !FOURIER EPWPSI=25. CDUAL=2.0 !END
  !DFT TYPE='PBE' !END
  !PSIDYN OPTIMIZER='DAV' SAFEORTHO=F STOP=T MIXER='PULAY' !END
!END
!EOB
```

```sh
vespaw foo.cntl
```

Output lands in `foo.prot` (the human-readable run log — the live
per-iteration energy/residual table plus a final `RUN-TIME REPORT`) and
`foo.rstrt` (the restart/wavefunction file for continuing the run). Add
`!HDF5OUT READ=T FILE='results.h5' !END` under `!CONTROL` to also get a
structured HDF5 results file readable from Python (§5).

## 4. The YAML input frontend

For the common case (a plane-wave DFT calculation: functional, cutoff,
k-mesh, SCF settings), a much shorter YAML input lowers to the
`.cntl`/`.strc` pair above automatically:

```yaml
vespaw:
  schema: 1

system:
  structure: POSCAR

electronic:
  xc:
    functional: PBE          # any name of the TYPE list, or PBE0R / PBE-SL0R / SCAN0R -- see FEATURES.md
  basis:
    wavefunction_cutoff: 25 Ry
  kpoints:
    mesh: [1, 1, 1]

scf:
  max_iterations: 10
  diagonalizer:
    method: cg                # cg | davidson | rmm-diis

ions:
  mode: static                 # or relax
```

```sh
vespaw run foo.yaml
```

This covers the common structured fields only (functional, basis,
k-points, SCF/mixing, ionic mode). Anything else — any `.cntl` keyword
without a dedicated YAML field — is reachable through a `raw:` section,
a literal mirror of the `.cntl` block tree:

```yaml
raw:
  FOURIER:
    CDUAL: 2.5
  MERMIN:
    READ: false
```

See `python/vespaw/yaml_input.py`'s module docstring for the exact
schema and what's deliberately not supported yet (it fails with a clear
message naming the offending key rather than silently ignoring
anything), and `testsuite/yaml_frontend/`, `testsuite/yaml_raw/` for
two more complete, runnable examples.

## 5. From Python

```sh
pip install -e python/
```

```python
from vespaw import Calculation

calc = Calculation("foo.cntl", env={"OPENBLAS_NUM_THREADS": "1"})
result = calc.run()

result.energy          # Hartree
result.forces          # (nat, 3) ndarray, Hartree/bohr
result.gap              # LUMO - HOMO, Hartree (None for dynamical-occupation runs)
result.provenance["git_revision"]
```

Requires the deck to already contain a `!HDF5OUT READ=T
FILE='results.h5' !END` block (this API reads results back, it does not
rewrite your input for you). `Calculation.run()` also writes a
`run.json` provenance record next to the results file; the human-readable
report is the `<root>.prot` file that VESPAW writes for every run. For YAML input, use
`vespaw.yaml_input.prepare_run(yaml_path)` instead, which validates and
lowers the YAML for you and hands back a ready-to-run `Calculation`.

See `python/vespaw/__init__.py`'s module docstring for the full API and
units (everything is VESPAW's own internal atomic units: bohr, Hartree,
Hartree/bohr).

## 6. Regression suite

```sh
python3 testsuite/run.py            # every checked-in case
python3 testsuite/run.py pbe_atom   # just one
```

A small, fast, checked-in set of known-good energies (distinct from the
gitignored exploratory `tests/` tree) — see `testsuite/run.py`'s own
module docstring.
