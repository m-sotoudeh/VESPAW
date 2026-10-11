# Using VESPAW

This page covers running a calculation from the command line and from
Python. The input reference (every block and keyword of `.cntl` and
`.strc`) is `doc/manual.pdf`; the methods are described in
`doc/theory.pdf`, and `doc/tutorial.pdf` works through examples.

## 1. The `vespaw` command

After the installation (see [INSTALL.md](INSTALL.md)), `vespaw` runs a
calculation and dispatches by file extension. Under `mpirun`, `mpiexec`
or `srun` it selects the parallel build.

```
usage: vespaw [--np] FILE.<cntl|dcntl|dpcntl|bcntl|pcntl|gcntl> [args...]
       vespaw run FILE.yaml [args...]   (or just: vespaw FILE.yaml)
```

| Extension | Runs |
|---|---|
| `.cntl` | the SCF/MD driver (the main calculation) |
| `.dcntl` | density-of-states post-processing |
| `.dpcntl` | DOS plotting |
| `.bcntl` | band structure |
| `.pcntl` | pre-optimizer |
| `.gcntl` | "grab" post-processor |
| `.yaml`/`.yml`, or `vespaw run ...` | the YAML input frontend (§4) |
| `vespaw convert ...` | structure-file importer (§2) |
| `vespaw doctor` | check the installation |

```sh
vespaw foo.cntl                                                      # serial
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpirun -np 4 vespaw foo.cntl # MPI
```

The other CP-PAW tools (`strc`, `tostrc`, `toxyz`, `tra`, `wave`,
`murnaghan`, …) are called by their executable names
`vespaw_<name>.x` in `bin/<arch>/`.

## 2. Building a structure file

`.strc` files can be written by hand (see `doc/manual.pdf` and the
annotated template `examples/installation/template_full.cntl`) or
generated from any structure file that ASE can read (CIF, POSCAR, …)
plus species definitions from the [setup library](../potentials/README.md):

```sh
vespaw convert POSCAR mgo.strc \
    --species-from potentials/mg_ndlss_v0.strc-species potentials/o_ndlss_v0.strc-species \
    --kdiv 4 4 4 --empty 4
```

## 3. A minimal calculation

The decks in `examples/tutorial/` and `testsuite/` are complete, checked
inputs and the best starting point. A minimal control file for a PBE
calculation with the default variational SCF:

```
!CONTROL
  !GENERIC NSTEP=300 DT=5.0 START=T !END
  !FOURIER EPWPSI=25. CDUAL=2.0 !END
  !DFT TYPE='PBE' !END
  !PSIDYN OPTIMIZER='DAV' MIXER='VARIATIONAL' SAFEORTHO=T STOP=T !END
!END
!EOB
```

```sh
vespaw foo.cntl
```

Output files:

- `foo.prot`: the report (input echo, system, model, SCF history,
  energies, k-points and occupations, charges and moments, forces,
  stress, timing) in eV, Å and GPa;
- `foo.log`: the full protocol of the run;
- `foo.rstrt`: the restart file (`START=F` continues from it).

`!HDF5OUT READ=T FILE='results.h5' !END` in `!CONTROL` writes the results
also to an HDF5 file, which the Python API reads (§5).

## 4. YAML input

For a standard plane-wave calculation (functional, cutoff, k-points,
SCF settings), a short YAML file is converted to `.cntl` and `.strc`:

```yaml
vespaw:
  schema: 1

system:
  structure: POSCAR

electronic:
  xc:
    functional: PBE          # any TYPE name, or PBE0R / PBE-SL0R / SCAN0R
  basis:
    wavefunction_cutoff: 25 Ry
  kpoints:
    mesh: [4, 4, 4]

scf:
  max_iterations: 100
  diagonalizer:
    method: davidson          # davidson | rmm-diis | cg

ions:
  mode: static                # or relax
```

```sh
vespaw run foo.yaml
```

Any `.cntl` keyword without a YAML field is reached through a `raw:`
section, which mirrors the `.cntl` block tree:

```yaml
raw:
  FOURIER:
    CDUAL: 2.5
  MERMIN:
    READ: false
```

The schema is documented in the docstring of
`python/vespaw/yaml_input.py`; unknown keys are rejected with a message
naming them. Runnable examples: `testsuite/yaml_frontend/`,
`testsuite/yaml_raw/`.

## 5. Python

```sh
pip install -e python/
```

```python
from vespaw import Calculation

calc = Calculation("foo.cntl", env={"OPENBLAS_NUM_THREADS": "1"})
result = calc.run()

result.energy          # eV
result.forces          # (nat, 3) ndarray, eV/Å
result.gap             # LUMO - HOMO in eV (None without eigenvalue data; see result.is_metal)
result.provenance["git_revision"]
```

The deck must contain `!HDF5OUT READ=T FILE='results.h5' !END`; the API
reads the results back and does not modify the input. `Calculation.run()`
also writes a provenance record `run.json` next to the results file. For
YAML input, `vespaw.yaml_input.prepare_run(yaml_path)` validates and
converts the input and returns a ready-to-run `Calculation`. The full API
is described in the docstring of `python/vespaw/__init__.py`.

## 6. Regression suite

```sh
python3 testsuite/run.py            # all cases
python3 testsuite/run.py pbe_atom   # one case
python3 testsuite/run.py si_pbe:dav_var   # one run of a case
VESPAW_REGRESSION_NP=4 python3 testsuite/run.py   # with 4 MPI tasks
```

Each case solves one system with several optimizers and checks energies,
gaps, forces and stresses against stored references and against each
other; see the docstring of `testsuite/run.py`.
