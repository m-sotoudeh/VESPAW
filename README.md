<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="doc/img/vespaw-logo-dark.png">
    <img src="doc/img/vespaw-logo.png" alt="VESPAW" width="480">
  </picture>
</p>

<p align="center">
  <b>Variational electronic structure with projector augmented waves</b><br>
  plane waves &middot; MPI and OpenMP &middot; GPL-3.0
</p>

<p align="center">
  <a href="https://github.com/m-sotoudeh/VESPAW/actions/workflows/ci.yml"><img src="https://github.com/m-sotoudeh/VESPAW/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://github.com/m-sotoudeh/VESPAW/releases/latest"><img src="https://img.shields.io/github/v/release/m-sotoudeh/VESPAW" alt="Release"></a>
  <a href="https://doi.org/10.5281/zenodo.23295633"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23295633.svg" alt="DOI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue" alt="License: GPL-3.0"></a>
  <img src="https://img.shields.io/badge/Fortran-OpenMP%20%7C%20MPI-734f96" alt="Fortran, OpenMP, MPI">
  <img src="https://img.shields.io/badge/platforms-Linux%20%7C%20macOS-lightgrey" alt="Linux, macOS">
</p>

---

**VESPAW** is a plane-wave electronic-structure code based on the
projector augmented-wave (PAW) method. It is a development of
[CP-PAW](https://github.com/cp-paw/cp-paw), the original PAW code of
Peter E. Blöchl, and keeps its input format, its PAW engine and its
Car–Parrinello molecular dynamics.

Its main addition is **variational self-consistency**: the Kohn–Sham
Hamiltonian, not the density, is the variable of the self-consistency
loop. The Mermin free energy is evaluated for the actual output state of
every trial Hamiltonian, and a step is accepted only if this energy does
not increase. Every energy reported during the iteration is therefore an
upper bound of the converged result, and a restart that lowers the
energy shows that a calculation had not converged.

## Highlights

- **Variational SCF**: Hamiltonian minimization (HM, the default,
  `MIXER='VARIATIONAL'`) with the Davidson, RMM-DIIS or CG eigensolver;
  direct minimization (DM, `OPTIMIZER='EDFT'`) over orbitals and
  occupations as an independent variational reference.
- **Functionals**, checked against libxc (to machine precision) and VASP:
  LDA (PZ, PW92), GGA (PBE, RPBE, PW91, revPBE, PBEsol, BLYP, BP86),
  meta-GGA (SCAN and r²SCAN with orbital τ, SCAN-L), the local hybrids
  PBE0r, PBE-SL0r, SCAN0r and r²SCAN0r, DFT+U. Selected by name (`TYPE='PBEsol'`) or by
  the VASP tags (`'PS'`, `'RP'`, …).
- **Physics**: collinear and non-collinear spin, spin–orbit coupling,
  finite electron temperature, dipole correction, forces and stress,
  Car–Parrinello molecular dynamics, space-group symmetry.
- **Geometry optimization**: damped dynamics or L-BFGS
  (`!RDYN OPTIMIZER='LBFGS'`) for the atoms, jointly with the unit cell
  under an external pressure or a general stress tensor.
- **Output**: a structured report `<name>.prot` (input, system, model,
  SCF history, energies, k-points and occupations, atomic charges and
  moments, forces, stress, timing) in eV, Å and GPa; results in
  `results.h5`; the full protocol in `<name>.log`.
- **Parallel**: MPI over k-points and plane waves, OpenMP threads;
  native builds for Apple Silicon and x86_64.
- **Tested**: a regression suite checks every optimizer against
  Car–Parrinello dynamics, serial and with MPI.

Feature status (stable / experimental): [doc/FEATURES.md](doc/FEATURES.md).
Known problems: [doc/KNOWN_ISSUES.md](doc/KNOWN_ISSUES.md).

## Quick start

```sh
git clone https://github.com/m-sotoudeh/VESPAW.git && cd VESPAW
./install.sh                      # finds or installs what is needed, builds, checks
export PATH="$HOME/.vespaw/bin:$PATH"
vespaw doctor                     # verify the installation
```

`install.sh` needs no administrator rights. On macOS with Homebrew it uses
the Homebrew compilers and libraries; otherwise it creates a private
conda-forge environment in `~/.vespaw`. On clusters, load the compiler,
MPI, OpenBLAS, FFTW and HDF5 modules and run `./install.sh --system`.
Details: [doc/INSTALL.md](doc/INSTALL.md).

A calculation is started with `vespaw name.cntl` (control file
`name.cntl`, structure `name.strc`) or from YAML input with
`vespaw run input.yaml`; see [doc/USAGE.md](doc/USAGE.md) and the
tutorial.

## Documentation

| | |
|---|---|
| [doc/theory.pdf](doc/theory.pdf) | *VESPAW: Theory and Methods* — methods, functionals, implementation |
| [doc/tutorial.pdf](doc/tutorial.pdf) | Tutorial with checked examples ([examples/tutorial](examples/tutorial)) |
| [doc/manual.pdf](doc/manual.pdf) | Input reference (all blocks and keywords) |
| [doc/INSTALL.md](doc/INSTALL.md) | Installation and build options |
| [doc/USAGE.md](doc/USAGE.md) | Running calculations, YAML input, Python API, test suite |
| [potentials/README.md](potentials/README.md) | The PAW setup library |

## Layout

```
arch/         compiler and library settings (parameter files for configure)
build/        makefile templates
src/          Fortran sources
tools/        post-processing programs and scripts
python/       Python package (API and YAML input)
potentials/   PAW setup library
data/         force-field and tessellation data
testsuite/    regression tests        (python3 testsuite/run.py)
examples/     input examples
doc/          documentation (PDFs, guides; LaTeX sources in doc/src)
```

## Built on CP-PAW

VESPAW grew out of [CP-PAW](https://github.com/cp-paw/cp-paw), developed
by Peter E. Blöchl and co-workers at IBM Zurich and Clausthal University
of Technology since the early 1990s. The PAW engine, the input format,
the setup construction and the Car–Parrinello dynamics are CP-PAW's. The
history of CP-PAW, its coding conventions and its original documentation
are maintained in the [CP-PAW repository](https://github.com/cp-paw/cp-paw).
Upstream corrections are ported to VESPAW (for example the 2024
correction of the PBE reduced gradient).

## Citing

If you use VESPAW, please cite

- the software: M. Sotoudeh, *VESPAW: variational electronic structure
  with projector augmented waves*, Zenodo,
  [doi:10.5281/zenodo.23295633](https://doi.org/10.5281/zenodo.23295633);
- the variational self-consistency method: M. Sotoudeh, *The Kohn–Sham
  Hamiltonian as the variational object* (manuscript in preparation, 2026);
- the PAW method: P. E. Blöchl, Phys. Rev. B **50**, 17953 (1994);
- for the local hybrids: M. Sotoudeh *et al.*, Phys. Rev. B **95**,
  235150 (2017); P. E. Blöchl, C. F. J. Walther and Th. Pruschke,
  Phys. Rev. B **84**, 205101 (2011).

See [CITATION.cff](CITATION.cff).

## Contributing

Bug reports, test cases and code are welcome; see
[CONTRIBUTING.md](CONTRIBUTING.md). Please follow the
[code of conduct](CODE_OF_CONDUCT.md).

## License

VESPAW is free software under the GNU General Public License, version 3
([LICENSE](LICENSE)). It is derived from CP-PAW, distributed under the
same license; see [NOTICE](NOTICE).

## Contact

Mohsen Sotoudeh, Helmholtz Institute Ulm (HIU) / Karlsruhe Institute of
Technology (KIT) — mohsen.sotoudeh@kit.edu. Questions and bug reports:
please open an issue.
