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
  meta-GGA (SCAN and r²SCAN with orbital τ, SCAN-L), DFT+U. Selected by
  name (`TYPE='PBEsol'`) or by the VASP tags (`'PS'`, `'RP'`, …).
- **Hybrid functionals**: plane-wave PAW Fock exchange for PBE0, HSE06,
  SCAN0 and r²SCAN0 (`!DFT!FOCK AEXX=0.25 HFSCREEN=0.2 !END`), with
  adaptively compressed exchange, k-points, symmetry, spin–orbit
  spinors, forces, stress and MPI; band-gap shifts within 2 % of VASP
  for Si, C, MgO and MnO.
- **On-site hybrid functionals**: PBE0r, PBE-SL0r, SCAN0r and r²SCAN0r
  (`!CONTROL!DFT!NTBO MODUS='HYBRID' HFWEIGHT=0.25 !END`) add exact
  exchange only inside the atomic spheres, evaluated in a basis of
  natural tight-binding orbitals, at nearly semilocal cost. Suited to
  localized d electrons in transition-metal oxides and to cells of about
  a hundred atoms. The operator acts in all eigensolvers, with forces,
  stress and MPI.
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

## System requirements and demo

**Operating systems.** Linux and macOS, on x86_64 or Apple Silicon. No
special hardware is needed. Tested on:

| System | Compiler and libraries |
|---|---|
| macOS 26.6, Apple M3 laptop (8 cores) | gfortran 15.1, OpenBLAS 0.3.29, FFTW 3.3.10, HDF5 1.14.5, Open MPI 5.0.7 (Homebrew) |
| macOS 15, arm64 (CI) | gfortran, OpenBLAS, FFTW, HDF5, Open MPI (Homebrew) |
| Ubuntu 24.04, x86_64 (CI) | gfortran, OpenBLAS, FFTW3, HDF5, Open MPI (system packages) |
| Linux cluster (bwHPC JUSTUS 2) | environment modules |

**Dependencies.**
- Required: a Fortran compiler (gfortran 13 to 15 are tested), BLAS/LAPACK
  (OpenBLAS), FFTW3, and HDF5 with the Fortran interface.
- Optional: an MPI library for the parallel build; Python 3 for the test
  suite, the YAML input and the Python API.

`install.sh` installs whatever is missing ([doc/INSTALL.md](doc/INSTALL.md)).

**Install time.** About 7 minutes on the M3 laptop above, with the
libraries already installed. Building the serial and MPI versions takes
about 5.5 minutes, and the demo below takes the rest. Installing the
libraries first adds time that depends on the network.

**Demo.** `install.sh` ends with a demo that checks the build. It runs the
regression cases `pbe_atom` (a PBE atom) and `si8_gamma` (an 8-atom Si cell,
solved with the Davidson, RMM-DIIS and CG eigensolvers and with
Car–Parrinello dynamics), first serially and then on two MPI tasks. To run
it again at any time:

```sh
python3 testsuite/run.py pbe_atom si8_gamma
```

*Expected output:*
1. one line per run, giving its energy, iteration count and band gap;
2. `[PASS]` lines comparing the energies and gaps with stored references
   and with each other;
3. a final line, `2/2 regression cases passed`.

*Expected run time* on the laptop above: about 1 minute serially, and
about 40 s for `si8_gamma` on two MPI tasks.

The full suite (`python3 testsuite/run.py`) has 23 cases: functionals,
DFT+U, local and plane-wave hybrids, spin–orbit coupling, forces,
relaxation and symmetry. It takes about 2 hours (2 h 4 min measured) on
the same laptop, without the long runs marked optional.

**Reproducing the results of the method paper** (arXiv:2611.20789).
The input files, scripts and per-system results of every calculation in
the paper are in the archive `inputs.zip`, with a README for each
calculation. The archive comes with the paper's Supplementary
Information and will also be archived on Zenodo.

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
- the variational self-consistency method and the C1/C2 certificates:
  M. Sotoudeh, *Detecting false self-consistency in density-functional
  calculations*, arXiv:2611.20789 (2026),
  [doi:10.48550/arXiv.2611.20789](https://doi.org/10.48550/arXiv.2611.20789);
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
