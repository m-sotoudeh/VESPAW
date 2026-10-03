# Changelog

## 1.0.0 — 2026-10-12, first public release

VESPAW is derived from CP-PAW (see README). Compared with CP-PAW:

- Variational self-consistency: the Kohn–Sham Hamiltonian as the
  variable (route 1, default `MIXER='VARIATIONAL'`) with Davidson,
  RMM-DIIS and CG eigensolvers; direct minimization over orbitals and
  occupations (route 2, `OPTIMIZER='EDFT'`).
- Exchange–correlation functionals checked against libxc and VASP;
  PBEsol, true PW91, SCAN with orbital τ and spin-polarized correlation,
  SCAN-L; local hybrids PBE0r, PBE-SL0r and SCAN0r (experimental).
- DFT+U, spin–orbit coupling, dipole correction, space-group symmetry.
- MPI for all new solvers and mixers.
- Structured report `<name>.prot` and output in eV/Å/GPa; results in
  HDF5; YAML input and a Python API.
- PAW setup library (`potentials/`), one-step installer (`install.sh`),
  regression test suite (`testsuite/`).
