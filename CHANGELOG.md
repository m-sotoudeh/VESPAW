# Changelog

## 1.0.0 — 2026-10-11, first public release

VESPAW is derived from CP-PAW (see README). Compared with CP-PAW:

- Variational self-consistency: Hamiltonian minimization (HM), with the
  Kohn–Sham Hamiltonian as the variable (default `MIXER='VARIATIONAL'`),
  with Davidson, RMM-DIIS and CG eigensolvers; direct minimization (DM)
  over orbitals and occupations (`OPTIMIZER='EDFT'`).
- Exchange–correlation functionals checked against libxc and VASP;
  PBEsol, true PW91, SCAN and r²SCAN with orbital τ and spin-polarized
  correlation, SCAN-L; local hybrids PBE0r, PBE-SL0r, SCAN0r and r²SCAN0r.
- Meta-GGA PAW treatment: setups built with PBE, pseudized core τ in the
  plane-wave part, one-center τ checked against the kinetic energy;
  SCAN and r²SCAN band gaps agree with VASP within 0.06 eV.
- Density mixers (Pulay, Anderson): residual measured on the density
  cutoff sphere; convergence requires small residuals of the density and
  of the one-center density matrices and a small change of the energy.
- DFT+U, spin–orbit coupling, dipole correction, space-group symmetry.
- MPI for all new solvers and mixers.
- Structured report `<name>.prot` and output in eV/Å/GPa; results in
  HDF5; YAML input and a Python API.
- PAW setup library (`potentials/`), one-step installer (`install.sh`),
  regression test suite (`testsuite/`).
- L-BFGS geometry optimization: `!RDYN OPTIMIZER='LBFGS'` for the atoms,
  jointly with the unit cell (`!CELL MOVE=T`), or the cell alone
  (`!CELL OPTIMIZER='LBFGS'`). Safeguards against energy rises, bad
  curvature and too long steps; SCF tolerances tightened automatically
  near convergence. Validated against friction dynamics (`si_relax`).
- `!CELL STRESS=` (general external stress tensor) is available again,
  with `STRESS[GPA]=` and `P[GPA]=`; at equilibrium the printed stress
  equals P·1 + STRESS, also for Parrinello–Rahman dynamics.
- `!CELL FRIC` is optional (default 0) as documented.
