# Changelog

## Unreleased

- Plane-wave hybrid functionals: nonlocal PAW Fock exchange for PBE0,
  HSE06, SCAN0 and r²SCAN0, block `!CONTROL!DFT!FOCK` (`AEXX`,
  `HFSCREEN` or `MU`, `CORE`, `ACE`, `RCOMP`). Pair densities in plane
  waves with PAW compensation charges, one-center and frozen-core
  core–valence exchange, erfc screening with the short-range PBE
  exchange hole for HSE06, the Born–von-Kármán treatment of the Coulomb
  singularity, adaptively compressed exchange (19 times faster for
  Si₈), k-points, space-group symmetry, collinear and non-collinear
  spin, forces and stress, MPI over k-point groups. Works with all
  eigensolvers and the variational mixer.
- Validated against VASP 6.4.2: the band-gap shift from PBE to PBE0 and
  HSE06 agrees within 0.3–2 % for Si, C, MgO and MnO, the PBE0
  FM–AFM energy of MnO within 1 meV and of NiO within 10 meV. With the
  HBS Ni setup the NiO shift is 9 % smaller than with the VASP Ni
  POTCAR (see `doc/KNOWN_ISSUES.md`).
- New regression case `si_hybrid` (PBE0 full mesh, symmetry and RMM-DIIS,
  HSE06, finite-difference force) and tutorial example
  `examples/tutorial/07_si_hybrid`.
- MPI: fixed the creation of a second communicator split from the same
  parent, which dropped the first one.
- Not yet available with hybrids: `GAMMASUPER=T`, the onsite hybrids
  (`!NTBO`) in the same run, B3LYP.

## 1.0.3 — 2026-11-21

- Onsite hybrids: the swap term of PBE-SL0r is now applied also at
  `a_x=0`. Before, it was skipped when the admixture was zero, which
  made the PBE-SL0r total energy jump between `a_x=0` and small `a_x`
  (MnO by 12 eV).
  Results of PBE-SL0r at `a_x=0` change. All other results are unchanged.
- Onsite hybrids with MPI: the onsite energy loop is split over the MPI
  tasks and the transformation to the tight-binding orbitals is cached
  per k-point. An 80-atom CaMnO₃ supercell with SCAN0r runs 1.5 times
  faster per iteration on 48 cores, with energies equal to within
  5×10⁻¹⁰ Ha.
- Diagnostic environment variable `VESPAW_DIAG_NOSWAP` (any nonempty
  value) omits the swap term.
- Installation: Open MPI transport settings for MPI runs on cluster
  nodes (`OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm`), see `doc/INSTALL.md`.
- Theory book: foundations, line search per mixer, certificates C1 and
  C2, L-BFGS, DFT+U, symmetry, dipole correction, relativistic setups and
  spin-orbit coupling, MPI and the full PAW force formula.
- The method paper is cited as arXiv:2611.20789. CI runs on macos-15.

## 1.0.2 — 2026-10-11

- Test suite: the runs behind Table I of the method paper (direct
  minimization for most systems, CG for Si with TSCAN, NiO with PBE+U
  and PBE-SL0r, and MgO, RMM-DIIS and CG for the water layer, friction
  dynamics for Al at 1000 K, direct minimization with symmetry) are part
  of the regression cases. They are marked `"optional"` because they take long;
  they run when named (`python3 testsuite/run.py si8_gamma:edft_var`) or
  with `VESPAW_REGRESSION_OPTIONAL=1`.
- The code and all results are unchanged from 1.0.1.

## 1.0.1 — 2026-10-11

First release archived on Zenodo, which assigns it a DOI. The code and
all results are unchanged from 1.0.0.

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
