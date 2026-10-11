# What VESPAW adds to CP-PAW

VESPAW is a development of CP-PAW (Peter E. Blöchl, Clausthal University
of Technology). The PAW engine, the input format and most of the physics
are those of CP-PAW; see `doc/manual.pdf` for them. This page lists what
VESPAW adds, each with a status:

- **Stable** — validated against independent references, in the
  regression suite, safe to use.
- **Experimental** — working and validated on the listed cases, but
  still changing or not yet broadly tested.
- **Prototype** — a research vehicle, not yet a finished feature.

All numbers below are reproduced by `python3 testsuite/run.py`
(case names in parentheses). Energies refer to scalar-relativistic
setups without ZORA, the default since 2026-09-29.

## Self-consistency

| Feature | Status | Description |
|---|---|---|
| Hamiltonian minimization (HM) — `MIXER='VARIATIONAL'` (**default**) | Stable | The Hamiltonian parameters (plane-wave potential, meta-GGA τ-potential, PAW one-center Hamiltonian) are the variables. The free energy is evaluated for the output state of every trial Hamiltonian, and a Kerker-preconditioned Anderson step is accepted only if it does not increase; the uniform-shift gauge mode is projected out. Works with `OPTIMIZER='DAV'`, `'RMM'` and `'CG'`. All three agree with Car–Parrinello friction dynamics within 1 µHa for Si, Al at 1000 K, AFM NiO (PBE, PBE+U, PBE0r, PBE-SL0r, TSCAN), MgO/PBE-SL0r and TSCAN Si (`si_pbe`, `al_hot`, `nio_*`, `mgo_pbesl0r`, `si_tscan`). `MIXTOL` is the mean-square potential residual [H²]. |
| Direct minimization (DM) — `OPTIMIZER='EDFT' MIXER='VARIATIONAL'` | Experimental | Direct minimization of the Mermin free energy over orbitals and an auxiliary occupation Hamiltonian (ensemble DFT), with a Kerker-preconditioned orbital direction and the same energy-monotone line search. Robust, but slower than HM for large or antiferromagnetic cells. Collinear spin only. Runs under MPI. |
| RMM-DIIS warm-up (`RMMDELAY`, default 8) | Stable | RMM-DIIS started from random orbitals converges to wrong states; the first `RMMDELAY` steps therefore use the Davidson solver. Set `RMMDELAY=0` for a warm restart. |
| Pulay / Anderson density mixing (`MIXER='PULAY'`, `'ANDERSON'`) | Stable | Conventional density mixing, kept for comparison. The energy printed is that of the mixed density and is not an upper bound. |

## Functionals and physics

| Feature | Status | Description |
|---|---|---|
| Exchange-correlation functionals | Stable | LDA (PZ, PW92), GGA (PBE, RPBE, PW91, revPBE, PBEsol, BLYP, BP86, …) and meta-GGA (SCAN, r²SCAN, SCAN-L), selected by number or name (`TYPE='PBEsol'`, also the VASP tags `PE`, `RP`, `91`, `RE`, `PS`). Checked pointwise against libxc to machine precision and against VASP (audit of October 2026). Full list: `!DFT TYPE` in the manual and Table 2.1 of the theory book. |
| Local hybrids PBE0r and PBE-SL0r | Stable | `!NTBO MODUS='HYBRID' HFWEIGHT=...` with `TYPE='PBE'` (PBE0r); add `DCTYPE='SCANL'` for PBE-SL0r (PBE everywhere, onsite exchange SCAN-L mixed with exact exchange; called SCAN0r before October 2026). The onsite Fock term is applied as an operator in all eigensolvers (`nio_pbe0r`, `nio_pbesl0r`, `mgo_pbesl0r`). |
| Local hybrids SCAN0r and r²SCAN0r | Experimental | `TYPE='SCAN'` or `'R2SCAN'` with `TSCANDELOC=T` and `!NTBO MODUS='HYBRID' HFWEIGHT=...`: orbital-τ meta-GGA everywhere plus onsite exact exchange. HM, RMM-DIIS and friction dynamics agree for NiO (`nio_scan0r_tscan`). |
| SCAN and r²SCAN meta-GGA (orbital τ) | Experimental | `TYPE='SCAN'` or `TYPE='R2SCAN'` with `TSCANDELOC=T`. Including the τ-potential in the variational SCF; spin-polarized correlation on the grid and in the one-center terms; setups built with PBE and the pseudized core τ in the plane-wave part. SCAN and r²SCAN band gaps of MgO agree with VASP within 0.06 eV (`si_tscan`, `si_r2scan`, `nio_tscan`). |
| DFT+U (Dudarev) | Stable | `!STRUCTURE!SPECIES!DFTU L= U[EV]= J[EV]= !END`. NiO, U−J = 6.3 eV: gap 3.12 eV; Davidson, RMM and friction MD agree (`nio_dftu`). |
| Spin–orbit coupling | Experimental | `!CONTROL!DFT SPINORBIT=T` with non-collinear spinors (`NSPIN=3`). Si: Γ-point splitting 48 meV (exp. 44 meV) (`si_soc`). |
| Dipole correction | Experimental | `!CONTROL!DFT DIPOLECORR=T` (Bengtsson), normal along the third lattice vector; the stress is not corrected (`h2o_dipole`). |
| Exact angular one-center XC | Experimental | `!CONTROL!DFT ANGULARXC='EXACT'` integrates the non-spherical one-center XC on an angular grid instead of the Taylor expansion (about 1.5× the cost per iteration). |
| Explicit partial-wave reference energies | Experimental | `!AUGMENT DEREF[H]=` for HBS setups (replaces the box-state energies). |

## Geometry optimization

| Feature | Status | Description |
|---|---|---|
| L-BFGS relaxation of atoms — `!RDYN OPTIMIZER='LBFGS'` | Experimental | Quasi-Newton optimization (limited-memory BFGS) as an alternative to damped dynamics. A step is taken only after the SCF has converged at the current geometry; the SCF tolerances are tightened automatically as the forces become small. Safeguards: backtracking when the energy rises, curvature and descent checks, adaptive trust radius (`MAXSTEP[AA]`). Converges on `FMAX[EV/AA]`. Needs `OPTIMIZER='DAV'`, `'RMM'` or `'CG'`; not yet with `!CONSTRAINTS`. Si: same structure and energy as friction dynamics, about 10× faster (`si_relax`); H2O: 7 geometries to 0.004 eV/Å. |
| Joint relaxation of atoms and cell | Experimental | `!CELL MOVE=T` with `!RDYN OPTIMIZER='LBFGS'` (or `!CELL OPTIMIZER='LBFGS'` alone at fixed fractional coordinates); `CONSTRAINTTYPE` `'FREE'`, `'NOSHEAR'`, `'ISOTROPIC'`; converges on `SMAX[GPA]`. Si: volume within 10⁻⁵ of Parrinello–Rahman friction dynamics (`si_relax`). |
| External pressure and stress | Experimental | `!CELL P[GPA]=` and the general tensor `STRESS[GPA]=` (also in Hartree/bohr³ as `P`, `STRESS`). The load acts on the actual cell: at equilibrium the printed stress equals P·1 + STRESS, for L-BFGS and for Parrinello–Rahman dynamics (`si_relax`). |

## Symmetry and parallelization

| Feature | Status | Description |
|---|---|---|
| Space-group symmetry | Stable | `!KPOINTS SYMMETRY=T`: irreducible k-points; density, one-center density matrices, forces and stress symmetrized. Agrees with the full mesh within 10⁻⁹ H including forces and stress (`symmetry_al`, `symmetry_si`). |
| MPI parallelization of the new solvers | Stable | k-point and plane-wave parallelization for HM, the DAV/RMM/CG solvers and symmetry. The full regression suite, including DM and the density mixers, gives the serial results on 4 MPI tasks (`VESPAW_REGRESSION_NP=4`). |
| Native Apple Silicon build | Stable | parameter file `parms.arm64`, about 2.2× faster than x86_64 under Rosetta. |

## Tooling

| Feature | Status | Description |
|---|---|---|
| `vespaw` command-line dispatcher | Stable | One command, dispatches by file extension; see [USAGE.md](USAGE.md). |
| YAML input frontend | Experimental | `vespaw run foo.yaml`; covers the common case, with a `raw:` passthrough for any `.cntl` keyword. |
| Python API (`python/`) | Experimental | Runs a deck and reads HDF5 results as numpy arrays. |
| `potentials` setup library | Experimental | One PAW setup file per element (105 elements). The lattice-constant benchmark is in progress; see [KNOWN_ISSUES.md](KNOWN_ISSUES.md) for setups that are not yet validated. |
| Regression suite (`testsuite/`) | Stable | 22 cases checking energies, gaps, forces, stress and cross-optimizer agreement, in serial and MPI runs. |

## Not yet available

GPU acceleration and machine-learned initial guesses are planned but not
implemented.
