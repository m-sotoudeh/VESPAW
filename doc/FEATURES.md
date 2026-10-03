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
| Variational SCF, route 1 — `MIXER='VARIATIONAL'` (**default**) | Stable | The Hamiltonian parameters (plane-wave potential, meta-GGA τ-potential, PAW one-center Hamiltonian) are the variables. The free energy is evaluated for the output state of every trial Hamiltonian, and a Kerker-preconditioned Anderson step is accepted only if it does not increase; the uniform-shift gauge mode is projected out. Works with `OPTIMIZER='DAV'`, `'RMM'` and `'CG'`. All three agree with Car–Parrinello friction dynamics within 1 µHa for Si, Al at 1000 K, AFM NiO (PBE, PBE+U, PBE0r, PBE-SL0r, TSCAN), MgO/PBE-SL0r and TSCAN Si (`si_pbe`, `al_hot`, `nio_*`, `mgo_pbesl0r`, `si_tscan`). `MIXTOL` is the mean-square potential residual [H²]. |
| Variational SCF, route 2 — `OPTIMIZER='EDFT' MIXER='VARIATIONAL'` | Experimental | Direct minimization of the Mermin free energy over orbitals and an auxiliary occupation Hamiltonian (ensemble DFT), with a Kerker-preconditioned orbital direction and the same energy-monotone line search. Robust, but slower than route 1 for large or antiferromagnetic cells. Collinear spin only. Runs under MPI. |
| RMM-DIIS warm-up (`RMMDELAY`, default 8) | Stable | RMM-DIIS started from random orbitals converges to wrong states; the first `RMMDELAY` steps therefore use the Davidson solver. Set `RMMDELAY=0` for a warm restart. |
| Pulay / Anderson density mixing (`MIXER='PULAY'`, `'ANDERSON'`) | Stable | Conventional density mixing, kept for comparison. The energy printed is that of the mixed density and is not an upper bound. |

## Functionals and physics

| Feature | Status | Description |
|---|---|---|
| Exchange-correlation functionals | Stable | LDA (PZ, PW92), GGA (PBE, RPBE, PW91, revPBE, PBEsol, BLYP, BP86, …) and meta-GGA (SCAN, SCAN-L), selected by number or name (`TYPE='PBEsol'`, also the VASP tags `PE`, `RP`, `91`, `RE`, `PS`). All agree with libxc and the VASP implementation to machine precision (audit of October 2026). Full list: `!DFT TYPE` in `doc/manual.pdf` and Table 2.1 of `doc/theory/theory.pdf`. |
| Local hybrids PBE0r and PBE-SL0r | Stable | `!NTBO MODUS='HYBRID' HFWEIGHT=...` with `TYPE='PBE'` (PBE0r); add `DCTYPE='SCANL'` for PBE-SL0r (PBE everywhere, onsite exchange SCAN-L mixed with exact exchange; called SCAN0r before October 2026). The onsite Fock term is applied as an operator in all eigensolvers (`nio_pbe0r`, `nio_pbesl0r`, `mgo_pbesl0r`). |
| Local hybrid SCAN0r | Experimental | `TYPE='SCAN' TSCANDELOC=T` with `!NTBO MODUS='HYBRID' HFWEIGHT=...`: orbital-τ SCAN everywhere plus onsite exact exchange (`nio_scan0r_tscan`). |
| SCAN / TSCAN meta-GGA | Experimental | Including the τ-potential in the variational SCF; spin-polarized SCAN correlation on the grid and in the one-center terms (`si_tscan`, `nio_tscan`). |
| DFT+U (Dudarev) | Stable | `!STRUCTURE!SPECIES!DFTU L= U[EV]= J[EV]= !END`. NiO, U−J = 6.3 eV: gap 3.12 eV; Davidson, RMM and friction MD agree (`nio_dftu`). |
| Spin–orbit coupling | Experimental | `!CONTROL!DFT SPINORBIT=T` with non-collinear spinors (`NSPIN=3`). Si: Γ-point splitting 48 meV (exp. 44 meV) (`si_soc`). |
| Dipole correction | Experimental | `!CONTROL!DFT DIPOLECORR=T` (Bengtsson), normal along the third lattice vector; the stress is not corrected (`h2o_dipole`). |
| Exact angular one-center XC | Experimental | `!CONTROL!DFT ANGULARXC='EXACT'` integrates the non-spherical one-center XC on an angular grid instead of the Taylor expansion (about 1.5× the cost per iteration). |
| Explicit partial-wave reference energies | Experimental | `!AUGMENT DEREF[H]=` for HBS setups (replaces the box-state energies). |

## Symmetry and parallelization

| Feature | Status | Description |
|---|---|---|
| Space-group symmetry | Stable | `!KPOINTS SYMMETRY=T`: irreducible k-points; density, one-center density matrices, forces and stress symmetrized. Agrees with the full mesh within 10⁻⁹ H including forces and stress (`symmetry_al`, `symmetry_si`). |
| MPI parallelization of the new solvers | Stable | k-point and plane-wave parallelization for route 1, the DAV/RMM/CG solvers and symmetry. The full regression suite, including route 2 and the density mixers, gives the serial results on 4 MPI tasks (`VESPAW_REGRESSION_NP=4`). |
| Native Apple Silicon build | Stable | `parms.in_use` / `parms.arm64`, about 2.2× faster than x86_64 under Rosetta. |

## Tooling

| Feature | Status | Description |
|---|---|---|
| `vespaw` command-line dispatcher | Stable | One command, dispatches by file extension; see [USAGE.md](USAGE.md). |
| YAML input frontend | Experimental | `vespaw run foo.yaml`; covers the common case, with a `raw:` passthrough for any `.cntl` keyword. |
| Python API (`python/`) | Experimental | Runs a deck and reads HDF5 results as numpy arrays. |
| `potentials` setup library | Experimental | One PAW setup file per element (105 elements). The lattice-constant benchmark is in progress; see [KNOWN_ISSUES.md](KNOWN_ISSUES.md) for setups that are not yet validated. |
| Regression suite (`testsuite/`) | Stable | 18 cases checking energies, gaps, forces, stress and cross-optimizer agreement, in serial and MPI runs. |

## Not yet available

GPU acceleration and machine-learned initial guesses are planned but not
implemented.
