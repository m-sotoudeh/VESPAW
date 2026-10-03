# Known issues

Open problems and limitations of the current version, with workarounds
where they exist. Please report anything not listed here as an issue.

## Correctness and robustness

- **Cells with large vacuum.** The variational optimizer realigns the
  potential gauge and fixes the potential in empty space, which makes
  slabs and molecules converge. A residual floor from the gradient-
  corrected XC potential in the density tail remains: for a water layer
  in an 18 Å cell (`h2o_dipole`) the Davidson run stops about 25 µHa
  above friction dynamics. If a molecule or slab stalls, a softer Kerker
  filter (`MIXG0SQ=0.1`) helps; for µHa accuracy, finish with
  `OPTIMIZER='MD'`.
- **PBE-SL0r with RMM-DIIS can miss the lowest empty band.** For MgO the
  energy is correct, but the gap comes out 6.54 instead of 4.23 eV
  (`mgo_pbesl0r`). Check gaps from this combination with Davidson.
- **PBE-SL0r ghost states at unphysical densities.** Far from
  self-consistency, the one-center SCAN-L double counting can produce
  deep spurious eigenstates (−70 to −100 eV). Start PBE-SL0r from a
  converged PBE or PBE0r restart.
- **Cold-started antiferromagnets.** Without a symmetry-breaking seed,
  NiO can converge to the wrong spin order, and with some Ni setups to a
  state below the physical minimum. Seed the spin order with an orbital
  potential (`!ORBPOT`, then `!ORBPOT_OFF`), as in `testsuite/nio_afm`.
- **Delocalized meta-GGA (`TSCANDELOC=T`).** Cold starts can become
  numerically unstable, and the restart is only readable by another
  `TSCANDELOC=T` run. The onsite PBE-SL0r recipe is better validated.
- **`GAMMASUPER=T` is not correct yet.** It converges to a wrong
  energy; do not use it.
- **Dipole correction:** the stress is not corrected.
- **Bit-level reproducibility.** Repeated runs can differ by ~3×10⁻¹⁰ H
  (vectorization order). Compare results with a tolerance, not bitwise.

## PAW setup library (`potentials`)

A lattice-constant benchmark of 23 elemental solids (PBE, 60 Ry,
compared with all-electron values) is in progress:

- within 1%: Al, Ag, C, Ca, Cu (−0.7%, cause under investigation), Ir,
  K, Li, Mo, Na, Nb, Si, Sr, V;
- too large an error: Pd `pd_ndlss_v0` (+12.6%), Sn `sn_ndlss_sc_v0`
  (−7.5%), Ge `ge_ndlss_sc_v0` (+2.6%), Au `au_hbs_v1` (+1.3%);
- failing: Pt `pt_hbs_sc_v1` and Rh `rh_hbs_sc_v1` (setup construction
  stops), Ta and W `*_ndlss_sc_v0` (spurious deep states), Ba
  `ba_hbs_sc_v1` (SCF does not converge).

Validate a setup for your system before production use; the setups in
`testsuite/` are validated.

## Parallelization

- All SCF routes run under MPI (k-point and plane-wave parallelization):
  route 1 with DAV/RMM/CG, route 2 (`OPTIMIZER='EDFT' MIXER='VARIATIONAL'`)
  and the Pulay/Anderson density mixers. Only the older, non-variational
  EDFT prototype (`OPTIMIZER='EDFT'` with a density mixer) is serial only.
- With MPI, set `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`
  unless you deliberately combine MPI and threads.

## Performance

- **Plane waves at Γ are not halved.** Time-reversal symmetry is
  recognized but not used to store half the coefficients (as VASP's
  `vasp_gam` does). Together with other open items this makes VESPAW
  noticeably slower than VASP for comparable systems.
- **Route 2 is slower than route 1** for large or antiferromagnetic
  cells (64-atom Si: 253 vs 19 iterations). Prefer route 1 there.

## Design limitations

- **RMM-DIIS needs a warm start.** Started from random orbitals it
  converges to wrong eigenstates, as residual minimization does in
  general. The first `RMMDELAY=8` steps therefore use Davidson; set
  `RMMDELAY=0` only when restarting from converged orbitals.
- **Density-mixed energies are not upper bounds.** With
  `MIXER='PULAY'` or `'ANDERSON'` the printed energy belongs to the
  mixed density and can lie below the true minimum. Use the default
  variational mixer when energies must be compared at the µHa level.
