# potentials setup audit (2026-09-27)

Each `*.strc-species` file was built for an isolated atom (PBE, setup construction only) with the
setup diagnostics added on 2026-09-27. Status:

- **UNUSABLE** -- setup construction stops (mostly `OMAT NOT POSITIVE DEFINITE`: the PAW overlap
  operator is not positive definite, i.e. negative-norm states). Do not use.
- **GHOST** -- the atomic PAW Hamiltonian (diagonalized in a Bessel basis) has a bound state far
  below the valence level of that channel. Do not use.
- **INCOMPLETE** -- the local potential binds a state far below the valence level in a channel
  with a single projector. Atomic tests pass, but in solids the variational minimum can collapse
  into such states (NiO: 3-10 Ha below the physical state, growing with the cutoff). A second projector for that channel removed the collapse for Ni.
- **OK** -- no warning. (Not a validation of accuracy.)

**Rebuild (2026-09-27):** every flagged setup has a replacement `<el>_hbs_v1` /
`<el>_hbs_sc_v1` built from CP-PAW's HBS parameters (`parameters/stp.cntl`) with two
projectors per s/p/d channel. All 63 replacements build cleanly (re-audited). Ten needed
a different variant than the original (standard <-> semicore; see each file's header), and
U needed its matching radii scaled by 1.1. The replacements are checked for the isolated
atom only; solid-state validation (lattice constants, cutoff convergence) is listed under
"Validated in solids" below as it is done. Use the replacement, not the flagged setup.

### Validated in solids

PBE equilibrium lattice constants (2026-09-28): 5 points at +-3 %, cubic
E(V) fit, 12x12x12 k-points, Mermin 500 K, DAV+VARIATIONAL. Reference: AE PBE
(Delta-benchmark volumes).

| setup | solid | cutoff | a0 [AA] | AE PBE [AA] | error |
|---|---|---|---|---|---|
| `fe_hbs_v1` (FM) | bcc Fe | 45 Ry | 2.791 | 2.832 | -1.5 % |
| `cu_hbs_v1` | fcc Cu | 45 Ry | 3.594 | 3.632 | -1.0 % |
| `cu_hbs_v1` | fcc Cu | 70 Ry | 3.593 | 3.632 | -1.1 % |
| `cu_hbs_v1`, NPRO=3 3 3 | fcc Cu | 70 Ry | 3.591 | 3.632 | -1.1 % |
| `cu_ndlss_v0` (old) | fcc Cu | 70 Ry | 3.666 | 3.632 | +0.9 % |
| `mo_hbs_sc_v1` (4s4p valence) | bcc Mo | 70 Ry | 3.172 | 3.163 | +0.3 % |
| control: `SI_NDLSS_V0` (regression) | diamond Si | 40 Ry | 5.490 | 5.469 | +0.4 % |

Conclusions:
- The code is fine (Si, Mo within normal PAW accuracy).
- The frozen-semicore 3d setups (Fe, Cu with ZV=8/11) are about 1 % off, and
  the error depends on the construction (old NDLSS +0.9 %, HBS -1.1 %), not
  on cutoff or projector count. **Do not use them for lattice constants,
  phonons or elastic properties.** They are fine for electronic structure
  (no ghosts, no collapse).
- The semicore setup of Mo is accurate.

Semicore 3d setups (2026-09-28, after the construction fix 9672667; all ten
3d `*_HBS_SC` parameter sets now build without ghost/incompleteness
warnings):

| setup | RCL/RCOV | suggested cutoff | a0 [AA] | error |
|---|---|---|---|---|
| CU_HBS_SC (stock) | 0.5 | 137 Ry | 3.649 (110 Ry) | +0.46 % |
| CU_HBS_SC, radii 0.8 | 0.8 | 43 Ry | 3.697 (45/60 Ry) | +1.8 % |

The error grows with the augmentation radius (frozen-core HBS 0.75: -1.1 %,
old NDLSS 0.7: +0.9 %, semicore 0.8: +1.8 %, semicore 0.5: +0.46 %). VASP's
Cu_pv reaches ~0.2 % with a larger radius at 27 Ry, so the construction,
not PAW itself, is too sensitive to the radius. OPEN: sweep RCL, the local
potential radius and the core pseudization radius separately to find the
term that drives the error. Until then: CU_HBS_SC stock for structural
accuracy (expensive), `cu_hbs_v1` for electronic structure.

Summary: 42 OK, 16 incomplete, 5 ghost, 42 unusable (of 105).

| setup | status | detail | replacement (2026-09-27 rebuild) |
|---|---|---|---|
| `ba_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `ba_hbs_sc_v1` (BA_HBS_SC, NPRO=2 2 2 1) |
| `bi_ndlss_sc_v0` | UNUSABLE | inaccuracy while unscreening ps potential | `bi_hbs_v1` (BI_HBS, NPRO=2 2 2 1) |
| `bk_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `bk_hbs_sc_v1` (BK_HBS_SC, NPRO=2 2 2 1) |
| `ce_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `ce_hbs_sc_v1` (CE_HBS_SC, NPRO=2 2 2 1) |
| `cf_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `cf_hbs_sc_v1` (CF_HBS_SC, NPRO=2 2 2 1) |
| `cs_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `cs_hbs_sc_v1` (CS_HBS_SC, NPRO=2 2 2) |
| `dy_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `dy_hbs_sc_v1` (DY_HBS_SC, NPRO=2 2 2 1) |
| `er_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `er_hbs_sc_v1` (ER_HBS_SC, NPRO=2 2 2 1) |
| `es_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `es_hbs_sc_v1` (ES_HBS_SC, NPRO=2 2 2 1) |
| `eu_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `eu_hbs_sc_v1` (EU_HBS_SC, NPRO=2 2 2 1) |
| `f_ndlss_v0` | UNUSABLE | omat not positive definite | `f_hbs_v1` (F_HBS, NPRO=2 2 2) |
| `fm_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `fm_hbs_sc_v1` (FM_HBS_SC, NPRO=2 2 2 1) |
| `gd_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `gd_hbs_sc_v1` (GD_HBS_SC, NPRO=2 2 2 1) |
| `hf_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `hf_hbs_sc_v1` (HF_HBS_SC, NPRO=2 2 2 1) |
| `ho_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `ho_hbs_sc_v1` (HO_HBS_SC, NPRO=2 2 2 1) |
| `lr_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `lr_hbs_sc_v1` (LR_HBS_SC, NPRO=2 2 2 1) |
| `lu_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `lu_hbs_sc_v1` (LU_HBS_SC, NPRO=2 2 2 1) |
| `md_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `md_hbs_sc_v1` (MD_HBS_SC, NPRO=2 2 2 1) |
| `mo_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `mo_hbs_sc_v1` (MO_HBS_SC, NPRO=2 2 2) |
| `na_ndlss_sc_v0` | UNUSABLE | inaccuracy while unscreening ps potential | `na_hbs_sc_v1` (NA_HBS_SC, NPRO=2 2 2) |
| `nb_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `nb_hbs_sc_v1` (NB_HBS_SC, NPRO=2 2 2) |
| `nd_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `nd_hbs_sc_v1` (ND_HBS_SC, NPRO=2 2 2 1) |
| `np_ndlss_v0` | UNUSABLE | omat not positive definite | `np_hbs_sc_v1` (NP_HBS_SC, NPRO=2 2 2 1) |
| `os_ndlss_v0` | UNUSABLE | omat not positive definite | `os_hbs_v1` (OS_HBS, NPRO=2 2 2 1) |
| `pa_ndlss_v0` | UNUSABLE | omat not positive definite | `pa_hbs_sc_v1` (PA_HBS_SC, NPRO=2 2 2 1) |
| `pm_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `pm_hbs_sc_v1` (PM_HBS_SC, NPRO=2 2 2 1) |
| `po_ndlss_sc_v0` | UNUSABLE | bisection loop not converged | `po_hbs_v1` (PO_HBS, NPRO=2 2 2 1) |
| `pr_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `pr_hbs_sc_v1` (PR_HBS_SC, NPRO=2 2 2 1) |
| `pu_ndlss_v0` | UNUSABLE | omat not positive definite | `pu_hbs_sc_v1` (PU_HBS_SC, NPRO=2 2 2 1) |
| `ra_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `ra_hbs_sc_v1` (RA_HBS_SC, NPRO=2 2 2 1) |
| `re_ndlss_v0` | UNUSABLE | omat not positive definite | `re_hbs_v1` (RE_HBS, NPRO=2 2 2 1) |
| `rf_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `rf_hbs_sc_v1` (RF_HBS_SC, NPRO=2 2 2 1) |
| `rh_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `rh_hbs_sc_v1` (RH_HBS_SC, NPRO=2 2 2) |
| `ru_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `ru_hbs_sc_v1` (RU_HBS_SC, NPRO=2 2 2) |
| `sm_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `sm_hbs_sc_v1` (SM_HBS_SC, NPRO=2 2 2 1) |
| `tb_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `tb_hbs_sc_v1` (TB_HBS_SC, NPRO=2 2 2 1) |
| `tc_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `tc_hbs_sc_v1` (TC_HBS_SC, NPRO=2 2 2) |
| `th_ndlss_v0` | UNUSABLE | omat not positive definite | `th_hbs_sc_v1` (TH_HBS_SC, NPRO=2 2 2 1) |
| `tm_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `tm_hbs_sc_v1` (TM_HBS_SC, NPRO=2 2 2 1) |
| `v_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `v_hbs_sc_v1` (V_HBS_SC, NPRO=2 2 2) |
| `y_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `y_hbs_sc_v1` (Y_HBS_SC, NPRO=2 2 2) |
| `zr_ndlss_sc_v0` | UNUSABLE | omat not positive definite | `zr_hbs_sc_v1` (ZR_HBS_SC, NPRO=2 2 2) |
| `ac_ndlss_v0` | GHOST | ghost in l=0 at -1.1737 H (valence -0.1405 H) | `ac_hbs_sc_v1` (AC_HBS_SC, NPRO=2 2 1 1) |
| `cm_ndlss_sc_v0` | GHOST | ghost in l=0 at -0.5589 H (valence -0.1654 H) | `cm_hbs_sc_v1` (CM_HBS_SC, NPRO=2 2 2 1) |
| `ha_ndlss_sc_v0` | GHOST | ghost in l=0 at -0.7788 H (valence -0.2362 H) | `ha_hbs_sc_v1` (HA_HBS_SC, NPRO=2 2 2 1) |
| `la_ndlss_v0` | GHOST | ghost in l=0 at -0.4886 H (valence -0.1310 H) | `la_hbs_sc_v1` (LA_HBS_SC, NPRO=2 2 2 1) |
| `u_ndlss_v0` | GHOST | ghost in l=0 at -2.4643 H (valence -0.1536 H) | `u_hbs_sc_v1` (U_HBS_SC RCLx1.1, NPRO=2 2 2 1) |
| `ag_ndlss_v0` | INCOMPLETE | l=0: local-potential state -0.9571 H vs valence -0.1675 H | `ag_hbs_v1` (AG_HBS, NPRO=2 2 2) |
| `as_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.1187 H vs valence -0.5370 H | `as_hbs_v1` (AS_HBS, NPRO=2 2 2) |
| `at_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.3444 H vs valence -0.7152 H | `at_hbs_v1` (AT_HBS, NPRO=2 2 2 1) |
| `au_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.2174 H vs valence -0.2173 H | `au_hbs_v1` (AU_HBS, NPRO=2 2 2 1) |
| `br_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.4133 H vs valence -0.7449 H | `br_hbs_v1` (BR_HBS, NPRO=2 2 2) |
| `co_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.0110 H vs valence -0.2056 H | `co_hbs_v1` (CO_HBS, NPRO=2 2 2) |
| `cu_ndlss_v0` | INCOMPLETE | l=0: local-potential state -0.7858 H vs valence -0.1750 H | `cu_hbs_v1` (CU_HBS, NPRO=2 2 2) |
| `fe_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.0422 H vs valence -0.1986 H | `fe_hbs_v1` (FE_HBS, NPRO=2 2 2) |
| `hg_ndlss_v0` | INCOMPLETE | l=0: local-potential state -0.8877 H vs valence -0.2558 H | `hg_hbs_v1` (HG_HBS, NPRO=2 2 2 1) |
| `i_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.4231 H vs valence -0.6442 H | `i_hbs_v1` (I_HBS, NPRO=2 2 2) |
| `ir_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.6099 H vs valence -0.2372 H | `ir_hbs_v1` (IR_HBS, NPRO=2 2 2 1) |
| `ni_ndlss_v0` | INCOMPLETE | l=0: local-potential state -0.9812 H vs valence -0.2123 H | `ni_hbs_v1` (NI_HBS, NPRO=2 2 2) |
| `pt_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.3871 H vs valence -0.2133 H | `pt_hbs_sc_v1` (PT_HBS_SC, NPRO=2 2 2 1) |
| `sb_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.1410 H vs valence -0.4752 H | `sb_hbs_v1` (SB_HBS, NPRO=2 2 2) |
| `se_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.2917 H vs valence -0.6401 H | `se_hbs_v1` (SE_HBS, NPRO=2 2 2) |
| `te_ndlss_v0` | INCOMPLETE | l=0: local-potential state -1.2914 H vs valence -0.5594 H | `te_hbs_v1` (TE_HBS, NPRO=2 2 2) |
| `al_ndlss_v0` | OK |  |  |
| `am_ndlss_sc_v0` | OK |  |  |
| `ar_ndlss_v0` | OK |  |  |
| `b_ndlss_v0` | OK |  |  |
| `be_ndlss_v0` | OK |  |  |
| `c_ndlss_v0` | OK |  |  |
| `ca_ndlss_sc_v0` | OK |  |  |
| `cd_ndlss_v0` | OK |  |  |
| `cl_ndlss_v0` | OK |  |  |
| `cr_ndlss_sc_v0` | OK |  |  |
| `fr_ndlss_sc_v0` | OK |  |  |
| `ga_ndlss_sc_v0` | OK |  |  |
| `ge_ndlss_sc_v0` | OK |  |  |
| `h_ndlss_v0` | OK |  |  |
| `he_ndlss_v0` | OK |  |  |
| `in_ndlss_sc_v0` | OK |  |  |
| `k_ndlss_sc_v0` | OK |  |  |
| `kr_ndlss_v0` | OK |  |  |
| `li_ndlss_v0` | OK |  |  |
| `mg_ndlss_v0` | OK |  |  |
| `mn_ndlss_sc_v0` | OK |  |  |
| `n_ndlss_v0` | OK |  |  |
| `ne_ndlss_v0` | OK |  |  |
| `no_ndlss_sc_v0` | OK |  |  |
| `o_ndlss_v0` | OK |  |  |
| `p_ndlss_v0` | OK |  |  |
| `pb_ndlss_sc_v0` | OK |  |  |
| `pd_ndlss_v0` | OK |  |  |
| `rb_ndlss_sc_v0` | OK |  |  |
| `rn_ndlss_v0` | OK |  |  |
| `s_ndlss_v0` | OK |  |  |
| `sc_ndlss_sc_v0` | OK |  |  |
| `si_ndlss_v0` | OK |  |  |
| `sn_ndlss_sc_v0` | OK |  |  |
| `sr_ndlss_sc_v0` | OK |  |  |
| `ta_ndlss_sc_v0` | OK |  |  |
| `ti_ndlss_sc_v0` | OK |  |  |
| `tl_ndlss_sc_v0` | OK |  |  |
| `w_ndlss_sc_v0` | OK |  |  |
| `xe_ndlss_v0` | OK |  |  |
| `yb_ndlss_sc_v0` | OK |  |  |
| `zn_ndlss_v0` | OK |  |  |
