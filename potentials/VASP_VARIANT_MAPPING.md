# Standard-vs-semicore variant selection

Per-element choice of CP-PAW's Element0ktest standard vs. semicore (`_SC`)
variant, made to roughly match VASP's typical recommended PAW-PBE POTCAR
default (plain vs. `_sv`/`_pv`/`_d`). **This is a best-effort mapping from
general knowledge of VASP's recommended-potentials conventions, cross-checked
against which suffixes actually exist in
a local copy of the VASP PBE.54 POTCAR library where
possible — it is NOT an independently verified, authoritative lookup.**
Correct individual entries here if you find they're wrong; each `.strc-species`
file's header points back to this table.

`only` = element has just one CP-PAW variant (no choice made).
`†` = no non-GW VASP PBE.54 potential exists for this element at all, so the
choice has no VASP reference point (defaulted by periodic-table pattern).

| Z | Element | Variant | File |
|---|---------|---------|------|
| 1 | H | only | `h_ndlss_v0.strc-species` |
| 2 | He | only | `he_ndlss_v0.strc-species` |
| 3 | Li | semicore | `li_ndlss_sc_v0.strc-species` |
| 4 | Be | standard | `be_ndlss_v0.strc-species` |
| 5 | B | only | `b_ndlss_v0.strc-species` |
| 6 | C | only | `c_ndlss_v0.strc-species` |
| 7 | N | only | `n_ndlss_v0.strc-species` |
| 8 | O | only | `o_ndlss_v0.strc-species` |
| 9 | F | only | `f_ndlss_v0.strc-species` |
| 10 | Ne | only | `ne_ndlss_v0.strc-species` |
| 11 | Na | semicore | `na_ndlss_sc_v0.strc-species` |
| 12 | Mg | standard | `mg_ndlss_v0.strc-species` |
| 13 | Al | only | `al_ndlss_v0.strc-species` |
| 14 | Si | only | `si_ndlss_v0.strc-species` |
| 15 | P | only | `p_ndlss_v0.strc-species` |
| 16 | S | only | `s_ndlss_v0.strc-species` |
| 17 | Cl | only | `cl_ndlss_v0.strc-species` |
| 18 | Ar | only | `ar_ndlss_v0.strc-species` |
| 19 | K | semicore | `k_ndlss_sc_v0.strc-species` |
| 20 | Ca | semicore | `ca_ndlss_sc_v0.strc-species` |
| 21 | Sc | semicore | `sc_ndlss_sc_v0.strc-species` |
| 22 | Ti | semicore | `ti_ndlss_sc_v0.strc-species` |
| 23 | V | semicore | `v_ndlss_sc_v0.strc-species` |
| 24 | Cr | semicore | `cr_ndlss_sc_v0.strc-species` |
| 25 | Mn | semicore | `mn_ndlss_sc_v0.strc-species` |
| 26 | Fe | standard | `fe_ndlss_v0.strc-species` |
| 27 | Co | standard | `co_ndlss_v0.strc-species` |
| 28 | Ni | standard | `ni_ndlss_v0.strc-species` |
| 29 | Cu | standard | `cu_ndlss_v0.strc-species` |
| 30 | Zn | standard | `zn_ndlss_v0.strc-species` |
| 31 | Ga | semicore | `ga_ndlss_sc_v0.strc-species` |
| 32 | Ge | semicore | `ge_ndlss_sc_v0.strc-species` |
| 33 | As | standard | `as_ndlss_v0.strc-species` |
| 34 | Se | standard | `se_ndlss_v0.strc-species` |
| 35 | Br | standard | `br_ndlss_v0.strc-species` |
| 36 | Kr | only | `kr_ndlss_v0.strc-species` |
| 37 | Rb | semicore | `rb_ndlss_sc_v0.strc-species` |
| 38 | Sr | semicore | `sr_ndlss_sc_v0.strc-species` |
| 39 | Y | semicore | `y_ndlss_sc_v0.strc-species` |
| 40 | Zr | semicore | `zr_ndlss_sc_v0.strc-species` |
| 41 | Nb | semicore | `nb_ndlss_sc_v0.strc-species` |
| 42 | Mo | semicore | `mo_ndlss_sc_v0.strc-species` |
| 43 | Tc | semicore | `tc_ndlss_sc_v0.strc-species` |
| 44 | Ru | semicore | `ru_ndlss_sc_v0.strc-species` |
| 45 | Rh | semicore | `rh_ndlss_sc_v0.strc-species` |
| 46 | Pd | standard | `pd_ndlss_v0.strc-species` |
| 47 | Ag | standard | `ag_ndlss_v0.strc-species` |
| 48 | Cd | standard | `cd_ndlss_v0.strc-species` |
| 49 | In | semicore | `in_ndlss_sc_v0.strc-species` |
| 50 | Sn | semicore | `sn_ndlss_sc_v0.strc-species` |
| 51 | Sb | standard | `sb_ndlss_v0.strc-species` |
| 52 | Te | standard | `te_ndlss_v0.strc-species` |
| 53 | I | standard | `i_ndlss_v0.strc-species` |
| 54 | Xe | standard | `xe_ndlss_v0.strc-species` |
| 55 | Cs | semicore | `cs_ndlss_sc_v0.strc-species` |
| 56 | Ba | semicore | `ba_ndlss_sc_v0.strc-species` |
| 57 | La | standard | `la_ndlss_v0.strc-species` |
| 58 | Ce | semicore | `ce_ndlss_sc_v0.strc-species` |
| 59 | Pr | semicore | `pr_ndlss_sc_v0.strc-species` |
| 60 | Nd | semicore | `nd_ndlss_sc_v0.strc-species` |
| 61 | Pm | semicore | `pm_ndlss_sc_v0.strc-species` |
| 62 | Sm | semicore | `sm_ndlss_sc_v0.strc-species` |
| 63 | Eu | semicore | `eu_ndlss_sc_v0.strc-species` |
| 64 | Gd | semicore | `gd_ndlss_sc_v0.strc-species` |
| 65 | Tb | semicore | `tb_ndlss_sc_v0.strc-species` |
| 66 | Dy | semicore | `dy_ndlss_sc_v0.strc-species` |
| 67 | Ho | semicore | `ho_ndlss_sc_v0.strc-species` |
| 68 | Er | semicore | `er_ndlss_sc_v0.strc-species` |
| 69 | Tm | semicore | `tm_ndlss_sc_v0.strc-species` |
| 70 | Yb | semicore | `yb_ndlss_sc_v0.strc-species` |
| 71 | Lu | semicore | `lu_ndlss_sc_v0.strc-species` |
| 72 | Hf | semicore | `hf_ndlss_sc_v0.strc-species` |
| 73 | Ta | semicore | `ta_ndlss_sc_v0.strc-species` |
| 74 | W | semicore | `w_ndlss_sc_v0.strc-species` |
| 75 | Re | standard | `re_ndlss_v0.strc-species` |
| 76 | Os | standard | `os_ndlss_v0.strc-species` |
| 77 | Ir | standard | `ir_ndlss_v0.strc-species` |
| 78 | Pt | standard | `pt_ndlss_v0.strc-species` |
| 79 | Au | standard | `au_ndlss_v0.strc-species` |
| 80 | Hg | standard | `hg_ndlss_v0.strc-species` |
| 81 | Tl | semicore | `tl_ndlss_sc_v0.strc-species` |
| 82 | Pb | semicore | `pb_ndlss_sc_v0.strc-species` |
| 83 | Bi | semicore | `bi_ndlss_sc_v0.strc-species` |
| 84 | Po | semicore | `po_ndlss_sc_v0.strc-species` |
| 85 | At | standard† | `at_ndlss_v0.strc-species` |
| 86 | Rn | only | `rn_ndlss_v0.strc-species` |
| 87 | Fr | semicore | `fr_ndlss_sc_v0.strc-species` |
| 88 | Ra | semicore | `ra_ndlss_sc_v0.strc-species` |
| 89 | Ac | standard | `ac_ndlss_v0.strc-species` |
| 90 | Th | standard | `th_ndlss_v0.strc-species` |
| 91 | Pa | standard | `pa_ndlss_v0.strc-species` |
| 92 | U | standard | `u_ndlss_v0.strc-species` |
| 93 | Np | standard | `np_ndlss_v0.strc-species` |
| 94 | Pu | standard | `pu_ndlss_v0.strc-species` |
| 95 | Am | semicore | `am_ndlss_sc_v0.strc-species` |
| 96 | Cm | semicore | `cm_ndlss_sc_v0.strc-species` |
| 97 | Bk | semicore† | `bk_ndlss_sc_v0.strc-species` |
| 98 | Cf | semicore | `cf_ndlss_sc_v0.strc-species` |
| 99 | Es | semicore† | `es_ndlss_sc_v0.strc-species` |
| 100 | Fm | semicore† | `fm_ndlss_sc_v0.strc-species` |
| 101 | Md | semicore† | `md_ndlss_sc_v0.strc-species` |
| 102 | No | semicore† | `no_ndlss_sc_v0.strc-species` |
| 103 | Lr | semicore† | `lr_ndlss_sc_v0.strc-species` |
| 104 | Rf | semicore† | `rf_ndlss_sc_v0.strc-species` |
| 105 | Ha | semicore† | `ha_ndlss_sc_v0.strc-species` |
