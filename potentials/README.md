# potentials

A curated library of validated, reusable PAW setup ("potential") definitions
for VESPAW.

## Why this exists

Upstream CP-PAW has no separate pseudopotential/setup library. Each `!SPECIES`
block (`!NTBO` + `!AUGMENT` sub-blocks) is inlined directly into every
`.strc` file and regenerated from scratch (`SETUP_NEWPRO`, `SETUP_TEST_NEWPRO`,
the "GHOSTBUSTER" checks, etc.) the moment the calculation starts. That works,
but it means:

- there is no single place to find "the validated Ti setup" — it has to be
  copy-pasted out of whichever `.strc` file used it last (see
  `tio2_test/tio2.strc`, `nio_pbe0r/*.strc`, etc. in the scan0r working
  directory for examples currently only living inline like this);
- two `.strc` files can silently drift to slightly different parameters for
  the "same" element with no diff-able record of why;
- there's no VASP-`POTCAR`-style artifact to hand to someone else or to a
  future ML-based initial-guess model as a stable input feature.

## Convention

One file per element/augmentation-type combination, named
`<element>_<augtype>_<tag>.strc-species` (a bare `!SPECIES...!END` block,
directly copy-pasteable into a `.strc` file's `!STRUCTURE` block), with a
header comment recording: source, which standard-vs-semicore variant was
picked and why, and known limitations.

## Status (2026-09-04): populated, full periodic table, mostly NOT production-validated

105 elements (H through element 104/105, `0_` "ghost/vacuum" species
excluded), one file each. Element/variant *structure* (`NPRO`, `ID`, `ZV`) —
**extracted directly from Peter Bloechl's own CP-PAW `Element0ktest`
benchmark species list**
(`src/Benchmarks/Element0ktest/src/specieslist` of CP-PAW,
the 0K elemental-properties benchmark validated against Lejaeghere et al.,
Crit. Rev. Solid State Mater. Sci. 39, 1 (2014)). For elements with both a
standard and a semicore (`_SC`) CP-PAW variant (89 of the 105), the one
roughly matching VASP's typical recommended PAW-PBE POTCAR default (plain vs.
`_sv`/`_pv`/`_d`) was picked — see `VASP_VARIANT_MAPPING.md` for the full
per-element table and reasoning (**best-effort from general VASP-convention
knowledge, not an independently verified lookup** — see the **known
mismatch on Li** below, keep that in mind for any other element too).

**Two tiers of data quality here, not one:**

- **H, Li, O, Mn** (`h_ndlss_v0`, `li_ndlss_v0`, `o_ndlss_v0`,
  `mn_ndlss_sc_v0`): **real, published, production parameters**, transcribed
  directly from Peter Bloechl's own setup blocks as printed in the PBE0r
  paper (Bloechl group, Phys. Rev. B 95, 235150 (2017)) — `RCL/RCOV=0.7`,
  `CORE POW=2`, `LRHOX=4`, `RAD/RCOV=1.2`, plus the paper's own `!NTBO`
  onsite-hybrid block (`LHFWEIGHT`, `TAILLAMBDA`, `RAUG/RCOV`, `RTAIL/RCOV`)
  for PBE0r/SCAN0r use. Live-tested (single-atom `SETUP_NEWPRO`/
  `GHOSTBUSTER`, `TYPE='PBE'`) — clean, no ghost-state warning. **This also
  caught a real error**: the variant mapping above had picked *semicore* Li
  (`Li_sv`-style) to match VASP's convention, but Peter's actual published Li
  is the *standard* (non-semicore) variant — `VASP_VARIANT_MAPPING.md`'s
  heuristic is not reliable evidence on its own, treat every entry there as
  unverified until checked some other way, as was possible here for Li.
- **The other 101 elements**: still generic, but the *ratios* now follow
  Peter's own PRB2017 general-purpose "house style"
  (`RCL/RCOV=0.7`/`CORE POW=2`/`LRHOX=4`/`RAD/RCOV=1.2`, matching what's
  consistent across all four validated elements above) rather than
  `Element0ktest`'s slightly different defaults (`RCL/RCOV=0.8`,
  `CORE POW=3`, `LRHOX=2`, no `RAD/RCOV`) — a better-grounded generic
  default, live-tested clean on Fe as a spot check, but **still not
  individually validated for that specific element**, and none of them carry
  the `!NTBO` hybrid block (only known for H/Li/O/Mn). Treat these as a
  documented starting point for building a `.strc`, not as validated
  production setups — run `GHOSTBUSTER`/`EPWPSI` convergence/equilibrium-
  volume-vs-experiment before trusting one in production. Elements already
  proven in real VESPAW runs (or found in a future paper/validated source)
  should get the same paper-exact treatment as H/Li/O/Mn, superseding the
  generic entry for that element.
