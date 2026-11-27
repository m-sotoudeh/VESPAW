# Plane-wave Fock exchange in VESPAW: design note

Started 2026-11-23. Aligned with the VASP 6.4.2 source on 2026-11-23 (Sec. 8;
read in third-party/VASP-6.4.2/vasp.6.4.2/src, not copied). Goal: VASP-consistent PBE0, HSE06, SCAN0 and B3LYP with a
nonlocal plane-wave PAW Fock operator, coupled to the variational SCF, MPI-correct,
and laid out so that the GPU port (stage 2: GPU-resident H|psi>) can reuse it.

## 1. Functionals (VASP conventions)

E_xc = a E_x^HF[v] + (1-a) E_x^SL[v] + E_c^SL, where v is the Coulomb kernel.

| Functional | a | v | Semilocal part | VASP tags |
|---|---|---|---|---|
| PBE0  | 0.25 | 1/r | 0.75 PBE x + PBE c | LHFCALC, AEXX=0.25 |
| HSE06 | 0.25 | erfc(mu r)/r, mu=0.2 A^-1 | PBE x - 0.25 wPBE_SR x + PBE c | HFSCREEN=0.2 |
| SCAN0 | 0.25 | 1/r | 0.75 SCAN x + SCAN c | METAGGA=SCAN, AEXX=0.25 |
| B3LYP | 0.20 | 1/r | 0.08 LDA x + 0.72 B88 x + 0.19 VWN3 c + 0.81 LYP c | GGA=B3, AGGAX=0.72, AGGAC=0.81, ALDAC=0.19 |

XC work needed: (a) scale exchange only (DFT_MODULE%SCALEX already exists for GGA,
needs a SCAN-exchange scale in the SCAN routines, one-center and grid); (b) separate
LDA and B88 gradient scales; VWN3 (check, the audit note says it was added); (c) the
short-range wPBE exchange of Heyd-Scuseria-Ernzerhof (HSE03/06 papers, EP98 hole),
written from the papers and validated against libxc GGA_X_WPBEH (mu=0.106 bohr^-1
in libxc units, check VASP 0.2 A^-1 = 0.10583545 bohr^-1). All four checked
pointwise against libxc as in the 2026-10-02 XC audit.

## 2. The operator

For spin s, band n at k (all quantities pseudo (PS) unless marked):

  (V_x psi_nk)(r) = - sum_{m,q} f_{m,k-q} w_{k-q} phi_{m,k-q}(r) int dr' v(r-r') phi*_{m,k-q}(r') psi_nk(r')

E_x^HF = -1/2 sum_{nk,mq} f_nk f_mq w_k w_q (nk,mq|v|mq,nk), spin diagonal (collinear),
both spin channels for NDIM=2 (noncollinear: later, not in milestone 1-4).

Plane-wave part: pair density rho_mn(r) = u*_{m,k-q}(r) u_{nk}(r) on the wave-function
grid of the k-group, FFT, multiply by v(G+q), FFT back, multiply by u_{m,k-q}.
Cost per H|psi> call ~ NB_occ * NB * NK * 2 FFT, which is why we use ACE (Sec. 5).

Kernels v(G+q):
- 1/r: 4 pi/|G+q|^2 and erfc(mu r)/r: 4 pi/|G+q|^2 (1 - exp(-|G+q|^2/4mu^2)).
  The G+q=0 element is the VASP FSG value in both cases (Sec. 8.3), so the
  HSE kernel also needs it for VASP agreement.

## 3. PAW

Same structure as VASP (Paier et al. JCP 122, 234102 (2005)):
1. Plane-wave part with the PS pair density plus compensation charges
   hat-rho_mn = sum_{ij,L} <psi_m|p_i> <p_j|psi_n> Q_ij^L g_L(r-R) so that the
   multipoles of the pair density are correct in the interstitial. VESPAW already
   builds compensation densities for the Hartree term (WAVES$MOMENTS, AUGMENTATION);
   reuse the shape g_L and Q_ij^L from the setups.
2. One-center correction per atom: sum_ijkl D_ij D_kl [ (phi_i phi_k|v|phi_l phi_j)
   - (tphi_i tphi_k + hat|v|tphi_l tphi_j + hat) ], with the exchange (not Hartree)
   index pattern and spin-resolved density matrices. Four-index radial Slater
   integrals per species, computed once (radial Poisson in atomlib, reuse
   RADIALFOCK and the onsite-hybrid U-tensor code of vespaw_lmto/ldaplusu where it
   fits). For HSE the one-center kernel is erfc-screened (VASP also screens it).
3. Core-valence exchange (frozen core): from the atom, adds a fixed one-center
   operator sum_c (phi_i phi_c|v|phi_c phi_j) to DH (scaled by a). Needs core
   orbitals from the setup construction (atomlib has them).

All one-center pieces enter through DH (the existing |p>DH<p| path), so the MD,
CG, DAV, RMM and EDFT paths pick them up without extra hooks. Forces then follow
from the existing dE/d<p|psi> machinery, except the compensation term in the
plane-wave part (new force and stress terms, milestone 6).

## 4. Data layout and MPI

Current layout: k-points distributed over k-groups (KMAP), real-space grid split
over tasks within a group (NR1L planes), bands packed as NBH=(NB+1)/2 at TINV
k-points (two real orbitals per complex vector), GAMMASUPER half-sphere storage
(not trusted yet, see memory).

Fock needs the occupied orbitals at all k-q. Plan:
- After each outer step: every k-group broadcasts its occupied orbitals
  (G-space, NGL x NOCC per k) to all k-groups (MPE in the 'MONOMER' world),
  stored in a FOCK_MODULE buffer. Memory: NK * NOCC * NGL * 16 B; for 64-atom
  oxides at 4 k-points this is ~1 GB, so later band-distributed storage.
- Pair densities live on the k-group's own grid slice (same FFT as WAVES_VPSI).
  Orbitals at k-q, k-q outside the IBZ: milestone 1-4 run without SYMMETRY=T
  (full mesh, TINV allowed). Symmetry (rotate orbitals from IBZ to the star,
  including projections and phases) is milestone 7.
- The GPU layout will be the same: batched pair-density FFTs over bands.

## 5. ACE and coupling to the variational SCF

Adaptively compressed exchange (Lin, JCTC 12, 2242 (2016)): after each outer step
build W_n = V_x[phi] phi_n for the occupied+computed bands, then
V_ACE = -W (M)^-1 W^dagger with M = phi^dagger W. Applying V_ACE in H|psi> costs
two ZGEMMs instead of NOCC FFT pairs. The diagonalization (DAV/RMM) then uses a
fixed Fock operator during the step, like VASP's inner loop.

Coupling to MIXER='VARIATIONAL' follows the onsite-hybrid precedent
(WAVES_VARMIXHAMILTONIAN comment): the Fock operator is taken from the current
state, not line-searched; the local potential and DH (now including the one-center
Fock pieces) are mixed as before. Every accepted energy is that of an actual
state (E_x^HF evaluated with the orbitals that built the density), so the
descent property and the C1/C2 certification carry over. If convergence with
the current-state operator stalls (as tried for the NTBO block), the fallback
is a two-loop scheme: inner variational SCF at fixed ACE, outer ACE update.

The MD (Car-Parrinello friction) path: V_x applied directly (not ACE) to the
current bands in WAVES$HPSI, energy and forces consistent by construction.

## 6. Input

!CONTROL!DFT TYPE='PBE0' | 'HSE06' | 'SCAN0' | 'B3LYP' (names map to the
semilocal TYPE plus a FOCK block), and an optional

  !CONTROL!DFT!FOCK AEXX=0.25 MU[1/ANGSTROM]=0.2 SING='VASP'|'SPENCER'|'NONE'
                    ACE=T NBANDSFOCK=... !END

Pair densities on the wave-function grid (equivalent to VASP PRECFOCK=Accurate,
Sec. 8.4). SING='VASP' = FSG of Sec. 8.3.

## 7. Milestones and validation

1. Plane-wave PS Fock at Gamma, TINV packing, NSPIN 1/2, energy + operator in
   WAVES_HPSI and WAVES$HPSI, kernels 1/r (no G=0 term) and erfc.
   Tests: (a) finite difference dE_x/dpsi vs operator (directional derivative,
   1e-8 rel); (b) H2 / He in a box: PS exchange energy vs -1/2 J for a single
   doubly occupied orbital (exact identity E_x = -E_H/2 for one orbital pair);
   (c) Hermiticity <a|V_x|b> = <b|V_x|a>*.
2. PAW one-center exchange + compensation charges in the pair density.
   Test: isolated atom (Ne, Mg, O in a large box, Gamma) HF exchange vs the
   radial atomic HF code (RADIALFOCK) with the same frozen orbitals.
3. XC pieces: exchange scaling, wPBE_SR, B3LYP parts, SCAN exchange scaling;
   pointwise libxc comparison.
4. k-points (full mesh), singularity correction, ACE, variational SCF coupling.
   Tests: Si, C, MgO PBE0/HSE06 gaps and total energies vs VASP 6.4.2
   (local install, same POTCARs is not possible: compare gaps and energy
   differences, e.g. lattice constant and Si-C formation energy).
5. MPI: k-group broadcast, band distribution; bitwise serial/MPI agreement.
6. Forces and stress (compensation term, singularity term), L-BFGS relaxation
   test: Si bond, MgO a0 vs VASP.
7. Symmetry (IBZ -> full BZ orbital rotation), noncollinear, GAMMASUPER.
8. Validation set: Si, C, MgO, NiO (AFM-II), MnO, CaMnO3, TiO2: PBE0/HSE06/SCAN0/B3LYP
   gaps, moments, a0 vs VASP; timing vs VASP and vs the onsite hybrids.

Ground rules: written from the papers; VASP source may be read, never copied
(same rule as r2SCAN). Numbers quoted in papers come only from the current code.

## 8. Alignment with VASP 6.4.2 (source read 2026-11-23)

Files read: fock_glb.F, fock.F (XC_FOCK_READER, SETUP_FOCK, SET_GFAC_VKPT,
SET_FOCK_ALPHA, SET_FSG, SET_FSG_STORE, FOCK_ACC, RSPARAMETERS), fock_ace.F,
pawfock.F (COLOUMB_4TERM(_PS), RAD_POT_EX_HAR, CALC_DHARTREE, SETUP_PAWFOCK,
CALC_PAWFOCK), paw.F (SET_PAW_ATOM_POT_RHF, SET_PAW_CORE_FOCK), setex.F,
xclib_grad.F (wpbe module), metagga.F, electron.F (ACE call site).

### 8.1 Defaults VASP uses (and VESPAW will copy as defaults)
- LHFCALC=.TRUE. -> AEXX=0.25, ALDAX=AGGAX=1-AEXX, ALDAC=AGGAC=1 (AEXX=1 also
  drops correlation). Meta-GGA: AMGGAX=1-AEXX (SCAN exchange scaled, SCAN
  correlation full). B3LYP: AEXX=0.2, AGGAX=0.72, AGGAC=0.81, ALDAC=0.19 with
  VWN3 local correlation (LEXCH=11, "B3"); ALDAX is 1-AEXX=0.8 (0.08 = 0.8-0.72
  is the pure LDA exchange share).
- HSE06 = PBE + HFSCREEN=0.2 A^-1 (0.105835 bohr^-1). libxc HYB_GGA_XC_HSE06
  uses omega=0.11 bohr^-1: compare XC pointwise against libxc with omega set to
  0.105835, not the libxc default.
- ISYM -> 3 (HF symmetry: full k mesh, orbitals rotated from the IBZ),
  ISIF -> 0, IMIX=1 (Kerker, no Broyden), FOCKCORR=1, HFALPHA automatic,
  LMAX_FOCK=4, LMAX_FOCKAE=-1 (no AE restoration on the grid), NMAX_FOCKAE=1,
  PRECFOCK=normal, ENCUTFOCK off, ACE on (LFOCKACE, -Dfock_dblbuf).

### 8.2 Operator and energy conventions (checked against FOCK_ACC)
- Occupation weight per q-band is FERWE in [0,1] per spin channel (VASP stores
  per-spin occupations also for ISPIN=1 and multiplies energies by RSPIN=2).
  In VESPAW OCC includes the k weight and the spin degeneracy (f<=2 for
  NSPIN=1), so V_x psi = -s sum_m f_m phi_m (v*(phi_m^* psi)) with s=1/2
  (NSPIN=1) or 1 (NSPIN=2), times AEXX, k-weights from the full mesh.
- Energy: E_x = 1/2 sum_n f_n <psi_n|V_x|psi_n> (VASP CDCHF "double counting"
  form). Same in VESPAW, evaluated from the snapshot orbitals.
- The kernel acts on every G of the Fock FFT box (no sphere), except
  ENCUTFOCK. 1/r: 4 pi/|G+q|^2. HSE: 4 pi/|G+q|^2 (1-exp(-|G+q|^2/4 mu^2)).
  The G+q=0 element is set to FSG for every pair, and the whole kernel is
  multiplied by AEXX.

### 8.3 Singularity (FSG) -- corrects Sec. 2
- 1/r with a full Monkhorst-Pack or Gamma-centered mesh: FSG = 2 x the
  Madelung (Ewald monopole) energy of a unit point charge in the Born-von
  Karman supercell (cell vectors a_i * NKP_i) with neutralizing background
  (EWALD_MONO). Not Gygi-Baldereschi. VESPAW: reuse its Ewald code on the
  supercell lattice.
- HSE (and other k sets): auxiliary-function form, FSG = self-energy of a
  Gaussian (HFALPHA) charge minus its lattice sum over G+q != 0 of
  e^{-alpha |G+q|^2} v(G+q), with the screened self-energy subtracted for HSE.
  So HSE also gets a finite G+q=0 value (not pi/mu^2): needed to match VASP.
- HFALPHA default: 6/(ENMAX[Ry]) / (2 pi) in VASP length units (set in
  SETUP_FOCK, the older SET_FOCK_ALPHA is commented out).
- FSG is stored per k-point but is the same for all k on a uniform mesh.

### 8.4 FFT grid for the Fock part -- corrects Sec. 2/6
- PRECFOCK=normal with PREC=Normal uses the wave-function grid itself.
  PRECFOCK=accurate rebuilds a grid with 2 x Gmax(psi) span (no wrap-around),
  normal/low/fast a 1.5 x span with some aliasing.
- VESPAW's wave-function grid is set from EPWRHO=CDUAL*EPWPSI. With CDUAL>=4
  (the default) its sphere spans 2 Gmax(psi), so the natural VESPAW choice is
  equivalent to PRECFOCK=accurate. With CDUAL=2 (used in many test decks) the
  pair densities are cut at sqrt(2) Gmax(psi), consistently in energy and
  operator, which resembles the coarser VASP settings; a warning is printed. Implementation: a
  separate plane-wave set 'FOCK' on the same 'K' grid slice (same NR1START/NR1L,
  NR2, NR3), KVEC=q, TINV=F, cutoff EPWRHO. A product of two functions inside
  Gmax(psi) lies inside 2 Gmax(psi), so the 'FOCK' sphere is exact for the
  plane-wave part, and the compensation part is cut at EPWRHO as for the
  Hartree term. To compare with VASP use PRECFOCK=Accurate there.

### 8.5 PAW -- refines Sec. 3
- Plane-wave pair density includes compensation charges up to L=LMAX_FOCK=4,
  built from the soft augmentation shape (AUG_SOFT) and the PAW moments QPAW.
  Its back-projection gives D_ij contributions applied through the
  projectors (OVERL_FOCK). In VESPAW: use the same Gaussian compensation
  shapes as WAVES$MOMENTS/ the Hartree term, L<=4 cut, D_ij via |p>D<p|.
- One-center: Slater integrals S(ijkl,L) = S_AE(phi) - S_PS(tphi + aug),
  all L up to 8 (not cut at LMAX_FOCK), with exchange index pattern
  DFOCK_il = - sum_kj c_ij^LM S c_kl^LM rho_kj (CALC_DHARTREE). Spin: per
  spin channel. NSPIN=1: rho is the per-spin matrix, energy x2.
- The one-center radial kernel is screened for HSE (erfc radial Green
  function, BRSEXX=AEXX, BEXX=0) and scaled by AEXX for PBE0/B3LYP/SCAN0.
- Core-valence exchange: from the AE core orbitals of the atom (recomputed in
  the frozen AE potential), sum over core shells, added as a fixed D_ij
  (scaled and screened like the valence part). The DFT core-valence xc is
  removed by evaluating the semilocal xc with the scaled exchange on
  rho_v+rho_c (as now in VESPAW) and adding a*E_x^HF(core-valence).
  VESPAW: core orbitals come from the setup construction in atomlib; store
  them (or recompute in the frozen AE potential) at setup read.
- VASP also changes the atomic reference energy (EATOM) for the hybrid. VESPAW
  has no such reference (absolute frozen-core total energies), so energy
  DIFFERENCES are compared with VASP, never total energies.

### 8.6 Semilocal parts
- HSE: E_x^SL = E_x^PBE - AEXX E_x^{wPBE,SR}(mu), LDA part analogously with the
  short-range LDA exchange. VASP's wPBE: HSE03 (JCP 118, 8207) with the
  Ernzerhof-Perdew 1998 hole, Lieb-Oxford-enforced, tabulated in (rho,s) and
  spline-interpolated. VESPAW: analytic EP98 hole (as libxc GGA_X_WPBEH), check
  that the VASP table differs by less than 1e-5 Ha/atom on Si and MgO.
- SCAN0: (1-AEXX) SCAN exchange + full SCAN correlation, grid and one-center.

### 8.7 SCF coupling -- confirms Sec. 5
- electron.F rebuilds the ACE operator from the current orbitals in every
  electronic step (after NELMDL), i.e. a fixed operator during each
  diagonalization, refreshed each step. Same as the VESPAW plan (current-state
  operator, Sec. 5).

### 8.8 Changes to the milestones
- Milestone 1 tests stay (operator/energy consistency, Hermiticity).
- Milestone 4: singularity = Ewald supercell form for Coulomb on MP/Gamma
  meshes, auxiliary function for HSE; check FSG against VASP OUTCAR values
  for Si/MgO.
- Milestone 8: VASP runs with PRECFOCK=Accurate, ISYM=3, same k mesh;
  compare gaps, band energies relative to VBM, lattice constants and
  formation energies.

## 9. Status

### Milestone 1 done (2026-11-23)
Code: src/vespaw_fock.f90 (FOCK_MODULE), hooks in vespaw_waves1.f90
(WAVES$INITIALIZE, WAVES$ETOT -> WAVES_FOCKETOT, WAVES_HPSI, WAVES$HPSI,
WAVES_FOCKCHECK), input !CONTROL!DFT!FOCK AEXX= MU= | HFSCREEN= !END in
READIN_DFT, object added to Makefile_bare(.in) and the arm64 object Makefiles.
Fock is added on top of the semilocal functional (no exchange scaling yet).
Tests in regression_hard/fock_m1 (make.py, make_si8.py), EPWPSI 25-30 Ry, CDUAL=2:
- operator vs energy, VESPAW_FOCK_CHECK=1: E_x = 1/2 sum f <psi|V_x|psi> to 1e-15 H
  in every iteration (H2O, TINV-packed Gamma).
- Si8 Gamma (no vacuum): Fock DAV converges like PBE (16 vs 15 iterations, residual
  2e-9), E=-896.01099 eV; friction MD -896.01086 (autopilot stop, still descending).
- H2O NSPIN=2 identical to NSPIN=1 (-488.21486467 eV): spin factors correct.
- serial vs 2-rank MPI: identical to 1e-8 eV over 20 iterations.
- no-Fock regression: Si8 PBE identical to the pre-Fock binary in every iteration.
- H2O in an 8 A box: DAV/RMM stop 0.40 meV above friction MD; PBE alone shows the
  same pattern (0.18 meV), i.e. the known vacuum GGA residual floor, not Fock.
- HSE kernel (HFSCREEN=0.2) runs; G=0 term still omitted (milestone 4).
- Cost without ACE: Si8 Fock 45 s vs PBE 8 s, FOCK$APPLY 82 %.

### Milestone 2 (2026-11-24): PAW plane-wave + one-center, G=0 term
Code: src/vespaw_fock.f90 (FOCK_SPECIES, FOCK_ADDCOMP, FOCK_COMPPROJ,
FOCK$ONECENTER, FOCK_V0, FOCK_SHAPEFT, FOCK_SELFENERGY), calls in
vespaw_waves1.f90 (FOCK$APPLY before WAVES_ADDPRO with PROJ/HPROJ,
FOCK$ONECENTER after ONSITE$ETOT). Input RCOMP= (bohr).
Design as built:
- Compensation charges in every pair density, all L up to 2*lmax, shape
  r^L (1-(r/Rc)^2)^6 confined to Rc (default 0.6*RBOX, about the augmentation
  radius), Fourier transform tabulated numerically (4-point Lagrange).
- Back-projection of the pair potentials onto the compensation charges goes to
  HPROJ, added by WAVES_ADDPRO.
- One-center: Slater tensor U = (ij|kl)_AE - (ij|kl)_PS+comp from AEPHI/PSPHI,
  radial densities cut at RINT=2*RBOX (the AE-PS difference of s,p partial
  waves decays smoothly, core-orthogonality tails, 1e-3 at RBOX for Si).
  Enters DH like DFT+U; E=-AEXX/2 sum D_ik D_lj U(ijkl) per spin.
- On-site compensation self-energy: the grid holds (ghat|ghat) only for
  |G|<Gmax, the radial term is exact; F_PS is corrected by
  Q_ab Q_cd [(ghat|ghat)_trunc - (ghat|ghat)_exact], (ghat|ghat) =
  1/(2pi^2) int FG(q)^2 dq (checked against the radial Poisson value to 1e-8).
- G=0: Gamma-only BZ average of the kernel. Coulomb: exact Ewald form
  Omega/sqrt(pi a) - sum' 4pi e^{-aG^2}/G^2 + 4pi a - Omega sum_R erfc(R/2sqrt a)/R
  (= -2 Omega E_Madelung, VASP MP branch); alpha-independent to all digits,
  Si8 cubic 298.746 bohr^2 = Madelung 2.8373 Omega/L. erfc kernel: auxiliary form.
Failed attempts (kept for the record):
- Gaussian compensation shapes: leak out of the sphere, E_x drifted 1.7 mHa with
  the width (frozen Si8).
- Band-limited back-transform of the confined shape: Gibbs ripples destroy the
  L=3,4 moments (moment -8 for L=4).
- G=0 aux form without the 4pi*alpha term: 5.8 mHa CDUAL dependence (Si8).
Tests (regression_hard/fock_m1, frozen orbitals from a converged PBE Si8, CDUAL 4/8):
- operator/energy 1e-13, one-center E1 = 1/2 sum DH1 D exactly.
- CDUAL 4 vs 8: E_x equal to 2e-7 H (exact pair densities at CDUAL>=4).
- Rc 1.4 -> 2.0: split moves by 63 mH, total by 0.44 mH (NPRO 2 2 1) and
  3.1 mH (NPRO 1 1 1): the residual is the setup completeness error.

### Milestone 2b (2026-11-24): screened one-center kernel, core-valence, screened G=0
Code: FOCK_ERFPART, FOCK_ERFPARTCV, FOCK_BESSELTRANSFORM, FOCK_ERFQGRID in
src/vespaw_fock.f90, core block in FOCK_SPECIES, KCV term in FOCK$ONECENTER.
Input !FOCK CORE=T|F (default T). The HSE guard is removed.
- erfc one-center kernel: the Coulomb Slater integrals minus the long-range part,
  (rho1 Y_L|erf|rho2 Y_L) = 8 int dq exp(-q^2/4mu^2) rho1(q) rho2(q),
  rho(q) = int rho(r) j_L(qr) r^2 dr (Simpson, q up to 2 mu sqrt(37)). This is the
  same interaction as the plane-wave kernel 4pi/q^2 (1-exp(-q^2/4mu^2)). VASP
  builds it in real space (RS_GreensFunc).
- Compensation self-energy correction with screening: truncated and exact values
  both use the screened kernel (exact = radial Coulomb - erf part over all q).
- Core-valence: KCV(i,j) = sum_{c,mc} (phi_i phi_c|phi_c phi_j), AE only, from the
  setup core states (AEPSI 1..NC, FB=2(2l+1)); E_cv = -AEXX tr(D_tot KCV),
  spin independent, enters DH(:,:,1). Same kernel (Coulomb or erfc) as the rest.
  AEPSI is the scalar-relativistic large component (norm 0.9975 for Si 1s, the
  missing (Z alpha)^2/4), renormalized to one. Core-core exchange is a constant and
  is omitted (as in VASP).
- Screened G=0 term: the auxiliary form is not alpha independent for the erfc
  kernel (Si8: 167.34, 168.55, 169.16 for alpha=.25,.125,.0625) and gave a 4.5 mH
  CDUAL dependence. Replaced by the same definition as the Coulomb Ewald term,
  V0 = pi/mu^2 - Omega sum_{R/=0} erfc(mu R)/R (Si8: 169.765), which is the alpha->0
  limit of the auxiliary form and tends to the Coulomb value for mu->0
  (298.75 - 2 Omega mu/sqrt(pi), checked numerically).
Tests (frozen PBE Si8 orbitals, regression_hard/fock_m1/run_m2b.sh), totals in H:
- Coulomb, CORE=F: CDUAL 4 -2.24530709, CDUAL 8 -2.24530708.
- HSE (HFSCREEN=0.2), CORE=F, CDUAL 4 / 8:
  Rc 1.4: -1.76726181 / -1.76726164, Rc 1.7: -1.76712309 / -1.76712309,
  Rc 2.0: -1.76694165 / -1.76694165. The split moves 63 mH, the total 0.32 mH
  (setup completeness, as for Coulomb).
- Core-valence: E_cv = -0.53438 (Coulomb), -0.53427 (HSE) for 8 atoms.
- operator/energy 1e-13, one-center E1 = 1/2 sum DH1 D to 1e-16.
SCF (regression_hard/fock_m1/run_m2c.sh, Si8, HFSCREEN=0.2, CORE=T, CDUAL 4):
- DAV converges in 15 iterations (as Coulomb), E=-921.45842362 eV, residual 1.7e-9,
  E_cv=-0.57789 H at self-consistency; operator/energy 2e-14 in every iteration.
- serial vs 2-rank MPI: energies and eigenvalues identical over 8 iterations.
- FOCK$APPLY is still 98 % of the run time (580 s), ACE is milestone 4.
Missing: semilocal exchange scaling (milestone 3).

### Milestone 3a (2026-11-24): PBE0 exchange scaling
Code: DFT_MODULE SCALEXHYB (vespaw_dft.f90; factor XFAC=SCALEX*SCALEXHYB on the
exchange in DFT, DFT2, DFT3, i.e. grid and one-center AE/PS), DFT$SETR8/GETR8,
DFT$SETL4('HYBRIDOFF'), READIN_FOCKXSCALE (vespaw_ioroutines.f90), SETUP$READ.
- E_xc = E_xc^GGA - AEXX E_x^GGA + AEXX E_x^HF (VASP ALDAX=AGGAX=1-AEXX), for the
  LDA/GGA types 1, 2, 10-14. SCAN0 and B3LYP still stop with an error.
- The atomic setups are built with HYBRIDOFF=T (no scaling, no short-range term),
  as VASP uses its PBE datasets for hybrids. psphi.dat identical to PBE (checked).
- Fock and the onsite hybrid (!NTBO) cannot be combined (error).
- VESPAW_FOCK_XSCALE=0 restores the milestone 1/2 test mode (Fock on top of the
  full semilocal exchange).
Tests (regression_hard/fock_m3/run_m3a.sh, frozen PBE Si8 orbitals):
- AEXX=0 identical to PBE (-857.70157758 eV).
- E(AEXX=0, .25, .5) exactly linear (second difference 1e-13 eV).
- PBE0 DAV SCF: 16 iterations, -858.46980451 eV, residual 7e-10.

### Milestone 3b (2026-11-24): HSE06 short-range wPBE
Code: WPBE_MODULE, WPBE_DIRECT, WPBE_INITIALIZE, WPBE_INTERPOLATE, WPBE_UNPOL,
WPBE$SR, WPBE_SELFTEST at the end of vespaw_dft.f90; HYBSRA/HYBMU in DFT_MODULE.
- E_xc^HSE = E_xc^PBE - AEXX E_x^wPBE,SR(mu) + AEXX E_x^HF,SR(mu), mu from HFSCREEN.
- Short-range wPBE from the Ernzerhof-Perdew hole with the HSE parameters and the
  hole sum rule (closed form with e^x erfc(sqrt x)), integrated with the exact erfc
  by a double-exponential rule (1e-13), tabulated on s in [0,8.58] (ds=0.01) and
  t=nu/(1+nu) in [0,1] (dt=0.002), Hermite bicubic interpolation (2e-10), energy
  and potential from the same interpolant. Large s as VASP (s>8.3 saturated).
- References: VASP 6.4.2 EXCHWPBE_R compiled as a standalone oracle (scratch only,
  not part of VESPAW): agreement 3e-8 relative in eps_x^SR (VASP uses a fitted
  erfc). libxc 7 GGA_X_WPBEH has a 1.3e-4 jump at s=1 (a branch absent in VASP and
  in the hole integral), so it is not used as the reference.
- The table has first derivatives only: DFT2/DFT3 stop when it is active, and the
  one-center XC uses the angular quadrature (ANGULARXC='EXACT', set automatically,
  also for spherical-only atoms).
Tests (regression_hard/fock_m3/run_m3b.sh, frozen PBE Si8 orbitals):
- HSE at AEXX=0 identical to PBE with ANGULARXC='EXACT' (-857.70168472 eV).
- E(AEXX=0, .25, .5) linear (to the 1e-8 eV print precision).
- VESPAW_WPBE_TEST=1: derivatives vs finite differences 1.6e-9.
- HSE06 DAV SCF: 15 iterations, -857.82942638 eV.
- Gamma-only Si8 gaps: PBE 0.605, HSE06 1.791, PBE0 2.586 eV.
- no-Fock regression identical.
VASP 6.4.2 comparison (regression_hard/fock_m3/vasp_si8/run_vasp.sh: same cell,
Gamma only, ENCUT 340 eV = 25 Ry, PREC=Accurate, PRECFOCK=Accurate, ISYM=0, PBE.54
Si POTCAR vs the VESPAW NDLSS setup NPRO 2 2 1), gaps at Gamma in eV:
            VESPAW   VASP    diff
  PBE       0.6047   0.6053  -0.6 meV
  HSE06     1.7912   1.7953  -4.1 meV
  PBE0      2.5857   2.5949  -9.2 meV
  hybrid shift of the gap: HSE06 1.1865 vs 1.1900 (0.3 %), PBE0 1.9810 vs 1.9896
  (0.4 %); HOMO shift PBE0 -1.167 vs -1.173, LUMO +0.814 vs +0.817.

### Milestone 3c (2026-11-24): SCAN0 / r2SCAN0
Code: SCAN_X_MODULE XSCALE and SCAN$SETXSCALE (vespaw_augmentation.f90),
READIN_FOCKXSCALE (types 21/22, requires TSCANDELOC=T, Coulomb kernel only).
- E_xc = E_xc^SCAN - AEXX E_x^SCAN + AEXX E_x^HF (VASP AMGGAX=1-AEXX). SCAN_X_FX is
  the single entry of the (r2)SCAN exchange (grid, one-center, nonspherical, all
  potentials by differences of it), so the factor 1-AEXX is applied to F_x there.
- The global SCAN uses complex wave functions, so it cannot restart from a TINV
  (PBE) restart file; the frozen tests start from converged SCAN orbitals.
Tests (regression_hard/fock_m3/run_m3c.sh, run_m3c_frozen.sh, Si8 Gamma):
- frozen SCAN orbitals: AEXX=0 identical to SCAN (-854.25071258 eV), E(AEXX=0,
  .25, .5) linear to the 1e-8 eV print precision, psphi.dat unchanged.
- SCF DAV: SCAN 21 iterations, SCAN0 24 iterations (-855.85455584 eV).
- VASP 6.4.2 (METAGGA=SCAN, LASPH=.TRUE., same settings as above), gaps in eV:
            VESPAW   VASP    diff
  SCAN      0.7985   0.7923  +6.2 meV
  SCAN0     2.6486   2.6531  -4.5 meV
  hybrid shift 1.8501 vs 1.8608 (0.6 %); HOMO -1.062 vs -1.067, LUMO +0.789 vs
  +0.794.
B3LYP: skipped at the user's request (2026-11-24).
Next: milestone 4 (k-points, ACE).

### Milestone 4a (2026-11-24): ACE (adaptively compressed exchange)
Code: FOCK$ACEBUILD, FOCK$APPLY (ACE), FOCK_APPLYEXACT (the former FOCK$APPLY),
FOCK_UNPACK in vespaw_fock.f90; build call in WAVES_FOCKETOT after FOCK$ENERGY;
input !FOCK ACE=T|F (default T).
- Once per energy evaluation, for every set: W = V_x psi for ALL NB bands (one exact
  application), M = <psi|W> (plane-wave part summed over 'K', projector part
  replicated), -M = L L^H (ZPOTRF), xi = W L^{-H} (ZTRSM) for the plane-wave and the
  projector part. Then V_x psi = -sum_k xi_k <xi_k|psi> (two ZGEMMs + one combine),
  HPROJ -= xi_P <xi|psi>. Exact on the span of the bands (as VASP), so residuals
  and the operator/energy identity at the current orbitals are exact; the energy is
  still the exact snapshot energy.
- TINV-packed pairs are split in real space (FFT), not with the -G map, which could
  cross 'K' tasks; at Gamma the xi are real functions, so the packed application
  works as for the exact operator.
- The snapshot reset invalidates the ACE vectors (fallback to the exact operator
  until the next build).
Tests (regression_hard/fock_m4/run_m4a.sh, Si8 PBE0, CDUAL 4, EPWPSI 25):
- DAV: -858.46980451 eV in 16 iterations, identical to the exact operator
  (fock_m3/pbe0_dav), wall time 32.5 s vs 616 s (19x).
- operator/energy check 3e-14 in every iteration (ACE application).
- RMM -858.46980450 (16 it.), friction MD -858.46972025 (autopilot stop).
- serial vs 2-rank MPI identical over 8 iterations.
- The build (one exact application per iteration) is now 57 % of the run time.

### Milestone 4b (2026-11-25): k-points
Code: FOCK$SETKMESH, FOCK_SELECTQ, FOCK_FINDK, FOCKQ_TYPE (per-Q pair-density
sets), FOCK_APPLYEXACT, FOCK$SNAPSHOT, FOCK$WBUILD (replaces FOCK$ACEBUILD) in
vespaw_fock.f90; WAVES_FOCKCHECK / WAVES_FOCKETOT in vespaw_waves1.f90.
- Full mesh = stored k plus time-reversal partners u_{-k} = conj(u_k), each with
  half the weight. FOCK$SETKMESH stops unless the full mesh has prod(NKDIV) points
  (no symmetry reduction yet).
- One pair-density plane-wave set per distinct Q = k - q mod 1 ('FOCK' for Q=0,
  'FOCKQnnnnn' otherwise). Kernels, Ylm, compensation shape transforms and
  structure factors are evaluated at K = Q+G. The folding vector G0 enters as the
  real-space phase exp(iG0 r) (PLANEWAVE_EIKR).
- G=0 term: V0 for the Born-von Karman supercell RBAS x NKDIV (Coulomb: Ewald form,
  screened: image sum). K=0 in a Q /= 0 set is an error.
- The energy is now taken from the build, E_x = 1/2 sum_n f_n Re M(n,n), for every
  mesh. FOCK$ENERGY (pair formula) is kept as a Gamma-only diagnostic and agrees
  with it to 3e-14.
- Restrictions: NKPTL = NKPT (no k-group distribution), SYMMETRY off, NDIM = 1.
Tests (regression_hard/fock_m4/run_m4b*.sh, Si8 PBE0, EPWPSI 22 Ry so that the FFT
grids 32/64 are commensurate, CDUAL 4):
- Gamma regression unchanged (-858.46980451 eV).
- Si8 1x1x2 vs Si16 Gamma (per Si8): PBE 1.4e-7 eV, PBE0 9.8e-7 eV, V0 identical
  (380.28315944). Checks the supercell V0, the Q sets, G0 folding and the k = X
  unpacking.
- Si8 1x1x4 vs Si16 1x1x2 (per Si8): PBE 2.5e-7 eV, PBE0 -8.0e-7 eV, V0 identical
  (-121.51818202). Checks the time-reversal partners (k = +-1/4).
- operator/energy check 3e-15 in every iteration of the k-mesh runs.
- Cost: Si8 PBE0 with ACE, 1x1x2 73 s, 1x1x3 146 s.
- Side finding (not Fock): the variational mixer stalls on the elongated Si16 and
  Si24 cells even for plain PBE (Si16 Gamma: 254 of 297 steps rejected). Pulay
  converges Si16 but charge-sloshes on Si24 Gamma (gap 0.003 eV), which is why the
  1x1x3 case is not compared against a supercell.

### Milestone 5 (2026-11-25): k-points over MPI k-groups
Code: WAVES_KDISTRIBUTEFOCK and the 'KX' communicator in WAVES$INITIALIZE,
FOCK$SNAPSHOTGATHER, FOCK$WBUILDDONE (vespaw_fock.f90), MPE$NEW list fix
(vespaw_mpelib.f90).
- With Fock on, all k-groups have the same number of tasks: the number of groups
  is the largest divisor of NTASKS not larger than NKPT, and the k-points go
  greedily to the least loaded group (a TINV point counts half). Equal groups have
  the same real-space slices, so the tasks of equal rank (communicator 'KX', one
  task per group) hold the same grid points.
- The snapshot is first filled with the local k-points, then FOCK$SNAPSHOTGATHER
  broadcasts every group's entries over 'KX' (orbitals on the local slice,
  projections, weights, mesh index). Each task then holds all occupied orbitals
  of the full mesh on its slice, as before with one group.
- The energy from FOCK$WBUILD (already summed inside a group over 'K') is summed
  over 'KX' in FOCK$WBUILDDONE. The one-center part needs nothing new (the
  density matrix is complete, as for DFT+U).
- MPE$NEW appended a new communicator directly after its parent, which dropped a
  communicator split earlier from the same parent ('K' was lost when 'KX' was
  split from 'MONOMER'). It now appends at the end of the list.
- Still on every task: the pair-density sets of all Q (also those only other
  groups need) and the full snapshot. Fine for the sizes so far.
Tests (regression_hard/fock_m4/run_m5.sh, run_m5b.sh, Si8 1x1x4 PBE0, 3 stored
k-points, 10 iterations, OMP_NUM_THREADS=1):
- serial vs 2 tasks (2 groups of 1), 3 tasks (3 groups of 1) and 4 tasks
  (2 groups of 2): energies identical in every iteration (to 1e-8 eV), operator/
  energy check 5e-15.
- spin-polarized (NSPIN=2) PBE0 with 4 tasks identical to serial.
- plain PBE with 4 tasks identical to serial (MPE$NEW fix), Gamma-only PBE0 with
  2 tasks identical to milestone 4a.
- wall time 160 s serial, 88 s with 2 or 3 tasks (the complex k-point dominates
  the load), 52 s with 4 tasks.

### Milestone 6 (2026-11-25): forces and stress
Code: FOCK$FORCE, FOCK$NEWCELL, FOCK_QSHAPES, FOCK_PAIRSET, FOCK_PAIRQ,
FOCK_ADDQL, FOCK_SYLMCHECK and the dFG/d|K| table in FOCK_SHAPEFT
(vespaw_fock.f90). The Fock projector term in WAVES$FORCE, and the calls of
FOCK$NEWCELL and FOCK$FORCE in WAVES$ETOT (vespaw_waves1.f90).
- The energy is E = sum_pairs lambda_mn t_mn with t_mn = Omega sum_K v(K)|rho_mn(K)|^2,
  lambda = -(AEXX/2) w_m w_n, doubled for m < n.
  rho_mn = conj(u_m) u_n exp(iG0 r) + rho^_mn, and the compensation part is built
  from Q^aL = sum_ij conj(<p_i|m>) <p_j|n> T^L_ij.
- Force, part 1 (projector derivatives): FOCK$APPLY gives the projector part of
  the Fock operator on each band (HPJ, exact on the band span through ACE).
  WAVES$FORCE adds f_n HPJ to DEDPROJ through WAVES_FORCE_ADDHTBC, the same
  routine as the onsite-hybrid term.
- Force, part 2 (explicit compensation charges): FOCK$FORCE loops over the pairs,
  rebuilds rho_mn(K), and accumulates Z^aL(K) = sum lambda conj(v rho) Q^aL with
  one ZGEMM per block of pairs. Then
  F_a = -2 Re sum_K (-iK) exp(-iK R_a) (-i)^L F_L(|K|) Y_L(K) Z^aL(K).
- Stress (+dE/deps at fixed sqrt(Omega) psi(G)), the sum of:
  - volume: -delta_ij E
  - kernel: sum lambda Omega (-K dv/dK) K_i K_j/K^2 |rho|^2, where
    -K dv/dK = 2v (Coulomb), or 2v - (2pi/mu^2) exp(-K^2/4mu^2) for erfc
  - G=0: sum lambda Omega |rho(0)|^2 dV0/deps_ij. dV0 is a central difference
    (h=1e-4) of FOCK_V0 on the strained Born-von Karman supercell.
  - compensation shapes: -2 Re sum exp(-iKR)(-i)^L [F'_L K_i K_j/K Y_L + F_L sym(K_i d_j Y_L)] Z
    with dF_L/dK = 4pi int g_L r^3 j'_L(Kr) dr (second table in FOCK_SHAPEFT, with
    the same interpolation as F_L).
- FOCK$NEWCELL recomputes V0, the G vectors, the kernels, Ylm, F_L, F'_L of every
  Q set after a cell change (it is called from WAVES$ETOT after WAVES_UPDATEGSET).
- MPI: pairs n split over 'KX', G vectors over 'K', sums over 'MONOMER'.
- Cost: computed in every SCF iteration when FORCES=T (as are the other force
  terms). Si2 PBE0 2x2x2: FOCK$FORCE 3.5 s vs 10.4 s for the build, about 19 % of
  the run. Could later be fused into FOCK$WBUILD.
- VESPAW_FOCK_CHECK=1 prints the pair energy vs E_x, the kernel, shape and G=0
  stress parts, and the FD checks of sym(K dY) and dF_L/dK.
Tests (regression_hard/fock_m6, Si2 diamond, SI2 displaced 0.05 AA along x,
FD steps dx=0.01 AA and eps=0.003, CDUAL 4):
- self-checks: E(pairs)-E_x 1e-13 to 6e-13 H, Ylm strain derivative vs FD 5e-8,
  dF_L/dK vs FD 5e-9.
- forces (analytic minus FD, eV/AA): PBE0 2x2x2 -4.8e-5, PBE -5.9e-5 (FD noise,
  the same for both). PBE0 4x1x1 at 50 Ry -6.0e-5, HSE06 2x2x2 -6.1e-5,
  spin-polarized PBE0 -4.8e-5. Newton's third law holds to 4e-6 (2x2x2).
- stress (analytic minus FD, eV):
  - 2x2x2 at 50 Ry, e23: PBE0 1.5e-6, HSE06 9.8e-6, PBE 2.1e-5.
  - e11 is off by 5.8e-4 (PBE0, HSE06) vs 6.6e-4 (PBE). This is the Pulay stress
    of the fixed-cutoff FD, and the Fock part matches to 7e-5.
  - at 25 Ry the PBE0 shear is off by 1.3e-3 vs 7e-4 for PBE. This difference is
    gone at 50 Ry, so it is a basis effect.
  - 4x1x1 at 50 Ry, PBE0 minus PBE: e11 -1.0e-4, e23 3.4e-4, e12 5e-4. Plain PBE
    is itself off by 1.0e-2 in e12 on this low-symmetry mesh (Pulay, not Fock).
- L-BFGS (PBE0 2x2x2, 25 Ry): atoms only converges in 44 iterations (51 s), with
  max force 6e-4 eV/AA. Atoms and cell together, from a cell stretched by 2 %,
  converges in 83 iterations (115 s) to V = 43.11 AA^3 (PBE 44.21) and P = 3e-4 GPa.
  A fresh static run at the returned structure differs by 1.3e-4 eV and 0.004 GPa,
  the same size as for PBE (2.0e-4 eV, 0.007 GPa). That gap is the basis held
  fixed during the cell relaxation, so FOCK$NEWCELL is consistent.
- MPI: Si2 PBE0 4x1x1 strained, serial vs 4 tasks (2 groups of 2): energy, forces
  and stress identical to all printed digits. Wall time 16.3 s vs 11.0 s.
- Restrictions as in milestone 4b/5: NDIM = 1, SYMMETRY off.

### Milestone 7 (2026-11-25): space-group symmetry, noncollinear spinors, GAMMASUPER
Code: FOCK$SETKMESH (stars), FOCK_ROTINIT, FOCK_ORBTOBOX, FOCK_ROTORB,
FOCK_ROTPROJ, FOCK_SYMCHECK, FOCK_ONECENTERNC, the NDIM argument of
FOCK$SNAPSHOT/APPLY/WBUILD and the star-restricted pair loop of FOCK$FORCE
(vespaw_fock.f90). SYMMETRY$GETOPS (vespaw_kpoints.f90). Call sites and
WAVES_FOCKCHECK (vespaw_waves1.f90).

Symmetry (SYMMETRY=T):
- The full mesh is the union of the stars q = s M^-T k of the irreducible
  points, over the operations x' = M x + tau and s = +1 or -1 (time reversal).
  A listed point is the first member of its own star. Without symmetry the only
  operation is the identity, which reproduces the earlier TR-partner logic.
  FOCK$SETKMESH stops if two listed points are related, if the stars do not fill
  the regular mesh, or if a weight differs from the star size / mesh size.
- Every occupied orbital is stored for all members of its star with the weight
  divided by the star size. The rotated orbital psi'(x) = psi(M^-1(x - tau)) has
  Bloch vector k' = M^-T k and periodic part
    u'(x) = sum_h c_h exp(-2 pi i (k'+h').tau) exp(2 pi i h'.x),  h' = M^-T h,
  built in G space on the whole FFT box (assembled over 'K'), so no grid
  commensurability with tau is needed. Time reversal conjugates.
- Projections: atom J = MAPAT(I) at x_J = M x_I + tau + L receives
  <p_J|psi'> = exp(2 pi i k'.L) U <p_I|psi>, with U the real-harmonic rotation
  of SYMMETRY$DENMAT ((-1)^l for improper operations).
- ACE and FOCK$APPLY act only on the irreducible points. FOCK$FORCE takes n only
  at the irreducible points (all partners m, weight = star size). Forces and
  stress are then symmetrized by the driver, as all other terms.
- VESPAW_FOCK_CHECK=1: FOCK_SYMCHECK applies every operation of the little group
  of k to the occupied orbitals of k and checks U' = U A and P' = P A with one
  matrix A (relative residuals 1e-8 to 1e-11 for Si on 2x2x2 and 3x3x3, limited
  by the SCF convergence).
Spinors (NSPIN=3, NDIM=2):
- Spinor components are stacked (orbitals NRLF*2, projections NPROF*2). The pair
  density sums over the components, the operator acts on each component, the
  ACE vectors are spinors. Occupations are per spinor (no factor 1/2).
- One center: E = -AEXX/2 sum_ab D^ba_ik D^ab_lj U(i,j,k,l) with the 2x2 spin
  density matrix from (total, x, y, z) as LDAPLUSU_SPINDENMAT, DH by the chain
  rule (FOCK_ONECENTERNC).
- Noncollinear meshes are never reduced by time reversal, and the symmetry
  module refuses noncollinear spin, so no spinor rotation is needed.
GAMMASUPER: still refused with Fock. Plain PBE with GAMMASUPER=T is itself wrong
(Si2 Gamma: -81.19 eV after 100 steps vs -199.95590804 eV, run_m7d.sh), see the
open GAMMASUPER bug. The TINV pair packing already gives the Gamma savings.

Tests (regression_hard/fock_m7, Si2, 25 Ry, CDUAL 4):
- symmetry vs full mesh (run_m7a.sh), energies in eV:
  - displaced 0.05 AA, 2x2x2 PBE0, 8 operations: -214.35272168 both
  - pristine 3x3x3 PBE0, 48 operations, 4 of 27 k: -216.10871060 both
  - displaced 3x3x3 HSE06, 8 operations, 8 of 27 k: -215.96589795 / -215.96589796
  - forces 2e-5 eV/AA and stress 1e-4 GPa, about the noise of the full-mesh runs
- star-restricted FOCK$FORCE (run_m7c.sh): same forces and stress to 1e-6.
  Fock force time 34.5 -> 20.8 s for the HSE06 3x3x3 case.
- spinors vs collinear, displaced 2x2x2 PBE0 (run_m7c.sh): -214.35272168 both,
  forces 1e-5 eV/AA, stress 1e-4 GPa, operator/energy 4e-16, one-center
  E1 vs 1/2 sum DH D 1e-17, pair energy vs E_x 5e-13.
- spin-orbit + PBE0 (run_m7b.sh): converges, Gamma valence band splits by
  49.5 meV (PBE 48 meV).
- MPI, 4 tasks vs serial (run_m7b.sh): symmetry 3x3x3 and spinors 2x2x2,
  energy, forces and stress identical to all printed digits.

### Milestone 8 (2026-11-25 to 2026-11-27): validation against VASP 6.4.2
Test decks and scripts in regression_hard/fock_m8 (mk.py, mkb.py, run_a.sh,
run_b.sh/run_b2.sh, diag/). VASP: PRECFOCK=Accurate, ISYM=3 (hybrids), full
mesh for the magnetic oxides. Compared quantity: the hybrid shift of the gap,
gap(hybrid) - gap(PBE) in the same code, so that setup differences of the
PBE gap drop out. Differences VESPAW - VASP:

| system, mesh | PBE gap VESPAW / VASP | PBE0 shift | HSE06 shift |
|---|---|---|---|
| Si, 4x4x4 | 0.693 / 0.693 | -12 meV (-0.9%) | -11 meV (-1.7%) |
| C, 4x4x4 | 4.444 / 4.431 | -7 meV (-0.3%) | -9 meV (-0.7%) |
| MgO, 4x4x4 | 4.731 / 4.725 | -33 meV (-1.3%) | -34 meV (-1.9%) |
| MnO AFM-II, 2x2x2 | 0.684 / 0.786 | -10 meV (-0.3%) | -10 meV (-0.4%) |
| NiO AFM-II, 2x2x2 | 1.528 / 1.320 | -356 meV (-8.6%) | -360 meV (-10.6%) |

- Magnetic energies: PBE0 E(FM) - E(AFM-II) NiO 0.406 / 0.416 eV, MnO
  0.353 / 0.352 eV.
- The VASP MnO hybrids (ALGO=Damped from scratch) stop at NELM=200 with a
  drift of 1-2e-4 eV per step; their gaps may still move by a few meV.
- NiO: the deficit is the same for PBE0 and HSE06, so it is not the long-range
  part. Checks that exclude the Fock implementation:
  - VESPAW cutoff 40 -> 60 Ry: shift 3.795 -> 3.810 eV.
  - VASP Ni -> Ni_pv (3p in the valence): shift 4.151 -> 4.172 eV.
  - Core-valence exchange: CORE=F lowers the VESPAW NiO PBE0 gap by 0.35 eV,
    but frozen-core MnO (MN_HBS, 3s3p in the core, KCV) gives a shift of
    3.050 eV against 3.022 eV with the semicore setup, so the frozen-core term
    is consistent.
  - MnO with the semicore MY_MN setup agrees to 10 meV.
  The remaining difference is attributed to the Ni setup (the PBE gaps
  already differ by 0.21 eV): the hybrid shift of a localized d shell is
  sensitive to the d partial waves. Not yet tested: a Ni setup with a third
  d projector, VASP Ni_GW.
- Regression case testsuite/si_hybrid (PBE0 full mesh / symmetry / RMM,
  HSE06, finite-difference force), about 3 minutes.
- Cost (regression_hard/onsite_hybrids/timing/results_hse.txt, Apple M3, 8 MPI
  ranks, 4x4x4 = 36 k, CDUAL 2, ACE): HSE06/PBE per SCF iteration NiO 102,
  MnO 174 (VASP 4 ranks, PRECFOCK=Normal: 131 and 339).
