!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  NON-SPHERICAL (LM>1) ONE-CENTER SCAN EXCHANGE, VIA DIRECT ANGULAR **
!**  QUADRATURE (NOT A TAYLOR EXPANSION). SEE SCAN0r notes        **
!**  SCAN0r notes FOR THE DESIGN RATIONALE (WHY THIS, NOT   **
!**  MIRRORING LMTO_RADXC_WITHCUT's QUADRATIC-TAYLOR-EXPANSION/        **
!**  analytic-2nd-3rd-derivative MACHINERY -- AUGMENTATION_SCAN_       **
!**  POINTWISE HAS NO SUCH DERIVATIVES, DERIVING THEM WOULD BE A LARGE,**
!**  SEPARATE UNDERTAKING).                                            **
!**                                                                    **
!**  PORTED FROM SCAN0r notes                **
!**  nsph_scan_module.f90, VALIDATED THERE TWO WAYS BEFORE PORTING:    **
!**  (1) REDUCES EXACTLY (~1E-9) TO THE ALREADY-VALIDATED SPHERICAL-   **
!**  ONLY AUGMENTATION_SCAN_POINTWISE WHEN FED A PURELY SPHERICAL      **
!**  INPUT; (2) THE GENUINELY NEW PIECE (FULL 3D GRADIENT FROM AN LM   **
!**  EXPANSION) MATCHES THE INDEPENDENT TEXTBOOK SPHERICAL-COORDINATE  **
!**  GRADIENT FORMULA TO MACHINE PRECISION (6.4E-14) FOR A GENUINELY   **
!**  ANISOTROPIC TEST DENSITY.                                         **
!**                                                                    **
!**  ADDED RELATIVE TO THE PROTOTYPE: CUT(R) IS NOW A DIRECT INPUT     **
!**  (SAME PURELY-RADIAL, LM=1-DERIVED CUTOFF FUNCTION USED THROUGHOUT **
!**  THIS PROJECT'S ONSITE DOUBLE-COUNTING MACHINERY -- IT MULTIPLIES  **
!**  UNIFORMLY ACROSS ALL LM CHANNELS SINCE IT IS ITSELF ANGLE-        **
!**  INDEPENDENT), AND VCUT(R) (THE UNCUT, ANGULARLY-INTEGRATED RADIAL **
!**  ENERGY PROFILE, BEFORE CUT-WEIGHTING) IS RETURNED SEPARATELY --   **
!**  BOTH MATCHING LMTO_RADXCTAU_WITHCUT's OWN CONVENTION EXACTLY, SO  **
!**  THE TWO ROUTINES' VCUT OUTPUTS ARE DIRECTLY INTERCHANGEABLE (E.G. **
!**  FOR THE FROZEN-CORE SUBTRACTION, WHICH STILL USES THE SPHERICAL-  **
!**  ONLY ROUTINE SINCE THE FROZEN CORE HAS NO REAL ANGULAR STRUCTURE).**
!***********************************************************************
MODULE NSPH_SCAN_MODULE
IMPLICIT NONE
LOGICAL(4)             :: TINI=.FALSE.
INTEGER(4)             :: LXSAVE=-1
INTEGER(4)             :: NTHSAVE=-1,NPHSAVE=-1
REAL(8)   ,ALLOCATABLE :: QUADX(:),QUADY(:),QUADZ(:),QUADW(:)  !(NQUAD) UNIT VECTORS + WEIGHTS
REAL(8)   ,ALLOCATABLE :: QUADYLM(:,:)      !(NQUAD,LMX) Y_LM AT EACH QUADRATURE POINT
REAL(8)   ,ALLOCATABLE :: QUADDYLMX(:,:)    !(NQUAD,LMX) TANGENTIAL GRADIENT (X COMPONENT)
REAL(8)   ,ALLOCATABLE :: QUADDYLMY(:,:)
REAL(8)   ,ALLOCATABLE :: QUADDYLMZ(:,:)
END MODULE NSPH_SCAN_MODULE
!
!     ..................................................................
      SUBROUTINE NSPH_SCAN_BUILDGRID(LX,NTH,NPH)
!     ******************************************************************
!     **  BUILDS/CACHES THE ANGULAR QUADRATURE GRID (GAUSS-LEGENDRE IN **
!     **  COS(THETA), UNIFORM IN PHI -- SAME CONVENTION AS             **
!     **  SCANTAU_BUILDTABLES, paw_scantau.f90) AND, AT EACH GRID      **
!     **  POINT, THE REAL SOLID-HARMONIC VALUE Y_LM AND TANGENTIAL     **
!     **  GRADIENT [GRAD_OMEGA Y_LM] FOR LM=1..(LX+1)**2. ONE-TIME COST**
!     **  (TINI GUARD), MIRRORING SCANTAU_BUILDTABLES'S OWN PATTERN.   **
!     ******************************************************************
      USE NSPH_SCAN_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LX
      INTEGER(4),INTENT(IN) :: NTH,NPH
      INTEGER(4)              :: LMX,INDX
      REAL(8)   ,ALLOCATABLE  :: YLMPOL(:,:)
      REAL(8)                  :: COSTH(NTH),WTH(NTH)
      REAL(8)                  :: PI,DELTAPHI,PHI,SINTH
      REAL(8)                  :: XPT,YPT,ZPT,WEIGHT
      INTEGER(4)                :: ITH,IPH,IQ,NQUAD,LM,L
      REAL(8)                    :: SVAL,DSX,DSY,DSZ
!     ******************************************************************
      IF(TINI.AND.LX.LE.LXSAVE.AND.NTH.LE.NTHSAVE.AND.NPH.LE.NPHSAVE) RETURN
      IF(ALLOCATED(QUADX)) THEN
        DEALLOCATE(QUADX,QUADY,QUADZ,QUADW,QUADYLM,QUADDYLMX,QUADDYLMY,QUADDYLMZ)
      END IF
      LMX=(LX+1)**2
      NQUAD=NTH*NPH
      ALLOCATE(QUADX(NQUAD),QUADY(NQUAD),QUADZ(NQUAD),QUADW(NQUAD))
      ALLOCATE(QUADYLM(NQUAD,LMX),QUADDYLMX(NQUAD,LMX))
      ALLOCATE(QUADDYLMY(NQUAD,LMX),QUADDYLMZ(NQUAD,LMX))
!
      INDX=(LX+1)*(LX+2)*(LX+3)/6
      ALLOCATE(YLMPOL(INDX,LMX))
      CALL SPHERICAL$YLMPOLYNOMIALS(LX,INDX,LMX,YLMPOL)
!
      PI=4.D0*ATAN(1.D0)
      CALL SCANTAU_GAULEG(-1.D0,1.D0,NTH,COSTH,WTH)
      DELTAPHI=2.D0*PI/REAL(NPH,KIND=8)
      IQ=0
      DO ITH=1,NTH
        SINTH=SQRT(MAX(0.D0,1.D0-COSTH(ITH)**2))
        DO IPH=1,NPH
          IQ=IQ+1
          PHI=REAL(IPH-1,KIND=8)*DELTAPHI
          XPT=SINTH*COS(PHI); YPT=SINTH*SIN(PHI); ZPT=COSTH(ITH)
          WEIGHT=WTH(ITH)*DELTAPHI
          QUADX(IQ)=XPT; QUADY(IQ)=YPT; QUADZ(IQ)=ZPT; QUADW(IQ)=WEIGHT
          DO LM=1,LMX
            CALL SCANTAU_SOLIDYLM(LX,YLMPOL,INDX,LMX,LM,XPT,YPT,ZPT,SVAL,DSX,DSY,DSZ)
            CALL SCANTAU_LVALOFLM(LM,L)
            QUADYLM(IQ,LM)=SVAL
            QUADDYLMX(IQ,LM)=DSX-REAL(L,KIND=8)*SVAL*XPT
            QUADDYLMY(IQ,LM)=DSY-REAL(L,KIND=8)*SVAL*YPT
            QUADDYLMZ(IQ,LM)=DSZ-REAL(L,KIND=8)*SVAL*ZPT
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(YLMPOL)
      LXSAVE=LX; NTHSAVE=NTH; NPHSAVE=NPH
      TINI=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_SCAN_NSPH(GID,NR,LX,NTH,NPH,RHOIN,TAUIN,CUT &
     &                                  ,TCORR,EXC,VXCOUT,VTAUOUT,VCUT)
!     ******************************************************************
!     **  FULL (NOT SPHERICAL-ONLY) ONE-CENTER SCAN EXCHANGE+CORRELATION,**
!     **  VIA REAL ANGULAR QUADRATURE (SEE MODULE HEADER). RHOIN/TAUIN  **
!     **  AND THE OUTPUT POTENTIALS ALL USE THE SAME RAW LM-COEFFICIENT **
!     **  CONVENTION AS TAUL/RHO_ALL ELSEWHERE IN THIS PROJECT (VALUE   **
!     **  AT (R,OMEGA) = SUM_LM COEF(R,LM)*Y_LM(OMEGA)).               **
!     **  CUT(R) IS THE SAME PURELY-RADIAL CUTOFF FUNCTION EVERY OTHER **
!     **  ONSITE DOUBLE-COUNTING ROUTINE IN THIS PROJECT USES --       **
!     **  APPLIED UNIFORMLY ACROSS ALL LM CHANNELS (IT HAS NO ANGULAR  **
!     **  DEPENDENCE). VCUT(R) IS THE UNCUT, ANGULARLY-INTEGRATED      **
!     **  RADIAL ENERGY PROFILE (BEFORE CUT-WEIGHTING) -- SAME ROLE,   **
!     **  SAME CONVENTION, AS LMTO_RADXCTAU_WITHCUT's OWN VCUT.        **
!     **                                                               **
!     **  CORRELATION (ADDED 2026-08-28, ONSITE NON-SPHERICITY PASS):  **
!     **  TCORR IS CALLER-SUPPLIED (SAME CONVENTION AS THE SPHERICAL   **
!     **  AUGMENTATION_SCAN_POINTWISE, paw_scantau.f90 -- NOT          **
!     **  DFT$GETL4('CORRELATION',...), SINCE THE HYBRID DOUBLE-       **
!     **  COUNTING CALL FORCES THAT GLOBAL FLAG OFF VIA XCONLY FOR     **
!     **  EVERY DCTYPE). RS/T2 USE THE SAME CONSTANTS/FORMULAS AS      **
!     **  AUGMENTATION_SCAN_POINTWISE'S OWN CORRELATION BLOCK, NOW     **
!     **  EVALUATED POINTWISE AT EACH (R,ANGULAR QUADRATURE POINT)     **
!     **  PAIR, EXACTLY LIKE THE EXCHANGE PIECE ALREADY WAS.           **
!     ******************************************************************
      USE NSPH_SCAN_MODULE
      USE SCAN_X_MODULE, ONLY: SCAN_X_FX
      USE SCAN_C_MODULE, ONLY: SCAN_C_EPSC
      USE DEORB_ALPHA_MODULE, ONLY: DEORB_SQ
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID,NR,LX,NTH,NPH
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,(LX+1)**2)
      REAL(8)   ,INTENT(IN) :: TAUIN(NR,(LX+1)**2)
      REAL(8)   ,INTENT(IN) :: CUT(NR)
      LOGICAL(4),INTENT(IN) :: TCORR
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXCOUT(NR,(LX+1)**2)
      REAL(8)   ,INTENT(OUT):: VTAUOUT(NR,(LX+1)**2)
      REAL(8)   ,INTENT(OUT):: VCUT(NR)
      REAL(8)   ,PARAMETER  :: EXFAC=-0.75D0*(3.D0/3.14159265358979D0)**(1.D0/3.D0)
      REAL(8)   ,PARAMETER  :: RHOMIN=1.D-8
      INTEGER(4)              :: LMX,NQUAD,IQ,LM,IR
      REAL(8)                  :: PI,CTF,RSFAC,KFFAC,KSFAC,TFAC
      REAL(8)                  :: R(NR)
      REAL(8)                  :: RHOT(NR),TAUT(NR),DRHODR(NR)
      REAL(8)                  :: GX(NR),GY(NR),GZ(NR),GRHO2(NR)
      REAL(8)                  :: S(NR),ALPHA(NR),EPSX(NR),FX(NR),EPSC(NR),RS(NR),T2(NR)
      REAL(8)                  :: TAUUNIF(NR),TAUW(NR)
      REAL(8)                  :: RHOP(NR),RHOM(NR),TAUP(NR),TAUM(NR)
      REAL(8)                  :: SP(NR),SM(NR),ALPHAP(NR),ALPHAM(NR),FXP(NR),FXM(NR)
      REAL(8)                  :: EPSCP(NR),EPSCM(NR),RSP(NR),RSM(NR),T2P(NR),T2M(NR)
      REAL(8)                  :: EXCDENSP(NR),EXCDENSM(NR)
      REAL(8)                  :: VRHOQ(NR),VTAUQ(NR),EXCDENS(NR)
      REAL(8)                  :: DRHO(NR),DTAU(NR)
      REAL(8)                  :: EG2Q(NR),AUX(NR),AUX2(NR)
      REAL(8)                  :: ARAD(NR,(LX+1)**2),BANG(NR,(LX+1)**2)
      REAL(8)                  :: DRHOIN(NR,(LX+1)**2) ! RADIAL DERIVATIVES OF RHOIN
      REAL(8)                   :: SVAR,G2P,G2M,DG2,EP,EM
!     ******************************************************************
      LMX=(LX+1)**2
      CALL NSPH_SCAN_BUILDGRID(LX,NTH,NPH)
      NQUAD=NTH*NPH
!
      ARAD(:,:)=0.D0
      BANG(:,:)=0.D0
      CALL RADIAL$R(GID,NR,R)
      PI=4.D0*ATAN(1.D0)
      CTF=0.3D0*(3.D0*PI**2)**(2.D0/3.D0)
      RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0)
      KFFAC=(3.D0*PI**2)**(1.D0/3.D0)
      KSFAC=SQRT(4.D0*KFFAC/PI)
      TFAC=PI/(16.D0*KFFAC) ! T2 OF PBE/SCAN: KF=KFFAC*RHO**(1/3) (WAS 0.62x TOO SMALL)
!
      VXCOUT(:,:)=0.D0
      VTAUOUT(:,:)=0.D0
      VCUT(:)=0.D0
!
!     == THE ANGULAR POINTS ARE INDEPENDENT AND RUN IN PARALLEL (OPENMP).==
!     == RADIAL$DERIVE IS NOT THREAD SAFE; THE DERIVATIVE IS LINEAR, SO  ==
!     == THE RADIAL DERIVATIVES OF THE LM COEFFICIENTS ARE TAKEN ONCE    ==
!     == HERE. ONE SERIAL CALL OF EACH KERNEL BUILDS THEIR LOOKUP TABLES ==
!     == (TINI) BEFORE THE THREADS START.                                ==
      DO LM=1,LMX
        CALL RADIAL$DERIVE(GID,NR,RHOIN(:,LM),DRHOIN(:,LM))
      ENDDO
      CALL DEORB_SQ(1.D0,1.D0,0.D0,SVAR,EP)
      CALL SCAN_X_FX(1.D0,1.D0,EP)
      CALL SCAN_C_EPSC(1.D0,1.D0,1.D0,1.D0,EP)
!$omp parallel do SCHEDULE(DYNAMIC) DEFAULT(SHARED) &
!$omp PRIVATE(IQ,LM,IR,RHOT,TAUT,DRHODR,GX,GY,GZ,GRHO2,S,ALPHA,EPSX,FX,EPSC,RS,T2) &
!$omp PRIVATE(TAUUNIF,TAUW,RHOP,RHOM,TAUP,TAUM,SP,SM,ALPHAP,ALPHAM,FXP,FXM) &
!$omp PRIVATE(EPSCP,EPSCM,RSP,RSM,T2P,T2M,EXCDENSP,EXCDENSM,VRHOQ,VTAUQ,EXCDENS) &
!$omp PRIVATE(DRHO,DTAU,EG2Q,SVAR,G2P,G2M,DG2,EP,EM) &
!$omp REDUCTION(+:VXCOUT,VTAUOUT,VCUT,ARAD,BANG)
      DO IQ=1,NQUAD
!       == RECONSTRUCT RHO,TAU,GRAD-RHO AT THIS ANGULAR POINT, ALL R ====
        RHOT(:)=0.D0; TAUT(:)=0.D0; DRHODR(:)=0.D0
        GX(:)=0.D0; GY(:)=0.D0; GZ(:)=0.D0
        DO LM=1,LMX
          RHOT(:)=RHOT(:)+RHOIN(:,LM)*QUADYLM(IQ,LM)
          TAUT(:)=TAUT(:)+TAUIN(:,LM)*QUADYLM(IQ,LM)
          DRHODR(:)=DRHODR(:)+DRHOIN(:,LM)*QUADYLM(IQ,LM)
          GX(:)=GX(:)+RHOIN(:,LM)*QUADDYLMX(IQ,LM)
          GY(:)=GY(:)+RHOIN(:,LM)*QUADDYLMY(IQ,LM)
          GZ(:)=GZ(:)+RHOIN(:,LM)*QUADDYLMZ(IQ,LM)
        ENDDO
        DO IR=1,NR
          GX(IR)=DRHODR(IR)*QUADX(IQ)+GX(IR)/MAX(R(IR),1.D-12)
          GY(IR)=DRHODR(IR)*QUADY(IQ)+GY(IR)/MAX(R(IR),1.D-12)
          GZ(IR)=DRHODR(IR)*QUADZ(IQ)+GZ(IR)/MAX(R(IR),1.D-12)
        ENDDO
        GRHO2(:)=GX(:)**2+GY(:)**2+GZ(:)**2
!
!       == POINTWISE SCAN EXCHANGE+CORRELATION AT THIS ANGULAR POINT, ==
!       == ALL R =========================================================
        DO IR=1,NR
          IF(RHOT(IR).LT.RHOMIN) THEN
            EPSX(IR)=0.D0; EPSC(IR)=0.D0; VRHOQ(IR)=0.D0; VTAUQ(IR)=0.D0
            EXCDENS(IR)=0.D0
            CYCLE
          END IF
          TAUUNIF(IR)=CTF*RHOT(IR)**(5.D0/3.D0)
          TAUW(IR)=GRHO2(IR)/(8.D0*RHOT(IR))
          ALPHA(IR)=MAX((TAUT(IR)-TAUW(IR))/MAX(TAUUNIF(IR),1.D-16),0.D0)
          CALL DEORB_SQ(RHOT(IR),GRHO2(IR),0.D0,S(IR),SVAR)
          CALL SCAN_X_FX(S(IR),ALPHA(IR),FX(IR))
          EPSX(IR)=EXFAC*RHOT(IR)**(1.D0/3.D0)*FX(IR)
          EPSC(IR)=0.D0
          IF(TCORR) THEN
            RS(IR)=RSFAC/RHOT(IR)**(1.D0/3.D0)
            T2(IR)=TFAC*GRHO2(IR)*RHOT(IR)**(-7.D0/3.D0)
            CALL SCAN_C_EPSC(RS(IR),S(IR),T2(IR),ALPHA(IR),EPSC(IR))
          END IF
          EXCDENS(IR)=RHOT(IR)*(EPSX(IR)+EPSC(IR))
        ENDDO
!
!       == LOCAL (VALUE) POTENTIAL BY FD AT FIXED GRAD-RHO,TAU ==========
        DO IR=1,NR
          IF(RHOT(IR).LT.RHOMIN) CYCLE
          DRHO(IR)=MAX(1.D-6*RHOT(IR),1.D-10)
          RHOP(IR)=RHOT(IR)+DRHO(IR)
          RHOM(IR)=MAX(RHOT(IR)-DRHO(IR),RHOMIN)
          CALL DEORB_SQ(RHOP(IR),GRHO2(IR),0.D0,SP(IR),SVAR)
          ALPHAP(IR)=MAX((TAUT(IR)-GRHO2(IR)/(8.D0*RHOP(IR))) &
     &                    /MAX(CTF*RHOP(IR)**(5.D0/3.D0),1.D-16),0.D0)
          CALL SCAN_X_FX(SP(IR),ALPHAP(IR),FXP(IR))
          EXCDENSP(IR)=RHOP(IR)*EXFAC*RHOP(IR)**(1.D0/3.D0)*FXP(IR)
!
          CALL DEORB_SQ(RHOM(IR),GRHO2(IR),0.D0,SM(IR),SVAR)
          ALPHAM(IR)=MAX((TAUT(IR)-GRHO2(IR)/(8.D0*RHOM(IR))) &
     &                    /MAX(CTF*RHOM(IR)**(5.D0/3.D0),1.D-16),0.D0)
          CALL SCAN_X_FX(SM(IR),ALPHAM(IR),FXM(IR))
          EXCDENSM(IR)=RHOM(IR)*EXFAC*RHOM(IR)**(1.D0/3.D0)*FXM(IR)
!
          EPSCP(IR)=0.D0; EPSCM(IR)=0.D0
          IF(TCORR) THEN
            RSP(IR)=RSFAC/RHOP(IR)**(1.D0/3.D0)
            T2P(IR)=TFAC*GRHO2(IR)*RHOP(IR)**(-7.D0/3.D0)
            CALL SCAN_C_EPSC(RSP(IR),SP(IR),T2P(IR),ALPHAP(IR),EPSCP(IR))
            RSM(IR)=RSFAC/RHOM(IR)**(1.D0/3.D0)
            T2M(IR)=TFAC*GRHO2(IR)*RHOM(IR)**(-7.D0/3.D0)
            CALL SCAN_C_EPSC(RSM(IR),SM(IR),T2M(IR),ALPHAM(IR),EPSCM(IR))
            EXCDENSP(IR)=EXCDENSP(IR)+RHOP(IR)*EPSCP(IR)
            EXCDENSM(IR)=EXCDENSM(IR)+RHOM(IR)*EPSCM(IR)
          END IF
!
          VRHOQ(IR)=(EXCDENSP(IR)-EXCDENSM(IR))/(RHOP(IR)-RHOM(IR))
        ENDDO
!
!       == TAU POTENTIAL BY FD AT FIXED RHO,GRAD-RHO ====================
        DO IR=1,NR
          IF(RHOT(IR).LT.RHOMIN) CYCLE
          DTAU(IR)=MAX(1.D-6*ABS(TAUT(IR)),1.D-10)
          TAUP(IR)=TAUT(IR)+DTAU(IR)
          TAUM(IR)=TAUT(IR)-DTAU(IR)
          ALPHAP(IR)=MAX((TAUP(IR)-TAUW(IR))/MAX(TAUUNIF(IR),1.D-16),0.D0)
          CALL SCAN_X_FX(S(IR),ALPHAP(IR),FXP(IR))
          EXCDENSP(IR)=RHOT(IR)*EXFAC*RHOT(IR)**(1.D0/3.D0)*FXP(IR)
          ALPHAM(IR)=MAX((TAUM(IR)-TAUW(IR))/MAX(TAUUNIF(IR),1.D-16),0.D0)
          CALL SCAN_X_FX(S(IR),ALPHAM(IR),FXM(IR))
          EXCDENSM(IR)=RHOT(IR)*EXFAC*RHOT(IR)**(1.D0/3.D0)*FXM(IR)
!
          IF(TCORR) THEN
!           == RS,T2 DEPEND ONLY ON RHOT,GRHO2 -- HELD FIXED, REUSE THE ==
!           == UNPERTURBED RS(IR),T2(IR) FROM THE ENERGY BLOCK ABOVE   ==
            CALL SCAN_C_EPSC(RS(IR),S(IR),T2(IR),ALPHAP(IR),EPSCP(IR))
            CALL SCAN_C_EPSC(RS(IR),S(IR),T2(IR),ALPHAM(IR),EPSCM(IR))
            EXCDENSP(IR)=EXCDENSP(IR)+RHOT(IR)*EPSCP(IR)
            EXCDENSM(IR)=EXCDENSM(IR)+RHOT(IR)*EPSCM(IR)
          END IF
!
          VTAUQ(IR)=(EXCDENSP(IR)-EXCDENSM(IR))/(TAUP(IR)-TAUM(IR))
        ENDDO
!
!       == DERIVATIVE WITH RESPECT TO |GRAD RHO|**2 AT FIXED RHO,TAU (IT ==
!       == ENTERS THROUGH S, TAU_W AND T2). DE/D(GRAD RHO)=2*EG2*GRAD RHO ==
!       == GIVES THE GRADIENT TERM OF THE POTENTIAL, WHICH WAS MISSING   ==
!       == BEFORE 2026-09-27 (THE POTENTIAL WAS NOT THE DERIVATIVE OF    ==
!       == THE ENERGY).                                                  ==
        EG2Q(:)=0.D0
        DO IR=1,NR
          IF(RHOT(IR).LT.RHOMIN) CYCLE
          DG2=MAX(1.D-6*GRHO2(IR),1.D-10*RHOT(IR)**(8.D0/3.D0),1.D-30)
          G2P=GRHO2(IR)+DG2
          G2M=MAX(GRHO2(IR)-DG2,0.D0)
          CALL NSPH_SCAN_EDENS(RHOT(IR),G2P,TAUT(IR),TCORR,EP)
          CALL NSPH_SCAN_EDENS(RHOT(IR),G2M,TAUT(IR),TCORR,EM)
          EG2Q(IR)=(EP-EM)/(G2P-G2M)
        ENDDO
        DO LM=1,LMX
          ARAD(:,LM)=ARAD(:,LM)+QUADW(IQ)*CUT(:)*2.D0*EG2Q(:) &
     &              *(GX(:)*QUADX(IQ)+GY(:)*QUADY(IQ)+GZ(:)*QUADZ(IQ))*QUADYLM(IQ,LM)
          BANG(:,LM)=BANG(:,LM)+QUADW(IQ)*CUT(:)*2.D0*EG2Q(:) &
     &              *(GX(:)*QUADDYLMX(IQ,LM)+GY(:)*QUADDYLMY(IQ,LM)+GZ(:)*QUADDYLMZ(IQ,LM))
        ENDDO
!
!       == ACCUMULATE UNCUT ENERGY PROFILE AND PROJECT POTENTIALS ONTO LM ==
        VCUT(:)=VCUT(:)+QUADW(IQ)*EXCDENS(:)
        DO LM=1,LMX
          VXCOUT(:,LM)=VXCOUT(:,LM)+QUADW(IQ)*VRHOQ(:)*QUADYLM(IQ,LM)
          VTAUOUT(:,LM)=VTAUOUT(:,LM)+QUADW(IQ)*VTAUQ(:)*QUADYLM(IQ,LM)
        ENDDO
      ENDDO
!$omp end parallel do
!
!     == APPLY THE (PURELY RADIAL, ANGLE-INDEPENDENT) CUT WEIGHTING =====
!     == UNIFORMLY ACROSS ALL LM CHANNELS -- SAME PATTERN AS            ==
!     == LMTO_RADXCTAU_WITHCUT's OWN "VCUT=FXC; FXC=FXC*CUT;            ==
!     == VXC=VXCLOC*CUT" (paw_lmto.f90).                                ==
      DO LM=1,LMX
        VXCOUT(:,LM)=VXCOUT(:,LM)*CUT(:)
        VTAUOUT(:,LM)=VTAUOUT(:,LM)*CUT(:)
      ENDDO
!
!     == GRADIENT TERM (WEAK FORM): FOR DRHO=F(R)Y_LM,                   ==
!     == GRAD DRHO=F'Y_LM RHAT+(F/R) GRAD_OMEGA Y_LM, SO                 ==
!     == V_LM += -(1/R**2)(R**2 ARAD_LM)' + BANG_LM/R  (CUT INCLUDED)    ==
      DO LM=1,LMX
        AUX(:)=ARAD(:,LM)*R(:)**2
        CALL RADIAL$DERIVE(GID,NR,AUX,AUX2)
        VXCOUT(:,LM)=VXCOUT(:,LM)-AUX2(:)/R(:)**2+BANG(:,LM)/R(:)
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,VCUT(:)*CUT(:)*R(:)**2,EXC)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE NSPH_SCAN_EDENS(RHO,G2,TAU,TCORR,EXCDENS)
!     ******************************************************************
!     **  SCAN (ORBITAL-TAU) ENERGY DENSITY RHO*(EPS_X+EPS_C) AT ONE   **
!     **  POINT, SAME FORMULAS AS AUGMENTATION_SCAN_NSPH ABOVE.        **
!     ******************************************************************
      USE SCAN_X_MODULE, ONLY: SCAN_X_FX
      USE SCAN_C_MODULE, ONLY: SCAN_C_EPSC
      USE DEORB_ALPHA_MODULE, ONLY: DEORB_SQ
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RHO,G2,TAU
      LOGICAL(4),INTENT(IN) :: TCORR
      REAL(8)   ,INTENT(OUT):: EXCDENS
      REAL(8)   ,PARAMETER  :: EXFAC=-0.75D0*(3.D0/3.14159265358979D0)**(1.D0/3.D0)
      REAL(8)               :: PI,CTF,RSFAC,KFFAC,KSFAC,TFAC
      REAL(8)               :: S,SVAR,ALPHA,FX,EPSC,RS,T2
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      CTF=0.3D0*(3.D0*PI**2)**(2.D0/3.D0)
      RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0)
      KFFAC=(3.D0*PI**2)**(1.D0/3.D0)
      KSFAC=SQRT(4.D0*KFFAC/PI)
      TFAC=PI/(16.D0*KFFAC) ! T2 OF PBE/SCAN: KF=KFFAC*RHO**(1/3) (WAS 0.62x TOO SMALL)
      ALPHA=MAX((TAU-G2/(8.D0*RHO))/MAX(CTF*RHO**(5.D0/3.D0),1.D-16),0.D0)
      CALL DEORB_SQ(RHO,G2,0.D0,S,SVAR)
      CALL SCAN_X_FX(S,ALPHA,FX)
      EXCDENS=RHO*EXFAC*RHO**(1.D0/3.D0)*FX
      IF(TCORR) THEN
        RS=RSFAC/RHO**(1.D0/3.D0)
        T2=TFAC*G2*RHO**(-7.D0/3.D0)
        CALL SCAN_C_EPSC(RS,S,T2,ALPHA,EPSC)
        EXCDENS=EXCDENS+RHO*EPSC
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_SCANC_NSPH(GID,NR,LX,NTH,NPH,RHOT,RHOS,TAUT,CUT &
     &                                  ,EXC,VT,VS,VTAUT,VCUT)
!     ******************************************************************
!     **  ONE-CENTER SPIN-POLARIZED SCAN CORRELATION (SCAN_C_EDENS) OF **
!     **  A COLLINEAR SPIN DENSITY BY ANGULAR QUADRATURE.              **
!     **  INPUT AS RAW LM COEFFICIENTS: RHOT=RHOUP+RHODN,              **
!     **  RHOS=RHOUP-RHODN, TAUT=TAUUP+TAUDN. THE CORRELATION DEPENDS  **
!     **  ON RHOUP,RHODN, |GRAD RHOT|**2 AND TAUT.                      **
!     **  OUTPUT: VT=D(EXC)/D(RHOT), VS=D(EXC)/D(RHOS),                 **
!     **  VTAUT=D(EXC)/D(TAUT), CUT-WEIGHTED, AS LM COEFFICIENTS; THE   **
!     **  GRADIENT TERM IS INCLUDED IN VT. VCUT IS THE UNCUT,           **
!     **  ANGULARLY INTEGRATED ENERGY DENSITY, EXC=INT(VCUT*CUT R**2).  **
!     **  CONVENTIONS AS AUGMENTATION_SCAN_NSPH.                        **
!     ******************************************************************
      USE NSPH_SCAN_MODULE
      USE SCAN_C_MODULE, ONLY: SCAN_C_EDENS
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID,NR,LX,NTH,NPH
      REAL(8)   ,INTENT(IN) :: RHOT(NR,(LX+1)**2)
      REAL(8)   ,INTENT(IN) :: RHOS(NR,(LX+1)**2)
      REAL(8)   ,INTENT(IN) :: TAUT(NR,(LX+1)**2)
      REAL(8)   ,INTENT(IN) :: CUT(NR)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VT(NR,(LX+1)**2)
      REAL(8)   ,INTENT(OUT):: VS(NR,(LX+1)**2)
      REAL(8)   ,INTENT(OUT):: VTAUT(NR,(LX+1)**2)
      REAL(8)   ,INTENT(OUT):: VCUT(NR)
      REAL(8)   ,PARAMETER  :: RHOMIN=1.D-8
      INTEGER(4)            :: LMX,NQUAD,IQ,LM,IR
      REAL(8)               :: R(NR),NT(NR),NS(NR),TQ(NR),DNDR(NR)
      REAL(8)               :: GX(NR),GY(NR),GZ(NR),G2(NR)
      REAL(8)               :: EQ(NR),VUPQ(NR),VDNQ(NR),VTQ(NR),EG2Q(NR)
      REAL(8)               :: ARAD(NR,(LX+1)**2),BANG(NR,(LX+1)**2)
      REAL(8)               :: AUX(NR),AUX2(NR)
      REAL(8)               :: DRHOT(NR,(LX+1)**2) ! RADIAL DERIVATIVES OF RHOT
      REAL(8)               :: NUP,NDN,D,EP,EM,G2P,G2M
!     ******************************************************************
      LMX=(LX+1)**2
      CALL NSPH_SCAN_BUILDGRID(LX,NTH,NPH)
      NQUAD=NTH*NPH
      CALL RADIAL$R(GID,NR,R)
      VT(:,:)=0.D0
      VS(:,:)=0.D0
      VTAUT(:,:)=0.D0
      VCUT(:)=0.D0
      ARAD(:,:)=0.D0
      BANG(:,:)=0.D0
!     == PARALLEL OVER ANGULAR POINTS, AS IN AUGMENTATION_SCAN_NSPH ======
      DO LM=1,LMX
        CALL RADIAL$DERIVE(GID,NR,RHOT(:,LM),DRHOT(:,LM))
      ENDDO
      CALL SCAN_C_EDENS(1.D0,1.D0,1.D0,1.D0,EP)
!$omp parallel do SCHEDULE(DYNAMIC) DEFAULT(SHARED) &
!$omp PRIVATE(IQ,LM,IR,NT,NS,TQ,DNDR,GX,GY,GZ,G2,EQ,VUPQ,VDNQ,VTQ,EG2Q) &
!$omp PRIVATE(NUP,NDN,D,EP,EM,G2P,G2M) &
!$omp REDUCTION(+:VT,VS,VTAUT,VCUT,ARAD,BANG)
      DO IQ=1,NQUAD
!       == DENSITIES, TAU AND GRAD RHOT AT THIS ANGULAR POINT, ALL R ===
        NT(:)=0.D0; NS(:)=0.D0; TQ(:)=0.D0; DNDR(:)=0.D0
        GX(:)=0.D0; GY(:)=0.D0; GZ(:)=0.D0
        DO LM=1,LMX
          NT(:)=NT(:)+RHOT(:,LM)*QUADYLM(IQ,LM)
          NS(:)=NS(:)+RHOS(:,LM)*QUADYLM(IQ,LM)
          TQ(:)=TQ(:)+TAUT(:,LM)*QUADYLM(IQ,LM)
          DNDR(:)=DNDR(:)+DRHOT(:,LM)*QUADYLM(IQ,LM)
          GX(:)=GX(:)+RHOT(:,LM)*QUADDYLMX(IQ,LM)
          GY(:)=GY(:)+RHOT(:,LM)*QUADDYLMY(IQ,LM)
          GZ(:)=GZ(:)+RHOT(:,LM)*QUADDYLMZ(IQ,LM)
        ENDDO
        DO IR=1,NR
          GX(IR)=DNDR(IR)*QUADX(IQ)+GX(IR)/MAX(R(IR),1.D-12)
          GY(IR)=DNDR(IR)*QUADY(IQ)+GY(IR)/MAX(R(IR),1.D-12)
          GZ(IR)=DNDR(IR)*QUADZ(IQ)+GZ(IR)/MAX(R(IR),1.D-12)
        ENDDO
        G2(:)=GX(:)**2+GY(:)**2+GZ(:)**2
!
!       == ENERGY DENSITY AND ITS PARTIAL DERIVATIVES ==================
        DO IR=1,NR
          EQ(IR)=0.D0; VUPQ(IR)=0.D0; VDNQ(IR)=0.D0; VTQ(IR)=0.D0; EG2Q(IR)=0.D0
          IF(NT(IR).LT.RHOMIN) CYCLE
          NUP=MAX(0.5D0*(NT(IR)+NS(IR)),0.D0)
          NDN=MAX(0.5D0*(NT(IR)-NS(IR)),0.D0)
          CALL SCAN_C_EDENS(NUP,NDN,G2(IR),TQ(IR),EQ(IR))
          D=MAX(1.D-6*NUP,1.D-10)
          CALL SCAN_C_EDENS(NUP+D,NDN,G2(IR),TQ(IR),EP)
          CALL SCAN_C_EDENS(MAX(NUP-D,0.D0),NDN,G2(IR),TQ(IR),EM)
          VUPQ(IR)=(EP-EM)/(NUP+D-MAX(NUP-D,0.D0))
          D=MAX(1.D-6*NDN,1.D-10)
          CALL SCAN_C_EDENS(NUP,NDN+D,G2(IR),TQ(IR),EP)
          CALL SCAN_C_EDENS(NUP,MAX(NDN-D,0.D0),G2(IR),TQ(IR),EM)
          VDNQ(IR)=(EP-EM)/(NDN+D-MAX(NDN-D,0.D0))
          D=MAX(1.D-6*ABS(TQ(IR)),1.D-10)
          CALL SCAN_C_EDENS(NUP,NDN,G2(IR),TQ(IR)+D,EP)
          CALL SCAN_C_EDENS(NUP,NDN,G2(IR),TQ(IR)-D,EM)
          VTQ(IR)=(EP-EM)/(2.D0*D)
          D=MAX(1.D-6*G2(IR),1.D-10*NT(IR)**(8.D0/3.D0),1.D-30)
          G2P=G2(IR)+D
          G2M=MAX(G2(IR)-D,0.D0)
          CALL SCAN_C_EDENS(NUP,NDN,G2P,TQ(IR),EP)
          CALL SCAN_C_EDENS(NUP,NDN,G2M,TQ(IR),EM)
          EG2Q(IR)=(EP-EM)/(G2P-G2M)
        ENDDO
!
!       == PROJECTIONS ONTO Y_LM ========================================
        VCUT(:)=VCUT(:)+QUADW(IQ)*EQ(:)
        DO LM=1,LMX
          VT(:,LM)=VT(:,LM)+QUADW(IQ)*0.5D0*(VUPQ(:)+VDNQ(:))*QUADYLM(IQ,LM)
          VS(:,LM)=VS(:,LM)+QUADW(IQ)*0.5D0*(VUPQ(:)-VDNQ(:))*QUADYLM(IQ,LM)
          VTAUT(:,LM)=VTAUT(:,LM)+QUADW(IQ)*VTQ(:)*QUADYLM(IQ,LM)
          ARAD(:,LM)=ARAD(:,LM)+QUADW(IQ)*CUT(:)*2.D0*EG2Q(:) &
     &              *(GX(:)*QUADX(IQ)+GY(:)*QUADY(IQ)+GZ(:)*QUADZ(IQ))*QUADYLM(IQ,LM)
          BANG(:,LM)=BANG(:,LM)+QUADW(IQ)*CUT(:)*2.D0*EG2Q(:) &
     &              *(GX(:)*QUADDYLMX(IQ,LM)+GY(:)*QUADDYLMY(IQ,LM)+GZ(:)*QUADDYLMZ(IQ,LM))
        ENDDO
      ENDDO
!$omp end parallel do
!
!     == CUT WEIGHTING AND GRADIENT TERM (AS AUGMENTATION_SCAN_NSPH) ====
      DO LM=1,LMX
        VT(:,LM)=VT(:,LM)*CUT(:)
        VS(:,LM)=VS(:,LM)*CUT(:)
        VTAUT(:,LM)=VTAUT(:,LM)*CUT(:)
        AUX(:)=ARAD(:,LM)*R(:)**2
        CALL RADIAL$DERIVE(GID,NR,AUX,AUX2)
        VT(:,LM)=VT(:,LM)-AUX2(:)/R(:)**2+BANG(:,LM)/R(:)
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,VCUT(:)*CUT(:)*R(:)**2,EXC)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SCANC_ONECENTER(GID,NR,LMX,NDIMD,RHOIN,TAUIN,CUT,EXC,VXC,VTAU,VCUT)
!     ******************************************************************
!     **  SPIN-POLARIZED SCAN CORRELATION OF A ONE-CENTER DENSITY IN   **
!     **  THE (TOTAL,SPIN) CONVENTION OF THE CALLERS (NDIMD=1 OR 2),   **
!     **  RAW LM COEFFICIENTS. ADDED BY THE TSCAN ROUTINES TO THEIR    **
!     **  SPIN-SCALED EXCHANGE. VXC,VTAU ARE CUT-WEIGHTED, VCUT UNCUT. **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID,NR,LMX,NDIMD
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMX,NDIMD)
      REAL(8)   ,INTENT(IN) :: TAUIN(NR,LMX,NDIMD)
      REAL(8)   ,INTENT(IN) :: CUT(NR)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMX,NDIMD)
      REAL(8)   ,INTENT(OUT):: VTAU(NR,LMX,NDIMD)
      REAL(8)   ,INTENT(OUT):: VCUT(NR)
      INTEGER(4),PARAMETER  :: NTH=14,NPH=24
      INTEGER(4)            :: LX
      REAL(8)               :: RHOS(NR,LMX),VS(NR,LMX)
!     ******************************************************************
      IF(NDIMD.NE.1.AND.NDIMD.NE.2) THEN
        CALL ERROR$MSG('ONLY NDIMD=1,2 ARE SUPPORTED')
        CALL ERROR$I4VAL('NDIMD',NDIMD)
        CALL ERROR$STOP('SCANC_ONECENTER')
      END IF
      LX=NINT(SQRT(DBLE(LMX)))-1
      RHOS(:,:)=0.D0
      IF(NDIMD.EQ.2) RHOS(:,:)=RHOIN(:,:,2)
      VXC(:,:,:)=0.D0
      VTAU(:,:,:)=0.D0
      CALL AUGMENTATION_SCANC_NSPH(GID,NR,LX,NTH,NPH,RHOIN(:,:,1),RHOS,TAUIN(:,:,1) &
     &                            ,CUT,EXC,VXC(:,:,1),VS,VTAU(:,:,1),VCUT)
      IF(NDIMD.EQ.2) VXC(:,:,2)=VS(:,:)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_XC_NSPH(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
!     ******************************************************************
!     **  ONE-CENTER LDA/GGA EXCHANGE AND CORRELATION WITHOUT THE     **
!     **  SECOND-ORDER TAYLOR EXPANSION IN THE NON-SPHERICAL DENSITY: **
!     **  THE DENSITY AND ITS FULL 3-D GRADIENT ARE EVALUATED ON AN   **
!     **  ANGULAR QUADRATURE GRID (GAUSS-LEGENDRE IN COS(THETA) TIMES **
!     **  A UNIFORM PHI GRID, AS IN AUGMENTATION_SCAN_NSPH) AND THE   **
!     **  POTENTIAL IS THE EXACT DERIVATIVE OF THE QUADRATURE ENERGY: **
!     **    V_LM = <Y_LM DE/DRHO> + <A.GRAD_OMEGA Y_LM>/R              **
!     **         - 1/R**2 D/DR [ R**2 <A_R Y_LM> ]                     **
!     **  WITH A=DE/D(GRAD RHO). COLLINEAR SPIN ONLY (NDIMD=1,2).     **
!     **  CONVENTION AS AUGMENTATION_XC: RHOIN(:,LM,1)=TOTAL,         **
!     **  RHOIN(:,LM,2)=SPIN DENSITY, VALUE=SUM_LM COEF*Y_LM.         **
!     **  SYNTHESIS AND PROJECTION ARE MATRIX PRODUCTS OVER ALL RADII.**
!     ******************************************************************
      USE NSPH_SCAN_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      INTEGER(4),SAVE       :: NTH=14,NPH=24
      LOGICAL(4),SAVE       :: TINIGRID=.FALSE.
      REAL(8)   ,PARAMETER  :: RHOMIN=1.D-10
      CHARACTER(16)         :: CENV
      LOGICAL(4)            :: TGRA
      INTEGER(4)            :: NSPIN,LX,NQ,IQ,IR,LM,IS
      REAL(8)               :: R(NR),RINV(NR),EOFR(NR),WORK(NR),WORK1(NR)
      REAL(8)   ,ALLOCATABLE:: DRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: Y(:,:),DX(:,:),DY(:,:),DZ(:,:)  ! (LMRX,NQ)
      REAL(8)   ,ALLOCATABLE:: RQ(:,:,:),DRQ(:,:,:)            ! (NR,NQ,NSPIN)
      REAL(8)   ,ALLOCATABLE:: GXQ(:,:,:),GYQ(:,:,:),GZQ(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VQ(:,:,:),AR(:,:,:),AX(:,:,:),AY(:,:,:),AZ(:,:,:)
      REAL(8)               :: RHAT(3),G1(3),G2(3),E,VT,VS,GT2,GS2,GST,A1(3),A2(3)
!     ******************************************************************
      IF(.NOT.TINIGRID) THEN
!       == EXPERT OVERRIDE OF THE ANGULAR GRID (CONVERGENCE TESTS) =====
        CENV=' '
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_XCNSPH_GRID',CENV)
        IF(LEN_TRIM(CENV).GT.0) READ(CENV,*) NTH,NPH
        TINIGRID=.TRUE.
      END IF
      EXC=0.D0
      VXC=0.D0
      IF(NDIMD.GT.2) THEN
        CALL ERROR$MSG('EXACT ANGULAR ONE-CENTER XC ONLY FOR NDIMD=1,2')
        CALL ERROR$STOP('AUGMENTATION_XC_NSPH')
      END IF
      NSPIN=NDIMD
      CALL DFT$GETL4('GC',TGRA)
      LX=NINT(SQRT(REAL(LMRX,KIND=8)))-1
      CALL NSPH_SCAN_BUILDGRID(LX,NTH,NPH)
      NQ=NTHSAVE*NPHSAVE      ! THE CACHED GRID MAY BE FINER THAN REQUESTED
      CALL RADIAL$R(GID,NR,R)
      RINV=0.D0
      WHERE(R.GT.1.D-8) RINV=1.D0/R
!     == ANGULAR TABLES (WEIGHTED COPIES FOR THE PROJECTION) ===========
      ALLOCATE(Y(LMRX,NQ),DX(LMRX,NQ),DY(LMRX,NQ),DZ(LMRX,NQ))
      Y=TRANSPOSE(QUADYLM(1:NQ,1:LMRX))
      DX=TRANSPOSE(QUADDYLMX(1:NQ,1:LMRX))
      DY=TRANSPOSE(QUADDYLMY(1:NQ,1:LMRX))
      DZ=TRANSPOSE(QUADDYLMZ(1:NQ,1:LMRX))
!     == DENSITY, RADIAL DERIVATIVE AND TANGENTIAL GRADIENT ON (R,Q) ====
      ALLOCATE(DRHO(NR,LMRX,NSPIN))
      ALLOCATE(RQ(NR,NQ,NSPIN),DRQ(NR,NQ,NSPIN))
      ALLOCATE(GXQ(NR,NQ,NSPIN),GYQ(NR,NQ,NSPIN),GZQ(NR,NQ,NSPIN))
      DO IS=1,NSPIN
        DO LM=1,LMRX
          CALL RADIAL$DERIVE(GID,NR,RHOIN(:,LM,IS),DRHO(:,LM,IS))
        ENDDO
        RQ(:,:,IS)=MATMUL(RHOIN(:,:,IS),Y)
        DRQ(:,:,IS)=MATMUL(DRHO(:,:,IS),Y)
        GXQ(:,:,IS)=MATMUL(RHOIN(:,:,IS),DX)
        GYQ(:,:,IS)=MATMUL(RHOIN(:,:,IS),DY)
        GZQ(:,:,IS)=MATMUL(RHOIN(:,:,IS),DZ)
      ENDDO
!     == FUNCTIONAL ON THE GRID; STORE WEIGHTED DERIVATIVES ===========
      ALLOCATE(VQ(NR,NQ,NSPIN),AR(NR,NQ,NSPIN))
      ALLOCATE(AX(NR,NQ,NSPIN),AY(NR,NQ,NSPIN),AZ(NR,NQ,NSPIN))
      VQ=0.D0; AR=0.D0; AX=0.D0; AY=0.D0; AZ=0.D0
      EOFR=0.D0
      G2=0.D0
      DO IQ=1,NQ
        RHAT=(/QUADX(IQ),QUADY(IQ),QUADZ(IQ)/)
        DO IR=1,NR
          IF(R(IR).LT.1.D-8) CYCLE
          IF(RQ(IR,IQ,1).LE.RHOMIN) CYCLE
          G1=DRQ(IR,IQ,1)*RHAT+(/GXQ(IR,IQ,1),GYQ(IR,IQ,1),GZQ(IR,IQ,1)/)*RINV(IR)
          IF(NSPIN.EQ.2) G2=DRQ(IR,IQ,2)*RHAT &
     &                     +(/GXQ(IR,IQ,2),GYQ(IR,IQ,2),GZQ(IR,IQ,2)/)*RINV(IR)
          IF(.NOT.TGRA) THEN
            G1=0.D0; G2=0.D0
          END IF
          IF(NSPIN.EQ.2) THEN
            CALL DFT(RQ(IR,IQ,1),RQ(IR,IQ,2),DOT_PRODUCT(G1,G1),DOT_PRODUCT(G2,G2) &
     &              ,DOT_PRODUCT(G1,G2),E,VT,VS,GT2,GS2,GST)
          ELSE
            CALL DFT(RQ(IR,IQ,1),0.D0,DOT_PRODUCT(G1,G1),0.D0,0.D0 &
     &              ,E,VT,VS,GT2,GS2,GST)
          END IF
          EOFR(IR)=EOFR(IR)+QUADW(IQ)*E
          A1=2.D0*GT2*G1+GST*G2
          VQ(IR,IQ,1)=QUADW(IQ)*VT
          AR(IR,IQ,1)=QUADW(IQ)*DOT_PRODUCT(A1,RHAT)
          AX(IR,IQ,1)=QUADW(IQ)*A1(1)*RINV(IR)
          AY(IR,IQ,1)=QUADW(IQ)*A1(2)*RINV(IR)
          AZ(IR,IQ,1)=QUADW(IQ)*A1(3)*RINV(IR)
          IF(NSPIN.EQ.2) THEN
            A2=2.D0*GS2*G2+GST*G1
            VQ(IR,IQ,2)=QUADW(IQ)*VS
            AR(IR,IQ,2)=QUADW(IQ)*DOT_PRODUCT(A2,RHAT)
            AX(IR,IQ,2)=QUADW(IQ)*A2(1)*RINV(IR)
            AY(IR,IQ,2)=QUADW(IQ)*A2(2)*RINV(IR)
            AZ(IR,IQ,2)=QUADW(IQ)*A2(3)*RINV(IR)
          END IF
        ENDDO
      ENDDO
!     == PROJECTION BACK ONTO Y_LM ======================================
      DO IS=1,NSPIN
        VXC(:,:,IS)=MATMUL(VQ(:,:,IS),TRANSPOSE(Y))
        IF(TGRA) THEN
          VXC(:,:,IS)=VXC(:,:,IS)+MATMUL(AX(:,:,IS),TRANSPOSE(DX)) &
     &                           +MATMUL(AY(:,:,IS),TRANSPOSE(DY)) &
     &                           +MATMUL(AZ(:,:,IS),TRANSPOSE(DZ))
          DRHO(:,:,IS)=MATMUL(AR(:,:,IS),TRANSPOSE(Y))    ! REUSED: <A_R Y_LM>
          DO LM=1,LMRX
            WORK=DRHO(:,LM,IS)*R**2
            CALL RADIAL$DERIVE(GID,NR,WORK,WORK1)
            VXC(:,LM,IS)=VXC(:,LM,IS)-WORK1*RINV**2
          ENDDO
        END IF
      ENDDO
!     == THE ORIGIN (IF ON THE GRID) TAKES THE VALUE OF THE NEXT POINT =
      IF(R(1).LT.1.D-8) VXC(1,:,:)=VXC(2,:,:)
      WORK=EOFR*R**2
      CALL RADIAL$INTEGRAL(GID,NR,WORK,EXC)
      DEALLOCATE(DRHO,Y,DX,DY,DZ,RQ,DRQ,GXQ,GYQ,GZQ,VQ,AR,AX,AY,AZ)
      RETURN
      END
