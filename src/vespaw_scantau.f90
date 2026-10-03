!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  MODULE/SUBROUTINE FOR THE ONE-CENTER TAU (KINETIC-ENERGY-DENSITY) **
!**  CONSTRUCTION NEEDED FOR TRUE SCAN-IN-PAW. NOT YET CALLED FROM     **
!**  ANYWHERE IN THE PRODUCTION CODE -- THIS IS THE PORT-THE-          **
!**  VALIDATED-PROTOTYPE-INTO-A-REAL-COMPILABLE-PRODUCTION-FILE        **
!**  STEP (SCAN0r notes), NOT YET THE HOOK-IT-    **
!**  INTO-THE-SCF-LOOP STEP.                                           **
!**                                                                    **
!**  MIRRORS AUGMENTATION_RHO OWN STRUCTURE AND CONVENTIONS             **
!**  (paw_augmentation.f90:600-653) AS CLOSELY AS POSSIBLE: SAME        **
!**  LMN1/LMN2 EXPANSION LOOP, SAME DENMAT(LMNXX,LMNXX) COMPLEX(8)      **
!**  ARGUMENT (REAL PART USED, SAME ASSUMES REAL PARTIAL WAVES          **
!**  CAVEAT), SAME PHI(NR,LNX) RADIAL-FUNCTION CONVENTION (CALL ONCE    **
!**  WITH AEPHI, AGAIN WITH PSPHI, EXACTLY LIKE AUGMENTATION_RHO OWN    **
!**  AERHO/PSRHO CALLS) -- BUT USING GID+RADIAL$DERIVE TO GET THE       **
!**  RADIAL DERIVATIVE INTERNALLY (LIKE AUGMENTATION_XC DOES FOR GRHO), **
!**  SINCE TAU NEEDS THE RADIAL DERIVATIVE OF PHI, NOT JUST PHI.        **
!**                                                                    **
!**  THE ANGULAR-MOMENTUM COUPLING TABLES (SCAN0r notes     **
!**  coupling_prototype/ylm_grad_module.f90 XTAB/YTAB) ARE REBUILT      **
!**  HERE USING CP-PAW OWN PRODUCTION SPHERICAL$YLMPOLYNOMIALS          **
!**  (paw_spherical.f90:152) AND GAUSSIAN_GAUSSINDEX (paw_gaussian.f90: **
!**  58) DIRECTLY -- NOT RE-PORTED COPIES -- SO THIS FILE REUSES REAL,  **
!**  ALREADY-TESTED PRODUCTION CODE FOR THE POLYNOMIAL/GAUNT MACHINERY. **
!**  ONLY THE GAUSS-LEGENDRE QUADRATURE (GAULEG) IS INCLUDED LOCALLY,   **
!**  SINCE NO EQUIVALENT EXISTS ELSEWHERE IN CP-PAW.                    **
!**                                                                    **
!**  DERIVATION OF THE TAU FORMULA (SAME AS ylm_grad_module.f90        **
!**  BUILD_TAU_COUPLING_TABLES HEADER, RE-STATED HERE FOR A READER OF   **
!**  THIS FILE WHO HAS NOT SEEN THE PROTOTYPE): FOR PHI_A(R,OMEGA)=     **
!**  R_A(R)*Y_A(OMEGA), GRAD PHI_A = D(R_A)/DR*Y_A*RHAT+(R_A(R)/R)*     **
!**  GRAD_OM Y_A (GRAD_OM = TANGENTIAL SURFACE GRADIENT, RHAT.GRAD_OM=0 **
!**  IDENTICALLY), SO GRAD PHI_A . GRAD PHI_B = DRA*DRB*Y_A*Y_B +      **
!**  (R_A*R_B/R**2)*(GRAD_OM Y_A . GRAD_OM Y_B) -- NO VALUE-TIMES-      **
!**  DERIVATIVE CROSS TERM. XTAB(J,K,I)=INT Y_I*Y_J*Y_K DOMEGA,         **
!**  YTAB(J,K,I)=INT Y_I*(GRAD_OM Y_J.GRAD_OM Y_K) DOMEGA, THE LATTER   **
!**  OBTAINED BY PROJECTING THE RADIAL COMPONENT OUT OF THE SOLID       **
!**  HARMONIC ORDINARY GRADIENT (GRAD_OM Y_LM=GRAD S_LM-L*Y_LM*RHAT,    **
!**  EULER THEOREM). THE RESULT (I) INDEX MUST COVER UP TO 2x THE       **
!**  CHANNEL L-RANGE (A REAL BUG CAUGHT WHILE BUILDING THE PROTOTYPE,   **
!**  SEE SCAN0r notes SECTION 5) -- LMXTAU BELOW, NOT    **
!**  JUST LMXCH.                                                       **
!***********************************************************************
MODULE SCANTAU_MODULE
IMPLICIT NONE
LOGICAL(4)            :: TINI=.FALSE.
INTEGER(4)            :: LXCHSAVE=-1
INTEGER(4)            :: LXTAUSAVE=-1
REAL(8)   ,ALLOCATABLE:: XTAB(:,:,:)   !(LMXCH,LMXCH,LMXTAU)
REAL(8)   ,ALLOCATABLE:: YTAB(:,:,:)   !(LMXCH,LMXCH,LMXTAU)
END MODULE SCANTAU_MODULE
!
!     ..................................................................
      SUBROUTINE SCANTAU_LVALOFLM(LM,L)
!     ******************************************************************
!     **  L OF INDEX LM=1+L**2+L-M (SAME CONVENTION AS SPHERICAL$YLM). **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM
      INTEGER(4),INTENT(OUT):: L
!     ******************************************************************
      IF(LM.LT.1) THEN
        CALL ERROR$MSG('LM MUST BE POSITIVE')
        CALL ERROR$I4VAL('LM',LM)
        CALL ERROR$STOP('SCANTAU_LVALOFLM')
      END IF
      L=INT(SQRT(REAL(LM-1,KIND=8)+1.D-9))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SCANTAU_GAULEG(X1,X2,N,X,W)
!     ******************************************************************
!     **  STANDARD GAUSS-LEGENDRE QUADRATURE NODES/WEIGHTS ON [X1,X2]. **
!     **  NO EQUIVALENT EXISTS ELSEWHERE IN CP-PAW (CHECKED BEFORE      **
!     **  WRITING THIS -- grep FOR GAULEG/GAUSSLEGENDRE FOUND NOTHING). **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: X1,X2
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(OUT):: X(N),W(N)
      REAL(8)   ,PARAMETER  :: EPS=3.D-14
      REAL(8)                :: PI,XM,XL,Z,Z1,P1,P2,P3,PP
      INTEGER(4)              :: M,I,J
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO I=1,M
        Z=COS(PI*(REAL(I,KIND=8)-0.25D0)/(REAL(N,KIND=8)+0.5D0))
        DO
          P1=1.D0; P2=0.D0
          DO J=1,N
            P3=P2; P2=P1
            P1=((2.D0*REAL(J,KIND=8)-1.D0)*Z*P2-(REAL(J,KIND=8)-1.D0)*P3)/REAL(J,KIND=8)
          ENDDO
          PP=REAL(N,KIND=8)*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
          IF(ABS(Z-Z1).LE.EPS) EXIT
        ENDDO
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SCANTAU_SOLIDYLM(LX,YLMPOL,INDX,LMX,LM,X,Y,Z,S,DSDX,DSDY,DSDZ)
!     ******************************************************************
!     **  EVALUATES THE SOLID HARMONIC S_LM(X,Y,Z)=R**L*Y_LM(OMEGA)   **
!     **  AND ITS EXACT GRADIENT FROM THE POLYNOMIAL COEFFICIENTS      **
!     **  ALREADY COMPUTED BY SPHERICAL$YLMPOLYNOMIALS -- PLAIN        **
!     **  TERM-BY-TERM POLYNOMIAL CALCULUS, NO NEW DERIVATION.         **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LX,INDX,LMX,LM
      REAL(8)   ,INTENT(IN) :: YLMPOL(INDX,LMX)
      REAL(8)   ,INTENT(IN) :: X,Y,Z
      REAL(8)   ,INTENT(OUT):: S,DSDX,DSDY,DSDZ
      INTEGER(4)              :: IND,I,J,K
      REAL(8)                  :: COEF,XI,YJ,ZK
!     ******************************************************************
      S=0.D0; DSDX=0.D0; DSDY=0.D0; DSDZ=0.D0
      DO IND=1,INDX
        COEF=YLMPOL(IND,LM)
        IF(COEF.EQ.0.D0) CYCLE
        CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
        XI=X**I; YJ=Y**J; ZK=Z**K
        S=S+COEF*XI*YJ*ZK
        IF(I.GT.0) DSDX=DSDX+COEF*REAL(I,KIND=8)*X**(I-1)*YJ*ZK
        IF(J.GT.0) DSDY=DSDY+COEF*REAL(J,KIND=8)*XI*Y**(J-1)*ZK
        IF(K.GT.0) DSDZ=DSDZ+COEF*REAL(K,KIND=8)*XI*YJ*Z**(K-1)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SCANTAU_BUILDTABLES(LXCH,LXTAU)
!     ******************************************************************
!     **  BUILDS/CACHES XTAB,YTAB (SCANTAU_MODULE) FOR CHANNEL RANGE   **
!     **  L=0..LXCH AND RESULT RANGE L=0..LXTAU. ONE-TIME COST (TINI   **
!     **  GUARD), MATCHING THE PROTOTYPE BUILD_TAU_COUPLING_TABLES.    **
!     ******************************************************************
      USE SCANTAU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LXCH,LXTAU
      INTEGER(4),PARAMETER   :: THPTS=14,PHPTS=24
      INTEGER(4)              :: LMXCH,LMXTAU,INDX
      REAL(8)   ,ALLOCATABLE  :: YLMPOL(:,:)
      REAL(8)                 :: COSTH(THPTS),WTH(THPTS)
      REAL(8)                 :: PI,DELTAPHI,PHI,SINTH
      REAL(8)                 :: XPT,YPT,ZPT,WEIGHT
      REAL(8)   ,ALLOCATABLE  :: SVAL(:),DSX(:),DSY(:),DSZ(:)
      REAL(8)   ,ALLOCATABLE  :: TX(:),TY(:),TZ(:)
      INTEGER(4)              :: ITH,IPH,I,J,K,L
      REAL(8)                  :: SDUM,DDX,DDY,DDZ
!     ******************************************************************
      IF(TINI.AND.LXCH.LE.LXCHSAVE.AND.LXTAU.LE.LXTAUSAVE) RETURN
      IF(ALLOCATED(XTAB)) DEALLOCATE(XTAB)
      IF(ALLOCATED(YTAB)) DEALLOCATE(YTAB)
      LMXCH=(LXCH+1)**2
      LMXTAU=(LXTAU+1)**2
      ALLOCATE(XTAB(LMXCH,LMXCH,LMXTAU))
      ALLOCATE(YTAB(LMXCH,LMXCH,LMXTAU))
      XTAB=0.D0; YTAB=0.D0
!
      INDX=(LXTAU+1)*(LXTAU+2)*(LXTAU+3)/6
      ALLOCATE(YLMPOL(INDX,LMXTAU))
      CALL SPHERICAL$YLMPOLYNOMIALS(LXTAU,INDX,LMXTAU,YLMPOL)
!
      ALLOCATE(SVAL(LMXTAU),DSX(LMXCH),DSY(LMXCH),DSZ(LMXCH))
      ALLOCATE(TX(LMXCH),TY(LMXCH),TZ(LMXCH))
      PI=4.D0*ATAN(1.D0)
      CALL SCANTAU_GAULEG(-1.D0,1.D0,THPTS,COSTH,WTH)
      DELTAPHI=2.D0*PI/REAL(PHPTS,KIND=8)
      DO ITH=1,THPTS
        SINTH=SQRT(MAX(0.D0,1.D0-COSTH(ITH)**2))
        DO IPH=1,PHPTS
          PHI=REAL(IPH-1,KIND=8)*DELTAPHI
          XPT=SINTH*COS(PHI); YPT=SINTH*SIN(PHI); ZPT=COSTH(ITH)
          WEIGHT=WTH(ITH)*DELTAPHI
          DO I=1,LMXTAU
            CALL SCANTAU_SOLIDYLM(LXTAU,YLMPOL,INDX,LMXTAU,I,XPT,YPT,ZPT,SVAL(I),SDUM,DDY,DDZ)
          ENDDO
          DO K=1,LMXCH
            CALL SCANTAU_SOLIDYLM(LXTAU,YLMPOL,INDX,LMXTAU,K,XPT,YPT,ZPT,SDUM,DSX(K),DSY(K),DSZ(K))
            CALL SCANTAU_LVALOFLM(K,L)
            TX(K)=DSX(K)-REAL(L,KIND=8)*SVAL(K)*XPT
            TY(K)=DSY(K)-REAL(L,KIND=8)*SVAL(K)*YPT
            TZ(K)=DSZ(K)-REAL(L,KIND=8)*SVAL(K)*ZPT
          ENDDO
          DO I=1,LMXTAU
            DO J=1,LMXCH
              DO K=1,LMXCH
                XTAB(J,K,I)=XTAB(J,K,I)+WEIGHT*SVAL(I)*SVAL(J)*SVAL(K)
                YTAB(J,K,I)=YTAB(J,K,I)+WEIGHT*SVAL(I)*(TX(J)*TX(K)+TY(J)*TY(K)+TZ(J)*TZ(K))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(YLMPOL,SVAL,DSX,DSY,DSZ,TX,TY,TZ)
      LXCHSAVE=LXCH; LXTAUSAVE=LXTAU
      TINI=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_TAU(GID,NR,LNX,LOX,PHI,LMNXX,DENMAT,LXCH,LXTAU,TAUL)
!     ******************************************************************
!     **  ONE-CENTER KINETIC-ENERGY DENSITY TAUL(R,LM), MIRRORING      **
!     **  AUGMENTATION_RHO OWN LMN1/LMN2 EXPANSION EXACTLY (SEE        **
!     **  paw_augmentation.f90:628-651) BUT CONTRACTING GRADIENT       **
!     **  PRODUCTS (VIA XTAB/YTAB) INSTEAD OF PLAIN PRODUCTS (VIA A    **
!     **  BARE GAUNT COEFFICIENT). CALL ONCE WITH PHI=AEPHI TO GET     **
!     **  TAU_1, AGAIN WITH PHI=PSPHI TO GET TAUTILDE_1 -- SAME        **
!     **  CALLING PATTERN AS AUGMENTATION_RHO OWN AERHO/PSRHO.           **
!     **                                                                **
!     **  LXCH MUST BE >= MAXVAL(LOX). LXTAU MUST BE >= 2*LXCH (THE    **
!     **  RESULT ANGULAR MOMENTUM RANGES UP TO L_A+L_B FOR A            **
!     **  PRODUCT OF TWO GRADIENTS) -- SEE THE MODULE HEADER.          **
!     ******************************************************************
      USE SCANTAU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! RADIAL GRID ID
      INTEGER(4),INTENT(IN) :: NR            ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LNX           ! #(RADIAL PARTIAL-WAVE CHANNELS)
      INTEGER(4),INTENT(IN) :: LOX(LNX)      ! ANGULAR MOMENTUM PER CHANNEL
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNX)   ! RADIAL PARTIAL WAVES (AE OR PS)
      INTEGER(4),INTENT(IN) :: LMNXX         ! MAX #(LM,N) SLOTS
      COMPLEX(8),INTENT(IN) :: DENMAT(LMNXX,LMNXX) ! ONE-CENTER DENSITY MATRIX
      INTEGER(4),INTENT(IN) :: LXCH          ! CHANNEL ANGULAR-MOMENTUM CUTOFF
      INTEGER(4),INTENT(IN) :: LXTAU         ! RESULT ANGULAR-MOMENTUM CUTOFF
      REAL(8)   ,INTENT(OUT):: TAUL(NR,(LXTAU+1)**2) ! ONE-CENTER TAU(R,LM)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: DPHI(:,:)
      INTEGER(4)              :: LMXTAU
      INTEGER(4)              :: LN1,LN2,L1,L2,IM1,IM2,LMN1,LMN2,LM1,LM2,LM3
      REAL(8)                  :: SVAR
!     ******************************************************************
      IF(LXTAU.LT.2*LXCH) THEN
        CALL ERROR$MSG('LXTAU MUST BE AT LEAST 2*LXCH')
        CALL ERROR$I4VAL('LXCH',LXCH)
        CALL ERROR$I4VAL('LXTAU',LXTAU)
        CALL ERROR$STOP('AUGMENTATION_TAU')
      END IF
      IF(MAXVAL(LOX).GT.LXCH) THEN
        CALL ERROR$MSG('LXCH TOO SMALL FOR THE CHANNELS SUPPLIED')
        CALL ERROR$I4VAL('LXCH',LXCH)
        CALL ERROR$I4VAL('MAXVAL(LOX)',MAXVAL(LOX))
        CALL ERROR$STOP('AUGMENTATION_TAU')
      END IF
      CALL SCANTAU_BUILDTABLES(LXCH,LXTAU)
      LMXTAU=(LXTAU+1)**2
!
      ALLOCATE(R(NR),DPHI(NR,LNX))
      CALL RADIAL$R(GID,NR,R)
      DO LN1=1,LNX
        CALL RADIAL$DERIVE(GID,NR,PHI(:,LN1),DPHI(:,LN1))
      ENDDO
!
      TAUL(:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LM1=L1**2+IM1
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
              SVAR=REAL(DENMAT(LMN1,LMN2),KIND=8)
              IF(SVAR.EQ.0.D0) CYCLE
!             == GUARD AGAINST THE SHLOG GRID EXACT R(1)=0: THE TRUE       ==
!             == PHI*PHI/R**2 TERM IS MATHEMATICALLY FINITE EVERYWHERE     ==
!             == (ANY S-CHANNEL (L=0) FACTOR MAKES YTAB IDENTICALLY ZERO,  ==
!             == SINCE THE TANGENTIAL GRADIENT OF A CONSTANT IS ZERO, AND  ==
!             == OTHERWISE PHI~R**L KEEPS THE RATIO REGULAR), BUT R(1)=0   ==
!             == EXACTLY WOULD GIVE A NUMERIC 0/0=NAN OR X/0=INF BEFORE    ==
!             == THE (EXACTLY ZERO) YTAB FACTOR CAN CANCEL IT.             ==
              DO LM3=1,LMXTAU
                TAUL(:,LM3)=TAUL(:,LM3)+SVAR*( &
     &                 DPHI(:,LN1)*DPHI(:,LN2)*XTAB(LM1,LM2,LM3) &
     &                +PHI(:,LN1)*PHI(:,LN2)/MAX(R(:)**2,1.D-16)*YTAB(LM1,LM2,LM3) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(R,DPHI)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_TAU_EXPECT(GID,NR,LNX,LOX,PHI,LMNXX &
     &                                   ,LXCH,LXTAU,VTAU,DATH)
!     ******************************************************************
!     **  ADJOINT OF AUGMENTATION_TAU: GIVEN THE TAU-POTENTIAL         **
!     **  VTAU(R,LM)=D(EXC)/D(TAU(R,LM)) (SAME (R,LM) CONVENTION AS    **
!     **  TAUL ABOVE, I.E. TAU(R,OMEGA)=SUM_LM VTAU_LM(R)*Y_LM(OMEGA)  **
!     **  FOR THE POTENTIAL TOO), RETURNS THE ONE-CENTER HAMILTONIAN   **
!     **  MATRIX ELEMENT DATH(LMN1,LMN2)=D(E_TAU)/D(DENMAT(LMN1,LMN2)) **
!     **  WHERE E_TAU=INT V_TAU(R,OMEGA)*TAU(R,OMEGA) D3R.             **
!     **                                                                **
!     **  DERIVATION: TAU_LM3(R) IS LINEAR IN DENMAT (SEE               **
!     **  AUGMENTATION_TAU ABOVE):                                      **
!     **    TAU_LM3(R)=SUM_LMN1,LMN2 D(LMN1,LMN2)*K_LMN1,LMN2,LM3(R)    **
!     **    K(R)=DPHI1*DPHI2*XTAB(LM1,LM2,LM3)                          **
!     **        +PHI1*PHI2/R**2*YTAB(LM1,LM2,LM3)                       **
!     **  SO, USING THE SAME (R,LM)-ORTHONORMAL-Y EXPANSION CONVENTION  **
!     **  AS AUGMENTATION_XC'S OWN POTENTIAL ARRAYS (E_TAU=SUM_LM       **
!     **  INT DR R**2 VTAU_LM(R)*TAU_LM(R)):                            **
!     **    D(E_TAU)/D(DENMAT(LMN1,LMN2))                               **
!     **      =SUM_LM3 INT DR R**2 VTAU_LM3(R)*K_LMN1,LMN2,LM3(R)       **
!     **  I.E. THE *SAME* LMN1/LMN2/LM3/XTAB/YTAB LOOP AS               **
!     **  AUGMENTATION_TAU, BUT INTEGRATING VTAU AGAINST THE KERNEL     **
!     **  INSTEAD OF CONTRACTING THE KERNEL WITH DENMAT -- THE EXACT    **
!     **  GRADIENT-GRADIENT COUNTERPART OF WHAT LMTO_EXPECT             **
!     **  (paw_lmto.f90:6869) ALREADY DOES FOR A VALUE-VALUE            **
!     **  (MULTIPLICATIVE) POTENTIAL. NO SUCH GRADIENT-GRADIENT ROUTINE **
!     **  EXISTED ANYWHERE IN CP-PAW BEFORE THIS -- A GENUINE GAP, NOT  **
!     **  JUST A PORT (SEE SCAN0r notes).         **
!     ******************************************************************
      USE SCANTAU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNX)
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: LXCH
      INTEGER(4),INTENT(IN) :: LXTAU
      REAL(8)   ,INTENT(IN) :: VTAU(NR,(LXTAU+1)**2)
      COMPLEX(8),INTENT(OUT):: DATH(LMNXX,LMNXX)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: DPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      INTEGER(4)              :: LMXTAU
      INTEGER(4)              :: LN1,LN2,L1,L2,IM1,IM2,LMN1,LMN2,LM1,LM2,LM3
      REAL(8)                  :: SVAR
!     ******************************************************************
      IF(LXTAU.LT.2*LXCH) THEN
        CALL ERROR$MSG('LXTAU MUST BE AT LEAST 2*LXCH')
        CALL ERROR$I4VAL('LXCH',LXCH)
        CALL ERROR$I4VAL('LXTAU',LXTAU)
        CALL ERROR$STOP('AUGMENTATION_TAU_EXPECT')
      END IF
      IF(MAXVAL(LOX).GT.LXCH) THEN
        CALL ERROR$MSG('LXCH TOO SMALL FOR THE CHANNELS SUPPLIED')
        CALL ERROR$I4VAL('LXCH',LXCH)
        CALL ERROR$I4VAL('MAXVAL(LOX)',MAXVAL(LOX))
        CALL ERROR$STOP('AUGMENTATION_TAU_EXPECT')
      END IF
      CALL SCANTAU_BUILDTABLES(LXCH,LXTAU)
      LMXTAU=(LXTAU+1)**2
!
      ALLOCATE(R(NR),DPHI(NR,LNX),WORK(NR))
      CALL RADIAL$R(GID,NR,R)
      DO LN1=1,LNX
        CALL RADIAL$DERIVE(GID,NR,PHI(:,LN1),DPHI(:,LN1))
      ENDDO
!
      DATH(:,:)=(0.D0,0.D0)
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LM1=L1**2+IM1
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
              SVAR=0.D0
              DO LM3=1,LMXTAU
                IF(XTAB(LM1,LM2,LM3).EQ.0.D0.AND.YTAB(LM1,LM2,LM3).EQ.0.D0) CYCLE
                WORK(:)=VTAU(:,LM3)*R(:)**2*( &
     &                 DPHI(:,LN1)*DPHI(:,LN2)*XTAB(LM1,LM2,LM3) &
     &                +PHI(:,LN1)*PHI(:,LN2)/MAX(R(:)**2,1.D-16)*YTAB(LM1,LM2,LM3) )
                CALL RADIAL$INTEGRAL(GID,NR,WORK,SVAR)
                DATH(LMN1,LMN2)=DATH(LMN1,LMN2)+CMPLX(SVAR,0.D0,KIND=8)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(R,DPHI,WORK)
      RETURN
      END
!
!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**  SUBROUTINE AUGMENTATION_SCAN_POINTWISE                            **
!**                                                                    **
!**  REAL (ORBITAL-TAU-BASED, NOT DEORBITALIZED) SCAN EXCHANGE-        **
!**  CORRELATION FOR THE ONE-CENTER SPHERICAL (LM=1) DENSITY -- THE    **
!**  TRUE-SCAN COUNTERPART OF AUGMENTATION_SCANL_POINTWISE ABOVE       **
!**  (paw_augmentation.f90:3810).                                      **
!**                                                                    **
!**  REUSES SCAN_X_MODULE::SCAN_X_FX (paw_augmentation.f90:3360)       **
!**  UNCHANGED -- THE ENHANCEMENT FACTOR Fx(s,alpha) DOES NOT CARE HOW **
!**  ALPHA WAS OBTAINED. THE ONLY DIFFERENCE FROM THE SCAN-L PATH IS   **
!**  ALPHA ITSELF: HERE ALPHA=(TAU-TAU_W)/TAU_UNIF IS BUILT FROM THE   **
!**  REAL ONE-CENTER TAU(R) (AUGMENTATION_TAU ABOVE), NOT FROM THE     **
!**  DENSITY-LAPLACIAN-BASED DEORB_ALPHA APPROXIMATION.                **
!**                                                                    **
!**  CORRELATION: REUSES SCAN_C_MODULE::SCAN_C_EPSC UNCHANGED, GATED   **
!**  ON THE CALLER-SUPPLIED TCORR ARGUMENT (NOT DFT$GETL4('CORRELATION'**
!**  ,...): THE PRODUCTION DOUBLE-COUNTING CALL FORCES THAT GLOBAL     **
!**  FLAG TO .FALSE. VIA XCONLY FOR *EVERY* DCTYPE, SO IT CANNOT BE    **
!**  USED TO TURN CORRELATION ON HERE -- SEE LMTO_MODULE'S SEPARATE    **
!**  LTSCANCORR FLAG, DEFAULT .FALSE., PLUMBED IN VIA LMTO_RADXCTAU_   **
!**  WITHCUT). SAME RS/T2 FORMULAS AS AUGMENTATION_SCANL_POINTWISE --  **
!**  BUT EVALUATED AT THE REAL-TAU ALPHA COMPUTED HERE, NOT THE        **
!**  DEORBITALIZED ONE.                                                **
!**  THIS IS THE ONLY DIFFERENCE FROM THE EXCHANGE-ONLY VERSION OF     **
!**  THIS ROUTINE: CORRELATION DEPENDS ON TAU ONLY THROUGH ALPHA (RS   **
!**  AND T2 DEPEND ONLY ON RHOT,GRHOT), SO THE VTAU FINITE DIFFERENCE  **
!**  BELOW ONLY NEEDS EPSC RE-EVALUATED AT ALPHA(TAU+/-DTAU) WITH      **
!**  RS,T2 HELD FIXED -- EXACTLY MIRRORING HOW EPSX'S OWN TAU-         **
!**  DEPENDENCE (THROUGH ALPHA, AT FIXED S) IS ALREADY HANDLED.        **
!**                                                                    **
!**  SCAN_X_MODULE'S COEFFICIENTS ARE INDEPENDENTLY VERIFIED AGAINST   **
!**  REAL LIBXC 7.0.0 mgga_x_scan TO MACHINE PRECISION ACROSS 42       **
!**  (S,ALPHA) POINTS (SCAN0r notes                   **
!**  verify_scanx_vs_libxc.f90). SCAN_C_MODULE'S COEFFICIENTS ARE      **
!**  INHERITED FROM THE ALREADY-ACCEPTED AUGMENTATION_SCANL_POINTWISE  **
!**  PRODUCTION PATH, NOT SEPARATELY RE-VERIFIED HERE.                 **
!**                                                                    **
!**  S IS BUILT VIA DEORB_ALPHA_MODULE::DEORB_SQ (PASSING A DUMMY      **
!**  LAPRHO=0 -- S's FORMULA DOES NOT DEPEND ON THE LAPLACIAN, ONLY Q  **
!**  DOES, AND Q IS DISCARDED) SO THE STANDARD REDUCED-GRADIENT        **
!**  FORMULA IS NOT RE-DERIVED, ONLY REUSED.                           **
!**                                                                    **
!**  VXCLOC AND VTAU ARE BOTH BUILT BY FINITE DIFFERENCE (VXCLOC AT    **
!**  FIXED GRHOT,TAU; VTAU AT FIXED RHOT,GRHOT) -- SAME METHOD, SAME   **
!**  *Y0*FOURPI CONVENTION, AS AUGMENTATION_SCANL_POINTWISE'S VXCLOC,  **
!**  DELIBERATELY, FOR CONSISTENCY WITH THAT ALREADY-ACCEPTED          **
!**  PRODUCTION CONVENTION RATHER THAN INTRODUCING A NEW ONE.          **
!***********************************************************************
      SUBROUTINE AUGMENTATION_SCAN_POINTWISE(GID,NR,RHOT,TAU,TCORR,FXC,VXCLOC,VTAU)
      USE SCAN_X_MODULE, ONLY: SCAN_X_FX
      USE SCAN_C_MODULE, ONLY: SCAN_C_EPSC
      USE DEORB_ALPHA_MODULE, ONLY: DEORB_SQ
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: RHOT(NR)   ! ACTUAL (Y0-SCALED) SPHERICAL DENSITY
      REAL(8)   ,INTENT(IN) :: TAU(NR)    ! REAL ONE-CENTER TAU, LM=1 (SPHERICAL) PART
      LOGICAL(4),INTENT(IN) :: TCORR      ! INCLUDE SCAN CORRELATION? CALLER-SUPPLIED
!                                         ! (NOT DFT$GETL4('CORRELATION',...): THE
!                                         ! PRODUCTION DOUBLE-COUNTING CALL FORCES THAT
!                                         ! GLOBAL FLAG TO .FALSE. VIA XCONLY FOR EVERY
!                                         ! DCTYPE, SO AN EXPLICIT ARGUMENT, SOURCED
!                                         ! FROM THE SEPARATE LMTO_MODULE LTSCANCORR
!                                         ! FLAG, IS THE ONLY WAY TO ACTUALLY ENABLE IT.
      REAL(8)   ,INTENT(OUT):: FXC(NR)    ! POINTWISE ENERGY DENSITY (FOURPI*(EPSX+EPSC)*RHOT)
      REAL(8)   ,INTENT(OUT):: VXCLOC(NR) ! D(EXCDENS)/D(RHOT) AT FIXED GRHOT,TAU
      REAL(8)   ,INTENT(OUT):: VTAU(NR)   ! D(EXCDENS)/D(TAU) AT FIXED RHOT,GRHOT
      REAL(8)   ,PARAMETER  :: EXFAC=-0.75D0*(3.D0/3.14159265358979D0)**(1.D0/3.D0)
      REAL(8)   ,PARAMETER  :: RHOMIN=1.D-8
      REAL(8)                :: PI,Y0,FOURPI,CTF
      REAL(8)                :: RSFAC,KFFAC,KSFAC,TFAC
      REAL(8)                :: GRHOT(NR)
      REAL(8)                :: S,ALPHA,EPSX,EPSC,EXCDENS,RS,T2
      REAL(8)                :: TAUW,TAUUNIF,DUMQ
      REAL(8)                :: RHOP,RHOM,SP,SM,ALPHAP,ALPHAM,EPSXP,EPSXM
      REAL(8)                :: EPSCP,EPSCM,RSP,RSM,T2P,T2M
      REAL(8)                :: EXCDENSP,EXCDENSM,TAUWP,TAUWM,TAUUNIFP,TAUUNIFM
      REAL(8)                :: TAUP,TAUM,DRHO,DTAU
      INTEGER(4)              :: IR
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CTF=0.3D0*(3.D0*PI**2)**(2.D0/3.D0)   ! SAME CONSTANT/FORMULA AS
                                             ! DEORB_ALPHA_MODULE's CTF
!     == SAME CONSTANTS/FORMULAS AS AUGMENTATION_SCANL_POINTWISE'S RS/T2 ==
      RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0)
      KFFAC=(3.D0*PI**2)**(1.D0/3.D0)
      KSFAC=SQRT(4.D0*KFFAC/PI)
      TFAC=PI/(16.D0*KFFAC) ! T2 OF PBE/SCAN: KF=KFFAC*RHO**(1/3) (WAS 0.62x TOO SMALL)
      FXC(:)=0.D0
      VXCLOC(:)=0.D0
      VTAU(:)=0.D0
      CALL RADIAL$DERIVE(GID,NR,RHOT,GRHOT)
!
!     ==================================================================
!     ==  POINTWISE SCAN EXCHANGE(+CORRELATION) ENERGY DENSITY        ==
!     ==================================================================
      DO IR=1,NR
        IF(RHOT(IR).LT.RHOMIN) CYCLE
        TAUUNIF=CTF*RHOT(IR)**(5.D0/3.D0)
        TAUW=GRHOT(IR)**2/(8.D0*RHOT(IR))
        ALPHA=MAX((TAU(IR)-TAUW)/MAX(TAUUNIF,1.D-16),0.D0)
        CALL DEORB_SQ(RHOT(IR),GRHOT(IR)**2,0.D0,S,DUMQ)
        CALL SCAN_X_FX(S,ALPHA,EPSX)
        EPSX=EXFAC*RHOT(IR)**(1.D0/3.D0)*EPSX
        EPSC=0.D0
        IF(TCORR) THEN
          RS=RSFAC/RHOT(IR)**(1.D0/3.D0)
          T2=TFAC*GRHOT(IR)**2*RHOT(IR)**(-7.D0/3.D0)
          CALL SCAN_C_EPSC(RS,S,T2,ALPHA,EPSC)
        END IF
        EXCDENS=RHOT(IR)*(EPSX+EPSC)
        FXC(IR)=FOURPI*EXCDENS
      ENDDO
!
!     ==================================================================
!     ==  VXCLOC: D(EXCDENS)/D(RHOT) AT FIXED GRHOT,TAU               ==
!     ==================================================================
      DO IR=1,NR
        IF(RHOT(IR).LT.RHOMIN) CYCLE
        DRHO=MAX(1.D-6*RHOT(IR),1.D-10)
        RHOP=RHOT(IR)+DRHO
        RHOM=MAX(RHOT(IR)-DRHO,RHOMIN)
!
        TAUUNIFP=CTF*RHOP**(5.D0/3.D0)
        TAUWP=GRHOT(IR)**2/(8.D0*RHOP)
        ALPHAP=MAX((TAU(IR)-TAUWP)/MAX(TAUUNIFP,1.D-16),0.D0)
        CALL DEORB_SQ(RHOP,GRHOT(IR)**2,0.D0,SP,DUMQ)
        CALL SCAN_X_FX(SP,ALPHAP,EPSXP)
        EPSXP=EXFAC*RHOP**(1.D0/3.D0)*EPSXP
        EPSCP=0.D0
        IF(TCORR) THEN
          RSP=RSFAC/RHOP**(1.D0/3.D0)
          T2P=TFAC*GRHOT(IR)**2*RHOP**(-7.D0/3.D0)
          CALL SCAN_C_EPSC(RSP,SP,T2P,ALPHAP,EPSCP)
        END IF
        EXCDENSP=RHOP*(EPSXP+EPSCP)
!
        TAUUNIFM=CTF*RHOM**(5.D0/3.D0)
        TAUWM=GRHOT(IR)**2/(8.D0*RHOM)
        ALPHAM=MAX((TAU(IR)-TAUWM)/MAX(TAUUNIFM,1.D-16),0.D0)
        CALL DEORB_SQ(RHOM,GRHOT(IR)**2,0.D0,SM,DUMQ)
        CALL SCAN_X_FX(SM,ALPHAM,EPSXM)
        EPSXM=EXFAC*RHOM**(1.D0/3.D0)*EPSXM
        EPSCM=0.D0
        IF(TCORR) THEN
          RSM=RSFAC/RHOM**(1.D0/3.D0)
          T2M=TFAC*GRHOT(IR)**2*RHOM**(-7.D0/3.D0)
          CALL SCAN_C_EPSC(RSM,SM,T2M,ALPHAM,EPSCM)
        END IF
        EXCDENSM=RHOM*(EPSXM+EPSCM)
!
        VXCLOC(IR)=(EXCDENSP-EXCDENSM)/(RHOP-RHOM)*Y0*FOURPI
      ENDDO
!
!     ==================================================================
!     ==  VTAU: D(EXCDENS)/D(TAU) AT FIXED RHOT,GRHOT                 ==
!     ==  (RS,T2 DEPEND ONLY ON RHOT,GRHOT -- HELD FIXED, RECOMPUTED   ==
!     ==  ONCE, REUSED FOR BOTH THE + AND - EVALUATIONS BELOW, EXACTLY ==
!     ==  MIRRORING HOW S IS ALREADY HANDLED FOR EPSX.)                ==
!     ==================================================================
      DO IR=1,NR
        IF(RHOT(IR).LT.RHOMIN) CYCLE
        TAUUNIF=CTF*RHOT(IR)**(5.D0/3.D0)
        TAUW=GRHOT(IR)**2/(8.D0*RHOT(IR))
        DTAU=MAX(1.D-6*ABS(TAU(IR)),1.D-10)
        TAUP=TAU(IR)+DTAU
        TAUM=TAU(IR)-DTAU
        CALL DEORB_SQ(RHOT(IR),GRHOT(IR)**2,0.D0,S,DUMQ)
        IF(TCORR) THEN
          RS=RSFAC/RHOT(IR)**(1.D0/3.D0)
          T2=TFAC*GRHOT(IR)**2*RHOT(IR)**(-7.D0/3.D0)
        END IF
!
        ALPHAP=MAX((TAUP-TAUW)/MAX(TAUUNIF,1.D-16),0.D0)
        CALL SCAN_X_FX(S,ALPHAP,EPSXP)
        EPSXP=EXFAC*RHOT(IR)**(1.D0/3.D0)*EPSXP
        EPSCP=0.D0
        IF(TCORR) CALL SCAN_C_EPSC(RS,S,T2,ALPHAP,EPSCP)
        EXCDENSP=RHOT(IR)*(EPSXP+EPSCP)
!
        ALPHAM=MAX((TAUM-TAUW)/MAX(TAUUNIF,1.D-16),0.D0)
        CALL SCAN_X_FX(S,ALPHAM,EPSXM)
        EPSXM=EXFAC*RHOT(IR)**(1.D0/3.D0)*EPSXM
        EPSCM=0.D0
        IF(TCORR) CALL SCAN_C_EPSC(RS,S,T2,ALPHAM,EPSCM)
        EXCDENSM=RHOT(IR)*(EPSXM+EPSCM)
!
        VTAU(IR)=(EXCDENSP-EXCDENSM)/(TAUP-TAUM)*Y0*FOURPI
      ENDDO
      RETURN
      END
