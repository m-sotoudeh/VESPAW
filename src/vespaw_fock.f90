!
!.......................................................................
!**  PLANE-WAVE PAW FOCK EXCHANGE                                     **
!**                                                                   **
!**  STATUS (SEE docs/FOCK_DESIGN.md, ALSO FOR THE ALIGNMENT WITH     **
!**  VASP 6.4.2):                                                     **
!**    MILESTONE 1 (2026-11-23): PLANE-WAVE PART AT THE GAMMA POINT,  **
!**      COLLINEAR SPIN, BARE OR ERFC-SCREENED COULOMB KERNEL.        **
!**    MILESTONE 2 (2026-11-24): PAW. COMPENSATION CHARGES IN EVERY   **
!**      PAIR DENSITY, ONE-CENTER (AE - PS) EXCHANGE AS A DENSITY-    **
!**      MATRIX FUNCTIONAL (DH), G=0 TERM (AUXILIARY FUNCTION).       **
!**      CORE-VALENCE EXCHANGE AND THE SCREENED ONE-CENTER KERNEL     **
!**      ARE STILL MISSING.                                           **
!**                                                                   **
!**  PAW DECOMPOSITION                                                **
!**    E_X = E_X^PW[RHOT_MN + RHOHAT_MN]                               **
!**        + SUM_R ( E_X^1[RHO1_MN] - E_X^1[RHOT1_MN + RHOHAT_MN] )    **
!**    RHOHAT_MN = SUM_{R,L} Q_MN^{RL} GHAT_L(|R-R_R|) Y_L,            **
!**    Q_MN^{RL} = SUM_IJ <PSI_M|P_I>^* <P_J|PSI_N> T_IJ^L,            **
!**    T_IJ^L = GAUNT(I,J,L) INT DR R^(2+L) (PHI_I PHI_J - TPHI_I TPHI_J)**
!**    GHAT_L(R) = N_L R^L (1-(R/RC)^2)^NPOLY FOR R<RC, ZERO OUTSIDE,  **
!**    INT GHAT_L R^(L+2) DR = 1. RC=0.6*RBOX (OR RCOMP, IN BOHR).     **
!**    THE SHAPE IS CONFINED TO THE SPHERE (WHERE THE PARTIAL-WAVE     **
!**    EXPANSION HOLDS) AND SMOOTH ENOUGH FOR THE DENSITY PLANE-WAVE   **
!**    SET (THE ANALOGUE OF VASP AUG_SOFT); THE ONE-CENTER PS TERM USES**
!**    THE SAME SHAPE. A GAUSSIAN SHAPE (FIRST TRY) LEAKED OUT OF THE  **
!**    SPHERE AND MADE THE ENERGY DEPEND ON ITS WIDTH.                 **
!**                                                                   **
!**  CONVENTIONS                                                      **
!**    OCC(N) FROM DYNOCC CONTAINS THE K-POINT WEIGHT AND THE SPIN    **
!**    DEGENERACY (F<=2 FOR NSPIN=1). THE OPERATOR IS                 **
!**      V_X|PSI> = -AEXX SUM_M W_M |PHI_M> (PHI_M^* PSI | V)  + PROJ. **
!**    WITH W_M = S*F_M, S=1/2 FOR NSPIN=1 (COLLINEAR), ELSE S=1.     **
!**    AND THE ENERGY E_X = 1/2 SUM_N F_N <PSI_N|V_X|PSI_N>.          **
!**    DENMAT(I,J)=SUM_N F_N <P_I|PSI_N><PSI_N|P_J> = (TOTAL[,SPIN]), **
!**    E = SUM DH(I,J) DENMAT(J,I) (SAME AS WAVES_HPROJ, ONSITE$ETOT).**
!**                                                                   **
!**  INTERFACE                                                        **
!**    FOCK$SETL4/GETL4('ON'), FOCK$SETR8/GETR8('AEXX','MU','RCOMP', **
!**      'ENERGY')                                                    **
!**    FOCK$SETKMESH(NKPT,XK,WK,NKDIV), FOCK$GETI4('NQFULL')          **
!**    FOCK$INITIALIZE(RBAS,EPW,NR1START,NR1L,NR2,NR3)                **
!**    FOCK$SNAPSHOTRESET(NSPIN,NDIM,NOCCX,NAT,NPRO,R)                **
!**    FOCK$SNAPSHOT(GSETID,ISPIN,NGL,NDIM,NBH,NB,NPRO,TINV,F,PSI,PROJ)**
!**    FOCK$SNAPSHOTGATHER, FOCK$WBUILD(...), FOCK$WBUILDDONE         **
!**    FOCK$ENERGY(EX)                                                **
!**    FOCK$APPLY(GSETID,ISPIN,NGL,NDIM,NBH,NPRO,PSI,PROJ,HPSI,HPROJ) **
!**    FOCK$ONECENTER(LMNXX,NDIMD,NAT,DENMAT,DH,E1)                    **
!**    FOCK$NEWCELL(RBAS), FOCK$FORCE(NAT,TSTRESS,FORCE,STRESS)       **
!**    PSI(NGL,NDIM,NBH), PROJ(NDIM,NBH,NPRO) AS IN WAVES_MODULE.     **
!.......................................................................
MODULE FOCK_MODULE
IMPLICIT NONE
TYPE FOCKSP_TYPE
  LOGICAL(4)             :: TINI=.FALSE.
  INTEGER(4)             :: LNX=0
  INTEGER(4)             :: LMNX=0
  INTEGER(4)             :: LMAXC=0      ! MAX. L OF THE COMPENSATION CHARGES
  INTEGER(4)             :: LMXC=0       ! (LMAXC+1)**2
  INTEGER(4),ALLOCATABLE :: LOX(:)       ! (LNX)
  INTEGER(4),ALLOCATABLE :: LNOF(:)      ! (LMNX) LN OF LMN
  INTEGER(4),ALLOCATABLE :: LMOF(:)      ! (LMNX) LM=L**2+M OF LMN
  REAL(8)   ,ALLOCATABLE :: T(:,:,:)     ! (LMNX,LMNX,LMXC) GAUNT*MOMENT
  REAL(8)   ,ALLOCATABLE :: U(:,:,:,:)   ! (LMNX**4) (IJ|KL)_AE-(IJ|KL)_PS
  REAL(8)   ,ALLOCATABLE :: KCV(:,:)     ! (LMNX,LMNX) SUM_C (IC|CJ) CORE-VALENCE
  REAL(8)                :: RC=0.D0      ! RADIUS OF GHAT
  INTEGER(4)             :: GID=0        ! RADIAL GRID OF THE SHAPES
  INTEGER(4)             :: NR=0
  REAL(8)   ,ALLOCATABLE :: R(:)         ! (NR)
  REAL(8)   ,ALLOCATABLE :: GHAT(:,:)    ! (NR,0:LMAXC) COMPENSATION SHAPES
END TYPE FOCKSP_TYPE
!== PLANE-WAVE SET OF THE PAIR DENSITIES FOR ONE BLOCH VECTOR Q (=K-Q' FOLDED) =
!== K=Q+G ARE THE PHYSICAL WAVE VECTORS (GVEC OF A SET INCLUDES Q)             =
TYPE FOCKQ_TYPE
  CHARACTER(16)          :: ID=' '
  REAL(8)                :: XQ(3)=0.D0   ! Q IN RELATIVE COORDINATES, IN [0,1)
  INTEGER(4)             :: NGL=0
  REAL(8)   ,ALLOCATABLE :: GVEC(:,:)    ! (3,NGL) K=Q+G
  REAL(8)   ,ALLOCATABLE :: VKER(:)      ! (NGL) KERNEL V(|Q+G|), V0 AT K=0
  REAL(8)   ,ALLOCATABLE :: YLM(:,:)     ! (NGL,LMXCX) Y_L(K/|K|)
  REAL(8)   ,ALLOCATABLE :: FG(:,:,:)    ! (NGL,0:LMAXCX,NSPF) SHAPE FT AT |K|
  REAL(8)   ,ALLOCATABLE :: FGD(:,:,:)   ! (NGL,0:LMAXCX,NSPF) D(FG)/D|K|
  COMPLEX(8),ALLOCATABLE :: EIGR(:,:)    ! (NGL,NATF) EXP(-I K R)
END TYPE FOCKQ_TYPE
!== ACE VECTORS OF ONE (K,SPIN) SET ============================================
TYPE FOCKACE_TYPE
  INTEGER(4)             :: N=0          ! #(VECTORS), 0: NOT BUILT
  INTEGER(4)             :: NGL=0
  COMPLEX(8),ALLOCATABLE :: XIG(:,:)     ! (NGL,N) PLANE-WAVE PART
  COMPLEX(8),ALLOCATABLE :: XIP(:,:)     ! (NPROF,N) PROJECTOR PART
END TYPE FOCKACE_TYPE
LOGICAL(4)             :: TON=.FALSE.
LOGICAL(4)             :: TINI=.FALSE.
REAL(8)                :: AEXX=0.25D0   ! FRACTION OF EXACT EXCHANGE
REAL(8)                :: MU=0.D0       ! ERFC SCREENING [1/BOHR], 0: BARE COULOMB
REAL(8)                :: RCOMP=-1.D0   ! RADIUS OF GHAT [BOHR], <0: 0.6*RBOX
LOGICAL(4)             :: TCV=.TRUE.    ! CORE-VALENCE EXCHANGE (ONE-CENTER)
INTEGER(4),PARAMETER   :: NPOLY=6       ! SMOOTHNESS OF GHAT AT RC
CHARACTER(16),PARAMETER:: FOCKID='FOCK' ! PLANE-WAVE SET OF Q=0 (QS(1))
INTEGER(4)             :: NRLF=0        ! #(LOCAL R-POINTS)
INTEGER(4)             :: NR1F=0,NR2F=0,NR3F=0,NR1STARTF=0,NR1LF=0
REAL(8)                :: GWEIGHTF=0.D0 ! CELL VOLUME
REAL(8)                :: RWEIGHTF=0.D0 ! VOLUME PER GRID POINT
REAL(8)                :: V0=0.D0       ! KERNEL AT K=0
REAL(8)                :: ALPHAF=0.D0   ! EWALD PARAMETER OF V0
REAL(8)                :: RBASF(3,3)=0.D0 ! LATTICE VECTORS OF THE SETS
REAL(8)                :: GMAXF=0.D0    ! RADIUS OF THE FOCK G SPHERE
INTEGER(4)             :: LMAXCX=-1     ! MAX LMAXC OVER SPECIES
!== K-POINT MESH ===============================================================
INTEGER(4)             :: NKSETF=0      ! #(K-POINTS OF THE WAVE FUNCTIONS)
REAL(8)   ,ALLOCATABLE :: XKSETF(:,:)   ! (3,NKSETF) RELATIVE
INTEGER(4)             :: NQFULL=0      ! #(POINTS OF THE FULL MESH)
REAL(8)   ,ALLOCATABLE :: XQFULL(:,:)   ! (3,NQFULL) RELATIVE
INTEGER(4)             :: NKDIVF(3)=1
!== THE FULL MESH IS THE UNION OF THE STARS OF THE K-POINTS OF THE WAVE        =
!== FUNCTIONS: Q = S MINVT(IOP) K WITH S=+1, OR S=-1 (TIME REVERSAL).          =
!== WITHOUT SPACE-GROUP SYMMETRY THE ONLY OPERATION IS THE IDENTITY.           =
INTEGER(4)             :: NOPF=0
INTEGER(4),ALLOCATABLE :: IROTF(:,:,:)  ! (3,3,NOPF) X'=IROT X+TAU (RELATIVE)
INTEGER(4),ALLOCATABLE :: MINVTF(:,:,:) ! (3,3,NOPF) K'=MINVT K =(IROT^-1)^T K
REAL(8)   ,ALLOCATABLE :: TAUF(:,:)     ! (3,NOPF)
REAL(8)   ,ALLOCATABLE :: CROTF(:,:,:)  ! (3,3,NOPF) CARTESIAN ROTATION
INTEGER(4),ALLOCATABLE :: MAPATF(:,:)   ! (NAT,NOPF) ATOM I -> MAPATF(I,IOP)
INTEGER(4)             :: NMEMXF=0      ! MAX. #(STAR MEMBERS)
INTEGER(4),ALLOCATABLE :: NMEMF(:)      ! (NKSETF) #(MESH POINTS IN THE STAR)
INTEGER(4),ALLOCATABLE :: MEMQF(:,:)    ! (NMEMXF,NKSETF) INDEX IN XQFULL
INTEGER(4),ALLOCATABLE :: MEMOPF(:,:)   ! (NMEMXF,NKSETF) OPERATION
INTEGER(4),ALLOCATABLE :: MEMSF(:,:)    ! (NMEMXF,NKSETF) +1, -1 (TIME REVERSAL)
LOGICAL(4),ALLOCATABLE :: MEMIDF(:,:)   ! (NMEMXF,NKSETF) OPERATION IS IDENTITY
REAL(8)   ,ALLOCATABLE :: XPOSF(:,:)    ! (3,NATF) RELATIVE ATOMIC POSITIONS
REAL(8)   ,ALLOCATABLE :: UROTF(:,:,:,:)! (LMNXX,LMNXX,NSPF,NOPF) Y_LM ROTATION
!== PAIR-DENSITY SETS (ONE PER Q), AND THE SELECTED ONE ========================
INTEGER(4)             :: NQF=0
TYPE(FOCKQ_TYPE),ALLOCATABLE,TARGET :: QS(:)
INTEGER(4)             :: NGLF=0        ! #(LOCAL G-VECTORS OF THE SELECTED SET)
REAL(8)   ,POINTER     :: VKER(:)=>NULL()
REAL(8)   ,POINTER     :: YLMGF(:,:)=>NULL()
REAL(8)   ,POINTER     :: FGQ(:,:,:)=>NULL()
COMPLEX(8),POINTER     :: EIGRF(:,:)=>NULL()
!== SPECIES AND ATOMS ==========================================================
INTEGER(4)             :: NSPF=0
TYPE(FOCKSP_TYPE),ALLOCATABLE :: SPF(:)
INTEGER(4)             :: NATF=0
INTEGER(4)             :: NPROF=0
INTEGER(4),ALLOCATABLE :: ISPF(:)       ! (NATF) SPECIES
INTEGER(4),ALLOCATABLE :: IPRO0F(:)     ! (NATF) FIRST PROJECTOR INDEX-1
!== SNAPSHOT OF THE OCCUPIED ORBITALS (PERIODIC PARTS, REAL SPACE) OF THE FULL =
!== MESH (TIME-REVERSAL PARTNERS ADDED AS COMPLEX CONJUGATES) ==================
INTEGER(4)             :: NSPINF=0
INTEGER(4)             :: NDIMF=1       ! 1 OR 2 (SPINORS, NONCOLLINEAR)
INTEGER(4)             :: NRLD=0        ! NRLF*NDIMF (SPINOR COMPONENTS STACKED)
INTEGER(4)             :: NPROD=0       ! NPROF*NDIMF
INTEGER(4)             :: NOCCXF=0
INTEGER(4),ALLOCATABLE :: NOCC(:)       ! (NSPINF)
REAL(8)   ,ALLOCATABLE :: WOCC(:,:)     ! (NOCCXF,NSPINF) S*F
INTEGER(4),ALLOCATABLE :: IQOCC(:,:)    ! (NOCCXF,NSPINF) INDEX IN XQFULL
COMPLEX(8),ALLOCATABLE :: PHIR(:,:,:)   ! (NRLD,NOCCXF,NSPINF) COMPONENT-MAJOR
COMPLEX(8),ALLOCATABLE :: CPROJ(:,:,:)  ! (NPROD,NOCCXF,NSPINF) <P|PHI_S>, I+NPROF*(S-1)
REAL(8)                :: EXF=0.D0      ! LAST PLANE-WAVE FOCK ENERGY
REAL(8)                :: EX1F=0.D0     ! LAST ONE-CENTER FOCK ENERGY (INCL. CV)
REAL(8)                :: EXCVF=0.D0    ! LAST CORE-VALENCE FOCK ENERGY
INTEGER(4),PARAMETER   :: NBLK=16       ! BANDS PER BATCHED FFT
!== ADAPTIVELY COMPRESSED EXCHANGE (ACE, LIN, J.CHEM.THEORY COMPUT.12,2242(2016)) =
!== V_X ~ -SUM_K |XI_K><XI_K|, XI = W L^{-H}, W=V_X PSI, -<PSI|W> = L L^H,       =
!== BUILT ONCE PER ENERGY EVALUATION FROM ALL BANDS OF THE SET (EXACT ON THEIR   =
!== SPAN, AS VASP). XI HAS A PLANE-WAVE PART (WAVE G SET) AND A PROJECTOR PART.  =
LOGICAL(4)             :: TACE=.TRUE.   ! USE ACE IN FOCK$APPLY
TYPE(FOCKACE_TYPE),ALLOCATABLE :: ACES(:,:) ! (NKSETF,NSPINF)
CONTAINS
!       ........................................................................
        SUBROUTINE FOCK_SELECTQ(IQ)
!       == MAKES PAIR-DENSITY SET IQ THE CURRENT ONE ==========================
        INTEGER(4),INTENT(IN) :: IQ
        NGLF=QS(IQ)%NGL
        VKER=>QS(IQ)%VKER
        YLMGF=>QS(IQ)%YLM
        FGQ=>QS(IQ)%FG
        EIGRF=>QS(IQ)%EIGR
        CALL PLANEWAVE$SELECT(QS(IQ)%ID)
        END SUBROUTINE FOCK_SELECTQ
!       ........................................................................
        INTEGER(4) FUNCTION FOCK_FINDK(X,N,XLIST)
!       == INDEX OF X IN XLIST MODULO RECIPROCAL LATTICE VECTORS, 0 IF ABSENT ==
        INTEGER(4),INTENT(IN) :: N
        REAL(8)   ,INTENT(IN) :: X(3),XLIST(3,N)
        INTEGER(4)            :: I
        REAL(8)               :: D(3)
        FOCK_FINDK=0
        DO I=1,N
          D=X-XLIST(:,I)
          D=D-ANINT(D)
          IF(MAXVAL(ABS(D)).LT.1.D-8) THEN
            FOCK_FINDK=I
            RETURN
          END IF
        ENDDO
        END FUNCTION FOCK_FINDK
END MODULE FOCK_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$SETL4(ID,VAL)
      USE FOCK_MODULE, ONLY : TON,TCV,TACE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
      IF(ID.EQ.'ON') THEN
        TON=VAL
      ELSE IF(ID.EQ.'CORE') THEN
        TCV=VAL
      ELSE IF(ID.EQ.'ACE') THEN
        TACE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('FOCK$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$GETL4(ID,VAL)
      USE FOCK_MODULE, ONLY : TON
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
      IF(ID.EQ.'ON') THEN
        VAL=TON
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('FOCK$GETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$SETR8(ID,VAL)
      USE FOCK_MODULE, ONLY : AEXX,MU,RCOMP,TINI
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
      IF(TINI) THEN
        CALL ERROR$MSG('PARAMETERS CANNOT BE CHANGED AFTER INITIALIZATION')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('FOCK$SETR8')
      END IF
      IF(ID.EQ.'AEXX') THEN
        AEXX=VAL
      ELSE IF(ID.EQ.'MU') THEN
        IF(VAL.LT.0.D0) THEN
          CALL ERROR$MSG('SCREENING PARAMETER MU MUST NOT BE NEGATIVE')
          CALL ERROR$R8VAL('MU',VAL)
          CALL ERROR$STOP('FOCK$SETR8')
        END IF
        MU=VAL
      ELSE IF(ID.EQ.'RCOMP') THEN
        RCOMP=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('FOCK$SETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$GETR8(ID,VAL)
      USE FOCK_MODULE, ONLY : AEXX,MU,RCOMP,EXF,EX1F,EXCVF
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VAL
      IF(ID.EQ.'AEXX') THEN
        VAL=AEXX
      ELSE IF(ID.EQ.'MU') THEN
        VAL=MU
      ELSE IF(ID.EQ.'RCOMP') THEN
        VAL=RCOMP
      ELSE IF(ID.EQ.'ENERGY') THEN
        VAL=EXF+EX1F
      ELSE IF(ID.EQ.'EPW') THEN
        VAL=EXF
      ELSE IF(ID.EQ.'E1C') THEN
        VAL=EX1F
      ELSE IF(ID.EQ.'ECV') THEN
        VAL=EXCVF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('FOCK$GETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$GETI4(ID,VAL)
      USE FOCK_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(OUT) :: VAL
      IF(ID.EQ.'NQFULL') THEN    ! #(POINTS OF THE FULL K MESH)
        VAL=NQFULL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('FOCK$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$SETKMESH(NKPT,XK,WK,NKDIV)
!     **************************************************************************
!     **  K-POINTS OF THE WAVE FUNCTIONS (RELATIVE COORDINATES), THEIR        **
!     **  WEIGHTS AND THE MESH DIVISIONS. THE FULL MESH IS THE UNION OF THE   **
!     **  STARS Q = S MINVT K OF THE LISTED POINTS UNDER THE SPACE-GROUP      **
!     **  OPERATIONS (ONLY THE IDENTITY WITHOUT SYMMETRY) AND TIME REVERSAL   **
!     **  (S=-1, U_{-K} = CONJG(U_K)). A POINT OF THE LIST IS THE FIRST MEMBER **
!     **  OF ITS OWN STAR (IDENTITY). THE STARS MUST NOT OVERLAP, THEIR SIZES **
!     **  MUST MATCH THE WEIGHTS AND THEY MUST FILL THE MESH.                 **
!     **  CALL BEFORE FOCK$INITIALIZE.                                        **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NKPT
      REAL(8)   ,INTENT(IN) :: XK(3,NKPT)
      REAL(8)   ,INTENT(IN) :: WK(NKPT)
      INTEGER(4),INTENT(IN) :: NKDIV(3)
      INTEGER(4),PARAMETER  :: NOPXX=1536
      INTEGER(4),ALLOCATABLE:: IROTX(:,:,:),MAPX(:,:),IKOFQ(:)
      REAL(8)   ,ALLOCATABLE:: TAUX(:,:),CROTX(:,:,:)
      INTEGER(4)            :: NAT,NOP,IOP,IOPID,IK,IS,JQ,M,IAT,I
      REAL(8)               :: A(3,3),AINV(3,3),XP(3),WSUM
      LOGICAL(4)            :: TSYM
!     **************************************************************************
!
!     == SPACE-GROUP OPERATIONS (IDENTITY ONLY WITHOUT SYMMETRY) ==============
      CALL ATOMLIST$NATOM(NAT)
      CALL SYMMETRY$GETL4('ON',TSYM)
      ALLOCATE(IROTX(3,3,NOPXX),TAUX(3,NOPXX),CROTX(3,3,NOPXX),MAPX(NAT,NOPXX))
      NOP=0
      IF(TSYM) CALL SYMMETRY$GETOPS(NOPXX,NAT,NOP,IROTX,TAUX,CROTX,MAPX)
      IF(NOP.EQ.0) THEN
        NOP=1
        IROTX(:,:,1)=0
        CROTX(:,:,1)=0.D0
        DO I=1,3
          IROTX(I,I,1)=1
          CROTX(I,I,1)=1.D0
        ENDDO
        TAUX(:,1)=0.D0
        DO IAT=1,NAT
          MAPX(IAT,1)=IAT
        ENDDO
      END IF
      IF(ALLOCATED(IROTF)) DEALLOCATE(IROTF,MINVTF,TAUF,CROTF,MAPATF)
      NOPF=NOP
      ALLOCATE(IROTF(3,3,NOPF),MINVTF(3,3,NOPF),TAUF(3,NOPF),CROTF(3,3,NOPF))
      ALLOCATE(MAPATF(NAT,NOPF))
      IROTF=IROTX(:,:,1:NOPF)
      TAUF=TAUX(:,1:NOPF)
      CROTF=CROTX(:,:,1:NOPF)
      MAPATF=MAPX(:,1:NOPF)
      DEALLOCATE(IROTX,TAUX,CROTX,MAPX)
      IOPID=0
      DO IOP=1,NOPF
        A=REAL(IROTF(:,:,IOP),KIND=8)
        CALL LIB$INVERTR8(3,A,AINV)
        MINVTF(:,:,IOP)=NINT(TRANSPOSE(AINV))
        IF(IOPID.EQ.0) THEN
          A=REAL(IROTF(:,:,IOP),KIND=8)
          DO I=1,3
            A(I,I)=A(I,I)-1.D0
          ENDDO
          XP=TAUF(:,IOP)-ANINT(TAUF(:,IOP))
          IF(MAXVAL(ABS(A)).LT.1.D-8.AND.MAXVAL(ABS(XP)).LT.1.D-8) IOPID=IOP
        END IF
      ENDDO
      IF(IOPID.EQ.0) THEN
        CALL ERROR$MSG('IDENTITY NOT AMONG THE SYMMETRY OPERATIONS')
        CALL ERROR$STOP('FOCK$SETKMESH')
      END IF
!
!     == STARS ================================================================
      NKSETF=NKPT
      NMEMXF=2*NOPF
      IF(ALLOCATED(XKSETF)) DEALLOCATE(XKSETF)
      IF(ALLOCATED(XQFULL)) DEALLOCATE(XQFULL)
      IF(ALLOCATED(NMEMF)) DEALLOCATE(NMEMF,MEMQF,MEMOPF,MEMSF,MEMIDF)
      ALLOCATE(XKSETF(3,NKSETF))
      ALLOCATE(XQFULL(3,NKSETF*NMEMXF))
      ALLOCATE(IKOFQ(NKSETF*NMEMXF))
      ALLOCATE(NMEMF(NKSETF))
      ALLOCATE(MEMQF(NMEMXF,NKSETF))
      ALLOCATE(MEMOPF(NMEMXF,NKSETF))
      ALLOCATE(MEMSF(NMEMXF,NKSETF))
      ALLOCATE(MEMIDF(NMEMXF,NKSETF))
      XKSETF=XK
      NKDIVF=NKDIV
      NQFULL=0
      DO IK=1,NKSETF
        IF(FOCK_FINDK(XK(:,IK),NQFULL,XQFULL(:,1:NQFULL)).NE.0) THEN
          CALL ERROR$MSG('K-POINT LISTED TWICE')
          CALL ERROR$I4VAL('IK',IK)
          CALL ERROR$STOP('FOCK$SETKMESH')
        END IF
        NQFULL=NQFULL+1
        XQFULL(:,NQFULL)=XK(:,IK)
        IKOFQ(NQFULL)=IK
        NMEMF(IK)=1
        MEMQF(1,IK)=NQFULL
        MEMOPF(1,IK)=IOPID
        MEMSF(1,IK)=1
        MEMIDF(1,IK)=.TRUE.
      ENDDO
!     == ROTATIONS FIRST (THEY MUST NOT REACH ANOTHER LISTED STAR), THEN     ==
!     == TIME REVERSAL (ADDS ONLY POINTS NOT REACHED OTHERWISE)              ==
      DO IS=1,-1,-2
        DO IK=1,NKSETF
          DO IOP=1,NOPF
            XP=REAL(IS,KIND=8)*MATMUL(REAL(MINVTF(:,:,IOP),KIND=8),XK(:,IK))
            JQ=FOCK_FINDK(XP,NQFULL,XQFULL(:,1:NQFULL))
            IF(JQ.EQ.0) THEN
              NQFULL=NQFULL+1
              XQFULL(:,NQFULL)=XP
              IKOFQ(NQFULL)=IK
              M=NMEMF(IK)+1
              NMEMF(IK)=M
              MEMQF(M,IK)=NQFULL
              MEMOPF(M,IK)=IOP
              MEMSF(M,IK)=IS
              MEMIDF(M,IK)=(IOP.EQ.IOPID)
            ELSE IF(IS.EQ.1.AND.IKOFQ(JQ).NE.IK) THEN
              CALL ERROR$MSG('K-POINT LIST IS NOT IRREDUCIBLE:')
              CALL ERROR$MSG('TWO LISTED POINTS ARE RELATED BY SYMMETRY')
              CALL ERROR$I4VAL('IK',IK)
              CALL ERROR$I4VAL('OTHER',IKOFQ(JQ))
              CALL ERROR$STOP('FOCK$SETKMESH')
            END IF
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(IKOFQ)
      IF(NQFULL.NE.NKDIV(1)*NKDIV(2)*NKDIV(3)) THEN
        CALL ERROR$MSG('THE STARS OF THE K-POINTS DO NOT FILL THE REGULAR MESH')
        CALL ERROR$I4VAL('#(POINTS OF THE STARS)',NQFULL)
        CALL ERROR$I4VAL('NKDIV(1)',NKDIV(1))
        CALL ERROR$I4VAL('NKDIV(2)',NKDIV(2))
        CALL ERROR$I4VAL('NKDIV(3)',NKDIV(3))
        CALL ERROR$STOP('FOCK$SETKMESH')
      END IF
!     == STAR SIZES MUST MATCH THE K-POINT WEIGHTS ============================
      WSUM=SUM(WK)
      DO IK=1,NKSETF
        IF(ABS(WK(IK)/WSUM-REAL(NMEMF(IK),KIND=8)/REAL(NQFULL,KIND=8)).GT.1.D-6) THEN
          CALL ERROR$MSG('K-POINT WEIGHT DIFFERS FROM THE SIZE OF ITS STAR')
          CALL ERROR$I4VAL('IK',IK)
          CALL ERROR$R8VAL('WEIGHT',WK(IK)/WSUM)
          CALL ERROR$I4VAL('#(STAR)',NMEMF(IK))
          CALL ERROR$I4VAL('#(MESH)',NQFULL)
          CALL ERROR$STOP('FOCK$SETKMESH')
        END IF
      ENDDO
      IF(NOPF.GT.1) WRITE(*,FMT='("FOCK: ",I4," SYMMETRY OPERATIONS, ",I6' &
     &   //'," IRREDUCIBLE K-POINTS, ",I6," POINTS IN THE FULL MESH")')NOPF,NKSETF,NQFULL
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$INITIALIZE(RBAS,EPW,NR1START,NR1L,NR2,NR3)
!     **************************************************************************
!     **  DEFINES ONE PLANE-WAVE SET PER DISTINCT BLOCH VECTOR Q=K-Q' (FOLDED **
!     **  INTO [0,1)) OF THE PAIR DENSITIES, ON THE REAL-SPACE GRID SLICE OF  **
!     **  THE WAVE FUNCTIONS (COMMUNICATOR 'K'), WITH THE DENSITY CUTOFF EPW. **
!     **  FOR EPW>=4*EPWPSI (CDUAL>=4) A PRODUCT OF TWO ORBITALS IS EXACT     **
!     **  (VASP PRECFOCK=ACCURATE). SET 1 IS Q=0 (ID 'FOCK').                 **
!     **  KERNEL V(|Q+G|) = 4PI/K^2 [1-EXP(-K^2/4MU^2)]; AT K=0 THE VALUE V0   **
!     **  OF THE BORN-VON KARMAN SUPERCELL OF THE MESH (FOCK_V0), SO THAT A    **
!     **  N1xN2xN3 MESH IS EXACTLY EQUIVALENT TO THE GAMMA-ONLY SUPERCELL.    **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: EPW
      INTEGER(4),INTENT(IN) :: NR1START,NR1L,NR2,NR3
      REAL(8)   ,ALLOCATABLE:: XQ(:,:),G2(:)
      COMPLEX(8),ALLOCATABLE:: FR(:),FG(:)
      REAL(8)               :: RBASS(3,3),D(3),PI,FOURPI,SVAR,GMAX,ALPHA,SUMG
      INTEGER(4)            :: IK,IQ,IG,NGAMMA,I,NQX
      CHARACTER(32)         :: OLDID
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(TINI) RETURN
                              CALL TRACE$PUSH('FOCK$INITIALIZE')
      IF(NKSETF.EQ.0) THEN
        CALL ERROR$MSG('FOCK$SETKMESH MUST BE CALLED BEFORE FOCK$INITIALIZE')
        CALL ERROR$STOP('FOCK$INITIALIZE')
      END IF
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      CALL PLANEWAVE$GETCH('ID',OLDID)
      NR1STARTF=NR1START
      NR1LF=NR1L
      NR2F=NR2
      NR3F=NR3
!
!     == DISTINCT Q=K-Q' MODULO G (Q=0 FIRST) ==================================
      NQX=NKSETF*NQFULL+1
      ALLOCATE(XQ(3,NQX))
      NQF=1
      XQ(:,1)=0.D0
      DO IK=1,NKSETF
        DO IQ=1,NQFULL
          D=XKSETF(:,IK)-XQFULL(:,IQ)
          D=MODULO(D,1.D0)
          DO I=1,3
            IF(ABS(D(I)-1.D0).LT.1.D-8) D(I)=0.D0
          ENDDO
          IF(FOCK_FINDK(D,NQF,XQ(:,1:NQF)).EQ.0) THEN
            NQF=NQF+1
            XQ(:,NQF)=D
          END IF
        ENDDO
      ENDDO
      ALLOCATE(QS(NQF))
      GMAX=SQRT(2.D0*EPW)
      GMAXF=GMAX
!
!     == G=0 ELEMENT: EWALD FORM FOR THE BORN-VON KARMAN SUPERCELL ============
      DO I=1,3
        RBASS(:,I)=RBAS(:,I)*REAL(NKDIVF(I),KIND=8)
      ENDDO
      ALPHA=25.D0/GMAX**2
      CALL FOCK_V0(RBASS,ALPHA,V0)
      CALL FOCK_V0(RBASS,0.5D0*ALPHA,SUMG)
      ALPHAF=ALPHA
      RBASF=RBAS
      WRITE(*,FMT='("FOCK: V(K=0)=",F16.8," (ALPHA=",F8.4,")",F16.8," (ALPHA/2)' &
     &    //' MESH",3I3," #Q-SETS=",I4)')V0,ALPHA,SUMG,NKDIVF,NQF
!
      DO IQ=1,NQF
        QS(IQ)%XQ=XQ(:,IQ)
        IF(IQ.EQ.1) THEN
          QS(IQ)%ID=FOCKID
        ELSE
          WRITE(QS(IQ)%ID,FMT='("FOCKQ",I5.5)')IQ
        END IF
        CALL PLANEWAVE$INITIALIZE(QS(IQ)%ID,'K',RBAS,QS(IQ)%XQ,.FALSE.,EPW &
     &                           ,NR1START,NR1L,NR2,NR3)
        CALL PLANEWAVE$SELECT(QS(IQ)%ID)
        CALL PLANEWAVE$GETI4('NGL',QS(IQ)%NGL)
        IF(IQ.EQ.1) THEN
          CALL PLANEWAVE$GETI4('NRL',NRLF)
          CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHTF)
          CALL PLANEWAVE$GETR8('RWEIGHT',RWEIGHTF)
          CALL PLANEWAVE$GETI4('NR1',NR1F)
        END IF
        ALLOCATE(QS(IQ)%GVEC(3,QS(IQ)%NGL))
        CALL PLANEWAVE$GETR8A('GVEC',3*QS(IQ)%NGL,QS(IQ)%GVEC)
        ALLOCATE(G2(QS(IQ)%NGL))
        DO IG=1,QS(IQ)%NGL
          G2(IG)=SUM(QS(IQ)%GVEC(:,IG)**2)
        ENDDO
        ALLOCATE(QS(IQ)%VKER(QS(IQ)%NGL))
        NGAMMA=0
        DO IG=1,QS(IQ)%NGL
          IF(G2(IG).LT.1.D-10) THEN
            NGAMMA=IG
            QS(IQ)%VKER(IG)=V0
          ELSE IF(MU.GT.0.D0) THEN
            QS(IQ)%VKER(IG)=FOURPI/G2(IG)*(1.D0-EXP(-G2(IG)/(4.D0*MU**2)))
          ELSE
            QS(IQ)%VKER(IG)=FOURPI/G2(IG)
          END IF
        ENDDO
        IF(NGAMMA.GT.0.AND.IQ.GT.1) THEN
          CALL ERROR$MSG('K=0 FOUND IN A SET WITH Q/=0 (GVEC WITHOUT Q?)')
          CALL ERROR$STOP('FOCK$INITIALIZE')
        END IF
        DEALLOCATE(G2)
      ENDDO
      DEALLOCATE(XQ)
!
!     == SELF TEST OF THE FFT NORMALIZATION: A CONSTANT 1 MUST MAP ONTO A  ==
!     == UNIT G=0 COEFFICIENT (F(R)=SUM_G F(G) EXP(IGR)), WHICH THE        ==
!     == EXCHANGE ENERGY FORMULA OMEGA*SUM_G V(G)|RHO(G)|^2 ASSUMES.       ==
      CALL PLANEWAVE$SELECT(FOCKID)
      ALLOCATE(FR(NRLF))
      ALLOCATE(FG(QS(1)%NGL))
      FR(:)=(1.D0,0.D0)
      CALL PLANEWAVE$FFT('RTOG',1,QS(1)%NGL,FG,NRLF,FR)
      SVAR=0.D0
      DO IG=1,QS(1)%NGL
        IF(SUM(QS(1)%GVEC(:,IG)**2).LT.1.D-10) SVAR=ABS(FG(IG)-(1.D0,0.D0))
      ENDDO
      CALL MPE$COMBINE('K','+',SVAR)
      IF(SVAR.GT.1.D-10) THEN
        CALL ERROR$MSG('UNEXPECTED FFT NORMALIZATION OF THE FOCK SET')
        CALL ERROR$R8VAL('|F(G=0)-1|',SVAR)
        CALL ERROR$STOP('FOCK$INITIALIZE')
      END IF
      DEALLOCATE(FR)
      DEALLOCATE(FG)
      IF(LEN_TRIM(OLDID).GT.0) CALL PLANEWAVE$SELECT(TRIM(OLDID))
      TINI=.TRUE.
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$REPORT(NFIL)
      USE FOCK_MODULE, ONLY : TON,AEXX,MU,RCOMP,V0,QS,NQF,NKDIVF,TACE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      IF(.NOT.TON) RETURN
      WRITE(NFIL,FMT='("FOCK EXCHANGE (PLANE-WAVE PAW)")')
      WRITE(NFIL,FMT='(T5,"K MESH",T30,3I4)')NKDIVF
      WRITE(NFIL,FMT='(T5,"#(PAIR-DENSITY SETS Q)",T30,I10)')NQF
      WRITE(NFIL,FMT='(T5,"ACE",T30,L10)')TACE
      WRITE(NFIL,FMT='(T5,"AEXX",T30,F10.5)')AEXX
      WRITE(NFIL,FMT='(T5,"MU [1/BOHR]",T30,F10.5)')MU
      WRITE(NFIL,FMT='(T5,"RCOMP [BOHR]",T30,F10.5)')RCOMP
      WRITE(NFIL,FMT='(T5,"V(K=0) [BOHR^2]",T30,F14.5)')V0
      WRITE(NFIL,FMT='(T5,"#(LOCAL G) FOCK SET Q=0",T30,I10)')QS(1)%NGL
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_SPECIES(ISP)
!     **************************************************************************
!     **  PAW DATA OF ONE SPECIES: MOMENTS T, ONE-CENTER TENSOR U AND THE     **
!     **  FOURIER TRANSFORM OF THE COMPENSATION SHAPE ON THE FOCK SET.        **
!     **                                                                      **
!     **  U(I,J,K,L) = SUM_L GAUNT(I,J,L) GAUNT(K,L,L)                        **
!     **               * [ F_AE^L(AB,CD) - F_PS^L(AB,CD) ]                    **
!     **  F_X^L(AB,CD) = INT R^2 RHO^X_CD(R) V_L[RHO^X_AB](R) DR,              **
!     **  RHO^AE_AB = PHI_A PHI_B, RHO^PS_AB = TPHI_A TPHI_B + Q^L_AB GHAT_L. **
!     **  THE RADIAL DENSITIES ARE CUT AT RINT, BEYOND WHICH AE AND PS PARTIAL**
!     **  WAVES AGREE AND GHAT VANISHES; THE PARTS OUTSIDE CANCEL EXACTLY.    **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP
      INTEGER(4)            :: GID,NR,LNX,LMNX,LMAXC,LMXC
      INTEGER(4)            :: LN,LN1,LN2,LN3,LN4,M,LMN,L,LM,IR,IG
      INTEGER(4)            :: LMN1,LMN2,LMN3,LMN4,LM1,LM2,LM3,LM4,L1,L2,L3,L4
      INTEGER(4)            :: NPAIR,IP1,IP2
      INTEGER(4),ALLOCATABLE:: IPAIR(:,:)
      REAL(8)   ,ALLOCATABLE:: R(:),AEPHI(:,:),PSPHI(:,:),GHAT(:,:),WORK(:)
      REAL(8)   ,ALLOCATABLE:: QL(:,:,:)        ! (LNX,LNX,0:LMAXC)
      REAL(8)   ,ALLOCATABLE:: RHOAE(:,:),RHOPS(:,:),VAE(:,:),VPS(:,:)
      REAL(8)   ,ALLOCATABLE:: F(:,:,:)         ! (NPAIR,NPAIR,0:LMAXC)
      REAL(8)   ,ALLOCATABLE:: GAUNT(:,:,:)     ! (LMNX,LMNX,LMXC)
      REAL(8)   ,ALLOCATABLE:: G2(:)
      REAL(8)               :: RBOX,RINT,PI,SVAR,NORM,GAMMAL,CG,G
      REAL(8)               :: SEXACT,STRUNC,SINF,SINFS
      CHARACTER(32)         :: CENV
      INTEGER(4)            :: NC,NB,IC,LC,MC,LMC,NCV
      INTEGER(4),ALLOCATABLE:: LB(:)
      REAL(8)   ,ALLOCATABLE:: AEPSI(:,:),FB(:),FCV(:,:,:,:),RHOC(:,:)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL SETUP$GETI4('NR',NR)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(SPF(ISP)%LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,SPF(ISP)%LOX)
      ALLOCATE(AEPHI(NR,LNX),PSPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI)
      CALL SETUP$GETR8('RBOX',RBOX)
      CALL SETUP$GETI4('NC',NC)
      CALL SETUP$GETI4('NB',NB)
      ALLOCATE(LB(NB),FB(NB),AEPSI(NR,NB))
      CALL SETUP$GETI4A('LB',NB,LB)
      CALL SETUP$GETR8A('FB',NB,FB)
      CALL SETUP$GETR8A('AEPSI',NR*NB,AEPSI)
      CALL SETUP$UNSELECT()
      ALLOCATE(R(NR),WORK(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     == INDEXING ==============================================================
      LMNX=SUM(2*SPF(ISP)%LOX+1)
      LMAXC=2*MAXVAL(SPF(ISP)%LOX)
      LMXC=(LMAXC+1)**2
      SPF(ISP)%LNX=LNX
      SPF(ISP)%LMNX=LMNX
      SPF(ISP)%LMAXC=LMAXC
      SPF(ISP)%LMXC=LMXC
      ALLOCATE(SPF(ISP)%LNOF(LMNX),SPF(ISP)%LMOF(LMNX))
      LMN=0
      DO LN=1,LNX
        L=SPF(ISP)%LOX(LN)
        DO M=1,2*L+1
          LMN=LMN+1
          SPF(ISP)%LNOF(LMN)=LN
          SPF(ISP)%LMOF(LMN)=L**2+M
        ENDDO
      ENDDO
!
!     == RADIUS OF THE COMPENSATION CHARGES ====================================
!     == DEFAULT 0.6*RBOX (ABOUT THE AUGMENTATION RADIUS RCL OF THE SETUPS):  ==
!     == SMALL ENOUGH TO STAY INSIDE THE REGION WHERE THE PARTIAL-WAVE        ==
!     == EXPANSION HOLDS AND TO AVOID OVERLAP WITH NEIGHBORS (FROZEN-ORBITAL  ==
!     == SI8 TEST: E_X DRIFTS 0.44 MH BETWEEN RC=1.4 AND 2.0 WITH NPRO=2 2 1, ==
!     == 3.1 MH WITH NPRO=1 1 1, I.E. THE DRIFT IS THE COMPLETENESS ERROR).    ==
      SPF(ISP)%RC=0.6D0*RBOX
      IF(RCOMP.GT.0.D0) SPF(ISP)%RC=RCOMP
!
!     == RADIUS BEYOND WHICH ALL RADIAL DENSITIES ARE CUT. THE AE-PS DIFFERENCE==
!     == OF S AND P PARTIAL WAVES DECAYS SMOOTHLY (CORE-ORTHOGONALITY TAILS),  ==
!     == SO THE CUT MUST LIE WELL OUTSIDE RBOX AND MUST NOT DEPEND ON RC.      ==
      RINT=MIN(2.D0*RBOX,0.9D0*R(NR))
      CALL GET_ENVIRONMENT_VARIABLE('VESPAW_FOCK_RINT',CENV)
      IF(LEN_TRIM(CENV).GT.0) READ(CENV,*)RINT
      IF(SPF(ISP)%RC.GT.RINT) THEN
        CALL ERROR$MSG('COMPENSATION RADIUS LARGER THAN THE INTEGRATION RADIUS')
        CALL ERROR$R8VAL('RC',SPF(ISP)%RC)
        CALL ERROR$R8VAL('RINT',RINT)
        CALL ERROR$STOP('FOCK_SPECIES')
      END IF
      DO IR=1,NR
        IF(R(IR).GT.RINT) THEN
          AEPHI(IR,:)=0.D0
          PSPHI(IR,:)=0.D0
        END IF
      ENDDO
!
!     == COMPENSATION SHAPES GHAT_L(R)=N_L R^L (1-(R/RC)^2)^NPOLY, R<RC ========
      ALLOCATE(GHAT(NR,0:LMAXC))
      DO L=0,LMAXC
        DO IR=1,NR
          IF(R(IR).LT.SPF(ISP)%RC) THEN
            GHAT(IR,L)=R(IR)**L*(1.D0-(R(IR)/SPF(ISP)%RC)**2)**NPOLY
          ELSE
            GHAT(IR,L)=0.D0
          END IF
        ENDDO
        WORK(:)=GHAT(:,L)*R(:)**(L+2)
        CALL RADIAL$INTEGRAL(GID,NR,WORK,SVAR)
        GHAT(:,L)=GHAT(:,L)/SVAR
      ENDDO
!
!     == MOMENTS Q^L_AB ========================================================
      ALLOCATE(QL(LNX,LNX,0:LMAXC))
      QL=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          DO L=0,LMAXC
            WORK(:)=(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2)) &
     &             *R(:)**(L+2)
            CALL RADIAL$INTEGRAL(GID,NR,WORK,QL(LN1,LN2,L))
          ENDDO
        ENDDO
      ENDDO
!
!     == GAUNT COEFFICIENTS AND T ==============================================
      ALLOCATE(GAUNT(LMNX,LMNX,LMXC))
      ALLOCATE(SPF(ISP)%T(LMNX,LMNX,LMXC))
      DO LMN1=1,LMNX
        LM1=SPF(ISP)%LMOF(LMN1)
        LN1=SPF(ISP)%LNOF(LMN1)
        DO LMN2=1,LMNX
          LM2=SPF(ISP)%LMOF(LMN2)
          LN2=SPF(ISP)%LNOF(LMN2)
          DO LM=1,LMXC
            CALL SPHERICAL$GAUNT(LM1,LM2,LM,CG)
            GAUNT(LMN1,LMN2,LM)=CG
            L=INT(SQRT(REAL(LM-1,KIND=8)+1.D-6))
            SPF(ISP)%T(LMN1,LMN2,LM)=CG*QL(LN1,LN2,L)
          ENDDO
        ENDDO
      ENDDO
!
!     == RADIAL SLATER INTEGRALS F^L(AB,CD) ====================================
      NPAIR=LNX*LNX
      ALLOCATE(IPAIR(LNX,LNX))
      IP1=0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IP1=IP1+1
          IPAIR(LN1,LN2)=IP1
        ENDDO
      ENDDO
      ALLOCATE(RHOAE(NR,NPAIR),RHOPS(NR,NPAIR),VAE(NR,NPAIR),VPS(NR,NPAIR))
      ALLOCATE(F(NPAIR,NPAIR,0:LMAXC))
      F=0.D0
      DO L=0,LMAXC
        DO LN1=1,LNX
          DO LN2=1,LNX
            IP1=IPAIR(LN1,LN2)
            RHOAE(:,IP1)=AEPHI(:,LN1)*AEPHI(:,LN2)
            RHOPS(:,IP1)=PSPHI(:,LN1)*PSPHI(:,LN2)+QL(LN1,LN2,L)*GHAT(:,L)
            CALL RADIAL$POISSON(GID,NR,L,RHOAE(:,IP1),VAE(:,IP1))
            CALL RADIAL$POISSON(GID,NR,L,RHOPS(:,IP1),VPS(:,IP1))
          ENDDO
        ENDDO
        DO IP1=1,NPAIR
          DO IP2=1,NPAIR
            WORK(:)=R(:)**2*(RHOAE(:,IP2)*VAE(:,IP1)-RHOPS(:,IP2)*VPS(:,IP1))
            CALL RADIAL$INTEGRAL(GID,NR,WORK,F(IP1,IP2,L))
          ENDDO
        ENDDO
!       == ERFC KERNEL: SUBTRACT THE LONG-RANGE ERF PART (Q SPACE) =============
        IF(MU.GT.0.D0) CALL FOCK_ERFPART(GID,NR,R,L,NPAIR,RHOAE,RHOPS,MU,F(:,:,L))
      ENDDO
!
!     == TRUNCATION OF THE COMPENSATION SELF-INTERACTION ON THE GRID: THE     ==
!     == PLANE-WAVE SET CONTAINS (GHAT|GHAT) ONLY FOR |G|<GMAXF, WHILE THE    ==
!     == RADIAL PS TERM ABOVE HAS THE EXACT VALUE. THE PS PAIR DENSITIES ARE  ==
!     == BAND-LIMITED (CDUAL>=4), SO THEIR CROSS TERMS WITH GHAT ARE EXACT ON ==
!     == THE GRID AND ONLY THE ON-SITE GHAT-GHAT TERM MUST BE MATCHED:        ==
!     ==   F_PS(AB,CD) += Q_AB Q_CD [ (GHAT|GHAT)_TRUNC - (GHAT|GHAT)_EXACT ]  ==
!     ==   (GHAT|GHAT)_TRUNC = 1/(2PI^2) INT_0^GMAXF FG(Q)^2 DQ               ==
      DO L=0,LMAXC
        WORK(:)=GHAT(:,L)
        CALL RADIAL$POISSON(GID,NR,L,WORK,VAE(:,1))
        WORK(:)=R(:)**2*GHAT(:,L)*VAE(:,1)
        CALL RADIAL$INTEGRAL(GID,NR,WORK,SEXACT)
        CALL FOCK_SELFENERGY(GID,NR,R,L,GHAT(:,L),GMAXF,0.D0,STRUNC)
        CALL FOCK_SELFENERGY(GID,NR,R,L,GHAT(:,L),20.D0*GMAXF,0.D0,SINF)
        WRITE(*,FMT='("FOCK: GHAT SELF-ENERGY L=",I2," RADIAL=",F14.8' &
     &        //'," Q-SPACE(INF)=",F14.8," TRUNCATED=",F14.8)')L,SEXACT,SINF,STRUNC
!       == ERFC KERNEL: THE GRID HOLDS THE SCREENED KERNEL, THE RADIAL PS TERM ==
!       == ABOVE HAS THE ERF PART REMOVED EXACTLY (FOCK_ERFPART)             ==
        IF(MU.GT.0.D0) THEN
          CALL FOCK_SELFENERGY(GID,NR,R,L,GHAT(:,L),20.D0*GMAXF,MU,SINFS)
          SEXACT=SEXACT-(SINF-SINFS)
          CALL FOCK_SELFENERGY(GID,NR,R,L,GHAT(:,L),GMAXF,MU,STRUNC)
        END IF
        DO LN1=1,LNX
          DO LN2=1,LNX
            IP1=IPAIR(LN1,LN2)
            DO LN3=1,LNX
              DO LN4=1,LNX
                IP2=IPAIR(LN3,LN4)
                F(IP1,IP2,L)=F(IP1,IP2,L) &
     &                -QL(LN1,LN2,L)*QL(LN3,LN4,L)*(STRUNC-SEXACT)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     == ONE-CENTER TENSOR =====================================================
      ALLOCATE(SPF(ISP)%U(LMNX,LMNX,LMNX,LMNX))
      SPF(ISP)%U=0.D0
      DO LMN1=1,LMNX
        LN1=SPF(ISP)%LNOF(LMN1)
        DO LMN2=1,LMNX
          LN2=SPF(ISP)%LNOF(LMN2)
          IP1=IPAIR(LN1,LN2)
          DO LMN3=1,LMNX
            LN3=SPF(ISP)%LNOF(LMN3)
            DO LMN4=1,LMNX
              LN4=SPF(ISP)%LNOF(LMN4)
              IP2=IPAIR(LN3,LN4)
              SVAR=0.D0
              DO LM=1,LMXC
                IF(GAUNT(LMN1,LMN2,LM).EQ.0.D0) CYCLE
                IF(GAUNT(LMN3,LMN4,LM).EQ.0.D0) CYCLE
                L=INT(SQRT(REAL(LM-1,KIND=8)+1.D-6))
                SVAR=SVAR+GAUNT(LMN1,LMN2,LM)*GAUNT(LMN3,LMN4,LM)*F(IP1,IP2,L)
              ENDDO
              SPF(ISP)%U(LMN1,LMN2,LMN3,LMN4)=SVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     == CORE-VALENCE EXCHANGE (ALL-ELECTRON, ONE-CENTER, FROZEN CORE) =========
!     ==   KCV(I,J) = SUM_{C,MC} (PHI_I PHI_C | PHI_C PHI_J)                  ==
!     == EACH CORE SPIN ORBITAL HOLDS ONE ELECTRON, SO PER SPIN               ==
!     ==   E_CV^S = -AEXX SUM_IJ D^S_JI KCV(I,J), E_CV = -AEXX TR(D_TOT KCV). ==
!     == THERE IS NO PSEUDO COUNTERPART (NO CORE IN THE PS WAVE FUNCTIONS).   ==
      ALLOCATE(SPF(ISP)%KCV(LMNX,LMNX))
      SPF(ISP)%KCV=0.D0
      NCV=0
      DO IC=1,NC
        LC=LB(IC)
        WORK(:)=(AEPSI(:,IC)*R(:))**2
        CALL RADIAL$INTEGRAL(GID,NR,WORK,NORM)
        NCV=NCV+1
        WRITE(*,FMT='("FOCK: CORE STATE",I3," L=",I2," F=",F8.4' &
     &             //'," NORM=",F12.8)')IC,LC,FB(IC),NORM
!       == AEPSI IS THE LARGE COMPONENT OF THE SCALAR-RELATIVISTIC STATE: ITS ==
!       == NORM IS 1-O((Z ALPHA)^2/4), E.G. 0.9975 FOR SI 1S, 0.990 FOR NI 1S,==
!       == 0.92 FOR AU 1S. RENORMALIZED. THE CHECK ONLY CATCHES BROKEN STATES. ==
        IF(ABS(NORM-1.D0).GT.0.15D0) THEN
          CALL ERROR$MSG('CORE STATE NOT NORMALIZED')
          CALL ERROR$I4VAL('IC',IC)
          CALL ERROR$R8VAL('NORM',NORM)
          CALL ERROR$STOP('FOCK_SPECIES')
        END IF
        AEPSI(:,IC)=AEPSI(:,IC)/SQRT(NORM)
        IF(ABS(FB(IC)-REAL(2*(2*LC+1),KIND=8)).GT.1.D-6) THEN
          CALL ERROR$MSG('CORE SHELL NOT FULLY OCCUPIED')
          CALL ERROR$I4VAL('IC',IC)
          CALL ERROR$R8VAL('F',FB(IC))
          CALL ERROR$STOP('FOCK_SPECIES')
        END IF
!       == RADIAL INTEGRALS FCV(A,B,L)=INT PHI_A PHI_C V_L[PHI_C PHI_B] ========
        ALLOCATE(FCV(LNX,LNX,0:LMAXC+2*LC,1),RHOC(NR,LNX))
        DO L=0,MAXVAL(SPF(ISP)%LOX)+LC
          DO LN1=1,LNX
            RHOC(:,LN1)=AEPHI(:,LN1)*AEPSI(:,IC)
          ENDDO
          DO LN2=1,LNX
            CALL RADIAL$POISSON(GID,NR,L,RHOC(:,LN2),VAE(:,1))
            DO LN1=1,LNX
              WORK(:)=R(:)**2*RHOC(:,LN1)*VAE(:,1)
              CALL RADIAL$INTEGRAL(GID,NR,WORK,FCV(LN1,LN2,L,1))
            ENDDO
          ENDDO
          IF(MU.GT.0.D0) THEN
            CALL FOCK_ERFPARTCV(GID,NR,R,L,LNX,RHOC,MU,FCV(:,:,L,1))
          END IF
        ENDDO
!       == ANGULAR SUM OVER MC AND THE MULTIPOLE LM ============================
        DO LMN1=1,LMNX
          LM1=SPF(ISP)%LMOF(LMN1)
          LN1=SPF(ISP)%LNOF(LMN1)
          L1=SPF(ISP)%LOX(LN1)
          DO LMN2=1,LMNX
            LM2=SPF(ISP)%LMOF(LMN2)
            LN2=SPF(ISP)%LNOF(LMN2)
            L2=SPF(ISP)%LOX(LN2)
            SVAR=0.D0
            DO MC=1,2*LC+1
              LMC=LC**2+MC
              DO L=ABS(L1-LC),L1+LC
                DO LM=L**2+1,(L+1)**2
                  CALL SPHERICAL$GAUNT(LM1,LMC,LM,CG)
                  IF(CG.EQ.0.D0) CYCLE
                  CALL SPHERICAL$GAUNT(LM2,LMC,LM,G)
                  SVAR=SVAR+CG*G*FCV(LN1,LN2,L,1)
                ENDDO
              ENDDO
            ENDDO
            SPF(ISP)%KCV(LMN1,LMN2)=SPF(ISP)%KCV(LMN1,LMN2)+SVAR
          ENDDO
        ENDDO
        DEALLOCATE(FCV,RHOC)
      ENDDO
      WRITE(*,FMT='("FOCK: SPECIES",I3," CORE STATES=",I3' &
     &          //'," KCV(DIAG)=",*(F9.5))')ISP,NCV,(SPF(ISP)%KCV(LMN1,LMN1),LMN1=1,LMNX)
!
!     == KEEP THE SHAPES: THEIR FOURIER TRANSFORMS ARE EVALUATED FOR EVERY  ==
!     == PAIR-DENSITY SET (|Q+G|) IN FOCK_PAWINIT                           ==
      SPF(ISP)%GID=GID
      SPF(ISP)%NR=NR
      ALLOCATE(SPF(ISP)%R(NR))
      ALLOCATE(SPF(ISP)%GHAT(NR,0:LMAXC))
      SPF(ISP)%R=R
      SPF(ISP)%GHAT=GHAT
      SPF(ISP)%TINI=.TRUE.
      WRITE(*,FMT='("FOCK: SPECIES",I3," LMNX=",I3," LMAXC=",I2," RC=",F8.4' &
     &           //'," RINT=",F8.4," Q0(DIAG)=",*(F9.5))') &
     &      ISP,LMNX,LMAXC,SPF(ISP)%RC,RINT,(QL(LN1,LN1,0),LN1=1,LNX)
      DEALLOCATE(R,WORK,AEPHI,PSPHI,GHAT,QL,GAUNT,IPAIR,RHOAE,RHOPS,VAE,VPS,F)
      DEALLOCATE(LB,FB,AEPSI)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_PAWINIT(NAT)
!     **************************************************************************
!     **  SPECIES DATA AND Y_L(G) ON FIRST USE                                **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4)            :: ISP,IAT,LMXCX,IQ,NG
!     **************************************************************************
      NATF=NAT
      ALLOCATE(ISPF(NATF))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NATF,ISPF)
      NSPF=MAXVAL(ISPF)
      ALLOCATE(SPF(NSPF))
      DO ISP=1,NSPF
        IF(.NOT.ANY(ISPF.EQ.ISP)) CYCLE
        CALL FOCK_SPECIES(ISP)
        LMAXCX=MAX(LMAXCX,SPF(ISP)%LMAXC)
      ENDDO
      ALLOCATE(IPRO0F(NATF))
      NPROF=0
      DO IAT=1,NATF
        IPRO0F(IAT)=NPROF
        NPROF=NPROF+SPF(ISPF(IAT))%LMNX
      ENDDO
      LMXCX=(LMAXCX+1)**2
!     == PER PAIR-DENSITY SET: Y_L(K), SHAPE TRANSFORMS AT |K|, STORAGE FOR  ==
!     == THE STRUCTURE FACTORS (K=Q+G)                                      ==
      DO IQ=1,NQF
        NG=QS(IQ)%NGL
        ALLOCATE(QS(IQ)%YLM(NG,LMXCX))
        ALLOCATE(QS(IQ)%FG(NG,0:LMAXCX,NSPF))
        ALLOCATE(QS(IQ)%FGD(NG,0:LMAXCX,NSPF))
        ALLOCATE(QS(IQ)%EIGR(NG,NATF))
        CALL FOCK_QSHAPES(IQ)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_QSHAPES(IQ)
!     **************************************************************************
!     **  Y_L(K) AND THE SHAPE TRANSFORMS FG, FGD AT |K| OF THE SET IQ        **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IQ
      INTEGER(4)            :: ISP,IG,LMXCX,NG
      REAL(8)               :: G,GU(3),PI,Y0
      REAL(8)   ,ALLOCATABLE:: GABS(:)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      LMXCX=(LMAXCX+1)**2
      NG=QS(IQ)%NGL
      ALLOCATE(GABS(NG))
      QS(IQ)%FG=0.D0
      QS(IQ)%FGD=0.D0
      DO IG=1,NG
        G=SQRT(SUM(QS(IQ)%GVEC(:,IG)**2))
        GABS(IG)=G
        IF(G.LT.1.D-8) THEN
          QS(IQ)%YLM(IG,:)=0.D0
          QS(IQ)%YLM(IG,1)=Y0
        ELSE
          GU(:)=QS(IQ)%GVEC(:,IG)/G
          CALL SPHERICAL$YLM(LMXCX,GU,QS(IQ)%YLM(IG,:))
        END IF
      ENDDO
      DO ISP=1,NSPF
        IF(.NOT.SPF(ISP)%TINI) CYCLE
        CALL FOCK_SHAPEFT(SPF(ISP)%GID,SPF(ISP)%NR,SPF(ISP)%R,SPF(ISP)%LMAXC &
     &                   ,SPF(ISP)%GHAT,NG,GABS,QS(IQ)%FG(:,0:SPF(ISP)%LMAXC,ISP) &
     &                   ,QS(IQ)%FGD(:,0:SPF(ISP)%LMAXC,ISP))
      ENDDO
      DEALLOCATE(GABS)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$NEWCELL(RBAS)
!     **************************************************************************
!     **  MOVING CELL: UPDATES THE PAIR-DENSITY SETS (SAME G INDICES, NEW     **
!     **  K VECTORS), THE KERNELS, V0, Y_L AND THE SHAPE TRANSFORMS. THE      **
!     **  STRUCTURE FACTORS FOLLOW IN FOCK$SNAPSHOTRESET.                     **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)               :: RBASS(3,3),G2,PI,FOURPI
      INTEGER(4)            :: IQ,IG,I
      CHARACTER(32)         :: OLDID
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(.NOT.TINI) RETURN
      IF(MAXVAL(ABS(RBAS-RBASF)).LT.1.D-12) RETURN
                              CALL TRACE$PUSH('FOCK$NEWCELL')
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      CALL PLANEWAVE$GETCH('ID',OLDID)
      RBASF=RBAS
      DO I=1,3
        RBASS(:,I)=RBAS(:,I)*REAL(NKDIVF(I),KIND=8)
      ENDDO
      CALL FOCK_V0(RBASS,ALPHAF,V0)
      DO IQ=1,NQF
        CALL PLANEWAVE$SELECT(QS(IQ)%ID)
        CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
        CALL PLANEWAVE$GETR8A('GVEC',3*QS(IQ)%NGL,QS(IQ)%GVEC)
        IF(IQ.EQ.1) THEN
          CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHTF)
          CALL PLANEWAVE$GETR8('RWEIGHT',RWEIGHTF)
        END IF
        DO IG=1,QS(IQ)%NGL
          G2=SUM(QS(IQ)%GVEC(:,IG)**2)
          IF(G2.LT.1.D-10) THEN
            QS(IQ)%VKER(IG)=V0
          ELSE IF(MU.GT.0.D0) THEN
            QS(IQ)%VKER(IG)=FOURPI/G2*(1.D0-EXP(-G2/(4.D0*MU**2)))
          ELSE
            QS(IQ)%VKER(IG)=FOURPI/G2
          END IF
        ENDDO
        IF(ALLOCATED(SPF)) CALL FOCK_QSHAPES(IQ)
      ENDDO
      IF(LEN_TRIM(OLDID).GT.0) CALL PLANEWAVE$SELECT(TRIM(OLDID))
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$SNAPSHOTRESET(NSPIN,NDIM,NOCCX,NAT,NPRO,R)
!     **************************************************************************
!     **  CLEARS THE ORBITAL SNAPSHOT, RESERVES ROOM FOR NOCCX ORBITALS PER   **
!     **  SPIN AND UPDATES THE STRUCTURE FACTORS FOR THE ATOMIC POSITIONS R.  **
!     **  NDIM=2: TWO-COMPONENT SPINORS (NONCOLLINEAR, NSPIN=1), STORED WITH  **
!     **  THE COMPONENTS STACKED; THE PAIR DENSITY SUMS OVER THE COMPONENTS.  **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSPIN
      INTEGER(4),INTENT(IN) :: NDIM
      INTEGER(4),INTENT(IN) :: NOCCX
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: NPRO
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)            :: IAT,IQ,IG
      REAL(8)               :: RBINV(3,3)
!     **************************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('FOCK OBJECT NOT INITIALIZED')
        CALL ERROR$STOP('FOCK$SNAPSHOTRESET')
      END IF
      IF(.NOT.ALLOCATED(SPF)) CALL FOCK_PAWINIT(NAT)
      IF(.NOT.ALLOCATED(UROTF)) CALL FOCK_ROTINIT()
      IF(NAT.NE.NATF.OR.NPRO.NE.NPROF) THEN
        CALL ERROR$MSG('INCONSISTENT NUMBER OF ATOMS OR PROJECTORS')
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$I4VAL('NATF',NATF)
        CALL ERROR$I4VAL('NPRO',NPRO)
        CALL ERROR$I4VAL('NPROF',NPROF)
        CALL ERROR$STOP('FOCK$SNAPSHOTRESET')
      END IF
!     == RELATIVE POSITIONS (LATTICE TRANSLATIONS OF THE SYMMETRY MAP) ========
      IF(.NOT.ALLOCATED(XPOSF)) ALLOCATE(XPOSF(3,NATF))
      CALL LIB$INVERTR8(3,RBASF,RBINV)
      XPOSF=MATMUL(RBINV,R)
!     == STRUCTURE FACTORS EXP(-I(Q+G)R) OF ALL PAIR-DENSITY SETS ============
      DO IQ=1,NQF
        DO IAT=1,NATF
          DO IG=1,QS(IQ)%NGL
            QS(IQ)%EIGR(IG,IAT)=EXP(-CI*SUM(QS(IQ)%GVEC(:,IG)*R(:,IAT)))
          ENDDO
        ENDDO
      ENDDO
      IF(NDIM.NE.1.AND.(NDIM.NE.2.OR.NSPIN.NE.1)) THEN
        CALL ERROR$MSG('UNSUPPORTED SPIN SETTING')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('FOCK$SNAPSHOTRESET')
      END IF
      IF(NDIM.EQ.2.AND.NOPF.GT.1) THEN
        CALL ERROR$MSG('SPACE-GROUP SYMMETRY WITH SPINORS IS NOT IMPLEMENTED')
        CALL ERROR$STOP('FOCK$SNAPSHOTRESET')
      END IF
      IF(NSPINF.NE.NSPIN.OR.NOCCXF.NE.NOCCX.OR.NDIMF.NE.NDIM.OR.NRLD.EQ.0) THEN
        IF(ALLOCATED(NOCC)) DEALLOCATE(NOCC)
        IF(ALLOCATED(WOCC)) DEALLOCATE(WOCC)
        IF(ALLOCATED(PHIR)) DEALLOCATE(PHIR)
        IF(ALLOCATED(CPROJ)) DEALLOCATE(CPROJ)
        IF(ALLOCATED(IQOCC)) DEALLOCATE(IQOCC)
        NSPINF=NSPIN
        NDIMF=NDIM
        NRLD=NRLF*NDIMF
        NPROD=NPROF*NDIMF
        NOCCXF=NOCCX
        ALLOCATE(NOCC(NSPINF))
        ALLOCATE(WOCC(NOCCXF,NSPINF))
        ALLOCATE(PHIR(NRLD,NOCCXF,NSPINF))
        ALLOCATE(CPROJ(NPROD,NOCCXF,NSPINF))
        ALLOCATE(IQOCC(NOCCXF,NSPINF))
!       == ACE STORAGE IS SIZED WITH NSPINF ====================================
        IF(ALLOCATED(ACES)) DEALLOCATE(ACES)
        ALLOCATE(ACES(NKSETF,NSPINF))
      END IF
      NOCC(:)=0
      WOCC(:,:)=0.D0
      EXF=0.D0
!     == THE ACE OPERATOR BELONGS TO THE PREVIOUS SNAPSHOT: INVALIDATE =========
      IF(ALLOCATED(ACES)) ACES(:,:)%N=0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$SNAPSHOT(GSETID,ISPIN,NGL,NDIM,NBH,NB,NPRO,TINV,F,PSI,PROJ)
!     **************************************************************************
!     **  STORES THE OCCUPIED ORBITALS OF ONE WAVE-FUNCTION SET (K,SPIN) AS   **
!     **  PERIODIC PARTS U_K(R) TOGETHER WITH THEIR PROJECTIONS.              **
!     **  TINV (K=-K MOD G): TWO REAL BLOCH FUNCTIONS 2I-1, 2I ARE PACKED AS  **
!     **  RE AND IM OF EXP(IKR) U(R) (AS IN WAVES_DENSITY); THEIR PROJECTIONS **
!     **  ARE REAL.                                                           **
!     **  EVERY ORBITAL IS STORED FOR ALL MEMBERS Q=S MINVT K OF THE STAR OF **
!     **  K (FOCK$SETKMESH), EACH WITH THE WEIGHT DIVIDED BY THE STAR SIZE:   **
!     **  ROTATED BY FOCK_ROTORB/FOCK_ROTPROJ, AND COMPLEX CONJUGATED FOR     **
!     **  S=-1 (TIME REVERSAL).                                               **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: GSETID
      INTEGER(4)  ,INTENT(IN) :: ISPIN
      INTEGER(4)  ,INTENT(IN) :: NGL
      INTEGER(4)  ,INTENT(IN) :: NDIM
      INTEGER(4)  ,INTENT(IN) :: NBH
      INTEGER(4)  ,INTENT(IN) :: NB
      INTEGER(4)  ,INTENT(IN) :: NPRO
      LOGICAL(4)  ,INTENT(IN) :: TINV
      REAL(8)     ,INTENT(IN) :: F(NB)
      COMPLEX(8)  ,INTENT(IN) :: PSI(NGL,NDIM,NBH)
      COMPLEX(8)  ,INTENT(IN) :: PROJ(NDIM,NBH,NPRO)
      COMPLEX(8)  ,ALLOCATABLE:: FR(:,:),FRS(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: FG(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: EIKR(:),PHI(:),CB(:)
      COMPLEX(8)  ,ALLOCATABLE:: FBOX(:,:,:),UROT(:),CROT(:)
      REAL(8)                 :: SFAC,XK(3),HALF
      REAL(8)     ,PARAMETER  :: FMIN=1.D-10
      INTEGER(4)              :: IB0,NB1,K,IBH,IB,N,IK,IQK,N0,IDIM
      LOGICAL(4)              :: TROT
      CHARACTER(16)           :: CENV
!     **************************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('FOCK OBJECT NOT INITIALIZED')
        CALL ERROR$STOP('FOCK$SNAPSHOT')
      END IF
      IF(ISPIN.LT.1.OR.ISPIN.GT.NSPINF) THEN
        CALL ERROR$MSG('SPIN INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('ISPIN',ISPIN)
        CALL ERROR$STOP('FOCK$SNAPSHOT')
      END IF
      IF(NDIM.NE.NDIMF.OR.(TINV.AND.NDIM.NE.1)) THEN
        CALL ERROR$MSG('INCONSISTENT NDIM (OR TINV WITH SPINORS)')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$I4VAL('NDIMF',NDIMF)
        CALL ERROR$STOP('FOCK$SNAPSHOT')
      END IF
      SFAC=1.D0
      IF(NSPINF.EQ.1.AND.NDIMF.EQ.1) SFAC=0.5D0
      CALL PLANEWAVE$SELECT(GSETID)
      CALL PLANEWAVE$GETR8A('XK',3,XK)
      IK=FOCK_FINDK(XK,NKSETF,XKSETF)
      IQK=FOCK_FINDK(XK,NQFULL,XQFULL)
      IF(IK.EQ.0.OR.IQK.EQ.0) THEN
        CALL ERROR$MSG('K-POINT OF THE SET NOT IN THE FOCK K LIST')
        CALL ERROR$STOP('FOCK$SNAPSHOT')
      END IF
      HALF=1.D0/REAL(NMEMF(IK),KIND=8)
      TROT=.NOT.ALL(MEMIDF(1:NMEMF(IK),IK))
      IF(TROT) THEN
        ALLOCATE(FBOX(NR1F,NR2F,NR3F))
        ALLOCATE(UROT(NRLF))
        ALLOCATE(CROT(NPROF))
      END IF
      ALLOCATE(FG(NGL,NBLK))
      ALLOCATE(FR(NRLD,NBLK))
      ALLOCATE(FRS(NRLF,NBLK))
      ALLOCATE(PHI(NRLF))
      ALLOCATE(CB(NPROD))
      IF(TINV) THEN
        ALLOCATE(EIKR(NRLF))
        CALL PLANEWAVE$GETC8A('EIKR',NRLF,EIKR)
      END IF
      N=NOCC(ISPIN)
      N0=N
      DO IB0=1,NBH,NBLK
        NB1=MIN(NBLK,NBH-IB0+1)
        DO IDIM=1,NDIM
          FG(:,1:NB1)=PSI(:,IDIM,IB0:IB0+NB1-1)
          CALL PLANEWAVE$FFT('GTOR',NB1,NGL,FG(:,1:NB1),NRLF,FRS(:,1:NB1))
          FR((IDIM-1)*NRLF+1:IDIM*NRLF,1:NB1)=FRS(:,1:NB1)
        ENDDO
        DO K=1,NB1
          IBH=IB0+K-1
          IF(TINV) THEN
            FR(:,K)=FR(:,K)*EIKR(:)          ! PACKED REAL BLOCH FUNCTIONS
            IB=2*IBH-1
            IF(F(IB).GT.FMIN) THEN
              PHI=CMPLX(REAL(FR(:,K),KIND=8),0.D0,KIND=8)*CONJG(EIKR)
              CALL ADDORB(HALF*SFAC*F(IB),PHI &
     &                   ,CMPLX(REAL(PROJ(1,IBH,:),KIND=8),0.D0,KIND=8))
            END IF
            IB=2*IBH
            IF(IB.LE.NB) THEN
              IF(F(IB).GT.FMIN) THEN
                PHI=CMPLX(AIMAG(FR(:,K)),0.D0,KIND=8)*CONJG(EIKR)
                CALL ADDORB(HALF*SFAC*F(IB),PHI &
     &                     ,CMPLX(AIMAG(PROJ(1,IBH,:)),0.D0,KIND=8))
              END IF
            END IF
          ELSE
            IF(F(IBH).GT.FMIN) THEN
              DO IDIM=1,NDIM
                CB((IDIM-1)*NPROF+1:IDIM*NPROF)=PROJ(IDIM,IBH,:)
              ENDDO
              CALL ADDORB(HALF*SFAC*F(IBH),FR(:,K),CB)
            END IF
          END IF
        ENDDO
      ENDDO
      NOCC(ISPIN)=N
      DEALLOCATE(FG)
      DEALLOCATE(FR,FRS)
      DEALLOCATE(PHI,CB)
      IF(TINV) DEALLOCATE(EIKR)
      IF(TROT) DEALLOCATE(FBOX,UROT,CROT)
!     == DIAGNOSTIC: THE LITTLE GROUP OF K MAPS THE OCCUPIED ORBITALS ONTO ==
!     == THEIR OWN SPAN (CHECKS THE ROTATION OF ORBITALS AND PROJECTIONS)  ==
      IF(NOPF.GT.1) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_FOCK_CHECK',CENV)
        IF(TRIM(CENV).EQ.'1') CALL FOCK_SYMCHECK(ISPIN,IK,N0,N)
      END IF
      RETURN
      CONTAINS
!       ........................................................................
        SUBROUTINE ADDORB(W,U,C)
        REAL(8)   ,INTENT(IN) :: W
        COMPLEX(8),INTENT(IN) :: U(NRLD)
        COMPLEX(8),INTENT(IN) :: C(NPROD)
        INTEGER(4)            :: IM,IOP,IQ
        IF(TROT) CALL FOCK_ORBTOBOX(U,FBOX)
        DO IM=1,NMEMF(IK)
          IOP=MEMOPF(IM,IK)
          IQ=MEMQF(IM,IK)
          IF(MEMIDF(IM,IK)) THEN
            IF(MEMSF(IM,IK).EQ.1) THEN
              CALL STOREORB(W,U,C,IQ)
            ELSE
              CALL STOREORB(W,CONJG(U),CONJG(C),IQ)
            END IF
          ELSE
!           == XQFULL(IQ) = S MINVT K EXACTLY (NO FOLDING) =====================
            CALL FOCK_ROTORB(FBOX,IOP,REAL(MEMSF(IM,IK),KIND=8)*XQFULL(:,IQ) &
     &                      ,(/0,0,0/),UROT)
            CALL FOCK_ROTPROJ(C,IOP,REAL(MEMSF(IM,IK),KIND=8)*XQFULL(:,IQ),CROT)
            IF(MEMSF(IM,IK).EQ.1) THEN
              CALL STOREORB(W,UROT,CROT,IQ)
            ELSE
              CALL STOREORB(W,CONJG(UROT),CONJG(CROT),IQ)
            END IF
          END IF
        ENDDO
        END SUBROUTINE ADDORB
!       ........................................................................
        SUBROUTINE STOREORB(W,U,C,IQ)
        REAL(8)   ,INTENT(IN) :: W
        COMPLEX(8),INTENT(IN) :: U(NRLD)
        COMPLEX(8),INTENT(IN) :: C(NPROD)
        INTEGER(4),INTENT(IN) :: IQ
        N=N+1
        IF(N.GT.NOCCXF) THEN
          CALL ERROR$MSG('TOO MANY OCCUPIED ORBITALS FOR THE SNAPSHOT')
          CALL ERROR$I4VAL('NOCCXF',NOCCXF)
          CALL ERROR$STOP('FOCK$SNAPSHOT')
        END IF
        WOCC(N,ISPIN)=W
        IQOCC(N,ISPIN)=IQ
        PHIR(:,N,ISPIN)=U(:)
        CPROJ(:,N,ISPIN)=C(:)
        END SUBROUTINE STOREORB
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ROTINIT()
!     **************************************************************************
!     **  ROTATION MATRICES OF THE PROJECTOR FUNCTIONS FOR EVERY OPERATION    **
!     **  AND SPECIES: P_LM(C S) = SUM_M' U_MM' P_LM'(S) IN EACH PARTIAL-WAVE **
!     **  BLOCK, (-1)^L FOR IMPROPER OPERATIONS (AS SYMMETRY$DENMAT).         **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      REAL(8)   ,ALLOCATABLE :: YLMROT(:,:)
      REAL(8)                :: C(3,3),DET,SGN
      INTEGER(4)             :: LMNXX,LX,LMX,IOP,ISP,LMN,LN,L
!     **************************************************************************
      LMNXX=0
      LX=0
      DO ISP=1,NSPF
        LMNXX=MAX(LMNXX,SPF(ISP)%LMNX)
        LX=MAX(LX,MAXVAL(SPF(ISP)%LOX))
      ENDDO
      LMX=(LX+1)**2
      IF(ALLOCATED(UROTF)) DEALLOCATE(UROTF)
      ALLOCATE(UROTF(LMNXX,LMNXX,NSPF,NOPF))
      ALLOCATE(YLMROT(LMX,LMX))
      UROTF=0.D0
      DO IOP=1,NOPF
        C=CROTF(:,:,IOP)
        DET=C(1,1)*(C(2,2)*C(3,3)-C(3,2)*C(2,3))-C(1,2)*(C(2,1)*C(3,3)-C(3,1)*C(2,3)) &
     &     +C(1,3)*(C(2,1)*C(3,2)-C(3,1)*C(2,2))
        IF(DET.LT.0.D0) C=-C
        CALL SPHERICAL$ROTATEYLM(LMX,C,YLMROT)
        DO ISP=1,NSPF
          LMN=0
          DO LN=1,SPF(ISP)%LNX
            L=SPF(ISP)%LOX(LN)
            SGN=1.D0
            IF(DET.LT.0.D0.AND.MOD(L,2).EQ.1) SGN=-1.D0
            UROTF(LMN+1:LMN+2*L+1,LMN+1:LMN+2*L+1,ISP,IOP) &
     &             =SGN*YLMROT(L**2+1:(L+1)**2,L**2+1:(L+1)**2)
            LMN=LMN+2*L+1
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(YLMROT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ORBTOBOX(U,FB)
!     **************************************************************************
!     **  FOURIER COEFFICIENTS U(R)=SUM_H FB(H) EXP(2PI I H.X) OF A PERIODIC  **
!     **  FUNCTION GIVEN ON THE LOCAL GRID SLICE, ON THE WHOLE FFT BOX        **
!     **  (ASSEMBLED OVER 'K', SAME ON ALL TASKS OF THE GROUP)                **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN)  :: U(NRLF)
      COMPLEX(8),INTENT(OUT) :: FB(NR1F,NR2F,NR3F)
!     **************************************************************************
      FB=(0.D0,0.D0)
      FB(NR1STARTF:NR1STARTF+NR1LF-1,:,:)=RESHAPE(U,(/NR1LF,NR2F,NR3F/))
      CALL MPE$COMBINE('K','+',FB)
      CALL SYMMETRY_FFT3D('RTOG',(/NR1F,NR2F,NR3F/),FB)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ROTORB(FB,IOP,XKP,IGF,UOUT)
!     **************************************************************************
!     **  PERIODIC PART OF THE ROTATED ORBITAL PSI'(X)=PSI(M^-1(X-TAU)) WITH  **
!     **  BLOCH VECTOR K'=MINVT K (RELATIVE), SHIFTED BY THE RECIPROCAL       **
!     **  LATTICE VECTOR IGF (BLOCH VECTOR K'-IGF):                           **
!     **    U'(X) = SUM_H C_H EXP(-2PI I (K'+H').TAU) EXP(2PI I (H'+IGF).X)   **
!     **  WITH H'=MINVT H AND C_H THE COEFFICIENTS OF U_K (FOCK_ORBTOBOX).    **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN)  :: FB(NR1F,NR2F,NR3F)
      INTEGER(4),INTENT(IN)  :: IOP
      REAL(8)   ,INTENT(IN)  :: XKP(3)
      INTEGER(4),INTENT(IN)  :: IGF(3)
      COMPLEX(8),INTENT(OUT) :: UOUT(NRLF)
      REAL(8)   ,PARAMETER   :: PI=3.141592653589793D0
      COMPLEX(8),ALLOCATABLE :: FS(:,:,:)
      INTEGER(4)             :: NR(3),H(3),HP(3),I1,I2,I3,I
      REAL(8)                :: PH
!     **************************************************************************
      NR=(/NR1F,NR2F,NR3F/)
      ALLOCATE(FS(NR1F,NR2F,NR3F))
      FS=(0.D0,0.D0)
      DO I3=0,NR(3)-1
        DO I2=0,NR(2)-1
          DO I1=0,NR(1)-1
            IF(FB(I1+1,I2+1,I3+1).EQ.(0.D0,0.D0)) CYCLE
            H=(/I1,I2,I3/)
            DO I=1,3
              IF(2*H(I).GT.NR(I)) H(I)=H(I)-NR(I)
            ENDDO
            HP=MATMUL(MINVTF(:,:,IOP),H)
            PH=-2.D0*PI*DOT_PRODUCT(XKP+REAL(HP,KIND=8),TAUF(:,IOP))
            HP=MODULO(HP+IGF,NR)
            FS(HP(1)+1,HP(2)+1,HP(3)+1)=FS(HP(1)+1,HP(2)+1,HP(3)+1) &
     &                  +FB(I1+1,I2+1,I3+1)*CMPLX(COS(PH),SIN(PH),KIND=8)
          ENDDO
        ENDDO
      ENDDO
      CALL SYMMETRY_FFT3D('GTOR',NR,FS)
      UOUT=RESHAPE(FS(NR1STARTF:NR1STARTF+NR1LF-1,:,:),(/NRLF/))
      DEALLOCATE(FS)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ROTPROJ(C,IOP,XKP,COUT)
!     **************************************************************************
!     **  PROJECTIONS OF THE ROTATED ORBITAL (FOCK_ROTORB) WITH BLOCH VECTOR  **
!     **  K': ATOM J=MAPAT(I) AT X_J = M X_I + TAU + L RECEIVES                **
!     **    <P_J|PSI'> = EXP(2PI I K'.L) U <P_I|PSI>                          **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN)  :: C(NPROF)
      INTEGER(4),INTENT(IN)  :: IOP
      REAL(8)   ,INTENT(IN)  :: XKP(3)
      COMPLEX(8),INTENT(OUT) :: COUT(NPROF)
      REAL(8)   ,PARAMETER   :: PI=3.141592653589793D0
      REAL(8)                :: XL(3),PH
      INTEGER(4)             :: IAT,JAT,ISP,LMNX
!     **************************************************************************
      DO IAT=1,NATF
        JAT=MAPATF(IAT,IOP)
        ISP=ISPF(IAT)
        LMNX=SPF(ISP)%LMNX
        XL=XPOSF(:,JAT)-MATMUL(REAL(IROTF(:,:,IOP),KIND=8),XPOSF(:,IAT))-TAUF(:,IOP)
        IF(MAXVAL(ABS(XL-ANINT(XL))).GT.1.D-3) THEN
          CALL ERROR$MSG('ATOMIC POSITIONS HAVE LOST THE SYMMETRY')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$I4VAL('IOP',IOP)
          CALL ERROR$STOP('FOCK_ROTPROJ')
        END IF
        PH=2.D0*PI*DOT_PRODUCT(XKP,ANINT(XL))
        COUT(IPRO0F(JAT)+1:IPRO0F(JAT)+LMNX)=CMPLX(COS(PH),SIN(PH),KIND=8) &
     &       *MATMUL(UROTF(1:LMNX,1:LMNX,ISP,IOP),C(IPRO0F(IAT)+1:IPRO0F(IAT)+LMNX))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_SYMCHECK(ISPIN,IK,N0,N1)
!     **************************************************************************
!     **  DIAGNOSTIC (VESPAW_FOCK_CHECK=1): EVERY OPERATION OF THE LITTLE     **
!     **  GROUP OF K (MINVT K = K + G) MAPS THE OCCUPIED ORBITALS OF K ONTO   **
!     **  THEIR SPAN, U' = U A AND <P|PSI'> = <P|PSI> A WITH ONE MATRIX A     **
!     **  (LEAST SQUARES FROM THE PLANE-WAVE PART). PRINTS THE RESIDUALS.     **
!     **  MEANINGFUL FOR FULLY OCCUPIED (CLOSED) SUBSPACES ONLY.              **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISPIN,IK,N0,N1
      COMPLEX(8),ALLOCATABLE:: FB(:,:,:),UR(:,:),CR(:,:),U(:,:),C(:,:)
      COMPLEX(8),ALLOCATABLE:: S(:,:),B(:,:),A(:,:)
      INTEGER(4),ALLOCATABLE:: IDX(:)
      REAL(8)               :: XK(3),XP(3),D(3),RU,RP,RUX,RPX,NU,NP
      INTEGER(4)            :: NV,N,IV,IOP,NLG
!     **************************************************************************
      ALLOCATE(IDX(N1-N0))
      NV=0
      DO N=N0+1,N1
        IF(IQOCC(N,ISPIN).NE.MEMQF(1,IK)) CYCLE
        NV=NV+1
        IDX(NV)=N
      ENDDO
      IF(NV.EQ.0) RETURN
      ALLOCATE(FB(NR1F,NR2F,NR3F),UR(NRLF,NV),CR(NPROF,NV),U(NRLF,NV),C(NPROF,NV))
      ALLOCATE(S(NV,NV),B(NV,NV),A(NV,NV))
      DO IV=1,NV
        U(:,IV)=PHIR(:,IDX(IV),ISPIN)
        C(:,IV)=CPROJ(:,IDX(IV),ISPIN)
      ENDDO
      S=MATMUL(CONJG(TRANSPOSE(U)),U)
      CALL MPE$COMBINE('K','+',S)
      XK=XKSETF(:,IK)
      RUX=0.D0
      RPX=0.D0
      NLG=0
      DO IOP=1,NOPF
        XP=MATMUL(REAL(MINVTF(:,:,IOP),KIND=8),XK)
        D=XP-XK
        IF(MAXVAL(ABS(D-ANINT(D))).GT.1.D-8) CYCLE
        NLG=NLG+1
        DO IV=1,NV
          CALL FOCK_ORBTOBOX(U(:,IV),FB)
          CALL FOCK_ROTORB(FB,IOP,XP,NINT(D),UR(:,IV))
          CALL FOCK_ROTPROJ(C(:,IV),IOP,XP,CR(:,IV))
        ENDDO
        B=MATMUL(CONJG(TRANSPOSE(U)),UR)
        CALL MPE$COMBINE('K','+',B)
        CALL LIB$MATRIXSOLVEC8(NV,NV,NV,S,A,B)
        UR=UR-MATMUL(U,A)
        CR=CR-MATMUL(C,A)
        RU=SUM(ABS(UR)**2)
        NU=SUM(ABS(U)**2)
        CALL MPE$COMBINE('K','+',RU)
        CALL MPE$COMBINE('K','+',NU)
        RP=SUM(ABS(CR)**2)
        NP=SUM(ABS(C)**2)
        RUX=MAX(RUX,SQRT(RU/NU))
        RPX=MAX(RPX,SQRT(RP/NP))
      ENDDO
      WRITE(*,FMT='("FOCK SYM CHECK: K",I5," SPIN",I2," #BANDS",I5," #LITTLE GROUP",I4' &
     &   //'," MAX REL. RESIDUAL ORBITALS",ES10.2," PROJECTIONS",ES10.2)')IK,ISPIN,NV,NLG,RUX,RPX
      DEALLOCATE(FB,UR,CR,U,C,S,B,A,IDX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$SNAPSHOTGATHER()
!     **************************************************************************
!     **  COMPLETES THE SNAPSHOT WITH THE OCCUPIED ORBITALS OF THE OTHER      **
!     **  K-GROUPS. ALL K-GROUPS HAVE THE SAME SIZE AND REAL-SPACE LAYOUT     **
!     **  (WAVES_KDISTRIBUTEFOCK), SO THE TASKS OF EQUAL RANK (COMMUNICATOR   **
!     **  'KX') HOLD THE SAME GRID SLICE AND EXCHANGE THEIR ENTRIES DIRECTLY. **
!     **  THE ORDER OF THE ENTRIES IS THE ORDER OF THE GROUPS. NO-OP WITH ONE **
!     **  K-GROUP.                                                            **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: NTX,ITX,ISPIN,IT,I0,N,I
      INTEGER(4),ALLOCATABLE:: NOF(:)
      REAL(8)   ,ALLOCATABLE:: WBUF(:)
      INTEGER(4),ALLOCATABLE:: IBUF(:)
      COMPLEX(8),ALLOCATABLE:: RBUF(:,:),PBUF(:,:)
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL MPE$QUERY('KX',NTX,ITX)
      IF(NTX.EQ.1) RETURN
                              CALL TRACE$PUSH('FOCK$SNAPSHOTGATHER')
                              CALL TIMING$CLOCKON('FOCK$SNAPSHOTGATHER')
      ALLOCATE(NOF(NTX))
      DO ISPIN=1,NSPINF
        NOF(:)=0
        NOF(ITX)=NOCC(ISPIN)
        CALL MPE$COMBINE('KX','+',NOF)
        IF(SUM(NOF).GT.NOCCXF) THEN
          CALL ERROR$MSG('TOO MANY OCCUPIED ORBITALS FOR THE SNAPSHOT')
          CALL ERROR$I4VAL('NOCCXF',NOCCXF)
          CALL ERROR$I4VAL('SUM OVER K-GROUPS',SUM(NOF))
          CALL ERROR$STOP('FOCK$SNAPSHOTGATHER')
        END IF
!       == MOVE THE OWN ENTRIES TO THEIR FINAL POSITION (BACKWARDS) ==========
        I0=SUM(NOF(1:ITX-1))
        IF(I0.GT.0) THEN
          DO I=NOF(ITX),1,-1
            WOCC(I0+I,ISPIN)=WOCC(I,ISPIN)
            IQOCC(I0+I,ISPIN)=IQOCC(I,ISPIN)
            PHIR(:,I0+I,ISPIN)=PHIR(:,I,ISPIN)
            CPROJ(:,I0+I,ISPIN)=CPROJ(:,I,ISPIN)
          ENDDO
        END IF
!       == EVERY GROUP BROADCASTS ITS ENTRIES ===================================
        I0=0
        DO IT=1,NTX
          N=NOF(IT)
          IF(N.GT.0) THEN
            ALLOCATE(WBUF(N),IBUF(N),RBUF(NRLD,N),PBUF(NPROD,N))
            IF(IT.EQ.ITX) THEN
              WBUF=WOCC(I0+1:I0+N,ISPIN)
              IBUF=IQOCC(I0+1:I0+N,ISPIN)
              RBUF=PHIR(:,I0+1:I0+N,ISPIN)
              PBUF=CPROJ(:,I0+1:I0+N,ISPIN)
            END IF
            CALL MPE$BROADCAST('KX',IT,WBUF)
            CALL MPE$BROADCAST('KX',IT,IBUF)
            CALL MPE$BROADCAST('KX',IT,RBUF)
            CALL MPE$BROADCAST('KX',IT,PBUF)
            IF(IT.NE.ITX) THEN
              WOCC(I0+1:I0+N,ISPIN)=WBUF
              IQOCC(I0+1:I0+N,ISPIN)=IBUF
              PHIR(:,I0+1:I0+N,ISPIN)=RBUF
              CPROJ(:,I0+1:I0+N,ISPIN)=PBUF
            END IF
            DEALLOCATE(WBUF,IBUF,RBUF,PBUF)
          END IF
          I0=I0+N
        ENDDO
        NOCC(ISPIN)=I0
      ENDDO
      DEALLOCATE(NOF)
                              CALL TIMING$CLOCKOFF('FOCK$SNAPSHOTGATHER')
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$WBUILDDONE()
!     **************************************************************************
!     **  AFTER FOCK$WBUILD FOR ALL LOCAL SETS: SUMS THE PLANE-WAVE FOCK      **
!     **  ENERGY OVER THE K-GROUPS (ONE TASK PER GROUP IN 'KX'; INSIDE A      **
!     **  GROUP EVERY TASK ALREADY HOLDS THE GROUP SUM).                      **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL MPE$COMBINE('KX','+',EXF)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$FORCE(NAT,TSTRESS,FORCE,STRESS)
!     **************************************************************************
!     **  EXPLICIT GEOMETRY DEPENDENCE OF THE PLANE-WAVE FOCK ENERGY OF THE   **
!     **  SNAPSHOT (THE DEPENDENCE THROUGH THE PROJECTIONS <P|PSI> IS TAKEN   **
!     **  CARE OF BY WAVES$FORCE WITH THE PROJECTOR PART OF FOCK$APPLY):      **
!     **                                                                      **
!     **    E = SUM_PAIRS LAMBDA_MN T_MN,  T = OMEGA SUM_K V(K) |RHO_MN(K)|^2, **
!     **    LAMBDA_MN = -AEXX/2 W_M W_N / S  (TWICE FOR M<N, PAIRS M<=N)      **
!     **    RHO_MN = CONJG(U_M) U_N EXP(I G0 R) + RHOHAT_MN, RHOHAT FROM      **
!     **    Q^{RL}_MN = SUM_IJ CONJG(<P_I|M>) <P_J|N> T_IJ^L AND GHAT_L(R-R_R)  **
!     **                                                                      **
!     **  FORCE (=-DE/DR, ONLY THE COMPENSATION CHARGES MOVE WITH THE ATOM):  **
!     **    F_RI = -2 RE SUM_K (-I K_I) EXP(-IKR) SUM_L (-I)^L F_L Y_L Z^{RL}(K)**
!     **    Z^{RL}(K) = SUM_PAIRS LAMBDA CONJG(V(K) RHO(K)) Q^{RL}             **
!     **  STRESS (=+DE/DEPSILON_IJ) AT FIXED COEFFICIENTS SQRT(OMEGA) PSI(G)  **
!     **  AND FIXED PROJECTIONS (K -> (1-EPS)K, OMEGA -> (1+TR EPS) OMEGA):   **
!     **    -DELTA_IJ E                                   (VOLUME)            **
!     **    +SUM LAMBDA OMEGA SUM_K (-K DV/DK) K_IK_J/K^2 |RHO|^2  (KERNEL)    **
!     **    +SUM LAMBDA OMEGA |RHO(0)|^2 DV0/DEPS_IJ      (G=0 TERM, NUMERIC) **
!     **    +2 RE SUM_K EXP(-IKR)(-I)^L DS_L/DEPS_IJ Z^{RL}  (SHAPES)           **
!     **    DS_L/DEPS_IJ = -[F_L' K K_IK_J/K^2 Y_L + F_L SYM(K_I D_J Y_L)]     **
!     **  THE PAIRS ARE DISTRIBUTED OVER 'KX', THE G VECTORS OVER 'K'.        **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NAT
      LOGICAL(4),INTENT(IN)  :: TSTRESS
      REAL(8)   ,INTENT(OUT) :: FORCE(3,NAT)
      REAL(8)   ,INTENT(OUT) :: STRESS(3,3)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER   :: RSMALL=1.D-20
      COMPLEX(8),ALLOCATABLE :: Z(:,:,:)       ! (NGLQX,LMXCX,NATF)
      COMPLEX(8),ALLOCATABLE :: RHOR(:,:),RHOG(:,:),VGL(:,:)
      COMPLEX(8),ALLOCATABLE :: QL(:,:,:)      ! (NBLK,LMXCX,NATF)
      COMPLEX(8),ALLOCATABLE :: EG0(:),AUX(:)
      REAL(8)   ,ALLOCATABLE :: GIJ(:,:),DVK(:),GABS(:),SYLM(:,:,:)
      INTEGER(4),ALLOCATABLE :: MB(:)
      REAL(8)   ,ALLOCATABLE :: G0B(:,:),LAMB(:)
      REAL(8)                :: SFAC,EPAIR,SKER(6),SSHP(6),SV0,DV0(3,3)
      REAL(8)                :: XK(3),G0(3),SVAR,PI,FOURPI,C,FR(3)
      REAL(8)                :: RBASS(3,3),DEPS(3,3),VP,VM
      REAL(8)   ,PARAMETER   :: HEPS=1.D-4
      INTEGER(4)             :: NTX,ITX,LMXCX,NGLQX,IQ,ISPIN,N,M,NB1,K,IG,IJ
      INTEGER(4)             :: I,J,IAT,ISP,LMXC,LM,L,IQM,NPAIR
      COMPLEX(8)             :: CSVAR,CFAC
      LOGICAL(4)             :: TPHASE,TSTARF
      CHARACTER(16)          :: CENV
      INTEGER(4)             :: NFILO,MLAST
      REAL(8)                :: SMULT
!     **************************************************************************
      FORCE(:,:)=0.D0
      STRESS(:,:)=0.D0
      IF(.NOT.TON) RETURN
      IF(.NOT.ALLOCATED(NOCC)) RETURN
      IF(NAT.NE.NATF) THEN
        CALL ERROR$MSG('INCONSISTENT NUMBER OF ATOMS')
        CALL ERROR$STOP('FOCK$FORCE')
      END IF
                              CALL TRACE$PUSH('FOCK$FORCE')
                              CALL TIMING$CLOCKON('FOCK$FORCE')
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      CALL MPE$QUERY('KX',NTX,ITX)
      SFAC=1.D0
      IF(NSPINF.EQ.1.AND.NDIMF.EQ.1) SFAC=0.5D0
!     == WITH SPACE-GROUP SYMMETRY: N ONLY AT THE IRREDUCIBLE K-POINTS ========
      TSTARF=(NOPF.GT.1)
      SMULT=1.D0
      LMXCX=(LMAXCX+1)**2
      NGLQX=MAXVAL(QS(:)%NGL)
      ALLOCATE(Z(NGLQX,LMXCX,NATF))
      ALLOCATE(RHOR(NRLF,NBLK),RHOG(NGLQX,NBLK),VGL(NGLQX,NBLK))
      ALLOCATE(QL(NBLK,LMXCX,NATF))
      ALLOCATE(EG0(NRLF),AUX(NGLQX))
      ALLOCATE(GIJ(6,NGLQX),DVK(NGLQX),GABS(NGLQX))
      ALLOCATE(MB(NBLK),G0B(3,NBLK),LAMB(NBLK))
      IF(TSTRESS) ALLOCATE(SYLM(NGLQX,LMXCX,6))
      EPAIR=0.D0
      SKER(:)=0.D0
      SSHP(:)=0.D0
      SV0=0.D0
      NPAIR=0
      DO IQ=1,NQF
        CALL FOCK_SELECTQ(IQ)
!       == K_I K_J/K^2, -K DV/DK, |K| OF THIS SET ===============================
        DO IG=1,NGLF
          SVAR=SUM(QS(IQ)%GVEC(:,IG)**2)
          GABS(IG)=SQRT(SVAR)
          IJ=0
          DO I=1,3
            DO J=I,3
              IJ=IJ+1
              GIJ(IJ,IG)=QS(IQ)%GVEC(I,IG)*QS(IQ)%GVEC(J,IG)/(SVAR+RSMALL)
            ENDDO
          ENDDO
          IF(SVAR.LT.1.D-10) THEN
            DVK(IG)=0.D0
          ELSE IF(MU.GT.0.D0) THEN
            DVK(IG)=2.D0*VKER(IG)-2.D0*PI/MU**2*EXP(-SVAR/(4.D0*MU**2))
          ELSE
            DVK(IG)=2.D0*VKER(IG)
          END IF
        ENDDO
        Z(:,:,:)=(0.D0,0.D0)
        DO ISPIN=1,NSPINF
          DO N=1,NOCC(ISPIN)
            IF(MOD(N-1,NTX).NE.ITX-1) CYCLE
            IF(TSTARF) THEN
!             == ONLY ORBITALS AT THE LISTED K-POINTS (STAR MEMBER 1), WITH  ==
!             == ALL PARTNERS M AND THE STAR SIZE AS WEIGHT (FORCES AND      ==
!             == STRESS ARE SYMMETRIZED BY THE CALLER)                       ==
              IF(IQOCC(N,ISPIN).GT.NKSETF) CYCLE
              SMULT=REAL(NMEMF(IQOCC(N,ISPIN)),KIND=8)
              MLAST=NOCC(ISPIN)
            ELSE
              MLAST=N
            END IF
            XK=XQFULL(:,IQOCC(N,ISPIN))
            NB1=0
            DO M=1,MLAST
              CALL FOCK_PAIRSET(XK,XQFULL(:,IQOCC(M,ISPIN)),IQM,G0)
              IF(IQM.NE.IQ) CYCLE
              NB1=NB1+1
              MB(NB1)=M
              G0B(:,NB1)=G0
              LAMB(NB1)=-0.5D0*AEXX*WOCC(M,ISPIN)*WOCC(N,ISPIN)/SFAC
              IF(TSTARF) THEN
                LAMB(NB1)=SMULT*LAMB(NB1)
              ELSE IF(M.LT.N) THEN
                LAMB(NB1)=2.D0*LAMB(NB1)
              END IF
              IF(NB1.EQ.NBLK) THEN
                CALL PAIRBLOCK()
                NB1=0
              END IF
            ENDDO
            IF(NB1.GT.0) CALL PAIRBLOCK()
          ENDDO
        ENDDO
!
!       == FORCES AND SHAPE STRESS FROM Z OF THIS SET ===========================
        IF(TSTRESS) CALL WAVES_STRAINEDYLM(NGLF,LMXCX,QS(IQ)%GVEC,YLMGF &
     &                                    ,SYLM(1:NGLF,:,:))
        DO IAT=1,NATF
          ISP=ISPF(IAT)
          LMXC=SPF(ISP)%LMXC
          FR(:)=0.D0
          DO LM=1,LMXC
            L=INT(SQRT(REAL(LM-1,KIND=8)+1.D-6))
            CFAC=(-CI)**L
            AUX(1:NGLF)=CFAC*EIGRF(:,IAT)*Z(1:NGLF,LM,IAT)
            DO IG=1,NGLF
              CSVAR=AUX(IG)*FGQ(IG,L,ISP)*YLMGF(IG,LM)
!             == -2 RE[(-I K) CSVAR] = -2 K IM(CSVAR) ===========================
              FR(:)=FR(:)-2.D0*QS(IQ)%GVEC(:,IG)*AIMAG(CSVAR)
            ENDDO
            IF(TSTRESS) THEN
              DO IG=1,NGLF
                SVAR=REAL(AUX(IG),KIND=8)
                C=QS(IQ)%FGD(IG,L,ISP)*GABS(IG)*YLMGF(IG,LM)
                DO IJ=1,6
                  SSHP(IJ)=SSHP(IJ)-2.D0*SVAR*(C*GIJ(IJ,IG) &
     &                               +FGQ(IG,L,ISP)*SYLM(IG,LM,IJ))
                ENDDO
              ENDDO
            END IF
          ENDDO
          FORCE(:,IAT)=FORCE(:,IAT)+FR(:)
        ENDDO
      ENDDO
!
!     == SUMS OVER THE G VECTORS ('K') AND THE PAIRS ('KX') ===================
      CALL MPE$COMBINE('MONOMER','+',FORCE)
      CALL MPE$COMBINE('MONOMER','+',EPAIR)
      CALL MPE$COMBINE('MONOMER','+',NPAIR)
      IF(TSTRESS) THEN
        CALL MPE$COMBINE('MONOMER','+',SKER)
        CALL MPE$COMBINE('MONOMER','+',SSHP)
        CALL MPE$COMBINE('MONOMER','+',SV0)
!       == DV0/DEPS FOR THE STRAINED BORN-VON KARMAN SUPERCELL (NUMERIC) ======
        DO I=1,3
          DO J=I,3
            DEPS(:,:)=0.D0
            DEPS(I,J)=0.5D0*HEPS
            DEPS(J,I)=DEPS(J,I)+0.5D0*HEPS
            DO K=1,3
              RBASS(:,K)=(RBASF(:,K)+MATMUL(DEPS,RBASF(:,K)))*REAL(NKDIVF(K),KIND=8)
            ENDDO
            CALL FOCK_V0(RBASS,ALPHAF,VP)
            DO K=1,3
              RBASS(:,K)=(RBASF(:,K)-MATMUL(DEPS,RBASF(:,K)))*REAL(NKDIVF(K),KIND=8)
            ENDDO
            CALL FOCK_V0(RBASS,ALPHAF,VM)
            DV0(I,J)=(VP-VM)/(2.D0*HEPS)
            DV0(J,I)=DV0(I,J)
          ENDDO
        ENDDO
        IJ=0
        DO I=1,3
          DO J=I,3
            IJ=IJ+1
            STRESS(I,J)=SKER(IJ)+SSHP(IJ)+SV0*DV0(I,J)
            IF(I.EQ.J) STRESS(I,J)=STRESS(I,J)-EPAIR
            STRESS(J,I)=STRESS(I,J)
          ENDDO
        ENDDO
      END IF
!
!     == DIAGNOSTIC: PAIR ENERGY VS THE ENERGY OF FOCK$WBUILD ==================
      CALL GET_ENVIRONMENT_VARIABLE('VESPAW_FOCK_CHECK',CENV)
      IF(TRIM(CENV).EQ.'1') THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(*,FMT='("FOCK FORCE CHECK: E(PAIRS)=",ES22.14," E_X=",ES22.14' &
     &         //'," DIFF=",ES10.2," #PAIRS=",I8)')EPAIR,EXF,EPAIR-EXF,NPAIR
        WRITE(NFILO,FMT='("FOCK FORCE CHECK: E(PAIRS)=",ES22.14," E_X=",ES22.14' &
     &         //'," DIFF=",ES10.2," #PAIRS=",I8)')EPAIR,EXF,EPAIR-EXF,NPAIR
        IF(TSTRESS) THEN
          WRITE(*,FMT='("FOCK STRESS KERNEL ",6ES14.6)')SKER
          WRITE(*,FMT='("FOCK STRESS SHAPES ",6ES14.6)')SSHP
          WRITE(*,FMT='("FOCK STRESS G=0    ",ES14.6," DV0",6ES14.6)')SV0 &
     &          ,DV0(1,1),DV0(1,2),DV0(1,3),DV0(2,2),DV0(2,3),DV0(3,3)
          CALL FOCK_SYLMCHECK()
        END IF
      END IF
      DEALLOCATE(Z,RHOR,RHOG,VGL,QL,EG0,AUX,GIJ,DVK,GABS,MB,G0B,LAMB)
      IF(TSTRESS) DEALLOCATE(SYLM)
                              CALL TIMING$CLOCKOFF('FOCK$FORCE')
                              CALL TRACE$POP()
      RETURN
      CONTAINS
!       ........................................................................
        SUBROUTINE PAIRBLOCK()
!       == NB1 PAIRS (MB(K),N) OF THE SET IQ: PAIR DENSITIES, ENERGY, KERNEL  ==
!       == STRESS AND ACCUMULATION OF Z                                       ==
        INTEGER(4) :: KK,IG1,IS1,I1
        REAL(8)    :: RHO2,CL
        DO KK=1,NB1
          RHOR(:,KK)=CONJG(PHIR(1:NRLF,MB(KK),ISPIN))*PHIR(1:NRLF,N,ISPIN)
          DO IS1=2,NDIMF
            I1=(IS1-1)*NRLF
            RHOR(:,KK)=RHOR(:,KK)+CONJG(PHIR(I1+1:I1+NRLF,MB(KK),ISPIN)) &
     &                           *PHIR(I1+1:I1+NRLF,N,ISPIN)
          ENDDO
          TPHASE=(MAXVAL(ABS(G0B(:,KK))).GT.0.5D0)
          IF(TPHASE) THEN
            CALL PLANEWAVE_EIKR(NR1F,NR2F,NR3F,NR1STARTF,NR1STARTF+NR1LF-1 &
     &                         ,G0B(:,KK),EG0)
            RHOR(:,KK)=RHOR(:,KK)*EG0(:)
          END IF
        ENDDO
        CALL PLANEWAVE$FFT('RTOG',NB1,NGLF,RHOG(1:NGLF,1:NB1),NRLF,RHOR(:,1:NB1))
        CALL FOCK_PAIRQ(NB1,MB(1:NB1),ISPIN,N,QL(1:NB1,:,:))
        CALL FOCK_ADDQL(NB1,QL(1:NB1,:,:),RHOG(1:NGLF,1:NB1))
        DO KK=1,NB1
          CL=LAMB(KK)*GWEIGHTF
          DO IG1=1,NGLF
            RHO2=REAL(CONJG(RHOG(IG1,KK))*RHOG(IG1,KK),KIND=8)
            EPAIR=EPAIR+CL*VKER(IG1)*RHO2
            IF(TSTRESS) THEN
              IF(GABS(IG1).LT.1.D-5) THEN
                SV0=SV0+CL*RHO2
              ELSE
                SKER(:)=SKER(:)+CL*DVK(IG1)*RHO2*GIJ(:,IG1)
              END IF
            END IF
            VGL(IG1,KK)=LAMB(KK)*CONJG(VKER(IG1)*RHOG(IG1,KK))
          ENDDO
        ENDDO
        NPAIR=NPAIR+NB1
!       == Z(G,(L,R)) += SUM_K VGL(G,K) QL(K,(L,R)) ============================
        CALL ZGEMM('N','N',NGLF,LMXCX*NATF,NB1,(1.D0,0.D0),VGL,NGLQX &
     &            ,QL,NBLK,(1.D0,0.D0),Z,NGLQX)
        END SUBROUTINE PAIRBLOCK
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_PAIRSET(XK,XQM,IQ,G0)
!     **************************************************************************
!     **  PAIR-DENSITY SET IQ AND FOLDING VECTOR G0 (RELATIVE) WITH           **
!     **  XK-XQM = Q(IQ) + G0  (SEE FOCK_APPLYEXACT)                          **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: XK(3),XQM(3)
      INTEGER(4),INTENT(OUT) :: IQ
      REAL(8)   ,INTENT(OUT) :: G0(3)
      REAL(8)                :: D(3),XQ(3)
      INTEGER(4)             :: I,JQ
!     **************************************************************************
      D=XK-XQM
      XQ=MODULO(D,1.D0)
      DO I=1,3
        IF(ABS(XQ(I)-1.D0).LT.1.D-8) XQ(I)=0.D0
      ENDDO
      G0=ANINT(D-XQ)
      IQ=0
      DO JQ=1,NQF
        D=XQ-QS(JQ)%XQ
        IF(MAXVAL(ABS(D-ANINT(D))).LT.1.D-8) THEN
          IQ=JQ
          G0=G0+ANINT(D)
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('NO PAIR-DENSITY SET FOR K-Q')
      CALL ERROR$STOP('FOCK_PAIRSET')
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_PAIRQ(NV,MB,ISPIN,N,QL)
!     **************************************************************************
!     **  MULTIPOLES OF THE PAIRS (MB(K),N) OF THE SNAPSHOT:                  **
!     **    QL(K,L,R) = SUM_S SUM_IJ CONJG(<P_I|M_KS>) <P_J|N_S> T_IJ^L        **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NV
      INTEGER(4),INTENT(IN)  :: MB(NV)
      INTEGER(4),INTENT(IN)  :: ISPIN,N
      COMPLEX(8),INTENT(OUT) :: QL(NV,(LMAXCX+1)**2,NATF)
      COMPLEX(8),ALLOCATABLE :: TJ(:,:)
      INTEGER(4)             :: IAT,ISP,LMNX,LMXC,I0,I,J,LM,K,IS,IC
!     **************************************************************************
      QL(:,:,:)=(0.D0,0.D0)
      DO IAT=1,NATF
        ISP=ISPF(IAT)
        LMNX=SPF(ISP)%LMNX
        LMXC=SPF(ISP)%LMXC
        ALLOCATE(TJ(LMNX,LMXC))
        DO IS=1,NDIMF
          I0=IPRO0F(IAT)+NPROF*(IS-1)
!         == TJ(I,L) = SUM_J T_IJ^L <P_J|N> ==
          TJ(:,:)=(0.D0,0.D0)
          DO LM=1,LMXC
            DO J=1,LMNX
              DO I=1,LMNX
                TJ(I,LM)=TJ(I,LM)+SPF(ISP)%T(I,J,LM)*CPROJ(I0+J,N,ISPIN)
              ENDDO
            ENDDO
          ENDDO
          DO LM=1,LMXC
            DO K=1,NV
              QL(K,LM,IAT)=QL(K,LM,IAT) &
     &                    +SUM(CONJG(CPROJ(I0+1:I0+LMNX,MB(K),ISPIN))*TJ(:,LM))
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(TJ)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ADDQL(NV,QL,RHOG)
!     **************************************************************************
!     **  RHOG_K(G) += 1/OMEGA SUM_{R,L} QL(K,L,R) EXP(-IGR) (-I)^L Y_L F_L(G)**
!     **  (AS FOCK_ADDCOMP, WITH THE MULTIPOLES GIVEN)                        **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NV
      COMPLEX(8),INTENT(IN)    :: QL(NV,(LMAXCX+1)**2,NATF)
      COMPLEX(8),INTENT(INOUT) :: RHOG(NGLF,NV)
      COMPLEX(8),ALLOCATABLE   :: AUX(:,:)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)
      COMPLEX(8)               :: CFAC
      INTEGER(4)               :: IAT,ISP,LMXC,LM,L,K
!     **************************************************************************
      ALLOCATE(AUX(NGLF,NV))
      DO IAT=1,NATF
        ISP=ISPF(IAT)
        LMXC=SPF(ISP)%LMXC
        AUX=(0.D0,0.D0)
        DO LM=1,LMXC
          L=INT(SQRT(REAL(LM-1,KIND=8)+1.D-6))
          CFAC=(-CI)**L
          DO K=1,NV
            IF(QL(K,LM,IAT).EQ.(0.D0,0.D0)) CYCLE
            AUX(:,K)=AUX(:,K)+CFAC*QL(K,LM,IAT)*YLMGF(:,LM)*FGQ(:,L,ISP)
          ENDDO
        ENDDO
        DO K=1,NV
          RHOG(:,K)=RHOG(:,K)+EIGRF(:,IAT)*AUX(:,K)/GWEIGHTF
        ENDDO
      ENDDO
      DEALLOCATE(AUX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_SYLMCHECK()
!     **************************************************************************
!     **  DIAGNOSTIC: STRAINED SPHERICAL HARMONICS OF WAVES_STRAINEDYLM VS A  **
!     **  FINITE DIFFERENCE  -D Y_L((1-EPS)K)/DEPS_IJ = SYM(K_I D_J Y_L)      **
!     **  AND D FG/D|K| VS A FINITE DIFFERENCE, FOR THE FIRST SET             **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      REAL(8)   ,ALLOCATABLE :: SYLM(:,:,:),YP(:),YM(:)
      REAL(8)                :: E(3,3),KP(3),KM(3),H,DEV,G,SCALE
      REAL(8)   ,ALLOCATABLE :: GA(:),FA(:,:),FDA(:,:)
      INTEGER(4)             :: LMXCX,NG,IG,I,J,IJ,LM,NGT,ISP,L,LMAXC
!     **************************************************************************
      LMXCX=(LMAXCX+1)**2
      NG=QS(1)%NGL
      NGT=MIN(NG,200)
      ALLOCATE(SYLM(NG,LMXCX,6),YP(LMXCX),YM(LMXCX))
      CALL WAVES_STRAINEDYLM(NG,LMXCX,QS(1)%GVEC,QS(1)%YLM,SYLM)
      H=1.D-5
      DEV=0.D0
      DO IG=1,NGT
        G=SQRT(SUM(QS(1)%GVEC(:,IG)**2))
        IF(G.LT.1.D-8) CYCLE
        IJ=0
        DO I=1,3
          DO J=I,3
            IJ=IJ+1
            E(:,:)=0.D0
            E(I,J)=0.5D0*H
            E(J,I)=E(J,I)+0.5D0*H
            KP=QS(1)%GVEC(:,IG)-MATMUL(E,QS(1)%GVEC(:,IG))
            KM=QS(1)%GVEC(:,IG)+MATMUL(E,QS(1)%GVEC(:,IG))
            CALL SPHERICAL$YLM(LMXCX,KP/SQRT(SUM(KP**2)),YP)
            CALL SPHERICAL$YLM(LMXCX,KM/SQRT(SUM(KM**2)),YM)
            DO LM=1,LMXCX
              DEV=MAX(DEV,ABS(-(YP(LM)-YM(LM))/(2.D0*H)-SYLM(IG,LM,IJ)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(*,FMT='("FOCK SYLM CHECK: MAX |SYLM-FD|=",ES10.2," (",I4," G-VECTORS)")') &
     &      DEV,NGT
!     == D FG/D|K| VS A FINITE DIFFERENCE OF FG (FIRST INITIALIZED SPECIES) ====
      DO ISP=1,NSPF
        IF(SPF(ISP)%TINI) EXIT
      ENDDO
      LMAXC=SPF(ISP)%LMAXC
      ALLOCATE(GA(3),FA(3,0:LMAXC),FDA(3,0:LMAXC))
      DEV=0.D0
      SCALE=0.D0
      H=1.D-4
      DO IG=1,NGT
        G=SQRT(SUM(QS(1)%GVEC(:,IG)**2))
        GA=(/G-H,G+H,G/)
        GA(1)=ABS(GA(1))
        CALL FOCK_SHAPEFT(SPF(ISP)%GID,SPF(ISP)%NR,SPF(ISP)%R,LMAXC,SPF(ISP)%GHAT &
     &                   ,3,GA,FA,FDA)
        DO L=0,LMAXC
          SCALE=MAX(SCALE,ABS(FDA(3,L)))
          IF(G.LT.2.D0*H) CYCLE
          DEV=MAX(DEV,ABS((FA(2,L)-FA(1,L))/(2.D0*H)-FDA(3,L)))
        ENDDO
      ENDDO
      WRITE(*,FMT='("FOCK SHAPE CHECK: MAX |DFG/DK-FD|=",ES10.2," MAX |DFG/DK|=",ES10.2)') &
     &      DEV,SCALE
      DEALLOCATE(GA,FA,FDA)
      DEALLOCATE(SYLM,YP,YM)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ADDCOMP(NV,ND,CM,PJ,RHOG)
!     **************************************************************************
!     **  ADDS THE COMPENSATION CHARGES OF THE PAIR DENSITIES                 **
!     **  RHO_K = SUM_S PHI_MS^* PSI_KS, K=1..NV, IN G-SPACE:                 **
!     **    RHOG_K(G) += 1/OMEGA SUM_{R,L} Q_K^{RL} EXP(-IGR) (-I)^L Y_L F_L(G)**
!     **    Q_K^{RL} = SUM_S SUM_IJ CONJG(CM_IS) PJ_KJS T_IJ^L                 **
!     **  ND=2: SPINOR COMPONENTS S (CM INDEX I+NPROF*(S-1))                  **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NV
      INTEGER(4),INTENT(IN)    :: ND
      COMPLEX(8),INTENT(IN)    :: CM(NPROF*ND)       ! <P|PHI_M>
      COMPLEX(8),INTENT(IN)    :: PJ(NV,NPROF,ND)    ! <P|PSI_K>
      COMPLEX(8),INTENT(INOUT) :: RHOG(NGLF,NV)
      COMPLEX(8),ALLOCATABLE   :: MJL(:,:),QKL(:,:),AUX(:,:)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)
      COMPLEX(8)               :: CFAC
      INTEGER(4)               :: IAT,ISP,LMNX,LMXC,I0,I,J,LM,L,K,IS,IC
!     **************************************************************************
      DO IAT=1,NATF
        ISP=ISPF(IAT)
        LMNX=SPF(ISP)%LMNX
        LMXC=SPF(ISP)%LMXC
        I0=IPRO0F(IAT)
        ALLOCATE(MJL(LMNX,LMXC),QKL(NV,LMXC),AUX(NGLF,NV))
        QKL=(0.D0,0.D0)
        DO IS=1,ND
          IC=I0+NPROF*(IS-1)
!         == MJL(J,L)=SUM_I CONJG(CM_I) T_IJ^L ==
          MJL=(0.D0,0.D0)
          DO LM=1,LMXC
            DO J=1,LMNX
              DO I=1,LMNX
                MJL(J,LM)=MJL(J,LM)+CONJG(CM(IC+I))*SPF(ISP)%T(I,J,LM)
              ENDDO
            ENDDO
          ENDDO
          QKL=QKL+MATMUL(PJ(:,I0+1:I0+LMNX,IS),MJL)
        ENDDO
        AUX=(0.D0,0.D0)
        DO LM=1,LMXC
          L=INT(SQRT(REAL(LM-1,KIND=8)+1.D-6))
          CFAC=(-CI)**L
          DO K=1,NV
            IF(QKL(K,LM).EQ.(0.D0,0.D0)) CYCLE
            AUX(:,K)=AUX(:,K)+CFAC*QKL(K,LM)*YLMGF(:,LM)*FGQ(:,L,ISP)
          ENDDO
        ENDDO
        DO K=1,NV
          RHOG(:,K)=RHOG(:,K)+EIGRF(:,IAT)*AUX(:,K)/GWEIGHTF
        ENDDO
        DEALLOCATE(MJL,QKL,AUX)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_COMPPROJ(NV,ND,CM,WG,HPROJ)
!     **************************************************************************
!     **  BACK-PROJECTION OF THE POTENTIALS W_K OF THE PAIR DENSITIES ONTO    **
!     **  THE COMPENSATION CHARGES:                                          **
!     **    X_K^{RL} = INT W_K(R) GHAT_L(|R-R_R|) Y_L D^3R                     **
!     **             = SUM_G W_K(G) EXP(+IGR) (+I)^L Y_L F_L(G)               **
!     **    HPROJ_KJ += SUM_{L,I} CM_I T_IJ^L X_K^{RL}                        **
!     **  (HPROJ IS NOT YET MULTIPLIED BY THE OCCUPATION WEIGHT)             **
!     **  ND=2: SPINOR COMPONENT S OF HPROJ FROM COMPONENT S OF CM            **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NV
      INTEGER(4),INTENT(IN)    :: ND
      COMPLEX(8),INTENT(IN)    :: CM(NPROF*ND)
      COMPLEX(8),INTENT(IN)    :: WG(NGLF,NV)
      COMPLEX(8),INTENT(INOUT) :: HPROJ(NV,NPROF,ND)
      COMPLEX(8),ALLOCATABLE   :: X(:,:,:),AUX(:)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)
      INTEGER(4)               :: IAT,ISP,LMNX,LMXC,I0,I,J,LM,L,K,LMXCX,IS,IC
!     **************************************************************************
      LMXCX=(LMAXCX+1)**2
      ALLOCATE(X(NV,LMXCX,NATF))
      ALLOCATE(AUX(NGLF))
      X=(0.D0,0.D0)
      DO IAT=1,NATF
        ISP=ISPF(IAT)
        LMXC=SPF(ISP)%LMXC
        DO K=1,NV
          AUX(:)=WG(:,K)*CONJG(EIGRF(:,IAT))
          DO LM=1,LMXC
            L=INT(SQRT(REAL(LM-1,KIND=8)+1.D-6))
            X(K,LM,IAT)=CI**L*SUM(AUX(:)*YLMGF(:,LM)*FGQ(:,L,ISP))
          ENDDO
        ENDDO
      ENDDO
      CALL MPE$COMBINE('K','+',X)
      DO IAT=1,NATF
        ISP=ISPF(IAT)
        LMNX=SPF(ISP)%LMNX
        LMXC=SPF(ISP)%LMXC
        I0=IPRO0F(IAT)
        DO IS=1,ND
          IC=I0+NPROF*(IS-1)
          DO LM=1,LMXC
            DO J=1,LMNX
              DO I=1,LMNX
                IF(SPF(ISP)%T(I,J,LM).EQ.0.D0) CYCLE
                DO K=1,NV
                  HPROJ(K,I0+J,IS)=HPROJ(K,I0+J,IS) &
     &                     +CM(IC+I)*SPF(ISP)%T(I,J,LM)*X(K,LM,IAT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(X,AUX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$ENERGY(EX)
!     **************************************************************************
!     **  PLANE-WAVE PART OF THE EXCHANGE ENERGY OF THE SNAPSHOT              **
!     **    E_X = -AEXX/2 SUM_S SUM_MN W_M W_N/S K_MN                         **
!     **    K_MN = OMEGA SUM_G V(G) |RHO_MN(G)|^2,                            **
!     **    RHO_MN = PHI_M^* PHI_N + RHOHAT_MN                                **
!     **  W/S = F (THE OCCUPATION). K_MN = K_NM, SO ONLY M<=N IS EVALUATED.  **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: EX
      COMPLEX(8),ALLOCATABLE :: RHOR(:,:)
      COMPLEX(8),ALLOCATABLE :: RHOG(:,:)
      COMPLEX(8),ALLOCATABLE :: PJ(:,:)
      REAL(8)                :: SFAC,SVAR,FAC
      INTEGER(4)             :: ISPIN,M,N0,NB1,K,N,IG,IS,I1
!     **************************************************************************
                              CALL TRACE$PUSH('FOCK$ENERGY')
                              CALL TIMING$CLOCKON('FOCK$ENERGY')
      SFAC=1.D0
      IF(NSPINF.EQ.1.AND.NDIMF.EQ.1) SFAC=0.5D0
      EX=0.D0
      IF(NQF.NE.1) THEN
!       == PAIR FORMULA ONLY FOR GAMMA (DIAGNOSTIC); THE ENERGY IS TAKEN FROM ==
!       == 1/2 SUM F <PSI|V_X|PSI> IN FOCK$WBUILD                             ==
                              CALL TIMING$CLOCKOFF('FOCK$ENERGY')
                              CALL TRACE$POP()
        RETURN
      END IF
      CALL FOCK_SELECTQ(1)
      ALLOCATE(RHOR(NRLF,NBLK))
      ALLOCATE(RHOG(NGLF,NBLK))
      ALLOCATE(PJ(NBLK,NPROD))
      EX=0.D0
      DO ISPIN=1,NSPINF
        DO M=1,NOCC(ISPIN)
          DO N0=M,NOCC(ISPIN),NBLK
            NB1=MIN(NBLK,NOCC(ISPIN)-N0+1)
            DO K=1,NB1
              N=N0+K-1
              RHOR(:,K)=(0.D0,0.D0)
              DO IS=1,NDIMF
                I1=(IS-1)*NRLF
                RHOR(:,K)=RHOR(:,K)+CONJG(PHIR(I1+1:I1+NRLF,M,ISPIN)) &
     &                             *PHIR(I1+1:I1+NRLF,N,ISPIN)
              ENDDO
              PJ(K,:)=CPROJ(:,N,ISPIN)
            ENDDO
            CALL PLANEWAVE$FFT('RTOG',NB1,NGLF,RHOG(:,1:NB1),NRLF,RHOR(:,1:NB1))
            CALL FOCK_ADDCOMP(NB1,NDIMF,CPROJ(:,M,ISPIN),PJ(1:NB1,:),RHOG(:,1:NB1))
            DO K=1,NB1
              N=N0+K-1
              SVAR=0.D0
              DO IG=1,NGLF
                SVAR=SVAR+VKER(IG)*REAL(CONJG(RHOG(IG,K))*RHOG(IG,K),KIND=8)
              ENDDO
              FAC=2.D0
              IF(N.EQ.M) FAC=1.D0
              EX=EX-FAC*WOCC(M,ISPIN)*WOCC(N,ISPIN)/SFAC*SVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      EX=0.5D0*AEXX*GWEIGHTF*EX
      CALL MPE$COMBINE('K','+',EX)
      DEALLOCATE(RHOR)
      DEALLOCATE(RHOG)
      DEALLOCATE(PJ)
                              CALL TIMING$CLOCKOFF('FOCK$ENERGY')
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_APPLYEXACT(GSETID,ISPIN,NGL,NDIM,NBH,NPRO,PSI,PROJ,HPSI,HPROJ)
!     **************************************************************************
!     **  EXACT OPERATOR (2*NOCC FFTS PER BAND); FOCK$APPLY USES ACE ONCE BUILT**
!     **  HPSI = HPSI + V_X PSI (PLANE-WAVE PART) AND HPROJ = HPROJ + THE     **
!     **  PROJECTOR PART FROM THE COMPENSATION CHARGES, FOR NBH (POSSIBLY     **
!     **  TINV-PACKED) VECTORS OF THE SET GSETID (K-POINT K). THE CALLER ADDS **
!     **  SUM_J |P_J> HPROJ_J (WAVES_ADDPRO).                                 **
!     **  FOR AN OCCUPIED ORBITAL M AT Q': K-Q' = Q + G0 (Q IN A FOCK SET);   **
!     **  PAIR DENSITY (PERIODIC PART FOR Q) CONJG(U_M) U_K EXP(I G0 R),      **
!     **  POTENTIAL P_Q FROM THE SET Q, V_X U_K += -W U_M P_Q EXP(-I G0 R).   **
!     **  AT TINV-PACKED SETS V_X IS REAL IN THE BLOCH SENSE AND ACTS ON THE  **
!     **  TWO PACKED ORBITALS INDEPENDENTLY.                                  **
!     **  NDIM=2 (SPINORS): THE PAIR DENSITY SUMS OVER THE COMPONENTS, AND    **
!     **  EACH COMPONENT S RECEIVES -W U_MS P_Q.                              **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: GSETID
      INTEGER(4)  ,INTENT(IN)    :: ISPIN
      INTEGER(4)  ,INTENT(IN)    :: NGL
      INTEGER(4)  ,INTENT(IN)    :: NDIM
      INTEGER(4)  ,INTENT(IN)    :: NBH
      INTEGER(4)  ,INTENT(IN)    :: NPRO
      COMPLEX(8)  ,INTENT(IN)    :: PSI(NGL,NDIM,NBH)
      COMPLEX(8)  ,INTENT(IN)    :: PROJ(NDIM,NBH,NPRO)
      COMPLEX(8)  ,INTENT(INOUT) :: HPSI(NGL,NDIM,NBH)
      COMPLEX(8)  ,INTENT(INOUT) :: HPROJ(NDIM,NBH,NPRO)
      COMPLEX(8)  ,ALLOCATABLE   :: PSIR(:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: ACC(:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: RHOR(:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: RHOG(:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: FG(:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: FRS(:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: HPJ(:,:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: PJB(:,:,:)
      COMPLEX(8)  ,ALLOCATABLE   :: EG0(:)
      INTEGER(4)  ,ALLOCATABLE   :: IQOFM(:)
      REAL(8)     ,ALLOCATABLE   :: G0OFM(:,:)
      INTEGER(4)                 :: IB0,NB1,K,M,NGLQX,IQ,I,IS,I1
      REAL(8)                    :: W,XK(3),D(3),XQ(3),G0(3),G0LAST(3)
      LOGICAL(4)                 :: TPHASE
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(ISPIN.GT.NSPINF) RETURN
      IF(NOCC(ISPIN).EQ.0) RETURN
      IF(NPRO.NE.NPROF.OR.NDIM.NE.NDIMF) THEN
        CALL ERROR$MSG('INCONSISTENT NUMBER OF PROJECTORS OR SPINOR COMPONENTS')
        CALL ERROR$I4VAL('NPRO',NPRO)
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$STOP('FOCK_APPLYEXACT')
      END IF
                              CALL TRACE$PUSH('FOCK_APPLYEXACT')
                              CALL TIMING$CLOCKON('FOCK_APPLYEXACT')
      CALL PLANEWAVE$SELECT(GSETID)
      CALL PLANEWAVE$GETR8A('XK',3,XK)
!     == SET Q AND FOLDING VECTOR G0 FOR EVERY OCCUPIED ORBITAL ===============
      ALLOCATE(IQOFM(NOCC(ISPIN)),G0OFM(3,NOCC(ISPIN)))
      DO M=1,NOCC(ISPIN)
        D=XK-XQFULL(:,IQOCC(M,ISPIN))
        XQ=MODULO(D,1.D0)
        DO I=1,3
          IF(ABS(XQ(I)-1.D0).LT.1.D-8) XQ(I)=0.D0
        ENDDO
        G0OFM(:,M)=ANINT(D-XQ)
        IQOFM(M)=0
        DO IQ=1,NQF
          D=XQ-QS(IQ)%XQ
          IF(MAXVAL(ABS(D-ANINT(D))).LT.1.D-8) THEN
            IQOFM(M)=IQ
            G0OFM(:,M)=G0OFM(:,M)+ANINT(D)
            EXIT
          END IF
        ENDDO
        IF(IQOFM(M).EQ.0) THEN
          CALL ERROR$MSG('NO PAIR-DENSITY SET FOR K-Q')
          CALL ERROR$STOP('FOCK_APPLYEXACT')
        END IF
      ENDDO
      NGLQX=MAXVAL(QS(:)%NGL)
      ALLOCATE(PSIR(NRLD,NBLK))
      ALLOCATE(ACC(NRLD,NBLK))
      ALLOCATE(RHOR(NRLF,NBLK))
      ALLOCATE(RHOG(NGLQX,NBLK))
      ALLOCATE(FG(NGL,NBLK))
      ALLOCATE(FRS(NRLF,NBLK))
      ALLOCATE(HPJ(NBLK,NPROF,NDIM))
      ALLOCATE(PJB(NBLK,NPROF,NDIM))
      ALLOCATE(EG0(NRLF))
      DO IB0=1,NBH,NBLK
        NB1=MIN(NBLK,NBH-IB0+1)
        CALL PLANEWAVE$SELECT(GSETID)
        DO IS=1,NDIM
          I1=(IS-1)*NRLF
          FG(:,1:NB1)=PSI(:,IS,IB0:IB0+NB1-1)
          CALL PLANEWAVE$FFT('GTOR',NB1,NGL,FG(:,1:NB1),NRLF,FRS(:,1:NB1))
          PSIR(I1+1:I1+NRLF,1:NB1)=FRS(:,1:NB1)
          PJB(1:NB1,:,IS)=PROJ(IS,IB0:IB0+NB1-1,:)
        ENDDO
        ACC(:,1:NB1)=(0.D0,0.D0)
        HPJ(1:NB1,:,:)=(0.D0,0.D0)
        G0LAST=1.D9
        DO M=1,NOCC(ISPIN)
          W=-AEXX*WOCC(M,ISPIN)
          G0=G0OFM(:,M)
          TPHASE=(MAXVAL(ABS(G0)).GT.0.5D0)
          IF(TPHASE.AND.MAXVAL(ABS(G0-G0LAST)).GT.0.5D0) THEN
            CALL PLANEWAVE_EIKR(NR1F,NR2F,NR3F,NR1STARTF,NR1STARTF+NR1LF-1,G0,EG0)
            G0LAST=G0
          END IF
          DO K=1,NB1
            RHOR(:,K)=CONJG(PHIR(1:NRLF,M,ISPIN))*PSIR(1:NRLF,K)
            DO IS=2,NDIM
              I1=(IS-1)*NRLF
              RHOR(:,K)=RHOR(:,K)+CONJG(PHIR(I1+1:I1+NRLF,M,ISPIN))*PSIR(I1+1:I1+NRLF,K)
            ENDDO
            IF(TPHASE) RHOR(:,K)=RHOR(:,K)*EG0(:)
          ENDDO
!         == PAIR DENSITY WITH COMPENSATION CHARGES -> POTENTIAL W(Q+G) =======
          CALL FOCK_SELECTQ(IQOFM(M))
          CALL PLANEWAVE$FFT('RTOG',NB1,NGLF,RHOG(1:NGLF,1:NB1),NRLF,RHOR(:,1:NB1))
          CALL FOCK_ADDCOMP(NB1,NDIM,CPROJ(:,M,ISPIN),PJB(1:NB1,:,:) &
     &                     ,RHOG(1:NGLF,1:NB1))
          DO K=1,NB1
            RHOG(1:NGLF,K)=RHOG(1:NGLF,K)*VKER(:)
          ENDDO
!         == PROJECTOR PART (FROM THE COMPENSATION CHARGES) ====================
          CALL FOCK_COMPPROJ(NB1,NDIM,W*CPROJ(:,M,ISPIN),RHOG(1:NGLF,1:NB1) &
     &                      ,HPJ(1:NB1,:,:))
!         == PLANE-WAVE PART ===================================================
          CALL PLANEWAVE$FFT('GTOR',NB1,NGLF,RHOG(1:NGLF,1:NB1),NRLF,RHOR(:,1:NB1))
          IF(TPHASE) THEN
            DO K=1,NB1
              RHOR(:,K)=RHOR(:,K)*CONJG(EG0(:))
            ENDDO
          END IF
          DO K=1,NB1
            DO IS=1,NDIM
              I1=(IS-1)*NRLF
              ACC(I1+1:I1+NRLF,K)=ACC(I1+1:I1+NRLF,K) &
     &                           +W*PHIR(I1+1:I1+NRLF,M,ISPIN)*RHOR(:,K)
            ENDDO
          ENDDO
        ENDDO
        CALL PLANEWAVE$SELECT(GSETID)
        DO IS=1,NDIM
          I1=(IS-1)*NRLF
          FRS(:,1:NB1)=ACC(I1+1:I1+NRLF,1:NB1)
          CALL PLANEWAVE$FFT('RTOG',NB1,NGL,FG(:,1:NB1),NRLF,FRS(:,1:NB1))
          HPSI(:,IS,IB0:IB0+NB1-1)=HPSI(:,IS,IB0:IB0+NB1-1)+FG(:,1:NB1)
          HPROJ(IS,IB0:IB0+NB1-1,:)=HPROJ(IS,IB0:IB0+NB1-1,:)+HPJ(1:NB1,:,IS)
        ENDDO
      ENDDO
      DEALLOCATE(PSIR,ACC,RHOR,RHOG,FG,FRS,HPJ,PJB,EG0,IQOFM,G0OFM)
      CALL PLANEWAVE$SELECT(GSETID)
                              CALL TIMING$CLOCKOFF('FOCK_APPLYEXACT')
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$APPLY(GSETID,ISPIN,NGL,NDIM,NBH,NPRO,PSI,PROJ,HPSI,HPROJ)
!     **************************************************************************
!     **  HPSI = HPSI + V_X PSI, HPROJ = HPROJ + PROJECTOR PART (SEE           **
!     **  FOCK_APPLYEXACT). WITH ACE (DEFAULT, ONCE FOCK$WBUILD RAN FOR THIS   **
!     **  K AND SPIN):  V_X PSI = -SUM_K XI_K <XI_K|PSI>,                      **
!     **    <XI|PSI> = <XI_PW|PSI>_G + SUM_I CONJG(XI_P,I) <P_I|PSI>.          **
!     **  AT TINV SETS THE XI ARE REAL BLOCH FUNCTIONS, SO THE OPERATOR ACTS   **
!     **  ON THE TWO PACKED ORBITALS INDEPENDENTLY, AS THE EXACT ONE.          **
!     **  NDIM=2: XI ARE SPINORS, PLANE-WAVE PART (NGL,NDIM) AND PROJECTOR    **
!     **  PART (NPRO,NDIM) STACKED.                                           **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: GSETID
      INTEGER(4)  ,INTENT(IN)    :: ISPIN
      INTEGER(4)  ,INTENT(IN)    :: NGL
      INTEGER(4)  ,INTENT(IN)    :: NDIM
      INTEGER(4)  ,INTENT(IN)    :: NBH
      INTEGER(4)  ,INTENT(IN)    :: NPRO
      COMPLEX(8)  ,INTENT(IN)    :: PSI(NGL,NDIM,NBH)
      COMPLEX(8)  ,INTENT(IN)    :: PROJ(NDIM,NBH,NPRO)
      COMPLEX(8)  ,INTENT(INOUT) :: HPSI(NGL,NDIM,NBH)
      COMPLEX(8)  ,INTENT(INOUT) :: HPROJ(NDIM,NBH,NPRO)
      COMPLEX(8)  ,ALLOCATABLE   :: A(:,:),PT(:,:),HT(:,:)
      INTEGER(4)                 :: N,IK,IS,I
      REAL(8)                    :: GWEIGHT,XK(3)
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(ISPIN.GT.NSPINF) RETURN
      IF(NOCC(ISPIN).EQ.0) RETURN
      N=0
      IK=0
      IF(TACE.AND.ALLOCATED(ACES)) THEN
        CALL PLANEWAVE$SELECT(GSETID)
        CALL PLANEWAVE$GETR8A('XK',3,XK)
        IK=FOCK_FINDK(XK,NKSETF,XKSETF)
        IF(IK.GT.0) THEN
          N=ACES(IK,ISPIN)%N
          IF(ACES(IK,ISPIN)%NGL.NE.NGL) N=0
        END IF
      END IF
      IF(N.EQ.0) THEN
        CALL FOCK_APPLYEXACT(GSETID,ISPIN,NGL,NDIM,NBH,NPRO,PSI,PROJ,HPSI,HPROJ)
        RETURN
      END IF
                              CALL TRACE$PUSH('FOCK$APPLY')
                              CALL TIMING$CLOCKON('FOCK$APPLY(ACE)')
      CALL PLANEWAVE$SELECT(GSETID)
      CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHT)
      ALLOCATE(A(N,NBH))
      ALLOCATE(PT(NPRO*NDIM,NBH))
      ALLOCATE(HT(NPRO*NDIM,NBH))
      DO IS=1,NDIM
        DO I=1,NPRO
          PT(I+NPRO*(IS-1),:)=PROJ(IS,:,I)
        ENDDO
      ENDDO
!     == A(K,J) = <XI_K|PSI_J> =================================================
      CALL ZGEMM('C','N',N,NBH,NGL*NDIM,(1.D0,0.D0),ACES(IK,ISPIN)%XIG,NGL*NDIM &
     &          ,PSI,NGL*NDIM,(0.D0,0.D0),A,N)
      A=A*GWEIGHT
      CALL MPE$COMBINE('K','+',A)
      CALL ZGEMM('C','N',N,NBH,NPRO*NDIM,(1.D0,0.D0),ACES(IK,ISPIN)%XIP,NPRO*NDIM &
     &          ,PT,NPRO*NDIM,(1.D0,0.D0),A,N)
!     == HPSI -= XI_PW A, HPROJ -= XI_P A ========================================
      CALL ZGEMM('N','N',NGL*NDIM,NBH,N,(-1.D0,0.D0),ACES(IK,ISPIN)%XIG,NGL*NDIM &
     &          ,A,N,(1.D0,0.D0),HPSI,NGL*NDIM)
      CALL ZGEMM('N','N',NPRO*NDIM,NBH,N,(-1.D0,0.D0),ACES(IK,ISPIN)%XIP,NPRO*NDIM &
     &          ,A,N,(0.D0,0.D0),HT,NPRO*NDIM)
      DO IS=1,NDIM
        DO I=1,NPRO
          HPROJ(IS,:,I)=HPROJ(IS,:,I)+HT(I+NPRO*(IS-1),:)
        ENDDO
      ENDDO
      DEALLOCATE(A,PT,HT)
                              CALL TIMING$CLOCKOFF('FOCK$APPLY(ACE)')
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$WBUILD(GSETID,ISPIN,NGL,NDIM,NBH,NB,NPRO,TINV,F,PSI,PROJ)
!     **************************************************************************
!     **  FOR ONE SET (K,SPIN), AFTER THE SNAPSHOT OF ALL OCCUPIED ORBITALS:   **
!     **    W = V_X PSI (EXACT) FOR ALL NB BANDS, M = <PSI|W>,                 **
!     **    ENERGY  EXF += 1/2 SUM_N F_N RE M_NN  (F WITH THE K WEIGHT),        **
!     **    ACE     -M = L L^H, XI = W L^{-H}  (IF TACE).                       **
!     **  WITH TINV THE PACKED PAIRS ARE SPLIT INTO REAL BLOCH FUNCTIONS IN    **
!     **  REAL SPACE (NO -G MAP ACROSS TASKS).                                 **
!     **************************************************************************
      USE MPE_MODULE
      USE FOCK_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: GSETID
      INTEGER(4)  ,INTENT(IN) :: ISPIN
      INTEGER(4)  ,INTENT(IN) :: NGL
      INTEGER(4)  ,INTENT(IN) :: NDIM
      INTEGER(4)  ,INTENT(IN) :: NBH
      INTEGER(4)  ,INTENT(IN) :: NB
      INTEGER(4)  ,INTENT(IN) :: NPRO
      LOGICAL(4)  ,INTENT(IN) :: TINV
      REAL(8)     ,INTENT(IN) :: F(NB)
      COMPLEX(8)  ,INTENT(IN) :: PSI(NGL,NDIM,NBH)
      COMPLEX(8)  ,INTENT(IN) :: PROJ(NDIM,NBH,NPRO)
      COMPLEX(8)  ,ALLOCATABLE:: WG(:,:,:),WP(:,:,:)  ! PACKED V_X PSI
      COMPLEX(8)  ,ALLOCATABLE:: UG(:,:),UP(:,:)      ! UNPACKED PSI (NGL*NDIM,NB),(NPRO*NDIM,NB)
      COMPLEX(8)  ,ALLOCATABLE:: VG(:,:),VP(:,:)      ! UNPACKED W
      COMPLEX(8)  ,ALLOCATABLE:: M(:,:)
      REAL(8)                 :: GWEIGHT,ASYM,XK(3)
      INTEGER(4)              :: INFO,N,IK
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(NOCC(ISPIN).EQ.0) RETURN
                              CALL TRACE$PUSH('FOCK$WBUILD')
                              CALL TIMING$CLOCKON('FOCK$WBUILD')
!     == W = V_X PSI WITH THE EXACT OPERATOR ===================================
      ALLOCATE(WG(NGL,NDIM,NBH),WP(NDIM,NBH,NPRO))
      WG=(0.D0,0.D0)
      WP=(0.D0,0.D0)
      CALL FOCK_APPLYEXACT(GSETID,ISPIN,NGL,NDIM,NBH,NPRO,PSI,PROJ,WG,WP)
!     == UNPACK INTO NB INDIVIDUAL ORBITALS ====================================
      ALLOCATE(UG(NGL*NDIM,NB),UP(NPRO*NDIM,NB),VG(NGL*NDIM,NB),VP(NPRO*NDIM,NB))
      CALL FOCK_UNPACK(GSETID,TINV,NGL,NDIM,NBH,NB,NPRO,PSI,PROJ,UG,UP)
      CALL FOCK_UNPACK(GSETID,TINV,NGL,NDIM,NBH,NB,NPRO,WG,WP,VG,VP)
      DEALLOCATE(WG,WP)
!     == M = <PSI|W> (PLANE-WAVE PART SUMMED OVER 'K', PROJECTORS REPLICATED) ===
      CALL PLANEWAVE$SELECT(GSETID)
      CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHT)
      CALL PLANEWAVE$GETR8A('XK',3,XK)
      ALLOCATE(M(NB,NB))
      CALL ZGEMM('C','N',NB,NB,NGL*NDIM,(1.D0,0.D0),UG,NGL*NDIM,VG,NGL*NDIM &
     &          ,(0.D0,0.D0),M,NB)
      M=M*GWEIGHT
      CALL MPE$COMBINE('K','+',M)
      CALL ZGEMM('C','N',NB,NB,NPRO*NDIM,(1.D0,0.D0),UP,NPRO*NDIM,VP,NPRO*NDIM &
     &          ,(1.D0,0.D0),M,NB)
!     == ENERGY 1/2 SUM_N F_N <PSI_N|V_X|PSI_N> =================================
      DO N=1,NB
        EXF=EXF+0.5D0*F(N)*REAL(M(N,N),KIND=8)
      ENDDO
      IF(TACE) THEN
!       == HERMITIAN (NEGATIVE DEFINITE) UP TO ROUNDING =========================
        ASYM=MAXVAL(ABS(M-TRANSPOSE(CONJG(M))))
        M=-0.5D0*(M+TRANSPOSE(CONJG(M)))
        CALL ZPOTRF('L',NB,M,NB,INFO)
        IF(INFO.NE.0) THEN
          CALL ERROR$MSG('ACE: -<PSI|V_X|PSI> IS NOT POSITIVE DEFINITE')
          CALL ERROR$I4VAL('ZPOTRF INFO',INFO)
          CALL ERROR$STOP('FOCK$WBUILD')
        END IF
!       == XI = W L^{-H} (PLANE-WAVE AND PROJECTOR PART) =========================
        CALL ZTRSM('R','L','C','N',NGL*NDIM,NB,(1.D0,0.D0),M,NB,VG,NGL*NDIM)
        CALL ZTRSM('R','L','C','N',NPRO*NDIM,NB,(1.D0,0.D0),M,NB,VP,NPRO*NDIM)
        IK=FOCK_FINDK(XK,NKSETF,XKSETF)
        IF(ALLOCATED(ACES(IK,ISPIN)%XIG)) DEALLOCATE(ACES(IK,ISPIN)%XIG)
        IF(ALLOCATED(ACES(IK,ISPIN)%XIP)) DEALLOCATE(ACES(IK,ISPIN)%XIP)
        ALLOCATE(ACES(IK,ISPIN)%XIG(NGL*NDIM,NB))
        ALLOCATE(ACES(IK,ISPIN)%XIP(NPRO*NDIM,NB))
        ACES(IK,ISPIN)%XIG=VG
        ACES(IK,ISPIN)%XIP=VP
        ACES(IK,ISPIN)%NGL=NGL
        ACES(IK,ISPIN)%N=NB
        IF(ASYM.GT.1.D-8) THEN
          WRITE(*,FMT='("FOCK ACE WARNING: <PSI|V_X|PSI> ASYMMETRY",ES10.2)')ASYM
        END IF
      END IF
      DEALLOCATE(UG,UP,VG,VP,M)
                              CALL TIMING$CLOCKOFF('FOCK$WBUILD')
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_UNPACK(GSETID,TINV,NGL,NDIM,NBH,NB,NPRO,F,FP,U,UP)
!     **************************************************************************
!     **  SPLITS TINV-PACKED PAIRS INTO NB INDIVIDUAL G-SPACE VECTORS: THE     **
!     **  PACKED VECTOR IS EXP(-IKR)(PHI_A + I PHI_B) WITH REAL BLOCH          **
!     **  FUNCTIONS PHI_A,PHI_B; WITHOUT TINV A PLAIN COPY.                    **
!     **  PROJECTIONS: RE/IM PARTS, FP(NDIM,NBH,NPRO) -> UP(NPRO*NDIM,NB)      **
!     **  NDIM=2 (NO TINV): SPINOR COMPONENTS STACKED, UP INDEX I+NPRO*(S-1)  **
!     **************************************************************************
      USE FOCK_MODULE, ONLY : NBLK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: GSETID
      LOGICAL(4)  ,INTENT(IN) :: TINV
      INTEGER(4)  ,INTENT(IN) :: NGL,NDIM,NBH,NB,NPRO
      COMPLEX(8)  ,INTENT(IN) :: F(NGL,NDIM,NBH)
      COMPLEX(8)  ,INTENT(IN) :: FP(NDIM,NBH,NPRO)
      COMPLEX(8)  ,INTENT(OUT):: U(NGL*NDIM,NB)
      COMPLEX(8)  ,INTENT(OUT):: UP(NPRO*NDIM,NB)
      COMPLEX(8)  ,ALLOCATABLE:: FG(:,:),FR(:,:),GR(:,:),EIKR(:)
      INTEGER(4)              :: NRL,IB0,NB1,K,IBH,IB,IS,I
!     **************************************************************************
      IF(.NOT.TINV) THEN
        DO IB=1,NB
          DO IS=1,NDIM
            U(NGL*(IS-1)+1:NGL*IS,IB)=F(:,IS,IB)
            DO I=1,NPRO
              UP(I+NPRO*(IS-1),IB)=FP(IS,IB,I)
            ENDDO
          ENDDO
        ENDDO
        RETURN
      END IF
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('TINV WITH SPINORS')
        CALL ERROR$STOP('FOCK_UNPACK')
      END IF
      CALL PLANEWAVE$SELECT(GSETID)
      CALL PLANEWAVE$GETI4('NRL',NRL)
      ALLOCATE(FG(NGL,2*NBLK),FR(NRL,NBLK),GR(NRL,2*NBLK),EIKR(NRL))
      CALL PLANEWAVE$GETC8A('EIKR',NRL,EIKR)
      DO IB0=1,NBH,NBLK
        NB1=MIN(NBLK,NBH-IB0+1)
        FG(:,1:NB1)=F(:,1,IB0:IB0+NB1-1)
        CALL PLANEWAVE$FFT('GTOR',NB1,NGL,FG(:,1:NB1),NRL,FR(:,1:NB1))
        DO K=1,NB1
          FR(:,K)=FR(:,K)*EIKR(:)
          GR(:,2*K-1)=CMPLX(REAL(FR(:,K),KIND=8),0.D0,KIND=8)*CONJG(EIKR(:))
          GR(:,2*K)  =CMPLX(AIMAG(FR(:,K)),0.D0,KIND=8)*CONJG(EIKR(:))
        ENDDO
        CALL PLANEWAVE$FFT('RTOG',2*NB1,NGL,FG(:,1:2*NB1),NRL,GR(:,1:2*NB1))
        DO K=1,NB1
          IBH=IB0+K-1
          IB=2*IBH-1
          U(:,IB)=FG(:,2*K-1)
          UP(:,IB)=CMPLX(REAL(FP(1,IBH,:),KIND=8),0.D0,KIND=8)
          IB=2*IBH
          IF(IB.LE.NB) THEN
            U(:,IB)=FG(:,2*K)
            UP(:,IB)=CMPLX(AIMAG(FP(1,IBH,:)),0.D0,KIND=8)
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(FG,FR,GR,EIKR)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK$ONECENTER(LMNXX,NDIMD,NAT,DENMAT,DH,E1)
!     **************************************************************************
!     **  ONE-CENTER CORRECTION E_X^1[AE] - E_X^1[PS+COMP] AS A FUNCTIONAL OF **
!     **  THE SPIN-RESOLVED ONE-CENTER DENSITY MATRICES D^S:                  **
!     **    E^S = -AEXX/2 SUM_IJKL D^S_IK D^S_LJ U(I,J,K,L)                    **
!     **    DE/DD^S_AB = -AEXX SUM_JL D^S_LJ U(A,J,B,L)                        **
!     **  DENMAT=(TOTAL[,SPIN]). DH(B,A) += DE/DD(A,B), E=SUM DH(I,J)D(J,I).  **
!     **  VESPAW_FOCK_CHECK=1 PRINTS E1 AND 1/2 SUM DH1 D (MUST AGREE).      **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LMNXX,NDIMD,NAT
      COMPLEX(8),INTENT(IN)    :: DENMAT(LMNXX,LMNXX,NDIMD,NAT)
      COMPLEX(8),INTENT(INOUT) :: DH(LMNXX,LMNXX,NDIMD,NAT)
      REAL(8)   ,INTENT(OUT)   :: E1
      COMPLEX(8),ALLOCATABLE   :: D(:,:,:),V(:,:,:)
      COMPLEX(8)               :: CSVAR,ECHK,CSCV
      REAL(8)                  :: ECV
      INTEGER(4)               :: IAT,ISP,LMNX,IS,NS,I,J,K,L
      LOGICAL(4),SAVE          :: TCHECKREAD=.FALSE.
      LOGICAL(4),SAVE          :: TCHECK=.FALSE.
      CHARACTER(16)            :: CENV
!     **************************************************************************
      E1=0.D0
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('FOCK$ONECENTER')
      IF(.NOT.TCHECKREAD) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_FOCK_CHECK',CENV)
        TCHECK=(TRIM(CENV).EQ.'1')
        TCHECKREAD=.TRUE.
      END IF
      IF(.NOT.ALLOCATED(SPF)) THEN
        CALL ERROR$MSG('PAW DATA NOT INITIALIZED (SNAPSHOT MISSING)')
        CALL ERROR$STOP('FOCK$ONECENTER')
      END IF
      IF(NDIMD.EQ.3.OR.NDIMD.GT.4) THEN
        CALL ERROR$MSG('NDIMD OUT OF RANGE')
        CALL ERROR$I4VAL('NDIMD',NDIMD)
        CALL ERROR$STOP('FOCK$ONECENTER')
      END IF
      ECHK=(0.D0,0.D0)
      ECV=0.D0
      DO IAT=1,NAT
        ISP=ISPF(IAT)
        LMNX=SPF(ISP)%LMNX
!       == CORE-VALENCE: LINEAR IN THE TOTAL DENSITY MATRIX ====================
        IF(TCV) THEN
          CSCV=(0.D0,0.D0)
          DO I=1,LMNX
            DO K=1,LMNX
              DH(I,K,1,IAT)=DH(I,K,1,IAT)-AEXX*SPF(ISP)%KCV(I,K)
              CSCV=CSCV-AEXX*SPF(ISP)%KCV(I,K)*DENMAT(K,I,1,IAT)
            ENDDO
          ENDDO
          ECV=ECV+REAL(CSCV,KIND=8)
        END IF
        IF(NDIMD.EQ.4) THEN
          CALL FOCK_ONECENTERNC(LMNX,ISP,DENMAT(1:LMNX,1:LMNX,:,IAT) &
     &                         ,DH(1:LMNX,1:LMNX,:,IAT),E1,ECHK)
          CYCLE
        END IF
        ALLOCATE(D(LMNX,LMNX,2),V(LMNX,LMNX,2))
        IF(NDIMD.EQ.1) THEN
          NS=1
          D(:,:,1)=0.5D0*DENMAT(1:LMNX,1:LMNX,1,IAT)
        ELSE
          NS=2
          D(:,:,1)=0.5D0*(DENMAT(1:LMNX,1:LMNX,1,IAT)+DENMAT(1:LMNX,1:LMNX,2,IAT))
          D(:,:,2)=0.5D0*(DENMAT(1:LMNX,1:LMNX,1,IAT)-DENMAT(1:LMNX,1:LMNX,2,IAT))
        END IF
        V=(0.D0,0.D0)
        DO IS=1,NS
          DO I=1,LMNX
            DO K=1,LMNX
              CSVAR=(0.D0,0.D0)
              DO J=1,LMNX
                DO L=1,LMNX
                  CSVAR=CSVAR+D(L,J,IS)*SPF(ISP)%U(I,J,K,L)
                ENDDO
              ENDDO
              V(I,K,IS)=-AEXX*CSVAR        ! DE^S/DD^S_IK
            ENDDO
          ENDDO
!         == E^S = 1/2 SUM_IK V_IK D_IK (HOMOGENEOUS OF DEGREE 2) ==============
          CSVAR=(0.D0,0.D0)
          DO I=1,LMNX
            DO K=1,LMNX
              CSVAR=CSVAR+V(I,K,IS)*D(I,K,IS)
            ENDDO
          ENDDO
          IF(NS.EQ.1) CSVAR=2.D0*CSVAR   ! BOTH SPIN CHANNELS
          E1=E1+0.5D0*REAL(CSVAR,KIND=8)
        ENDDO
!       == DH(B,A) += DE/DD_TOTAL(A,B) [AND DE/DD_SPIN(A,B)] =====================
        DO I=1,LMNX
          DO K=1,LMNX
            IF(NS.EQ.1) THEN
              DH(K,I,1,IAT)=DH(K,I,1,IAT)+V(I,K,1)
              ECHK=ECHK+V(I,K,1)*DENMAT(I,K,1,IAT)
            ELSE
              DH(K,I,1,IAT)=DH(K,I,1,IAT)+0.5D0*(V(I,K,1)+V(I,K,2))
              DH(K,I,2,IAT)=DH(K,I,2,IAT)+0.5D0*(V(I,K,1)-V(I,K,2))
              ECHK=ECHK+0.5D0*(V(I,K,1)+V(I,K,2))*DENMAT(I,K,1,IAT) &
     &                 +0.5D0*(V(I,K,1)-V(I,K,2))*DENMAT(I,K,2,IAT)
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(D,V)
      ENDDO
      EXCVF=ECV
      E1=E1+ECV
      EX1F=E1
      IF(TCHECK) THEN
        WRITE(*,FMT='("FOCK CHECK CV: E_CV=",ES22.14)')ECV
        WRITE(*,FMT='("FOCK CHECK 1C: E1=",ES22.14," 1/2SUM DH1 D=",ES22.14' &
     &              //'," DIFF=",ES10.2)')E1-ECV,0.5D0*REAL(ECHK,KIND=8) &
     &              ,E1-ECV-0.5D0*REAL(ECHK,KIND=8)
      END IF
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ONECENTERNC(LMNX,ISP,DENMAT,DH,E1,ECHK)
!     **************************************************************************
!     **  ONE-CENTER EXCHANGE FOR SPINORS (NDIMD=4, DENMAT=(TOTAL,X,Y,Z)).    **
!     **  WITH D^AB_IJ = SUM_N F_N <P_I|PSI_NA><PSI_NB|P_J>:                  **
!     **    E = -AEXX/2 SUM_AB SUM_IJKL D^BA_IK D^AB_LJ U(I,J,K,L)             **
!     **    DE/DD^AB_IK = V^AB_IK = -AEXX SUM_JL D^BA_LJ U(I,J,K,L)            **
!     **  (COLLINEAR CASE: ONLY A=B). D^AB FROM DENMAT AS LDAPLUSU_SPINDENMAT,**
!     **  DH_C(K,I) += DE/DD_C(I,K) BY THE CHAIN RULE. E1, ECHK ACCUMULATE.   **
!     **************************************************************************
      USE FOCK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LMNX,ISP
      COMPLEX(8),INTENT(IN)    :: DENMAT(LMNX,LMNX,4)
      COMPLEX(8),INTENT(INOUT) :: DH(LMNX,LMNX,4)
      REAL(8)   ,INTENT(INOUT) :: E1
      COMPLEX(8),INTENT(INOUT) :: ECHK
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)
      COMPLEX(8)               :: D(LMNX,LMNX,2,2),V(LMNX,LMNX,2,2),G(LMNX,LMNX,4)
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: IA,IB,I,J,K,L,IC
!     **************************************************************************
      D(:,:,1,1)=0.5D0*(DENMAT(:,:,1)+DENMAT(:,:,4))
      D(:,:,2,2)=0.5D0*(DENMAT(:,:,1)-DENMAT(:,:,4))
      D(:,:,1,2)=0.5D0*(DENMAT(:,:,2)-CI*DENMAT(:,:,3))
      D(:,:,2,1)=0.5D0*(DENMAT(:,:,2)+CI*DENMAT(:,:,3))
      DO IA=1,2
        DO IB=1,2
          DO I=1,LMNX
            DO K=1,LMNX
              CSVAR=(0.D0,0.D0)
              DO J=1,LMNX
                DO L=1,LMNX
                  CSVAR=CSVAR+D(L,J,IB,IA)*SPF(ISP)%U(I,J,K,L)
                ENDDO
              ENDDO
              V(I,K,IA,IB)=-AEXX*CSVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CSVAR=(0.D0,0.D0)
      DO IA=1,2
        DO IB=1,2
          CSVAR=CSVAR+SUM(V(:,:,IA,IB)*D(:,:,IA,IB))
        ENDDO
      ENDDO
      E1=E1+0.5D0*REAL(CSVAR,KIND=8)
!     == CHAIN RULE: D11=(D1+D4)/2, D22=(D1-D4)/2, D12=(D2-ID3)/2, D21=(D2+ID3)/2
      G(:,:,1)=0.5D0*(V(:,:,1,1)+V(:,:,2,2))
      G(:,:,2)=0.5D0*(V(:,:,1,2)+V(:,:,2,1))
      G(:,:,3)=0.5D0*CI*(V(:,:,2,1)-V(:,:,1,2))
      G(:,:,4)=0.5D0*(V(:,:,1,1)-V(:,:,2,2))
      DO IC=1,4
        DH(:,:,IC)=DH(:,:,IC)+TRANSPOSE(G(:,:,IC))
        ECHK=ECHK+SUM(G(:,:,IC)*DENMAT(:,:,IC))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_SHAPEFT(GID,NR,R,LMAXC,GHAT,NG,GABS,FG,FGD)
!     **************************************************************************
!     **  FG(IG,L) = 4PI INT GHAT_L(R) J_L(|G| R) R^2 DR, TABULATED ON A      **
!     **  UNIFORM |G| GRID AND INTERPOLATED WITH 4-POINT LAGRANGE POLYNOMIALS **
!     **  FGD(IG,L) = D FG/D|G| = 4PI INT GHAT_L(R) J_L'(|G| R) R^3 DR,       **
!     **  J_L' = [L J_(L-1) - (L+1) J_(L+1)]/(2L+1), TABULATED IN THE SAME WAY**
!     **  (STRESS OF THE COMPENSATION CHARGES)                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: GID,NR,LMAXC,NG
      REAL(8)   ,INTENT(IN)  :: R(NR)
      REAL(8)   ,INTENT(IN)  :: GHAT(NR,0:LMAXC)
      REAL(8)   ,INTENT(IN)  :: GABS(NG)
      REAL(8)   ,INTENT(OUT) :: FG(NG,0:LMAXC)
      REAL(8)   ,INTENT(OUT) :: FGD(NG,0:LMAXC)
      REAL(8)   ,PARAMETER   :: DG=0.004D0
      REAL(8)   ,ALLOCATABLE :: TAB(:,:),TABD(:,:),WORK(:)
      REAL(8)                :: PI,G,X,FAC(4),T,TD,X1,DJ
      INTEGER(4)             :: NT,IT,L,IR,IG,I0,K,J
      REAL(8)   ,EXTERNAL    :: FOCK_SBESSEL
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      NT=INT(MAXVAL(GABS)/DG)+4
      ALLOCATE(TAB(NT,0:LMAXC))
      ALLOCATE(TABD(NT,0:LMAXC))
      ALLOCATE(WORK(NR))
      DO IT=1,NT
        G=DG*REAL(IT-2,KIND=8)       ! TABLE STARTS AT G=-DG (SYMMETRIC/ODD EXT.)
        DO L=0,LMAXC
          DO IR=1,NR
            IF(GHAT(IR,L).EQ.0.D0) THEN
              WORK(IR)=0.D0
            ELSE
              WORK(IR)=GHAT(IR,L)*FOCK_SBESSEL(L,ABS(G)*R(IR))*R(IR)**2
            END IF
          ENDDO
          CALL RADIAL$INTEGRAL(GID,NR,WORK,TAB(IT,L))
          TAB(IT,L)=4.D0*PI*TAB(IT,L)
          IF(G.LT.0.D0.AND.MOD(L,2).EQ.1) TAB(IT,L)=-TAB(IT,L) ! J_L(-X)=(-1)^L J_L(X)
!         == DERIVATIVE: J_L'(-X)=(-1)^(L+1) J_L'(X) ==============================
          DO IR=1,NR
            IF(GHAT(IR,L).EQ.0.D0) THEN
              WORK(IR)=0.D0
            ELSE
              X1=ABS(G)*R(IR)
              DJ=-REAL(L+1,KIND=8)*FOCK_SBESSEL(L+1,X1)
              IF(L.GT.0) DJ=DJ+REAL(L,KIND=8)*FOCK_SBESSEL(L-1,X1)
              DJ=DJ/REAL(2*L+1,KIND=8)
              WORK(IR)=GHAT(IR,L)*DJ*R(IR)**3
            END IF
          ENDDO
          CALL RADIAL$INTEGRAL(GID,NR,WORK,TABD(IT,L))
          TABD(IT,L)=4.D0*PI*TABD(IT,L)
          IF(G.LT.0.D0.AND.MOD(L,2).EQ.0) TABD(IT,L)=-TABD(IT,L)
        ENDDO
      ENDDO
      DO IG=1,NG
        X=GABS(IG)/DG+2.D0            ! TABLE COORDINATE (1-BASED)
        I0=INT(X)-1                   ! POINTS I0..I0+3 AROUND X
        I0=MAX(1,MIN(I0,NT-3))
        DO K=1,4
          FAC(K)=1.D0
          DO J=1,4
            IF(J.EQ.K) CYCLE
            FAC(K)=FAC(K)*(X-REAL(I0+J-1,KIND=8))/REAL(K-J,KIND=8)
          ENDDO
        ENDDO
        DO L=0,LMAXC
          T=0.D0
          TD=0.D0
          DO K=1,4
            T=T+FAC(K)*TAB(I0+K-1,L)
            TD=TD+FAC(K)*TABD(I0+K-1,L)
          ENDDO
          FG(IG,L)=T
          FGD(IG,L)=TD
        ENDDO
      ENDDO
      DEALLOCATE(TAB,TABD,WORK)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      REAL(8) FUNCTION FOCK_SBESSEL(L,X)
!     **************************************************************************
!     **  SPHERICAL BESSEL FUNCTION J_L(X), X>=0: POWER SERIES FOR X<L+2,     **
!     **  UPWARD RECURSION (STABLE FOR X>L) OTHERWISE                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: X
      REAL(8)               :: TERM,SUM,DFAC,J0,J1,J2
      INTEGER(4)            :: K
!     **************************************************************************
      IF(X.LT.REAL(L+2,KIND=8)) THEN
        DFAC=1.D0
        DO K=1,L
          DFAC=DFAC*REAL(2*K+1,KIND=8)
        ENDDO
        TERM=X**L/DFAC
        SUM=TERM
        DO K=1,60
          TERM=-TERM*X**2/(2.D0*REAL(K,KIND=8)*REAL(2*L+2*K+1,KIND=8))
          SUM=SUM+TERM
          IF(ABS(TERM).LT.1.D-17*ABS(SUM)) EXIT
        ENDDO
        FOCK_SBESSEL=SUM
      ELSE
        J0=SIN(X)/X
        IF(L.EQ.0) THEN
          FOCK_SBESSEL=J0
          RETURN
        END IF
        J1=SIN(X)/X**2-COS(X)/X
        DO K=1,L-1
          J2=REAL(2*K+1,KIND=8)/X*J1-J0
          J0=J1
          J1=J2
        ENDDO
        FOCK_SBESSEL=J1
      END IF
      RETURN
      END


!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_SELFENERGY(GID,NR,R,L,GHAT,QMAX,SMU,S)
!     **************************************************************************
!     **  S = 1/(2PI^2) INT_0^QMAX FG(Q)^2 W(Q) DQ,                           **
!     **  FG = 4PI INT GHAT J_L(QR) R^2 DR, W=1 (SMU=0) OR 1-EXP(-Q^2/4SMU^2) **
!     **  = (SCREENED) COULOMB SELF-ENERGY OF GHAT_L*Y_L FOR |Q|<QMAX          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: GID,NR,L
      REAL(8)   ,INTENT(IN)  :: R(NR),GHAT(NR),QMAX,SMU
      REAL(8)   ,INTENT(OUT) :: S
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      REAL(8)                :: PI,DQ,Q,FQ,W
      INTEGER(4)             :: NQ,IQ,IR
      REAL(8)   ,EXTERNAL    :: FOCK_SBESSEL
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      NQ=2*INT(QMAX/0.01D0)+1
      DQ=QMAX/REAL(NQ-1,KIND=8)
      ALLOCATE(WORK(NR))
      S=0.D0
      DO IQ=1,NQ
        Q=DQ*REAL(IQ-1,KIND=8)
        DO IR=1,NR
          IF(GHAT(IR).EQ.0.D0) THEN
            WORK(IR)=0.D0
          ELSE
            WORK(IR)=GHAT(IR)*FOCK_SBESSEL(L,Q*R(IR))*R(IR)**2
          END IF
        ENDDO
        CALL RADIAL$INTEGRAL(GID,NR,WORK,FQ)
        FQ=4.D0*PI*FQ
        W=2.D0
        IF(MOD(IQ,2).EQ.0) W=4.D0
        IF(IQ.EQ.1.OR.IQ.EQ.NQ) W=1.D0
        IF(SMU.GT.0.D0) W=W*(1.D0-EXP(-0.25D0*Q**2/SMU**2))
        S=S+W*FQ**2
      ENDDO
      S=S*DQ/3.D0/(2.D0*PI**2)
      DEALLOCATE(WORK)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_V0(RBAS,ALPHA,V0)
!     **************************************************************************
!     **  KERNEL AT G=0 FOR THE GAMMA-ONLY MESH = AVERAGE OF V(Q) OVER THE BZ **
!     **                                                                      **
!     **  COULOMB (MU=0): EWALD FORM OF -2*OMEGA*E_MADELUNG(UNIT CHARGE+BG),  **
!     **    V0 = OMEGA/SQRT(PI ALPHA) - SUM_{G/=0} 4PI EXP(-ALPHA G^2)/G^2     **
!     **       + 4PI ALPHA - OMEGA SUM_{R/=0} ERFC(R/(2 SQRT(ALPHA)))/R        **
!     **  (EXACT FOR ANY ALPHA; VASP SET_FSG, MONKHORST-PACK BRANCH)          **
!     **  ERFC-SCREENED (MU>0): THE SAME DEFINITION (KERNEL OF A POINT CHARGE **
!     **  MINUS ITS PERIODIC IMAGES), CONVERGENT WITHOUT EWALD SPLITTING:     **
!     **    V0 = PI/MU^2 - OMEGA SUM_{R/=0} ERFC(MU R)/R   (ALPHA UNUSED)     **
!     **  IT IS THE ALPHA->0 LIMIT OF THE AUXILIARY FORM                      **
!     **    OMEGA/SQRT(PI)[1/SQRT(A)-1/SQRT(A+B)] - SUM' 4PI E^{-AG^2}       **
!     **    (1-E^{-BG^2})/G^2, B=1/(4MU^2), WHICH DEPENDS ON ALPHA (SI8:      **
!     **    167.34/168.55/169.16 FOR A=.25,.125,.0625 VS 169.77) AND GAVE A    **
!     **    4.5 MH CDUAL DEPENDENCE; AND IT TENDS TO THE COULOMB VALUE FOR    **
!     **    MU->0 (298.75-2 OMEGA MU/SQRT(PI)).                               **
!     **************************************************************************
      USE FOCK_MODULE, ONLY : MU
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: RBAS(3,3)
      REAL(8),INTENT(IN)  :: ALPHA
      REAL(8),INTENT(OUT) :: V0
      REAL(8)             :: GBAS(3,3),OMEGA,PI,G(3),G2,GMAX2,R(3),RR,RMAX
      INTEGER(4)          :: N1,N2,N3,NX1,NX2,NX3,I
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL GBASS(RBAS,GBAS,OMEGA)
      IF(MU.GT.0.D0) THEN
!       == ERFC(X) < 1.E-18 FOR X > 6.2 ========================================
        RMAX=6.2D0/MU
        NX1=INT(RMAX/SQRT(SUM(RBAS(:,1)**2)))+1
        NX2=INT(RMAX/SQRT(SUM(RBAS(:,2)**2)))+1
        NX3=INT(RMAX/SQRT(SUM(RBAS(:,3)**2)))+1
        V0=PI/MU**2
        DO N1=-NX1,NX1
          DO N2=-NX2,NX2
            DO N3=-NX3,NX3
              IF(N1.EQ.0.AND.N2.EQ.0.AND.N3.EQ.0) CYCLE
              R(:)=RBAS(:,1)*N1+RBAS(:,2)*N2+RBAS(:,3)*N3
              RR=SQRT(SUM(R**2))
              IF(RR.GT.RMAX) CYCLE
              V0=V0-OMEGA*ERFC(MU*RR)/RR
            ENDDO
          ENDDO
        ENDDO
        RETURN
      END IF
!     == RECIPROCAL SPACE: EXP(-ALPHA G^2) < 1.E-17 ============================
      GMAX2=40.D0/ALPHA
      NX1=INT(SQRT(GMAX2)*SQRT(SUM(RBAS(:,1)**2))/(2.D0*PI))+1
      NX2=INT(SQRT(GMAX2)*SQRT(SUM(RBAS(:,2)**2))/(2.D0*PI))+1
      NX3=INT(SQRT(GMAX2)*SQRT(SUM(RBAS(:,3)**2))/(2.D0*PI))+1
      V0=0.D0
      DO N1=-NX1,NX1
        DO N2=-NX2,NX2
          DO N3=-NX3,NX3
            IF(N1.EQ.0.AND.N2.EQ.0.AND.N3.EQ.0) CYCLE
            DO I=1,3
              G(I)=GBAS(I,1)*N1+GBAS(I,2)*N2+GBAS(I,3)*N3
            ENDDO
            G2=SUM(G**2)
            IF(G2.GT.GMAX2) CYCLE
            V0=V0-4.D0*PI*EXP(-ALPHA*G2)/G2
          ENDDO
        ENDDO
      ENDDO
      V0=V0+OMEGA/SQRT(PI*ALPHA)+4.D0*PI*ALPHA
!     == REAL SPACE: ERFC(R/(2 SQRT(ALPHA))) < 1.E-17 ==========================
      RMAX=2.D0*SQRT(ALPHA)*6.2D0
      NX1=INT(RMAX/SQRT(SUM(RBAS(:,1)**2)))+1
      NX2=INT(RMAX/SQRT(SUM(RBAS(:,2)**2)))+1
      NX3=INT(RMAX/SQRT(SUM(RBAS(:,3)**2)))+1
      DO N1=-NX1,NX1
        DO N2=-NX2,NX2
          DO N3=-NX3,NX3
            IF(N1.EQ.0.AND.N2.EQ.0.AND.N3.EQ.0) CYCLE
            R(:)=RBAS(:,1)*N1+RBAS(:,2)*N2+RBAS(:,3)*N3
            RR=SQRT(SUM(R**2))
            IF(RR.GT.RMAX) CYCLE
            V0=V0-OMEGA*ERFC(RR/(2.D0*SQRT(ALPHA)))/RR
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_BESSELTRANSFORM(GID,NR,R,L,N,RHO,NQ,DQ,RQ)
!     **************************************************************************
!     **  RQ(IQ,I) = INT RHO(R,I) J_L(Q R) R^2 DR, Q=(IQ-1)*DQ                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: GID,NR,L,N,NQ
      REAL(8)   ,INTENT(IN)  :: R(NR),RHO(NR,N),DQ
      REAL(8)   ,INTENT(OUT) :: RQ(NQ,N)
      REAL(8)   ,ALLOCATABLE :: BJ(:),WORK(:)
      INTEGER(4)             :: IQ,IR,I
      REAL(8)   ,EXTERNAL    :: FOCK_SBESSEL
!     **************************************************************************
      ALLOCATE(BJ(NR),WORK(NR))
      DO IQ=1,NQ
        DO IR=1,NR
          BJ(IR)=FOCK_SBESSEL(L,DQ*REAL(IQ-1,KIND=8)*R(IR))*R(IR)**2
        ENDDO
        DO I=1,N
          WORK(:)=RHO(:,I)*BJ(:)
          CALL RADIAL$INTEGRAL(GID,NR,WORK,RQ(IQ,I))
        ENDDO
      ENDDO
      DEALLOCATE(BJ,WORK)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ERFQGRID(SMU,NQ,DQ)
!     **************************************************************************
!     **  Q GRID FOR THE ERF PART: EXP(-Q^2/4MU^2) < 1.E-16 BEYOND QMAX,      **
!     **  ODD NUMBER OF POINTS FOR SIMPSON                                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: SMU
      INTEGER(4),INTENT(OUT) :: NQ
      REAL(8)   ,INTENT(OUT) :: DQ
      REAL(8)                :: QMAX
!     **************************************************************************
      QMAX=2.D0*SMU*SQRT(37.D0)
      NQ=2*MAX(200,INT(QMAX/0.004D0))+1
      DQ=QMAX/REAL(NQ-1,KIND=8)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      REAL(8) FUNCTION FOCK_SIMPSONW(IQ,NQ)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IQ,NQ
      FOCK_SIMPSONW=2.D0
      IF(MOD(IQ,2).EQ.0) FOCK_SIMPSONW=4.D0
      IF(IQ.EQ.1.OR.IQ.EQ.NQ) FOCK_SIMPSONW=1.D0
      FOCK_SIMPSONW=FOCK_SIMPSONW/3.D0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ERFPART(GID,NR,R,L,NPAIR,RHOAE,RHOPS,SMU,F)
!     **************************************************************************
!     **  REMOVES THE LONG-RANGE ERF(MU R)/R PART FROM THE RADIAL SLATER      **
!     **  INTEGRALS F^L(AB,CD) (COULOMB -> ERFC KERNEL):                      **
!     **    (RHO_1 Y_L | ERF | RHO_2 Y_L) = 8 INT DQ EXP(-Q^2/4MU^2)           **
!     **                                    RHO_1(Q) RHO_2(Q),                **
!     **    RHO(Q) = INT RHO(R) J_L(QR) R^2 DR                                **
!     **  (FROM 4PI/Q^2 EXP(-Q^2/4MU^2), THE SAME FORM AS THE PLANE-WAVE      **
!     **  KERNEL, SO GRID AND ONE-CENTER TERMS USE THE SAME INTERACTION.      **
!     **  VASP BUILDS THE SAME KERNEL IN REAL SPACE, RS_GREENSFUNC.)          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: GID,NR,L,NPAIR
      REAL(8)   ,INTENT(IN)    :: R(NR),RHOAE(NR,NPAIR),RHOPS(NR,NPAIR),SMU
      REAL(8)   ,INTENT(INOUT) :: F(NPAIR,NPAIR)
      REAL(8)   ,ALLOCATABLE   :: QAE(:,:),QPS(:,:),W(:)
      REAL(8)                  :: DQ
      INTEGER(4)               :: NQ,IQ,I,J
      REAL(8)   ,EXTERNAL      :: FOCK_SIMPSONW
!     **************************************************************************
      CALL FOCK_ERFQGRID(SMU,NQ,DQ)
      ALLOCATE(QAE(NQ,NPAIR),QPS(NQ,NPAIR),W(NQ))
      CALL FOCK_BESSELTRANSFORM(GID,NR,R,L,NPAIR,RHOAE,NQ,DQ,QAE)
      CALL FOCK_BESSELTRANSFORM(GID,NR,R,L,NPAIR,RHOPS,NQ,DQ,QPS)
      DO IQ=1,NQ
        W(IQ)=8.D0*DQ*FOCK_SIMPSONW(IQ,NQ) &
     &       *EXP(-0.25D0*(DQ*REAL(IQ-1,KIND=8))**2/SMU**2)
      ENDDO
      DO I=1,NPAIR
        DO J=1,NPAIR
          F(I,J)=F(I,J)-SUM(W(:)*(QAE(:,I)*QAE(:,J)-QPS(:,I)*QPS(:,J)))
        ENDDO
      ENDDO
      DEALLOCATE(QAE,QPS,W)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOCK_ERFPARTCV(GID,NR,R,L,N,RHO,SMU,F)
!     **************************************************************************
!     **  ERF PART OF THE CORE-VALENCE INTEGRALS (SEE FOCK_ERFPART), AE ONLY  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: GID,NR,L,N
      REAL(8)   ,INTENT(IN)    :: R(NR),RHO(NR,N),SMU
      REAL(8)   ,INTENT(INOUT) :: F(N,N)
      REAL(8)   ,ALLOCATABLE   :: Q(:,:),W(:)
      REAL(8)                  :: DQ
      INTEGER(4)               :: NQ,IQ,I,J
      REAL(8)   ,EXTERNAL      :: FOCK_SIMPSONW
!     **************************************************************************
      CALL FOCK_ERFQGRID(SMU,NQ,DQ)
      ALLOCATE(Q(NQ,N),W(NQ))
      CALL FOCK_BESSELTRANSFORM(GID,NR,R,L,N,RHO,NQ,DQ,Q)
      DO IQ=1,NQ
        W(IQ)=8.D0*DQ*FOCK_SIMPSONW(IQ,NQ) &
     &       *EXP(-0.25D0*(DQ*REAL(IQ-1,KIND=8))**2/SMU**2)
      ENDDO
      DO I=1,N
        DO J=1,N
          F(I,J)=F(I,J)-SUM(W(:)*Q(:,I)*Q(:,J))
        ENDDO
      ENDDO
      DEALLOCATE(Q,W)
      RETURN
      END
