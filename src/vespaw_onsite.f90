!........1.........2.........3.........4.........5.........6.........7.........8
!**************************************************************************
!**                                                                      **
!**  ONE-CENTER ADDITIONS TO THE PAW HAMILTONIAN THAT DEPEND ONLY ON THE **
!**  ONE-CENTER DENSITY MATRICES:                                        **
!**                                                                      **
!**  DFT+U (DUDAREV, PRB 57, 1505 (1998); VASP LDAUTYPE=2):              **
!**    E_U = (U-J)/2 SUM TR[ N - N N ],                                  **
!**    N(M S,M' S') = SUM_{I,J IN SHELL L} Q_IJ D(I M S, J M' S')        **
!**    WITH Q_IJ=<AEPHI_I|AEPHI_J> INSIDE RCUT (THE PAW-SPHERE OCCUPATION**
!**    MATRIX, AS IN VASP). N IS THE FULL SPINOR OCCUPATION MATRIX, SO   **
!**    THE SAME CODE COVERS NON-SPIN, COLLINEAR AND NONCOLLINEAR RUNS.   **
!**                                                                      **
!**  SPIN-ORBIT COUPLING (NONCOLLINEAR RUNS, NSPIN=3):                   **
!**    H_SO = XI(R) L.S,  XI = 1/(2 M^2 C^2) 1/R DV/DR,                  **
!**    M = 1 - V/(2C^2) (SCALAR-RELATIVISTIC MASS, V THE SPHERICAL       **
!**    ALL-ELECTRON POTENTIAL OF THE ISOLATED ATOM), EVALUATED BETWEEN   **
!**    ALL-ELECTRON PARTIAL WAVES INSIDE RCUT.                           **
!**                                                                      **
!**  CONVENTIONS (SEE WAVES_DENMAT, WAVES_HPROJ): D(:,:,1..NDIMD) =      **
!**  (TOTAL[,SPIN]) OR (TOTAL,X,Y,Z); THE SPINOR HAMILTONIAN IS          **
!**  DH_T*1+SUM_A DH_A SIGMA_A AND E = SUM_A SUM_IJ DH_A(I,J) D_A(J,I).  **
!**                                                                      **
!**************************************************************************
MODULE ONSITE_MODULE
IMPLICIT NONE
TYPE ONSITESP_TYPE
  LOGICAL(4)          :: TU=.FALSE.
  INTEGER(4)          :: LU=-1           ! ANGULAR MOMENTUM OF THE U SHELL
  REAL(8)             :: UMJ=0.D0        ! U-J [HARTREE]
  REAL(8)             :: RCUT=1.D0       ! INTEGRATION RADIUS IN UNITS OF RCOV
  LOGICAL(4)          :: TQ=.FALSE.
  REAL(8),ALLOCATABLE :: Q(:,:)          ! (LNX,LNX) <AEPHI|AEPHI> INSIDE RCUT
  REAL(8),ALLOCATABLE :: K(:,:)          ! (LNX,LNX) SPIN-ORBIT RADIAL INTEGRALS
END TYPE ONSITESP_TYPE
TYPE(ONSITESP_TYPE),ALLOCATABLE :: SP(:)
INTEGER(4)          :: NSP=0
LOGICAL(4)          :: TSO=.FALSE.
INTEGER(4),PARAMETER:: LXSO=3
COMPLEX(8)          :: LMAT(16,16,3)    ! <Y_LM|L_A|Y_LM'>, L<=3
LOGICAL(4)          :: TLMAT=.FALSE.
END MODULE ONSITE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONSITE$SETDFTU(NSP_,ISP,L,UMJ,RCUT)
!     **  DEFINES A HUBBARD SHELL FOR SPECIES ISP (UMJ=U-J IN HARTREE,       **
!     **  RCUT IN UNITS OF THE COVALENT RADIUS)                              **
      USE ONSITE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSP_,ISP,L
      REAL(8)   ,INTENT(IN) :: UMJ,RCUT
      IF(.NOT.ALLOCATED(SP)) THEN
        NSP=NSP_
        ALLOCATE(SP(NSP))
      END IF
      SP(ISP)%TU=.TRUE.
      SP(ISP)%LU=L
      SP(ISP)%UMJ=UMJ
      SP(ISP)%RCUT=RCUT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONSITE$SETSO(NSP_,TSO_)
      USE ONSITE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSP_
      LOGICAL(4),INTENT(IN) :: TSO_
      TSO=TSO_          ! SPECIES DATA ARE ALLOCATED LATER (ONSITE$ETOT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONSITE$GETL4(ID,VAL)
      USE ONSITE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
      INTEGER(4)              :: ISP
      IF(ID.EQ.'ON') THEN
        VAL=TSO
        DO ISP=1,NSP
          VAL=VAL.OR.SP(ISP)%TU
        ENDDO
      ELSE IF(ID.EQ.'SO') THEN
        VAL=TSO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('ONSITE$GETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONSITE$ETOT(LMNXX,NDIMD,NAT,DENMAT,DH,EU,ESO)
!     **************************************************************************
!     **  ADDS DE/DD OF THE DFT+U AND SPIN-ORBIT TERMS TO DH; RETURNS ENERGIES **
!     **************************************************************************
      USE ONSITE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LMNXX,NDIMD,NAT
      COMPLEX(8),INTENT(IN)    :: DENMAT(LMNXX,LMNXX,NDIMD,NAT)
      COMPLEX(8),INTENT(INOUT) :: DH(LMNXX,LMNXX,NDIMD,NAT)
      REAL(8)   ,INTENT(OUT)   :: EU,ESO
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)
      INTEGER(4),ALLOCATABLE   :: ISPECIES(:),LOX(:),LNOFLMN(:),MOFLMN(:)
      COMPLEX(8),ALLOCATABLE   :: D1(:,:,:,:),G1(:,:,:,:),N(:,:),V(:,:)
      INTEGER(4)               :: IAT,ISP,LNX,LMNX,LN,LMN,I,J,IS,JS,M,L,NM
      INTEGER(4)               :: LMN1,LMN2,IA,IB,LN1,LN2
      COMPLEX(8)               :: CSVAR
!     **************************************************************************
      EU=0.D0
      ESO=0.D0
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      IF(.NOT.ALLOCATED(SP)) THEN
        NSP=MAXVAL(ISPECIES)
        ALLOCATE(SP(NSP))
      END IF
      IF(TSO.AND.NDIMD.NE.4) THEN
        CALL ERROR$MSG('SPIN-ORBIT COUPLING REQUIRES NONCOLLINEAR SPIN (NSPIN=3)')
        CALL ERROR$STOP('ONSITE$ETOT')
      END IF
      IF(TSO.AND..NOT.TLMAT) CALL ONSITE_LMATRICES()
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IF(.NOT.(SP(ISP)%TU.OR.TSO)) CYCLE
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        IF(.NOT.SP(ISP)%TQ) CALL ONSITE_RADIAL(ISP,LNX,LOX)
        CALL SETUP$UNSELECT()
        LMNX=SUM(2*LOX+1)
        ALLOCATE(LNOFLMN(LMNX),MOFLMN(LMNX))
        LMN=0
        DO LN=1,LNX
          DO M=1,2*LOX(LN)+1
            LMN=LMN+1
            LNOFLMN(LMN)=LN
            MOFLMN(LMN)=M
          ENDDO
        ENDDO
!       == SPINOR DENSITY MATRIX D1(I,J,S,S')=SUM F <P_I|PSI_S><PSI_S'|P_J> ====
        ALLOCATE(D1(LMNX,LMNX,2,2),G1(LMNX,LMNX,2,2))
        D1=(0.D0,0.D0)
        IF(NDIMD.EQ.1) THEN
          D1(:,:,1,1)=0.5D0*DENMAT(1:LMNX,1:LMNX,1,IAT)
          D1(:,:,2,2)=0.5D0*DENMAT(1:LMNX,1:LMNX,1,IAT)
        ELSE IF(NDIMD.EQ.2) THEN
          D1(:,:,1,1)=0.5D0*(DENMAT(1:LMNX,1:LMNX,1,IAT)+DENMAT(1:LMNX,1:LMNX,2,IAT))
          D1(:,:,2,2)=0.5D0*(DENMAT(1:LMNX,1:LMNX,1,IAT)-DENMAT(1:LMNX,1:LMNX,2,IAT))
        ELSE
          D1(:,:,1,1)=0.5D0*(DENMAT(1:LMNX,1:LMNX,1,IAT)+DENMAT(1:LMNX,1:LMNX,4,IAT))
          D1(:,:,2,2)=0.5D0*(DENMAT(1:LMNX,1:LMNX,1,IAT)-DENMAT(1:LMNX,1:LMNX,4,IAT))
          D1(:,:,1,2)=0.5D0*(DENMAT(1:LMNX,1:LMNX,2,IAT)-CI*DENMAT(1:LMNX,1:LMNX,3,IAT))
          D1(:,:,2,1)=0.5D0*(DENMAT(1:LMNX,1:LMNX,2,IAT)+CI*DENMAT(1:LMNX,1:LMNX,3,IAT))
        END IF
!       == G1=DE/DD1(I,J,S,S') ==================================================
        G1=(0.D0,0.D0)
!
!       == DFT+U ================================================================
        IF(SP(ISP)%TU) THEN
          L=SP(ISP)%LU
          NM=2*L+1
          ALLOCATE(N(2*NM,2*NM),V(2*NM,2*NM))
          N=(0.D0,0.D0)
          DO LMN1=1,LMNX
            LN1=LNOFLMN(LMN1)
            IF(LOX(LN1).NE.L) CYCLE
            DO LMN2=1,LMNX
              LN2=LNOFLMN(LMN2)
              IF(LOX(LN2).NE.L) CYCLE
              DO IS=1,2
                DO JS=1,2
                  IA=MOFLMN(LMN1)+NM*(IS-1)
                  IB=MOFLMN(LMN2)+NM*(JS-1)
                  N(IA,IB)=N(IA,IB)+SP(ISP)%Q(LN1,LN2)*D1(LMN1,LMN2,IS,JS)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CSVAR=(0.D0,0.D0)
          DO IA=1,2*NM
            CSVAR=CSVAR+N(IA,IA)
            DO IB=1,2*NM
              CSVAR=CSVAR-N(IA,IB)*N(IB,IA)
            ENDDO
          ENDDO
          EU=EU+0.5D0*SP(ISP)%UMJ*REAL(CSVAR,KIND=8)
!         == DE/DN(A,B)=(U-J)/2*(DELTA_AB-2*N(B,A)) ============================
          V=-SP(ISP)%UMJ*TRANSPOSE(N)
          DO IA=1,2*NM
            V(IA,IA)=V(IA,IA)+0.5D0*SP(ISP)%UMJ
          ENDDO
          DO LMN1=1,LMNX
            LN1=LNOFLMN(LMN1)
            IF(LOX(LN1).NE.L) CYCLE
            DO LMN2=1,LMNX
              LN2=LNOFLMN(LMN2)
              IF(LOX(LN2).NE.L) CYCLE
              DO IS=1,2
                DO JS=1,2
                  IA=MOFLMN(LMN1)+NM*(IS-1)
                  IB=MOFLMN(LMN2)+NM*(JS-1)
                  G1(LMN1,LMN2,IS,JS)=G1(LMN1,LMN2,IS,JS)+SP(ISP)%Q(LN1,LN2)*V(IA,IB)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(N,V)
        END IF
!
!       == SPIN-ORBIT: E=SUM K_IJ <I|L|J>.<S> -> DE/DD1(I,J,S,S')=K_JI/2 ======
!       == <J|L_A|I> SIGMA_A(S',S) =============================================
        IF(TSO) THEN
          DO LMN1=1,LMNX
            LN1=LNOFLMN(LMN1)
            DO LMN2=1,LMNX
              LN2=LNOFLMN(LMN2)
              IF(LOX(LN1).NE.LOX(LN2).OR.LOX(LN1).EQ.0.OR.LOX(LN1).GT.LXSO) CYCLE
              L=LOX(LN1)
              I=L**2+MOFLMN(LMN1)
              J=L**2+MOFLMN(LMN2)
!             == H_SO(J,S';I,S)=K/2 <J|L|I>.SIGMA(S',S);  DE/DD1(I,J,S,S')=H(J,S';I,S)
              CSVAR=0.5D0*SP(ISP)%K(LN2,LN1)
              G1(LMN1,LMN2,1,1)=G1(LMN1,LMN2,1,1)+CSVAR*LMAT(J,I,3)
              G1(LMN1,LMN2,2,2)=G1(LMN1,LMN2,2,2)-CSVAR*LMAT(J,I,3)
!             == SIGMA_X(S',S)+SIGMA_Y(S',S): (S',S)=(2,1) -> 1+I ; (1,2) -> 1-I
              G1(LMN1,LMN2,1,2)=G1(LMN1,LMN2,1,2)+CSVAR*(LMAT(J,I,1)+CI*LMAT(J,I,2))
              G1(LMN1,LMN2,2,1)=G1(LMN1,LMN2,2,1)+CSVAR*(LMAT(J,I,1)-CI*LMAT(J,I,2))
            ENDDO
          ENDDO
          CSVAR=(0.D0,0.D0)
          DO LMN1=1,LMNX
            DO LMN2=1,LMNX
              IF(SP(ISP)%TU) CYCLE
              DO IS=1,2
                DO JS=1,2
                  CSVAR=CSVAR+G1(LMN1,LMN2,IS,JS)*D1(LMN1,LMN2,IS,JS)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF(SP(ISP)%TU) CALL ONSITE_SOENERGY(LMNX,LNX,LOX,LNOFLMN,MOFLMN,ISP,D1,CSVAR)
          ESO=ESO+REAL(CSVAR,KIND=8)
        END IF
!
!       == DE/DD_A(J,I) -> DH_A(I,J) ============================================
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            IF(NDIMD.EQ.1) THEN
              DH(LMN2,LMN1,1,IAT)=DH(LMN2,LMN1,1,IAT) &
     &           +0.5D0*(G1(LMN1,LMN2,1,1)+G1(LMN1,LMN2,2,2))
            ELSE IF(NDIMD.EQ.2) THEN
              DH(LMN2,LMN1,1,IAT)=DH(LMN2,LMN1,1,IAT) &
     &           +0.5D0*(G1(LMN1,LMN2,1,1)+G1(LMN1,LMN2,2,2))
              DH(LMN2,LMN1,2,IAT)=DH(LMN2,LMN1,2,IAT) &
     &           +0.5D0*(G1(LMN1,LMN2,1,1)-G1(LMN1,LMN2,2,2))
            ELSE
              DH(LMN2,LMN1,1,IAT)=DH(LMN2,LMN1,1,IAT) &
     &           +0.5D0*(G1(LMN1,LMN2,1,1)+G1(LMN1,LMN2,2,2))
              DH(LMN2,LMN1,4,IAT)=DH(LMN2,LMN1,4,IAT) &
     &           +0.5D0*(G1(LMN1,LMN2,1,1)-G1(LMN1,LMN2,2,2))
              DH(LMN2,LMN1,2,IAT)=DH(LMN2,LMN1,2,IAT) &
     &           +0.5D0*(G1(LMN1,LMN2,1,2)+G1(LMN1,LMN2,2,1))
              DH(LMN2,LMN1,3,IAT)=DH(LMN2,LMN1,3,IAT) &
     &           +0.5D0*CI*(G1(LMN1,LMN2,2,1)-G1(LMN1,LMN2,1,2))
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(D1,G1,LOX,LNOFLMN,MOFLMN)
      ENDDO
      DEALLOCATE(ISPECIES)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONSITE_SOENERGY(LMNX,LNX,LOX,LNOFLMN,MOFLMN,ISP,D1,E)
!     **  SPIN-ORBIT ENERGY SUM K/2 <J|L|I>.SIGMA(S',S) D1(I,J,S,S') ALONE    **
!     **  (G1 ALSO HOLDS THE NONLINEAR DFT+U PART WHEN BOTH ARE ACTIVE)     **
      USE ONSITE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LMNX,LNX,ISP
      INTEGER(4),INTENT(IN)  :: LOX(LNX),LNOFLMN(LMNX),MOFLMN(LMNX)
      COMPLEX(8),INTENT(IN)  :: D1(LMNX,LMNX,2,2)
      COMPLEX(8),INTENT(OUT) :: E
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      INTEGER(4)             :: LMN1,LMN2,LN1,LN2,L,I,J
      COMPLEX(8)             :: C
      E=(0.D0,0.D0)
      DO LMN1=1,LMNX
        LN1=LNOFLMN(LMN1)
        DO LMN2=1,LMNX
          LN2=LNOFLMN(LMN2)
          IF(LOX(LN1).NE.LOX(LN2).OR.LOX(LN1).EQ.0.OR.LOX(LN1).GT.LXSO) CYCLE
          L=LOX(LN1)
          I=L**2+MOFLMN(LMN1)
          J=L**2+MOFLMN(LMN2)
          C=0.5D0*SP(ISP)%K(LN2,LN1)
          E=E+C*(LMAT(J,I,3)*(D1(LMN1,LMN2,1,1)-D1(LMN1,LMN2,2,2)) &
     &          +(LMAT(J,I,1)+CI*LMAT(J,I,2))*D1(LMN1,LMN2,1,2) &
     &          +(LMAT(J,I,1)-CI*LMAT(J,I,2))*D1(LMN1,LMN2,2,1))
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONSITE_RADIAL(ISP,LNX,LOX)
!     **************************************************************************
!     **  Q_IJ=INT_0^RCUT R^2 PHI_I PHI_J AND THE SPIN-ORBIT INTEGRALS K_IJ   **
!     **  (SETUP ISP MUST BE SELECTED)                                        **
!     **************************************************************************
      USE ONSITE_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP,LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,PARAMETER  :: C=137.035999084D0
      INTEGER(4)            :: GID,NR,LN1,LN2,IR
      REAL(8)               :: AEZ,RCOV,RCUT,Y0,PI
      REAL(8)   ,ALLOCATABLE:: R(:),PHI(:,:),POT(:),DPOT(:),XI(:),WORK(:)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL SETUP$GETI4('GID',GID)
      CALL SETUP$GETI4('NR',NR)
      CALL SETUP$GETR8('AEZ',AEZ)
      ALLOCATE(R(NR),PHI(NR,LNX),POT(NR),DPOT(NR),XI(NR),WORK(NR))
      CALL RADIAL$R(GID,NR,R)
      CALL SETUP$GETR8A('AEPHI',NR*LNX,PHI)
      CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RCOV)
      RCUT=SP(ISP)%RCUT*RCOV
      ALLOCATE(SP(ISP)%Q(LNX,LNX),SP(ISP)%K(LNX,LNX))
      SP(ISP)%Q=0.D0
      SP(ISP)%K=0.D0
!     == SPHERICAL AE POTENTIAL OF THE ATOM (STORED AS Y0-COEFFICIENT) ==========
      CALL SETUP$GETR8A('AEPOT',NR,POT)
      POT=POT*Y0
      CALL RADIAL$DERIVE(GID,NR,POT,DPOT)
      XI=0.D0
      DO IR=1,NR
        IF(R(IR).LT.1.D-10) CYCLE
        XI(IR)=DPOT(IR)/R(IR)/(2.D0*C**2*(1.D0-POT(IR)/(2.D0*C**2))**2)
      ENDDO
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          WORK=R**2*PHI(:,LN1)*PHI(:,LN2)
          WHERE(R.GT.RCUT) WORK=0.D0
          CALL RADIAL$INTEGRAL(GID,NR,WORK,SP(ISP)%Q(LN1,LN2))
          WORK=WORK*XI
          CALL RADIAL$INTEGRAL(GID,NR,WORK,SP(ISP)%K(LN1,LN2))
        ENDDO
      ENDDO
      SP(ISP)%TQ=.TRUE.
      WRITE(*,FMT='("ONSITE: SPECIES",I3," RCUT=",F8.4," Q(DIAG)=",10F8.4)') &
     &      ISP,RCUT,(SP(ISP)%Q(LN1,LN1),LN1=1,LNX)
      DEALLOCATE(R,PHI,POT,DPOT,XI,WORK)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONSITE_LMATRICES()
!     **************************************************************************
!     **  <Y_I|L_A|Y_J> FOR REAL SPHERICAL HARMONICS, L=-I R X GRAD,          **
!     **  BY ANGULAR QUADRATURE (EXACT FOR L<=3), CHECKED AGAINST [LX,LY]=ILZ **
!     **************************************************************************
      USE ONSITE_MODULE
      USE NSPH_SCAN_MODULE, ONLY : QUADX,QUADY,QUADZ,QUADW,QUADYLM &
     &                            ,QUADDYLMX,QUADDYLMY,QUADDYLMZ,NTHSAVE,NPHSAVE
      IMPLICIT NONE
      COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)
      INTEGER(4),PARAMETER :: NTH=14,NPH=24
      INTEGER(4)           :: IQ,I,J,L,LM1,LM2
      REAL(8)              :: RH(3),G(3),C(3),DEV
      COMPLEX(8)           :: COMM(16,16)
!     **************************************************************************
      CALL NSPH_SCAN_BUILDGRID(LXSO,NTH,NPH)
      LMAT=(0.D0,0.D0)
      DO IQ=1,NTHSAVE*NPHSAVE
        RH=(/QUADX(IQ),QUADY(IQ),QUADZ(IQ)/)
        DO J=1,(LXSO+1)**2
          G=(/QUADDYLMX(IQ,J),QUADDYLMY(IQ,J),QUADDYLMZ(IQ,J)/)
          C=(/RH(2)*G(3)-RH(3)*G(2),RH(3)*G(1)-RH(1)*G(3),RH(1)*G(2)-RH(2)*G(1)/)
          DO I=1,(LXSO+1)**2
            LMAT(I,J,:)=LMAT(I,J,:)-CI*QUADW(IQ)*QUADYLM(IQ,I)*C
          ENDDO
        ENDDO
      ENDDO
!     == ONLY BLOCKS WITH EQUAL L; CHECK THE ANGULAR-MOMENTUM ALGEBRA ===========
      DO L=0,LXSO
        DO LM1=L**2+1,(L+1)**2
          DO LM2=1,(LXSO+1)**2
            IF(LM2.LE.L**2.OR.LM2.GT.(L+1)**2) LMAT(LM1,LM2,:)=(0.D0,0.D0)
          ENDDO
        ENDDO
      ENDDO
      COMM=MATMUL(LMAT(:,:,1),LMAT(:,:,2))-MATMUL(LMAT(:,:,2),LMAT(:,:,1))-CI*LMAT(:,:,3)
      DEV=MAXVAL(ABS(COMM))
      IF(DEV.GT.1.D-8) THEN
        CALL ERROR$MSG('ANGULAR-MOMENTUM MATRICES VIOLATE [LX,LY]=ILZ')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('ONSITE_LMATRICES')
      END IF
      TLMAT=.TRUE.
      RETURN
      END
