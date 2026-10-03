!*******************************************************************************
!** CONTAINER FOR K-POINT-RELATED INFORMATION.                                **
!** DRIVER FOR THE PAW_BRILLOUIN OBJECT.                                      **
!** IS ALSO RESPONSIBLE FOR THE CASES                                         **
!** WHERE THE BRILLOUIN OBJECT IS NOT USED, SUCH AS ISOLATED MOLECULES        **
!*******************************************************************************
MODULE KPOINTS_MODULE
LOGICAL(4)          :: TINV           ! TIME-INVERSION SYMMETRY
INTEGER(4)          :: NKDIV(3)       ! DIVISION OF K-POINTS
INTEGER(4)          :: ISHIFT(3)      ! DISPLACEMENT OF GRID FROM GAMMA
INTEGER(4)          :: NKPT=0
REAL(8),ALLOCATABLE :: XK(:,:)  !(3,NKPT) K-POINTS IN RELATIVE COORDINATES
REAL(8),ALLOCATABLE :: WKPT(:)  ! GEOMETRIC K-POINT WEIGHTS
END MODULE KPOINTS_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$SETL4(ID,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : TINV
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'TINV') THEN
        TINV=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$GETL4(ID,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : TINV
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'TINV') THEN
        VAL=TINV
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$GETI4(ID,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : NKPT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'NKPT') THEN
        VAL=NKPT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$GETI4A(ID,LEN,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : NKDIV,ISHIFT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(OUT):: VAL(LEN)
!     **************************************************************************
      IF(LEN.NE.3) THEN
        CALL ERROR$MSG('INCONSISTENT SIZE')
        CALL ERROR$STOP('KPOINTS$GETI4A')
      END IF
      IF(ID.EQ.'DIV') THEN
        VAL=NKDIV
      ELSE IF(ID.EQ.'SHIFT') THEN
        VAL=ISHIFT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$GETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$SETI4A(ID,LEN,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : NKDIV,ISHIFT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: VAL(LEN)
!     **************************************************************************
      IF(ID.EQ.'DIV') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$SETI4A')
        END IF
        NKDIV=VAL
!
      ELSE IF(ID.EQ.'SHIFT') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$SETI4A')
        END IF
        ISHIFT=VAL
! 
     ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$GETR8A(ID,LEN,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : NKPT,WKPT,XK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
!     **************************************************************************
      IF(ID.EQ.'XK') THEN
        IF(LEN.NE.3*NKPT) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT 3*NKPT')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NKPT',NKPT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$GETR8A')
        END IF
        VAL=RESHAPE(XK,(/3*NKPT/))
!
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(LEN.NE.NKPT) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT NKPT')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NKPT',NKPT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$GETR8A')
        END IF
        VAL=WKPT
! 
     ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$GETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS_NKDIV(RBAS,RMAX,NKDIV)
!     **************************************************************************
!     **  HELPER ROUTINE                                                      **
!     **  CALCULATES THE DIVISION OF THE RECIPROCAL UNIT CELL                 **
!     **  CONSISTENT WITH A RADIUS RMAX IN REAL SPACE.                        **
!     **  THE K-POINT GRID IS CONSISTENT WITH A SUPERCELL WITH                **
!     **  THE LONGEST DIAMETER EQUAL TO RMAX. (CHECK THIS)                    **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: RMAX
      INTEGER(4),INTENT(OUT):: NKDIV(3)
      REAL(8)               :: GBAS(3,3)
      REAL(8)               :: VOL,SVAR
      INTEGER(4)            :: I
!     **************************************************************************
      CALL GBASS(RBAS,GBAS,VOL)
      DO I=1,3
        SVAR=RMAX*SQRT(DOT_PRODUCT(GBAS(:,I),GBAS(:,I))) &
     &           /DOT_PRODUCT(RBAS(:,I),GBAS(:,I))
        NKDIV(I)=INT(ABS(SVAR+1.D0))
      ENDDO
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE KPOINTS_KPOINTS(TINV,NKDIV,ISHIFT,NKPT,XK,WKPT)
!!$!     ******************************************************************
!!$!     **  SELECTS K-POINTS ON AN EQUISPACED GRID WITH THE PROVIDED    **
!!$!     **  DIVISIONS NKDIV. TIME INVERSION SYMMETRY IS CONSIDERED      **
!!$!     **  IF TINV=.TRUE.                                              **
!!$!     **                                                              **
!!$!     **  THE NUMBER OF K-POINTS MUST BE ESTIMATED USING KPOINTS_NKPT **
!!$!     **                                                              **
!!$!     ******************************************************************
!!$      IMPLICIT NONE
!!$      LOGICAL(4),INTENT(IN) :: TINV        ! TIME INVERSION SYMMETRY
!!$      INTEGER(4),INTENT(IN) :: NKDIV(3)    ! #(DIVISIONS)
!!$      INTEGER(4),INTENT(IN) :: ISHIFT(3)   ! DISPLACEMENT AWAY FROM GAMMA
!!$      INTEGER(4),INTENT(IN) :: NKPT        ! #(KPOINTS)
!!$      REAL(8)   ,INTENT(OUT):: XK(3,NKPT)  !KPOINT IN RELATIVE COORDINATES
!!$      REAL(8)   ,INTENT(OUT):: WKPT(NKPT)  ! GEOMETRIC K-POINT WEIGHT
!!$      REAL(8)               :: DUMMYRBAS(3,3)
!!$      INTEGER(4)            :: NKPT1
!!$!     ******************************************************************
!!$!     USE SIC UNIT CELL TO OBTAIN RELATIVE COORDINATES
!!$                               CALL TRACE$PUSH('KPOINTS_KPOINTS')
!!$      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
!!$      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
!!$      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
!!$      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
!!$      CALL BRILLOUIN$GETI4('NK',NKPT1)
!!$      IF(NKPT1.NE.NKPT) THEN
!!$        CALL ERROR$STOP('KPOINTS_KPOINTS')
!!$      END IF 
!!$      CALL BRILLOUIN$GETR8A('XK',3*NKPT,XK)
!!$      CALL BRILLOUIN$GETR8A('WKPT',NKPT,WKPT)
!!$                                CALL TRACE$POP
!!$      RETURN
!!$      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE KPOINTS_NKPT(TINV,NKDIV,ISHIFT,NKPT)
!!$!     ******************************************************************
!!$!     **                                                              **
!!$!     ******************************************************************
!!$      IMPLICIT NONE
!!$      LOGICAL(4),INTENT(IN) :: TINV
!!$      INTEGER(4),INTENT(IN) :: NKDIV(3)
!!$      INTEGER(4),INTENT(IN) :: ISHIFT(3)   ! DISPLACEMENT AWAY FROM GAMMA
!!$      INTEGER(4),INTENT(OUT):: NKPT
!!$      REAL(8)               :: DUMMYRBAS(3,3)
!!$!     ******************************************************************
!!$      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
!!$      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
!!$      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
!!$      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
!!$      CALL BRILLOUIN$GETI4('NK',NKPT)
!!$      RETURN
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$MAKEGRID()
!     ******************************************************************
!     **  SELECTS K-POINTS ON AN EQUISPACED GRID WITH THE PROVIDED    **
!     **  DIVISIONS NKDIV. TIME INVERSION SYMMETRY IS CONSIDERED      **
!     **  IF TINV=.TRUE.                                              **
!     **                                                              **
!     **  THE NUMBER OF K-POINTS MUST BE ESTIMATED USING KPOINTS_NKPT **
!     ******************************************************************
      USE KPOINTS_MODULE, ONLY : TINV,NKDIV,ISHIFT,NKPT,XK,WKPT
      IMPLICIT NONE
      REAL(8)               :: DUMMYRBAS(3,3)
      INTEGER(4) :: IKPT
!     ******************************************************************
                               CALL TRACE$PUSH('KPOINTS$MAKEGRID')
!     USE SIC UNIT CELL TO OBTAIN RELATIVE COORDINATES
      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
      CALL BRILLOUIN$GETI4('NK',NKPT)
      ALLOCATE(XK(3,NKPT))
      ALLOCATE(WKPT(NKPT))
      CALL BRILLOUIN$GETR8A('XK',3*NKPT,XK)
      CALL BRILLOUIN$GETR8A('WKPT',NKPT,WKPT)
                                CALL TRACE$POP
      RETURN
      END






!
!     ...1.........2.........3.........4.........5.........6.........7.........8
!*******************************************************************************
!**  SPACE-GROUP SYMMETRY (2026-09-28)                                        **
!**                                                                           **
!**  FINDS THE SPACE-GROUP OPERATIONS OF THE CRYSTAL FROM THE LATTICE AND THE **
!**  ATOMIC POSITIONS, REDUCES THE K-POINT MESH TO IRREDUCIBLE POINTS, AND    **
!**  SYMMETRIZES THE QUANTITIES THAT ARE SUMMED OVER THE IRREDUCIBLE WEDGE:   **
!**  PSEUDO DENSITY, ONE-CENTER DENSITY MATRICES, FORCES AND STRESS.          **
!**  AN OPERATION MAPS FRACTIONAL COORDINATES X -> IROT*X+TAU; ATOMS ARE      **
!**  EQUIVALENT WHEN THEY BELONG TO THE SAME SPECIES (USE DIFFERENT SPECIES   **
!**  TO DISTINGUISH MAGNETIC SUBLATTICES).                                    **
!*******************************************************************************
MODULE SYMMETRY_MODULE
LOGICAL(4)             :: TON=.FALSE.
INTEGER(4)             :: NOP=0
INTEGER(4)             :: NATSYM=0
INTEGER(4),ALLOCATABLE :: IROT(:,:,:)  !(3,3,NOP) FRACTIONAL ROTATION
REAL(8)   ,ALLOCATABLE :: TAU(:,:)     !(3,NOP)   FRACTIONAL TRANSLATION
REAL(8)   ,ALLOCATABLE :: CROT(:,:,:)  !(3,3,NOP) CARTESIAN ROTATION
INTEGER(4),ALLOCATABLE :: MAPAT(:,:)   !(NAT,NOP) ATOM I -> ATOM MAPAT(I,IOP)
END MODULE SYMMETRY_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY$SETL4(ID,VAL)
      USE SYMMETRY_MODULE, ONLY : TON
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
      IF(ID.EQ.'ON') THEN
        TON=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SYMMETRY$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY$GETL4(ID,VAL)
      USE SYMMETRY_MODULE, ONLY : TON,NOP
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
      IF(ID.EQ.'ON') THEN
        VAL=TON.AND.NOP.GT.1
      ELSE IF(ID.EQ.'REQUESTED') THEN
        VAL=TON
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SYMMETRY$GETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY$FIND(RBAS,NAT,R,ISPECIES)
!     **************************************************************************
!     **  FINDS ALL OPERATIONS (IROT,TAU) THAT MAP THE LATTICE ONTO ITSELF    **
!     **  (IROT^T G IROT = G FOR THE METRIC G) AND EVERY ATOM ONTO AN ATOM OF **
!     **  THE SAME SPECIES.                                                   **
!     **************************************************************************
      USE SYMMETRY_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)     ! LATTICE VECTORS (COLUMNS)
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: R(3,NAT)      ! CARTESIAN POSITIONS
      INTEGER(4),INTENT(IN) :: ISPECIES(NAT)
      INTEGER(4),PARAMETER  :: NOPX=1536     ! 48 X UP TO 32 PURE TRANSLATIONS
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      REAL(8)               :: RBINV(3,3),G(3,3),GM(3,3),X(3,NAT),XP(3),D(3)
      REAL(8)               :: T(3),TMAX
      INTEGER(4)            :: M(3,3),IROTX(3,3,NOPX),MAPX(NAT,NOPX),MAP1(NAT)
      REAL(8)               :: TAUX(3,NOPX)
      INTEGER(4)            :: I1,I2,I3,I4,I5,I6,I7,I8,I9,IAT,JAT,KAT,IOP,DET
      LOGICAL(4)            :: TOK,TFOUND
!     **************************************************************************
      CALL LIB$INVERTR8(3,RBAS,RBINV)
      X=MATMUL(RBINV,R)
      X=X-FLOOR(X+TOL)
      G=MATMUL(TRANSPOSE(RBAS),RBAS)
      TMAX=MAXVAL(ABS(G))
      NOP=0
      DO I1=-1,1; DO I2=-1,1; DO I3=-1,1
      DO I4=-1,1; DO I5=-1,1; DO I6=-1,1
      DO I7=-1,1; DO I8=-1,1; DO I9=-1,1
        M=RESHAPE((/I1,I2,I3,I4,I5,I6,I7,I8,I9/),(/3,3/))
        DET=M(1,1)*(M(2,2)*M(3,3)-M(3,2)*M(2,3)) &
     &     -M(1,2)*(M(2,1)*M(3,3)-M(3,1)*M(2,3)) &
     &     +M(1,3)*(M(2,1)*M(3,2)-M(3,1)*M(2,2))
        IF(ABS(DET).NE.1) CYCLE
        GM=MATMUL(TRANSPOSE(REAL(M,KIND=8)),MATMUL(G,REAL(M,KIND=8)))
        IF(MAXVAL(ABS(GM-G)).GT.TOL*TMAX) CYCLE
!       == TRANSLATIONS: MAP ATOM 1 ONTO EVERY ATOM OF ITS SPECIES ============
        DO JAT=1,NAT
          IF(ISPECIES(JAT).NE.ISPECIES(1)) CYCLE
          T=X(:,JAT)-MATMUL(REAL(M,KIND=8),X(:,1))
          T=T-FLOOR(T+TOL)
          TOK=.TRUE.
          DO IAT=1,NAT
            XP=MATMUL(REAL(M,KIND=8),X(:,IAT))+T
            TFOUND=.FALSE.
            DO KAT=1,NAT
              IF(ISPECIES(KAT).NE.ISPECIES(IAT)) CYCLE
              D=XP-X(:,KAT)
              D=D-NINT(D)
              IF(MAXVAL(ABS(D)).LT.TOL) THEN
                MAP1(IAT)=KAT
                TFOUND=.TRUE.
                EXIT
              END IF
            ENDDO
            IF(.NOT.TFOUND) THEN
              TOK=.FALSE.
              EXIT
            END IF
          ENDDO
          IF(.NOT.TOK) CYCLE
          NOP=NOP+1
          IF(NOP.GT.NOPX) THEN
            CALL ERROR$MSG('TOO MANY SYMMETRY OPERATIONS')
            CALL ERROR$STOP('SYMMETRY$FIND')
          END IF
          IROTX(:,:,NOP)=M
          TAUX(:,NOP)=T
          MAPX(:,NOP)=MAP1
        ENDDO
      ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
      IF(ALLOCATED(IROT)) DEALLOCATE(IROT,TAU,CROT,MAPAT)
      ALLOCATE(IROT(3,3,NOP),TAU(3,NOP),CROT(3,3,NOP),MAPAT(NAT,NOP))
      IROT=IROTX(:,:,1:NOP)
      TAU=TAUX(:,1:NOP)
      MAPAT=MAPX(:,1:NOP)
      DO IOP=1,NOP
        CROT(:,:,IOP)=MATMUL(RBAS,MATMUL(REAL(IROT(:,:,IOP),KIND=8),RBINV))
      ENDDO
      NATSYM=NAT
      WRITE(*,FMT='("SYMMETRY: ",I4," SPACE-GROUP OPERATIONS FOUND")')NOP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY$MAKEKGRID()
!     **************************************************************************
!     **  REPLACES THE K-POINT GRID BY ITS IRREDUCIBLE WEDGE. K IN RELATIVE   **
!     **  COORDINATES TRANSFORMS WITH (IROT^-1)^T.                            **
!     **************************************************************************
      USE SYMMETRY_MODULE, ONLY : NOP,IROT
      USE KPOINTS_MODULE, ONLY : TINV,NKDIV,ISHIFT,NKPT,XK,WKPT
      IMPLICIT NONE
      INTEGER(4),ALLOCATABLE :: IIO(:,:,:)
      REAL(8)                :: DUMMYRBAS(3,3),A(3,3),AINV(3,3)
      REAL(8)   ,ALLOCATABLE :: XKA(:,:),WKA(:)
      INTEGER(4)             :: IOP,NK0
!     **************************************************************************
      ALLOCATE(IIO(3,3,NOP))
      DO IOP=1,NOP
        A=REAL(IROT(:,:,IOP),KIND=8)
        CALL LIB$INVERTR8(3,A,AINV)
        IIO(:,:,IOP)=NINT(TRANSPOSE(AINV))
      ENDDO
      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
      CALL BRILLOUIN$MSHSYM(DUMMYRBAS,NKDIV,ISHIFT,TINV,NOP,IIO)
      DEALLOCATE(IIO)
      IF(ALLOCATED(XK)) DEALLOCATE(XK)
      IF(ALLOCATED(WKPT)) DEALLOCATE(WKPT)
      CALL BRILLOUIN$GETI4('NK',NKPT)
      ALLOCATE(XKA(3,NKPT))
      ALLOCATE(WKA(NKPT))
      CALL BRILLOUIN$GETR8A('XK',3*NKPT,XKA)
      CALL BRILLOUIN$GETR8A('WKPT',NKPT,WKA)
!     == DROP POINTS WITHOUT WEIGHT (ZONE-BOUNDARY DUPLICATES OF THE MESH) ====
      NK0=COUNT(WKA.GT.1.D-10)
      ALLOCATE(XK(3,NK0))
      ALLOCATE(WKPT(NK0))
      NK0=0
      DO IOP=1,NKPT
        IF(WKA(IOP).LE.1.D-10) CYCLE
        NK0=NK0+1
        XK(:,NK0)=XKA(:,IOP)
        WKPT(NK0)=WKA(IOP)
      ENDDO
      NKPT=NK0
      DEALLOCATE(XKA,WKA)
      DO IOP=1,NKPT
        WRITE(*,FMT='("SYMMETRY: K",I5,3F10.5," W=",F10.6)')IOP,XK(:,IOP),WKPT(IOP)
      ENDDO
      WRITE(*,FMT='("SYMMETRY: ",I6," IRREDUCIBLE K-POINTS, SUM OF WEIGHTS ",F12.8)') &
     &      NKPT,SUM(WKPT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY$RHO(NRL,NDIMD,RHO)
!     **************************************************************************
!     **  SYMMETRIZES THE PSEUDO DENSITY ON THE 'DENSITY' GRID. WITH MPI THE  **
!     **  GRID IS SPLIT INTO SHEETS ALONG THE FIRST LATTICE VECTOR; THE FULL  **
!     **  GRID IS ASSEMBLED ON EVERY TASK, SYMMETRIZED, AND EACH TASK KEEPS   **
!     **  ITS OWN SHEETS.                                                     **
!     **************************************************************************
      USE SYMMETRY_MODULE, ONLY : NOP
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NRL
      INTEGER(4),INTENT(IN)    :: NDIMD
      REAL(8)   ,INTENT(INOUT) :: RHO(NRL,NDIMD)
      REAL(8)   ,ALLOCATABLE   :: FULL(:,:,:,:)
      REAL(8)   ,ALLOCATABLE   :: LOC(:,:,:,:)
      INTEGER(4)               :: NR1,NR2,NR3,NR1L,NR1START,NTASKS,THISTASK
!     **************************************************************************
      IF(NOP.LE.1) RETURN
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$GETI4('NR1',NR1)
      CALL PLANEWAVE$GETI4('NR2',NR2)
      CALL PLANEWAVE$GETI4('NR3',NR3)
      CALL PLANEWAVE$GETI4('NR1L',NR1L)
      CALL PLANEWAVE$GETI4('NR1START',NR1START)
      IF(NR1L*NR2*NR3.NE.NRL) THEN
        CALL ERROR$MSG('GRID SIZE MISMATCH')
        CALL ERROR$STOP('SYMMETRY$RHO')
      END IF
      IF(NTASKS.EQ.1) THEN
        CALL SYMMETRY_RHOFULL(NRL,NDIMD,RHO)
        RETURN
      END IF
      ALLOCATE(FULL(NR1,NR2,NR3,NDIMD))
      ALLOCATE(LOC(NR1L,NR2,NR3,NDIMD))
      LOC=RESHAPE(RHO,(/NR1L,NR2,NR3,NDIMD/))
      FULL=0.D0
      FULL(NR1START:NR1START+NR1L-1,:,:,:)=LOC
      CALL MPE$COMBINE('MONOMER','+',FULL)
      CALL SYMMETRY_RHOFULL(NR1*NR2*NR3,NDIMD,FULL)
      LOC=FULL(NR1START:NR1START+NR1L-1,:,:,:)
      RHO=RESHAPE(LOC,(/NRL,NDIMD/))
      DEALLOCATE(FULL,LOC)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY_RHOFULL(NRL,NDIMD,RHO)
!     **************************************************************************
!     **  SYMMETRIZES THE PSEUDO DENSITY ON THE 'DENSITY' GRID:              **
!     **    RHO_SYM(X) = 1/N SUM_G RHO(M_G X + TAU_G)   (X IN RELATIVE COORD.)  **
!     **  EXACT ON THE REAL-SPACE GRID WHEN EVERY OPERATION MAPS GRID POINTS  **
!     **  ONTO GRID POINTS. OTHERWISE (E.G. TAU=1/4 ON AN 18^3 GRID) IN       **
!     **  RECIPROCAL SPACE: THE COEFFICIENT AT H'=M^T H IS                    **
!     **  RHO(H)*EXP(2*PI*I*H.TAU), WITH RHO(X)=SUM_H RHO(H)EXP(2*PI*I*H.X).  **
!     **************************************************************************
      USE SYMMETRY_MODULE, ONLY : NOP,IROT,TAU
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NRL
      INTEGER(4),INTENT(IN)    :: NDIMD
      REAL(8)   ,INTENT(INOUT) :: RHO(NRL,NDIMD)
      REAL(8)   ,PARAMETER     :: PI=3.141592653589793D0
      COMPLEX(8),ALLOCATABLE   :: FG(:,:,:),FS(:,:,:)
      INTEGER(4),ALLOCATABLE   :: MINVS(:,:,:)
      REAL(8)   ,ALLOCATABLE   :: DRHO(:,:)
      REAL(8)                  :: A(3,3),AINV(3,3),PH
      COMPLEX(8)               :: C
      INTEGER(4)               :: IOP,IDIM,I,J
      INTEGER(4)               :: H(3),HP(3)
      INTEGER(4)               :: NR(3),N0(3),NP(3),I1,I2,I3,IR0,IRP
      LOGICAL(4)               :: TGRID
      LOGICAL(4),SAVE          :: TWARNED=.FALSE.
!     **************************************************************************
      IF(NOP.LE.1) RETURN
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$GETI4('NR1',NR(1))
      CALL PLANEWAVE$GETI4('NR2',NR(2))
      CALL PLANEWAVE$GETI4('NR3',NR(3))
      IF(NR(1)*NR(2)*NR(3).NE.NRL) THEN
        CALL ERROR$MSG('GRID SIZE MISMATCH')
        CALL ERROR$STOP('SYMMETRY_RHOFULL')
      END IF
!
!     ==========================================================================
!     == EXACT SYMMETRIZATION ON THE REAL-SPACE GRID IF POSSIBLE              ==
!     ==========================================================================
      TGRID=.TRUE.
      DO IOP=1,NOP
        DO I=1,3
          IF(ABS(TAU(I,IOP)*NR(I)-NINT(TAU(I,IOP)*NR(I))).GT.1.D-6) TGRID=.FALSE.
          DO J=1,3
            IF(MOD(IROT(I,J,IOP)*NR(I),NR(J)).NE.0) TGRID=.FALSE.
          ENDDO
        ENDDO
      ENDDO
      IF(TGRID) THEN
        ALLOCATE(DRHO(NRL,NDIMD))
        DRHO=0.D0
        DO IOP=1,NOP
          DO I3=0,NR(3)-1
            DO I2=0,NR(2)-1
              DO I1=0,NR(1)-1
                N0=(/I1,I2,I3/)
!               == X'=M X+TAU WITH X=N/NR: N'_I=SUM_J M_IJ N_J NR_I/NR_J+TAU_I NR_I
                DO I=1,3
                  NP(I)=NINT(TAU(I,IOP)*NR(I))
                  DO J=1,3
                    NP(I)=NP(I)+(IROT(I,J,IOP)*NR(I)/NR(J))*N0(J)
                  ENDDO
                  NP(I)=MODULO(NP(I),NR(I))
                ENDDO
                IR0=1+I1+NR(1)*(I2+NR(2)*I3)
                IRP=1+NP(1)+NR(1)*(NP(2)+NR(2)*NP(3))
                DRHO(IR0,:)=DRHO(IR0,:)+RHO(IRP,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        RHO=DRHO/REAL(NOP,KIND=8)
        DEALLOCATE(DRHO)
        RETURN
      END IF
!
!     ==========================================================================
!     == OTHERWISE: SYMMETRIZE IN RECIPROCAL SPACE ON THE WHOLE FFT BOX.      ==
!     == A COMPONENT IS SYMMETRIZED IF ALL ITS IMAGES LIE IN THE BOX (NYQUIST ==
!     == PLANES EXCLUDED); INCOMPLETE STARS ARE LEFT UNCHANGED                ==
!     ==========================================================================
      IF(.NOT.TWARNED) THEN
        WRITE(*,FMT='("SYMMETRY: TRANSLATIONS NOT COMMENSURATE WITH THE ",3I4 &
     &        ," DENSITY GRID; SYMMETRIZING IN RECIPROCAL SPACE")')NR
        TWARNED=.TRUE.
      END IF
      ALLOCATE(FG(0:NR(1)-1,0:NR(2)-1,0:NR(3)-1))
      ALLOCATE(FS(0:NR(1)-1,0:NR(2)-1,0:NR(3)-1))
      ALLOCATE(MINVS(3,3,NOP))
      DO IOP=1,NOP
        A=REAL(IROT(:,:,IOP),KIND=8)
        CALL LIB$INVERTR8(3,A,AINV)
        MINVS(:,:,IOP)=NINT(TRANSPOSE(AINV))    ! H = (M^T)^-1 H'
      ENDDO
      DO IDIM=1,NDIMD
        FG=RESHAPE(CMPLX(RHO(:,IDIM),0.D0,KIND=8),NR)
        CALL SYMMETRY_FFT3D('RTOG',NR,FG)
        FS=FG
        DO I3=0,NR(3)-1
          DO I2=0,NR(2)-1
            DO I1=0,NR(1)-1
              HP=(/I1,I2,I3/)
              DO I=1,3
                HP(I)=MERGE(HP(I)-NR(I),HP(I),2*HP(I).GT.NR(I))
              ENDDO
              IF(ANY(2*ABS(HP).GE.NR)) CYCLE            ! NYQUIST PLANE
              C=(0.D0,0.D0)
              TGRID=.TRUE.
              DO IOP=1,NOP
                H=MATMUL(MINVS(:,:,IOP),HP)
                IF(ANY(2*ABS(H).GE.NR)) THEN
                  TGRID=.FALSE.
                  EXIT
                END IF
                PH=2.D0*PI*DOT_PRODUCT(REAL(H,KIND=8),TAU(:,IOP))
                C=C+FG(MODULO(H(1),NR(1)),MODULO(H(2),NR(2)),MODULO(H(3),NR(3))) &
     &             *CMPLX(COS(PH),SIN(PH),KIND=8)
              ENDDO
              IF(TGRID) FS(I1,I2,I3)=C/REAL(NOP,KIND=8)
            ENDDO
          ENDDO
        ENDDO
        CALL SYMMETRY_FFT3D('GTOR',NR,FS)
        RHO(:,IDIM)=RESHAPE(REAL(FS,KIND=8),(/NRL/))
      ENDDO
      DEALLOCATE(FG,FS,MINVS)
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY_FFT3D(DIR,NR,F)
!     **************************************************************************
!     **  3-D FFT OF A FULL GRID BY 1-D PASSES (LIB_FFTW3 CONVENTIONS):        **
!     **  'GTOR': F(R)=SUM_G F(G)EXP(+IGR); 'RTOG': F(G)=1/N SUM_R F(R)EXP(-IGR)**
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: DIR
      INTEGER(4)  ,INTENT(IN)    :: NR(3)
      COMPLEX(8)  ,INTENT(INOUT) :: F(NR(1),NR(2),NR(3))
      COMPLEX(8)  ,ALLOCATABLE   :: X(:,:),Y(:,:)
      INTEGER(4)                 :: I1,I2,I3
!     **************************************************************************
      ALLOCATE(X(NR(1),NR(2)),Y(NR(1),NR(2)))
      DO I3=1,NR(3)
        X=F(:,:,I3)
        CALL LIB_FFTW3(DIR,NR(1),NR(2),X,Y)
        F(:,:,I3)=Y
      ENDDO
      DEALLOCATE(X,Y)
      ALLOCATE(X(NR(2),NR(1)),Y(NR(2),NR(1)))
      DO I3=1,NR(3)
        X=TRANSPOSE(F(:,:,I3))
        CALL LIB_FFTW3(DIR,NR(2),NR(1),X,Y)
        F(:,:,I3)=TRANSPOSE(Y)
      ENDDO
      DEALLOCATE(X,Y)
      ALLOCATE(X(NR(3),NR(1)),Y(NR(3),NR(1)))
      DO I2=1,NR(2)
        DO I1=1,NR(1)
          X(:,I1)=F(I1,I2,:)
        ENDDO
        CALL LIB_FFTW3(DIR,NR(3),NR(1),X,Y)
        DO I1=1,NR(1)
          F(I1,I2,:)=Y(:,I1)
        ENDDO
      ENDDO
      DEALLOCATE(X,Y)
      RETURN
      END
!     ...1.........2.........3.........4.........5.........6.........7.........8

      SUBROUTINE SYMMETRY$DENMAT(LMNXX,NDIMD,NAT,DENMAT)
!     **************************************************************************
!     **  SYMMETRIZES ONE-CENTER DENSITY MATRICES: D(J)=1/N SUM_G U D(I) U^T  **
!     **  WITH J=MAPAT(I,G) AND U THE ROTATION OF REAL SPHERICAL HARMONICS    **
!     **  (TIMES (-1)^L FOR IMPROPER OPERATIONS) IN EACH PARTIAL-WAVE BLOCK.  **
!     **  TIME REVERSAL (USED BY THE K-POINT REDUCTION) MAKES THE RESULT REAL **
!     **  FOR COLLINEAR SPIN.                                                 **
!     **************************************************************************
      USE SYMMETRY_MODULE, ONLY : NOP,CROT,MAPAT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LMNXX,NDIMD,NAT
      COMPLEX(8),INTENT(INOUT) :: DENMAT(LMNXX,LMNXX,NDIMD,NAT)
      INTEGER(4),PARAMETER     :: LMX=16      ! L<=3
      COMPLEX(8),ALLOCATABLE   :: DSUM(:,:,:,:)
      REAL(8)   ,ALLOCATABLE   :: U(:,:)
      INTEGER(4),ALLOCATABLE   :: ISPECIES(:),LOX(:)
      REAL(8)                  :: YLMROT(LMX,LMX),C(3,3),DET,SGN
      INTEGER(4)               :: IOP,IAT,JAT,ISP,LNX,LMNX,LN,L,LMN,IDIM
!     **************************************************************************
      IF(NOP.LE.1) RETURN
      IF(NDIMD.GT.2) THEN
        CALL ERROR$MSG('SYMMETRIZATION NOT IMPLEMENTED FOR NONCOLLINEAR SPIN')
        CALL ERROR$STOP('SYMMETRY$DENMAT')
      END IF
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      ALLOCATE(DSUM(LMNXX,LMNXX,NDIMD,NAT))
      DSUM=(0.D0,0.D0)
      DO IOP=1,NOP
        C=CROT(:,:,IOP)
        DET=C(1,1)*(C(2,2)*C(3,3)-C(3,2)*C(2,3))-C(1,2)*(C(2,1)*C(3,3)-C(3,1)*C(2,3)) &
     &     +C(1,3)*(C(2,1)*C(3,2)-C(3,1)*C(2,2))
        IF(DET.LT.0.D0) C=-C
        CALL SPHERICAL$ROTATEYLM(LMX,C,YLMROT)
        DO IAT=1,NAT
          JAT=MAPAT(IAT,IOP)
          ISP=ISPECIES(IAT)
          CALL SETUP$ISELECT(ISP)
          CALL SETUP$GETI4('LNX',LNX)
          ALLOCATE(LOX(LNX))
          CALL SETUP$GETI4A('LOX',LNX,LOX)
          CALL SETUP$UNSELECT()
          LMNX=SUM(2*LOX+1)
          ALLOCATE(U(LMNX,LMNX))
          U=0.D0
          LMN=0
          DO LN=1,LNX
            L=LOX(LN)
            SGN=1.D0
            IF(DET.LT.0.D0.AND.MOD(L,2).EQ.1) SGN=-1.D0
            U(LMN+1:LMN+2*L+1,LMN+1:LMN+2*L+1)=SGN*YLMROT(L**2+1:(L+1)**2,L**2+1:(L+1)**2)
            LMN=LMN+2*L+1
          ENDDO
          DO IDIM=1,NDIMD
            DSUM(1:LMNX,1:LMNX,IDIM,JAT)=DSUM(1:LMNX,1:LMNX,IDIM,JAT) &
     &          +MATMUL(U,MATMUL(DENMAT(1:LMNX,1:LMNX,IDIM,IAT),TRANSPOSE(U)))
          ENDDO
          DEALLOCATE(U,LOX)
        ENDDO
      ENDDO
      DENMAT=CMPLX(REAL(DSUM,KIND=8)/REAL(NOP,KIND=8),0.D0,KIND=8)
      DEALLOCATE(DSUM,ISPECIES)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8

      SUBROUTINE SYMMETRY$FORCE(NAT,F)
      USE SYMMETRY_MODULE, ONLY : NOP,CROT,MAPAT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NAT
      REAL(8)   ,INTENT(INOUT) :: F(3,NAT)
      REAL(8)                  :: FS(3,NAT)
      INTEGER(4)               :: IOP,IAT
      IF(NOP.LE.1) RETURN
      FS=0.D0
      DO IOP=1,NOP
        DO IAT=1,NAT
          FS(:,MAPAT(IAT,IOP))=FS(:,MAPAT(IAT,IOP))+MATMUL(CROT(:,:,IOP),F(:,IAT))
        ENDDO
      ENDDO
      F=FS/REAL(NOP,KIND=8)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SYMMETRY$STRESS(S)
      USE SYMMETRY_MODULE, ONLY : NOP,CROT
      IMPLICIT NONE
      REAL(8),INTENT(INOUT) :: S(3,3)
      REAL(8)               :: SS(3,3)
      INTEGER(4)            :: IOP
      IF(NOP.LE.1) RETURN
      SS=0.D0
      DO IOP=1,NOP
        SS=SS+MATMUL(CROT(:,:,IOP),MATMUL(S,TRANSPOSE(CROT(:,:,IOP))))
      ENDDO
      S=SS/REAL(NOP,KIND=8)
      RETURN
      END
!
