!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  OBJECTNAME: AUGMENTATION                                         **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    MOMENTS                                                        **
!**    SPHERE                                                         **
!**    SYNC                                                           **
!**                                                                   **
!**  MODULES USED: MPE                                                **
!**  USES: MPE; SETUP; HYPERFINE, ENERGYLIST, ATOMLIST                **
!**    EXTERNAL1CPOT
!**    CLEBSCH, RADIAL, DFT, GAUSSN                                   **
!**    SELFTEST, TRACE, FILEHANDLER, ERROR                            **
!**                                                                   **
!***********************************************************************
!***********************************************************************
MODULE AUGMENTATION_MODULE
INTEGER(4)  ,PARAMETER :: NE=10
CHARACTER(32)          :: ID(NE)
REAL(8)                :: VAL(NE)
LOGICAL(4)             :: TINI=.FALSE.
!== SCAN0R: ONE-CENTER REAL (ORBITAL-TAU-BASED) KINETIC-ENERGY DENSITY, ===
!== SET BY THE CALLER (AUGMENTATION$SPHERE) IMMEDIATELY BEFORE EACH OF   ===
!== ITS THREE AUGMENTATION_XC CALLS (CORE-ONLY, AE, PS) VIA              ===
!== AUGMENTATION_XC$SETTAU, AND CONSUMED ONCE, IMMEDIATELY, BY           ===
!== AUGMENTATION_XC ITSELF WHEN TSCANDELOC IS ACTIVE -- SAME "SET RIGHT  ===
!== BEFORE THE CALL, CONSUME IMMEDIATELY" PATTERN AS POTENTIAL$SETTAU.   ===
!== RAW (NOT Y0-SCALED) Y_LM COEFFICIENTS, SAME CONVENTION AS RHOIN      ===
!== ITSELF -- LM DIMENSION MATCHES WHATEVER LMRX THE CALLER IS ABOUT TO  ===
!== PASS TO AUGMENTATION_XC FOR THAT SAME CALL (LM=1 FOR THE SPHERICAL-  ===
!== ONLY CASE; LMXCOM, THE PADDED COMMON RHO/TAU RANGE, WHEN NON-        ===
!== SPHERICAL -- SEE AUGMENTATION$SPHERE'S OWN LXCOM CONSTRUCTION).      ===
REAL(8)   ,ALLOCATABLE :: XCTAU(:,:,:)   !(NR,LM,NDIMD)
!== SCAN0R: D(EXCDENS)/D(TAU) FROM THE LAST AUGMENTATION_XC_TSCANDELOC   ===
!== CALL, SAME RAW CONVENTION AND SAME LM RANGE AS XCTAU ABOVE. SET BY   ===
!== AUGMENTATION_XC_TSCANDELOC, READ (AND COPIED OUT BEFORE THE NEXT     ===
!== CALL OVERWRITES IT) BY AUGMENTATION$SPHERE RIGHT AFTER EACH OF ITS   ===
!== AE/PS AUGMENTATION_XC CALLS.                                         ===
REAL(8)   ,ALLOCATABLE :: XCVTAU(:,:,:)  !(NR,LM,NDIMD)
!== ONE-CENTER LDA/GGA: EXACT ANGULAR QUADRATURE (T) OR SECOND-ORDER     ===
!== TAYLOR EXPANSION IN THE NON-SPHERICAL DENSITY (F). SET FROM          ===
!== !CONTROL!DFT:ANGULARXC VIA AUGMENTATION_XC$SETEXACT.                 ===
LOGICAL(4)             :: TXCEXACT=.FALSE.
END MODULE AUGMENTATION_MODULE
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_XC$SETEXACT(TEXACT)
      USE AUGMENTATION_MODULE, ONLY: TXCEXACT
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TEXACT
      TXCEXACT=TEXACT
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_XC$SETTAU(NR_,LMX_,NDIMD_,TAU_)
!     ******************************************************************
!     **  SCAN0R: STORES THE RAW (NOT Y0-SCALED) ONE-CENTER TAU FOR    **
!     **  THE NEXT AUGMENTATION_XC CALL. SEE AUGMENTATION_MODULE'S     **
!     **  XCTAU DECLARATION FOR THE CALLING CONVENTION (LMX_ MUST      **
!     **  MATCH THE LMRX THE CALLER IS ABOUT TO PASS TO AUGMENTATION_XC**
!     **  FOR THIS SAME CALL).                                         **
!     ******************************************************************
      USE AUGMENTATION_MODULE, ONLY: XCTAU
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR_
      INTEGER(4),INTENT(IN) :: LMX_
      INTEGER(4),INTENT(IN) :: NDIMD_
      REAL(8)   ,INTENT(IN) :: TAU_(NR_,LMX_,NDIMD_)
!     ******************************************************************
      IF(ALLOCATED(XCTAU)) DEALLOCATE(XCTAU)
      ALLOCATE(XCTAU(NR_,LMX_,NDIMD_))
      XCTAU=TAU_
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_INI
!     ******************************************************************
!     **  RESETS AUGMENTATION_MODULE                                  **
!     **  (CALLED BY AUGMENTATION_ADD)                                **
!     **                                                              **
!     ******************************************************************
      USE AUGMENTATION_MODULE
      IMPLICIT NONE
      INTEGER(4) :: I
!     ******************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!
      ID(1)='AE1 EXCHANGE-CORRELATION'
      ID(2)='PS1 EXCHANGE-CORRELATION'
      ID(3)='AE1 ELECTROSTATIC'       
      ID(4)='PS1 ELECTROSTATIC'       
      ID(5)='AE1 BACKGROUND'          
      ID(6)='PS1 BACKGROUND'          
      ID(7)='AE1-PS1 KINETIC'         
      ID(8)='LDA+U EXCHANGE'          
      ID(9)='CORE RELAXATION'          
      ID(10)='EXTERNAL 1CENTER POTENTIAL'          
      DO I=1,NE
        VAL(NE)=0.D0
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADD(ID_,VAL_)
!     ******************************************************************
!     **  ADDS DATA INTO THE TEMPORARY STACK FOR SYNCHRONIZATION      **
!     **  WITH THE ENERGYLIST                                         **
!     ******************************************************************
      USE AUGMENTATION_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_
      INTEGER(4)              :: I
!     ******************************************************************
      IF(.NOT.TINI) CALL AUGMENTATION_INI
      DO I=1,NE
        IF(ID_.EQ.ID(I)) THEN
          VAL(I)=VAL(I)+VAL_
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ID NOT RECOGNIZED')
      CALL ERROR$CHVAL('ID ON INPUT   ',ID_)
      DO I=1,NE
        CALL ERROR$CHVAL('POSSIBLE ID ',ID(I))
      ENDDO
      CALL ERROR$STOP('AUGMENTATION_ADD')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION$SYNC
!     ******************************************************************
!     **  SYNCHRONIZES THE TEMPORARY STACK WITH THE ENERGYLIST        **
!     **  AND RESETS THE TEMPORARY STACK                              **
!     **                                                              **
!     **  MUST BE CALLED AFTER AUGMENTATION$SPHERE HAS BEEN           **
!     **  PROCESSED FOR ALL ATOMS                                     **
!     ******************************************************************
      USE AUGMENTATION_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4) :: I
      REAL(8)    :: SVAR
!     ******************************************************************
      TINI=.FALSE.
      CALL MPE$COMBINE('MONOMER','+',VAL)
      DO I=1,NE
        CALL ENERGYLIST$SET(ID(I),VAL(I))
      ENDDO
!     ==  ID(1)='AE1 EXCHANGE-CORRELATION' =============================
!     ==  ID(2)='PS1 EXCHANGE-CORRELATION' =============================
!     ==  ID(3)='AE1 ELECTROSTATIC' ====================================
!     ==  ID(4)='PS1 ELECTROSTATIC' ====================================
!     ==  ID(5)='AE1 BACKGROUND' =======================================
!     ==  ID(6)='PS1 BACKGROUND' =======================================
!     ==  ID(7)='AE1-PS1 KINETIC' ======================================
!     ==  ID(8)='LDA+U EXCHANGE' =======================================
!     ==  ID(9)='CORE RELAXATION' ======================================
!     ==  ID(10)='EXTERNAL 1CENTER POTENTIAL' ===========================
      SVAR=VAL(7)+VAL(1)-VAL(2)+VAL(3)-VAL(4)+VAL(5)-VAL(6)+VAL(8)+VAL(9)+VAL(10)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',SVAR)
      CALL ENERGYLIST$ADD('AE  KINETIC',VAL(7))
      CALL ENERGYLIST$ADD('AE  EXCHANGE-CORRELATION',VAL(1)-VAL(2))
      CALL ENERGYLIST$ADD('    LOCAL CORRELATIONS',VAL(8))
      CALL ENERGYLIST$ADD('AE  ELECTROSTATIC',VAL(3)-VAL(4))
      CALL ENERGYLIST$ADD('BACKGROUND',VAL(5)-VAL(6))
      VAL(:)=0.D0
      RETURN
      END
!
!     .....................................................MOMNTS.......
      SUBROUTINE AUGMENTATION$MOMENTS(ISP,LMNX,DENMAT,LMRX,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ELECTROSTATIC MULTIPOLE MOMENTS FOR THE      **
!     **  COMPENSTATION CHARGE DENSITY                                **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      INTEGER(4),INTENT(IN)  :: ISP
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: LMRX
      COMPLEX(8),INTENT(IN)  :: DENMAT(LMNX,LMNX)
      REAL(8)   ,INTENT(OUT) :: QLM(LMRX)
      REAL(8)   ,ALLOCATABLE :: AERHO(:,:)   !(NR,LMRX)
      REAL(8)   ,ALLOCATABLE :: PSRHO(:,:)   !(NR,LMRX)
      REAL(8)   ,ALLOCATABLE :: AECORE(:)    !(NR)
      REAL(8)   ,ALLOCATABLE :: PSCORE(:)    !(NR)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)   !(NR,LNX)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)   !(NR,LNX)
      INTEGER(4),ALLOCATABLE :: LOX(:)       !(LN)
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)                :: AEZ
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMN1,LMN2
!     ******************************************************************
                          CALL TRACE$PUSH('AUGMENTATION$MOMENTS')
!
!     == PRINT DENSITY MATRIX FOR TEST ========== ======================
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(*,FMT='("TOTAL DENSITY MATRIX")')
        DO LMN1=1,LMNX
          WRITE(*,FMT='(9F10.6)')(DENMAT(LMN1,LMN2),LMN2=1,LMNX)
        ENDDO
      END IF
!     
!     ==================================================================
!     ==   RECEIVE PARTIAL WAVES FROM SETUP OBJECT                    ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
!
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI)
!
!     ==================================================================
!     ==  RECEIVE CORE DENSITY FROM SETUPS OBJECT                     ==
!     ==  RECEIVE ATOMIC NUMBER FROM ATOMLIST                         ==
!     ==================================================================
      ALLOCATE(PSCORE(NR))
      ALLOCATE(AECORE(NR))
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
      CALL SETUP$GETR8A('PSCORE',NR,PSCORE)
      CALL SETUP$GETR8('AEZ',AEZ)
!     
!     ==================================================================
!     ==   CALCULATE ONE CENTER CHARGE DENSITY                        ==
!     ==================================================================
      ALLOCATE(AERHO(NR,LMRX))
      ALLOCATE(PSRHO(NR,LMRX))
      CALL AUGMENTATION_RHO(NR,LNX,LOX,AEPHI,LMNX,DENMAT,LMRX,AERHO)
      CALL AUGMENTATION_RHO(NR,LNX,LOX,PSPHI,LMNX,DENMAT,LMRX,PSRHO)
      DEALLOCATE(AEPHI)
      DEALLOCATE(PSPHI)
      DEALLOCATE(LOX)
!     
!     ==================================================================
!     ==   CALCULATE ONE MULTIPOLE MOMENTS                            ==
!     ==================================================================
      CALL AUGMENTATION_QLM(GID,NR,LMRX,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
      DEALLOCATE(AERHO)
      DEALLOCATE(PSRHO)
      DEALLOCATE(PSCORE)
      DEALLOCATE(AECORE)
      IF(TPR) THEN
        WRITE(*,FMT='(A63)')'MULTIPOLE MOMENTS'
        WRITE(*,FMT='(9E12.5)')QLM(:)
      END IF
      CALL SETUP$UNSELECT()
                   CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION$SPHERE(ISP,IAT,LMNX,NDIMD,DENMAT,EDENMAT &
     &                              ,LMRX,VQLM,RHOB,POTB,DATH,DO)
      USE AUGMENTATION_MODULE, ONLY: XCVTAU
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ** THREE OPTIONS:                                               **
!     **   1 - NON-SPIN POLARIZED (NT)                                **
!     **   2 - SPIN-POLARIZED (NT,NS)                                 **
!     **   4 - NONCOLLINEAR SPIN  (NT,NSX,NSY,NSZ)                    **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: ISP
      INTEGER(4),INTENT(IN)   :: IAT
      INTEGER(4),INTENT(IN)   :: LMNX
      INTEGER(4),INTENT(IN)   :: NDIMD
      COMPLEX(8),INTENT(INOUT):: DENMAT(LMNX,LMNX,NDIMD)
      COMPLEX(8),INTENT(INOUT):: EDENMAT(LMNX,LMNX,NDIMD)
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(INOUT):: VQLM(LMRX)
      REAL(8)   ,INTENT(IN)   :: RHOB ! NEUTRALIZING BACKGROUND DENSITY
!                               RHOB MAY BE SET TO ZERO BY ISOLATE OBJECT
      REAL(8)   ,INTENT(OUT)  :: POTB ! NEG. AV. EL. AUGM. POT.
      COMPLEX(8),INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT)  :: DO(LMNX,LMNX,NDIMD)
      INTEGER(4)              :: GID   ! GRID ID
      REAL(8)   ,ALLOCATABLE  :: R(:)  ! RADIAL GRID
      REAL(8)                 :: AEZ
      INTEGER(4)              :: NR
      INTEGER(4)              :: LNX
      INTEGER(4),ALLOCATABLE  :: LOX(:)
      REAL(8)   ,ALLOCATABLE  :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: AECORE(:)
      REAL(8)   ,ALLOCATABLE  :: PSCORE(:)
      REAL(8)   ,ALLOCATABLE  :: AETAUCORE(:) !(NR) TABULATED CORE KIN.-EN. DENSITY
      REAL(8)   ,ALLOCATABLE  :: PSTAUCORE(:) !(NR) PSEUDIZED COUNTERPART
      REAL(8)                 :: RCSM
      REAL(8)   ,ALLOCATABLE  :: VADD(:)
      REAL(8)   ,ALLOCATABLE  :: DOVER(:,:)
      REAL(8)   ,ALLOCATABLE  :: DTKIN(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DTKIN1(:,:)
      REAL(8)                 :: DELTAH(LMNX,LMNX,NDIMD)  ! SOFTCORE CORRECTION TO DTKIN
      REAL(8)                 :: DELTAO(LMNX,LMNX,NDIMD)  ! SOFTCORE CORRECTION TO DO
      REAL(8)   ,ALLOCATABLE  :: DELTARHO(:,:,:)    ! SOFCORE CORRECTION TO DAERHO
      REAL(8)                 :: DECORE
      REAL(8)   ,ALLOCATABLE  :: AERHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DATP(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DWORK1(:)
      INTEGER(4)              :: IDIM,IR
      INTEGER(4)              :: LM,LMN1
      INTEGER(4)              :: NFILO
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=1
      LOGICAL(4)              :: TBACK
      REAL(8)                 :: DETOT,PSEHARTREE,AEEHARTREE,COREEXC
      REAL(8)                 :: EKINNL,ENL,AEEXC,PSEXC
      CHARACTER(32)           :: ATOM
      REAL(8)                 :: VQLM1(LMRX)
      REAL(8)                 :: QLM(LMRX)
      REAL(8)                 :: PI
      REAL(8)   ,ALLOCATABLE  :: AEHPOT(:,:)
      REAL(8)   ,ALLOCATABLE  :: AEXCPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSHPOT(:,:)
      REAL(8)   ,ALLOCATABLE  :: PSXCPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: AETOTPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSTOTPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: RHO(:,:,:)
      CHARACTER(32)           :: SOFTCORETYPE
!     == SCAN0R: ONE-CENTER REAL TAU FOR AUGMENTATION_XC_TSCANDELOC, SEE  ==
!     == AUGMENTATION_XC$SETTAU / AUGMENTATION_MODULE'S XCTAU DECLARATION ==
      LOGICAL(4)              :: TTSCANDELOC
      INTEGER(4)              :: LXCH,LXTAU,LMXTAU
      INTEGER(4)              :: LRX,LXCOM,LMXCOM
      REAL(8)                 :: Y0
      REAL(8)   ,ALLOCATABLE  :: AETAUL(:,:),PSTAUL(:,:)  !(NR,LMXTAU)
      REAL(8)   ,ALLOCATABLE  :: AETAU(:,:,:),PSTAU(:,:,:)  !(NR,LMXTAU,NDIMD)
      REAL(8)   ,ALLOCATABLE  :: VTAU_AE(:,:,:),VTAU_PS(:,:,:) !(NR,LMXTAU,NDIMD), D(EXCDENS)/D(TAU)
      REAL(8)   ,ALLOCATABLE  :: RHO_PAD(:,:,:),TAU_PAD(:,:,:) !(NR,LMXCOM,NDIMD)
      REAL(8)   ,ALLOCATABLE  :: VXC_PAD(:,:,:)               !(NR,LMXCOM,NDIMD)
      COMPLEX(8),ALLOCATABLE  :: DATH_TAU_AE(:,:),DATH_TAU_PS(:,:) !(LMNX,LMNX)
!     ******************************************************************
                            CALL TRACE$PUSH('AUGMENTATION$SPHERE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
!
!     ==================================================================
!     ==  COLLECT ATOM-TYPE SPECIFIC INFORMATION FROM SETUP OBJECT    ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      CALL SETUP$GETR8('AEZ',AEZ)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI)
      ALLOCATE(AECORE(NR))
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
      ALLOCATE(PSCORE(NR))
      CALL SETUP$GETR8A('PSCORE',NR,PSCORE)
      ALLOCATE(AETAUCORE(NR))
      CALL SETUP$GETR8A('AETAUCORE',NR,AETAUCORE)
      ALLOCATE(PSTAUCORE(NR))
      CALL SETUP$GETR8A('PSTAUCORE',NR,PSTAUCORE)
      CALL SETUP$GETR8('RCSM',RCSM)
      ALLOCATE(VADD(NR))
      CALL SETUP$GETR8A('VADD',NR,VADD)
!
!     ==================================================================
!     ==  EXPAND OVERLAP AND KINETIC ENERGY MATRIX ELEMENTS           ==
!     ==================================================================
      ALLOCATE(DOVER(LNX,LNX))
      CALL SETUP$GETR8A('DO',LNX*LNX,DOVER)
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DOVER,DO)
      DEALLOCATE(DOVER)
!
      ALLOCATE(DTKIN(LMNX,LMNX,NDIMD))
      ALLOCATE(DTKIN1(LNX,LNX))
      CALL SETUP$GETR8A('DEKIN',LNX*LNX,DTKIN1)
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DTKIN1,DTKIN)
      DEALLOCATE(DTKIN1)
      CALL SETUP$UNSELECT()
!
!     ==================================================================
!     ==  SELF TEST                                                   ==
!     ==================================================================
 1000 CONTINUE
      IF(TTEST) THEN
        IF(ITEST.EQ.1) THEN 
          CALL SELFTEST$START('SPHERE',LMNX*LMNX*NDIMD,DENMAT,1.D-3)
          VQLM(:)=0.D0
          DO IDIM=1,NDIMD
            DENMAT(:,:,IDIM)=0.5D0*(DENMAT(:,:,IDIM)+TRANSPOSE(CONJG(DENMAT(:,:,IDIM))))
          ENDDO
        END IF
      ENDIF
!     == PRINT DENSITY MATRIX FOR TEST
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,FMT='("DENMAT FOR IDIM= ",I2)') IDIM
          DO LMN1=1,LMNX
            WRITE(NFILO,FMT='(10F15.6)')DENMAT(LMN1,:,IDIM)
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     ==  CALCULATE 1-CENTER KINETIC ENERGY                           ==
!     ==================================================================
!     == DENSITY MATRIX IS HERMITEAN FOR EACH SPIN DIRECTION
      EKINNL=REAL(SUM(CONJG(DENMAT)*DTKIN),KIND=8)
      CALL AUGMENTATION_ADD('AE1-PS1 KINETIC',EKINNL)
!
!     ==================================================================
!     ==  CALCULATE 1-CENTER CHARGE DENSITY                           ==
!     ==================================================================
      ALLOCATE(AERHO(NR,LMRX,NDIMD))
      ALLOCATE(PSRHO(NR,LMRX,NDIMD))
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,AEPHI &
     &                  ,LMNX,DENMAT(1,1,IDIM),LMRX,AERHO(1,1,IDIM))
        CALL AUGMENTATION_RHO(NR,LNX,LOX,PSPHI &
     &                  ,LMNX,DENMAT(1,1,IDIM),LMRX,PSRHO(1,1,IDIM))
      ENDDO
!
!     ==================================================================
!     ==  SCAN0R: ONE-CENTER REAL TAU FOR AUGMENTATION_XC_TSCANDELOC,  ==
!     ==  IF TSCANDELOC IS ACTIVE. SAME LNX/LOX/AEPHI/PSPHI/DENMAT     ==
!     ==  USED FOR AERHO/PSRHO ABOVE, VIA AUGMENTATION_TAU INSTEAD OF  ==
!     ==  AUGMENTATION_RHO (paw_scantau.f90). FULL LMXTAU=(2*LXCH+1)**2==
!     ==  RANGE IS KEPT (2026-08-28: NO LONGER TRUNCATED TO LM=1) --   ==
!     ==  SEE THE LXCOM PADDING BELOW FOR HOW THIS RECONCILES WITH     ==
!     ==  RHO'S OWN (GENERALLY DIFFERENT) LMRX RANGE, MIRRORING        ==
!     ==  LMTO_SIMPLEDC_NEW_NEW'S OWN LNSPH BRANCH (paw_lmto.f90).     ==
!     ==  THOMAS-FERMI CORE TAU ADDED TO THE T-CHANNEL, LM=1 ONLY,     ==
!     ==  MATCHING HOW AECORE/PSCORE THEMSELVES ARE ADDED ONLY TO      ==
!     ==  RHO(:,1,1) BELOW (SAME "CRITICAL" BUG AVOIDED HERE AS        ==
!     ==  DOCUMENTED FOR DCTYPE=21'S OWN TAU_ALL CONSTRUCTION IN       ==
!     ==  LMTO_SIMPLEDC_NEW_NEW, paw_lmto.f90: A VALENCE-ONLY TAU      ==
!     ==  PAIRED WITH A CORE-INCLUSIVE RHO VIOLATES TAU>=TAU_W AND     ==
!     ==  SILENTLY CLAMPS ALPHA TO 0 NEAR THE NUCLEUS).                ==
!     ==================================================================
      CALL POTENTIAL$GETL4('TSCANDELOC',TTSCANDELOC)
      IF(TTSCANDELOC) THEN
        LXCH=MAXVAL(LOX)
        LXTAU=2*LXCH
        LMXTAU=(LXTAU+1)**2
        LRX=NINT(SQRT(DBLE(LMRX)))-1
        LXCOM=MAX(LRX,LXTAU)
        LMXCOM=(LXCOM+1)**2
        ALLOCATE(AETAUL(NR,LMXTAU))
        ALLOCATE(PSTAUL(NR,LMXTAU))
        ALLOCATE(AETAU(NR,LMXTAU,NDIMD))
        ALLOCATE(PSTAU(NR,LMXTAU,NDIMD))
        DO IDIM=1,NDIMD
          CALL AUGMENTATION_TAU(GID,NR,LNX,LOX,AEPHI,LMNX,DENMAT(1,1,IDIM) &
     &                          ,LXCH,LXTAU,AETAUL)
          AETAU(:,:,IDIM)=AETAUL(:,:)
          CALL AUGMENTATION_TAU(GID,NR,LNX,LOX,PSPHI,LMNX,DENMAT(1,1,IDIM) &
     &                          ,LXCH,LXTAU,PSTAUL)
          PSTAU(:,:,IDIM)=PSTAUL(:,:)
        ENDDO
        DEALLOCATE(AETAUL)
        DEALLOCATE(PSTAUL)
!       == TABULATED (ORBITAL-BASED) CORE TAU (2026-08-28), REPLACING THE  ==
!       == FORMER ON-THE-FLY THOMAS-FERMI STAND-IN -- AETAUCORE/PSTAUCORE  ==
!       == ARE PRECOMPUTED ONCE DURING SETUP CONSTRUCTION (paw_setups.f90,==
!       == SAME RAW Y_LM=1-COEFFICIENT CONVENTION AS AECORE/PSCORE, SO NO ==
!       == Y0 CONVERSION IS NEEDED HERE -- JUST ADD DIRECTLY, MIRRORING   ==
!       == HOW AECORE/PSCORE THEMSELVES ARE ADDED TO RHO(:,1,1) ABOVE).   ==
        AETAU(:,1,1)=AETAU(:,1,1)+AETAUCORE(:)
        PSTAU(:,1,1)=PSTAU(:,1,1)+PSTAUCORE(:)
      END IF
!
!     ================================================================
!     ==  EVALUATE MULTIPOLE MOMENTS                                ==
!     ================================================================
      CALL AUGMENTATION_QLM(GID,NR,LMRX,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
      ALLOCATE(RHO(NR,LMRX,1))
      RHO(:,:,1)=PSRHO(:,:,1)
      RHO(:,1,1)=RHO(:,1,1)+PSCORE(:)   
      CALL AUGMENTATION_ADJUSTVQLM(GID,NR,LMRX,RHO,RCSM,QLM,VQLM1)
      VQLM1(:)=VQLM1(:)+VQLM(:)
      DEALLOCATE(RHO)
!     
!     ================================================================
!     ==  NEW SOFT CORE                                             ==
!     ================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETCH('SOFTCORETYPE',SOFTCORETYPE)
      CALL SETUP$UNSELECT()
      IF(SOFTCORETYPE.EQ.'NONE') THEN
        DECORE=0.D0
      ELSE 
        ALLOCATE(DELTARHO(NR,LMRX,NDIMD))
        DELTAH=0.D0
        DELTAO=0.D0
!        CALL AUGMENTATION_NEWSOFTCORE(SOFTCORETYPE,GID,NR,LMRX,NDIMD,AEZ,LMNX &
!     &          ,DENMAT,EDENMAT,VQLM1,RHOB,DELTAH,DELTAO,DELTARHO,DECORE)
        DTKIN=DTKIN+DELTAH
        DO=DO+DELTAO
        DEALLOCATE(DELTARHO)
      END IF
!
!     =================================================================
!     == HARTREE ENERGY AND POTENTIAL                                ==
!     =================================================================
      ALLOCATE(PSHPOT(NR,LMRX))
      ALLOCATE(AEHPOT(NR,LMRX))
      CALL AUGMENTATION_PSHARTREE(GID,NR,LMRX,PSCORE,PSRHO(:,:,1) &
     &                 ,VADD,RCSM,QLM,VQLM1,RHOB,PSHPOT,PSEHARTREE)
      CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AEZ,AECORE,AERHO(:,:,1) &
     &                               ,VQLM1,RHOB,AEHPOT,AEEHARTREE)
!
!     ================================================================
!     ==   ADD EXCHANGE AND CORRELATION POTENTIAL                   ==
!     ================================================================
      ALLOCATE(AEXCPOT(NR,LMRX,NDIMD))
      ALLOCATE(PSXCPOT(NR,LMRX,NDIMD))
      ALLOCATE(RHO(NR,LMRX,NDIMD))
!     == CORE ONLY EXCHANGE ENERGY (ALWAYS SPHERICAL -- LMRX=1, NO       ==
!     == LXCOM PADDING NEEDED) ========================================
      IF(TTSCANDELOC) CALL AUGMENTATION_XC$SETTAU(NR,1,1,RESHAPE(AETAUCORE,(/NR,1,1/)))
      CALL AUGMENTATION_XC(GID,NR,1,1,AECORE,COREEXC,AEXCPOT(:,1,1))
      AEXCPOT(:,:,:)=0.D0
!     == AE-EXCHANGE ENERGY AND POTENTIAL. WHEN NON-SPHERICAL           ==
!     == (LMXCOM>1), ZERO-PAD RHO/TAU TO THE COMMON LXCOM RANGE BEFORE  ==
!     == THE CALL (MIRRORING LMTO_SIMPLEDC_NEW_NEW'S OWN LNSPH BRANCH), ==
!     == THEN TRUNCATE THE RETURNED POTENTIALS BACK DOWN TO WHAT EACH   ==
!     == CONSUMER NEEDS: AEXCPOT AT THE ORIGINAL LMRX (FOR THE VALUE-   ==
!     == TYPE HAMILTONIAN FOLD VIA AUGMENTATION_EXPECT), VTAU_AE AT     ==
!     == LMXTAU (FOR THE TAU-VARIATIONAL FOLD VIA AUGMENTATION_TAU_     ==
!     == EXPECT, WHICH EXPECTS EXACTLY THAT RANGE).                     ==
!     ==================================================================
      RHO(:,:,:)=AERHO(:,:,:)
      RHO(:,1,1)=RHO(:,1,1)+AECORE(:)
      IF(TTSCANDELOC) THEN
        ALLOCATE(RHO_PAD(NR,LMXCOM,NDIMD),TAU_PAD(NR,LMXCOM,NDIMD))
        ALLOCATE(VXC_PAD(NR,LMXCOM,NDIMD))
        RHO_PAD(:,:,:)=0.D0
        RHO_PAD(:,1:LMRX,:)=RHO(:,:,:)
        TAU_PAD(:,:,:)=0.D0
        TAU_PAD(:,1:LMXTAU,:)=AETAU(:,:,:)
        CALL AUGMENTATION_XC$SETTAU(NR,LMXCOM,NDIMD,TAU_PAD)
        CALL AUGMENTATION_XC(GID,NR,LMXCOM,NDIMD,RHO_PAD,AEEXC,VXC_PAD)
        AEXCPOT(:,:,:)=VXC_PAD(:,1:LMRX,:)
        ALLOCATE(VTAU_AE(NR,LMXTAU,NDIMD))
        VTAU_AE(:,:,:)=XCVTAU(:,1:LMXTAU,:)
        DEALLOCATE(RHO_PAD,TAU_PAD,VXC_PAD)
      ELSE
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,AEEXC,AEXCPOT)
      END IF
      AEEXC=AEEXC-COREEXC
!     == PS-EXCHANGE ENERGY AND POTENTIAL -- SAME LXCOM PADDING PATTERN ==
      RHO(:,:,:)=PSRHO(:,:,:)
      RHO(:,1,1)=RHO(:,1,1)+PSCORE(:)
      IF(TTSCANDELOC) THEN
        ALLOCATE(RHO_PAD(NR,LMXCOM,NDIMD),TAU_PAD(NR,LMXCOM,NDIMD))
        ALLOCATE(VXC_PAD(NR,LMXCOM,NDIMD))
        RHO_PAD(:,:,:)=0.D0
        RHO_PAD(:,1:LMRX,:)=RHO(:,:,:)
        TAU_PAD(:,:,:)=0.D0
        TAU_PAD(:,1:LMXTAU,:)=PSTAU(:,:,:)
        CALL AUGMENTATION_XC$SETTAU(NR,LMXCOM,NDIMD,TAU_PAD)
        CALL AUGMENTATION_XC(GID,NR,LMXCOM,NDIMD,RHO_PAD,PSEXC,VXC_PAD)
        PSXCPOT(:,:,:)=VXC_PAD(:,1:LMRX,:)
        ALLOCATE(VTAU_PS(NR,LMXTAU,NDIMD))
        VTAU_PS(:,:,:)=XCVTAU(:,1:LMXTAU,:)
        DEALLOCATE(RHO_PAD,TAU_PAD,VXC_PAD)
      ELSE
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,PSEXC,PSXCPOT)
      END IF
      DEALLOCATE(RHO)
      IF(TTSCANDELOC) THEN
        DEALLOCATE(AETAU)
        DEALLOCATE(PSTAU)
      END IF
!
!     ================================================================
!     ==   ADD UP POTENTIALS                                        ==
!     ================================================================
      ALLOCATE(AETOTPOT(NR,LMRX,NDIMD))
      ALLOCATE(PSTOTPOT(NR,LMRX,NDIMD))
      AETOTPOT=AEXCPOT
      PSTOTPOT=PSXCPOT
      AETOTPOT(:,:,1)=AETOTPOT(:,:,1)+AEHPOT(:,:)
      PSTOTPOT(:,:,1)=PSTOTPOT(:,:,1)+PSHPOT(:,:)
!
!     ================================================================
!     ==   REPORT ENERGIES                                          ==
!     ================================================================
      CALL AUGMENTATION_ADD('AE1 BACKGROUND',0.D0)
      CALL AUGMENTATION_ADD('PS1 BACKGROUND',0.D0)
      CALL AUGMENTATION_ADD('AE1 EXCHANGE-CORRELATION',AEEXC)
      CALL AUGMENTATION_ADD('PS1 EXCHANGE-CORRELATION',PSEXC)
      CALL AUGMENTATION_ADD('AE1 ELECTROSTATIC',AEEHARTREE)
      CALL AUGMENTATION_ADD('PS1 ELECTROSTATIC',PSEHARTREE)
      CALL AUGMENTATION_ADD('CORE RELAXATION',DECORE)
!
!     =================================================================
!     == AVERAGE ELECTROSTATIC ONE-CENTER POTENTIAL                  ==
!     =================================================================
      ALLOCATE(DWORK1(NR))
      DWORK1(:)=AEHPOT(:,1)-(PSHPOT(:,1)-VADD(:))
      DWORK1(:)=-DWORK1(:)*R(:)**2*SQRT(4.D0*PI)  !SQRT(4*PI)=4*PI*Y_0
      CALL RADIAL$INTEGRAL(GID,NR,DWORK1,POTB)
      DEALLOCATE(DWORK1)
!
!     ==========================================================================
!     ==  SEND DENSITIES TO HYPERFINE-PARAMETER OBJECT                        ==
!     ==========================================================================
      CALL AUGMENTATION_FEEDHYPERFINE(GID,NR,LMRX,NDIMD &
     &         ,IAT,AERHO,AECORE,PSRHO,PSCORE,AEHPOT,PSHPOT)
!     
!     ==========================================================================
!     ==  SEND DENSITIES AND POTENTIALS TO GRAPHICS OBJECT                    ==
!     ==========================================================================
      CALL GRAPHICS$SET1CPOT('HARTREE','AE',IAT,GID,NR,NR,LMRX,AEHPOT)
      CALL GRAPHICS$SET1CPOT('HARTREE','PS',IAT,GID,NR,NR,LMRX,PSHPOT)
      CALL GRAPHICS$SET1CPOT('TOT','AE',IAT,GID,NR,NR,LMRX,AETOTPOT(:,:,1))
      CALL GRAPHICS$SET1CPOT('TOT','PS',IAT,GID,NR,NR,LMRX,PSTOTPOT(:,:,1))
!     
!     ==========================================================================
!     ==  EVALUATE CORE SHIFTS                                                ==
!     ==========================================================================
      CALL CORE_CORESHIFTS(IAT,ISP,GID,NR,LMRX,AETOTPOT)
!
      DEALLOCATE(AEHPOT)
      DEALLOCATE(PSHPOT)
      DEALLOCATE(AEXCPOT)
      DEALLOCATE(PSXCPOT)
!
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,*)'AE POTENTIAL FOR SPIN ',IDIM
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(AETOTPOT(IR,LM,IDIM),LM=1,LMRX)
          ENDDO
        ENDDO
        DO IDIM=1,NDIMD
          WRITE(NFILO,*)'PS POTENTIAL FOR ATOM ',IAT,IDIM
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(PSTOTPOT(IR,LM,IDIM),LM=1,LMRX)
          ENDDO
        ENDDO
      ENDIF
!     
!     ==================================================================
!     ==  SELF TEST FOR HAMILTONIAN AND OVERLAP                       ==
!     ==================================================================
 1001 CONTINUE
      IF(TTEST) THEN
        IF(ITEST.EQ.2) THEN
          CALL SELFTEST$START('SPHERE-EKI',LMNX*LMNX*NDIMD,DENMAT,1.D-3)
        END IF
      END IF
!     
!     ==================================================================
!     ==  EVALUATE KINETIC HAMILTONIAN AND OVERLAP                    ==
!     ==================================================================
      DATH(:,:,:)=CMPLX(DTKIN(:,:,:),0.D0,KIND=8)
!
      IF(TTEST) THEN
        IF(ITEST.EQ.2) THEN
          CALL SELFTEST$END('SPHERE-EKI',LMNX*LMNX*NDIMD,DATH,EKINNL,TBACK)
          IF(TBACK) GOTO 1001
          CALL ERROR$MSG('STOP AFTER SELFTEST')
          CALL ERROR$STOP('SPHERE')
        END IF
      END IF
!     
!     
!     ================================================================
!     ==   ADD POTENTIAL ENERGY TO THE ONE-CENTER HAMILTONIAN       ==
!     ==   DATH =DATH + <AEPHI|AEPOT|AEPHI>-<PSPHI|PSPOT|PSPHI>     ==
!     ================================================================
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL AUGMENTATION_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                        ,AETOTPOT,PSTOTPOT,AEPHI,PSPHI,DATP)
      DATH(:,:,:)=DATH(:,:,:)+CMPLX(DATP(:,:,:),0.D0,KIND=8)
      DEALLOCATE(DATP)
!
!     ================================================================
!     ==  SCAN0R: TAU-VARIATIONAL CONTRIBUTION TO THE ONE-CENTER    ==
!     ==  HAMILTONIAN, IF TSCANDELOC IS ACTIVE -- THE TAU-CHANNEL   ==
!     ==  COUNTERPART OF THE VALUE-TYPE FOLDING JUST ABOVE:         ==
!     ==  DATH=DATH + <AEPHI|VTAU_AE|AEPHI>-<PSPHI|VTAU_PS|PSPHI>   ==
!     ==  (GRADIENT-GRADIENT KERNEL, NOT A PLAIN MULTIPLICATIVE     ==
!     ==  ONE -- SEE AUGMENTATION_TAU_EXPECT, paw_scantau.f90).     ==
!     ==  PER (T,Z) SPIN CHANNEL -- THE KERNEL (XTAB/YTAB, BUILT FROM  ==
!     ==  AEPHI/PSPHI GEOMETRY ONLY) IS SPIN-INDEPENDENT, SO EACH      ==
!     ==  CHANNEL'S CONTRIBUTION IS INDEPENDENT OF THE OTHER, EXACTLY  ==
!     ==  LIKE AUGMENTATION_EXPECT'S OWN PER-ISPIN LOOP ABOVE.         ==
!     ==  VTAU_AE/VTAU_PS ALREADY SPAN THE FULL LMXTAU RANGE (NO       ==
!     ==  ZERO-PADDING NEEDED HERE ANY MORE, UNLIKE THE ORIGINAL       ==
!     ==  SPHERICAL-ONLY FIRST PASS -- SEE THEIR LXCOM-PADDED          ==
!     ==  CONSTRUCTION ABOVE).                                         ==
!     ================================================================
      IF(TTSCANDELOC) THEN
        ALLOCATE(DATH_TAU_AE(LMNX,LMNX))
        ALLOCATE(DATH_TAU_PS(LMNX,LMNX))
        DO IDIM=1,NDIMD
          CALL AUGMENTATION_TAU_EXPECT(GID,NR,LNX,LOX,AEPHI,LMNX,LXCH,LXTAU &
     &                                 ,VTAU_AE(:,:,IDIM),DATH_TAU_AE)
          CALL AUGMENTATION_TAU_EXPECT(GID,NR,LNX,LOX,PSPHI,LMNX,LXCH,LXTAU &
     &                                 ,VTAU_PS(:,:,IDIM),DATH_TAU_PS)
          DATH(:,:,IDIM)=DATH(:,:,IDIM)+DATH_TAU_AE(:,:)-DATH_TAU_PS(:,:)
        ENDDO
        DEALLOCATE(DATH_TAU_AE)
        DEALLOCATE(DATH_TAU_PS)
        DEALLOCATE(VTAU_AE)
        DEALLOCATE(VTAU_PS)
      END IF
!     WRITE(TESTSTRING,FMT='("R8DATH",I2,12(" "))')IAT
!     CALL STOREIT(TESTSTRING,8*LMNX*LMNX*NSPIN,DATH)
!
!
!     ================================================================
!     ==  CI                                                        ==
!     ================================================================
!      CALL PAW_CI(ISP,LMNX,NDIMD,DENMAT)
!
!     ================================================================
!     ==  APPLY EXTERNAL POTENTIAL                                  ==
!     ================================================================
      CALL ATOMLIST$GETCH('NAME',IAT,ATOM)
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATP,DETOT)
      DATH(:,:,:)=DATH(:,:,:)+CMPLX(DATP(:,:,:),0.D0,KIND=8)
      DEALLOCATE(DATP)
      CALL AUGMENTATION_ADD('EXTERNAL 1CENTER POTENTIAL',DETOT)
!     
!     ================================================================
!     ==  SELF-TEST                                                 ==
!     ================================================================
      IF(TTEST) THEN
        IF(ITEST.EQ.1) THEN
          ENL=AEEXC-PSEXC+AEEHARTREE-PSEHARTREE+EKINNL+DETOT
          CALL SELFTEST$END('SPHERE',LMNX*LMNX*NDIMD,DATH,ENL,TBACK)
          IF(TBACK) GOTO 1000
          CALL ERROR$MSG('STOP AFTER SELFTEST')
          CALL ERROR$STOP('SPHERE')
        ENDIF
      END IF
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,FMT='("DATH FOR IDIM= ",I2)') IDIM
          DO LMN1=1,LMNX
            WRITE(NFILO,FMT='(10F15.6)')REAL(DATH(LMN1,:,IDIM))
          ENDDO
        ENDDO
STOP
      END IF
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(AETOTPOT)
      DEALLOCATE(PSTOTPOT)
      DEALLOCATE(AERHO)
      DEALLOCATE(PSRHO)
      DEALLOCATE(LOX)
      DEALLOCATE(AEPHI)
      DEALLOCATE(PSPHI)
      DEALLOCATE(AECORE)
      DEALLOCATE(PSCORE)
      DEALLOCATE(AETAUCORE)
      DEALLOCATE(PSTAUCORE)
      DEALLOCATE(VADD)
      DEALLOCATE(DTKIN)
      DEALLOCATE(R)
                               CALL TRACE$POP()
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_RHO(NR,LNX,LOX,PHI,LMNXX,DENMAT,LMRX,RHOL)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: LMRX
      COMPLEX(8),INTENT(IN) :: DENMAT(LMNXX,LMNXX)
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: RHOL(NR,LMRX)
      INTEGER(4)             :: LMN1,LN1,L1,IM1,LMN2,LN2,L2,IM2
      INTEGER(4)             :: LM1,LM2,LM3
      REAL(8)                :: CG     !CLEBSCH GORDAN COEFFICIENT
      REAL(8)                :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     ==   CALCULATE ONE CENTER CHARGE DENSITY                        ==
!     ==================================================================
!
!     ==================================================================
!     ==  ADD VALENCE CHARGE DENSITY                                  ==
!     ==================================================================
      RHOL(:,:)=0.D0
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
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
!                 == ASSUMES REAL PARTIAL WAVES ==================
                  SVAR=CG*REAL(DENMAT(LMN1,LMN2))
                  RHOL(:,LM3)=RHOL(:,LM3)+SVAR*PHI(:,LN1)*PHI(:,LN2)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_LAPRHO(GID,NR,LMRX,NSPIN,RHO,LAPRHO)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE LAPLACIAN OF A ONE-CENTER DENSITY GIVEN ON   **
!     **  A RADIAL GRID TIMES REAL SPHERICAL HARMONICS, PER (L,M)     **
!     **  CHANNEL:                                                    **
!     **    LAPRHO(:,LM,ISPIN)=D2RHO/DR2+(2/R)*DRHO/DR-L(L+1)/R**2*RHO**
!     **  UNLIKE THE GRADIENT (SEE GRHO IN AUGMENTATION_XC, WHICH     **
!     **  COUPLES DIFFERENT LM CHANNELS THROUGH GAUNT COEFFICIENTS),  **
!     **  THE LAPLACIAN IS DIAGONAL IN (L,M): REAL SPHERICAL          **
!     **  HARMONICS ARE EIGENFUNCTIONS OF THE ANGULAR LAPLACIAN, SO   **
!     **  EACH CHANNEL'S RESULT DEPENDS ONLY ON ITS OWN RADIAL        **
!     **  FUNCTION. NO CROSS-CHANNEL COUPLING IS NEEDED HERE.         **
!     **                                                              **
!     **  ADDED FOR THE SCAN0R / DEORBITALIZED-META-GGA (SCAN-L)      **
!     **  PROJECT: FEEDS THE REDUCED LAPLACIAN "Q" USED BY THE        **
!     **  DEORBITALIZED ISO-ORBITAL INDICATOR (SEE                   **
!     **  DEORB_ALPHA_MODULE, PAW_SCANL.F90).                        **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NSPIN
      REAL(8)   ,INTENT(IN) :: RHO(NR,LMRX,NSPIN)
      REAL(8)   ,INTENT(OUT):: LAPRHO(NR,LMRX,NSPIN)
      REAL(8)                :: R(NR)
      REAL(8)                :: DFDR(NR)
      REAL(8)                :: D2FDR2(NR)
      INTEGER(4)              :: ISPIN,LM,L
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      LAPRHO(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
          CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),DFDR)
          CALL RADIAL$VERLETD2(GID,NR,RHO(:,LM,ISPIN),D2FDR2)
          LAPRHO(:,LM,ISPIN)=D2FDR2(:)+2.D0/R(:)*DFDR(:)
          IF(L.GT.0) THEN
            LAPRHO(2:,LM,ISPIN)=LAPRHO(2:,LM,ISPIN) &
     &                         -DBLE(L*(L+1))/R(2:)**2*RHO(2:,LM,ISPIN)
!           == SAME R(1)-DIVIDE-BY-ZERO GUARD PATTERN USED FOR VRHO IN ==
!           == AUGMENTATION_XC (BELOW, "CORRECT FOR DIVERGENCE AT THE ==
!           == ORIGIN") ========================================================
            LAPRHO(1,LM,ISPIN)=LAPRHO(2,LM,ISPIN)
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_QLM(GID,NR,LMRX &
     &                           ,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE MULTIPOLE MOMENTS OF THE DEVIATION OF             **
!     **  ALL-ELECTRON AND THE PSEUDO DENSITY FROM ONE CENTER         **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)   ,INTENT(IN) :: PSRHO(NR,LMRX)
      REAL(8)   ,INTENT(OUT):: QLM(LMRX)
      REAL(8)               :: DWORK(NR)
      REAL(8)               :: PI,Y0
      REAL(8)               :: RES
      INTEGER(4)            :: LM,L
      REAL(8)               :: R(NR)
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      QLM(:)=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        DWORK(:)=(AERHO(:,LM)-PSRHO(:,LM))*R(:)**(L+2)
        CALL RADIAL$INTEGRAL(GID,NR,DWORK,QLM(LM))
      ENDDO
!
      DWORK(:)=(AECORE(:)-PSCORE(:))*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,DWORK,RES)
      QLM(1)=QLM(1)+RES-AEZ*Y0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
      USE AUGMENTATION_MODULE, ONLY: XCTAU,TXCEXACT
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE                          **
!     **  SPHERICAL CONTRIBUTION OF THE DENSITY UP TO QUADRATIC               **
!     **  ORDER IN THE NON-SPHERICAL CONTRIBUTIONS                            **
!     **                                                                      **
!     **  EXC = EXC(XVAL(L=0)*Y0)                                             **
!     **      + 0.5 * D2[EXC]/D[XVAL(L=0)*Y0]**2 * XVAL(L>0)**2               **
!     **                                                                      **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)              **
!     **  IS AN SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.             **
!     **                                                                      **
!     **  DEPENDECIES:                                                        **
!     **    DFT                                                               **
!     **    TIMING                                                            **
!     **    TRACE                                                             **
!     **                                                                      **
!     **  REMARKS: THE GRADIENTS ARE CORRECT ONLY IF DFT SUPPORTS             **
!     **    THIRD DERIVATIVES OF THE XC-ENERGY                                **
!     **   - WHEN USING SELFTEST ON THIS ROUTINE, THEN                        **
!     **     D(EXC)/DRHO(I)=POT(I)*DEX*R(I)**3                                **
!     **     AND THE VALUES AT LARGE RADII MUST BE SURPRESSED                 **
!     **                                                                      **
!     **  REMARK: FOR A COLLINEAR DENSITY THE ROUTINE GIVES DIFFERENT RESULTS **
!     **          WITH NDIMD=2 AND NDIMD=4 DUE TO THE TAYLOR EXPANSION IN     **
!     **          ANGULAR MOMENTUM EXPANSIONS                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1996 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD      ! CAN BE 1,2,4
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      LOGICAL(4)            :: TSCANL ! SCAN0R: USE SCAN-L EXCHANGE INSTEAD (DFT TYPE=20)
      LOGICAL(4)            :: TTSCANDELOC ! SCAN0R: REAL TAU-SCAN ONE-CENTER CORRECTION
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: Y0
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM
      INTEGER(4)            :: IMAX
      REAL(8)               :: FAC
      REAL(8)               :: CG0LL
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
      REAL(8)   ,PARAMETER  :: XX=1.D0  !BUGFIX 160509 XX=0.5->1.0
!     **************************************************************************
      CALL TRACE$PUSH('AUGMENTATION_XC')
      EXC=0.D0
      VXC(:,:,:)=0.D0
!
!     ==========================================================================
!     ==  SCAN-L (DEORBITALIZED META-GGA) EXCHANGE, IF SELECTED               ==
!     ==  (SCAN0R PROJECT). COMPLETELY BYPASSES THE GGA MACHINERY BELOW --    ==
!     ==  SEE AUGMENTATION_XC_SCANL FOR SCOPE/CAVEATS (ONSITE-SPHERE ONLY,    ==
!     ==  SPHERICAL (LM=1) PART ONLY, APPROXIMATE POTENTIAL FOR THIS FIRST    ==
!     ==  VERSION). NO EFFECT WHATSOEVER UNLESS TSCANL IS SET (DFT TYPE=20).  ==
!     ==========================================================================
      CALL DFT$GETL4('TSCANL',TSCANL)
      IF(TSCANL) THEN
        CALL AUGMENTATION_XC_SCANL(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     ==  REAL (ORBITAL-TAU-BASED) SCAN, ONE-CENTER "AE MINUS PS" CORRECTION, ==
!     ==  IF TSCANDELOC IS ACTIVE (SCAN0R PROJECT). NOTE: TSCANDELOC IS A     ==
!     ==  FLAG ON THE *POTENTIAL* OBJECT (SET VIA !DFT TSCANDELOC=T), NOT    ==
!     ==  THE *DFT* OBJECT -- SEE POTENTIAL_MODULE (paw_potential.f90).      ==
!     ==  WITHOUT THIS BRANCH, THE STANDARD PAW ONE-CENTER RECONSTRUCTION    ==
!     ==  (THIS ROUTINE, CALLED FOR EVERY ATOM/STEP REGARDLESS OF HYBRID/    ==
!     ==  NTBO/DCTYPE SETTINGS) WOULD SILENTLY FALL THROUGH TO PLAIN PBE     ==
!     ==  EVEN THOUGH THE INTERSTITIAL TERM (POTENTIAL_VOFRHO) USES REAL     ==
!     ==  SCAN -- A GENUINE ARCHITECTURAL GAP, NOT A DELIBERATE SCOPE CUT,   ==
!     ==  FOUND BY COMPARING WHICH FILES REFERENCE TSCANDELOC (ONLY          **
!     ==  paw_ioroutines.f90/paw_potential.f90/paw_waves1.f90 DID, NEVER     ==
!     ==  paw_augmentation.f90) -- SEE scan0r NOTES ON THE ONE-CENTER/       ==
!     ==  INTERSTITIAL COMPOSITION GAP.                                      ==
!     ==  FULLY SPIN-RESOLVED (NDIMD=1,2), INCLUDES THE TAU-VARIATIONAL       ==
!     ==  HAMILTONIAN CONTRIBUTION (VIA AUGMENTATION_TAU_EXPECT IN THE       ==
!     ==  CALLER), AND SUPPORTS NON-SPHERICAL (LM>1) ONE-CENTER EVALUATION   ==
!     ==  VIA AUGMENTATION_SCAN_NSPH (paw_scannsph.f90, DIRECT ANGULAR       ==
!     ==  QUADRATURE) WHEN LMRX>1, MIRRORING LMTO_RADXCTAU_NSPH_WITHCUT'S    ==
!     ==  OWN DISPATCH PATTERN (paw_lmto.f90). FALLS BACK TO THE FASTER      ==
!     ==  SPHERICAL-ONLY PATH (AUGMENTATION_SCAN_POINTWISE) WHEN LMRX=1.     ==
!     ==  TAU IS SUPPLIED BY THE CALLER (AUGMENTATION$SPHERE) VIA            ==
!     ==  AUGMENTATION_XC$SETTAU IMMEDIATELY BEFORE EACH CALL -- SEE         ==
!     ==  AUGMENTATION_MODULE'S XCTAU DECLARATION.                          ==
!     ==========================================================================
      CALL POTENTIAL$GETL4('TSCANDELOC',TTSCANDELOC)
      IF(TTSCANDELOC) THEN
        CALL AUGMENTATION_XC_TSCANDELOC(GID,NR,LMRX,NDIMD,RHOIN,XCTAU,EXC,VXC)
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     ==  EXACT ANGULAR QUADRATURE INSTEAD OF THE TAYLOR EXPANSION BELOW      ==
!     ==========================================================================
      IF(TXCEXACT.AND.LMRX.GT.1.AND.NDIMD.LE.2) THEN
        CALL AUGMENTATION_XC_NSPH(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CG0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  OBTAIN SPIN DENSITY                                                 ==
!     ==========================================================================
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      ALLOCATE(RHO(NR,LMRX,NSPIN))
      ALLOCATE(GRHO(NR,LMRX,NSPIN))
      ALLOCATE(VRHO(NR,LMRX,NSPIN))
      RHO(:,:,1)=RHOIN(:,:,1)
      IF(NDIMD.EQ.2) THEN
        RHO(:,:,2)=RHOIN(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == HERE WE NEED TO CALCULATE THE ABSOLUTE VALUE OF THE SPIN DENSITY ====
!       == IN AN ANGULAR MOMENTUM EXPANSION. THIS IS ONLY POSSIBLE APPROXIMATELY
!       == USING A TAYLOR EXPANSION ABOUT THE SPHERICAL PART OF THE SQUARE =====
!       == OF THE SPIN DENSITY. ================================================
        VRHO(:,:,:)=0.D0
        CALL AUGMENTATION_NCOLLTRANS(GID,'RHO',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF
!
!     == IMAX ALLOWS TO RESTRICT SOME LOOPS (1:5) TO (1:IMAX)
      IF(TGRA) THEN
        IF(NSPIN.EQ.2) THEN; IMAX=5; ELSE; IMAX=3; END IF
      ELSE 
        IF(NSPIN.EQ.2) THEN; IMAX=2; ELSE; IMAX=1; END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE RADIAL GRADIENT OF THE DENSITY                            ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE GRADIENTS')
      IF(TGRA) THEN
        GRHO(:,:,:)=0.D0
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
          ENDDO
        ENDDO
      ELSE
        GRHO(:,:,:)=0.D0
      END IF
!
!     ==========================================================================
!     ==  DEFINE VECTOR (RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS)             ==
!     ==========================================================================
      XVAL(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
          XVAL(:,ISPIN,LM)=RHO(:,LM,ISPIN)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN          ! THIS LOOP PUTS T,T->3; S,S->4 ;T,S->5
          DO ISPIN2=ISPIN1,1,-1    ! AND ASSURES CONSISTENCY WITH NSPIN
            II=II+1
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              XVAL(:,II,1)=XVAL(:,II,1) &
        &         +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
        &         +XX*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                      +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY                 ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE DFT')
      WORK1(:)=0.D0
      XDER(:,:,:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ================================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ================================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)
!       == NOW CALCULATE ENERGY DENSITY AND DERIAVTIVES ========================
        WORK1(IR)=FOURPI*EXC1
        XDER(IR,:,1)  =FOURPI*VXC5(:)*Y0
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              WORK1(IR)=WORK1(IR) &
       &              +0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &              +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,WORK1(:)*R(:)**2,EXC)
!
!     ==========================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                             ==
!     ==========================================================================
      ALLOCATE(VGRHO(NR,LMRX,NSPIN))
      VRHO(:,:,:)=0.D0
      VGRHO(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          VRHO(:,LM,ISPIN)=XDER(:,ISPIN,LM)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN
          DO ISPIN2=ISPIN1,1,-1
            II=II+1
!           == FIRST RESOLVE XVAL(:,II,1) ======================================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
!             == THE FOLLOWING LINES DIVIDE BY ZERO IF R=0. ====================
!             == THE VALUES ARE OVERWRITTEN AT THE END OF THE ROUTINE... =======
              VRHO(:,LM,ISPIN1)  =VRHO(:,LM,ISPIN1) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN2)
              VRHO(:,LM,ISPIN2)  =VRHO(:,LM,ISPIN2) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN1) =VGRHO(:,LM,ISPIN1) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN2)
              VGRHO(:,LM,ISPIN2) =VGRHO(:,LM,ISPIN2) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN1)
            ENDDO
!           == NOW RESOLVE XVAL(:,II,LM) =======================================
            DO LM=2,LMRX
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              VGRHO(:,1,ISPIN1) =VGRHO(:,1,ISPIN1) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2) =VGRHO(:,1,ISPIN2) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN2)=VGRHO(:,LM,ISPIN2) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
              VGRHO(:,LM,ISPIN1)=VGRHO(:,LM,ISPIN1) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
            ENDDO
          ENDDO
        ENDDO               
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==  V = V -1/R**2 D/DR [ R**2 VGRHO ]                                   ==
!     ==  V = V -[2/R VGRHO+ D/DR VGRHO ]                                     ==
!     ==========================================================================
      IF(TGRA) THEN
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
!           == FIRST ALTERNATIVE
!           CALL RADIAL$DERIVE(GID,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(GID,NR,WORK2,WORK1)
            WORK1(2:)=WORK1(2:)/R(2:)**2
            WORK(1)=WORK(2)  ! AVOID DIVIDE BY ZERO
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(VGRHO)
!
!     ==========================================================================
!     ==   CORRECT FOR DIVERGENCE AT THE ORIGIN:                              ==
!     ==   IF A SHIFTED LOGARITHMIC GRID IS USED THE FIRST GRID POINT         ==
!     ==   IS MESSED UP BECAUSE OF FACTORS 1/R                                ==
!     ==========================================================================
      IF(R(1).LT.1.D-5) THEN
        VRHO(1,:,:)=VRHO(2,:,:)
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==========================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=VRHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        CALL AUGMENTATION_NCOLLTRANS(GID,'POT',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF     
      DEALLOCATE(RHO)
      DEALLOCATE(GRHO)
      DEALLOCATE(VRHO)
                      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_XC_PREV(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE                          **
!     **  SPHERICAL CONTRIBUTION OF THE DENSITY UP TO QUADRATIC               **
!     **  ORDER IN THE NON-SPHERICAL CONTRIBUTIONS                            **
!     **                                                                      **
!     **  EXC = EXC(XVAL(L=0)*Y0)                                             **
!     **      + 0.5 * D2[EXC]/D[XVAL(L=0)*Y0]**2 * XVAL(L>0)**2               **
!     **                                                                      **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)              **
!     **  IS AN SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.             **
!     **                                                                      **
!     **  DEPENDECIES:                                                        **
!     **    DFT                                                               **
!     **    TIMING                                                            **
!     **    TRACE                                                             **
!     **                                                                      **
!     **  REMARKS: THE GRADIENTS ARE CORRECT ONLY IF DFT SUPPORTS             **
!     **    THIRD DERIVATIVES OF THE XC-ENERGY                                **
!     **   - WHEN USING SELFTEST ON THIS ROUTINE, THEN                        **
!     **     D(EXC)/DRHO(I)=POT(I)*DEX*R(I)**3                                **
!     **     AND THE VALUES AT LARGE RADII MUST BE SURPRESSED                 **
!     **                                                                      **
!     **  REMARK: FOR A COLLINEAR DENSITY THE ROUTINE GIVES DIFFERENT RESULTS **
!     **          WITH NDIMD=2 AND NDIMD=4 DUE TO THE TAYLOR EXPANSION IN     **
!     **          ANGULAR MOMENTUM EXPANSIONS                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1996 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD      ! CAN BE 1,2,4
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: B(:,:)    
      REAL(8)               :: B0INV(NR) 
      REAL(8)   ,ALLOCATABLE:: C(:,:)    
      REAL(8)   ,ALLOCATABLE:: VB(:,:)    
      REAL(8)   ,ALLOCATABLE:: VC(:,:)    
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: Y0
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM,LM1,LM2,LM3
      INTEGER(4)            :: IMAX
      REAL(8)               :: FAC
      REAL(8)               :: CG0LL
      REAL(8)               :: CG
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
      REAL(8)               :: WORK3(NR)
      REAL(8)   ,PARAMETER  :: R8SMALL=1.D-20
!     **************************************************************************
      CALL TRACE$PUSH('AUGMENTATION_XC')
      EXC=0.D0
      VXC(:,:,:)=0.D0
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CG0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  OBTAIN SPIN DENSITY                                                 ==
!     ==========================================================================
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      ALLOCATE(RHO(NR,LMRX,NSPIN))
      ALLOCATE(GRHO(NR,LMRX,NSPIN))
      ALLOCATE(VRHO(NR,LMRX,NSPIN))
      RHO(:,:,1)=RHOIN(:,:,1)
      IF(NDIMD.EQ.2) THEN
        RHO(:,:,2)=RHOIN(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == HERE WE NEED TO CALCULATE THE ABSOLUTE VALUE OF THE SPIN DENSITY ====
!       == IN AN ANGULAR MOMENTUM EXPANSION. THIS IS ONLY POSSIBLE APPROXIMATELY
!       == USING A TAYLOR EXPANSION ABOUT THE SPHERICAL PART OF THE SQUARE =====
!       == OF THE SPIN DENSITY. ================================================
        ALLOCATE(B(NR,LMRX))
!       == EVALUATE B, THE SQUARE OF THE SPIN DENSITY ==========================
        B(:,:)=0.D0
        DO LM2=1,LMRX
          DO LM3=LM2,LMRX
            WORK(:)=RHOIN(:,LM2,2)*RHOIN(:,LM3,2) &
     &             +RHOIN(:,LM2,3)*RHOIN(:,LM3,3) &
     &             +RHOIN(:,LM2,4)*RHOIN(:,LM3,4)
            IF(LM2.NE.LM3)WORK(:)=2.D0*WORK(:)
            DO LM1=1,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              B(:,LM1)=B(:,LM1)+CG*WORK(:)
            ENDDO
          ENDDO
        ENDDO
        B0INV=1.D0/(B(:,1)+R8SMALL)
!       ==  TRANSFORM B INTO C=B/(2*B0*Y0) =====================================
!       ==  THE SPHERICAL PART OF C REAMAINS ZERO ==============================
        ALLOCATE(C(NR,LMRX))
        C(:,:)=0.D0
        WORK(:)=0.5D0*B0INV(:)/Y0
        DO LM1=2,LMRX
          C(:,LM1)=B(:,LM1)*WORK(:)
        ENDDO
!       == CALCULATE SPIN DENSITY ==============================================
        RHO(:,:,2)=0.D0
        DO LM2=2,LMRX
          DO LM3=LM2,LMRX
            WORK(:)=C(:,LM2)*C(:,LM3)
            IF(LM2.NE.LM3) WORK(:)=2.D0*WORK(:)
            DO LM1=1,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              RHO(:,LM1,2)=RHO(:,LM1,2)+CG*WORK(:)
            ENDDO
          ENDDO
        ENDDO
        RHO(:,:,2)=-0.5D0*RHO(:,:,2)
        RHO(:,1,2)=RHO(:,1,2)+1.D0/Y0
        DO LM1=2,LMRX
          RHO(:,LM1,2)=RHO(:,LM1,2)+C(:,LM1)
        ENDDO
        WORK(:)=SQRT(B(:,1)*Y0)
        DO LM1=1,LMRX
          RHO(:,LM1,2)=WORK(:)*RHO(:,LM1,2)
        ENDDO
      END IF
!
!     == IMAX ALLOWS TO RESTRICT SOME LOOPS (1:5) TO (1:IMAX)
      IF(TGRA) THEN
        IF(NSPIN.EQ.2) THEN; IMAX=5; ELSE; IMAX=3; END IF
      ELSE 
        IF(NSPIN.EQ.2) THEN; IMAX=2; ELSE; IMAX=1; END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE RADIAL GRADIENT OF THE DENSITY                            ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE GRADIENTS')
      IF(TGRA) THEN
        GRHO(:,:,:)=0.D0
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
          ENDDO
        ENDDO
      ELSE
        GRHO(:,:,:)=0.D0
      END IF
!
!     ==========================================================================
!     ==  DEFINE VECTOR (RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS)             ==
!     ==========================================================================
      XVAL(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
          XVAL(:,ISPIN,LM)=RHO(:,LM,ISPIN)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN          ! THIS LOOP PUTS T,T->3; S,S->4 ;T,S->5
          DO ISPIN2=ISPIN1,1,-1    ! AND ASSURES CONSISTENCY WITH NSPIN
            II=II+1
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              XVAL(:,II,1)=XVAL(:,II,1) &
        &         +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
        &         +0.5D0*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                      +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY                 ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE DFT')
      WORK1(:)=0.D0
      XDER(:,:,:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ================================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ================================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)
!       == NOW CALCULATE ENERGY DENSITY AND DERIAVTIVES ========================
        WORK1(IR)=FOURPI*EXC1
        XDER(IR,:,1)  =FOURPI*VXC5(:)*Y0
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              WORK1(IR)=WORK1(IR) &
       &              +0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &              +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,WORK1(:)*R(:)**2,EXC)
!
!     ==========================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                             ==
!     ==========================================================================
      ALLOCATE(VGRHO(NR,LMRX,NSPIN))
      VRHO(:,:,:)=0.D0
      VGRHO(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          VRHO(:,LM,ISPIN)=XDER(:,ISPIN,LM)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN
          DO ISPIN2=ISPIN1,1,-1
            II=II+1
!           == FIRST RESOLVE XVAL(:,II,1) ======================================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              VRHO(:,LM,ISPIN1)  =VRHO(:,LM,ISPIN1) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN2)
              VRHO(:,LM,ISPIN2)  =VRHO(:,LM,ISPIN2) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN1) =VGRHO(:,LM,ISPIN1) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN2)
              VGRHO(:,LM,ISPIN2) =VGRHO(:,LM,ISPIN2) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN1)
            ENDDO
!           == NOW RESOLVE XVAL(:,II,LM) =======================================
            DO LM=2,LMRX
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              VGRHO(:,1,ISPIN1) =VGRHO(:,1,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2) =VGRHO(:,1,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN2)=VGRHO(:,LM,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
              VGRHO(:,LM,ISPIN1)=VGRHO(:,LM,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
            ENDDO
          ENDDO
        ENDDO               
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==  V = V -1/R**2 D/DR [ R**2 VGRHO ]                                   ==
!     ==  V = V -[2/R VGRHO+ D/DR VGRHO ]                                     ==
!     ==========================================================================
      IF(TGRA) THEN
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
!           == FIRST ALTERNATIVE
!           CALL RADIAL$DERIVE(GID,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(GID,NR,WORK2,WORK1)
            WORK1(:)=WORK1(:)/R(:)**2
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(VGRHO)
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==========================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=VRHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == CALCULATE VC=DE/DC===================================================
        ALLOCATE(VC(NR,LMRX))
        VC(:,:)=0.D0
!       ---- L.NEQ.0 -----------------------------------------------------------
        DO LM2=1,LMRX
          DO LM3=2,LMRX
            WORK(:)=VRHO(:,LM2,2)*C(:,LM3)
            DO LM1=2,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              VC(:,LM1)=VC(:,LM1)+CG*WORK(:)
            ENDDO
          ENDDO
        ENDDO
        WORK(:)=SQRT(B(:,1)*Y0)
        DO LM1=2,LMRX
          VC(:,LM1)=WORK(:)*(VRHO(:,LM1,2)-VC(:,LM1))
        ENDDO
!       ----- L=0 --------------------------------------------------------------
        WORK(:)=0.5D0*B0INV
        DO LM1=1,LMRX
          VC(:,1)=VC(:,1)+WORK(:)*VRHO(:,LM1,2)*RHO(:,LM1,2)
        ENDDO 
!       == CALCULATE VB=DE/DB ==================================================
        ALLOCATE(VB(NR,LMRX))
        VB(:,:)=0.D0
!       ----- L=0 --------------------------------------------------------------
        DO LM1=2,LMRX
          VB(:,1)=VB(:,1)-VC(:,LM1)*C(:,LM1)
        ENDDO
        VB(:,1)=VC(:,1)-B0INV(:)*VB(:,1)
!       ---- L.NEQ.0 -----------------------------------------------------------
        WORK(:)=0.5D0*B0INV(:)/Y0
        DO LM1=2,LMRX
          VB(:,LM1)=VC(:,LM1)*WORK(:)
        ENDDO      
!       == TRANSFORM VB=DE/DB INTO POTENTIAL FOR SPIN DENSITY ==================
        DO LM1=1,LMRX
          VB(:,LM1)=2.D0*VB(:,LM1)
          VXC(:,LM1,2)=0.D0
          VXC(:,LM1,3)=0.D0
          VXC(:,LM1,4)=0.D0
        ENDDO
        DO LM2=1,LMRX
          DO LM3=1,LMRX
            WORK1(:)=VB(:,LM2)*RHOIN(:,LM3,2)
            WORK2(:)=VB(:,LM2)*RHOIN(:,LM3,3)
            WORK3(:)=VB(:,LM2)*RHOIN(:,LM3,4)
            DO LM1=1,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              VXC(:,LM1,2)=VXC(:,LM1,2)+CG*WORK1(:)
              VXC(:,LM1,3)=VXC(:,LM1,3)+CG*WORK2(:)
              VXC(:,LM1,4)=VXC(:,LM1,4)+CG*WORK3(:)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(B)
        DEALLOCATE(VB)
        DEALLOCATE(C)
        DEALLOCATE(VC)
      END IF     
      DEALLOCATE(RHO)
      DEALLOCATE(GRHO)
      DEALLOCATE(VRHO)
!
!     ==========================================================================
!     ==   CORRECT FOR DIVERGENCE AT THE ORIGIN:                              ==
!     ==   IF A SHIFTED LOGARITHMIC GRID IS USED THE FIRST GRID POINT         ==
!     ==   IS MESSED UP BECAUSE OF FACTORS 1/R                                ==
!     ==========================================================================
      IF(R(1).LT.1.D-5) THEN
        VXC(1,:,:)=VXC(2,:,:)
      END IF
                      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS(GID,ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO SECOND ORDER. THUS THE RESULTS      **
!     **  FOR COLLINEAR DENSITIES DO NOT AGREE WITH THOSE OF COLLINEAR        **
!     **  TREATED AS NONCOLLINEAR ONES                                        **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)    :: GID
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)                    :: Q(NR)
      REAL(8)                    :: VQ(NR)
      REAL(8)                    :: A(NR,LMRX,3)
      REAL(8)                    :: VA(NR,LMRX,3)
      REAL(8)                    :: P(NR,LMRX)
      REAL(8)                    :: VP(NR,LMRX)
      REAL(8)                    :: S(NR,LMRX)
      REAL(8)                    :: VS(NR,LMRX)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: CG   ! GAUNT COEFFICIENT
      REAL(8)     ,PARAMETER     :: SMALL=(1.D-2)**2
      INTEGER(4)                 :: ISIG,LM,LM1,LM2,LM3
      LOGICAL(4)  ,PARAMETER     :: TOTHER=.TRUE.
!     **************************************************************************
      IF(TOTHER) THEN
        CALL AUGMENTATION_NCOLLTRANS_OTHER3(GID,ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z_OTHER.DAT',GID,NR,LMRX,RHO4(:,:,4))
!!$CALL AUGMENTATION_WRITEPHI('RHO2_Z_OTHER.DAT',GID,NR,LMRX,RHO2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT2_Z_OTHER.DAT',GID,NR,LMRX,POT2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT4_Z_OTHER.DAT',GID,NR,LMRX,POT4(:,:,4))
        RETURN
      END IF
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT AUXILIARY VARIABLES                                         ==
!     ==========================================================================
      Q(:)=SQRT(SMALL+(RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)*Y0**2)
!
!     == CONSTRUCT A=SIGMA/Q ===================================================
      DO ISIG=1,3
        DO LM=1,LMRX
          A(:,LM,ISIG)=RHO4(:,LM,ISIG+1)/Q(:)
        ENDDO
      ENDDO
!
!     ==  CONSTRUCT P=A^2 ======================================================
      P(:,:)=0.D0
      DO LM1=1,LMRX
        DO LM2=LM1,LMRX
          AUX(:)=A(:,LM1,1)*A(:,LM2,1) &
     &          +A(:,LM1,2)*A(:,LM2,2) &
     &          +A(:,LM1,3)*A(:,LM2,3)
          IF(LM1.NE.LM2)AUX(:)=2.D0*AUX(:)
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            P(:,LM3)=P(:,LM3)+CG*AUX(:)
          ENDDO
        ENDDO
      ENDDO
!
!     == CONSTRUCT S ===========================================================
      S(:,:)=0.D0
      DO LM1=1,LMRX
        S(:,LM1)=S(:,LM1)+0.75D0*P(:,LM1)
        DO LM2=LM1,LMRX
          AUX(:)=P(:,LM1)*P(:,LM2)
          IF(LM1.NE.LM2)AUX(:)=2.D0*AUX(:)
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            S(:,LM3)=S(:,LM3)-0.125D0*CG*AUX(:)
          ENDDO
        ENDDO
      ENDDO
      S(:,1)=S(:,1)+0.375D0/Y0
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
      RHO2(:,:,1)=RHO4(:,:,1)
      RHO2(:,1,2)=0.D0
      DO LM=1,LMRX
        RHO2(:,LM,2)=RHO2(:,LM,2)+Q(:)*S(:,LM)
      ENDDO
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
!
!     == CONSTRUCT VQ ==========================================================
      VQ(:)=0.D0
      DO LM=1,LMRX
        VQ(:)=VQ(:)+POT2(:,LM,2)*S(:,LM)
        VS(:,LM)=POT2(:,LM,2)*Q(:)
      ENDDO
!
!     == CONSTRUCT VP ==========================================================
      DO LM1=1,LMRX
        VP(:,LM1)=0.75D0*VS(:,LM1)
        DO LM2=1,LMRX
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            VP(:,LM1)=VP(:,LM1)-0.25D0*CG*VS(:,LM2)*P(:,LM3)
          ENDDO
        ENDDO
      ENDDO
!
!     == CONSTRUCT VA ==========================================================
      VA(:,:,:)=0.D0
      DO ISIG=1,3
        DO LM1=1,LMRX
          DO LM2=1,LMRX
            DO LM3=1,LMRX
              CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
              VA(:,LM1,ISIG)=VA(:,LM1,ISIG)+2.D0*CG*VP(:,LM2)*A(:,LM3,ISIG)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     == CONSTRUCT POTENTIAL FOR NON-COLLINEAR SPIN DENSITY ====================
      AUX(:)=0.D0
      DO ISIG=1,3
        DO LM=1,LMRX
          AUX(:)=AUX(:)+VA(:,LM,ISIG)*A(:,LM,ISIG)
        ENDDO
      ENDDO
      AUX(:)=(VQ(:)-AUX(:)/Q(:))*Y0**2
!
      POT4(:,:,:)=0.D0
      DO ISIG=1,3
        POT4(:,1,ISIG+1)=POT4(:,1,ISIG)+AUX(:)*A(:,1,ISIG)
        DO LM1=1,LMRX
          POT4(:,LM1,ISIG+1)=POT4(:,LM1,ISIG+1)+VA(:,LM1,ISIG)/Q(:)
        ENDDO
      ENDDO
!
!     == ADD TOTAL POTENTIAL ===================================================
      POT4(:,:,1)=POT2(:,:,1)
!!$CALL AUGMENTATION_WRITEPHI('Q.DAT',GID,NR,1,Q)
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z.DAT',GID,NR,LMRX,RHO4(:,:,4))
!!$CALL AUGMENTATION_WRITEPHI('RHO2_Z.DAT',GID,NR,LMRX,RHO2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT2_Z.DAT',GID,NR,LMRX,POT2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT4_Z.DAT',GID,NR,LMRX,POT4(:,:,4))
!!$STOP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS_OTHER1(ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO FIRST ORDER ONLY. IT IS CONSTRUCTED **
!     **  SUCH THAT THE RESULTS FOR COLLINEAR DENSITIES DO AGREE WITH         **
!     **  THOSE OF COLLINEAR DENSITIES TREATED AS NONCOLLINEAR ONES           **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)                    :: A(NR,LMRX,3)
      REAL(8)                    :: VA(NR,LMRX,3)
      REAL(8)                    :: Q(NR)
      REAL(8)                    :: VQ(NR)
      REAL(8)                    :: PI,Y0
      REAL(8)     ,PARAMETER     :: SMALL=1.D-4 ! DO NOT CHOOSE TOO SMALL!
      INTEGER(4)                 :: ISIG,LM
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT AUXILIARY VARIABLES                                         ==
!     ==========================================================================
      Q(:)=SQRT(SMALL+RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)
      DO ISIG=1,3
        DO LM=1,LMRX
          A(:,LM,ISIG)=RHO4(:,LM,ISIG+1)/Q(:)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
      RHO2(:,:,1)=RHO4(:,:,1)
      RHO2(:,1,2)=Q(:)
      DO LM=2,LMRX
        RHO2(:,LM,2)=0.D0
        DO ISIG=1,3
          RHO2(:,LM,2)=RHO2(:,LM,2)+A(:,1,ISIG)*A(:,LM,ISIG)
        ENDDO
        RHO2(:,LM,2)=RHO2(:,LM,2)*Q(:)
      ENDDO
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
!     -- CALCULATE VQ ----------------------------------------------------------
      VQ(:)=0.D0
      DO LM=1,LMRX
        VQ(:)=VQ(:)+POT2(:,LM,2)*RHO2(:,LM,2)/Q(:)
      ENDDO
!     -- CALCULATE VA ----------------------------------------------------------
      DO ISIG=1,3
        VA(:,1,ISIG)=0.D0
        DO LM=2,LMRX
          VA(:,1,ISIG)=VA(:,1,ISIG)+POT2(:,LM,2)*A(:,LM,ISIG)
          VA(:,LM,ISIG)=POT2(:,LM,2)*A(:,1,ISIG)
        ENDDO
      ENDDO
!     -- NEW VQ ----------------------------------------------------------------
      DO ISIG=1,3
        DO LM=1,LMRX
          VQ(:)=VQ(:)-VA(:,LM,ISIG)*A(:,LM,ISIG)
        ENDDO
      ENDDO
!     --------------------------------------------------------------------------
      POT4(:,:,1)=POT2(:,:,1)
      DO ISIG=1,3
        DO LM=1,LMRX
          POT4(:,LM,ISIG+1)=VA(:,LM,ISIG)
        ENDDO
        POT4(:,1,ISIG+1)=POT4(:,1,ISIG+1)+VQ(:)*A(:,1,ISIG)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS_OTHER2(ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO FIRST ORDER ONLY. IT IS CONSTRUCTED **
!     **  SUCH THAT THE RESULTS FOR COLLINEAR DENSITIES DO AGREE WITH         **
!     **  THOSE OF COLLINEAR DENSITIES TREATED AS NONCOLLINEAR ONES           **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)                    :: PI,Y0
      REAL(8)     ,PARAMETER     :: SMALL=1.D-4 ! DO NOT CHOOSE TOO SMALL!
      INTEGER(4)                 :: ISIG,LM
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
      RHO2(:,:,1)=RHO4(:,:,1)
      RHO2(:,1,2)=SQRT(RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)
      DO LM=2,LMRX
        RHO2(:,LM,2)=Y0*(RHO4(:,1,2)*RHO4(:,LM,2) &
     &                  +RHO4(:,1,3)*RHO4(:,LM,3) &
     &                  +RHO4(:,1,4)*RHO4(:,LM,4)) 
      ENDDO
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
      POT4(:,:,1)=POT2(:,:,1)
      DO ISIG=2,4
        POT4(:,1,ISIG)=POT2(:,1,2)*RHO4(:,1,ISIG)/(SMALL+RHO2(:,1,2))
      ENDDO
      DO LM=2,LMRX
        DO ISIG=2,4
           POT4(:,LM,ISIG)=POT2(:,LM,2)*RHO4(:,1,ISIG)*Y0
           POT4(:,1,ISIG)=POT4(:,1,ISIG)+POT2(:,LM,2)*RHO4(:,LM,ISIG)*Y0
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS_OTHER3(GID,ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO FIRST ORDER ONLY. IT IS CONSTRUCTED **
!     **  SUCH THAT THE RESULTS FOR COLLINEAR DENSITIES DO AGREE WITH         **
!     **  THOSE OF COLLINEAR DENSITIES TREATED AS NONCOLLINEAR ONES           **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: GID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)                    :: ARRAY(NR)
      REAL(8)                    :: A(NR,LMRX,3)
      REAL(8)                    :: B(NR,LMRX,3)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: PI,Y0,C0LL,CG
      REAL(8)     ,PARAMETER     :: SMALL=1.D-6 ! DO NOT CHOOSE TOO SMALL!
      INTEGER(4)                 :: ISIG,LM,LM1,LM2,LM3
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=Y0
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
!     == TOTAL DENSITY =========================================================
      RHO2(:,:,1)=RHO4(:,:,1)
!
!     == SQUARE OF THE SPIN DENSITY=============================================
      DO LM1=1,LMRX
        AUX(:)=0.D0
        DO LM2=1,LMRX
          DO LM3=LM2,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            IF(LM2.NE.LM3)CG=2.D0*CG
            AUX(:)=AUX(:)+CG*(RHO4(:,LM2,2)*RHO4(:,LM3,2) &
     &                       +RHO4(:,LM2,3)*RHO4(:,LM3,3) &
     &                       +RHO4(:,LM2,4)*RHO4(:,LM3,4)) 
          ENDDO
        ENDDO
        RHO2(:,LM1,2)=AUX(:)
      ENDDO
!
!     == ABSOLUTE VALUE OF THE SPIN DENSITY=====================================
      ARRAY(:)=1.D0/SQRT(SMALL+RHO2(:,1,2)/Y0)
      AUX(:)=ARRAY(:)/(2.D0*Y0)
      DO LM=1,LMRX
        RHO2(:,LM,2)=RHO2(:,LM,2)*AUX(:)
      ENDDO
      RHO2(:,1,2)=2.D0*RHO2(:,1,2)
!!$CALL AUGMENTATION_WRITEPHI('RHO4_X',GID,NR,LMRX,RHO4(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Y',GID,NR,LMRX,RHO4(:,:,3))
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z',GID,NR,LMRX,RHO4(:,:,4))
!!$CALL AUGMENTATION_WRITEPHI('RHO2',GID,NR,LMRX,RHO2(:,:,2))
!!$AUX(:)=SQRT(RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)
!!$CALL AUGMENTATION_WRITEPHI('AUX',GID,NR,1,AUX)
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
!     == TOTAL POTENTIAL =======================================================
      POT4(:,:,1)=POT2(:,:,1)
!
      DO ISIG=2,4
        DO LM=1,LMRX
          A(:,LM,ISIG-1)=RHO4(:,LM,ISIG)*ARRAY(:)
        ENDDO
      ENDDO
!
      B(:,:,:)=0.D0
      DO LM1=1,LMRX
        DO LM2=1,LMRX
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            IF(LM2.EQ.1)CG=2.D0*CG
            DO ISIG=2,4
              B(:,LM1,ISIG-1)=B(:,LM1,ISIG-1)+CG*POT2(:,LM2,2)*A(:,LM3,ISIG-1)        
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      AUX(:)=0.D0
      DO ISIG=2,4
        DO LM=1,LMRX
          AUX(:)=AUX(:)+B(:,LM,ISIG-1)*A(:,LM,ISIG-1)
        ENDDO
      ENDDO
!
      DO ISIG=2,4
        DO LM=1,LMRX
          POT4(:,LM,ISIG)=-0.5D0*AUX(:)*A(:,LM,ISIG-1)+B(:,LM,ISIG-1)
        ENDDO
      ENDDO
      POT4(:,:,2:)=POT4(:,:,2:)/Y0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_PSHARTREE(GID,NR,LMRX,PSRHOC,PSRHO &
     &            ,VADD,RCSM,QLM,VQLM,RHOB,PSPOT,PSEH)
!     **************************************************************************
!     **                                                                      **
!     ** CALCULATES HARTREE ENERGY                                            **
!     ** ADDS COMPENSATION DENSITY TO THE PSEUDO DENSITY                      **
!     **                                                                      **
!     **  1) AEE=AERHO*AEZ/R                                                  **
!     **     PSE=PSRHO*VADD                                                   **
!     **  2) PSRHO=PSRHO+GAUSS(R)*QLM                                         **
!     **  3) AEE=AEE + 0.5*AERHO(R)*AERHO(R')/|R-R'|                          **
!     **     PSE=PSE + 0.5*PSRHO(R)*PSRHO(R')/|R-R'|                          **
!     **     VQLM=VQLM - GAUSS(R)*PSRHO(R')/|R-R'|                            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: PSRHOC(NR)
      REAL(8)    ,INTENT(IN) :: PSRHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: VADD(NR)
      REAL(8)    ,INTENT(IN) :: RCSM
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)    ,INTENT(IN) :: RHOB
      REAL(8)    ,INTENT(OUT):: PSPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: PSEH
      REAL(8)                :: R(NR)
      REAL(8)                :: RHO1(NR)
      REAL(8)                :: POT(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: G(NR)
      REAL(8)                :: RHOHAT(NR,LMRX)
      REAL(8)                :: PI
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      INTEGER(4)             :: LM,L,M,LX
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  INITIALIZE ENERGIES AND POTENTIALS TO  ZERO                 ==
!     ==================================================================
      PSE(:)=0.D0
      PSPOT(:,:)=0.D0
!
!     ==================================================================
!     ==  CONSTRUCT COMPENSATION DENSITY                              ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      LX=INT(SQRT(REAL(LMRX))-1.D0)   ! LMX=(LX+1)**2
      IF((LX+1)**2.NE.LMRX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMRX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMRX',LMRX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('AUGMENTATION_PSHARTREE')
      END IF
      LM=0
      DO L=0,LX
        CALL GAUSSN(L,ALPHA,CL)
        G(:)=CL*EXP(-ALPHA*R(:)**2)*R(:)**L
        DO M=1,2*L+1
          LM=LM+1
          RHOHAT(:,LM)=QLM(LM)*G(:)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==   ADD UNSCREENING POTENTIAL VADD                             ==
!     ==================================================================
      PSE(:)    =VADD(:)*(PSRHO(:,1)+PSRHOC(:))
      PSPOT(:,1)=PSPOT(:,1)+VADD(:)
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL GAUSSN(L,ALPHA,CL)
        SVAR=QLM(LM)*CL
        RHO1(:)=PSRHO(:,LM)+RHOHAT(:,LM)
        IF(L.EQ.0)RHO1(:)=RHO1(:)+PSRHOC(:)
        CALL RADIAL$POISSON(GID,NR,L,RHO1,POT)
        PSE(:)=PSE(:)+0.5D0*POT(:)*RHO1(:)
        PSPOT(:,LM)=PSPOT(:,LM)+POT(:)
      ENDDO
!
!     ==================================================================
!     ==  ADD EXTERNAL POTENTIAL                                      ==
!     ==================================================================
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        POT(:)=VQLM(LM)*R(:)**L
        PSPOT(:,LM)=PSPOT(:,LM)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
        RHO1(:)=PSRHO(:,LM)+RHOHAT(:,LM)
        IF(L.EQ.0) RHO1(:)=RHO1(:)+PSRHOC(:)
        PSE(:)=PSE(:)+RHO1(:)*POT(:)
      ENDDO
!
!     ==================================================================
!     ==  ADD POTENTIAL FROM THE BACKGROUND                           ==
!     ==================================================================
      SVAR=-2.D0*PI*RHOB/3.D0*SQRT(4.D0*PI)
      POT(:)=SVAR*R(:)**2
      PSPOT(:,1)=PSPOT(:,1)+POT(:) ! POTENTIAL OF THE BACKGROUND
      PSE(:)=PSE(:)+(PSRHO(:,1)+RHOHAT(:,1)+PSRHOC(:))*POT(:)
!
!     ==================================================================
!     ==  CALCULATE TOTAL ENERGY                                      ==
!     ==================================================================
      PSE(:)=PSE(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,PSE,PSEH)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_AEHARTREE(GID,NR,LMRX,AEZ,RHOC,AERHO &
     &                                 ,VQLM,RHOB,AEPOT,AEEH)
!     ******************************************************************
!     **  ELECTROSTATIC ENERGY OF THE ALL-ELECTRON ONE-CENTER DENSITY **
!     **  INCLUDING THE EXTERNAL POTENTIAL AND THE POTENTIAL OF THE   **
!     **  COMPENSATING CHARGE BACKGROUND FOR CHARGED SYSTEMS          **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: AEZ        ! ATOMIC NUMBER
      REAL(8)    ,INTENT(IN) :: RHOC(NR)   ! CORE DENSITY
      REAL(8)    ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)    ,INTENT(IN) :: RHOB       ! COMPENSATING BACKGROUND
      REAL(8)    ,INTENT(OUT):: AEPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: AEEH       ! ENERGY
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX(NR)
      REAL(8)                :: POT(NR)
      REAL(8)                :: AEE(NR)  ! ENERGY DENSITY
      REAL(8)                :: SVAR
      REAL(8)                :: EEXT
      REAL(8)                :: PI,Y0
      INTEGER(4)             :: LM
      INTEGER(4)             :: L
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  INITIALIZE ENERGIES AND POTENTIALS TO  ZERO                 ==
!     ==================================================================
      AEE(:)=0.D0
      AEPOT(:,:)=0.D0
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC VALENCE-CORE INTERACTION           ==
!     ==   AND VALENCE-NUCLEUS INTERACTION                            ==
!     ==================================================================
      CALL RADIAL$POISSON(GID,NR,0,RHOC,POT)
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX)
      POT(:)  = POT(:)+AUX(:)
      AEE(:)    = POT(:)* AERHO(:,1)
      AEPOT(:,1)= AEPOT(:,1)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
!      AEE(:)=AE(:)+RHOC(:)*VQLM(1)
!AUX(:)=AEE(:)*R(:)**2
!CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!PRINT*,'EL C-V INT', SVAR
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
!AUX(:)=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL RADIAL$POISSON(GID,NR,L,AERHO(:,LM),POT)
        AEE(:)=AEE(:)+0.5D0*POT(:)*AERHO(:,LM)
        AEPOT(:,LM)=AEPOT(:,LM)+POT(:)
!AUX(:)=AUX(:)+0.5D0*POT(:)*AERHO(:,LM)*R(:)**2
      ENDDO
!CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!PRINT*,'EL V-V INT',SVAR
!
!AUX(:)=AEE(:)*R(:)**2
!CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!PRINT*,'EL HARTREE INT ',SVAR
!
!     ==================================================================
!     ==  ADD EXTERNAL POTENTIAL                                      ==
!     ==================================================================
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        POT(:)=VQLM(LM)*R(:)**L
        AEPOT(:,LM)=AEPOT(:,LM)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
        AEE(:)=AEE(:)+AERHO(:,LM)*POT(:)
        IF(L.EQ.0)AEE(:)=AEE(:)+RHOC(:)*POT(:)
      ENDDO
      EEXT=-VQLM(1)*Y0*AEZ  ! WILL BE ADDED TO THE TOTAL ENERGY
!
!     ==================================================================
!     ==  ADD POTENTIAL FROM THE BACKGROUND                           ==
!     ==================================================================
      SVAR=-2.D0*PI*RHOB/3.D0*SQRT(4.D0*PI)
      POT(:)=SVAR*R(:)**2
      AEPOT(:,1)=AEPOT(:,1)+POT(:) ! POTENTIAL OF THE BACKGROUND
      AEE(:)=AEE(:)+(AERHO(:,1)+RHOC(:))*POT(:)
!CALL RADIAL$INTEGRAL(GID,NR,(AERHO(:,1)+RHOC(:))*POT(:)*R**2,SVAR)
!PRINT*,'EL BACKGROUND',SVAR
!
!     ==================================================================
!     ==  CALCULATE ENERGY                                            ==
!     ==================================================================
      AEE(:)=AEE(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AEE,AEEH)
      AEEH=AEEH+EEXT
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADJUSTVQLM(GID,NR,LMRX,PSRHO,RCSM,QLM,VQLM)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      REAL(8)    ,INTENT(IN) :: RCSM
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: PSRHO(NR,LMRX) ! INCLUDES PSEUDO CORE DENSITY
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(OUT):: VQLM(LMRX)
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX(NR)
      REAL(8)                :: RHO(NR)
      REAL(8)                :: G(NR)
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      INTEGER(4)             :: L,M,LM
      INTEGER(4)             :: LX
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      LX=INT(SQRT(REAL(LMRX))-1.D0) ! LMX=(LX+1)**2
      IF((LX+1)**2.NE.LMRX) THEN   
        CALL ERROR$MSG('LMRX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMRX',LMRX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('AUGMENTATION_ADJUSTVQLM')
      END IF
      LM=0
      DO L=0,LX
        CALL GAUSSN(L,ALPHA,CL)
        G(:)=CL*EXP(-ALPHA*R(:)**2)*R(:)**L
        DO M=1,2*L+1
          LM=LM+1
          RHO=PSRHO(:,LM)+QLM(LM)*G(:)
          CALL RADIAL$POISSON(GID,NR,L,RHO,AUX)
          AUX(:)=AUX(:)*G(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          VQLM(LM)=-SVAR
        ENDDO
      ENDDO
      RETURN
      END
!     ................................................................
      SUBROUTINE AUGMENTATION_FEEDHYPERFINE(GID,NR,LMRX,NDIMD &
     &         ,IAT,AERHO,AECORE,PSRHO,PSCORE,AEHPOT,PSHPOT)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: PSRHO(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(IN) :: AEHPOT(NR,LMRX)
      REAL(8)   ,INTENT(IN) :: PSHPOT(NR,LMRX)
      REAL(8)               :: RHO(NR,LMRX)
!     ****************************************************************
!
!     ================================================================
!     ==  TOTAL DENSITY                                             ==
!     ================================================================
      CALL HYPERFINE$SET1CRHO('PS','TOT',IAT,GID,NR,NR,LMRX,PSRHO)
      RHO(:,:)=AERHO(:,:,1)
      RHO(:,1)=RHO(:,1)+AECORE(:)
      CALL HYPERFINE$SET1CRHO('AE','TOT',IAT,GID,NR,NR,LMRX,RHO)
!
!     ================================================================
!     ==  SPIN DENSITY                                              ==
!     ================================================================
      IF(NDIMD.GT.1) THEN
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          RHO(:,:)=PSRHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          RHO(:,:)=SQRT(PSRHO(:,:,2)**2+PSRHO(:,:,3)**2+PSRHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('PS','SPIN',IAT,GID,NR,NR,LMRX,RHO)
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          RHO(:,:)=AERHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          RHO(:,:)=SQRT(AERHO(:,:,2)**2+AERHO(:,:,3)**2+AERHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('AE','SPIN',IAT,GID,NR,NR,LMRX,RHO)
      END IF
!
!     ================================================================
!     ==  TOTAL POTENTIAL FOR ELECTRIC FIELD GRADIENTS              ==
!     ================================================================
      CALL HYPERFINE$SET1CPOT('AE',IAT,GID,NR,NR,LMRX,AEHPOT)
      CALL HYPERFINE$SET1CPOT('PS',IAT,GID,NR,NR,LMRX,PSHPOT)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &           ,AEPOT,PSPOT,AEPHI,PSPHI,DATH)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE EXPECTATION VALUE OF                         **
!     **  THE ONE-CENTER POTENTIALS WITH THE PARTIAL WAVES            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: PSPOT(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: DATH(LMNX,LMNX,NDIMD)
      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LM1,LM2,LM3
      INTEGER(4)            :: L1,L2
      INTEGER(4)            :: IM1,IM2
      INTEGER(4)            :: ISPIN
      REAL(8)               :: AEDMU(NR,NDIMD)
      REAL(8)               :: PSDMU(NR,NDIMD)
      REAL(8)               :: DWORK1(NR)
      REAL(8)               :: CG
      REAL(8)               :: SVAR
      REAL(8)               :: R(NR)
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      DATH(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          LM1=L1**2+IM1
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
!     
!             ==========================================================
!             ==  SUM ALL POTENTIALS THAT ACT ON THE GIVEN PAIR       ==
!             ==  OF PARTIAL WAVES                                    ==
!             ==========================================================
              AEDMU(:,:)=0.D0
              PSDMU(:,:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  DO ISPIN=1,NDIMD
                    AEDMU(:,ISPIN)=AEDMU(:,ISPIN)+CG*AEPOT(:,LM3,ISPIN)
                    PSDMU(:,ISPIN)=PSDMU(:,ISPIN)+CG*PSPOT(:,LM3,ISPIN)
                  ENDDO
                END IF
              ENDDO
!     
!             ==========================================================
!             ==  PERFORM NOW THE INTEGRATION                         ==
!             ==========================================================
              DO ISPIN=1,NDIMD
                DWORK1(:)= &
     &               (AEDMU(:,ISPIN)*AEPHI(:,LN1)*AEPHI(:,LN2) &
     &               -PSDMU(:,ISPIN)*PSPHI(:,LN1)*PSPHI(:,LN2))*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,SVAR)
                DATH(LMN1,LMN2,ISPIN)=DATH(LMN1,LMN2,ISPIN)+SVAR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE AUGMENTATION_CORELEVELS(NDIMD,R1,DEX,NR,LMRX,AEPOT)
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: LMRX
!!$      INTEGER(4),INTENT(IN) :: NDIMD
!!$      REAL(8)   ,INTENT(IN) :: R1
!!$      REAL(8)   ,INTENT(IN) :: DEX
!!$      INTEGER(4),INTENT(IN) :: NR
!!$      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
!!$      INTEGER(4)            :: NC         ! NUMBER OF CORE STATES
!!$      INTEGER(4),ALLOCATABLE:: LOFC(:)    !(NC)
!!$      REAL(8)   ,ALLOCATABLE:: EOFC(:)    !(NC)
!!$      REAL(8)   ,ALLOCATABLE:: VATOM(NR)  !(NC)
!!$      REAL(8)   ,ALLOCATABLE:: PHIC(:,:)  !(NR,NC)
!!$      REAL(8)   ,ALLOCATABLE:: R2(:)   ! R**2
!!$      REAL(8)   ,ALLOCATABLE:: H(:,:)  ! HAMILTONIAN
!!$      COMPLEX(8),ALLOCATABLE:: UC(:,:) ! CORE EIGENVECTORS
!!$      REAL(8)   ,ALLOCATABLE:: EIGC(:) ! CORE EIGENVALUES
!!$      INTEGER(4)            :: NE ! DIMENSIONA OF THE HAMILTONIAN
!!$      INTEGER(4)            :: I1,I2,M1,M2,L1,L2,LM1,LM2,LM3
!!$!     ******************************************************************
!!$      IF(NDIMD.NE.1) THEN
!!$        CALL ERROR$STOP('AUGMENTATION_CORELEVELS')
!!$      END IF
!!$      PI=4.D0*ATAN(1.D0)
!!$      C000=1/SQRT(4.D0*PI)
!!$!
!!$!     ==================================================================
!!$!     ==  READ CORE LEVELS                                            ==
!!$!     ==================================================================
!!$!===  READ NC 
!!$      ALLOCATE(LOFC(NC))
!!$      ALLOCATE(PHIC(NR,NC))
!!$!
!!$!     ==================================================================
!!$!     ==                                                              ==
!!$!     ==================================================================
!!$      NE=0      
!!$      DO I=1,NC
!!$        NE=NE+2*LC(I)+1
!!$      ENDDO
!!$      ALLOCATE(H(NC,NC))
!!$      ALLOCATE(EIGC(NC))
!!$      ALLOCATE(UCC(NC,NC))
!!$!
!!$!     ==================================================================
!!$!     ==  RADIAL GRID R2=R**2                                         ==
!!$!     ==================================================================
!!$      XEXP=EXP(DEX)
!!$      RI=R1/XEXP
!!$      DO IR=1,NR
!!$        RI=RI*XEXP
!!$        R2(IR)=RI**2
!!$      ENDDO
!!$!
!!$!     ==================================================================
!!$!     ==  CALCULATE HAMILTONIAN                                       ==
!!$!     ==================================================================
!!$      DO I1=1,NC
!!$        L1=LOFC(I1)
!!$        DO I2=I1,NC
!!$          L2=LOFC(I2)
!!$          AUX(:)=PHI(:,I1)*PHI(:,I2)*R2(:)
!!$          DO M1=1,2*L1+1
!!$            LM1=L1**2+M1
!!$            DO M2=1,2*L2+1
!!$              LM2=L2**2+M2
!!$              SUM(:)=0.D0
!!$              DO LM3=1,LMRX
!!$                CALL CLEBSCH(LM1,LM2,LM3,CG)
!!$                IF(CG.EQ.0.D0) CYCLE
!!$                SUM(:)=SUM(:)+CG*AEPOT(:,LM3,1)
!!$              ENDDO
!!$              IF(LM1.EQ.LM2) SUM(:)=SUM(:)-VC(:)*C000
!!$              SUM(:)=SUM(:)*AUX(:)
!!$              CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,SUM,H(IE1,IE2))
!!$              H(IE2,IE1)=H(IE1,IE2)
!!$! ADD EIGENVALUES
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==================================================================
!!$!     ==  DIAGONALIZE HAMILTONIAN                                     ==
!!$!     ==================================================================
!!$      CALL LIB$DIAGR8(NE,H,EIGC,UC)
!!$!
!!$!     ==================================================================
!!$!     ==  PRINT EIGENVALUES                                           ==
!!$!     ==================================================================
!!$      CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$      WRITE(NFILO,FMT='("CORE EIGENVALUES"//"================")')
!!$      WRITE(NFILO,FMT='(10F10.3)')EIGC
!!$      RETURN
!!$      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE EXPERTNAL1CPOT_MODULE
!*******************************************************************************
!**                                                                           **
!**  APPLIES AN EXTERNAL POTENTIAL ACTING ON THE ELECTRONS                    **
!**  WITHIN THE AUGMENTATION REGION                                           **
!**                                                                           **
!**    IDIMD EXTERNALLY THIS INDEX REFERS TO                                  **
!**             0: TOTAL  FOR NDIMD=1,2,4                                     **
!**             1: SZ FOR NDIMD=2                                             **
!**             1: SX FOR NDIMD=4                                             **
!**             2: SY FPR NDIMD=4                                             **
!**             3: SZ FPR NDIMD=4                                             **
!**          INTERNALLY THE VALUE IS RAISED BY ONE SO THAT IT                 **
!**          CORRESPONDS TO THE INDEXING                                      **
!**                                                                           **
!*******************************************************************************
TYPE EXTPOT
  REAL(8)          :: VALUE
  CHARACTER(LEN=32):: ATOM
  REAL(8)          :: RC  
  REAL(8)          :: PWR
  CHARACTER(LEN=32):: TYPE   ! CAN BE 'S','P','D','F','ALL'
  INTEGER(4)       :: IDIMD  ! IDIMD=0 REFERES TO ALL SPIN DIRECTIONS
END TYPE EXTPOT
INTEGER(4)             :: NPOT=0
INTEGER(4)             :: NPOTX=0
INTEGER(4),PARAMETER   :: DNPOT=5
TYPE (EXTPOT), ALLOCATABLE :: POT(:)
!*******************************************************************************
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CREATE
      TYPE (EXTPOT),ALLOCATABLE :: TMPPOT(:)
!     **************************************************************************
      IF(NPOT.LT.NPOTX) RETURN
      IF(ALLOCATED(POT)) THEN
        ALLOCATE(TMPPOT(NPOTX))
        TMPPOT=POT
        DEALLOCATE(POT)
        NPOTX=NPOTX+DNPOT
        ALLOCATE(POT(NPOTX))
        POT(1:NPOT)=TMPPOT(1:NPOT)
        DEALLOCATE(TMPPOT)
      ELSE
        NPOTX=DNPOT
        ALLOCATE(POT(NPOTX))
      END IF          
      END SUBROUTINE CREATE
END MODULE EXPERTNAL1CPOT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTERNAL1CPOT$SETPOT(ATOM,TYPE,IDIMD,VALUE,RC,PWR)
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      REAL(8)          ,INTENT(IN) :: VALUE
      CHARACTER(LEN=32),INTENT(IN) :: ATOM
      CHARACTER(LEN=32),INTENT(IN) :: TYPE
      INTEGER(4)       ,INTENT(IN) :: IDIMD !SPIN DIRECTION OR 0
      REAL(8)          ,INTENT(IN) :: RC
      REAL(8)          ,INTENT(IN) :: PWR
!     **************************************************************************
      CALL CREATE
      NPOT=NPOT+1
      POT(NPOT)%ATOM =ATOM
      POT(NPOT)%VALUE=VALUE
      POT(NPOT)%RC=RC
      POT(NPOT)%PWR=PWR
      POT(NPOT)%IDIMD=IDIMD+1
      POT(NPOT)%TYPE =TYPE
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTERNAL1CPOT$REPORT(NFIL)
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IPOT
!     **************************************************************************
      IF(NPOT.EQ.0) RETURN
!     CALL REPORT$TITLE(NFIL,"EXTERNAL POTENTIALS ON ORBITALS")
      WRITE(NFIL,*)
      WRITE(NFIL,FMT='("EXTERNAL POTENTIALS ON ORBITALS"/30("="))')
      DO IPOT=1,NPOT
!       CALL REPORT$I4VAL(NFIL,"POTENTIAL NR:",IPOT,' ')
!       CALL REPORT$CHVAL(NFIL,'ATOM',POT(IPOT)%ATOM)
!       CALL REPORT$CHVAL(NFIL,'TYPE',POT(IPOT)%TYPE)
!       CALL REPORT$R8VAL(NFIL,'VALUE',POT(IPOT)%VALUE,'H')
!       CALL REPORT$I4VAL(NFIL,'SPIN',POT(IPOT)%NSPIN,' ')
        WRITE(NFIL,FMT='("POTENTIAL NR: ",I5)')IPOT
        CALL WRITECH(NFIL,'ATOM',POT(IPOT)%ATOM)
        CALL WRITECH(NFIL,'TYPE',POT(IPOT)%TYPE)
        CALL WRITER8(NFIL,'VALUE',POT(IPOT)%VALUE,'H')
        CALL WRITER8(NFIL,'PWR',POT(IPOT)%PWR,' ')
        CALL WRITER8(NFIL,'RC',POT(IPOT)%RC,'ABOHR')
        CALL WRITEI4(NFIL,'IDIMD ([0=NT],[0=NT,1=NS],[0=NT,1=NX,2=NY,3=NZ])' &
     &                   ,POT(IPOT)%IDIMD-1,' ')
      ENDDO
      RETURN
      CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE WRITECH(NFIL,NAME,VALUE)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        CHARACTER(LEN=*),INTENT(IN) :: VALUE
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,A)')NAME,TRIM(VALUE)
        RETURN
        END SUBROUTINE WRITECH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE WRITER8(NFIL,NAME,VALUE,UNIT)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        REAL(8)          ,INTENT(IN) :: VALUE
        CHARACTER(LEN=*),INTENT(IN) :: UNIT
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,F10.5," ",A)')NAME,VALUE,UNIT
        RETURN
        END SUBROUTINE WRITER8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE WRITEI4(NFIL,NAME,VALUE,UNIT)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        INTEGER(4)       ,INTENT(IN) :: VALUE
        CHARACTER(LEN=*),INTENT(IN) :: UNIT
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,I10," ",A)')NAME,VALUE,UNIT
        RETURN
        END SUBROUTINE WRITEI4
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATH,ETOT)
!     **************************************************************************
!     **  TOTAL ENERGY AND ONE-CENTER HAMILTONIAN CONTRIBUTIONS FROM          **
!     **  EXTERNAL POTENTIAL                                                  **
!     **************************************************************************
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: ATOM
      INTEGER(4)  ,INTENT(IN)   :: NDIMD
      INTEGER(4)  ,INTENT(IN)   :: LMNX
      COMPLEX(8)  ,INTENT(IN)   :: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(OUT)  :: ETOT
      TYPE(EXTPOT)              :: POT1
      INTEGER(4)                :: LANG
      INTEGER(4)                :: MANG
      INTEGER(4)                :: IAT    ! ATOM INDEX
      INTEGER(4)                :: ISP    ! ATOM TYPE INDEX
      INTEGER(4)                :: NR     ! #(GRID RADIAL POINTS)
      INTEGER(4)                :: LNX    ! #(PARTIAL WAVES (L,N))
      INTEGER(4)   ,ALLOCATABLE :: LOX(:) !(LNX) #(GRID RADIAL POINTS)
      REAL(8)      ,ALLOCATABLE :: AEPHI(:,:) !(NR,LNX) AE PARTIAL WAVES
      REAL(8)      ,ALLOCATABLE :: UONE(:,:)  !(LNX,LNX)
      REAL(8)      ,ALLOCATABLE :: RDEP(:)    !(NR)
      REAL(8)      ,ALLOCATABLE :: AUX(:)    !(NR)
      CHARACTER(32)             :: SPECIES
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: LN1,LN2,LMN1,LMN2,IPOT,IDIMD,L1,L2,I
      INTEGER(4)                :: GID     ! GRID ID
      REAL(8)      ,ALLOCATABLE :: R(:)    ! RADIAL GRID
      LOGICAL(4)                :: TACTIVE
      INTEGER(4)                :: NFILTRACE
!     **************************************************************************
      DATH(:,:,:)=0.D0
      ETOT=0.D0
!
!     ==========================================================================
!     ==  SELECT POTENTIAL                                                    ==
!     ==========================================================================
!     == DETERMINE ID OF THE ATOM TYPE (SPECIES)
      CALL ATOMLIST$INDEX(ATOM,IAT)
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
!
      TACTIVE=.FALSE.
      DO IPOT=1,NPOT
        POT1=POT(IPOT)
        TCHK=(ATOM.EQ.POT1%ATOM.OR.SPECIES.EQ.POT1%ATOM)
        IF(.NOT.TCHK) CYCLE
        TACTIVE=.TRUE.
!
!       ========================================================================
!       ==  UNSCRAMBLE TYPE                                                   ==  
!       ========================================================================
!       == LANG=MAIN ANGULAR MOMENTUM  LANG=-1 ALL ANGULAR MOMENTA
!       == MANG="MAGNETIC" QUANTUM NUMBER  MANG=-1 ALL "MAGNETIC" QUANTUM NUMBERS
!       == CAUTION: MANG REFERS TO REAL SPHERICAL HARMONICS
        MANG=-1
        IF(TRIM(POT1%TYPE).EQ.'S') THEN
          LANG=0
!
        ELSE IF(TRIM(POT1%TYPE(1:1)).EQ.'P') THEN
          LANG=1
          IF(TRIM(POT1%TYPE(2:)).EQ.'X') THEN
            MANG=1  
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'Z') THEN
            MANG=2  
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'Y') THEN
            MANG=3
          ELSE IF(LEN_TRIM(POT1%TYPE(2:)).GT.0) THEN
            CALL ERROR$MSG('VALUE OF TYPE NOT RECOGNIZED')
            CALL ERROR$CHVAL('POT1%TYPE',POT1%TYPE)
            CALL ERROR$STOP('EXTERNAL1CPOT$APPLY')
          END IF
!
        ELSE IF(TRIM(POT1%TYPE(1:1)).EQ.'D') THEN
          LANG=2
          IF(TRIM(POT1%TYPE(2:)).EQ.'X2-Y2') THEN
            MANG=1  
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'XZ') THEN
            MANG=2  
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'3Z2-R2') THEN
            MANG=3
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'YZ') THEN
            MANG=4
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'XY') THEN
            MANG=5
          ELSE IF(LEN_TRIM(POT1%TYPE(2:)).GT.0) THEN
            CALL ERROR$MSG('VALUE OF TYPE NOT RECOGNIZED')
            CALL ERROR$CHVAL('POT1%TYPE',POT1%TYPE)
            CALL ERROR$STOP('EXTERNAL1CPOT$APPLY')
          END IF
!
        ELSE IF(TRIM(POT1%TYPE(1:1)).EQ.'F') THEN
          LANG=3
          IF(TRIM(POT1%TYPE(2:)).EQ.'1') THEN
            MANG=1  
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'2') THEN
            MANG=2  
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'3') THEN
            MANG=3
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'4') THEN
            MANG=4
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'5') THEN
            MANG=5
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'6') THEN
            MANG=6
          ELSE IF(TRIM(POT1%TYPE(2:)).EQ.'7') THEN
            MANG=7
          ELSE IF(LEN_TRIM(POT1%TYPE(2:)).GT.0) THEN
            CALL ERROR$MSG('VALUE OF TYPE NOT RECOGNIZED')
            CALL ERROR$CHVAL('POT1%TYPE',POT1%TYPE)
            CALL ERROR$STOP('EXTERNAL1CPOT$APPLY')
          END IF
!
        ELSE IF(TRIM(POT1%TYPE).EQ.'ALL') THEN
          LANG=-1
!
        ELSE
          CALL ERROR$MSG('VALUE OF TYPE NOT RECOGNIZED')
          CALL ERROR$CHVAL('POT1%TYPE',POT1%TYPE)
          CALL ERROR$STOP('EXTERNAL1CPOT$APPLY')
        END IF
!
!       ========================================================================
!       == COLLECT PARTIAL WAVES                                              ==
!       ========================================================================
        CALL ATOMLIST$INDEX(ATOM,IAT)
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETCH('ID',SPECIES)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        ALLOCATE(AEPHI(NR,LNX))
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
        CALL SETUP$UNSELECT()
!
!       ========================================================================
!       == SPECIFY SHAPE OF THE EXTERNAL POTENTIAL                            ==
!       ========================================================================
        ALLOCATE(RDEP(NR))
        ALLOCATE(AUX(NR))
        RDEP(:)=POT1%VALUE*EXP(-(R(:)/POT1%RC)**POT1%PWR)
!
!       ========================================================================
!       == DETERMINE SCALAR PRODUCTS WITH PARTIAL WAVES                       ==
!       == UONE(LN1,LN1)=<AEPHI(LN1)|U(R)|AEPHI(LN2)>                         ==
!       ========================================================================
        ALLOCATE(UONE(LNX,LNX))
        UONE(:,:)=0.D0
        LMN1=0
        DO LN1=1,LNX
          IF(.NOT.(LANG.EQ.LOX(LN1).OR.LANG.EQ.-1)) CYCLE
          DO LN2=1,LNX
            IF(.NOT.(LANG.EQ.LOX(LN2).OR.LANG.EQ.-1)) CYCLE
            IF(LOX(LN2).NE.LOX(LN1)) CYCLE
            IF(LN2.LT.LN1) THEN
              UONE(LN1,LN2)=UONE(LN2,LN1)
              CYCLE
            END IF
            AUX(:)=RDEP(:)*AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,UONE(LN1,LN2))
          ENDDO
        ENDDO
        DEALLOCATE(R)
        DEALLOCATE(RDEP)
        DEALLOCATE(AUX)
        DEALLOCATE(AEPHI)
!
!       ========================================================================
!       ==  CALCULATE ONE-CENTER HAMILTONIAN                                  ==
!       ========================================================================
        IF(POT1%IDIMD.LT.1.OR.POT1%IDIMD.GT.NDIMD) THEN
          CALL ERROR$MSG('VALUE OF POT1%IDIMD NOT RECOGNIZED')
          CALL ERROR$I4VAL('POT1%IDIMD',POT1%IDIMD)
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$STOP('EXTERNAL1CPOT$APPLY')
        END IF
        DO IDIMD=1,NDIMD
          IF(IDIMD.EQ.POT1%IDIMD) THEN
            LMN1=0
            DO LN1=1,LNX
              L1=LOX(LN1)
              LMN2=0
              DO LN2=1,LNX
                L2=LOX(LN2)
                IF(L1.EQ.L2) THEN
                  DO I=1,2*L1+1
                    IF(.NOT.(MANG.EQ.I.OR.MANG.EQ.-1)) CYCLE
                    DATH(LMN1+I,LMN2+I,IDIMD)=DATH(LMN1+I,LMN2+I,IDIMD) &
       &                                     +UONE(LN1,LN2)
                  ENDDO
                END IF
                LMN2=LMN2+2*L2+1
              ENDDO
              LMN1=LMN1+2*L1+1
            ENDDO
          END IF
        ENDDO
        DEALLOCATE(UONE)
        DEALLOCATE(LOX)
      ENDDO
!
!     ==========================================================================
!     ==  RETURN IF NO POTENTIAL FOR THIS ATOM FOUND                          ==
!     ==========================================================================
      IF(.NOT.TACTIVE) RETURN

CALL TRACE$GETL4('TRACEFILE',TCHK)
IF(TCHK) THEN
  CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
  WRITE(NFILTRACE,*)'DENSITY MATRIX IN EXTERNAL1CPOT. ATOM ',POT1%ATOM 
  WRITE(NFILTRACE,*)TRIM(ATOM),TRIM(SPECIES),LNX,LMNX
  DO IDIMD=1,NDIMD
    WRITE(NFILTRACE,*)'IDIMD',IDIMD
    DO LMN1=1,LMNX
      WRITE(NFILTRACE,FMT='(20F10.3)')REAL(DENMAT(LMN1,:,IDIMD))  
    ENDDO
  ENDDO
  WRITE(NFILTRACE,*)'DELTA HAMILTONIAN FROM EXTERNAL1CPOT'
  DO IDIMD=1,NDIMD
    WRITE(NFILTRACE,*)'IDIMD',IDIMD
    DO LMN1=1,LMNX
      WRITE(NFILTRACE,FMT='(20F10.3)')DATH(LMN1,:,IDIMD) 
    ENDDO
  ENDDO
END IF
!
!     ==========================================================================
!     ==  CALCULATE TOTAL ENERGY                                              ==
!     ==========================================================================
      ETOT=0.D0
      DO IDIMD=1,NDIMD
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            ETOT=ETOT+REAL(DENMAT(LMN1,LMN2,IDIMD)*DATH(LMN2,LMN1,IDIMD),KIND=8)
          ENDDO
        ENDDO
      ENDDO
!     PRINT*,'ETOT ',ETOT
      RETURN 
      END
!
!     ..................................................................
      SUBROUTINE GAUSSPOT(RC,L,R,V)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE HARTREE POTENTIAL OF A GENERALIZED           **
!     **  GAUSSIAN  R**L*EXP(-(R/RC)**2)*YLM                          **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8),   INTENT(IN) :: RC
      REAL(8),   INTENT(IN) :: R
      INTEGER(4),INTENT(IN) :: L
      REAL(8),   INTENT(OUT):: V
      REAL(8)               :: PI
      REAL(8)               :: FAC,X,SVAR1
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      FAC=4.D0*PI/DBLE(2*L+1)*RC**DBLE(L+2)
      X=R/RC
      CALL DG2NDX(L+1,X,SVAR1)
      V=FAC*(SVAR1*X**DBLE(-L-1)+0.5D0*X**DBLE(L)*EXP(-X**2))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE DG2NDX(N,X,V)
!     ******************************************************************
!     **                                                              **
!     **  V = INT{DR(0,R)| X**(2*N) *  EXP(-(X)**2) }                 **
!     **                                                              **
!     **  START FROM D/DX{X**(2*N-1)*EXP(-X**2)}                      **
!     **  BUILD RECURSION FOR THE INTEGRAL                            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8),   INTENT(IN) :: X
      REAL(8),   INTENT(OUT):: V
      REAL(8)               :: PI
      REAL(8)               :: X2,GAUSS,FAC1,FAC2
      REAL(8)               :: ERFX
      INTEGER(4)            :: I
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      X2=X*X
      CALL LIB$ERFR8(X,ERFX)
      V=0.5D0*SQRT(PI)*ERFX
      GAUSS=EXP(-X2)
      FAC1=0.5D0
      FAC2=0.5D0*X
      DO I=1,N
        V=FAC1*V - FAC2*GAUSS 
        FAC1=FAC1+1.D0
        FAC2=FAC2*X2
      ENDDO
      RETURN
      END
!!$!
!!$!     .....................................................................
!!$      SUBROUTINE AUGMENTATION_SOFTCORE(ISP,R1,DEX,NR,AERHO1,AEPOT1)
!!$      USE SCHRGL_INTERFACE_MODULE, ONLY : OUTBOUNDARY
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: ISP
!!$      REAL(8)   ,INTENT(IN) :: R1
!!$      REAL(8)   ,INTENT(IN) :: DEX
!!$      INTEGER(4),INTENT(IN) :: NR
!!$      INTEGER(4),INTENT(IN) :: LMRX
!!$      REAL(8)   ,INTENT(IN) :: AERHO1(NR)
!!$      REAL(8)   ,INTENT(IN) :: AEPOT1(NR)
!!$      REAL(8)               :: AECORE0(NR)
!!$      REAL(8)               :: PHI(NR,3)
!!$      REAL(8)   ,PARAMETER  :: TOLB=1.D-6
!!$      INTEGER(4),PARAMETER  :: IBIX=100
!!$      INTEGER(4)            :: NB
!!$      INTEGER(4),ALLOCATABLE:: LOFB(:)  ! MAIN ANGULAR MOMENTUM
!!$      INTEGER(4),ALLOCATABLE:: NNOFB(:) ! #(NODES)
!!$      REAL(8)   ,ALLOCATABLE:: FOFB(:)  ! OCCUPATION
!!$      REAL(8)   ,ALLOCATABLE:: EOFB(:)  ! ENERGY
!!$      REAL(8)   ,ALLOCATABLE:: RCLOFB(:) 
!!$      REAL(8)               :: AEZ
!!$      INTEGER(4)            :: NITER
!!$      INTEGER(4)            :: ITER,IB
!!$!     **********************************************************************
!!$      CALL SETUP$AECORE(ISP,NR,AECORE0)
!!$      CALL SETUP$AEZ(ISP,AEZ)
!!$!     ======================================================================
!!$!     ==  HARDWIRE PARAMETERS FOR IRON                                    ==
!!$!     ======================================================================
!!$      ALLOCATE(LOFB(NB))
!!$      ALLOCATE(FOFB(NB))
!!$      ALLOCATE(EOFB(NB))
!!$      ALLOCATE(RCLOFB(NB))
!!$      NB=5
!!$      LOFB=(/0,0,1,0,1/)
!!$      NNOFB=(/0,1,0,2,1/)
!!$      DO IB=1,NB
!!$        FOFB(IB)=REAL(2*L+1)
!!$        ZEFF=AEZ-SUM(FOFB(1:IB-1))
!!$        RCLOFB(IB)=1.D0/ZEFF/REAL(LOFB(IB)+NNOFB(IB)+1)
!!$        IR=1+NINT(LOG(RCLOFB(IB)/R1)/DEX)
!!$        EOFC(IB)=AEPOT(IR)
!!$      ENDDO
!!$
!!$!     ======================================================================
!!$!     ==  START LOOP                                                      ==
!!$!     ======================================================================
!!$      NITER=2
!!$      DO ITER=1,NITER
!!$!
!!$!       ==============================================================
!!$!       ==  FIND CORE STATES                                        ==
!!$!       ==============================================================
!!$        AECORE(:)=0.D0
!!$        DO IB=1,NLB
!!$          X0=EOFB(IB)
!!$          ISTART=1
!!$          TCONV=.FALSE.
!!$          CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
!!$          DO IBI=1,IBIX
!!$            TFORWARD=.TRUE.
!!$            IR=1+NINT(LOG(RCLOFB(IB)/R1)/DEX)
!!$            CALL SCHROEDER(TFORWARD,.FALSE.,R1,DEX,NR,LOFB(IB),EOFB(IB) &
!!$     &                    ,AEZ,IR,4.D0 &
!!$     &                     ,AEPOT,PHI,GINH,DLG)
!!$            DIN=RCLOFB(IB)*PHI(IR,2)/PHI(IR,1)
!!$            DIN=NNOFB(IB)+0.5D0+ATAN(DIN)/PI            
!!$            TFORWARD=.FALSE.
!!$            CALL SCHROEDER(TFORWARD,.FALSE.,R1,DEX,NR,LOFB(IB),EOFB(IB) &
!!$     &                    ,AEZ,IR,4.D0 &
!!$     &                     ,AEPOT,PHI,GINH,DLG)
!!$            DOUT=RCLOFB(IB)*PHI(IR,2)/PHI(IR,1)
!!$            DIN=0.5D0+ATAN(DOUT)/PI            
!!$            Y0=DIN-DOUT
!!$            TCONV=ABS(DX).LT.TOLB
!!$            IF(TCONV) EXIT
!!$            CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
!!$            EOFB(IB)=X0
!!$          ENDDO
!!$          IF(.NOT.TCONV) THEN
!!$            CALL ERROR$MSG('LOOP NOT CONVERGED')
!!$            CALL ERROR$STOP('AUGMENTATION_SOFTCORE')
!!$          END IF
!!$          AECORE(:)=AECORE(:)+FB(IB)*PHI(:)**2
!!$        ENDDO
!!$!
!!$!       ==============================================================
!!$!       ==  CONSTRUCT DENSITY                                       ==
!!$!       ==============================================================
!!$        RHO(:)=RHO0(:)+AECORE(:)-AECORE0(:)
!!$!
!!$!       ==============================================================
!!$!       ==  CONSTRUCT POTENTIAL                                     ==
!!$!       ==============================================================
!!$
!!$
!!$
!!$      ENDDO
!!$      RETURN
!!$      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADDVQLM(GID,NR,LMRX,VQLM,AEPOT,PSPOT)
!     ******************************************************************
!     **                                                             **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: GID
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(IN)   :: VQLM(LMRX)
      REAL(8)   ,INTENT(INOUT):: AEPOT(NR,LMRX)
      REAL(8)   ,INTENT(INOUT):: PSPOT(NR,LMRX)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: SVAR
      INTEGER(4)              :: LM,L
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        SVAR=VQLM(LM)
        AEPOT(:,LM)=AEPOT(:,LM)+SVAR*R(:)**L
        PSPOT(:,LM)=PSPOT(:,LM)+SVAR*R(:)**L
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADDBACKGROUND(GID,NR,RHOB &
     &                        ,RHO,EB,POT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ONE-CENTER CONTRIBUTION OF THE COMPENSATING  **
!     **  CHARGE BACKGROUND TO TOTAL ENERGY AND POTENTIAL             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: GID
      INTEGER(4),INTENT(IN)   :: NR
      REAL(8)   ,INTENT(IN)   :: RHOB
      REAL(8)   ,INTENT(IN)   :: RHO(NR)
      REAL(8)   ,INTENT(OUT)  :: EB
      REAL(8)   ,INTENT(INOUT):: POT(NR)
      REAL(8)                 :: PI
      REAL(8)                 :: SVAR
      REAL(8)                 :: WORK(NR)
      REAL(8)                 :: R(NR)
!     ******************************************************************
      IF(RHOB.EQ.0.D0) THEN
        EB=0.D0
        RETURN
      END IF
      PI=4.D0*ATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      SVAR=-2.D0*PI*RHOB/3.D0*SQRT(4.D0*PI)
      WORK(:)=SVAR*RHO(:)*R(:)**4
      POT(:)=POT(:)+SVAR*R(:)**2 ! POTENTIAL OF THE BACKGROUND
      CALL RADIAL$INTEGRAL(GID,NR,WORK,EB)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_HARTREE(IAT,GID,NR,AEZ,RHOCOR &
     &            ,VADD,RCSM,LMRX,AERHO,PSRHO,QLM,AEPOT,PSPOT,VQLM &
     &            ,AEEHARTREE,PSEHARTREE)
!     ******************************************************************
!     **                                                              **
!     ** CALCULATES HARTREE ENERGY                                    **
!     ** ADDS COMPENSATION DENSITY TO THE PSEUDO DENSITY              **
!     **                                                              **
!     **  1) AEE=AERHO*AEZ/R                                          **
!     **     PSE=PSRHO*VADD                                           **
!     **  2) PSRHO=PSRHO+GAUSS(R)*QLM                                 **
!     **  3) AEE=AEE + 0.5*AERHO(R)*AERHO(R')/|R-R'|                  **
!     **     PSE=PSE + 0.5*PSRHO(R)*PSRHO(R')/|R-R'|                  **
!     **     VQLM=VQLM - GAUSS(R)*PSRHO(R')/|R-R'|                    **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: IAT
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      REAL(8)    ,INTENT(IN) :: AEZ
      REAL(8)    ,INTENT(IN) :: RHOCOR(NR)
      REAL(8)    ,INTENT(IN) :: VADD(NR)
      REAL(8)    ,INTENT(IN) :: RCSM
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)    ,INTENT(INOUT) :: PSRHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(OUT):: VQLM(LMRX)
      REAL(8)    ,INTENT(OUT):: AEPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: PSPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: AEEHARTREE
      REAL(8)    ,INTENT(OUT):: PSEHARTREE
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX1(NR)
      REAL(8)                :: AEDMU(NR)
      REAL(8)                :: PSDMU(NR)
      REAL(8)                :: AEE(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      REAL(8)                :: GAUSSIAN(NR)
      REAL(8)                :: AEH,PSH
      INTEGER(4)             :: LM,IR
      INTEGER(4)             :: L
!     == ARRAYS NEEDED DETAILED REPORT
      LOGICAL(4) ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4)             :: NFILO
!     ******************************************************************
      CALL ERROR$MSG('THIS ROUTINE IS APPARENTLY NOT USED')
      CALL ERROR$MSG('IT IS MARKED FOR DELETION')
      CALL ERROR$STOP('AUGMENTATION_HARTREE')

      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  INITIALIZE ENERGIES AND POTENTIALS TO  ZERO                 ==
!     ==================================================================
      AEEHARTREE=0.D0
      PSEHARTREE=0.D0
      VQLM(:)=0.D0
      AEPOT(:,:)=0.D0
      PSPOT(:,:)=0.D0
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC VALENCE-CORE INTERACTION           ==
!     ==   AND VALENCE-NUCLEUS INTERACTION                            ==
!     ==================================================================
      AEH=0.D0
      PSH=0.D0
      AEE(:)=0.D0
      PSE(:)=0.D0
      AEDMU(:)=0.D0
      CALL RADIAL$POISSON(GID,NR,0,RHOCOR,AEDMU)
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX1)
      AEDMU(:)  = AEDMU(:)+AUX1(:)
      AEE(:)    = AEDMU(:)* AERHO(:,1) *R(:)**2
      AEPOT(:,1)= AEPOT(:,1)+AEDMU(:)
      PSE(:)    = VADD(:) * PSRHO(:,1) *R(:)**2
      PSPOT(:,1)= PSPOT(:,1)+VADD(:)
      CALL RADIAL$INTEGRAL(GID,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(GID,NR,PSE,PSH)
!PRINT*,'AEH C-V INT',AEH
      AEEHARTREE=AEEHARTREE+AEH
      PSEHARTREE=PSEHARTREE+PSH
!
!     ==================================================================
!     ==   ADD COMPENSATION CHARGE DENSITY                            ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL GAUSSN(L,ALPHA,CL)
        SVAR=QLM(LM)*CL
        PSRHO(:,LM)=PSRHO(:,LM)+SVAR*EXP(-ALPHA*R(:)**2)*R(:)**L
      ENDDO
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      AEE(:)=0.D0
      PSE(:)=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL RADIAL$POISSON(GID,NR,L,AERHO(1,LM),AEDMU)
        CALL RADIAL$POISSON(GID,NR,L,PSRHO(1,LM),PSDMU)
        ALPHA=1.D0/RCSM**2
        CALL GAUSSN(L,ALPHA,CL)
!       == ADD COMPENSATION DENSITY AND ITS POTENTIAL ================
        GAUSSIAN(:)=CL*EXP(-ALPHA*R(:)**2)*R(:)**L
!
!       == CALCULATE TOTAL ENERGY AND POTENTIAL ======================
        AEE(:)=AEE(:)+0.5D0*AEDMU(:)*AERHO(:,LM)*R(:)**2
        PSE(:)=PSE(:)+0.5D0*PSDMU(:)*PSRHO(:,LM)*R(:)**2
        AEPOT(:,LM)=AEPOT(:,LM)+AEDMU(:)
        PSPOT(:,LM)=PSPOT(:,LM)+PSDMU(:)
        AUX1(:)=PSDMU(:)*GAUSSIAN*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX1,SVAR)
        VQLM(LM)=VQLM(LM)-SVAR
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(GID,NR,PSE,PSH)
      AEEHARTREE=AEEHARTREE+AEH
      PSEHARTREE=PSEHARTREE+PSH
!PRINT*,'AEH V-V INT',AEH
!
!     ==================================================================
!     ==  PRINTOUT FOR TESTING                                        ==
!     ==================================================================
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'AE POTENTIAL FOR ATOM ',IAT
        DO IR=1,NR,49
          WRITE(NFILO,FMT='(9F10.5)')(AEPOT(IR,LM),LM=1,LMRX)
        ENDDO
        WRITE(NFILO,*)'PS POTENTIAL FOR ATOM ',IAT
        DO IR=1,NR,49
          WRITE(NFILO,FMT='(9F10.5)')(PSPOT(IR,LM),LM=1,LMRX)
        ENDDO
        WRITE(NFILO,*)'PS DENSITY FOR ATOM ',IAT
        DO IR=1,NR,49
          WRITE(NFILO,FMT='(9F10.5)')(PSRHO(IR,LM),LM=1,LMRX)
        ENDDO
      ENDIF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DA1,DA)
!     **                                                              **
!     **  CONVERTS THE MATRIX ELEMENTS OF A SPHERICALLY SYMMETRIC     **
!     **  OPERATOR INTO THE FORM OF THE DENSITY MATRIX                **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LMNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: DA1(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: DA(LMNX,LMNX,NDIMD)
      INTEGER(4)            :: LMN1,LMN2,LN1,LN2,L1,L2,IM
!     ******************************************************************      
      DA(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            DO IM=1,2*L1+1
              DA(LMN1+IM,LMN2+IM,1)=DA1(LN1,LN2)
            ENDDO
          END IF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI(IR,:))).GT.1.D+3) EXIT
        WRITE(100,FMT='(F15.10,2X,20(E25.10,2X))')R(IR),PHI(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!***********************************************************************
!***********************************************************************
!**                                                                    **
!**  SCAN0R PROJECT: SCAN-L (DEORBITALIZED META-GGA) EXCHANGE          **
!**                                                                    **
!**  MODULES SCAN_X_MODULE AND DEORB_ALPHA_MODULE, AND SUBROUTINE      **
!**  AUGMENTATION_XC_SCANL BELOW, ARE APPENDED HERE (RATHER THAN IN A  **
!**  NEW FILE) SO NO CHANGES TO THE BUILD SYSTEM'S FILE LIST/          **
!**  DEPENDENCY GRAPH ARE NEEDED -- SEE THE SCAN0R DESIGN NOTES FOR THE**
!**  FULL RATIONALE AND FOR THE STANDALONE-VALIDATED PROTOTYPES THIS   **
!**  IS COPIED FROM (01_scan_exchange_prototype/,                      **
!**  04_deorbitalized_alpha_prototype/).                                **
!**                                                                    **
!**  ACTIVATED BY DFT TYPE=20 (see paw_dft.f90) VIA AUGMENTATION_XC'S  **
!**  EARLY-RETURN TSCANL BRANCH ABOVE. NO EFFECT ON ANY OTHER          **
!**  FUNCTIONAL TYPE.                                                  **
!**                                                                    **
!***********************************************************************
!***********************************************************************
!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  MODULE SCAN_X_MODULE                                              **
!**                                                                    **
!**  PURPOSE: EXCHANGE ENHANCEMENT FACTOR OF THE SCAN META-GGA         **
!**    SUN, RUZSINSZKY, PERDEW, PHYS. REV. LETT. 115, 036402 (2015)    **
!**                                                                    **
!**  STATUS: VALUE-ONLY (NO ANALYTIC DERIVATIVES). AUGMENTATION_XC_    **
!**    SCANL BELOW USES FINITE DIFFERENCES FOR THE POTENTIAL.          **
!**                                                                    **
!**  PROVENANCE OF THE COEFFICIENTS BELOW:                             **
!**    C1X, C2X, DX, K1 CROSS-CHECKED DIRECTLY AGAINST THE LIBXC       **
!**    SOURCE (src/mgga_x_scan.c, SCAN PARAMETER SET).                 **
!**    A1, MU_GE, B1..B4 TAKEN FROM SEARCH-ENGINE-SUMMARIZED LITERATURE**
!**    (CONSISTENT ACROSS TWO INDEPENDENT QUERIES) AND CROSS-CHECKED   **
!**    ALGEBRAICALLY AGAINST THE EXACT UNIFORM-ELECTRON-GAS LIMIT      **
!**    FX(S=0,ALPHA=1)=1 -- SEE SCAN0r notes     **
!**    scan_x_test.f90.                                                **
!**                                                                    **
!**    INDEPENDENTLY VERIFIED (SCAN0r notes            **
!**    verify_scanx_vs_libxc.f90) AGAINST THE REAL, INSTALLED LIBXC     **
!**    7.0.0 mgga_x_scan (func_id=263): 42 (S,ALPHA) POINTS SPANNING    **
!**    S=0..3, ALPHA=0..2.5 (SINGLE-ORBITAL THROUGH WELL PAST THE       **
!**    UNIFORM-GAS LIMIT), CONSTRUCTING A CONSISTENT (RHO,SIGMA,TAU)    **
!**    TRIPLE FOR EACH POINT AND CALLING LIBXC'S xc_mgga_exc DIRECTLY   **
!**    (NOT JUST COMPARING A SHARED (S,ALPHA) DEFINITION IN ISOLATION). **
!**    AGREEMENT: MACHINE PRECISION (~1E-15 ABSOLUTE) AT EVERY POINT.   **
!**                                                                    **
!***********************************************************************
MODULE SCAN_X_MODULE
IMPLICIT NONE
LOGICAL(4)         :: TINI=.FALSE.
REAL(8)            :: A1        ! GX(S) PARAMETER
REAL(8)            :: MUGE      ! SECOND-ORDER GRADIENT-EXPANSION COEFFICIENT =10/81
REAL(8)            :: K1        ! H1X INTERPOLATION PARAMETER
REAL(8)            :: B1,B2,B3,B4 ! H1X(S,ALPHA) "X" FUNCTION PARAMETERS
REAL(8)            :: C1X,C2X,DX  ! FX(ALPHA) SWITCHING-FUNCTION PARAMETERS
REAL(8),PARAMETER  :: FX1=1.174D0 ! SINGLE-ORBITAL-LIMIT PREFACTOR (SUN ET AL. EQ.12)
!     == R2SCAN (FURNESS ET AL., J.PHYS.CHEM.LETT 11, 8208 (2020) AND ERRATUM ==
!     == 11, 9248 (2020)). TR2SCAN IS SET BY SCAN$SETR2SCAN (DFT TYPE 22). ===
!     == ALL CALLERS PASS S AND THE ALPHA OF SCAN; ALPHA-BAR FOLLOWS FROM ====
!     == TAU_W/TAU_UNIF=(5/3)S**2:  ALPHA-BAR=ALPHA/(1+ETA*(5/3)S**2) =======
LOGICAL(4)         :: TR2SCAN=.FALSE.
REAL(8),PARAMETER  :: R2ETA=1.D-3    ! REGULARIZATION OF ALPHA-BAR
REAL(8),PARAMETER  :: R2DP2=0.361D0  ! DAMPING OF THE GRADIENT-EXPANSION TERMS
!     == INTERPOLATION POLYNOMIAL FOR 0<=ALPHA-BAR<=2.5, EXCHANGE =============
REAL(8),PARAMETER  :: R2CX(0:7)=(/1.D0,-0.667D0,-0.4445555D0,-0.663086601049D0 &
     &                  ,1.451297044490D0,-0.887998041597D0,0.234528941479D0 &
     &                  ,-0.023185843322D0/)
CONTAINS
!      .................................................................
       SUBROUTINE SCAN_X_INITIALIZE
       IMPLICIT NONE
!      *****************************************************************
       C1X=0.667D0
       C2X=0.8D0
       DX=1.24D0
       K1=0.065D0
       A1=4.9479D0
       MUGE=10.D0/81.D0
       B3=0.5D0
       B2=SQRT(5913.D0/405000.D0)
       B1=(511.D0/13500.D0)/(2.D0*B2)
       B4=MUGE**2/K1-1606.D0/18225.D0-B1**2
       TINI=.TRUE.
       RETURN
       END SUBROUTINE SCAN_X_INITIALIZE
!
!      .................................................................
       SUBROUTINE SCAN_X_GX(S,G)
!      *****************************************************************
!      **  GX(S)=1-EXP(-A1/SQRT(S))                                   **
!      **  AS S->0+, A1/SQRT(S)->+INFINITY SO GX->1. THE BRANCH BELOW **
!      **  ONLY AVOIDS THE SQRT(0) DIVISION, IT DOES NOT CHANGE THE   **
!      **  LIMIT (GX->1, NOT 0).                                      **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: S
       REAL(8),INTENT(OUT):: G
!      *****************************************************************
       IF(.NOT.TINI) CALL SCAN_X_INITIALIZE
       IF(S.LT.1.D-12) THEN
         G=1.D0
       ELSE
         G=1.D0-EXP(-A1/SQRT(S))
       END IF
       RETURN
       END SUBROUTINE SCAN_X_GX
!
!      .................................................................
       SUBROUTINE SCAN_X_H1X(S,ALPHA,H)
!      *****************************************************************
!      **  H1X(S,ALPHA)=1+K1*X/(K1+X)                                 **
!      **  APPROXIMATE RESUMMATION OF THE 4TH-ORDER GRADIENT EXPANSION**
!      **  FOR EXCHANGE (SUN ET AL. EQ.9-11)                          **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: S
       REAL(8),INTENT(IN) :: ALPHA
       REAL(8),INTENT(OUT):: H
       REAL(8)            :: S2,X,TERM1,TERM2,OMA
!      *****************************************************************
       IF(.NOT.TINI) CALL SCAN_X_INITIALIZE
       S2=S*S
       OMA=1.D0-ALPHA
       TERM1=MUGE*S2*(1.D0+B4*S2/MUGE*EXP(-B4*S2/MUGE))
       TERM2=(B1*S2+B2*OMA*EXP(-B3*OMA**2))**2
       X=TERM1+TERM2
       H=1.D0+K1*X/(K1+X)
       RETURN
       END SUBROUTINE SCAN_X_H1X
!
!      .................................................................
       SUBROUTINE SCAN_X_FALPHA(ALPHA,F)
!      *****************************************************************
!      **  ALPHA-INTERPOLATION/SWITCHING FUNCTION (SUN ET AL. EQ.13)  **
!      **  F(ALPHA)=EXP(-C1X*ALPHA/(1-ALPHA))            FOR ALPHA<1  **
!      **  F(ALPHA)=-DX*EXP( C2X/(1-ALPHA))               FOR ALPHA>1 **
!      **  F(1)=0 (CONTINUOUS FROM BOTH SIDES)                        **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: ALPHA
       REAL(8),INTENT(OUT):: F
       REAL(8),PARAMETER  :: TOL=1.D-8
!      *****************************************************************
       IF(ALPHA.LT.1.D0-TOL) THEN
         F=EXP(-C1X*ALPHA/(1.D0-ALPHA))
       ELSE IF(ALPHA.GT.1.D0+TOL) THEN
         F=-DX*EXP(C2X/(1.D0-ALPHA))
       ELSE
         F=0.D0   ! LIMIT FROM BOTH SIDES AT ALPHA=1
       END IF
       RETURN
       END SUBROUTINE SCAN_X_FALPHA
!
!      .................................................................
       SUBROUTINE SCAN_X_FX(S,ALPHA,FX)
!      *****************************************************************
!      **  FULL SCAN EXCHANGE ENHANCEMENT FACTOR                      **
!      **  FX(S,ALPHA)={H1X(S,ALPHA)+F(ALPHA)*[1.174-H1X(S,ALPHA)]}   **
!      **               *GX(S)                                       **
!      **  E_X^SCAN=INTEGRAL{ N * EPS_X^UNIF(N) * FX(S,ALPHA) }       **
!      **  (EPS_X^UNIF ALREADY IMPLEMENTED IN EXCHANGE_X1,             **
!      **   PAW_DFT.F90:1867 -- REUSABLE UNCHANGED, SEE ITS USE VIA    **
!      **   EXFAC BELOW)                                              **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: S      ! REDUCED DENSITY GRADIENT
       REAL(8),INTENT(IN) :: ALPHA  ! ISO-ORBITAL INDICATOR (TAU-TAU_W)/TAU_UNIF
       REAL(8),INTENT(OUT):: FX
       REAL(8)            :: G,H,F
!      *****************************************************************
       IF(TR2SCAN) THEN
         CALL R2SCAN_X_FX(S,ALPHA,FX)
         RETURN
       END IF
       CALL SCAN_X_GX(S,G)
       CALL SCAN_X_H1X(S,ALPHA,H)
       CALL SCAN_X_FALPHA(ALPHA,F)
       FX=(H+F*(FX1-H))*G
       RETURN
       END SUBROUTINE SCAN_X_FX
!
!      .................................................................
       SUBROUTINE R2SCAN_FALPHA(A,C,C2,D,F)
!      *****************************************************************
!      **  INTERPOLATION FUNCTION OF RSCAN/R2SCAN IN ALPHA-BAR:       **
!      **    F=EXP(-C1*A/(1-A)), C1=-C(1)    FOR A<0                  **
!      **    F=SUM_I C(I)*A**I               FOR 0<=A<=2.5            **
!      **    F=-D*EXP(C2/(1-A))              FOR A>2.5                **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: A
       REAL(8),INTENT(IN) :: C(0:7)
       REAL(8),INTENT(IN) :: C2,D
       REAL(8),INTENT(OUT):: F
       INTEGER(4)         :: I
!      *****************************************************************
       IF(A.LT.0.D0) THEN
         F=EXP(C(1)*A/(1.D0-A))
       ELSE IF(A.LE.2.5D0) THEN
         F=C(7)
         DO I=6,0,-1
           F=F*A+C(I)
         ENDDO
       ELSE
         F=-D*EXP(C2/(1.D0-A))
       END IF
       RETURN
       END SUBROUTINE R2SCAN_FALPHA
!
!      .................................................................
       SUBROUTINE R2SCAN_X_FX(S,ALPHA,FX)
!      *****************************************************************
!      **  R2SCAN EXCHANGE ENHANCEMENT FACTOR                         **
!      **    FX={H1X(P)+F(ALPHA-BAR)*[1.174-H1X(P)]}*GX(P), P=S**2    **
!      **    H1X(P)=1+K1-K1/(1+P*(MU+C_ETA*C2*EXP(-P**2/DP2**4))/K1)  **
!      **    C2=-SUM_{I=1..7} I*C(I)*(1-1.174), C_ETA=20/27+ETA*5/3   **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: S      ! REDUCED DENSITY GRADIENT
       REAL(8),INTENT(IN) :: ALPHA  ! (TAU-TAU_W)/TAU_UNIF AS FOR SCAN
       REAL(8),INTENT(OUT):: FX
       REAL(8)            :: P,ABAR,DELF2,C2,CETA,DAMP,H,F,G
       INTEGER(4)         :: I
!      *****************************************************************
       IF(.NOT.TINI) CALL SCAN_X_INITIALIZE
       P=S*S
       ABAR=ALPHA/(1.D0+R2ETA*5.D0/3.D0*P)
       DELF2=0.D0
       DO I=1,7
         DELF2=DELF2+DBLE(I)*R2CX(I)
       ENDDO
       C2=-DELF2*(1.D0-FX1)
       CETA=20.D0/27.D0+R2ETA*5.D0/3.D0
       DAMP=EXP(-P**2/R2DP2**4)
       H=1.D0+K1-K1/(1.D0+P*(MUGE+CETA*C2*DAMP)/K1)
       CALL R2SCAN_FALPHA(ABAR,R2CX,C2X,DX,F)
       CALL SCAN_X_GX(S,G)
       FX=(H+F*(FX1-H))*G
       RETURN
       END SUBROUTINE R2SCAN_X_FX
!
END MODULE SCAN_X_MODULE
!
!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  MODULE DEORB_ALPHA_MODULE                                         **
!**                                                                    **
!**  PURPOSE: DEORBITALIZED (LAPLACIAN-LEVEL) APPROXIMATION TO THE     **
!**    SCAN ISO-ORBITAL INDICATOR ALPHA=(TAU-TAU_W)/TAU_UNIF, BUILT    **
!**    FROM N, |GRAD N|^2, LAPLACIAN(N) ONLY -- NO ORBITALS NEEDED.    **
!**                                                                    **
!**  METHOD: PERDEW-CONSTANTIN (PC) INTERPOLATION BETWEEN THE          **
!**    MODIFIED 4TH-ORDER GRADIENT EXPANSION OF TAU (MGE4) AND THE     **
!**    VON WEIZSAECKER LIMIT, AS USED FOR DEORBITALIZED META-GGAS      **
!**    (SCAN-L) BY MEJIA-RODRIGUEZ & TRICKEY,                          **
!**    PHYS. REV. A 96, 052512 (2017), ARXIV:1710.06032, EQS. 6-8,     **
!**    14, 24-29. ORIGINAL SWITCHING FUNCTION: PERDEW & CONSTANTIN,    **
!**    PHYS. REV. B 75, 155109 (2007).                                 **
!**                                                                    **
!**  PROVENANCE / VERIFICATION NOTES:                                  **
!**    - S, Q, T_TF, T_W, F_T^(0,2,4), F_T^MGE4, THE PC INTERPOLATION  **
!**      STRUCTURE (Z^PC, THETA_PC, A=0.5389, B=3) ARE TAKEN FROM A    **
!**      LITERATURE-SEARCH EXTRACTION OF ARXIV:1710.06032 (NOT READ    **
!**      DIRECTLY FROM THE PDF -- EXTRACTION TOOLING COULD NOT PARSE   **
!**      IT, SO THIS WENT THROUGH AN INTERMEDIATE SUMMARIZATION STEP). **
!**    - THE THETA_PC FORMULA WAS EXTRACTED TWICE AND DISAGREED ON     **
!**      BRACKET PLACEMENT BETWEEN ATTEMPTS; THE SECOND (USED HERE)    **
!**      WAS KEPT BECAUSE IT IS THE ONE THAT SATISFIES CONTINUITY AT   **
!**      Z=0 AND Z=A, VERIFIED BY HAND -- NOT BECAUSE THE SECOND       **
!**      FETCH IS INHERENTLY MORE TRUSTWORTHY.                         **
!**    - F_T^W=(5/3)S^2 (NO ADDITIVE CONSTANT) WAS CROSS-CHECKED BY    **
!**      REQUIRING THE UNIFORM-GAS LIMIT ALPHA(S=0,Q=0)=1 TO COME OUT  **
!**      RIGHT -- IT ONLY DOES WITHOUT THE "+1" A CANDIDATE EXTRACTION **
!**      INCLUDED. SEE SCAN0r notes        **
!**      alpha_test.f90.                                               **
!**    - VERIFIED 2026-10-02 AGAINST LIBXC 7.1.2 (MGGA_X_SCANL): THE   **
!**      SCAN-L EXCHANGE ENERGY DENSITY AGREES TO 1E-15 (RELATIVE).    **
!**                                                                    **
!***********************************************************************
MODULE DEORB_ALPHA_MODULE
IMPLICIT NONE
LOGICAL(4)         :: TINI=.FALSE.
REAL(8)            :: PI
REAL(8)            :: CTF        ! T_TF=CTF*RHO**(5/3); CTF=(3/10)(3PI^2)**(2/3)
REAL(8)            :: SFAC       ! S=SFAC*SQRT(GRHO2)/RHO**(4/3)
REAL(8)            :: QFAC       ! Q=QFAC*LAPRHO/RHO**(5/3)
!== PC07OPT PARAMETERS OF MEJIA-RODRIGUEZ & TRICKEY (2017) FOR SCAN-L, AS IN
!== LIBXC (MGGA_X_SCANL, MGGA_C_SCANL) AND VASP (METAGGA=SCANL). THE ORIGINAL
!== PERDEW-CONSTANTIN VALUES A=0.5389, B=3 WERE USED BEFORE 2026-10-02.
REAL(8),PARAMETER  :: A_PC=1.784720D0  ! PC SWITCHING FUNCTION PARAMETER (EQ.27)
REAL(8),PARAMETER  :: B_PC=0.258304D0  ! PC SWITCHING FUNCTION EXPONENT   (EQ.27)
REAL(8),PARAMETER  :: EXPARGMAX=700.D0 ! OVERFLOW GUARD FOR EXP()
CONTAINS
!      .................................................................
       SUBROUTINE DEORB_INITIALIZE
       IMPLICIT NONE
!      *****************************************************************
       PI=4.D0*ATAN(1.D0)
       CTF=0.3D0*(3.D0*PI**2)**(2.D0/3.D0)
       SFAC=1.D0/(2.D0*(3.D0*PI**2)**(1.D0/3.D0))
       QFAC=1.D0/(4.D0*(3.D0*PI**2)**(2.D0/3.D0))
       TINI=.TRUE.
       RETURN
       END SUBROUTINE DEORB_INITIALIZE
!
!      .................................................................
       SUBROUTINE DEORB_SQ(RHO,GRHO2,LAPRHO,S,Q)
!      *****************************************************************
!      **  REDUCED DENSITY GRADIENT S (EQ.8) AND REDUCED LAPLACIAN Q  **
!      **  (EQ.14) OF MEJIA-RODRIGUEZ & TRICKEY (2017).               **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHO     ! DENSITY (SPIN-UNPOLARIZED)
       REAL(8),INTENT(IN) :: GRHO2   ! |GRAD RHO|**2
       REAL(8),INTENT(IN) :: LAPRHO  ! LAPLACIAN(RHO)
       REAL(8),INTENT(OUT):: S
       REAL(8),INTENT(OUT):: Q
       REAL(8)            :: RHO1
!      *****************************************************************
       IF(.NOT.TINI) CALL DEORB_INITIALIZE
       RHO1=MAX(RHO,1.D-12)
       S=SFAC*SQRT(MAX(GRHO2,0.D0))/RHO1**(4.D0/3.D0)
       Q=QFAC*LAPRHO/RHO1**(5.D0/3.D0)
       RETURN
       END SUBROUTINE DEORB_SQ
!
!      .................................................................
       SUBROUTINE DEORB_THETAPC(Z,THETA)
!      *****************************************************************
!      **  PERDEW-CONSTANTIN SWITCHING FUNCTION (EQ.27):              **
!      **  THETA_PC(Z)=0                                 Z<=0         **
!      **  THETA_PC(Z)=[(1+EXP(A/(A-Z)))/(EXP(A/Z)+EXP(A/(A-Z)))]**B  **
!      **                                                 0<Z<A       **
!      **  THETA_PC(Z)=1                                  Z>=A        **
!      **  CONTINUITY CHECKED BY HAND AT BOTH ENDPOINTS (SEE HEADER). **
!      **  OVERFLOW OF EXP() NEAR THE ENDPOINTS IS GUARDED EXPLICITLY **
!      **  USING THE SAME LIMITING VALUES.                            **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: Z
       REAL(8),INTENT(OUT):: THETA
       REAL(8)            :: ARG1,ARG2,NUM,DEN,RATIO
!      *****************************************************************
       IF(Z.LE.0.D0) THEN
         THETA=0.D0
       ELSE IF(Z.GE.A_PC) THEN
         THETA=1.D0
       ELSE
         ARG1=A_PC/Z          ! -> +INFTY AS Z->0+
         ARG2=A_PC/(A_PC-Z)   ! -> +INFTY AS Z->A-
         IF(ARG1.GT.EXPARGMAX) THEN
           THETA=0.D0         ! DENOMINATOR DOMINATED BY EXP(ARG1)->INFTY
         ELSE IF(ARG2.GT.EXPARGMAX) THEN
           THETA=1.D0         ! NUM. AND DEN. BOTH DOMINATED BY EXP(ARG2)
         ELSE
           NUM=1.D0+EXP(ARG2)
           DEN=EXP(ARG1)+EXP(ARG2)
           RATIO=NUM/DEN
           THETA=RATIO**B_PC
         END IF
       END IF
       RETURN
       END SUBROUTINE DEORB_THETAPC
!
!      .................................................................
       SUBROUTINE DEORB_ALPHA(S,Q,ALPHA)
!      *****************************************************************
!      **  DEORBITALIZED ISO-ORBITAL INDICATOR ALPHA_PC(S,Q).         **
!      **  DERIVATION (NOT STATED EXPLICITLY IN THE SOURCE PAPER,     **
!      **  WORKED OUT HERE FROM ITS DEFINITIONS):                     **
!      **    F_T := TAU/T_TF ,  F_T^W := TAU_W/T_TF = (5/3)*S**2      **
!      **    ALPHA := (TAU-TAU_W)/T_TF = F_T - F_T^W                  **
!      **    F_T^PC = F_T^W + Z^PC*THETA_PC(Z^PC)  (EQ.28)            **
!      **    => ALPHA_PC = F_T^PC-F_T^W = Z^PC*THETA_PC(Z^PC)         **
!      **  (T_TF IS THE SAME "TAU_UNIF" USED IN THE STANDARD SCAN     **
!      **  DEFINITION OF ALPHA, SO THIS ALPHA_PC DROPS STRAIGHT INTO  **
!      **  SCAN_X_MODULE'S SCAN_X_FX(S,ALPHA).)                       **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: S
       REAL(8),INTENT(IN) :: Q
       REAL(8),INTENT(OUT):: ALPHA
       REAL(8)            :: S2,Q2,S4
       REAL(8)            :: FTW,FT2,FT4,FTMGE4,Z,THETA
!      *****************************************************************
       S2=S*S
       S4=S2*S2
       Q2=Q*Q
       FTW=(5.D0/3.D0)*S2                                     ! EQ.7/T_TF
       FT2=(5.D0/27.D0)*S2+(20.D0/9.D0)*Q                     ! EQ.25
       FT4=(8.D0/81.D0)*Q2-(1.D0/9.D0)*S2*Q+(8.D0/243.D0)*S4  ! EQ.26
       FTMGE4=(1.D0+FT2+FT4)/SQRT(1.D0+(FT4/(1.D0+FTW))**2)   ! EQ.24
       Z=FTMGE4-FTW                                           ! EQ.29
       CALL DEORB_THETAPC(Z,THETA)
       ALPHA=Z*THETA                                          ! SEE HEADER DERIVATION
       RETURN
       END SUBROUTINE DEORB_ALPHA
!
END MODULE DEORB_ALPHA_MODULE
!
!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  MODULE SCAN_C_MODULE                                              **
!**                                                                    **
!**  PURPOSE: SCAN META-GGA CORRELATION ENERGY PER PARTICLE, SPIN-     **
!**    UNPOLARIZED (ZETA=0) ONLY -- MATCHES THE SPIN SCOPE CUT ALREADY **
!**    MADE FOR SCAN_X_MODULE/AUGMENTATION_XC_SCANL (SEE THE SCAN0R    **
!**    DESIGN DOC). ORIGINAL 2015 SCAN (NOT RSCAN/R2SCAN -- IDELEC=0,  **
!**    IINTERP=0, IALPHA=0 IN THE SOURCE THIS WAS TRANSCRIBED FROM).   **
!**    SUN, RUZSINSZKY, PERDEW, PHYS. REV. LETT. 115, 036402 (2015).   **
!**                                                                    **
!**  STATUS: VALUE-ONLY (NO ANALYTIC DERIVATIVES), SAME AS             **
!**    SCAN_X_MODULE.                                                  **
!**                                                                    **
!**  PROVENANCE (DIFFERENT SOURCING ROUTE THAN SCAN_X_MODULE/          **
!**  DEORB_ALPHA_MODULE -- THOSE WENT THROUGH A WEB-SEARCH SUMMARIZER  **
!**  OF THE PDF, WHICH FAILED REPEATEDLY FOR CORRELATION. THIS WENT    **
!**  THROUGH XCFUN (github.com/dftlibs/xcfun), A HAND-WRITTEN,         **
!**  HUMAN-READABLE SYMBOLIC C++ XC LIBRARY -- PLAIN SOURCE CODE,      **
!**  NOT A PDF, FETCHED AND READ DIRECTLY (src/functionals/            **
!**  SCAN_like_eps.hpp, pw92eps.hpp). CROSS-CHECKS THAT INCREASE       **
!**  CONFIDENCE RELATIVE TO THE EXCHANGE MODULE'S SOURCING:            **
!**    - GAMMA=0.031090690869655 MATCHES THE STANDARD PBE CONSTANT     **
!**      (1-LN2)/PI**2 TO THE FULL 12 DIGITS QUOTED.                   **
!**    - THE PW92 UNPOLARIZED PARAMETERS (A,ALPHA1,BETA1..4) MATCH     **
!**      INDEPENDENTLY-RECALLED STANDARD PW92 VALUES EXACTLY.          **
!**    - CFC1=0.64, CFC2=1.5, CFDC1=0.7 (THE ALPHA-SWITCHING FUNCTION, **
!**      SAME FORM AS SCAN_X_MODULE'S FALPHA) MATCHED A SEPARATE,      **
!**      EARLIER WEB-SEARCH RESULT INDEPENDENTLY.                      **
!**    *** STILL NOT VERIFIED AGAINST A REFERENCE NUMERICAL            **
!**    *** IMPLEMENTATION (E.G. RUNNING ACTUAL LIBXC/XCFUN AND         **
!**    *** COMPARING POINTWISE). TWO EXACT-LIMIT CHECKS ARE DONE       **
!**    *** BELOW (scan_c_test.f90) BUT THAT IS NOT THE SAME AS         **
!**    *** POINTWISE VALIDATION ACROSS THE FULL (S,T,ALPHA) DOMAIN.    **
!**                                                                    **
!**  INPUTS: RS (WIGNER-SEITZ RADIUS), S (STANDARD REDUCED DENSITY     **
!**    GRADIENT, SAME CONVENTION AS SCAN_X_MODULE/DEORB_ALPHA_MODULE), **
!**    T2 (SQUARED PBE-STYLE CORRELATION REDUCED GRADIENT -- SAME      **
!**    CONVENTION AND CONSTANTS AS THE EXISTING, ALREADY-VALIDATED     **
!**    PBE_MODULE IN PAW_DFT.F90, T2=TFAC*GRHO2*RHOT**(-7/3),          **
!**    TFAC=0.25*RSFAC/KSFAC**2 -- SEE paw_dft.f90 AROUND LINE 5616),  **
!**    ALPHA (SAME DEORBITALIZED ALPHA AS DEORB_ALPHA_MODULE).         **
!**                                                                    **
!***********************************************************************
MODULE SCAN_C_MODULE
USE SCAN_X_MODULE, ONLY: TR2SCAN,R2ETA,R2DP2,R2SCAN_FALPHA
IMPLICIT NONE
LOGICAL(4)         :: TINI=.FALSE.
!     == PW92 UNPOLARIZED LSDA CORRELATION PARAMETERS (STANDARD) ==============
REAL(8),PARAMETER  :: PW92_A =0.03109070D0
REAL(8),PARAMETER  :: PW92_A1=0.21370D0
REAL(8),PARAMETER  :: PW92_B1=7.59570D0
REAL(8),PARAMETER  :: PW92_B2=3.5876D0
REAL(8),PARAMETER  :: PW92_B3=1.63820D0
REAL(8),PARAMETER  :: PW92_B4=0.49294D0
!     == SCAN ALPHA=0 (SINGLE-ORBITAL) BRANCH PARAMETERS =======================
REAL(8),PARAMETER  :: B1C=0.0285764D0
REAL(8),PARAMETER  :: B2C=0.0889D0
REAL(8),PARAMETER  :: B3C=0.125541D0
REAL(8),PARAMETER  :: CHI_LD=0.12802585262625815D0
!     == SCAN ALPHA=1 (GGA) BRANCH: PBE-STYLE H1 CONSTRUCTION WITH AN =========
!     == RS-DEPENDENT BETA (SCAN'S MODIFICATION VS. PLAIN PBE, WHICH USES ======
!     == A CONSTANT BETA) =======================================================
REAL(8),PARAMETER  :: GAMMA_PBE=0.031090690869655D0
REAL(8),PARAMETER  :: BETA_MB=0.06672455060314922D0 ! AS VASP, LIBXC (WAS 0.066725 BEFORE 2026-10-02)
REAL(8),PARAMETER  :: AFACTOR=0.1D0
REAL(8),PARAMETER  :: BFACTOR=0.1778D0
!     == ALPHA-INTERPOLATION SWITCHING FUNCTION (SAME FORM AS SCAN_X_MODULE'S =
!     == FALPHA, DIFFERENT CONSTANTS) ===========================================
REAL(8),PARAMETER  :: CFC1=0.64D0
REAL(8),PARAMETER  :: CFC2=1.5D0
REAL(8),PARAMETER  :: CFDC1=0.7D0
!     == PW92 FULLY POLARIZED (P) AND SPIN-STIFFNESS (A) PARAMETERS ===========
REAL(8),PARAMETER  :: PW92P(6)=(/0.01554535D0,0.20548D0,14.1189D0,6.1977D0 &
     &                          ,3.3662D0,0.62517D0/)
REAL(8),PARAMETER  :: PW92S(6)=(/0.0168869D0,0.11125D0,10.357D0,3.6231D0 &
     &                          ,0.88026D0,0.49671D0/)
REAL(8),PARAMETER  :: FZZ=1.709920934161365617563962776245D0 ! F''(0)
REAL(8),PARAMETER  :: DXC=2.363D0   ! G_C(ZETA) PARAMETER: 2.363 AS VASP (VSCANC) AND LIBXC;
!                                      ! THE SCAN SUPPLEMENT QUOTES 2.3631
!     == R2SCAN CORRELATION: INTERPOLATION POLYNOMIAL FOR 0<=ALPHA-BAR<=2.5 ===
REAL(8),PARAMETER  :: R2CC(0:7)=(/1.D0,-0.64D0,-0.4352D0,-1.535685604549D0 &
     &                  ,3.061560252175D0,-1.915710236206D0,0.516884468372D0 &
     &                  ,-0.051848879792D0/)
REAL(8),PARAMETER  :: PW92U(6)=(/0.0310907D0,0.21370D0,7.5957D0,3.5876D0 &
     &                          ,1.63820D0,0.49294D0/)
CONTAINS
!      .................................................................
       SUBROUTINE SCAN_C_PW92(RS,EU)
!      *****************************************************************
!      **  STANDARD PERDEW-WANG 1992 UNPOLARIZED LSDA CORRELATION     **
!      **  ENERGY PER PARTICLE.                                       **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RS
       REAL(8),INTENT(OUT):: EU
       REAL(8)            :: SQRTRS,POLY
!      *****************************************************************
       SQRTRS=SQRT(RS)
       POLY=SQRTRS*(PW92_B1+SQRTRS*(PW92_B2+SQRTRS*(PW92_B3+PW92_B4*SQRTRS)))
       EU=-2.D0*PW92_A*(1.D0+PW92_A1*RS)*LOG(1.D0+0.5D0/(PW92_A*POLY))
       RETURN
       END SUBROUTINE SCAN_C_PW92
!
!      .................................................................
       SUBROUTINE SCAN_C_EC0(RS,S,EC0)
!      *****************************************************************
!      **  SCAN'S ALPHA=0 (SINGLE-ORBITAL) CORRELATION BRANCH.        **
!      **  USES ITS OWN SIMPLIFIED LDA FORM (B1C,B2C,B3C), NOT PW92.  **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RS
       REAL(8),INTENT(IN) :: S
       REAL(8),INTENT(OUT):: EC0
       REAL(8)            :: ECLDA,W0,GINF,H0
!      *****************************************************************
       ECLDA=-B1C/(1.D0+B2C*SQRT(RS)+B3C*RS)
       W0=EXP(-ECLDA/B1C)-1.D0
       GINF=1.D0/(1.D0+4.D0*CHI_LD*S*S)**0.25D0
       H0=B1C*LOG(1.D0+W0*(1.D0-GINF))
       EC0=ECLDA+H0
       RETURN
       END SUBROUTINE SCAN_C_EC0
!
!      .................................................................
       SUBROUTINE SCAN_C_EC1(RS,T2,EC1)
!      *****************************************************************
!      **  SCAN'S ALPHA=1 (SLOWLY-VARYING/GGA) CORRELATION BRANCH:    **
!      **  PBE-STYLE H1 CONSTRUCTION ON TOP OF PW92, WITH AN          **
!      **  RS-DEPENDENT BETA (SCAN'S MODIFICATION OF PBE, WHICH USES  **
!      **  A CONSTANT BETA).                                          **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RS
       REAL(8),INTENT(IN) :: T2
       REAL(8),INTENT(OUT):: EC1
       REAL(8)            :: EU,W1,BETA,Y,GY,H1
!      *****************************************************************
       CALL SCAN_C_PW92(RS,EU)
       W1=EXP(-EU/GAMMA_PBE)-1.D0
       BETA=BETA_MB*(1.D0+AFACTOR*RS)/(1.D0+BFACTOR*RS)
       Y=BETA/(GAMMA_PBE*W1)*T2
       GY=1.D0/(1.D0+4.D0*Y)**0.25D0
       H1=GAMMA_PBE*LOG(1.D0+W1*(1.D0-GY))
       EC1=EU+H1
       RETURN
       END SUBROUTINE SCAN_C_EC1
!
!      .................................................................
       SUBROUTINE SCAN_C_FALPHA(ALPHA,IEF)
!      *****************************************************************
!      **  ALPHA-INTERPOLATION SWITCHING FUNCTION FOR CORRELATION.    **
!      **  SAME FORM AS SCAN_X_MODULE'S SCAN_X_FALPHA, DIFFERENT      **
!      **  CONSTANTS (CFC1,CFC2,CFDC1 VS. C1X,C2X,DX).                **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: ALPHA
       REAL(8),INTENT(OUT):: IEF
       REAL(8),PARAMETER  :: TOL=1.D-8
!      *****************************************************************
       IF(ALPHA.LT.1.D0-TOL) THEN
         IEF=EXP(-CFC1*ALPHA/(1.D0-ALPHA))
       ELSE IF(ALPHA.GT.1.D0+TOL) THEN
         IEF=-CFDC1*EXP(CFC2/(1.D0-ALPHA))
       ELSE
         IEF=0.D0   ! LIMIT FROM BOTH SIDES AT ALPHA=1 (SAME AS EXCHANGE)
       END IF
       RETURN
       END SUBROUTINE SCAN_C_FALPHA
!
!      .................................................................
       SUBROUTINE SCAN_C_EPSC(RS,S,T2,ALPHA,EPSC)
!      *****************************************************************
!      **  SCAN CORRELATION ENERGY PER PARTICLE OF AN UNPOLARIZED     **
!      **  DENSITY: EPSC = EC1 + IEF(ALPHA)*(EC0-EC1)                 **
!      **  E_C^SCAN = INTEGRAL{ N * EPSC }                            **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RS
       REAL(8),INTENT(IN) :: S
       REAL(8),INTENT(IN) :: T2
       REAL(8),INTENT(IN) :: ALPHA
       REAL(8),INTENT(OUT):: EPSC
!      *****************************************************************
       CALL SCAN_C_EPSC_SPIN(RS,0.D0,S,T2,ALPHA,EPSC)
       RETURN
       END SUBROUTINE SCAN_C_EPSC
!
!      .................................................................
       SUBROUTINE SCAN_C_PW92G(P,RS,G)
!      *****************************************************************
!      **  PW92 INTERPOLATION G(RS)=-2A(1+A1 RS)LN(1+1/(2A POLY))     **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: P(6)   ! A,ALPHA1,BETA1..BETA4
       REAL(8),INTENT(IN) :: RS
       REAL(8),INTENT(OUT):: G
       REAL(8)            :: SQRTRS,POLY
!      *****************************************************************
       SQRTRS=SQRT(RS)
       POLY=SQRTRS*(P(3)+SQRTRS*(P(4)+SQRTRS*(P(5)+P(6)*SQRTRS)))
       G=-2.D0*P(1)*(1.D0+P(2)*RS)*LOG(1.D0+0.5D0/(P(1)*POLY))
       RETURN
       END SUBROUTINE SCAN_C_PW92G
!
!      .................................................................
       SUBROUTINE SCAN_C_PW92GD(P,RS,G,DG)
!      *****************************************************************
!      **  PW92 INTERPOLATION G(RS) AND ITS DERIVATIVE DG/DRS          **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: P(6)   ! A,ALPHA1,BETA1..BETA4
       REAL(8),INTENT(IN) :: RS
       REAL(8),INTENT(OUT):: G,DG
       REAL(8)            :: SQRTRS,Q1,DQ1,LG
!      *****************************************************************
       SQRTRS=SQRT(RS)
       Q1=2.D0*P(1)*SQRTRS*(P(3)+SQRTRS*(P(4)+SQRTRS*(P(5)+P(6)*SQRTRS)))
       DQ1=P(1)*(P(3)/SQRTRS+2.D0*P(4)+3.D0*P(5)*SQRTRS+4.D0*P(6)*RS)
       LG=LOG(1.D0+1.D0/Q1)
       G=-2.D0*P(1)*(1.D0+P(2)*RS)*LG
       DG=-2.D0*P(1)*P(2)*LG+2.D0*P(1)*(1.D0+P(2)*RS)*DQ1/(Q1*(1.D0+Q1))
       RETURN
       END SUBROUTINE SCAN_C_PW92GD
!
!      .................................................................
       SUBROUTINE SCAN_C_EPSC_SPIN(RS,ZETA,S,T2,ALPHA0,EPSC)
!      *****************************************************************
!      **  SCAN CORRELATION ENERGY PER PARTICLE FOR SPIN POLARIZATION **
!      **  ZETA (SUN, RUZSINSZKY, PERDEW, PRL 115, 036402 (2015),     **
!      **  SUPPLEMENT), AS IN VASP AND LIBXC (MGGA_C_SCAN):           **
!      **    EPSC=EC1+F_C(ALPHA)*(EC0-EC1)                            **
!      **    EC1=EPS_PW92(RS,ZETA)+GAMMA PHI**3 LN[1+W1(1-G(A T**2))]  **
!      **    EC0=[EPS_LDA0(RS)+H0(S)]*G_C(ZETA)                        **
!      **  CONVENTIONS OF THE INPUT (ALL REFER TO THE TOTAL DENSITY): **
!      **    S     =|GRAD N|/(2 KF N)                                 **
!      **    T2    =|GRAD N|**2/(2 KS N)**2, I.E. WITHOUT PHI; THE    **
!      **           FACTOR PHI**2 IS APPLIED HERE                     **
!      **    ALPHA0=(TAU-TAU_W)/TAU_UNIF(N), I.E. WITHOUT D_S(ZETA);  **
!      **           THE FACTOR 1/D_S(ZETA) IS APPLIED HERE            **
!      **  FOR ZETA=0 ALL SPIN FACTORS ARE ONE.                       **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RS
       REAL(8),INTENT(IN) :: ZETA
       REAL(8),INTENT(IN) :: S
       REAL(8),INTENT(IN) :: T2
       REAL(8),INTENT(IN) :: ALPHA0
       REAL(8),INTENT(OUT):: EPSC
       REAL(8),PARAMETER  :: ZMAX=1.D0-1.D-12
       REAL(8)            :: Z,OPZ,OMZ,PHI,PHI3,DX,DS,FZ,Z4,GCZ
       REAL(8)            :: EU,EP,EA,ELSDA,W1,BETA,Y,GY,H1,EC1
       REAL(8)            :: ECLDA0,W0,GINF,H0,EC0,ALPHA,IEF
       REAL(8)            :: P,DELF2,DEU,DEP,DEA,DELSDA,DECLDA0,DELY
       INTEGER(4)         :: I
!      *****************************************************************
       Z=MIN(MAX(ZETA,-ZMAX),ZMAX)
       OPZ=1.D0+Z
       OMZ=1.D0-Z
       PHI=0.5D0*(OPZ**(2.D0/3.D0)+OMZ**(2.D0/3.D0))
       DX =0.5D0*(OPZ**(4.D0/3.D0)+OMZ**(4.D0/3.D0))
       DS =0.5D0*(OPZ**(5.D0/3.D0)+OMZ**(5.D0/3.D0))
       PHI3=PHI**3
!
!      == SPIN-POLARIZED PW92 ==========================================
       CALL SCAN_C_PW92(RS,EU)
       FZ=(OPZ**(4.D0/3.D0)+OMZ**(4.D0/3.D0)-2.D0)/(2.D0**(4.D0/3.D0)-2.D0)
       Z4=Z**4
       ELSDA=EU
       IF(Z.NE.0.D0) THEN
         CALL SCAN_C_PW92G(PW92P,RS,EP)
         CALL SCAN_C_PW92G(PW92S,RS,EA)   ! =-ALPHA_C
         ELSDA=EU*(1.D0-FZ*Z4)+EP*FZ*Z4-EA*FZ*(1.D0-Z4)/FZZ
       END IF
!
!      == ALPHA=1 BRANCH ===============================================
       W1=EXP(-ELSDA/(GAMMA_PBE*PHI3))-1.D0
       BETA=BETA_MB*(1.D0+AFACTOR*RS)/(1.D0+BFACTOR*RS)
       Y=BETA/(GAMMA_PBE*W1)*T2/PHI**2
       GY=1.D0/(1.D0+4.D0*Y)**0.25D0
       H1=GAMMA_PBE*PHI3*LOG(1.D0+W1*(1.D0-GY))
       EC1=ELSDA+H1
!
!      == ALPHA=0 BRANCH ===============================================
       ECLDA0=-B1C/(1.D0+B2C*SQRT(RS)+B3C*RS)
       W0=EXP(-ECLDA0/B1C)-1.D0
       GINF=1.D0/(1.D0+4.D0*CHI_LD*S*S)**0.25D0
       H0=B1C*LOG(1.D0+W0*(1.D0-GINF))
       GCZ=(1.D0-DXC*(DX-1.D0))*(1.D0-Z**12)
       EC0=(ECLDA0+H0)*GCZ
!
!      == R2SCAN: H1 WITH THE GRADIENT-EXPANSION CORRECTION DELTA-Y ====
!      == DELY=DF2/(27 GAMMA DS PHI3 W1)*{20 RS [D(GC*ELDA0)/DRS        ==
!      ==   -D(ELSDA)/DRS]-45 ETA [GC*ELDA0-ELSDA]}*P*EXP(-P**2/DP2**4)  ==
       IF(TR2SCAN) THEN
         P=S*S
         CALL SCAN_C_PW92GD(PW92U,RS,EU,DEU)
         DELSDA=DEU
         IF(Z.NE.0.D0) THEN
           CALL SCAN_C_PW92GD(PW92P,RS,EP,DEP)
           CALL SCAN_C_PW92GD(PW92S,RS,EA,DEA)
           DELSDA=DEU*(1.D0-FZ*Z4)+DEP*FZ*Z4-DEA*FZ*(1.D0-Z4)/FZZ
         END IF
         DECLDA0=B1C*(0.5D0*B2C/SQRT(RS)+B3C)/(1.D0+B2C*SQRT(RS)+B3C*RS)**2
         DELF2=0.D0
         DO I=1,7
           DELF2=DELF2+DBLE(I)*R2CC(I)
         ENDDO
         DELY=DELF2/(27.D0*GAMMA_PBE*DS*PHI3*W1) &
     &       *(20.D0*RS*(GCZ*DECLDA0-DELSDA)-45.D0*R2ETA*(GCZ*ECLDA0-ELSDA)) &
     &       *P*EXP(-P**2/R2DP2**4)
         GY=1.D0/MAX(1.D0+4.D0*(Y-DELY),1.D-12)**0.25D0
         H1=GAMMA_PBE*PHI3*LOG(1.D0+W1*(1.D0-GY))
         EC1=ELSDA+H1
       END IF
!
!      == INTERPOLATION ================================================
       IF(TR2SCAN) THEN
!        == ALPHA-BAR=(TAU-TAU_W)/(TAU_UNIF*DS+ETA*TAU_W) ==============
         ALPHA=ALPHA0/(DS+R2ETA*5.D0/3.D0*S*S)
         CALL R2SCAN_FALPHA(ALPHA,R2CC,CFC2,CFDC1,IEF)
       ELSE
         ALPHA=ALPHA0/DS
         CALL SCAN_C_FALPHA(ALPHA,IEF)
       END IF
       EPSC=EC1+IEF*(EC0-EC1)
       RETURN
       END SUBROUTINE SCAN_C_EPSC_SPIN
!
!      .................................................................
       SUBROUTINE SCAN_C_EDENS(RHOUP,RHODN,G2,TAU,EC)
!      *****************************************************************
!      **  SCAN CORRELATION ENERGY PER VOLUME N*EPSC AT ONE POINT     **
!      **  FROM THE SPIN DENSITIES, THE SQUARED GRADIENT OF THE TOTAL **
!      **  DENSITY AND THE TOTAL KINETIC-ENERGY DENSITY.              **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOUP,RHODN
       REAL(8),INTENT(IN) :: G2      ! |GRAD(RHOUP+RHODN)|**2
       REAL(8),INTENT(IN) :: TAU     ! TAUUP+TAUDN
       REAL(8),INTENT(OUT):: EC
       REAL(8),PARAMETER  :: RHOMIN=1.D-10
       REAL(8)            :: PI,N,ZETA,KF,RS,S,T2,ALPHA0
!      *****************************************************************
       EC=0.D0
       N=MAX(RHOUP,0.D0)+MAX(RHODN,0.D0)
       IF(N.LT.RHOMIN) RETURN
       PI=4.D0*ATAN(1.D0)
       ZETA=(MAX(RHOUP,0.D0)-MAX(RHODN,0.D0))/N
       KF=(3.D0*PI**2*N)**(1.D0/3.D0)
       RS=(3.D0/(4.D0*PI*N))**(1.D0/3.D0)
       S=SQRT(MAX(G2,0.D0))/(2.D0*KF*N)
       T2=MAX(G2,0.D0)*PI/(16.D0*KF*N**2)   ! |GRAD N|**2/(2 KS N)**2
       ALPHA0=MAX((TAU-G2/(8.D0*N))/(0.3D0*KF**2*N),0.D0)
       CALL SCAN_C_EPSC_SPIN(RS,ZETA,S,T2,ALPHA0,EC)
       EC=N*EC
       RETURN
       END SUBROUTINE SCAN_C_EDENS
!
END MODULE SCAN_C_MODULE
!
!     ..................................................................
      SUBROUTINE SCAN$SETR2SCAN(TVAL)
!     ******************************************************************
!     **  SELECTS THE R2SCAN FORMS IN SCAN_X_MODULE AND SCAN_C_MODULE **
!     **  (CALLED BY DFT_INITIALIZE FOR EVERY CHANGE OF THE DFT TYPE) **
!     ******************************************************************
      USE SCAN_X_MODULE, ONLY: TR2SCAN
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TVAL
!     ******************************************************************
      TR2SCAN=TVAL
      RETURN
      END
!
!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  SUBROUTINE AUGMENTATION_XC_SCANL                                  **
!**                                                                    **
!**  SCAN-L EXCHANGE (+CORRELATION, IF DFT$SETL4('XCONLY',...) HAS NOT **
!**  DISABLED IT) FOR THE ONE-CENTER (AUGMENTATION-SPHERE) DENSITY,    **
!**  SPHERICAL (LM=1) PART ONLY, TOTAL-DENSITY (SPIN-UNPOLARIZED-      **
!**  TREATED) APPROXIMATION. FIRST WIRING PASS -- SCOPE CUTS BELOW ARE **
!**  DELIBERATE AND DOCUMENTED, NOT OVERSIGHTS:                        **
!**                                                                    **
!**   (1) NON-SPHERICAL (LM>1) DENSITY COMPONENTS ARE IGNORED. THE     **
!**       PBE0R PAPER'S OWN TARGET SYSTEMS (TRANSITION-METAL OXIDES)   **
!**       ARE EXACTLY WHERE THIS MATTERS MOST -- SEE                   **
!**       SCAN0r notes'S DISCUSSION OF VASP'S LASPH     **
!**       FLAG FOR THE SAME CAVEAT IN THE ANALOGOUS TRUE-SCAN CASE.    **
!**   (2) NO EXPLICIT SPIN-POLARIZATION: USES THE TOTAL DENSITY ONLY   **
!**       (RHOIN(:,1,1)), NOT THE EXACT SPIN-SCALING RELATION.         **
!**   (3) THE POTENTIAL (VXC) ONLY INCLUDES THE LOCAL                  **
!**       (D(EXDENS)/D(RHO) AT FIXED GRAD-RHO, LAPLACIAN-RHO) PIECE,   **
!**       BUILT BY FINITE DIFFERENCE. THE GRADIENT- AND LAPLACIAN-     **
!**       DEPENDENT CONTRIBUTIONS TO THE VARIATIONALLY EXACT POTENTIAL **
!**       (WHICH FOR THE LAPLACIAN TERM NEEDS A SECOND INTEGRATION BY  **
!**       PARTS, I.E. NABLA**2[D(EXDENS)/D(LAPRHO)] -- SEE             **
!**       SCAN0r notes SECTION 2, ITEM 3) ARE NOT INCLUDED. THIS   **
!**       MEANS THE ENERGY IS THE INTENDED SCAN-L VALUE, BUT SCF       **
!**       CALCULATIONS USING TYPE=20 WILL NOT EXACTLY SELF-            **
!**       CONSISTENTLY MINIMIZE IT.                                    **
!**   (4) CORRELATION IS GATED BY DFT$GETL4('CORRELATION',...), THE    **
!**       SAME TCORRELATION FLAG (SET VIA DFT$SETL4('XCONLY',...))     **
!**       THAT GATES EVERY OTHER CORRELATION MODULE IN DFT/DFT2/DFT3.  **
!**       THIS MATTERS FOR THE ACTUAL SCAN0R USE CASE: PBE0R'S OWN     **
!**       ONSITE DOUBLE-COUNTING CORRECTION (paw_ldaplusu.f90,         **
!**       paw_lmto.f90) ALWAYS CALLS AUGMENTATION_XC WITH XCONLY=T,    **
!**       I.E. EXCHANGE-ONLY -- SEE SCAN0r notes SECTION 1.        **
!**                                                                    **
!***********************************************************************
!      .................................................................
       SUBROUTINE AUGMENTATION_XC_SCANL(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
!      *****************************************************************
!      **  ONE-CENTER SCAN-L (SPHERICAL PART). SINCE 2026-10-02 SPIN   **
!      **  RESOLVED: EXACT SPIN SCALING FOR EXCHANGE AND SPIN-         **
!      **  POLARIZED CORRELATION (SCANL_EDENS_SPIN); BEFORE, THE       **
!      **  SPIN-SUMMED DENSITY WAS TREATED AS UNPOLARIZED.             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: GID
       INTEGER(4),INTENT(IN) :: NR
       INTEGER(4),INTENT(IN) :: LMRX
       INTEGER(4),INTENT(IN) :: NDIMD
       REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
       REAL(8)   ,INTENT(OUT):: EXC
       REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
       REAL(8)                :: PI,Y0
       REAL(8)                :: R(NR)
       REAL(8)                :: RHOUP(NR),RHODN(NR) ! ACTUAL (Y0-SCALED) SPHERICAL DENSITIES
       REAL(8)                :: FXC(NR),VUP(NR),VDN(NR)
!      *****************************************************************
       IF(NDIMD.EQ.4) THEN
         CALL ERROR$MSG('SCAN-L ONE-CENTER TERM NOT IMPLEMENTED FOR NONCOLLINEAR SPIN')
         CALL ERROR$STOP('AUGMENTATION_XC_SCANL')
       END IF
       PI=4.D0*ATAN(1.D0)
       Y0=1.D0/SQRT(4.D0*PI)
       EXC=0.D0
       VXC(:,:,:)=0.D0
       CALL RADIAL$R(GID,NR,R)
       IF(NDIMD.EQ.2) THEN
         RHOUP(:)=0.5D0*(RHOIN(:,1,1)+RHOIN(:,1,2))*Y0
         RHODN(:)=0.5D0*(RHOIN(:,1,1)-RHOIN(:,1,2))*Y0
       ELSE
         RHOUP(:)=0.5D0*RHOIN(:,1,1)*Y0
         RHODN(:)=RHOUP(:)
       END IF
       CALL AUGMENTATION_SCANL_SPIN(GID,NR,RHOUP,RHODN,FXC,VUP,VDN)
       CALL RADIAL$INTEGRAL(GID,NR,FXC(:)*R(:)**2,EXC)
       VXC(:,1,1)=0.5D0*(VUP(:)+VDN(:))
       IF(NDIMD.EQ.2) VXC(:,1,2)=0.5D0*(VUP(:)-VDN(:))
       RETURN
       END SUBROUTINE AUGMENTATION_XC_SCANL
!
!      .................................................................
       SUBROUTINE AUGMENTATION_SCANL_SPIN(GID,NR,RHOUP,RHODN,FXC,VUP,VDN)
!      *****************************************************************
!      **  SPHERICAL SPIN-RESOLVED SCAN-L ENERGY DENSITY FXC=4PI*E(R)  **
!      **  AND POTENTIALS (AS LM=1 COEFFICIENTS) FOR EACH SPIN:        **
!      **    V_S = E_N(S) - (1/R**2)(R**2 E_G(S))' + LAP(E_L(S))       **
!      **  WITH E_G(S),E_L(S) THE DERIVATIVES WITH RESPECT TO THE      **
!      **  RADIAL DERIVATIVE AND THE LAPLACIAN OF THE SPIN DENSITY,    **
!      **  ALL BY CENTRAL FINITE DIFFERENCES OF SCANL_EDENS_SPIN.      **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: GID
       INTEGER(4),INTENT(IN) :: NR
       REAL(8)   ,INTENT(IN) :: RHOUP(NR),RHODN(NR)
       REAL(8)   ,INTENT(OUT):: FXC(NR)
       REAL(8)   ,INTENT(OUT):: VUP(NR),VDN(NR)
       REAL(8)   ,PARAMETER  :: RHOMIN=1.D-8
       REAL(8)                :: PI,Y0,FOURPI
       REAL(8)                :: R(NR)
       REAL(8)                :: G(NR,2),L(NR,2),RHO1(NR,1,1),LAP1(NR,1,1)
       REAL(8)                :: EN(NR,2),EG(NR,2),EL(NR,2),V(NR,2)
       REAL(8)                :: AUX(NR),AUX2(NR)
       REAL(8)                :: X(6),XP(6),XM(6),E0,EP,EM,D
       LOGICAL(4)             :: TCORR
       INTEGER(4)             :: IR,IS,K
!      *****************************************************************
       PI=4.D0*ATAN(1.D0)
       FOURPI=4.D0*PI
       Y0=1.D0/SQRT(FOURPI)
       CALL DFT$GETL4('CORRELATION',TCORR)
       CALL RADIAL$R(GID,NR,R)
       CALL RADIAL$DERIVE(GID,NR,RHOUP,G(:,1))
       CALL RADIAL$DERIVE(GID,NR,RHODN,G(:,2))
       RHO1(:,1,1)=RHOUP(:)/Y0
       CALL AUGMENTATION_LAPRHO(GID,NR,1,1,RHO1,LAP1)
       L(:,1)=LAP1(:,1,1)*Y0
       RHO1(:,1,1)=RHODN(:)/Y0
       CALL AUGMENTATION_LAPRHO(GID,NR,1,1,RHO1,LAP1)
       L(:,2)=LAP1(:,1,1)*Y0
!
       FXC(:)=0.D0
       EN(:,:)=0.D0
       EG(:,:)=0.D0
       EL(:,:)=0.D0
       DO IR=1,NR
         IF(RHOUP(IR)+RHODN(IR).LT.RHOMIN) CYCLE
!        == X=(NUP,NDN,GUP,GDN,LUP,LDN) ===============================
         X=(/MAX(RHOUP(IR),0.D0),MAX(RHODN(IR),0.D0),G(IR,1),G(IR,2),L(IR,1),L(IR,2)/)
         CALL SCANL_EDENS_SPIN(X(1),X(2),X(3),X(4),X(5),X(6),TCORR,E0)
         FXC(IR)=FOURPI*E0
         DO K=1,6
           IF(K.LE.2) THEN
             D=MAX(1.D-6*X(K),1.D-12)
           ELSE IF(K.LE.4) THEN
             D=MAX(1.D-6*ABS(X(K)),1.D-10*(X(1)+X(2))**(4.D0/3.D0),1.D-30)
           ELSE
             D=MAX(1.D-6*ABS(X(K)),1.D-10*(X(1)+X(2))**(5.D0/3.D0),1.D-30)
           END IF
           XP=X; XM=X
           XP(K)=X(K)+D
           XM(K)=X(K)-D
           IF(K.LE.2) XM(K)=MAX(XM(K),0.D0)
           CALL SCANL_EDENS_SPIN(XP(1),XP(2),XP(3),XP(4),XP(5),XP(6),TCORR,EP)
           CALL SCANL_EDENS_SPIN(XM(1),XM(2),XM(3),XM(4),XM(5),XM(6),TCORR,EM)
           IS=2-MOD(K,2)
           IF(K.LE.2) THEN
             EN(IR,IS)=(EP-EM)/(XP(K)-XM(K))
           ELSE IF(K.LE.4) THEN
             EG(IR,IS)=(EP-EM)/(XP(K)-XM(K))
           ELSE
             EL(IR,IS)=(EP-EM)/(XP(K)-XM(K))
           END IF
         ENDDO
       ENDDO
!
!      == V = E_N - DIV(E_G) + LAP(E_L), AS LM=1 COEFFICIENTS ============
       DO IS=1,2
         AUX(:)=EG(:,IS)*R(:)**2
         CALL RADIAL$DERIVE(GID,NR,AUX,AUX2)
         V(:,IS)=EN(:,IS)-AUX2(:)/R(:)**2
         CALL RADIAL$DERIVE(GID,NR,EL(:,IS),AUX)
         AUX(:)=AUX(:)*R(:)**2
         CALL RADIAL$DERIVE(GID,NR,AUX,AUX2)
         V(:,IS)=(V(:,IS)+AUX2(:)/R(:)**2)*FOURPI*Y0
       ENDDO
       VUP(:)=V(:,1)
       VDN(:)=V(:,2)
       RETURN
       END SUBROUTINE AUGMENTATION_SCANL_SPIN
!
!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  SUBROUTINE AUGMENTATION_XC_TSCANDELOC                             **
!**                                                                    **
!**  REAL (ORBITAL-TAU-BASED, NOT DEORBITALIZED) SCAN EXCHANGE         **
!**  (+CORRELATION, GATED THE SAME WAY AS EVERY OTHER CORRELATION      **
!**  MODULE VIA DFT$GETL4('CORRELATION',...)) FOR THE ONE-CENTER       **
!**  (AUGMENTATION-SPHERE) DENSITY -- THE STANDARD PAW "AE MINUS PS"   **
!**  RECONSTRUCTION COUNTERPART TO TSCANDELOC'S OWN INTERSTITIAL/      **
!**  PLANE-WAVE TERM (POTENTIAL_VOFRHO, paw_potential.f90) AND         **
!**  HAMILTONIAN ACTION (WAVES_TAUHPSI, paw_scankinden.f90).           **
!**                                                                    **
!**  BEFORE THIS ROUTINE EXISTED, AUGMENTATION_XC HAD NO TSCANDELOC    **
!**  BRANCH AT ALL (UNLIKE TSCANL, WHICH HAS AUGMENTATION_XC_SCANL) -- **
!**  A GENUINE ARCHITECTURAL GAP, NOT A DELIBERATE SCOPE CUT: THE      **
!**  STANDARD ONE-CENTER RECONSTRUCTION, CALLED FOR EVERY ATOM/STEP    **
!**  REGARDLESS OF HYBRID/NTBO/DCTYPE SETTINGS, SILENTLY FELL THROUGH  **
!**  TO PLAIN PBE EVEN WHEN THE INTERSTITIAL TERM WAS REAL SCAN. SEE   **
!**  AUGMENTATION_XC'S OWN TSCANDELOC DISPATCH COMMENT FOR HOW THIS    **
!**  WAS FOUND (COMPARING WHICH FILES REFERENCE TSCANDELOC AT ALL).    **
!**                                                                    **
!**  NON-SPHERICAL (LM>1) ONE-CENTER TREATMENT (2026-08-28, THE        **
!**  LASPH-EQUIVALENT EXTENSION): WHEN LMRX>1, DISPATCHES TO           **
!**  AUGMENTATION_SCAN_NSPH (paw_scannsph.f90, DIRECT ANGULAR          **
!**  QUADRATURE, NOW INCLUDING CORRELATION), WITH CUT(:)=1 EVERYWHERE  **
!**  (NO LOCAL/TOTAL SPLIT HERE, UNLIKE THE LMTO ONSITE DOUBLE-        **
!**  COUNTING SCHEME -- THE WHOLE SPHERE IS TREATED UNIFORMLY, SAME AS **
!**  THE LM=1 FALLBACK BELOW). MIRRORS LMTO_RADXCTAU_NSPH_WITHCUT'S    **
!**  OWN SPIN DISPATCH (paw_lmto.f90) EXACTLY: NDIMD=1 DIRECT CALL;    **
!**  NDIMD=2 VIA THE SAME EXACT SPIN-SCALING RELATION USED BELOW FOR   **
!**  THE SPHERICAL PATH, NOW APPLIED TO THE RAW LM-COEFFICIENT ARRAYS  **
!**  DIRECTLY (NO Y0 CONVERSION NEEDED -- AUGMENTATION_SCAN_NSPH TAKES **
!**  RAW LM COEFFICIENTS, NOT A Y0-PREMULTIPLIED SCALAR DENSITY).      **
!**  REMAINING FIRST-PASS SCOPE CUTS, DELIBERATE AND DOCUMENTED:       **
!**   (1) (REMOVED -- SEE ABOVE.)                                      **
!**   (2) NDIMD=1 (SPIN-RESTRICTED) OR NDIMD=2 (COLLINEAR SPIN) ONLY,  **
!**       USING THE EXACT SPIN-SCALING RELATION FOR EXCHANGE,          **
!**       E_X[N_UP,N_DN]=(1/2)*{E_X_UNPOL[2*N_UP]+E_X_UNPOL[2*N_DN]},  **
!**       SAME PATTERN AS LMTO_RADXCTAU_WITHCUT (paw_lmto.f90) --      **
!**       "2*N_UP"=(RHOT+RHOZ) WITH TAU=(TAUT+TAUZ) (SAME REASONING,   **
!**       TAU IS ADDITIVE UNDER THE SAME "DOUBLE THE SYSTEM"           **
!**       CONSTRUCTION AS RHO), CALLING THE *SAME*, UNCHANGED,         **
!**       AUGMENTATION_SCAN_POINTWISE TWICE AND AVERAGING/HALF-        **
!**       DIFFERENCING INTO THE (T,Z) CHANNELS. NDIMD=4 (NONCOLLINEAR) **
!**       IS NOT SUPPORTED -- ERROR$STOP, MATCHING HOW TSCAN (DCTYPE=  **
!**       21) ITSELF REFUSES NDIMD=4 IN LMTO_SIMPLEDC_NEW_NEW.         **
!**   (3) THE VALUE-TYPE POTENTIAL (VXC) ONLY INCLUDES THE LOCAL       **
!**       (D(EXCDENS)/D(RHOT) AT FIXED TAU) PIECE, MATCHING SCANL'S    **
!**       OWN SCOPE. THE TAU-VARIATIONAL CONTRIBUTION IS NOW INCLUDED  **
!**       SEPARATELY: VTAU (D(EXCDENS)/D(TAU) AT FIXED RHOT) IS        **
!**       STASHED INTO AUGMENTATION_MODULE'S XCVTAU FOR THE CALLER     **
!**       (AUGMENTATION$SPHERE) TO FOLD BACK INTO DATH VIA             **
!**       AUGMENTATION_TAU_EXPECT (paw_scantau.f90), THE SAME WAY      **
!**       DCTYPE=21'S OWN LMTOAUGMENTATION$ADDC8A('DH',...) CALL DOES  **
!**       IT FOR THE HYBRID DOUBLE-COUNTING TERM.                      **
!**   (4) TAU IS SUPPLIED BY THE CALLER (AUGMENTATION$SPHERE) VIA      **
!**       AUGMENTATION_XC$SETTAU IMMEDIATELY BEFORE EACH OF ITS THREE  **
!**       AUGMENTATION_XC CALLS (CORE-ONLY, AE, PS) -- SEE             **
!**       AUGMENTATION_MODULE'S XCTAU DECLARATION FOR THE CONVENTION.  **
!***********************************************************************
!      .................................................................
       SUBROUTINE AUGMENTATION_XC_TSCANDELOC(GID,NR,LMRX,NDIMD,RHOIN,TAUIN,EXC,VXC)
       USE AUGMENTATION_MODULE, ONLY: XCVTAU
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: GID
       INTEGER(4),INTENT(IN) :: NR
       INTEGER(4),INTENT(IN) :: LMRX
       INTEGER(4),INTENT(IN) :: NDIMD
       REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
       REAL(8)   ,INTENT(IN) :: TAUIN(NR,LMRX,NDIMD) ! RAW (NOT Y0-SCALED), SAME LM RANGE AS RHOIN
       REAL(8)   ,INTENT(OUT):: EXC
       REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
       LOGICAL(4)              :: TCORR
       REAL(8)                :: PI,Y0
       REAL(8)                :: R(NR)
       REAL(8)                :: RHOT_UP(NR),RHOT_DN(NR)  ! ACTUAL (Y0-SCALED), =2*RHO(UP/DN)
       REAL(8)                :: TAUT_UP(NR),TAUT_DN(NR)  ! ACTUAL (Y0-SCALED), =2*TAU(UP/DN)
       REAL(8)                :: FXC(NR),FXC_UP(NR),FXC_DN(NR)
       REAL(8)                :: VXCLOC_UP(NR),VXCLOC_DN(NR)
       REAL(8)                :: VTAU_UP(NR),VTAU_DN(NR)
       REAL(8)                :: WORK(NR)
       INTEGER(4)              :: LX
       INTEGER(4),PARAMETER    :: NTH=14,NPH=24  ! SAME DEFAULT AS LMTO_RADXCTAU_NSPH_WITHCUT -- CONFIRMED
!      ADEQUATE VS. A 20,32 (640-POINT) DIAGNOSTIC RUN ON NIO (2026-08-29): GAP
!      CHANGED BY <0.01 EV, SEE scan0r ANOMALY STATUS MEMORY.
       REAL(8)   ,ALLOCATABLE  :: CUT(:)
       REAL(8)   ,ALLOCATABLE  :: RHO_UP(:,:),RHO_DN(:,:),TAU_UP(:,:),TAU_DN(:,:)
       REAL(8)   ,ALLOCATABLE  :: VXC_UP(:,:),VXC_DN(:,:),VTAU_UP_A(:,:),VTAU_DN_A(:,:)
       REAL(8)   ,ALLOCATABLE  :: VCUT(:),VCUT_UP(:),VCUT_DN(:)
       REAL(8)                 :: EXC_UP,EXC_DN
!      *****************************************************************
       PI=4.D0*ATAN(1.D0)
       Y0=1.D0/SQRT(4.D0*PI)
       EXC=0.D0
       VXC(:,:,:)=0.D0
       CALL RADIAL$R(GID,NR,R)
       CALL DFT$GETL4('CORRELATION',TCORR)
       IF(ALLOCATED(XCVTAU)) DEALLOCATE(XCVTAU)
!
!      =================================================================
!      ==  NON-SPHERICAL (LM>1) PATH -- AUGMENTATION_SCAN_NSPH VIA     ==
!      ==  ANGULAR QUADRATURE, MIRRORING LMTO_RADXCTAU_NSPH_WITHCUT'S  ==
!      ==  SPIN DISPATCH EXACTLY (paw_lmto.f90).                       ==
!      =================================================================
       IF(LMRX.GT.1) THEN
         LX=NINT(SQRT(DBLE(LMRX)))-1
         ALLOCATE(CUT(NR)); CUT(:)=1.D0   ! NO LOCAL/TOTAL SPLIT HERE
         ALLOCATE(VCUT(NR))
         ALLOCATE(RHO_UP(NR,LMRX),TAU_UP(NR,LMRX))
         ALLOCATE(VXC_UP(NR,LMRX),VTAU_UP_A(NR,LMRX))
         IF(NDIMD.EQ.1) THEN
           RHO_UP(:,:)=RHOIN(:,:,1)
           TAU_UP(:,:)=TAUIN(:,:,1)
           CALL AUGMENTATION_SCAN_NSPH(GID,NR,LX,NTH,NPH,RHO_UP,TAU_UP,CUT &
     &                                 ,.FALSE.,EXC,VXC(:,:,1),VTAU_UP_A,VCUT)
           ALLOCATE(XCVTAU(NR,LMRX,1))
           XCVTAU(:,:,1)=VTAU_UP_A(:,:)
         ELSE IF(NDIMD.EQ.2) THEN
           ALLOCATE(RHO_DN(NR,LMRX),TAU_DN(NR,LMRX))
           ALLOCATE(VXC_DN(NR,LMRX),VTAU_DN_A(NR,LMRX))
           ALLOCATE(VCUT_UP(NR),VCUT_DN(NR))
           RHO_UP(:,:)=RHOIN(:,:,1)+RHOIN(:,:,2)   ! =2*RHO(SPIN UP), PER LM
           RHO_DN(:,:)=RHOIN(:,:,1)-RHOIN(:,:,2)   ! =2*RHO(SPIN DOWN), PER LM
           TAU_UP(:,:)=TAUIN(:,:,1)+TAUIN(:,:,2)
           TAU_DN(:,:)=TAUIN(:,:,1)-TAUIN(:,:,2)
           CALL AUGMENTATION_SCAN_NSPH(GID,NR,LX,NTH,NPH,RHO_UP,TAU_UP,CUT &
     &                                 ,.FALSE.,EXC_UP,VXC_UP,VTAU_UP_A,VCUT_UP)
           CALL AUGMENTATION_SCAN_NSPH(GID,NR,LX,NTH,NPH,RHO_DN,TAU_DN,CUT &
     &                                 ,.FALSE.,EXC_DN,VXC_DN,VTAU_DN_A,VCUT_DN)
           EXC=0.5D0*(EXC_UP+EXC_DN)
           VXC(:,:,1)=0.5D0*(VXC_UP(:,:)+VXC_DN(:,:))
           VXC(:,:,2)=0.5D0*(VXC_UP(:,:)-VXC_DN(:,:))
           ALLOCATE(XCVTAU(NR,LMRX,2))
           XCVTAU(:,:,1)=0.5D0*(VTAU_UP_A(:,:)+VTAU_DN_A(:,:))
           XCVTAU(:,:,2)=0.5D0*(VTAU_UP_A(:,:)-VTAU_DN_A(:,:))
           DEALLOCATE(RHO_DN,TAU_DN,VXC_DN,VTAU_DN_A,VCUT_UP,VCUT_DN)
         ELSE
           CALL ERROR$MSG('TSCANDELOC ONE-CENTER RECONSTRUCTION SUPPORTS')
           CALL ERROR$MSG('NDIMD=1 (SPIN-RESTRICTED) OR NDIMD=2 (COLLINEAR')
           CALL ERROR$MSG('SPIN) ONLY -- NOT NDIMD=4 (NONCOLLINEAR)')
           CALL ERROR$I4VAL('NDIMD',NDIMD)
           CALL ERROR$STOP('AUGMENTATION_XC_TSCANDELOC')
         END IF
         IF(TCORR) CALL TSCANDELOC_ADDC(GID,NR,LMRX,NDIMD,RHOIN,TAUIN,EXC,VXC)
         DEALLOCATE(CUT,VCUT,RHO_UP,TAU_UP,VXC_UP,VTAU_UP_A)
         RETURN
       END IF
!
!      =================================================================
!      ==  SPHERICAL-ONLY (LM=1) FALLBACK -- THE ORIGINAL FIRST-PASS   ==
!      ==  IMPLEMENTATION, VIA AUGMENTATION_SCAN_POINTWISE.            ==
!      =================================================================
       IF(NDIMD.EQ.1) THEN
         RHOT_UP(:)=RHOIN(:,1,1)*Y0
         TAUT_UP(:)=TAUIN(:,1,1)*Y0
         CALL AUGMENTATION_SCAN_POINTWISE(GID,NR,RHOT_UP,TAUT_UP,.FALSE. &
     &                                    ,FXC,VXCLOC_UP,VTAU_UP)
         VXC(:,1,1)=VXCLOC_UP(:)
         ALLOCATE(XCVTAU(NR,1,1))
         XCVTAU(:,1,1)=VTAU_UP(:)
       ELSE IF(NDIMD.EQ.2) THEN
         RHOT_UP(:)=(RHOIN(:,1,1)+RHOIN(:,1,2))*Y0   ! =2*RHO(SPIN UP)
         RHOT_DN(:)=(RHOIN(:,1,1)-RHOIN(:,1,2))*Y0   ! =2*RHO(SPIN DOWN)
         TAUT_UP(:)=(TAUIN(:,1,1)+TAUIN(:,1,2))*Y0
         TAUT_DN(:)=(TAUIN(:,1,1)-TAUIN(:,1,2))*Y0
         CALL AUGMENTATION_SCAN_POINTWISE(GID,NR,RHOT_UP,TAUT_UP,.FALSE. &
     &                                    ,FXC_UP,VXCLOC_UP,VTAU_UP)
         CALL AUGMENTATION_SCAN_POINTWISE(GID,NR,RHOT_DN,TAUT_DN,.FALSE. &
     &                                    ,FXC_DN,VXCLOC_DN,VTAU_DN)
         FXC(:)=0.5D0*(FXC_UP(:)+FXC_DN(:))
         VXC(:,1,1)=0.5D0*(VXCLOC_UP(:)+VXCLOC_DN(:))
         VXC(:,1,2)=0.5D0*(VXCLOC_UP(:)-VXCLOC_DN(:))
         ALLOCATE(XCVTAU(NR,1,2))
         XCVTAU(:,1,1)=0.5D0*(VTAU_UP(:)+VTAU_DN(:))
         XCVTAU(:,1,2)=0.5D0*(VTAU_UP(:)-VTAU_DN(:))
       ELSE
         CALL ERROR$MSG('TSCANDELOC ONE-CENTER RECONSTRUCTION SUPPORTS')
         CALL ERROR$MSG('NDIMD=1 (SPIN-RESTRICTED) OR NDIMD=2 (COLLINEAR')
         CALL ERROR$MSG('SPIN) ONLY -- NOT NDIMD=4 (NONCOLLINEAR)')
         CALL ERROR$I4VAL('NDIMD',NDIMD)
         CALL ERROR$STOP('AUGMENTATION_XC_TSCANDELOC')
       END IF
       WORK(:)=FXC(:)*R(:)**2
       CALL RADIAL$INTEGRAL(GID,NR,WORK,EXC)
       IF(TCORR) CALL TSCANDELOC_ADDC(GID,NR,LMRX,NDIMD,RHOIN,TAUIN,EXC,VXC)
       RETURN
       END SUBROUTINE AUGMENTATION_XC_TSCANDELOC
!
!      .................................................................
       SUBROUTINE TSCANDELOC_ADDC(GID,NR,LMRX,NDIMD,RHOIN,TAUIN,EXC,VXC)
!      *****************************************************************
!      **  ADDS THE SPIN-POLARIZED SCAN CORRELATION (SCANC_ONECENTER)  **
!      **  TO THE SPIN-SCALED ONE-CENTER EXCHANGE OF                  **
!      **  AUGMENTATION_XC_TSCANDELOC: ENERGY, POTENTIAL AND XCVTAU.  **
!      *****************************************************************
       USE AUGMENTATION_MODULE, ONLY: XCVTAU
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: GID,NR,LMRX,NDIMD
       REAL(8)   ,INTENT(IN)   :: RHOIN(NR,LMRX,NDIMD)
       REAL(8)   ,INTENT(IN)   :: TAUIN(NR,LMRX,NDIMD)
       REAL(8)   ,INTENT(INOUT):: EXC
       REAL(8)   ,INTENT(INOUT):: VXC(NR,LMRX,NDIMD)
       REAL(8)                 :: CUT(NR),VCUT(NR),EC
       REAL(8)                 :: VC(NR,LMRX,NDIMD),VTC(NR,LMRX,NDIMD)
!      *****************************************************************
       CUT(:)=1.D0
       CALL SCANC_ONECENTER(GID,NR,LMRX,NDIMD,RHOIN,TAUIN,CUT,EC,VC,VTC,VCUT)
       EXC=EXC+EC
       VXC(:,:,:)=VXC(:,:,:)+VC(:,:,:)
       XCVTAU(:,:,1)=XCVTAU(:,:,1)+VTC(:,:,1)
       RETURN
       END SUBROUTINE TSCANDELOC_ADDC
!
!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**  SUBROUTINE AUGMENTATION_SCANL_POINTWISE                          **
!**                                                                    **
!**  SCAN-L (DEORBITALIZED SCAN) EXCHANGE(+CORRELATION) OF THE ONE-    **
!**  CENTER SPHERICAL, SPIN-SUMMED DENSITY, SHARED BY                  **
!**  AUGMENTATION_XC_SCANL AND LMTO_RADXC_WITHCUT (SCAN0R DOUBLE       **
!**  COUNTING).                                                        **
!**                                                                    **
!**  THE ENERGY IS E=INT W(R)*FXC(R) R**2 DR WITH A CALLER-SUPPLIED    **
!**  WEIGHT W (1, OR THE CUT FUNCTION OF THE DOUBLE COUNTING), AND     **
!**  FXC=FOURPI*RHO*EPS_XC(RHO,|GRAD RHO|**2,LAP RHO). VXC IS ITS FULL **
!**  FUNCTIONAL DERIVATIVE WITH RESPECT TO THE LM=1 DENSITY COMPONENT: **
!**     V = W E_RHO - (1/R**2)(R**2 W 2RHO' E_G2)'                     **
!**           + (1/R**2)(R**2 (W E_LAP)')'                             **
!**  (E=RHO*EPS, G2=|GRAD RHO|**2), TIMES FOURPI*Y0. THE GRADIENT AND  **
!**  LAPLACIAN TERMS WERE MISSING BEFORE 2026-09-27, SO THE SCAN0R     **
!**  POTENTIAL WAS NOT THE DERIVATIVE OF ITS ENERGY (SELF-CONSISTENT   **
!**  SOLUTIONS WERE NOT STATIONARY POINTS OF THE ENERGY).              **
!**                                                                    **
!**  INPUT RHOT(NR) IS THE ACTUAL (Y0-SCALED) SPHERICAL DENSITY.       **
!**  OUTPUT FXC IS UNWEIGHTED (W NOT APPLIED); VXC INCLUDES W.         **
!***********************************************************************
       SUBROUTINE AUGMENTATION_SCANL_POINTWISE(GID,NR,RHOT,W,FXC,VXC)
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: GID
       INTEGER(4),INTENT(IN) :: NR
       REAL(8)   ,INTENT(IN) :: RHOT(NR)
       REAL(8)   ,INTENT(IN) :: W(NR)
       REAL(8)   ,INTENT(OUT):: FXC(NR)
       REAL(8)   ,INTENT(OUT):: VXC(NR)
       REAL(8)   ,PARAMETER  :: RHOMIN=1.D-8
       REAL(8)                :: PI,Y0,FOURPI
       REAL(8)                :: R(NR),GRHOT(NR),LAPRHOT(NR)
       REAL(8)                :: RHO1(NR,1,1),LAPRHO1(NR,1,1)
       REAL(8)                :: ERHO(NR),EG2(NR),ELAP(NR)
       REAL(8)                :: A(NR),B(NR),AUX(NR),AUX2(NR)
       REAL(8)                :: RHO,G2,XLAP,E0,EP,EM,D
       LOGICAL(4)             :: TCORR
       INTEGER(4)             :: IR
!      *****************************************************************
       PI=4.D0*ATAN(1.D0)
       FOURPI=4.D0*PI
       Y0=1.D0/SQRT(FOURPI)
       CALL DFT$GETL4('CORRELATION',TCORR)
       CALL RADIAL$R(GID,NR,R)
       CALL RADIAL$DERIVE(GID,NR,RHOT,GRHOT)
       RHO1(:,1,1)=RHOT(:)/Y0
       CALL AUGMENTATION_LAPRHO(GID,NR,1,1,RHO1,LAPRHO1)
       LAPRHOT(:)=LAPRHO1(:,1,1)*Y0
!
!      == ENERGY DENSITY AND ITS PARTIAL DERIVATIVES (CENTRAL DIFFERENCES)
       FXC=0.D0
       ERHO=0.D0
       EG2=0.D0
       ELAP=0.D0
       DO IR=1,NR
         RHO=RHOT(IR)
         IF(RHO.LT.RHOMIN) CYCLE
         G2=GRHOT(IR)**2
         XLAP=LAPRHOT(IR)
         CALL AUGMENTATION_SCANL_EDENS(RHO,G2,XLAP,TCORR,E0)
         FXC(IR)=FOURPI*E0
         D=MAX(1.D-6*RHO,1.D-12)
         CALL AUGMENTATION_SCANL_EDENS(RHO+D,G2,XLAP,TCORR,EP)
         CALL AUGMENTATION_SCANL_EDENS(MAX(RHO-D,RHOMIN),G2,XLAP,TCORR,EM)
         ERHO(IR)=(EP-EM)/(RHO+D-MAX(RHO-D,RHOMIN))
         D=MAX(1.D-6*G2,1.D-10*RHO**(8.D0/3.D0),1.D-30)
         CALL AUGMENTATION_SCANL_EDENS(RHO,G2+D,XLAP,TCORR,EP)
         CALL AUGMENTATION_SCANL_EDENS(RHO,MAX(G2-D,0.D0),XLAP,TCORR,EM)
         EG2(IR)=(EP-EM)/(G2+D-MAX(G2-D,0.D0))
         D=MAX(1.D-6*ABS(XLAP),1.D-10*RHO**(5.D0/3.D0),1.D-30)
         CALL AUGMENTATION_SCANL_EDENS(RHO,G2,XLAP+D,TCORR,EP)
         CALL AUGMENTATION_SCANL_EDENS(RHO,G2,XLAP-D,TCORR,EM)
         ELAP(IR)=(EP-EM)/(2.D0*D)
       ENDDO
!
!      == V = W E_RHO - DIV(W DE/DGRAD) + LAP(W DE/DLAP) ==================
       A(:)=W(:)*2.D0*GRHOT(:)*EG2(:)*R(:)**2
       CALL RADIAL$DERIVE(GID,NR,A,AUX)
       VXC(:)=W(:)*ERHO(:)-AUX(:)/R(:)**2
       B(:)=W(:)*ELAP(:)
       CALL RADIAL$DERIVE(GID,NR,B,AUX)
       AUX(:)=AUX(:)*R(:)**2
       CALL RADIAL$DERIVE(GID,NR,AUX,AUX2)
       VXC(:)=VXC(:)+AUX2(:)/R(:)**2
       VXC(:)=VXC(:)*FOURPI*Y0
       RETURN
       END SUBROUTINE AUGMENTATION_SCANL_POINTWISE
!
!.........1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE AUGMENTATION_SCANL_EDENS(RHO,G2,XLAP,TCORR,EXCDENS)
!      *****************************************************************
!      **  SCAN-L ENERGY DENSITY RHO*(EPS_X+EPS_C) AT ONE POINT OF A   **
!      **  SPIN-UNPOLARIZED DENSITY WITH |GRAD RHO|**2=G2 AND          **
!      **  LAPLACIAN XLAP.                                             **
!      *****************************************************************
       USE SCAN_X_MODULE, ONLY: SCAN_X_FX
       USE SCAN_C_MODULE, ONLY: SCAN_C_EPSC
       USE DEORB_ALPHA_MODULE, ONLY: DEORB_SQ,DEORB_ALPHA
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: RHO,G2,XLAP
       LOGICAL(4),INTENT(IN) :: TCORR
       REAL(8)   ,INTENT(OUT):: EXCDENS
       REAL(8)   ,PARAMETER  :: EXFAC=-0.75D0*(3.D0/3.14159265358979D0)**(1.D0/3.D0)
       REAL(8)               :: PI,RSFAC,KFFAC,KSFAC,TFAC
       REAL(8)               :: S,Q,ALPHA,EPSX,EPSC,RS,T2
!      *****************************************************************
       PI=4.D0*ATAN(1.D0)
       RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0)
       KFFAC=(3.D0*PI**2)**(1.D0/3.D0)
       KSFAC=SQRT(4.D0*KFFAC/PI)
       TFAC=PI/(16.D0*KFFAC) ! T2 OF PBE/SCAN: KF=KFFAC*RHO**(1/3) (WAS 0.62x TOO SMALL)
       CALL DEORB_SQ(RHO,G2,XLAP,S,Q)
       CALL DEORB_ALPHA(S,Q,ALPHA)
       CALL SCAN_X_FX(S,ALPHA,EPSX)
       EPSX=EXFAC*RHO**(1.D0/3.D0)*EPSX
       EPSC=0.D0
       IF(TCORR) THEN
         RS=RSFAC/RHO**(1.D0/3.D0)
         T2=TFAC*G2*RHO**(-7.D0/3.D0)
         CALL SCAN_C_EPSC(RS,S,T2,ALPHA,EPSC)
       END IF
       EXCDENS=RHO*(EPSX+EPSC)
       RETURN
       END SUBROUTINE AUGMENTATION_SCANL_EDENS
!
!.........1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SCAN_XSIG_TAU(RHOSIG,G2SIG,TAUSIG,EX)
!      *****************************************************************
!      **  SCAN EXCHANGE ENERGY PER VOLUME OF ONE SPIN CHANNEL BY THE  **
!      **  EXACT SPIN-SCALING RELATION                                **
!      **    EX[NUP,NDN]=(EX[2NUP]+EX[2NDN])/2                         **
!      **  I.E. EX=RHOSIG*EPSX(2*RHOSIG,4*G2SIG,2*TAUSIG), WITH THE   **
!      **  ORBITAL KINETIC-ENERGY DENSITY TAUSIG OF THIS SPIN.        **
!      *****************************************************************
       USE SCAN_X_MODULE, ONLY: SCAN_X_FX
       USE DEORB_ALPHA_MODULE, ONLY: DEORB_SQ
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOSIG,G2SIG,TAUSIG
       REAL(8),INTENT(OUT):: EX
       REAL(8),PARAMETER  :: RHOMIN=1.D-10
       REAL(8)            :: PI,N,G2,TAU,S,Q,TAUUNIF,TAUW,ALPHA,FX
!      *****************************************************************
       EX=0.D0
       IF(RHOSIG.LT.RHOMIN) RETURN
       PI=4.D0*ATAN(1.D0)
       N=2.D0*RHOSIG
       G2=4.D0*MAX(G2SIG,0.D0)
       TAU=2.D0*TAUSIG
       CALL DEORB_SQ(N,G2,0.D0,S,Q)
       TAUUNIF=0.3D0*(3.D0*PI**2)**(2.D0/3.D0)*N**(5.D0/3.D0)
       TAUW=G2/(8.D0*N)
       ALPHA=MAX((TAU-TAUW)/MAX(TAUUNIF,1.D-16),0.D0)
       CALL SCAN_X_FX(S,ALPHA,FX)
       EX=RHOSIG*(-0.75D0)*(3.D0/PI)**(1.D0/3.D0)*N**(1.D0/3.D0)*FX
       RETURN
       END SUBROUTINE SCAN_XSIG_TAU
!
!.........1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SCANL_TAUSIG(RHOSIG,G2SIG,LAPSIG,TAUSIG)
!      *****************************************************************
!      **  DEORBITALIZED (SCAN-L) KINETIC-ENERGY DENSITY OF ONE SPIN  **
!      **  CHANNEL, USING THE SPIN SCALING OF THE KINETIC ENERGY      **
!      **    TAU[NUP,NDN]=(TAU[2NUP]+TAU[2NDN])/2                      **
!      **  AS IN LIBXC: TAU_PC(N)=TAU_W(N)+TAU_UNIF(N)*ALPHA_PC(S,Q). **
!      *****************************************************************
       USE DEORB_ALPHA_MODULE, ONLY: DEORB_SQ,DEORB_ALPHA
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOSIG,G2SIG,LAPSIG
       REAL(8),INTENT(OUT):: TAUSIG
       REAL(8),PARAMETER  :: RHOMIN=1.D-10
       REAL(8)            :: PI,N,G2,S,Q,ALPHA
!      *****************************************************************
       TAUSIG=0.D0
       IF(RHOSIG.LT.RHOMIN) RETURN
       PI=4.D0*ATAN(1.D0)
       N=2.D0*RHOSIG
       G2=4.D0*MAX(G2SIG,0.D0)
       CALL DEORB_SQ(N,G2,2.D0*LAPSIG,S,Q)
       CALL DEORB_ALPHA(S,Q,ALPHA)
       TAUSIG=0.5D0*(G2/(8.D0*N)+0.3D0*(3.D0*PI**2)**(2.D0/3.D0)*N**(5.D0/3.D0)*ALPHA)
       RETURN
       END SUBROUTINE SCANL_TAUSIG
!
!.........1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SCANL_EDENS_SPIN(RHOUP,RHODN,GUP,GDN,LUP,LDN,TCORR,EXC)
!      *****************************************************************
!      **  SCAN-L EXCHANGE-CORRELATION ENERGY PER VOLUME AT ONE POINT **
!      **  OF A COLLINEAR SPIN DENSITY WITH PARALLEL GRADIENTS (RADIAL**
!      **  ONE-CENTER CASE): GUP,GDN ARE THE SIGNED RADIAL DERIVATIVES **
!      **  OF THE SPIN DENSITIES, LUP,LDN THEIR LAPLACIANS.           **
!      **  EXCHANGE: EXACT SPIN SCALING. CORRELATION: SPIN-POLARIZED   **
!      **  SCAN WITH TAU=TAUUP+TAUDN FROM THE DEORBITALIZATION.       **
!      *****************************************************************
       USE SCAN_X_MODULE, ONLY: SCAN_X_FX
       USE SCAN_C_MODULE, ONLY: SCAN_C_EDENS
       USE DEORB_ALPHA_MODULE, ONLY: DEORB_SQ,DEORB_ALPHA
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: RHOUP,RHODN,GUP,GDN,LUP,LDN
       LOGICAL(4),INTENT(IN) :: TCORR
       REAL(8)   ,INTENT(OUT):: EXC
       REAL(8)   ,PARAMETER  :: RHOMIN=1.D-10
       REAL(8)               :: PI,EXFAC,S,Q,ALPHA,FX,TUP,TDN,EC
!      *****************************************************************
       PI=4.D0*ATAN(1.D0)
       EXFAC=-0.75D0*(3.D0/PI)**(1.D0/3.D0)
       EXC=0.D0
       IF(RHOUP.GE.RHOMIN) THEN
         CALL DEORB_SQ(2.D0*RHOUP,4.D0*GUP**2,2.D0*LUP,S,Q)
         CALL DEORB_ALPHA(S,Q,ALPHA)
         CALL SCAN_X_FX(S,ALPHA,FX)
         EXC=EXC+RHOUP*EXFAC*(2.D0*RHOUP)**(1.D0/3.D0)*FX
       END IF
       IF(RHODN.GE.RHOMIN) THEN
         CALL DEORB_SQ(2.D0*RHODN,4.D0*GDN**2,2.D0*LDN,S,Q)
         CALL DEORB_ALPHA(S,Q,ALPHA)
         CALL SCAN_X_FX(S,ALPHA,FX)
         EXC=EXC+RHODN*EXFAC*(2.D0*RHODN)**(1.D0/3.D0)*FX
       END IF
       IF(TCORR) THEN
         CALL SCANL_TAUSIG(RHOUP,GUP**2,LUP,TUP)
         CALL SCANL_TAUSIG(RHODN,GDN**2,LDN,TDN)
         CALL SCAN_C_EDENS(RHOUP,RHODN,(GUP+GDN)**2,TUP+TDN,EC)
         EXC=EXC+EC
       END IF
       RETURN
       END SUBROUTINE SCANL_EDENS_SPIN



