!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAW
!     **************************************************************************
!     **  THIS IS THE MAIN PAW SUBROUTINE                                     **
!     **    1) READ AND INITIALIZE                                            **
!     **    2) ITERATE                                                        **
!     **    3) CLOSE DOWN                                                     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4)   :: TNWSTR
      LOGICAL(4)   :: TSTOP
      LOGICAL(4)   :: TFIRST
      LOGICAL(4)   :: TLAST
      LOGICAL(4)   :: TPRINT
      LOGICAL(4)   :: TMERMN
      INTEGER(4)   :: NFILO     ! FILE UNIT OF PROTOCOL
      INTEGER(4)   :: NFI0      ! TIME STEP COUNTER AT START
      INTEGER(4)   :: NFI       ! TIME STEP COUNTER
      INTEGER(4)   :: NOMORE    ! 
      INTEGER(4)   :: IPRINT    ! 
      INTEGER(4)   :: NBEG      ! DECIDES IF RESTART FILE IS READ
      REAL(8)      :: DELT
      LOGICAL(4)   :: TCHK
!     **************************************************************************
                              CALL TRACE$PUSH('PAW')
                              CALL TIMING$CLOCKON('INITIALIZATION')
      CALL STOPIT$SETSTARTTIME
!
!     ==================================================================
!     ====  READ CONTROL INPUT DATA FILE "CNTL"                     ====
!     ==================================================================
      CALL READIN(NBEG,NOMORE,IPRINT,DELT,TMERMN,TNWSTR)
!
!     ==================================================================
!     ==  READ STRUCTURAL DATA FROM FILE "STRC"                       ==
!     ==================================================================
      CALL STRCIN

!     ==================================================================
!     ==  SET DIMER DIMENSION                                         ==
!     ==================================================================
      CALL DIMER$INITDIM()

!
!     ==================================================================
!     ==  INITIALIZE ATOMIC SETUPS                                    ==
!     ==================================================================
!     CALL SETUP$READ()   !NOW CALLED FROM STRCIN
!
!     ==================================================================
!     ==  INITIALIZE ATOMS OBJECT                                     ==
!     ==================================================================
      CALL ATOMS$INITIALIZE
!
!     ==================================================================
!     ==  GENERATE G-VECTORS                                          ==
!     ==================================================================
                              CALL TIMING$CLOCKON('WAVES$INITIALIZE')
      CALL WAVES$INITIALIZE
                              CALL TIMING$CLOCKOFF('WAVES$INITIALIZE')
!
!     ==================================================================
!     ==  READ RESTART FILE (WAVE FUNCTION COEFFICIENTS, ETC)         ==
!     ==================================================================
      IF(NBEG.GE.0) THEN
                              CALL TIMING$CLOCKON('RESTART I/O')
        CALL READRESTART
                              CALL TIMING$CLOCKOFF('RESTART I/O')
      END IF
!
!     ================================================================
!     ==  REPORT INPUT DATA                                         ==
!     ================================================================
      CALL IO$REPORT
      CALL PROTREPORT$HEADER()
!
!     ==================================================================
!     ==  RUN TIME STATISTICS FOR THE INITIALIZATION                  ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
                              CALL TIMING$CLOCKOFF('INITIALIZATION')
                              CALL TIMING$PRINT('MONOMER',NFILO)
!     == CLOCK WILL BE RESTARTED AGAIN AT THE END OF THE FIRST ITERATION TO
!     == ACCOUNT FOR THE FACT THAT MANY ROUTINES INIITIALIZE THEMSELFES IN THE 
!     == FIRST ITERATION. THE TIME TO THIS RESET IS NOT PRINTED 
                              CALL TIMING$START
!
!     ==================================================================
!     ==  RETURN IF NO ITERATIONS ARE REQUIRED                        ==
!     ==================================================================
      IF(NOMORE.EQ.0) THEN
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==================================================================
!     ==================================================================
!     == THE BASIC LOOP FOR MOLECULAR DYNAMICS STARTS HERE            ==
!     ==================================================================
!     ==================================================================
      CALL TIMESTEP$GETI4('ISTEP',NFI)
      NFI0=NFI
      TFIRST=.TRUE.
      TLAST=.FALSE.
      TSTOP=.FALSE.
1000  CONTINUE
!
!     ==================================================================
!     ==  ITERATION CONTROL (PROPER STOP ETC. )                       ==
!     ==================================================================
      NFI=NFI+1
      CALL STOPIT$GETL4('STOP',TSTOP)
      IF(TSTOP) TLAST=.TRUE.
      IF(NFI.GE.NFI0+NOMORE) TLAST=.TRUE.
      TPRINT=(MOD(NFI,IPRINT).EQ.0.OR.TFIRST.OR.TLAST)

      CALL HYPERFINE$SETL4('WAKE',TPRINT)
      CALL GRAPHICS$SETL4('WAKE',TPRINT.AND.TLAST)
      CALL OPTEELS$SETL4('ON',TLAST)
      CALL CORE$SETL4('ON',TLAST)
!
!     ==================================================================
!     ==   WRITE RESTART_OUT                                          ==
!     ==================================================================
      IF(TLAST.OR.(TPRINT.AND.(.NOT.TFIRST))) THEN
                              CALL TIMING$CLOCKON('RESTART I/O') 
        CALL WRITERESTART
        CALL WAVES$GETL4('WRITERHO',TCHK)
        IF(TCHK) CALL WAVES$FIXRHOWRITE()
                              CALL TIMING$CLOCKOFF('RESTART I/O')
      END IF
!
!     ==================================================================
!     ==   PERFORM ONE TIME STEP                                      ==
!     ==================================================================
                              CALL TIMING$CLOCKON('TIMESTEP')
!     ==USE THE LINE WITH "NOT TSTOP" TO AVOID AN ADDITIONAL LAST TIME STEP
!     IF(.NOT.TSTOP) CALL TIMESTEP(DELT,TPRINT,NFI,TSTOP)
      CALL DYNOCC$GETL4('DYN',TMERMN)
      CALL TIMESTEP(DELT,TPRINT,NFI,TSTOP)
                              CALL TIMING$CLOCKOFF('TIMESTEP')
!
!     ==================================================================
!     ==  LIVE SCF ITERATION TABLE (DEFAULT OUTPUT, NOT VERBOSITY-    ==
!     ==  GATED -- SEE SCFREPORT_MODULE, BELOW)                       ==
!     ==================================================================
      CALL SCFREPORT$ITERATION(NFI)
      IF(TLAST) CALL SCFREPORT$FINISH(NFI,TSTOP)
      IF(TLAST) CALL FORCESTRESS$WRITE()
!
!     ==================================================================
!     ==   WRITE INFORMATION AND TRAJECTORIES                         ==
!     ==================================================================
!     __ADD TO TRAJECTORIES (TEMPORARY BUFFER)__________________________
      IF(.NOT.TLAST) THEN
        CALL WRITETRAJECTORY(NFI,DELT)
      END IF
!     __ WRITE TRAJECTORY FROM TEMPORARY BUFFER TO FILE_________________
      IF(TPRINT.OR.TLAST) THEN
        CALL TRAJECTORYIO$FLUSHALL
      END IF
!
!     ==================================================================
!     == RESET CLOCKS FOR TIMING                                      ==
!     == THIS IS DONE AFTER THE FIRST ITERATION TO AVOID COUNTING     ==
!     == THE TIME FOR SELF-INITIALIZATION OF OBJECTS                  ==
!     ==================================================================
      IF(TFIRST) THEN
        CALL TIMING$START
      ELSE
        CALL TIMING$COUNT
      END IF
!
!     ==================================================================
!     ==   STOP OR CONTINUE LOOP                                      ==
!     ==================================================================
      TFIRST=.FALSE.
      IF(.NOT.TLAST) GO TO 1000
!     ==================================================================
!     ==================================================================
!     ====   END OF TIME-STEP ITERATION                               ==
!     ==================================================================
!     ==================================================================
!     == PRINT THE PER-CLOCK TIMING BREAKDOWN (SEE VESPAW_TIMING.F90,   ==
!     == TIMING$PRINT) TO THE .PROT FILE, FOR THE MAIN ITERATION LOOP.  ==
!     == TIMING$PRINT ITSELF ALREADY EXISTED AND IS ALREADY CALLED ONCE ==
!     == ABOVE (RIGHT AFTER INITIALIZATION, REPORTING ONLY THE ONE-TIME ==
!     == SETUP COST) -- BUT NOTHING EVER CALLED IT AGAIN TO REPORT THE  ==
!     == ACTUAL PER-ITERATION BREAKDOWN OF THE EXTENSIVE EXISTING       ==
!     == CLOCKON/CLOCKOFF INSTRUMENTATION THROUGHOUT SRC/ (PLANEWAVE$   ==
!     == FFT, W:HPSI AND ITS SUB-TIMERS, VOFRHO:*, CG-M*, ETC.) -- ADDED==
!     == 2026-09-04 WHILE DIGGING INTO WHY VESPAW IS SLOWER THAN VASP   ==
!     == (SEE THE DEV NOTES).                                        ==
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL TIMING$PRINT('MONOMER',NFILO)
      CALL FORCESTRESS$PROTREPORT(NFI,TSTOP)
      CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      MODULE TIMESTEP_MODULE
      INTEGER(4) :: ISTEPNUMBER
      REAL(8)    :: DELTAT
      LOGICAL(4) :: TNEWTHERMOSTAT=.TRUE. !SWITCHES WAVES THERMOSTAT
      END MODULE TIMESTEP_MODULE
!
!     ..................................................................
      SUBROUTINE TIMESTEP$SETR8(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DELTAT') THEN
       DELTAT=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$GETR8(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DELTAT') THEN
        VAL_=DELTAT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$SETI4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'ISTEP') THEN
        ISTEPNUMBER=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$GETI4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'ISTEP') THEN
        VAL_=ISTEPNUMBER
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$SETL4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'NEWTHERMOSTAT') THEN
        TNEWTHERMOSTAT=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$GETL4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'NEWTHERMOSTAT') THEN
        VAL_=TNEWTHERMOSTAT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$GETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP(DELT,TPRINT,NFI,TSTOP)
!     ******************************************************************      
!     ******************************************************************      
      USE TIMESTEP_MODULE ,ONLY : DELTAT,ISTEPNUMBER,TNEWTHERMOSTAT
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: DELT   ! TIME STEP
      LOGICAL(4),INTENT(IN)   :: TPRINT ! FLAG FOR LONG PRINTOUT
      LOGICAL(4),INTENT(IN)   :: TSTOP  ! FLAG FOR LAST TIME STEP
      INTEGER(4),INTENT(INOUT):: NFI    ! TIME STEP COUNTER
      INTEGER(4)              :: NFILO
      LOGICAL(4)              :: TFOR   ! ON/OFF SWITCH FOR ATOMIC MOTION
      LOGICAL(4)              :: TGRA   ! ON/OFF SWITCH FOR GRADIENT CORRECTION
      REAL(8)                 :: ANNEE  ! FRICTION COEFFICIENT FOR WAVE FUNCTIONS
      REAL(8)                 :: ANNER  ! FRICTION COEFFICIENT FOR NUCLEI
      LOGICAL(4)              :: TCHK
      REAL(8)                 :: TEMPINST
      REAL(8)                 :: ENOSE  ! GENERIC THERMOSTAT ENERGY 
      REAL(8)                 :: EKIN   ! GENERIC KINETIC ENERGY
      REAL(8)                 :: SVAR
      LOGICAL(4)              :: TCHK1,TCHK2
      REAL(8)                 :: FAV,FMAX
      INTEGER(4)              :: NAT
      REAL(8)   ,ALLOCATABLE  :: FSYM(:,:)
      REAL(8)                 :: SSYM(3,3)
      CHARACTER(8)            :: OPTIMIZERVAL
      CHARACTER(16)           :: MIXERVAL
      REAL(8)                 :: EDIFFVAL
      REAL(8)                 :: ETOTNOW
      REAL(8)     ,SAVE       :: EPREVEDIFF=1.D+20
      LOGICAL(4)  ,SAVE       :: TFIRSTEDIFF=.TRUE.
!     == HAUX-RESIDUAL-AWARE EDIFF STOPPING (LIVE OCCUPATIONS ONLY): SEE THE  ==
!     == GATE JUST BELOW, AND THE TWIN GATE ON THE MIXTOL/TCONV STOP SIGNAL   ==
!     == IN VESPAW_WAVES1.F90 (WAVES$ETOT). BOTH ARE NEEDED -- EDIFF (HERE)   ==
!     == AND MIXTOL/TCONV (THERE) ARE INDEPENDENT STOP SIGNALS, AND WHICH ONE ==
!     == ACTUALLY FIRES FIRST FOR A GIVEN DECK/MIXER/TOLERANCE COMBINATION    ==
!     == VARIES -- SEE THE DEV NOTES FOR HOW THIS WAS FOUND (INITIALLY    ==
!     == ONLY THE MIXTOL/TCONV GATE WAS ADDED, AND TESTED AS A NO-OP BECAUSE  ==
!     == THIS DECK'S ACTUAL STOP SIGNAL TURNED OUT TO BE EDIFF, NOT MIXTOL).  ==
      LOGICAL(4)  ,SAVE       :: TEDFTLIVEOCCCHECKED2=.FALSE.
      LOGICAL(4)  ,SAVE       :: TEDFTLIVEOCC2=.FALSE.
      CHARACTER(16)           :: CENVVALEDFTOCC2
      LOGICAL(4)              :: TEDFTINI2
      REAL(8)                 :: EDFTTEMPCHECK2
      REAL(8)                 :: EDFTHAUXRESNOW2
      LOGICAL(4)  ,SAVE       :: TEDFTHAUXTOLCHECKED2=.FALSE.
      REAL(8)     ,SAVE       :: EDFTHAUXCONVTOL2=1.D-5
      CHARACTER(32)           :: CENVVALHAUXTOL2
      LOGICAL(4)              :: TEDFTOKTOSTOP
      LOGICAL(4)              :: TRELAX   ! L-BFGS GEOMETRY OPTIMIZATION
      LOGICAL(4)              :: TSCFCONV ! SCF CONVERGED AT THIS GEOMETRY
      LOGICAL(4)              :: TMOVED   ! RELAX HAS PROPOSED A NEW GEOMETRY
!     ******************************************************************
                              CALL TRACE$PUSH('TIMESTEP')
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!      
      DELTAT=DELT       !-> TIMESTEP_MODULE
      ISTEPNUMBER=NFI   !-> TIMESTEP_MODULE
      CALL ENERGYLIST$RESET
!
      CALL ATOMS$GETL4('MOVE',TFOR)
      CALL DFT$GETL4('GC',TGRA)
      CALL WAVES$SETL4('HAMILTON',TPRINT)
!
!     ==================================================================
!     ==  COMMUNICATE TIME STEP TO OBJECTS                            ==
!     ==================================================================
      CALL WAVES$SETR8('TIMESTEP',DELT)
      CALL ATOMS$SETR8('TIMESTEP',DELT)
      CALL CELL$SETR8('DT',DELT)
      CALL DYNOCC$SETR8('TIMESTEP',DELT)
      CALL QMMM$SETR8('TIMESTEP',DELT)
      CALL THERMOSTAT$SELECT('ATOMS') 
      CALL THERMOSTAT$SETR8('TIMESTEP',DELT)
      CALL THERMOSTAT$SELECT('WAVES') 
      CALL THERMOSTAT$SETR8('TIMESTEP',DELT)
      CALL DIALS$SETR8('TIMESTEP',DELT)
!
!     ==================================================================
!     ==================================================================
!     ==  CALCULATE TOTAL ENERGY AND FORCES                           ==
!     ==================================================================
!     ==================================================================
!     == LDA TOTAL ENERGY ==============================================
      CALL WAVES$ETOT
!     == EXTERNAL POTENTIAL ACTING ON ATOMS ============================
      CALL VEXT$APPLY
!     == SPACE-GROUP SYMMETRY: SYMMETRIZE FORCES AND STRESS =============
      CALL SYMMETRY$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL ATOMLIST$NATOM(NAT)
        ALLOCATE(FSYM(3,NAT))
        CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FSYM)
        CALL SYMMETRY$FORCE(NAT,FSYM)
        CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FSYM)
        DEALLOCATE(FSYM)
        CALL CELL$GETR8A('STRESS_I',9,SSYM)
        CALL SYMMETRY$STRESS(SSYM)
        CALL CELL$SETR8A('STRESS_I',9,SSYM)
      END IF
!     == OCCUPATIONS ===================================================
      CALL DYNOCC$GETR8('EPOT',SVAR)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',SVAR)
      CALL ENERGYLIST$SET('OCCUPATIONAL ENTROPY TERM (-TS)',SVAR)
!
!     ==================================================================
!     ==================================================================
!     ==  SET FRICTION FOR ANNEALING SCHEDULE                         ==
!     ==================================================================
!     ==================================================================
      CALL AUTOPI(TCHK)
      IF(.NOT.TSTOP.AND.TCHK) THEN
        CALL STOPIT$SETL4('STOP',TCHK)
        CALL STOPIT$SETCH('REASON','CONVERGED (AUTOPILOT)')
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'STOP SIGNAL FROM AUTOPILOT'
      END IF
!
!     ==================================================================
!     ==  EDIFF: STOP ONCE THE TOTAL ENERGY HAS CONVERGED (CG/DAV/RMM  ==
!     ==  OPTIMIZERS ONLY -- VASP'S OWN EDIFF). COMPLEMENTS THE        ==
!     ==  DENSITY-RESIDUAL MIXTOL CRITERION (WAVES$ETOT/WAVES_MIXRHO_  ==
!     ==  PUL/_ANDERSON): THE RESIDUAL CAN STILL OSCILLATE JUST ABOVE  ==
!     ==  MIXTOL (E.G. A RESTARTED-DIIS MIXER REPEATEDLY RESTARTING)   ==
!     ==  LONG AFTER THE ENERGY ITSELF HAS STOPPED CHANGING, WHICH     ==
!     ==  OTHERWISE FORCES THE RUN ALL THE WAY TO NSTEP FOR NO         ==
!     ==  PHYSICAL BENEFIT.                                            ==
!     ==================================================================
      IF(.NOT.TEDFTLIVEOCCCHECKED2) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_EDFT_LIVE_OCC',CENVVALEDFTOCC2)
        TEDFTLIVEOCC2=(LEN_TRIM(CENVVALEDFTOCC2).GT.0)
        TEDFTLIVEOCCCHECKED2=.TRUE.
      END IF
      IF(.NOT.TEDFTHAUXTOLCHECKED2) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_EDFT_HAUXCONVTOL',CENVVALHAUXTOL2)
        IF(LEN_TRIM(CENVVALHAUXTOL2).GT.0) READ(CENVVALHAUXTOL2,*) EDFTHAUXCONVTOL2
        TEDFTHAUXTOLCHECKED2=.TRUE.
      END IF
      CALL WAVES$GETCH('OPTIMIZER',OPTIMIZERVAL)
!     == MIXER='VARIATIONAL' INTERLEAVES REJECTED LINE-SEARCH TRIALS WITH   ==
!     == ACCEPTED STEPS, SO CONSECUTIVE ENERGIES SAY NOTHING ABOUT         ==
!     == CONVERGENCE -- THAT MIXER APPLIES EDIFF TO ACCEPTED ENERGIES ITSELF==
!     == (WAVES_MIXPOT_VAR) AND THIS CHECK IS SKIPPED.                     ==
      CALL WAVES$GETCH('MIXER',MIXERVAL)
      IF(TRIM(OPTIMIZERVAL).NE.'MD'.AND.TRIM(MIXERVAL).NE.'VARIATIONAL') THEN
        CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOTNOW)
        IF(.NOT.TFIRSTEDIFF) THEN
          CALL WAVES$GETR8('EDIFF',EDIFFVAL)
          IF(.NOT.TSTOP.AND.ABS(ETOTNOW-EPREVEDIFF).LT.EDIFFVAL) THEN
            TEDFTOKTOSTOP=.TRUE.
            IF(TEDFTLIVEOCC2) THEN
              CALL EDFT$ISINI(TEDFTINI2)
              CALL DYNOCC$GETR8('TEMP',EDFTTEMPCHECK2)
              IF(TEDFTINI2.AND.EDFTTEMPCHECK2.GT.0.D0) THEN
                CALL EDFT$GETHAUXRESMAX(EDFTHAUXRESNOW2)
                TEDFTOKTOSTOP=(EDFTHAUXRESNOW2.LT.EDFTHAUXCONVTOL2)
              END IF
            END IF
            IF(TEDFTOKTOSTOP) THEN
              CALL STOPIT$SETL4('STOP',.TRUE.)
              CALL STOPIT$SETCH('REASON','CONVERGED (ENERGY CHANGE BELOW EDIFF)')
              CALL FILEHANDLER$UNIT('PROT',NFILO)
              WRITE(NFILO,*)'STOP SIGNAL: TOTAL ENERGY CONVERGED (EDIFF)'
            END IF
          END IF
        END IF
        EPREVEDIFF=ETOTNOW
        TFIRSTEDIFF=.FALSE.
      END IF
!
!     ==================================================================
!     ==  GEOMETRY OPTIMIZATION (!RDYN/!CELL OPTIMIZER='LBFGS'): THE   ==
!     ==  SCF STOP SIGNAL ONLY MEANS "CONVERGED AT THIS GEOMETRY".     ==
!     ==  RELAX$STEP THEN PROPOSES THE NEXT GEOMETRY, WHICH            ==
!     ==  ATOMS$PROPAGATE AND CELL$PROPAGATE BELOW TAKE OVER, OR SETS  ==
!     ==  THE STOP ITSELF WHEN FORCES AND STRESS ARE CONVERGED.        ==
!     ==================================================================
      CALL RELAX$GETL4('ON',TRELAX)
      IF(TRELAX) THEN
        TSCFCONV=.FALSE.
        IF(.NOT.TSTOP) CALL STOPIT$TAKESCFCONVERGED(TSCFCONV)
        CALL RELAX$STEP(TSCFCONV,TSTOP,TMOVED)
        IF(TMOVED) THEN
!         __ NEW ENERGY SURFACE: RESTART THE SCF HISTORIES________________
          CALL MIXER$NEWGEOMETRY()
          TFIRSTEDIFF=.TRUE.
        END IF
      END IF
!
!     ==================================================================
!     ==  OBTAIN INSTANTANEOUS FRICTION VALUES FROM THERMOSTATS       ==
!     ==  AND COMMUNICATE TO WAVES AND ATOMS OBJECTS                  ==
!     ==================================================================
!
!     == SET THERMOSTAT FRICTION FOR ATOMS =============================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK1)
      IF(TCHK1) THEN
        CALL THERMOSTAT$GETR8('COOLING',ANNER)
        CALL ATOMS$SETR8('FRICTION',ANNER)
      END IF
!
!     == SET THERMOSTAT FRICTION FOR WAVE FUNCTIONS ====================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK2)
      IF(TCHK2) THEN
        CALL THERMOSTAT$SELECT('WAVES')
        CALL THERMOSTAT$GETR8('COOLING',ANNEE)
        CALL WAVES$SETR8('FRICTION',ANNEE)
        CALL ATOMS$SETR8('ANNEE',ANNEE)
      ELSE
        CALL ATOMS$SETR8('ANNEE',0.D0)
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  PROPAGATE:                                                  ==
!     ==================================================================
!     ==================================================================
! 
!     ==================================================================
!     ==  PROPAGATE NUCLEI (CONSTRAINTS FOLLLOW LATER...)             ==
!     ==================================================================
!---DIMERMERGE FIX      
      CALL DIMER$GETL4('DIMER',TCHK)
      IF(TCHK) THEN
         CALL ATOMS$PROPAGATE_DIMER()
      ELSE
         CALL ATOMS$PROPAGATE()
      END IF
!---END DIMERMERGE FIX      
! 
!     ==================================================================
!     ==  PROPAGATE UNIT CELL                                         ==
!     ==================================================================
      CALL CELL$PROPAGATE()
      CALL CELL$GETR8('EPOT',SVAR)
      CALL ENERGYLIST$SET('CELLOSTAT POTENTIAL',SVAR)     
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',SVAR)
! 
!     ==================================================================
!     ==  APPLY CONSTRAINTS TO ATOMIC COORDINATES                     ==
!     ==================================================================
      CALL ATOMS$CONSTRAINTS()
! 
!     ==================================================================
!     ==  PROPAGATE WAVE FUNCTIONS                                    ==
!     ==================================================================
      CALL WAVES$PROPAGATE()
      CALL WAVES$ORTHOGONALIZE()
!
!     ==================================================================
!     ==  PROPAGATE OCCUPATIONS                                       ==
!     ==================================================================
      CALL DYNOCC$PROPAGATE()
! 
!     ==================================================================
!     ==  PROPAGATE THERMOSTAT FOR THE NUCLEI                         ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
!       __EXTRACT KINETIC ENERGY (USING THE TRUE, EFFECTIVE MASSES)_____
        CALL ATOMS$EKIN(EKIN)
        CALL THERMOSTAT$SETR8('EKIN(SYSTEM)',EKIN)
!       __PROPAGATE THERMOSTAT AND RECORD ITS ENERGY____________________
        CALL THERMOSTAT$PROPAGATE()
      END IF
!
!     ==================================================================
!     ==  PROPAGATE THERMOSTAT FOR THE WAVE FUNCTIONS                 ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        IF(TNEWTHERMOSTAT) THEN
          CALL ATOMS$EFFEKIN(EKIN)
          CALL THERMOSTAT$SETR8('TARGET',EKIN)
        END IF
        CALL WAVES$GETR8('EKIN(PSI)',EKIN)
        CALL THERMOSTAT$SETR8('EKIN(SYSTEM)',EKIN)
        CALL THERMOSTAT$PROPAGATE()
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  COLLECT KINETIC ENERGIES                                    ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     __WAVE FUNCTIONS__________________________________________________
!     ==================================================================
      CALL WAVES$GETR8('EKIN(PSI)',EKIN)
      CALL ENERGYLIST$SET('WAVEFUNCTION KINETIC ENERGY',EKIN)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKIN)
!
!     ==================================================================
!     __ NUCLEI__________________________________________________________
!     ==================================================================
!     __CALCULATE KINETIC ENERGY OF THE ATOMS_ _________________________
      CALL ATOMS$EKIN(EKIN)
      CALL ENERGYLIST$SET('IONIC KINETIC ENERGY',EKIN)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKIN)
!     __SUBTRACT THE EFFECTIVE WAVE FUNCTION KINETIC ENERGY_____________
      CALL ATOMS$EFFEKIN(EKIN)
      CALL ENERGYLIST$SET('BO-WAVEFUNCTION KINETIC ENERGY',EKIN)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',-EKIN)
!     __STORE INSTANTANEOUS ATOMIC TEMPERATURE__________________________
      CALL ATOMS$TEMPERATURE(TEMPINST)
      CALL ENERGYLIST$SET('IONIC TEMPERATURE',TEMPINST)
!
!     ==================================================================
!     == UNIT CELL                                                    ==
!     ==================================================================
      CALL CELL$GETR8('EKIN',EKIN)
      CALL ENERGYLIST$SET('CELLOSTAT KINETIC',EKIN)     
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKIN)
!
!     ==================================================================
!     == OCCUPATIONS                                                  ==
!     ==================================================================
      CALL DYNOCC$GETR8('EKIN',SVAR)
      CALL ENERGYLIST$SET('OCCUPATION KINETIC ENERGY',SVAR)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',SVAR)
! 
!     ==================================================================
!     ==  THERMOSTAT ACTING ON THE NUCLEI                             ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL THERMOSTAT$GETR8('ENERGY',ENOSE)
        CALL ENERGYLIST$SET('ATOM THERMOSTAT',ENOSE)
        CALL ENERGYLIST$ADD('CONSTANT ENERGY',ENOSE)
      END IF
!
!     ==================================================================
!     ==  COLLECT ENERGY OF WAVE FUNCTION THERMOSTAT                  ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        IF(TNEWTHERMOSTAT) THEN
          CALL THERMOSTAT$GETR8('EKIN',ENOSE)
        ELSE
          CALL THERMOSTAT$GETR8('ENERGY',ENOSE)
        END IF
        CALL ENERGYLIST$SET('ELECTRON THERMOSTAT',ENOSE)
        CALL ENERGYLIST$ADD('CONSTANT ENERGY',ENOSE)
      END IF
!     EFLUXE=-EKINC*(2.D0*ANNEE/DELT)/EMASS
!
!     ==================================================================
!     ==================================================================
!     == CHECK CONVERGENCE                                            ==
!     ==================================================================
!     ==================================================================
      CALL ATOMS$FORCECRITERION(FAV,FMAX)
!PRINT*,'FORCE FAV=',FAV,' FMAX=',FMAX
!
!     ==================================================================
!     ==================================================================
!     == WRITE INFORMATION                                            ==
!     ==================================================================
!     ==================================================================
!     __KEEP FORCES AND STRESS BEFORE THE SWITCHES BELOW ZERO THEM______
      CALL FORCESTRESS$CAPTURE()
!     __WRITE PROTOCOLL_________________________________________________
      CALL PRINFO(TPRINT,TSTOP,NFI,DELT)
!
!     ==================================================================
!     ==================================================================
!     == UPDATE DYNAMICAL VARIABLES                                   ==
!     ==================================================================
!     ==================================================================
!     __ELECTRONIC WAVE FUNCTIONS_______________________________________
      CALL WAVES$SWITCH()
!     __ATOMIC POSITIONS________________________________________________
      CALL ATOMS$SWITCH()
!     __UNIT CELL_______________________________________________________
      CALL CELL$SWITCH()
!     __OCCUPATIONS_____________________________________________________
      CALL DYNOCC$SWITCH() 
!     __WAVE FUNCTION THERMOSTAT________________________________________
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$SWITCH()
!     __ATOM THERMOSTAT_________________________________________________
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$SWITCH()
!
!     ==================================================================
!     == TURN DIALS                                                   ==
!     ==================================================================
      CALL DIALS$APPLY
                              CALL TRACE$POP()
      RETURN
      END SUBROUTINE TIMESTEP
!
!     ...................................................STOPIT.........
MODULE STOPIT_MODULE
!**                                                                   ** 
!**  SET SWITCH TO INITIATE SOFT-KILL                                 ** 
!**  SWITCH CAN BE SET BY                                             ** 
!**   A) CREATE A PREDEFINED EXITFILE                                 ** 
!**   B) BY THE CODE BY SETTING TSTOP EXPLICITELY                     ** 
!**                                                                   ** 
USE CLOCK_MODULE 
!USE CLOCK_MODULE , ONLY : DATE_TIME,(.LATER.)
LOGICAL(4)     :: TSTOP=.FALSE.   ! INITIATE SOFT-KILL
LOGICAL(4)     :: DISTRIBUTED=.FALSE. ! TSTOP=TRUE ON ALL NODES
LOGICAL(4)     :: TNOTIFY=.TRUE.  ! NOTIFICATION ABOUT STOP REQUIRED
LOGICAL(4)     :: EXITFILEREMOVED=.FALSE.
CHARACTER(64)  :: STOPREASON=' ' ! WHY THE RUN STOPPED (FOR THE .PROT REPORT)
INTEGER(4)     :: RUNTIME=-1         ! RUNTIME IN SECONDS
TYPE(DATE_TIME) :: STARTTIME
END MODULE STOPIT_MODULE
!
!     ..................................................................
      SUBROUTINE STOPIT$SETL4(ID_,VAL_)
!     ****************************************************************** 
!     **  STOPIT$SET                                                  ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     ****************************************************************** 
      IF(TRIM(ID_).EQ.'STOP') THEN
        TSTOP=VAL_
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STOPIT$SETI4(ID_,VAL_)
!     ****************************************************************** 
!     **  STOPIT$SET                                                  ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ****************************************************************** 
      IF(TRIM(ID_).EQ.'RUNTIME') THEN 
        RUNTIME=VAL_    ! RUNTIME IN SECONDS
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STOPIT$GETL4(ID_,VAL_)
!     ****************************************************************** 
!     **  STOPIT$GET                                                  ** 
!     **  USE STOPIT$GET('TSTOP',4,TSTOP) TO SEE IF STOP IS SET       ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(OUT):: VAL_
!     ****************************************************************** 
      IF(TRIM(ID_).EQ.'STOP') THEN
        CALL  STOPIT_UPDATE
        VAL_=TSTOP
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$GETL4')
      END IF
      RETURN
      END
!
!     ...................................................STOPIT.........
      SUBROUTINE STOPIT$SETCH(ID_,VAL_)
!     ** 'REASON': WHY THE RUN STOPS; THE FIRST REASON SET IS KEPT **
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(*),INTENT(IN) :: VAL_
      IF(TRIM(ID_).EQ.'REASON') THEN
        IF(STOPREASON.EQ.' ') STOPREASON=VAL_
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$SETCH')
      END IF
      RETURN
      END
      SUBROUTINE STOPIT$TAKESCFCONVERGED(TCONV)
!     ** GEOMETRY OPTIMIZATION (RELAX OBJECT): AN SCF STOP SIGNAL -- THE   **
!     ** DENSITY RESIDUAL, THE EDIFF OR THE VARIATIONAL-MIXER CRITERION -- **
!     ** MEANS "ELECTRONS CONVERGED AT THIS GEOMETRY", NOT "END THE RUN".  **
!     ** TAKE IT BACK AND REPORT IT AS TCONV. ANY OTHER STOP REQUEST (EXIT **
!     ** FILE, RUN-TIME LIMIT, AUTOPILOT, L-BFGS CONVERGED) STAYS.         **
      USE STOPIT_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(OUT) :: TCONV
      INTEGER(4)             :: NVAL
!     ******************************************************************
      TCONV=TSTOP.AND.(STOPREASON(1:19).EQ.'CONVERGED (DENSITY' &
     &             .OR.STOPREASON(1:25).EQ.'CONVERGED (ENERGY CHANGE' &
     &             .OR.STOPREASON(1:23).EQ.'CONVERGED (VARIATIONAL')
!     == ALL TASKS MUST TAKE THE SAME DECISION ===========================
      NVAL=0
      IF(TCONV) NVAL=1
      CALL MPE$COMBINE('~','+',NVAL)
      TCONV=(NVAL.NE.0)
      IF(TCONV) THEN
        TSTOP=.FALSE.
        STOPREASON=' '
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STOPIT$GETCH(ID_,VAL_)
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(*),INTENT(OUT):: VAL_
      IF(TRIM(ID_).EQ.'REASON') THEN
        VAL_=STOPREASON
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$GETCH')
      END IF
      RETURN
      END
      SUBROUTINE STOPIT$SETSTARTTIME
!     ****************************************************************** 
!     **  STOPIT_UPDATE                                               ** 
!     **                                                              ** 
!     **  TESTS WHETHER TSTOP HAS BECOME TRUE IN THE MEANWHILE        ** 
!     **                                                              ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CALL CLOCK$NOW(STARTTIME)
      RETURN
      END

!
!     ...................................................STOPIT.........
      SUBROUTINE STOPIT_UPDATE
!     ****************************************************************** 
!     **  STOPIT_UPDATE                                               ** 
!     **                                                              ** 
!     **  TESTS WHETHER TSTOP HAS BECOME TRUE IN THE MEANWHILE        ** 
!     **                                                              ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)     :: NVAL
      INTEGER(4)     :: NFILO
      CHARACTER(256) :: EXITFILE=' '
      INTEGER(4)     :: NTASKS,THISTASK
      TYPE(DATE_TIME) :: ENDTIME,NOW
      INTEGER(4)      :: ISVAR
!     ****************************************************************** 
!
!     ==================================================================
!     ==  REMOVE EXIT FILE IN THE FIRST REQUEST                       ==
!     ==================================================================
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        IF(.NOT.EXITFILEREMOVED) THEN
          CALL FILEHANDLER$DELETE('EXIT')
          EXITFILEREMOVED=.TRUE.
        END IF
      ELSE
        EXITFILEREMOVED=.TRUE.
      END IF

!     ==================================================================
!     ==  CHECK IF TIME EXCEEDS LIMIT                                 ==
!     ==================================================================
      IF(RUNTIME.GT.0) THEN
        ENDTIME=STARTTIME
        ISVAR=RUNTIME
        ENDTIME%SECOND=ENDTIME%SECOND+ISVAR
!
        ISVAR=INT(ENDTIME%SECOND/60)
        ENDTIME%SECOND=ENDTIME%SECOND-ISVAR*60
        ENDTIME%MINUTE=ENDTIME%MINUTE+ISVAR
!
        ISVAR=INT(ENDTIME%MINUTE/60)
        ENDTIME%MINUTE=ENDTIME%MINUTE-ISVAR*60
        ENDTIME%HOUR=ENDTIME%HOUR+ISVAR
!
        ISVAR=INT(ENDTIME%HOUR/24)
        ENDTIME%HOUR=ENDTIME%HOUR-ISVAR*24
        ENDTIME%DAY=ENDTIME%DAY+ISVAR
!
!       WARNING! THIS CHOICE STOPS AT THE LAST SECOND OF THE CURRENT MONTH!!!!!!!
        CALL CLOCK$NOW(NOW)
        IF(NOW.LATER.ENDTIME) THEN
          TSTOP=.TRUE.
          IF(STOPREASON.EQ.' ') STOPREASON='RUN-TIME LIMIT REACHED'
        END IF
      END IF
!
!     ==================================================================
!     ==  CHECK IF EXITFILE EXISTS                                    ==
!     ==================================================================
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      IF(.NOT.TSTOP.AND.THISTASK.EQ.1) THEN
        CALL FILEHANDLER$FILENAME('EXIT',EXITFILE)
        INQUIRE(FILE=EXITFILE,EXIST=TSTOP)
        IF(TSTOP.AND.STOPREASON.EQ.' ') STOPREASON='EXIT FILE FOUND'
      END IF
!
!     ==================================================================
!     ==  CHECK WHETHER TSTOP=T FOR ANY TASK                          ==
!     ==================================================================
      IF(.NOT.DISTRIBUTED) THEN
        IF(TSTOP) THEN 
          NVAL=1
        ELSE
          NVAL=0
        ENDIF
        CALL MPE$COMBINE('~','+',NVAL)
        TSTOP=(NVAL.NE.0) 
        DISTRIBUTED=TSTOP
      END IF
!
!     ==================================================================
!     == WRITE MESSAGE IF NOT DONE ALREADY =============================
!     ==================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(TSTOP.AND.TNOTIFY.AND.THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'STOP SIGNAL RECEIVED'
        CALL LIB$FLUSHFILE(NFILO)
        TNOTIFY=.FALSE.
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PRINFO(TPRINT,TLAST,NFI,DELT)
!     **************************************************************************
!     **  REPORTS ON THE PROCESS OF THE SIMULATION AND INVOKES                **
!     **  ANALYSIS ROUTINES                                                   **
!     **                                                                      **
!     **  PRINFO IS CALLED ONCE PER TIMESTEP                                  **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TPRINT  ! FLAG FOR EXTENSIVE PRINTOUT
      LOGICAL(4),INTENT(IN) :: TLAST   ! FLAG FOR ONE-TIME EXECUTION
      INTEGER(4),INTENT(IN) :: NFI     ! TIME STEP NUMBER
      REAL(8)   ,INTENT(IN) :: DELT    ! RIME STEP
      LOGICAL(4),SAVE       :: TFIRST=.TRUE.
      INTEGER(4)            :: NAT
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)               :: PICO
      REAL(8)               :: SECOND
      REAL(8)               :: TIME     ! ACTUAL TIME
      REAL(8)               :: TME1
      REAL(8)               :: ECONS
      REAL(8)               :: ETOT
      REAL(8)               :: EKINP
      REAL(8)               :: EKINC
      REAL(8)               :: ENOSEE
      REAL(8)               :: ENOSEP
      REAL(8)               :: EKINFC
      REAL(8)               :: HEAT
      REAL(8)               :: EEXT
      REAL(8)               :: OCCKIN
      REAL(8)               :: TEMPINST
      REAL(8)               :: CELVIN
      INTEGER(4)            :: ITEMP
      REAL(8)               :: ANNEE
      REAL(8)               :: ANNER
      INTEGER(4)            :: ISVAR
      REAL(8)               :: EFFEKIN
      REAL(8)               :: ECELLPOT
      REAL(8)               :: ECELLKIN
      LOGICAL(4)            :: TCHK,TCHK1
      LOGICAL(4)            :: TQMMM=.FALSE.
      REAL(8)               :: QMMMKIN   ! EKIN OF QM-MM ENVIRONMENT
      REAL(8)               :: QMMMPOT   ! EPOT OF QM-MM ENVIRONMENT
      REAL(8)               :: QMMMTHERM ! ENERGY OF THE QM-MM THERMOSTAT
      REAL(8)               :: QMMMTEMP  ! TEMPERATURE OF THE QM-MM
      LOGICAL(4)            :: TCOSMO
      REAL(8)               :: EKINCOSMO,EPOTCOSMO
      CHARACTER(8)          :: FFTYPE
!     **************************************************************************
                              CALL TRACE$PUSH('PRINFO')
      TIME=DBLE(NFI)*DELT
      CALL ATOMLIST$NATOM(NAT)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!   
!     ==========================================================================
!     ==   HYPERFINE PARAMETERS                                               ==
!     ==========================================================================
      IF(TPRINT) CALL HYPERFINE$PRINT
!   
!     ==========================================================================
!     ==  PLOT WAVE FUNCTIONS, DENSITIES                                      ==
!     ==========================================================================
!     == GRAPHICS OBJECT IS SET INTO THE WAKE STATE ONLY IF TLAST.AND.TPRINT
                             CALL TRACE$PASS('BEFORE GRAPHICS$PLOT')
      CALL GRAPHICS$PLOT()
                             CALL TRACE$PASS('AFTER GRAPHICS$PLOT')
!
!     ==========================================================================
      IF(TPRINT) THEN
!       == OPTEELS OBJECT IS SWITCHED ON ONLY IN THE LAST ITERATION ============
        CALL OPTEELS$PLOT()
      END IF
!
!     ==========================================================================
!     ==   WRITE HEADER FOR PROTOCOLL FOR EACH TIME STEP                      ==
!     ==========================================================================
                             CALL TRACE$PASS('WRITE HEADER')
      IF(THISTASK.EQ.1.AND.TFIRST) THEN
        CALL COSMO$GETL4('ON',TCOSMO)
        IF(TQMMM) THEN
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A10,2A13,2A6,3A10)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER','T(ENV)','EPOT(ENV)','ECNS(ENV)'
        ELSE IF (TCOSMO) THEN
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A10,2A13,2A8,2A11)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER','COSMO-EPOT','COSMO-EKIN'
        ELSE
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A10,2A13,2A8)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER'
        END IF
        TFIRST=.FALSE.
      END IF
      IF(TPRINT) TFIRST=.TRUE.
!   
!     ==========================================================================
!     ==   WRITE PROTOCOLL FOR EACH TIME STEP                                 ==
!     ==========================================================================
                             CALL TRACE$PASS('WRITE PROTOCOLL')
      IF(THISTASK.EQ.1) THEN

        CALL CONSTANTS('PICO',PICO)
        CALL CONSTANTS('SECOND',SECOND)
        CALL CONSTANTS('KB',CELVIN)
        TME1=TIME/(PICO*SECOND)
!       ========================================================================
!       == ADD UP CONSERVED ENERGY                                            ==
!       ========================================================================
        ECONS=0.D0
!
!       == BASIC LDA + ATOMIC AND FICTITIOUS ELECTRONIC KINETIC ENERGY =
        CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)     
        CALL ENERGYLIST$RETURN('IONIC KINETIC ENERGY',EKINP)
        CALL ENERGYLIST$RETURN('WAVEFUNCTION KINETIC ENERGY',EKINC)     
        CALL ENERGYLIST$RETURN('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
        ECONS=ECONS+EKINC-EFFEKIN+EKINP+ETOT
!
!       == ELECTRON AND ATOM THERMOSTATS =======================================
        CALL ENERGYLIST$RETURN('CELLOSTAT KINETIC',ECELLKIN)     
        CALL ENERGYLIST$RETURN('CELLOSTAT POTENTIAL',ECELLPOT)     
        ECONS=ECONS+ECELLKIN+ECELLPOT
!PRINT*,'ECELLKIN/POT ',ECELLKIN,ECELLPOT,ECELLKIN+ECELLPOT
!
!       == ELECTRON AND ATOM THERMOSTATS =======================================
        CALL ENERGYLIST$RETURN('ATOM THERMOSTAT',ENOSEP)     
        CALL ENERGYLIST$RETURN('ELECTRON THERMOSTAT',ENOSEE)     
        ECONS=ECONS+ENOSEP+ENOSEE
!
        CALL ENERGYLIST$RETURN('CONSTRAINT KINETIC ENERGY',EKINFC)     
        ECONS=ECONS+EKINFC
!
!       == OCCUPATIONS =========================================================
        CALL ENERGYLIST$RETURN('EPOT',HEAT) ! -T*S_MERMIN-MU*N-B*S
        CALL ENERGYLIST$RETURN('OCCUPATION KINETIC ENERGY',OCCKIN)
        ECONS=ECONS+HEAT+OCCKIN
!
!       == QM-MM ENVIRONMENT ===================================================
        CALL QMMM$GETL4('ON',TQMMM)
        IF(TQMMM) THEN
          CALL QMMM$GETR8('EKIN',QMMMKIN)
          CALL QMMM$GETR8('EPOT',QMMMPOT)
          ECONS=ECONS+QMMMKIN   !POTENTIAL ENERGY ALREADY INCLUDED IN ETOT
!
!         == THERMOSTAT FOR THE ENVIRONMENT ====================================
          CALL QMMM$GETR8('ETHERM',QMMMTHERM)
          ECONS=ECONS+QMMMTHERM
        END IF
!
!       == COSMO ===============================================================
        CALL COSMO$GETL4('ON',TCOSMO)
        IF(TCOSMO) THEN
          CALL ENERGYLIST$RETURN('COSMO KINETIC ENERGY',EKINCOSMO)
          CALL ENERGYLIST$RETURN('COSMO POTENTIAL ENERGY',EPOTCOSMO)
!IS ALREADY CONTAINED IN ECONS AND ETOT
!          ECONS=ECONS+EKINCOSMO+EPOTCOSMO
!          ETOT=ETOT+EPOTCOSMO
        END IF
!
!       == EXTERNAL POTENTIAL===================================================
        CALL ENERGYLIST$RETURN('EXTERNAL POTENTIAL',EEXT)
        ECONS=ECONS+EEXT
!
!       == SOME OTHER STUFF ====================================================
        CALL ENERGYLIST$RETURN('IONIC TEMPERATURE',TEMPINST)
        ITEMP=NINT(TEMPINST/CELVIN)
        CALL WAVES$GETR8('FRICTION',ANNEE)
        CALL ATOMS$GETR8('FRICTION',ANNER)
!
        IF (TCOSMO) THEN
          WRITE(NFILO,FMT='("!>",I6,F10.5,1X,I5,F10.5,2F13.5,2F8.5' &
     &                     //',2F11.5)') &
     &               NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER &
     &              ,EPOTCOSMO,EKINCOSMO
        ELSE IF(TQMMM) THEN
          CALL CONSTANTS('KB',CELVIN)
          CALL QMMM$GETI4('NAT:ENV',ISVAR)
          QMMMTEMP=2.D0*QMMMKIN/REAL(3*ISVAR,KIND=8)/CELVIN
          WRITE(NFILO,FMT='("!>",I5,F10.5,1X,I5,F10.5,2F13.5,2F6.3' &
     &                //',I10,2F10.5)') &
     &                NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER &
     &               ,NINT(QMMMTEMP),QMMMPOT,QMMMKIN+QMMMPOT
        ELSE
          WRITE(NFILO,FMT='("!>",I6,F10.5,1X,I5,F10.5,2F13.5,2F8.5)') &
     &                NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER
        ENDIF
!
!       ========================================================================
!       ==   WRITE ENERGIES TO PROTOCOLL                                      ==
!       ========================================================================
                             CALL TRACE$PASS('WRITE ENERGIES')
        IF(TPRINT) THEN
          CALL ENERGYLIST$PRINTHEADER(NFILO)     
!         ==  BASIC LDA ========================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'TOTAL ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  KINETIC')     
!         == AE ELECTROSTATIC ENERGY DOES NOT INCLUDE ISOLATE ENERGY  ==========
!         == (PREVIOUSLY IT WAS ADDED)                             =============
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  EXCHANGE-CORRELATION')     
          CALL ENERGYLIST$PRINTONE(NFILO,'    LOCAL CORRELATIONS')     
          CALL ENERGYLIST$PRINTONE(NFILO,'BACKGROUND')     
          CALL ENERGYLIST$PRINTONE(NFILO,'ISOLATE ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  KINETIC')     
!         == PS ELECTROSTATIC ENERGY DOES NOT INCLUDE ISOLATE ENERGY  ==========
!         == (PREVIOUSLY IT WAS ADDED)                             =============
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  EXCHANGE-CORRELATION')     
          CALL ENERGYLIST$PRINTONE(NFILO,'CORE RELAXATION')     
!         == ATOM KINETIC ENERGY ===============================================
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC TEMPERATURE')
          CALL ENERGYLIST$PRINTONE(NFILO,'ATOM THERMOSTAT')     
!
          CALL ENERGYLIST$PRINTONE(NFILO,'WAVEFUNCTION KINETIC ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'BO-WAVEFUNCTION KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'ELECTRON THERMOSTAT')     
          CALL ENERGYLIST$PRINTONE(NFILO,'EXTERNAL 1CENTER POTENTIAL')     
          CALL ENERGYLIST$PRINTONE(NFILO,'LMTO INTERFACE')     
          CALL ENERGYLIST$PRINTONE(NFILO,'DMFT INTERFACE')     
!
!         == VAN DER WAALS ENERGY===============================================
          CALL VDW$GETL4('ON',TCHK)
          IF(TCHK) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'VAN DER WAALS ENERGY')
          END IF
!
!         == CELLOSTAT   =======================================================
          CALL CELL$GETL4('ON',TCHK)
          CALL CELL$GETL4('MOVE',TCHK1)
          IF(TCHK.AND.TCHK1) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'CELLOSTAT KINETIC')
            CALL ENERGYLIST$PRINTONE(NFILO,'CELLOSTAT POTENTIAL')
          END IF
!          
!         == OCCUPATIONS =======================================================
          CALL DYNOCC$GETL4('DYN',TCHK)
          IF(TCHK) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'OCCUPATIONAL ENTROPY TERM (-TS)')
            CALL ENERGYLIST$PRINTONE(NFILO,'OCCUPATION KINETIC ENERGY')
          END IF
!
!         == EXTERNAL POTENTIAL ================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'EXTERNAL POTENTIAL')
!
!         == QM-MM =============================================================
          CALL QMMM$GETL4('ON',TQMMM)
          IF(TQMMM) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'QMMM KINETIC ENERGY')
            CALL ENERGYLIST$PRINTONE(NFILO,'QMMM POTENTIAL ENERGY')
          END IF
!
!         == COSMO =============================================================
          IF(TCOSMO) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'COSMO KINETIC ENERGY')
            CALL ENERGYLIST$PRINTONE(NFILO,'COSMO POTENTIAL ENERGY')
          END IF
        END IF
      END IF
!   
!     ==========================================================================
!     ==   PROJECTED DENSITY OF STATES                                        ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE STATEANALYSIS')
      IF(TPRINT) THEN
        CALL WAVES$WRITEPDOS
        CALL ATOMLIST$REPORT(NFILO)
        CALL QMMM$REPORT(NFILO)
      ENDIF
                              CALL TRACE$PASS('AFTER STATEANALYSIS')
!   
!     ==========================================================================
!     ==   BANDDATA                                                           ==
!     ==========================================================================
      IF(TPRINT)THEN
        CALL BANDDATA$WRITEFILE
      ENDIF
!   
!     ==========================================================================
!     ==   WRITE FILE FOR COSMOTHERM                                          ==
!     ==========================================================================
      CALL COSMO$PRINTOUT()
!
!     ==========================================================================
!     ==   WRITE OCCUPATIONS AND ENERGIES TO PROTOCOLL                        ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE OCCUPATIONS')
      IF(TPRINT) THEN
        WRITE(NFILO,FMT='()')
        CALL DYNOCC$GETL4('DYN',TCHK)
        IF(TCHK) THEN
          CALL DYNOCC$REPORT(NFILO)
        ELSE
          CALL WAVES$REPORTEIG(NFILO)
        END IF
        CALL CORE$REPORT(NFILO)
!
!       ==================================================================
!       ==  OPTIONAL HDF5 RESULTS CONTAINER (!CONTROL!HDF5OUT) -- A NO-OP ==
!       ==  UNLESS THAT BLOCK REQUESTED IT (src/vespaw_ioroutines.f90).  ==
!       ==  DELIBERATELY UNGUARDED BY "IS THIS THE LAST STEP" HERE --    ==
!       ==  PRINFO's OWN TLAST ARGUMENT IS BOUND TO TSTOP (CONVERGENCE   ==
!       ==  ONLY, SEE THE CALL IN TIMESTEP BELOW), MISSING THE "REACHED  ==
!       ==  NSTEP WITHOUT CONVERGING" CASE THAT ONLY THE OUTER PAW LOOP  ==
!       ==  KNOWS ABOUT (VESPAW_DRIVER.F90, SUBROUTINE PAW) -- AND BY    ==
!       ==  THE TIME PAW's LOOP LEARNS THAT, WAVES$SWITCH() (CALLED      ==
!       ==  LATER IN THIS SAME TIMESTEP CALL, AFTER THIS BLOCK) HAS      ==
!       ==  ALREADY INVALIDATED THIS%EIGVAL -- CONFIRMED BY TWO REAL     ==
!       ==  FAILED ATTEMPTS: GUARDING BY PRINFO's TLAST ARGUMENT SILENTLY==
!       ==  NEVER WROTE ANYTHING ON A RUN THAT REACHED NSTEP NATURALLY;  ==
!       ==  MOVING THE CALL TO AFTER TIMESTEP RETURNS IN PAW CRASHED     ==
!       ==  WITH "EIGENVALUES NOT PRESENT" EVEN ON A RUN THAT DID        ==
!       ==  CONVERGE. SO: FIRE ON *EVERY* TPRINT INSTEAD (H5FCREATE_F    ==
!       ==  TRUNCATES EACH TIME) -- TPRINT IS GUARANTEED TRUE WHENEVER   ==
!       ==  PAW's TLAST IS (TPRINT=(...OR.TLAST) THERE), SO THE ACTUAL   ==
!       ==  FINAL CALL ALWAYS INCLUDES THE TRUE FINAL STATE; ANY EARLIER ==
!       ==  NWRITE-CHECKPOINT CALLS ARE JUST REDUNDANT, HARMLESS         ==
!       ==  OVERWRITES, NOT A CORRECTNESS ISSUE. UNCONDITIONAL ON        ==
!       ==  THISTASK, NOT JUST TASK 1 -- HDF5IO$WRITE's OWN EIGENVALUE   ==
!       ==  GATHER IS AN MPI-COLLECTIVE CALL, SAME AS WAVES$REPORTEIG    ==
!       ==  JUST ABOVE.                                                  ==
!       ==================================================================
        CALL HDF5IO$WRITE()
!       == EIGENVALUES/OCCUPATIONS FOR THE .PROT K-POINT TABLE (COLLECTIVE) ==
        CALL PROTREPORT$CAPTUREEIG()
      END IF
!
!     ==========================================================================
!     ==   WRITE CONSTRAINT INFORMATION TO PROTOCOLL                          ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE CONSTRAINTS')
      IF(TPRINT) THEN
        CALL CONSTRAINTS$REPORT(NFILO,'FULL')
      END IF
                              CALL TRACE$PASS('BEFORE FLUSH')
!   
!     ==========================================================================
!     ==   WRITE CONSTRAINT INFORMATION                                       ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CONSTRAINTS',NFIL)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL,FMT='(20("="),2X,"TIMESTEP: ",I10,2X,20("="))')NFI
        CALL CONSTRAINTS$REPORT(NFIL,'SHORT')
        CALL LIB$FLUSHFILE(NFIL)
      END IF
!   
!     ==========================================================================
!     ==   WRITE FILE STRC_OUT                                                ==
!     ==========================================================================
      CALL STRCOUT
      CALL QMMM$GETL4('ON',TQMMM)
      IF(TQMMM) THEN
        CALL FORCEFIELD$GETCH('FORCEFIELD',FFTYPE)
        IF(FFTYPE.EQ.'AMBER') THEN
           CALL FORCEFIELD$WRITE_MMSTRC
        ELSE IF(FFTYPE.EQ.'UFF') THEN
           CALL FORCEFIELD$WRITE_UFFSTRC
        END IF
      END IF
!
!     ==========================================================================
!     ==   CLOSE PROTOCOL FILE                                                ==
!     ==========================================================================
      CALL FILEHANDLER$CLOSE('PROT')
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITETRAJECTORY(NFI,DELT)
!     **************************************************************************
!     **  transfers trajectory data into the trajectory object                **
!     **  data will be accumulated an be written in bigu chunks               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFI
      REAL(8)   ,INTENT(IN)  :: DELT
      REAL(8)                :: CELVIN
      REAL(8)                :: GLIB
      REAL(8)                :: PICO
      REAL(8)                :: SECOND
      INTEGER(4)             :: NTASKS,THISTASK
      REAL(8)   ,ALLOCATABLE :: DWORK(:) 
      REAL(8)                :: TIME
      REAL(8)                :: SVAR
      INTEGER(4)             :: NAT,NATM
      REAL(8)                :: ECONS
      REAL(8)                :: ETOT
      REAL(8)                :: EKINP
      REAL(8)                :: EKINC
      REAL(8)                :: EFFEKIN
      REAL(8)                :: ENOSEP
      REAL(8)                :: ENOSEE
      REAL(8)                :: EKINFC
      REAL(8)                :: HEAT
      REAL(8)                :: OCCKIN
      LOGICAL(4)             :: TQMMM,TCHK
      REAL(8)                :: QMMMKIN,QMMMPOT,QMMMTHERM
      REAL(8)                :: EEXT
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
                                 CALL TRACE$PUSH('WRITETRAJECTORY')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL CONSTANTS('PICO',PICO)
      CALL CONSTANTS('SECOND',SECOND)
      CALL CONSTANTS('KB',CELVIN)
      TIME=DBLE(NFI)*DELT
      CALL ATOMLIST$NATOM(NAT)
!
!     ==========================================================================
!     ==  COLLECT ENERGIES                                                    ==
!     ==========================================================================
!     ==========================================================================
!     == ADD UP CONSERVED ENERGY                                              ==
!     ==========================================================================
      ECONS=0.D0
!     
!     == BASIC LDA + ATOMIC AND FICTITIOUS ELECTRONIC KINETIC ENERGY ===========
      CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)     
      CALL ENERGYLIST$RETURN('IONIC KINETIC ENERGY',EKINP)
      CALL ENERGYLIST$RETURN('WAVEFUNCTION KINETIC ENERGY',EKINC)     
      CALL ENERGYLIST$RETURN('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
      ECONS=ECONS+EKINC-EFFEKIN+EKINP+ETOT
!     
!     == ELECTRON AND ATOM THERMOSTATS =========================================
      CALL ENERGYLIST$RETURN('ATOM THERMOSTAT',ENOSEP)     
      CALL ENERGYLIST$RETURN('ELECTRON THERMOSTAT',ENOSEE)     
      ECONS=ECONS+ENOSEP+ENOSEE
!     
      CALL ENERGYLIST$RETURN('CONSTRAINT KINETIC ENERGY',EKINFC)     
      ECONS=ECONS+EKINFC
!     
!     == OCCUPATIONS ===========================================================
      CALL ENERGYLIST$RETURN('OCCUPATIONAL ENTROPY TERM (-TS)',HEAT) ! -T*S_MERMIN
      CALL ENERGYLIST$RETURN('OCCUPATION KINETIC ENERGY',OCCKIN)
      ECONS=ECONS+HEAT+OCCKIN
!     
!     == QM-MM ENVIRONMENT =====================================================
      CALL QMMM$GETL4('ON',TQMMM)
      IF(TQMMM) THEN
        CALL QMMM$GETR8('EKIN',QMMMKIN)
        CALL QMMM$GETR8('EPOT',QMMMPOT)
        ECONS=ECONS+QMMMKIN   !POTENTIAL ENERGY ALREADY INCLUDED IN ETOT
!     
!       == THERMOSTAT FOR THE ENVIRONMENT ======================================
        CALL QMMM$GETR8('ETHERM',QMMMTHERM)
        ECONS=ECONS+QMMMTHERM
      END IF
!
!     == EXTERNAL POTENTIAL=====================================================
      CALL ENERGYLIST$RETURN('EXTERNAL POTENTIAL',EEXT)
      ECONS=ECONS+EEXT
!   
!     ==========================================================================
!     ==   WRITE ENERGY TRAJECTORY                                            ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE E-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('ENERGY-TRAJECTORY')
      CALL TRAJECTORYIO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        ALLOCATE(DWORK(MAX(4*NAT,8)))
        CALL CONSTANTS('KB',CELVIN)
        GLIB=3*NAT-3
        IF(GLIB.NE.0.D0) THEN 
          SVAR=NINT(EKINP/(.5D0*GLIB*CELVIN))
        ELSE
          SVAR=0.D0
        END IF
        DWORK(1)=SVAR
        DWORK(2)=EKINC
        DWORK(3)=EKINP
        DWORK(4)=ETOT
        DWORK(5)=ECONS
        DWORK(6)=ENOSEE
        DWORK(7)=ENOSEP
        DWORK(8)=HEAT
        CALL TRAJECTORYIO$ADD(NFI,TIME,8,DWORK)
        DEALLOCATE(DWORK)
      END IF
      CALL TRAJECTORYIO$SELECT('NONE')
!   
!     ==========================================================================
!     ==   WRITE POSITION TRAJECTORY                                          ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE R-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('POSITION-TRAJECTORY')
      CALL TRAJECTORYIO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        ALLOCATE(DWORK(9+8*NAT))
        CALL CELL$GETR8A('T(0)',9,DWORK(1:9))
        CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,DWORK(9+1:9+3*NAT))
        CALL ATOMLIST$GETR8A('Q',0,NAT,DWORK(9+3*NAT+1:9+3*NAT+NAT))
        CALL ATOMLIST$GETR8A('CHARGEANDMOMENTS',0,4*NAT &
     &                                  ,DWORK(9+3*NAT+NAT+1:9+3*NAT+NAT+4*NAT))
        CALL TRAJECTORYIO$ADD(NFI,TIME,9+8*NAT,DWORK)
        DEALLOCATE(DWORK)
      END IF
      CALL TRAJECTORYIO$SELECT('NONE')
!   
!     ==========================================================================
!     ==   WRITE FORCE TRAJECTORY                                             ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE F-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('FORCE-TRAJECTORY')
      CALL TRAJECTORYIO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        ALLOCATE(DWORK(4*NAT))
        CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,DWORK)
        DWORK(3*NAT+1:4*NAT)=0.D0 ! SHALL CONTAIN IN FUTURE THE POTENTIALS
        CALL TRAJECTORYIO$ADD(NFI,TIME,4*NAT,DWORK)
        DEALLOCATE(DWORK)
      END IF
      CALL TRAJECTORYIO$SELECT('NONE')
!   
!     ==========================================================================
!     ==   WRITE POSITION TRAJECTORY FOR QMMM                                 ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE QM-MM R-TRAJECTORY')
      CALL QMMM$GETL4('ON',TQMMM)
      IF(TQMMM) THEN
         CALL TRAJECTORYIO$SELECT('QMMM-POS-TRA')
         CALL TRAJECTORYIO$GETL4('ON',TCHK)
         IF(TCHK) THEN      
           CALL CLASSICAL$SELECT('QMMM')
           CALL CLASSICAL$GETI4('NAT',NATM)
           ALLOCATE(DWORK(9+4*NATM))
           CALL CELL$GETR8A('T(0)',9,DWORK(1:9))
           CALL CLASSICAL$GETR8A('R(0)',3*NATM,DWORK(10:9+3*NATM))
           CALL CLASSICAL$GETR8A('QEL',NATM,DWORK(10+3*NATM:))

           CALL TRAJECTORYIO$ADD(NFI,TIME,9+4*NATM,DWORK)
           DEALLOCATE(DWORK)
         END IF
         CALL TRAJECTORYIO$SELECT('NONE')
      END IF
!
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE SCFREPORT_MODULE
!     **************************************************************************
!     ** STATE FOR THE COMPACT LIVE SCF-ITERATION TABLE. THIS TABLE IS THE     **
!     ** DEFAULT LIVE VIEW OF A RUNNING CALCULATION (SEE THE DEV NOTES'S   **
!     ** LIVE-SCF-OUTPUT DESIGN) -- IT IS *NOT* GATED BY VERBOSITY_MODULE,     **
!     ** UNLIKE THE OLD PER-BAND/PER-CYCLE PRINT* DIAGNOSTICS IN               **
!     ** VESPAW_MIXER.F90 AND VESPAW_CG.F90, WHICH ARE NOW HIDDEN BEHIND       **
!     ** VERBOSITY='DEBUG' AND REPLACED BY THIS TABLE AS THE DEFAULT OUTPUT.   **
!     **************************************************************************
      LOGICAL(4) ,SAVE :: TFIRSTCALL=.TRUE.
      REAL(8)    ,SAVE :: EOLD=0.D0
      REAL(8)    ,SAVE :: TIMEOLD=0.D0
      END MODULE SCFREPORT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCFREPORT$ITERATION(NFI)
!     **************************************************************************
!     ** PRINTS ONE ROW OF THE LIVE SCF TABLE. CALLED ONCE PER OUTER ITERATION **
!     ** FROM SUBROUTINE PAW, RIGHT AFTER CALL TIMESTEP RETURNS -- UNLIKE      **
!     ** HDF5IO$WRITE (SEE PRINFO, ABOVE) THIS HAS NO WAVE-FUNCTION-LIFETIME   **
!     ** CONSTRAINT SINCE IT ONLY NEEDS THE ENERGYLIST/MIXER RESIDUAL, SO IT   **
!     ** CAN SAFELY RUN ON EVERY ITERATION RATHER THAN BEING TIED TO TPRINT.   **
!     **************************************************************************
      USE SCFREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFI
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)               :: ETOT,RESIDUAL,TNOW,DT,EV
      CHARACTER(12)         :: DESTR
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)
      CALL MIXER$GETR8('RESIDUAL',RESIDUAL)
      CALL PROTREPORT$ITERATION(NFI,ETOT,RESIDUAL)
      CALL CONSTANTS('EV',EV)
      CALL TIMING_CLOCK(TNOW)
      IF(TFIRSTCALL) THEN
        WRITE(*,FMT='()')
        WRITE(*,FMT='("SCF ITERATION")')
        WRITE(*,FMT='(72("-"))')
        WRITE(*,FMT='(A5,A20,A16,A14,A11)') &
     &        'ITER','ENERGY [EV]','DE [EV]','RESIDUAL','TIME'
        WRITE(*,FMT='(72("-"))')
        DESTR='--'
        DT=0.D0
      ELSE
        WRITE(DESTR,FMT='(ES12.3)') (ETOT-EOLD)/EV
        DT=TNOW-TIMEOLD
      END IF
      WRITE(*,FMT='(I5,F20.8,2X,A12,2X,ES12.3,F9.2," S")') &
     &      NFI,ETOT/EV,ADJUSTL(DESTR),RESIDUAL,DT
      EOLD=ETOT
      TIMEOLD=TNOW
      TFIRSTCALL=.FALSE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCFREPORT$FINISH(NFI,TSTOPPEDEARLY)
!     **************************************************************************
!     ** PRINTS THE ONE-LINE SUMMARY THAT CLOSES THE LIVE SCF TABLE. CALLED    **
!     ** ONCE FROM SUBROUTINE PAW WHEN ITS OWN (CORRECTLY-COMPUTED) TLAST      **
!     ** BECOMES TRUE. TSTOPPEDEARLY IS STOPIT$GETL4('STOP',...) -- TRUE IF    **
!     ** THE RUN ENDED BEFORE REACHING NSTEP, FALSE IF IT SIMPLY RAN OUT OF    **
!     ** ITERATIONS. CORRECTED 2026-09-04: THIS FLAG DOES **NOT** MEAN ONLY    **
!     ** AN EXTERNAL SOFT-KILL (STOP FILE/TIME LIMIT) -- STOPIT$SETL4('STOP',  **
!     ** ...) IS ALSO SET BY THE AUTOPILOT (!AUTO/MIXTOL) AND BY THE EDIFF     **
!     ** TOTAL-ENERGY-CONVERGENCE CHECK (BOTH IN SUBROUTINE TIMESTEP, ABOVE)   **
!     ** -- ALL FOUR REASONS SHARE THE SAME FLAG, SO THIS ROUTINE CANNOT TELL  **
!     ** THEM APART (THE .PROT FILE CAN: EACH SITE WRITES ITS OWN "STOP        **
!     ** SIGNAL: ..." LINE THERE). NEITHER CASE HERE IMPLIES ELECTRONIC/IONIC  **
!     ** CONVERGENCE WAS REACHED -- THAT IS A SEPARATE, INNER CHECK (SEE THE   **
!     ** ETOTNOW TEST INSIDE SUBROUTINE TIMESTEP), NOT EXPOSED HERE EITHER.    **
!     **************************************************************************
      USE SCFREPORT_MODULE, ONLY : TFIRSTCALL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFI
      LOGICAL(4),INTENT(IN) :: TSTOPPEDEARLY
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)               :: ETOT,EV
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      IF(TFIRSTCALL) RETURN   ! SCFREPORT$ITERATION NEVER RAN (NOMORE=0 ETC.)
      CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)
      CALL CONSTANTS('EV',EV)
      ETOT=ETOT/EV
      WRITE(*,FMT='(72("-"))')
      IF(TSTOPPEDEARLY) THEN
        WRITE(*,FMT='("STOPPED EARLY (see .prot for the reason -- ",'// &
     &        '"autopilot/EDIFF convergence or an external stop) AFTER ",'// &
     &        'I0," ITERATIONS | E = ",F20.8," EV")') NFI,ETOT
      ELSE
        WRITE(*,FMT='("FINISHED AFTER ",I0," ITERATIONS | E = ",'// &
     &        'F20.8," EV")') NFI,ETOT
      END IF
      WRITE(*,FMT='()')
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE FORCESTRESSREPORT_MODULE
!     **************************************************************************
!     ** FINAL FORCES AND STRESS AT FULL PRECISION, ALSO FOR FIXED ATOMS AND  **
!     ** A FIXED CELL (!CONTROL!GENERIC FORCES=T STRESS=T). FORCES AND THE    **
!     ** INTERNAL STRESS ARE ZEROED BY ATOMS$SWITCH/CELL$SWITCH AFTER EVERY   **
!     ** STEP, SO THEY ARE CAPTURED BEFORE THE SWITCH AND WRITTEN AT THE END. **
!     **************************************************************************
      LOGICAL(4)            ,SAVE :: TFORCES=.FALSE.
      LOGICAL(4)            ,SAVE :: TSTRESS=.FALSE.
      LOGICAL(4)            ,SAVE :: TCAPTURED=.FALSE.
      REAL(8)   ,ALLOCATABLE,SAVE :: FORCE(:,:)  !(3,NAT) [H/ABOHR]
      REAL(8)               ,SAVE :: STRESS(3,3)=0.D0 ! -DE/DEPSILON [H]
      REAL(8)               ,SAVE :: RBAS(3,3)=0.D0
      END MODULE FORCESTRESSREPORT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FORCESTRESS$SETL4(ID,VAL)
      USE FORCESTRESSREPORT_MODULE, ONLY : TFORCES,TSTRESS
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'FORCES') THEN
        TFORCES=VAL
      ELSE IF(ID.EQ.'STRESS') THEN
        TSTRESS=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('FORCESTRESS$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FORCESTRESS$CAPTURE()
!     **************************************************************************
!     ** KEEPS THE FORCES AND STRESS OF THE CURRENT STEP. CALLED FROM TIMESTEP**
!     ** AFTER THE ENERGY EVALUATION AND BEFORE ATOMS$SWITCH/CELL$SWITCH.    **
!     **************************************************************************
      USE FORCESTRESSREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NAT
!     **************************************************************************
      CALL ATOMLIST$NATOM(NAT)
      IF(.NOT.ALLOCATED(FORCE)) ALLOCATE(FORCE(3,NAT))
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCE)
      CALL CELL$GETR8A('STRESS_I',9,STRESS)
      CALL CELL$GETR8A('T(0)',9,RBAS)
      TCAPTURED=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FORCESTRESS$PROTREPORT(NFI,TSTOP)
!     **************************************************************************
!     ** HANDS THE LAST CAPTURED FORCES AND STRESS TO THE .PROT REPORT       **
!     **************************************************************************
      USE FORCESTRESSREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFI
      LOGICAL(4),INTENT(IN) :: TSTOP
      REAL(8)               :: DUMMY(3,1)
!     **************************************************************************
      IF(TCAPTURED) THEN
        CALL PROTREPORT$FINISH(NFI,TSTOP,SIZE(FORCE,2),FORCE,.TRUE.,STRESS,TSTRESS,RBAS)
      ELSE
        DUMMY=0.D0
        CALL PROTREPORT$FINISH(NFI,TSTOP,1,DUMMY,.FALSE.,STRESS,.FALSE.,RBAS)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FORCESTRESS$WRITE()
!     **************************************************************************
!     ** WRITES THE LAST CAPTURED FORCES AND STRESS TO STDOUT AND THE .PROT   **
!     ** FILE. THE STRESS SIGMA=-(1/V) DE/DEPSILON IS POSITIVE WHEN THE CELL  **
!     ** WANTS TO EXPAND; THE PRESSURE IS TR(SIGMA)/3.                        **
!     **************************************************************************
      USE FORCESTRESSREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4)    :: NTASKS,THISTASK
      INTEGER(4)    :: NFILO,NFIL,IU,NAT,IAT,I
      REAL(8)       :: ANGSTROM,EV,GPA,VOL
      CHARACTER(32) :: NAME
!     **************************************************************************
      IF(.NOT.TCAPTURED) RETURN
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL CONSTANTS('EV',EV)
      CALL CONSTANTS('PASCAL',GPA)
      GPA=1.D+9*GPA
      VOL=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3))
      VOL=ABS(VOL)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      DO IU=1,2
        NFIL=6
        IF(IU.EQ.2) NFIL=NFILO
        IF(TFORCES) THEN
          NAT=SIZE(FORCE,2)
          WRITE(NFIL,FMT='(/"FINAL FORCES [EV/AA]")')
          DO IAT=1,NAT
            CALL ATOMLIST$GETCH('NAME',IAT,NAME)
            WRITE(NFIL,FMT='("FORCE ",A10,3F18.10)') &
     &            NAME(1:10),FORCE(:,IAT)/(EV/ANGSTROM)
          ENDDO
        END IF
        IF(TSTRESS) THEN
          WRITE(NFIL,FMT='(/"FINAL STRESS  SIGMA=-(1/V)DE/DEPSILON [GPA]",'// &
     &              '"  AND DE/DEPSILON [EV]")')
          DO I=1,3
            WRITE(NFIL,FMT='("STRESS ",3F14.7,2X,3F18.10)') &
     &            STRESS(I,:)/VOL/GPA,-STRESS(I,:)/EV
          ENDDO
          WRITE(NFIL,FMT='("PRESSURE[GPA] ",F14.7,"  VOLUME[AA^3] ",F16.8)') &
     &          (STRESS(1,1)+STRESS(2,2)+STRESS(3,3))/3.D0/VOL/GPA,VOL/ANGSTROM**3
        END IF
      ENDDO
      RETURN
      END
