!
!*******************************************************************************
!**  EDFT_MODULE -- MINIMAL PROTOTYPE OF ENSEMBLE-DFT-STYLE (MARZARI-        **
!**  VANDERBILT-PAYNE) OCCUPATION HANDLING VIA AN AUXILIARY HAMILTONIAN.     **
!**                                                                          **
!**  MOTIVATION (SEE THE DEV NOTES): REPLACE THE "DIAGONALIZE AT A       **
!**  FROZEN POTENTIAL, FERMI-FILL THE RESULTING EIGENVALUES, MIX THE         **
!**  DENSITY WITH HISTORY, REPEAT" TWO-LEVEL SCF LOOP WITH A SINGLE JOINT    **
!**  VARIATIONAL DESCENT ON THE MERMIN FREE ENERGY, FOLLOWING THE ENSEMBLE-  **
!**  DFT SCHEME OF MARZARI, VANDERBILT AND PAYNE (PRL 79, 1337 (1997)): AN   **
!**  AUXILIARY, UNITARILY-INVARIANT HAMILTONIAN HAUX (NOT THE PHYSICAL ONE)  **
!**  IS DIAGONALIZED TO GENERATE FERMI-DIRAC OCCUPATIONS, AND HAUX ITSELF    **
!**  IS THE OBJECT THAT GETS DRIVEN TOWARD SELF-CONSISTENCY (VIA A           **
!**  GENERALIZED FORCE, NOT IMPLEMENTED YET), RATHER THAN A SEPARATELY-      **
!**  MIXED REAL-SPACE DENSITY.                                               **
!**                                                                          **
!**  ROUTINES IN THIS FILE:                                                  **
!**    EDFT$INI          -- INITIALIZE HAUX=DIAG(EIG) FOR EACH (IKPT,ISPIN)  **
!**    EDFT$INIBLOCK      -- LAZY, PER-BLOCK VARIANT OF EDFT$INI FOR A       **
!**                         CALLER THAT VISITS ONE (IKPT,ISPIN) AT A TIME.   **
!**    EDFT$ISINI         -- WHETHER HAUX HAS BEEN INITIALIZED YET.          **
!**    EDFT$SETHAUXRESMAX -- RECORD THIS ITERATION'S WORST FERMI-WINDOW      **
!**    EDFT$GETHAUXRESMAX    HAUX-TRACKING RESIDUAL, FOR THE LIVE-OCCUPATION **
!**                         STOPPING GATE IN VESPAW_WAVES1.F90.              **
!**    EDFT$HAUXEIGVEC    -- LIKE EDFT$HAUXEIG, ALSO RETURNS HAUX'S OWN      **
!**                         EIGENVECTORS -- THE NATURAL-ORBITAL ROTATION     **
!**                         EDFT$STATE_BY_BLOCK (VESPAW_CG.F90) NEEDS FOR    **
!**                         THE REAL JOINT ORBITAL+OCCUPATION OPTIMIZER.     **
!**    EDFT$OCCUPATIONS  -- DIAGONALIZE HAUX AND FERMI-FILL ITS EIGENVALUES, **
!**                         REUSING DYNOCC_MERMIN'S BISECTION-ON-CHEMICAL-   **
!**                         POTENTIAL LOGIC UNCHANGED.                       **
!**    EDFT$HAUXFORCE    -- GENERALIZED FORCE ON HAUX: THE MISMATCH BETWEEN  **
!**                         THE REAL HAMILTONIAN (PROJECTED INTO THE SAME    **
!**                         BASIS, SUPPLIED BY THE CALLER) AND HAUX ITSELF.  **
!**    EDFT$HAUXSTEP     -- DAMPED LINEAR-MIXING STEP OF HAUX ALONG THAT     **
!**                         FORCE (UNIFORM WEIGHT -- SEE EDFT$HAUXSTEP_MVP   **
!**                         BELOW FOR THE ENTROPY-DERIVATIVE-WEIGHTED ONE).  **
!**    EDFT$HAUXSTEP_MVP -- LIKE EDFT$HAUXSTEP, BUT WEIGHTS THE RESIDUAL BY  **
!**                         THE FERMI-DIRAC DERIVATIVE F_I(1-F_I)/T IN       **
!**                         HAUX'S OWN EIGENBASIS -- THE PIECE OF THE REAL   **
!**                         MARZARI-VANDERBILT-PAYNE GENERALIZED FORCE THAT  **
!**                         EDFT$HAUXSTEP OMITS (SEE ITS OWN HEADER COMMENT  **
!**                         AND THE DEV NOTES FOR THE DERIVATION). FOUND **
!**                         2026-09-05 TO TRACK FERMI-ACTIVE BANDS ~110X     **
!**                         WORSE THAN PLAIN EDFT$HAUXSTEP UNDER A FIXED     **
!**                         LAMBDA -- KEPT FOR PROVENANCE/COMPARISON, NOT    **
!**                         THE DEFAULT SHADOW-MODE STEP ANY MORE.          **
!**    EDFT$HAUXSTEP_ADAPTIVE -- PER-BAND-ADAPTIVE FIX FOR THAT REGRESSION:  **
!**                         SAME WEIGHT SHAPE AS EDFT$HAUXSTEP_MVP, BUT      **
!**                         GATED (MIN(1,W_I/WREF)) RATHER THAN USED AS A    **
!**                         RAW MULTIPLIER, SO ANY MEANINGFULLY FRACTIONAL   **
!**                         BAND GETS THE SAME FULL-STRENGTH LAMBDA AS THE   **
!**                         UNIFORM STEP INSTEAD OF BEING SYSTEMATICALLY     **
!**                         SLOWED DOWN JUST FOR NOT SITTING EXACTLY AT      **
!**                         F=0.5 (SEE ITS OWN HEADER COMMENT). THE CURRENT  **
!**                         DEFAULT SHADOW-MODE STEP WHEN TEMP>0.            **
!**    EDFT$HAUXSTEP_MVPEXACT -- THE ACTUAL, VERIFIED dF/dHAUX (NOT JUST A   **
!**                         WEIGHTED RESIDUAL LIKE EVERY STEP ABOVE):        **
!**                         FOUND 2026-09-05, BY AN INDEPENDENT FINITE-      **
!**                         DIFFERENCE CHECK (TESTS/UNITTESTS/              **
!**                         EDFT_GRADIENT_CHECK), THAT EDFT$HAUXSTEP_        **
!**                         ADAPTIVE'S RESIDUAL MISMATCHES THE TRUE          **
!**                         GRADIENT OF THE MERMIN FREE ENERGY BY 28-107%,   **
!**                         DX-INDEPENDENTLY (A REAL FORMULA BUG, NOT        **
!**                         TRUNCATION NOISE), BOTH NEAR AND FAR FROM        **
!**                         SELF-CONSISTENCY. THIS ROUTINE IMPLEMENTS THE    **
!**                         CORRECTED FORMULA -- A FIXED-ELECTRON-COUNT      **
!**                         LAGRANGE CORRECTION ON THE DIAGONAL, AND MVP'S   **
!**                         OWN DEGENERATE-PERTURBATION-THEORY OFF-DIAGONAL  **
!**                         FORMULA IN PLACE OF THE GEOMETRIC-MEAN GATE      **
!**                         WEIGHT -- AND MATCHES THE FINITE-DIFFERENCE      **
!**                         GRADIENT TO 0.00-0.22% (SEE ITS OWN HEADER       **
!**                         COMMENT AND THE DEV NOTES FOR THE FULL       **
!**                         DERIVATION). NOT YET THE DEFAULT ANYWHERE --     **
!**                         ITS STABILITY AS A LIVE, FED-BACK-INTO-PSI STEP  **
!**                         (VS. EDFT$HAUXSTEP_ADAPTIVE'S ALREADY-           **
!**                         CHARACTERIZED LAMBDA_HAUX BEHAVIOR) IS NOT YET   **
!**                         TESTED.                                         **
!**    EDFT$HAUXSTEP_ANDERSON -- REPLACES THE HAND-TUNED LAMBDA_HAUX        **
!**                         DAMPING CONSTANT ENTIRELY WITH THE SAME         **
!**                         ADAPTIVE, RESTARTED-DIIS ANDERSON MIXING        **
!**                         ALREADY VALIDATED FOR THE DENSITY SELF-         **
!**                         CONSISTENCY LOOP (WAVES_MIXRHO_ANDERSON,        **
!**                         VESPAW_MIXER.F90). TAKES NO LAMBDA ARGUMENT --  **
!**                         THE DAMPING IS ENTIRELY SELF-TUNING, PER        **
!**                         (K-POINT,SPIN) BLOCK (EDFTMIXANDERSON_MODULE,   **
!**                         BELOW). DIRECT RESPONSE TO THE FINDING (SEE     **
!**                         THE DEV NOTES) THAT NO SINGLE LAMBDA_HAUX   **
!**                         VALUE WORKED ACROSS SYSTEMS FOR                 **
!**                         EDFT$HAUXSTEP_MVPEXACT.                         **
!**    EDFT$HAUXSTEP_RESGATE -- FIXES A DIFFERENT BUG IN                    **
!**                         EDFT$HAUXSTEP_ADAPTIVE'S GATE, FOUND 2026-09-06 **
!**                         BY DIRECT INSTRUMENTATION OF HAUX'S OWN         **
!**                         EIGENVALUES ON THE SI TEST CASE (SEE            **
!**                         THE DEV NOTES): THE GATE (MIN(1,W_I/WREF),  **
!**                         W_I FROM THE BAND'S *CURRENT* EIGENVALUE'S      **
!**                         DISTANCE FROM THE CHEMICAL POTENTIAL) HAS NO    **
!**                         WAY TO OPEN FOR A BAND WHOSE CURRENT EIGENVALUE **
!**                         SITS FAR FROM THE FERMI LEVEL BUT WHOSE TRUE    **
!**                         RESIDUAL (HPHYS-HAUX) IS NONETHELESS LARGE --   **
!**                         A SELF-REINFORCING LOCK (CONFIRMED NUMERICALLY: **
!**                         GATE~6.D-5, RESIDUAL~0.17 H, NEITHER MOVING).   **
!**                         FIX: OR THE EXISTING FERMI-WINDOW GATE WITH A   **
!**                         SECOND, INDEPENDENT GATE ON THE RESIDUAL'S OWN  **
!**                         SIZE (MIN(1,|RESIDUAL_I|/RESREF)) -- EITHER     **
!**                         SIGNAL CAN OPEN THE GATE, SO A GENUINELY        **
!**                         SETTLED, SMALL-RESIDUAL BAND IS STILL LEFT      **
!**                         ALONE (THE ORIGINAL DESIGN INTENT), BUT A BAND  **
!**                         WITH A LARGE RESIDUAL NO LONGER GETS STUCK JUST **
!**                         BECAUSE ITS OWN (WRONG) EIGENVALUE HAPPENS TO   **
!**                         SIT AWAY FROM THE FERMI LEVEL.                  **
!**                                                                          **
!**  BY CONSTRUCTION, WITH HAUX INITIALIZED TO DIAG(EIG), EDFT$OCCUPATIONS   **
!**  MUST REPRODUCE DYNOCC_MERMIN(...,EIG,...) BIT FOR BIT (THE EIGENVALUES **
!**  OF A DIAGONAL MATRIX ARE ITS DIAGONAL ENTRIES) -- THE VALIDATION        **
!**  ANCHOR FOR EVERYTHING ELSE IN THIS FILE.                                **
!**                                                                          **
!**  ON EDFT$HAUXFORCE/STEP: THE PHYSICALLY CORRECT SELF-CONSISTENCY         **
!**  CONDITION FOR THE AUXILIARY HAMILTONIAN IS HAUX = HPHYS[RHO(HAUX)],     **
!**  WHERE HPHYS[RHO] IS THE REAL PHYSICAL HAMILTONIAN (BUILT FROM THE       **
!**  CURRENT DENSITY) PROJECTED INTO THE SAME BASIS THAT HAUX LIVES IN --    **
!**  EXACTLY THE MANUSCRIPT'S OWN SELF-CONSISTENCY RESIDUAL R=Q-Q[H(R,Q)]    **
!**  (SEC. "SPECIALIZATION TO SELF-CONSISTENT TIGHT BINDING"), PHRASED FOR   **
!**  A HAMILTONIAN-VALUED VARIABLE INSTEAD OF A CHARGE VECTOR. EDFT$HAUXFORCE**
!**  COMPUTES THAT RESIDUAL (HPHYS-HAUX, SUPPLIED BY THE CALLER, WHICH OWNS  **
!**  THE MACHINERY TO BUILD HPHYS -- NOT IMPLEMENTED HERE, SEE               **
!**  THE DEV NOTES FOR THE WIRING-IN PLAN) AND EDFT$HAUXSTEP APPLIES A   **
!**  SIMPLE DAMPED-MIXING STEP ALONG IT. THIS IS A STEPPING STONE TOWARD,    **
!**  NOT YET, THE FULL MARZARI-VANDERBILT-PAYNE GENERALIZED-FORCE            **
!**  EXPRESSION (WHICH ADDITIONALLY WEIGHTS BY THE ENTROPY-DERIVATIVE/       **
!**  OCCUPATION-NUMBER STRUCTURE): IT HAS THE CORRECT FIXED POINT AND IS     **
!**  SIMPLE ENOUGH TO VALIDATE NOW, BUT ITS CONVERGENCE PROPERTIES (RATE,    **
!**  STABILITY FOR A GENUINE SELF-CONSISTENT, NONLINEAR HPHYS[RHO(HAUX)])    **
!**  ARE NOT YET CHARACTERIZED -- ONLY THE LINEAR, FIXED-TARGET TOY PROBLEM  **
!**  (HPHYS HELD CONSTANT ACROSS ITERATIONS) HAS BEEN VALIDATED SO FAR.      **
!**  AN OCCUPATION-WEIGHTED ORBITAL GRADIENT IN VESPAW_CG.F90 IS ALSO NOT    **
!**  IMPLEMENTED YET.                                                        **
!*******************************************************************************
      MODULE EDFT_MODULE
      COMPLEX(8),ALLOCATABLE :: HAUX(:,:,:,:)   ! (NB,NB,NKPT,NSPIN)
      LOGICAL(4),ALLOCATABLE :: TBLOCKINI(:,:)  ! (NKPT,NSPIN) -- SEE EDFT$INIBLOCK
      INTEGER(4)             :: NB=0
      INTEGER(4)             :: NKPT=0
      INTEGER(4)             :: NSPIN=0
      LOGICAL(4)             :: TINI=.FALSE.
!     == LARGEST FERMI-WINDOW |HAUXEIG-EIGVAL| SEEN OVER THE OUTER ITERATION  ==
!     == JUST FINISHED (SEE EDFT$SETHAUXRESMAX/EDFT$GETHAUXRESMAX BELOW) --  ==
!     == INITIALIZED TO A LARGE SENTINEL SO A CALLER CHECKING IT BEFORE ANY  ==
!     == BLOCK HAS EVER REPORTED A REAL VALUE CORRECTLY SEES "NOT CONVERGED" ==
!     == RATHER THAN A STALE ZERO.                                           ==
      REAL(8)                :: HAUXRESMAX=1.D+20
      END MODULE EDFT_MODULE
!
!*******************************************************************************
!**  EDFTMIXANDERSON_MODULE -- PER-(K-POINT,SPIN)-BLOCK ANDERSON-MIXING STATE **
!**  FOR EDFT$HAUXSTEP_ANDERSON, BELOW. KEPT SEPARATE FROM EDFT_MODULE ITSELF **
!**  (SAME SEPARATION PRINCIPLE AS MIXANDERSON_MODULE VS. THE DENSITY-MIXING  **
!**  MACHINERY IT IS MODELED ON, VESPAW_MIXER.F90) SINCE THIS IS BOOKKEEPING  **
!**  FOR THE ITERATIVE MIXER'S OWN RESIDUAL/INPUT HISTORY, NOT HAUX ITSELF.   **
!**  ONE INDEPENDENT HISTORY PER (K-POINT,SPIN) BLOCK, UNLIKE DENSITY         **
!**  MIXING'S SINGLE GLOBAL HISTORY -- SEE EDFT$HAUXSTEP_ANDERSON'S OWN       **
!**  HEADER COMMENT FOR WHY.                                                  **
!*******************************************************************************
      MODULE EDFTMIXANDERSON_MODULE
      INTEGER(4),PARAMETER    :: NHISTORY=4
      REAL(8)   ,PARAMETER    :: BETAMIN=0.02D0,BETAMAX=0.5D0
      LOGICAL(4),ALLOCATABLE,SAVE :: TFIRST(:,:)     ! (NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE,SAVE :: BETA(:,:)       ! (NKPT,NSPIN) -- ADAPTIVE DAMPING
      REAL(8)   ,ALLOCATABLE,SAVE :: RESBEST(:,:)    ! (NKPT,NSPIN)
      INTEGER(4),ALLOCATABLE,SAVE :: NVALID(:,:)     ! (NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE,SAVE :: XOLD(:,:,:,:)   ! (NTOT,NHISTORY+1,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE,SAVE :: FOLD(:,:,:,:)   ! (NTOT,NHISTORY+1,NKPT,NSPIN)
      LOGICAL(4)             ,SAVE :: TALLOCATED=.FALSE.
      END MODULE EDFTMIXANDERSON_MODULE
!
!     ..................................................................
      SUBROUTINE EDFT$INI(NB_,NKPT_,NSPIN_,EIG)
!     ******************************************************************
!     **  INITIALIZE THE AUXILIARY HAMILTONIAN HAUX=DIAG(EIG), ONE     **
!     **  (NB,NB) HERMITIAN BLOCK PER (K-POINT,SPIN).                  **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB_
      INTEGER(4),INTENT(IN) :: NKPT_
      INTEGER(4),INTENT(IN) :: NSPIN_
      REAL(8)   ,INTENT(IN) :: EIG(NB_,NKPT_,NSPIN_)
      INTEGER(4)            :: IB,IKPT,ISPIN
!     ******************************************************************
      IF(ALLOCATED(HAUX)) DEALLOCATE(HAUX)
      NB=NB_
      NKPT=NKPT_
      NSPIN=NSPIN_
      ALLOCATE(HAUX(NB,NB,NKPT,NSPIN))
      HAUX(:,:,:,:)=(0.D0,0.D0)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NB
            HAUX(IB,IB,IKPT,ISPIN)=CMPLX(EIG(IB,IKPT,ISPIN),0.D0,KIND=8)
          ENDDO
        ENDDO
      ENDDO
      TINI=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$INIBLOCK(IKPT_,ISPIN_,NB_,NKPT_,NSPIN_,EIGBLOCK)
!     ******************************************************************
!     **  LAZY, PER-(K-POINT,SPIN)-BLOCK VARIANT OF EDFT$INI, FOR A     **
!     **  CALLER (E.G. DAVIDSON$STATE_BY_BLOCK) THAT VISITS ONE BLOCK   **
!     **  AT A TIME AND DOES NOT HAVE THE FULL (NB,NKPT,NSPIN) EIG      **
!     **  ARRAY ASSEMBLED AT ONCE. ALLOCATES HAUX (ALL BLOCKS, ZEROED)  **
!     **  ON THE FIRST EVER CALL, THEN SEEDS HAUX(:,:,IKPT_,ISPIN_)=    **
!     **  DIAG(EIGBLOCK) THE FIRST TIME *THAT* BLOCK IS VISITED ONLY --**
!     **  SUBSEQUENT CALLS FOR AN ALREADY-SEEDED BLOCK ARE A NO-OP, SO  **
!     **  THIS IS SAFE TO CALL ON EVERY OUTER ITERATION.                **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      INTEGER(4),INTENT(IN) :: NB_
      INTEGER(4),INTENT(IN) :: NKPT_
      INTEGER(4),INTENT(IN) :: NSPIN_
      REAL(8)   ,INTENT(IN) :: EIGBLOCK(NB_)
      INTEGER(4)            :: IB
!     ******************************************************************
      IF(.NOT.TINI) THEN
        NB=NB_
        NKPT=NKPT_
        NSPIN=NSPIN_
        ALLOCATE(HAUX(NB,NB,NKPT,NSPIN))
        ALLOCATE(TBLOCKINI(NKPT,NSPIN))
        HAUX(:,:,:,:)=(0.D0,0.D0)
        TBLOCKINI(:,:)=.FALSE.
        TINI=.TRUE.
      END IF
      IF(NB_.NE.NB.OR.NKPT_.NE.NKPT.OR.NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('DIMENSIONS INCONSISTENT WITH A PRIOR EDFT$INIBLOCK CALL')
        CALL ERROR$MSG('(E.G. NB VARIES ACROSS K-POINTS -- NOT SUPPORTED)')
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$I4VAL('NB (FROM EARLIER CALL)',NB)
        CALL ERROR$I4VAL('NKPT_',NKPT_)
        CALL ERROR$I4VAL('NKPT (FROM EARLIER CALL)',NKPT)
        CALL ERROR$I4VAL('NSPIN_',NSPIN_)
        CALL ERROR$I4VAL('NSPIN (FROM EARLIER CALL)',NSPIN)
        CALL ERROR$STOP('EDFT$INIBLOCK')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$INIBLOCK')
      END IF
      IF(.NOT.TBLOCKINI(IKPT_,ISPIN_)) THEN
        HAUX(:,:,IKPT_,ISPIN_)=(0.D0,0.D0)
        DO IB=1,NB
          HAUX(IB,IB,IKPT_,ISPIN_)=CMPLX(EIGBLOCK(IB),0.D0,KIND=8)
        ENDDO
        TBLOCKINI(IKPT_,ISPIN_)=.TRUE.
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$ISINI(TINI_)
!     ******************************************************************
!     **  RETURNS WHETHER HAUX HAS BEEN INITIALIZED YET (TRUE FROM THE  **
!     **  FIRST EVER EDFT$INI/EDFT$INIBLOCK CALL ONWARD -- SEE THEIR    **
!     **  OWN HEADER COMMENTS). LETS A CALLER (E.G. THE LIVE-OCCUPATION **
!     **  FEED IN VESPAW_WAVES1.F90) CHECK READINESS BEFORE CALLING     **
!     **  ANYTHING ELSE IN THIS MODULE, WITHOUT POKING TINI DIRECTLY.   **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(OUT):: TINI_
!     ******************************************************************
      TINI_=TINI
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$SETHAUXRESMAX(RES_)
!     ******************************************************************
!     **  RECORDS THE LARGEST FERMI-WINDOW |HAUXEIG-EIGVAL| SEEN OVER    **
!     **  ONE OUTER ITERATION (ACROSS ALL (K-POINT,SPIN) BLOCKS), FOR    **
!     **  THE HAUX-RESIDUAL-AWARE STOPPING CHECK IN VESPAW_WAVES1.F90 --**
!     **  SEE ITS OWN HEADER COMMENT AND THE DEV NOTES. CALLED FROM  **
!     **  DAVIDSON$STATE_BY_BLOCK (VESPAW_CG.F90) ONCE PER ITERATION,    **
!     **  AFTER THE PER-BLOCK LOOP, WHEN LIVE OCCUPATIONS ARE ACTIVE.    **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: RES_
!     ******************************************************************
      HAUXRESMAX=RES_
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$GETHAUXRESMAX(RES_)
!     ******************************************************************
!     **  RETURNS THE VALUE LAST SET BY EDFT$SETHAUXRESMAX (OR THE       **
!     **  MODULE'S LARGE SENTINEL DEFAULT IF NEVER SET) -- SEE THAT       **
!     **  ROUTINE'S HEADER COMMENT.                                       **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: RES_
!     ******************************************************************
      RES_=HAUXRESMAX
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXEIG(IKPT_,ISPIN_,NB_,E)
!     ******************************************************************
!     **  RETURNS THE CURRENT EIGENVALUES OF HAUX FOR A SINGLE           **
!     **  (K-POINT,SPIN) BLOCK, WITHOUT FERMI-FILLING. TWO USES:         **
!     **  (1) DIAGNOSTIC -- COMPARING HAUX'S OWN SPECTRUM AGAINST THE    **
!     **  REAL HAMILTONIAN (E.G. EIGVAL) DURING SHADOW-MODE VALIDATION;  **
!     **  (2) LIVE -- THE ACTUAL VALUES SUBSTITUTED FOR THE REAL,        **
!     **  INSTANTANEOUS <PSI|H|PSI> EXPECTATION VALUE BEFORE THEY REACH  **
!     **  DYNOCC'S OCCUPATION-DETERMINATION MACHINERY, WHEN              **
!     **  VESPAW_EDFT_LIVE_OCC IS SET (SEE VESPAW_WAVES1.F90, WAVES$ETOT,**
!     **  AND THE DEV NOTES). NOT ITSELF THE FERMI-FILL -- SEE       **
!     **  EDFT$OCCUPATIONS FOR THAT.                                     **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      INTEGER(4),INTENT(IN) :: NB_
      REAL(8)   ,INTENT(OUT):: E(NB_)
      COMPLEX(8)             :: U(NB_,NB_)
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXEIG CALLED BEFORE EDFT$INI/EDFT$INIBLOCK')
        CALL ERROR$STOP('EDFT$HAUXEIG')
      END IF
      IF(NB_.NE.NB.OR.IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('DIMENSIONS OR INDICES INCONSISTENT')
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXEIG')
      END IF
      CALL LIB$DIAGC8(NB,HAUX(:,:,IKPT_,ISPIN_),E,U)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXEIGVEC(IKPT_,ISPIN_,NB_,E,U)
!     ******************************************************************
!     **  LIKE EDFT$HAUXEIG, BUT ALSO RETURNS HAUX'S OWN EIGENVECTOR      **
!     **  MATRIX U(NB_,NB_) -- THE UNITARY "NATURAL ORBITAL" ROTATION     **
!     **  EDFT$STATE_BY_BLOCK (VESPAW_CG.F90) APPLIES TO THE TRIAL        **
!     **  ORBITALS BEFORE WRITING THEM BACK TO THIS%PSI0, PER THE REAL    **
!     **  MARZARI-VANDERBILT-PAYNE JOINT ORBITAL+OCCUPATION SCHEME (SEE   **
!     **  THE DEV NOTES): DIAGONALIZING HAUX GIVES BOTH THE FERMI-    **
!     **  FILLED OCCUPATIONS (EDFT$OCCUPATIONS) AND THE NATURAL-ORBITAL   **
!     **  BASIS THE TRIAL SUBSPACE MUST BE ROTATED INTO, SINCE ROTATING   **
!     **  BETWEEN TWO TRIAL ORBITALS WITH DIFFERENT OCCUPATION WEIGHT     **
!     **  CHANGES THE TOTAL ENERGY (UNLIKE THE INSULATOR CASE, WHERE ANY  **
!     **  ORTHONORMAL BASIS OF THE OCCUPIED SUBSPACE IS EQUIVALENT).      **
!     **  E MUST BE BIT-IDENTICAL TO EDFT$HAUXEIG'S OUTPUT ON THE SAME    **
!     **  HAUX STATE (SAME UNDERLYING LIB$DIAGC8 CALL) -- A CHEAP,        **
!     **  BUILT-IN CONSISTENCY CHECK BETWEEN THE TWO ACCESSORS.           **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      INTEGER(4),INTENT(IN) :: NB_
      REAL(8)   ,INTENT(OUT):: E(NB_)
      COMPLEX(8),INTENT(OUT):: U(NB_,NB_)
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXEIGVEC CALLED BEFORE EDFT$INI/EDFT$INIBLOCK')
        CALL ERROR$STOP('EDFT$HAUXEIGVEC')
      END IF
      IF(NB_.NE.NB.OR.IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('DIMENSIONS OR INDICES INCONSISTENT')
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXEIGVEC')
      END IF
      CALL LIB$DIAGC8(NB,HAUX(:,:,IKPT_,ISPIN_),E,U)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$GETHAUX(IKPT_,ISPIN_,NB_,H)
!     ******************************************************************
!     **  RETURNS ONE (K-POINT,SPIN) BLOCK OF HAUX. USED BY THE         **
!     **  VARIATIONAL EDFT STEP (EDFT$VARSTEP, VESPAW_CG.F90), WHICH    **
!     **  TREATS HAUX AS AN EXPLICIT MINIMIZATION VARIABLE.             **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      INTEGER(4),INTENT(IN) :: NB_
      COMPLEX(8),INTENT(OUT):: H(NB_,NB_)
!     ******************************************************************
      IF(.NOT.TINI.OR.NB_.NE.NB) THEN
        CALL ERROR$MSG('HAUX NOT INITIALIZED OR NB INCONSISTENT')
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$STOP('EDFT$GETHAUX')
      END IF
      H=HAUX(:,:,IKPT_,ISPIN_)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$SETHAUX(IKPT_,ISPIN_,NB_,H)
!     ******************************************************************
!     **  OVERWRITES ONE (K-POINT,SPIN) BLOCK OF HAUX (HERMITIANIZED).  **
!     **  COUNTERPART OF EDFT$GETHAUX.                                  **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      INTEGER(4),INTENT(IN) :: NB_
      COMPLEX(8),INTENT(IN) :: H(NB_,NB_)
!     ******************************************************************
      IF(.NOT.TINI.OR.NB_.NE.NB) THEN
        CALL ERROR$MSG('HAUX NOT INITIALIZED OR NB INCONSISTENT')
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$STOP('EDFT$SETHAUX')
      END IF
      HAUX(:,:,IKPT_,ISPIN_)=0.5D0*(H+CONJG(TRANSPOSE(H)))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$OCCUPATIONS(NB_,NKPT_,NSPIN_,TOTCHA,ISPINDEG,TEMP &
     &                            ,WKPT,F,CHMPOT,EMERMN)
!     ******************************************************************
!     **  DIAGONALIZES THE CURRENT AUXILIARY HAMILTONIAN HAUX AND      **
!     **  FERMI-FILLS ITS EIGENVALUES -- I.E. THE SAME BISECTION-ON-   **
!     **  CHEMICAL-POTENTIAL PROCEDURE AS DYNOCC_MERMIN, APPLIED TO    **
!     **  HAUX'S EIGENVALUES INSTEAD OF THE PHYSICAL HAMILTONIAN'S.    **
!     **  WITH HAUX=DIAG(EIG) (SEE EDFT$INI) THIS REPRODUCES           **
!     **  DYNOCC_MERMIN(...,EIG,...) EXACTLY, SINCE THE EIGENVALUES OF **
!     **  A DIAGONAL MATRIX ARE ITS DIAGONAL ENTRIES.                  **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB_
      INTEGER(4),INTENT(IN) :: NKPT_
      INTEGER(4),INTENT(IN) :: NSPIN_
      REAL(8)   ,INTENT(IN) :: TOTCHA
      INTEGER(4),INTENT(IN) :: ISPINDEG
      REAL(8)   ,INTENT(IN) :: TEMP
      REAL(8)   ,INTENT(IN) :: WKPT(NKPT_)
      REAL(8)   ,INTENT(OUT):: F(NB_,NKPT_,NSPIN_)
      REAL(8)   ,INTENT(OUT):: CHMPOT
      REAL(8)   ,INTENT(OUT):: EMERMN
      REAL(8)                :: EPSTILDE(NB_,NKPT_,NSPIN_)
      REAL(8)                :: E(NB_)
      COMPLEX(8)             :: U(NB_,NB_)
      INTEGER(4)             :: IKPT,ISPIN
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$OCCUPATIONS CALLED BEFORE EDFT$INI')
        CALL ERROR$STOP('EDFT$OCCUPATIONS')
      END IF
      IF(NB_.NE.NB.OR.NKPT_.NE.NKPT.OR.NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('DIMENSIONS INCONSISTENT WITH EDFT$INI')
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$I4VAL('NB (FROM EDFT$INI)',NB)
        CALL ERROR$I4VAL('NKPT_',NKPT_)
        CALL ERROR$I4VAL('NKPT (FROM EDFT$INI)',NKPT)
        CALL ERROR$I4VAL('NSPIN_',NSPIN_)
        CALL ERROR$I4VAL('NSPIN (FROM EDFT$INI)',NSPIN)
        CALL ERROR$STOP('EDFT$OCCUPATIONS')
      END IF
!
!     ==========================================================================
!     ==  DIAGONALIZE THE AUXILIARY HAMILTONIAN, ONE BLOCK PER (K-POINT,SPIN) ==
!     ==========================================================================
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL LIB$DIAGC8(NB,HAUX(:,:,IKPT,ISPIN),E,U)
          EPSTILDE(:,IKPT,ISPIN)=E(:)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  FERMI-FILL EXACTLY AS DYNOCC_MERMIN DOES, USING THE AUXILIARY       ==
!     ==  EIGENVALUES IN PLACE OF THE PHYSICAL HAMILTONIAN'S EIGENVALUES      ==
!     ==========================================================================
      CALL DYNOCC_MERMIN(NB,NKPT,NSPIN,TOTCHA,ISPINDEG,TEMP,WKPT &
     &                   ,EPSTILDE,F,CHMPOT,EMERMN)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXFORCE(IKPT_,ISPIN_,HPHYS,RESIDUAL)
!     ******************************************************************
!     **  GENERALIZED FORCE ON THE AUXILIARY HAMILTONIAN, FOR A SINGLE **
!     **  (K-POINT,SPIN) BLOCK: THE MISMATCH BETWEEN THE REAL           **
!     **  HAMILTONIAN HPHYS (PROJECTED INTO THE SAME BASIS AS HAUX,     **
!     **  SUPPLIED BY THE CALLER) AND HAUX ITSELF. VANISHES EXACTLY AT  **
!     **  THE SELF-CONSISTENT FIXED POINT HAUX=HPHYS[RHO(HAUX)]. SEE    **
!     **  THE FILE HEADER FOR THE FULL JUSTIFICATION.                   **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      COMPLEX(8),INTENT(IN) :: HPHYS(NB,NB)
      COMPLEX(8),INTENT(OUT):: RESIDUAL(NB,NB)
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXFORCE CALLED BEFORE EDFT$INI')
        CALL ERROR$STOP('EDFT$HAUXFORCE')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXFORCE')
      END IF
      RESIDUAL(:,:)=HPHYS(:,:)-HAUX(:,:,IKPT_,ISPIN_)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXSTEP(LAMBDA,IKPT_,ISPIN_,RESIDUAL)
!     ******************************************************************
!     **  DAMPED LINEAR-MIXING STEP OF HAUX ALONG THE GENERALIZED      **
!     **  FORCE FROM EDFT$HAUXFORCE, FOR A SINGLE (K-POINT,SPIN)        **
!     **  BLOCK:  HAUX <- HAUX + LAMBDA*RESIDUAL.  RE-HERMITIANIZES     **
!     **  AFTERWARD TO GUARD AGAINST ROUNDOFF ASYMMETRY ACCUMULATING    **
!     **  OVER MANY STEPS.                                              **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: LAMBDA
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      COMPLEX(8),INTENT(IN) :: RESIDUAL(NB,NB)
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXSTEP CALLED BEFORE EDFT$INI')
        CALL ERROR$STOP('EDFT$HAUXSTEP')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXSTEP')
      END IF
      HAUX(:,:,IKPT_,ISPIN_)=HAUX(:,:,IKPT_,ISPIN_)+CMPLX(LAMBDA,0.D0,KIND=8)*RESIDUAL(:,:)
      HAUX(:,:,IKPT_,ISPIN_)=0.5D0*(HAUX(:,:,IKPT_,ISPIN_) &
     &                             +TRANSPOSE(CONJG(HAUX(:,:,IKPT_,ISPIN_))))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXSTEP_MVP(LAMBDA,IKPT_,ISPIN_,CHMPOT,TEMP,RESIDUAL)
!     ******************************************************************
!     **  ENTROPY-DERIVATIVE-WEIGHTED VERSION OF EDFT$HAUXSTEP -- THE   **
!     **  PIECE OF THE REAL MARZARI-VANDERBILT-PAYNE GENERALIZED FORCE  **
!     **  THAT THE UNIFORM-WEIGHT EDFT$HAUXSTEP OMITS.                  **
!     **                                                                **
!     **  DERIVATION SKETCH (SEE THE DEV NOTES FOR THE FULL         **
!     **  ARGUMENT): FOR THE MERMIN FREE ENERGY A_T=SUM_I F_I*EPS_I-T*S,**
!     **  WITH HAUX'S OWN EIGENVALUES LAMBDA_I GENERATING THE           **
!     **  OCCUPATIONS F_I=FERMI(LAMBDA_I) AND EPS_I=<PSI_I|H|PSI_I> THE  **
!     **  REAL HAMILTONIAN'S EXPECTATION VALUE, THE GRADIENT OF A_T     **
!     **  W.R.T. LAMBDA_I (IN HAUX'S OWN EIGENBASIS) IS                 **
!     **  (DF_I/DLAMBDA_I)*(EPS_I-LAMBDA_I) -- THE SAME RESIDUAL         **
!     **  EDFT$HAUXFORCE ALREADY COMPUTES, BUT WEIGHTED BY THE          **
!     **  FERMI-DIRAC DERIVATIVE DF_I/DLAMBDA_I=-F_I*(1-F_I)/T. THIS     **
!     **  WEIGHT VANISHES FOR FULLY OCCUPIED/EMPTY STATES (F NEAR 0 OR  **
!     **  1) AND PEAKS FOR STATES NEAR THE FERMI LEVEL (F NEAR 0.5) --  **
!     **  I.E. IT AUTOMATICALLY LEAVES SETTLED STATES ALONE AND         **
!     **  CONCENTRATES THE UPDATE ON STATES WHOSE OCCUPATION IS         **
!     **  GENUINELY STILL BEING DETERMINED. THE WEIGHT ACTUALLY USED IS **
!     **  THIS DERIVATIVE'S *SHAPE*, RENORMALIZED TO PEAK AT 1.0 (NOT   **
!     **  THE RAW 1/(4*TEMP) PEAK VALUE, WHICH IS >>1 IN HARTREE UNITS  **
!     **  AT ORDINARY SMEARING TEMPERATURES AND MAKES THE STEP DIVERGE  **
!     **  -- SEE THE WEIGHT-COMPUTATION LOOP BELOW FOR THE FULL STORY). **
!     **  THIS SHAPE-ONLY WEIGHT IS THE MECHANISM                       **
!     **  MISSING FROM VESPAW_CG.F90'S OCCUPATION-WEIGHTED LOCKING      **
!     **  TOLERANCE (WHICH USES A CRUDE, UNIFORM-IN-ENERGY PROXY        **
!     **  INSTEAD), AND A CANDIDATE EXPLANATION FOR WHY THAT MECHANISM  **
!     **  DID NOT GENERALIZE FROM SI TO A REAL METAL.                   **
!     **                                                                **
!     **  OFF-DIAGONAL ELEMENTS OF THE RESIDUAL (COUPLING BETWEEN TWO   **
!     **  DIFFERENT HAUX EIGENSTATES I,J) ARE WEIGHTED BY THE GEOMETRIC  **
!     **  MEAN SQRT(W_I*W_J), WHICH REDUCES TO W_I ON THE DIAGONAL --   **
!     **  A SMOOTH, WELL-DEFINED CHOICE, NOT A LITERAL REPRODUCTION OF  **
!     **  MVP'S OWN DEGENERATE-PERTURBATION-THEORY OFF-DIAGONAL FORMULA.**
!     **  ONLY THE DIAGONAL PART IS DERIVED HERE WITH THE SAME RIGOR AS **
!     **  THE MANUSCRIPT'S OWN LEMMA 1 (GRAND-POTENTIAL GRADIENT).      **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: LAMBDA
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      REAL(8)   ,INTENT(IN) :: CHMPOT
      REAL(8)   ,INTENT(IN) :: TEMP
      COMPLEX(8),INTENT(IN) :: RESIDUAL(NB,NB)
      REAL(8)                :: EIGVAL(NB)
      COMPLEX(8)             :: U(NB,NB)
      COMPLEX(8)             :: RESEIG(NB,NB)
      REAL(8)                :: WEIGHT(NB)
      REAL(8)                :: SVAR,EX
      INTEGER(4)             :: I,J
      LOGICAL(4),SAVE        :: TMVPFULLOFFDIAGCHECKED=.FALSE.
      LOGICAL(4),SAVE        :: TMVPFULLOFFDIAG=.FALSE.
      CHARACTER(16)          :: CENVVALMVP
!     ******************************************************************
      IF(.NOT.TMVPFULLOFFDIAGCHECKED) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_EDFT_MVP_FULL_OFFDIAG',CENVVALMVP)
        TMVPFULLOFFDIAG=(LEN_TRIM(CENVVALMVP).GT.0)
        TMVPFULLOFFDIAGCHECKED=.TRUE.
        IF(TMVPFULLOFFDIAG) THEN
          PRINT*,'VESPAW_EDFT_MVP_FULL_OFFDIAG SET -- USING MAX(W_I,W_J)' &
     &          ,' INSTEAD OF SQRT(W_I*W_J) FOR OFF-DIAGONAL SUPPRESSION'
        END IF
      END IF
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXSTEP_MVP CALLED BEFORE EDFT$INI/INIBLOCK')
        CALL ERROR$STOP('EDFT$HAUXSTEP_MVP')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXSTEP_MVP')
      END IF
      IF(TEMP.LE.0.D0) THEN
        CALL ERROR$MSG('TEMP MUST BE POSITIVE')
        CALL ERROR$R8VAL('TEMP',TEMP)
        CALL ERROR$STOP('EDFT$HAUXSTEP_MVP')
      END IF
!
!     ==========================================================================
!     ==  DIAGONALIZE THE CURRENT HAUX TO GET ITS OWN EIGENBASIS             ==
!     ==========================================================================
      CALL LIB$DIAGC8(NB,HAUX(:,:,IKPT_,ISPIN_),EIGVAL,U)
!
!     ==========================================================================
!     ==  FERMI-DIRAC-DERIVATIVE *SHAPE*, NORMALIZED TO PEAK AT 1.0 AT THE    ==
!     ==  FERMI LEVEL: W_I=4*F_I*(1-F_I) (=4*EXP(X)/(1+EXP(X))**2, X=         ==
!     ==  (EIGVAL-CHMPOT)/TEMP). TEMP STILL SETS THE *WIDTH* OF THE ACTIVE    ==
!     ==  WINDOW THROUGH X, AS INTENDED -- BUT THE RAW DERIVATIVE             ==
!     ==  F_I*(1-F_I)/TEMP (WITHOUT THE 4*TEMP RESCALING) IS *NOT* USED       ==
!     ==  HERE, BECAUSE ITS PEAK VALUE 1/(4*TEMP) DIVERGES AS TEMP->0 AND IS  ==
!     ==  ALREADY >>1 IN HARTREE UNITS AT ORDINARY SMEARING TEMPERATURES      ==
!     ==  (E.G. ~79 AT T=1000K) -- MULTIPLYING AN ENERGY-SCALE RESIDUAL BY   ==
!     ==  THAT AND THEN TAKING A LAMBDA~0.3 STEP OF IT IS WILDLY UNSTABLE.   ==
!     ==  FOUND 2026-09-05: THE UNNORMALIZED VERSION MADE HAUX'S OWN         ==
!     ==  SELF-CONSISTENCY DIAGNOSTIC DIVERGE (MAX|HAUXEIG-EIGVAL| RUNNING   ==
!     ==  AWAY TO ~0.24 H ON A REAL METAL TEST) INSTEAD OF CONVERGING --      ==
!     ==  HARMLESS TO ANY REAL RUN SINCE THIS IS STILL SHADOW-MODE-ONLY, BUT ==
!     ==  A REAL BUG IN THE STEP ITSELF, NOT JUST THE DIAGNOSTIC. THE        ==
!     ==  NORMALIZED VERSION KEEPS THE SAME SHAPE (VANISHES AT F=0/1, PEAKS  ==
!     ==  AT F=0.5) WITHOUT THE DIVERGENT ABSOLUTE SCALE.                    ==
!     ==========================================================================
      DO I=1,NB
        SVAR=(EIGVAL(I)-CHMPOT)/TEMP
        IF(ABS(SVAR).GT.50.D0) THEN
          WEIGHT(I)=0.D0
        ELSE
          EX=EXP(SVAR)
          WEIGHT(I)=4.D0*EX/((1.D0+EX)**2)
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  TRANSFORM RESIDUAL INTO HAUX'S EIGENBASIS, WEIGHT, TRANSFORM BACK.  ==
!     ==  OFF-DIAGONAL (I.NE.J) COUPLING BETWEEN TWO HAUX EIGENSTATES: THE   ==
!     ==  GEOMETRIC MEAN SQRT(W_I*W_J) SUPPRESSES ANY COUPLING INVOLVING AN  ==
!     ==  INACTIVE STATE (W~0) EVEN WHEN THE OTHER STATE IS ACTIVE -- BUT    ==
!     ==  THAT COUPLING IS PRECISELY WHAT ROTATES HAUX'S EIGENBASIS TOWARD   ==
!     ==  THE PHYSICAL HAMILTONIAN'S TRUE NEAR-FERMI EIGENSTATES. DIAGNOSTIC ==
!     ==  ENV VAR VESPAW_EDFT_MVP_FULL_OFFDIAG (SEE FIRST CALL, CHECKED      ==
!     ==  ONCE) SWITCHES TO MAX(W_I,W_J) INSTEAD -- STILL SUPPRESSES         ==
!     ==  COUPLING BETWEEN TWO MUTUALLY INACTIVE STATES, BUT NO LONGER       ==
!     ==  SUPPRESSES COUPLING TOUCHING AN ACTIVE ONE -- TO TEST WHETHER THIS ==
!     ==  IS WHY THE GEOMETRIC-MEAN VERSION TRACKS WORSE THAN THE UNIFORM    ==
!     ==  STEP NEAR THE FERMI LEVEL (FOUND 2026-09-05, SEE DESIGN).         ==
!     ==========================================================================
      RESEIG=MATMUL(TRANSPOSE(CONJG(U)),MATMUL(RESIDUAL,U))
      IF(TMVPFULLOFFDIAG) THEN
        DO I=1,NB
          DO J=1,NB
            RESEIG(I,J)=RESEIG(I,J)*MAX(WEIGHT(I),WEIGHT(J))
          ENDDO
        ENDDO
      ELSE
        DO I=1,NB
          DO J=1,NB
            RESEIG(I,J)=RESEIG(I,J)*SQRT(MAX(WEIGHT(I)*WEIGHT(J),0.D0))
          ENDDO
        ENDDO
      END IF
      RESEIG=MATMUL(U,MATMUL(RESEIG,TRANSPOSE(CONJG(U))))
!
!     ==========================================================================
!     ==  DAMPED STEP, RE-HERMITIANIZE (SAME AS EDFT$HAUXSTEP)                ==
!     ==========================================================================
      HAUX(:,:,IKPT_,ISPIN_)=HAUX(:,:,IKPT_,ISPIN_)+CMPLX(LAMBDA,0.D0,KIND=8)*RESEIG(:,:)
      HAUX(:,:,IKPT_,ISPIN_)=0.5D0*(HAUX(:,:,IKPT_,ISPIN_) &
     &                             +TRANSPOSE(CONJG(HAUX(:,:,IKPT_,ISPIN_))))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXSTEP_ADAPTIVE(LAMBDA,IKPT_,ISPIN_,CHMPOT,TEMP,RESIDUAL)
!     ******************************************************************
!     **  PER-BAND-ADAPTIVE VERSION OF EDFT$HAUXSTEP_MVP -- FIXES THE   **
!     **  MECHANISM FOUND 2026-09-05 (SEE THE DEV NOTES) BY WHICH   **
!     **  EDFT$HAUXSTEP_MVP TRACKED FERMI-ACTIVE BANDS ~110X *WORSE*    **
!     **  THAN THE TRIVIAL UNIFORM STEP ON A REAL METAL TEST.           **
!     **                                                                **
!     **  ROOT CAUSE OF THAT REGRESSION: THE FERMI-DIRAC-DERIVATIVE     **
!     **  WEIGHT (NORMALIZED TO PEAK AT 1.0 AT F=0.5, SEE                **
!     **  EDFT$HAUXSTEP_MVP) IS STRICTLY <1 FOR *ANY* BAND NOT EXACTLY  **
!     **  AT HALF OCCUPATION -- WHICH IS MOST "ACTIVE" BANDS. MULTIPLYING**
!     **  THE RESIDUAL BY THIS WEIGHT AND THEN TAKING A *FIXED* LAMBDA  **
!     **  STEP (TUNED FOR THE OLD, UNWEIGHTED SCHEME) SYSTEMATICALLY     **
!     **  SHRINKS THE EFFECTIVE PER-ITERATION CORRECTION FOR EVERY       **
!     **  WEIGHTED BAND -- A STRUCTURAL MISMATCH BETWEEN A GRADIENT      **
!     **  WEIGHT (MEANT FOR A PROPER ADAPTIVE-STEP OPTIMIZER) AND A      **
!     **  FIXED-STEP DAMPED-MIXING SCHEME, CONFIRMED NOT TO BE AN        **
!     **  OFF-DIAGONAL-FORMULA ARTIFACT (VESPAW_EDFT_MVP_FULL_OFFDIAG    **
!     **  GAVE NUMERICALLY IDENTICAL RESULTS).                           **
!     **                                                                **
!     **  FIX: REPLACE THE WEIGHT ITSELF WITH A *GATE*,                  **
!     **  GATE_I=MIN(1,W_I/WREF) -- STILL 0 FOR SETTLED (F NEAR 0/1)     **
!     **  STATES, BUT SATURATES AT EXACTLY 1.0 (FULL LAMBDA, SAME SPEED  **
!     **  AS THE UNIFORM STEP) FOR ANY BAND WHOSE OCCUPATION IS          **
!     **  MEANINGFULLY FRACTIONAL (W_I ABOVE THE SMALL REFERENCE WREF),  **
!     **  INSTEAD OF SCALING DOWN EVERY ACTIVE BAND IN PROPORTION TO HOW **
!     **  CLOSE IT HAPPENS TO SIT TO EXACTLY F=0.5. THE GATE IS A        **
!     **  CONTINUOUS (CLAMPED-LINEAR) FUNCTION OF W_I -- NO HARD         **
!     **  DISCONTINUITY, THE SAME "SMOOTH RAMP, NOT A HARD CLAMP"        **
!     **  PRINCIPLE AS A3'S FIX TO THE LOCKING-TOLERANCE MECHANISM IN     **
!     **  VESPAW_CG.F90.                                                 **
!     **                                                                **
!     **  WREF IS RUNTIME-CONFIGURABLE VIA VESPAW_EDFT_ADAPTIVE_WREF     **
!     **  (DEFAULT 1.D-6, EMPIRICALLY TUNED 2026-09-05 ON THE REAL AL    **
!     **  METAL TEST -- SEE THE DEV NOTES. AN INITIAL GUESS OF 4.D-2 **
!     **  (GATE=1 ONCE F*(1-F)>1.D-2) STILL LEFT MANY GENUINELY ACTIVE   **
!     **  BANDS BELOW THE GATE'S SATURATION POINT (A BAND CAN SIT AT,    **
!     **  SAY, F=0.999996 AND STILL COUNT AS "ACTIVE" BY THE FERMI-      **
!     **  WINDOW DIAGNOSTIC'S F*(1-F)>1.D-4 CRITERION, YET HAVE A RAW    **
!     **  WEIGHT AS SMALL AS ~3.D-5); 1.D-6 (GATE=1 WITHIN ROUGHLY        **
!     **  15*TEMP OF THE FERMI LEVEL) WAS THE VALUE THAT ACTUALLY CLOSED **
!     **  THE GAP -- THE FERMI-WINDOW-RESTRICTED TRACKING RESIDUAL AT    **
!     **  THIS SETTING MATCHED THE PLAIN UNIFORM STEP'S (MEAN 3.76.D-5   **
!     **  VS. 3.80.D-5 OVER THE ACTIVE BLOCKS AT THE FINAL ITERATION),   **
!     **  VS. 4.D-2'S MEAN 1.66.D-3 (STILL ~44X WORSE THAN UNIFORM) AND  **
!     **  THE UNGATED RAW WEIGHT'S 2.93.D-3 (~77X WORSE). GOING SMALLER  **
!     **  THAN THIS EVENTUALLY DEFEATS THE POINT (AN ASTRONOMICALLY      **
!     **  TINY WREF SATURATES GATE=1 EVEN FOR GENUINELY SETTLED BANDS,   **
!     **  RECOVERING THE PLAIN UNIFORM STEP FOR EVERYTHING) -- 1.D-6 IS  **
!     **  THE VALUE ACTUALLY MEASURED TO WORK, NOT A THEORETICAL LIMIT.  **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: LAMBDA
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      REAL(8)   ,INTENT(IN) :: CHMPOT
      REAL(8)   ,INTENT(IN) :: TEMP
      COMPLEX(8),INTENT(IN) :: RESIDUAL(NB,NB)
      REAL(8)                :: EIGVAL(NB)
      COMPLEX(8)             :: U(NB,NB)
      COMPLEX(8)             :: RESEIG(NB,NB)
      REAL(8)                :: WEIGHT(NB)
      REAL(8)                :: GATE(NB)
      REAL(8)                :: SVAR,EX
      INTEGER(4)             :: I,J
      LOGICAL(4),SAVE        :: TWREFCHECKED=.FALSE.
      REAL(8)   ,SAVE        :: WREF=1.D-6
      CHARACTER(32)          :: CENVVALWREF
!     ******************************************************************
      IF(.NOT.TWREFCHECKED) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_EDFT_ADAPTIVE_WREF',CENVVALWREF)
        IF(LEN_TRIM(CENVVALWREF).GT.0) THEN
          READ(CENVVALWREF,*) WREF
          PRINT*,'VESPAW_EDFT_ADAPTIVE_WREF SET -- USING WREF=',WREF
        END IF
        TWREFCHECKED=.TRUE.
      END IF
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXSTEP_ADAPTIVE CALLED BEFORE EDFT$INI/INIBLOCK')
        CALL ERROR$STOP('EDFT$HAUXSTEP_ADAPTIVE')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXSTEP_ADAPTIVE')
      END IF
      IF(TEMP.LE.0.D0) THEN
        CALL ERROR$MSG('TEMP MUST BE POSITIVE')
        CALL ERROR$R8VAL('TEMP',TEMP)
        CALL ERROR$STOP('EDFT$HAUXSTEP_ADAPTIVE')
      END IF
!
!     ==========================================================================
!     ==  DIAGONALIZE THE CURRENT HAUX TO GET ITS OWN EIGENBASIS             ==
!     ==========================================================================
      CALL LIB$DIAGC8(NB,HAUX(:,:,IKPT_,ISPIN_),EIGVAL,U)
!
!     ==========================================================================
!     ==  SAME NORMALIZED FERMI-DIRAC-DERIVATIVE SHAPE AS EDFT$HAUXSTEP_MVP, ==
!     ==  W_I=4*F_I*(1-F_I), THEN GATED: GATE_I=MIN(1,W_I/WREF).             ==
!     ==========================================================================
      DO I=1,NB
        SVAR=(EIGVAL(I)-CHMPOT)/TEMP
        IF(ABS(SVAR).GT.50.D0) THEN
          WEIGHT(I)=0.D0
        ELSE
          EX=EXP(SVAR)
          WEIGHT(I)=4.D0*EX/((1.D0+EX)**2)
        END IF
        GATE(I)=MIN(1.D0,WEIGHT(I)/WREF)
      ENDDO
!
!     ==========================================================================
!     ==  TRANSFORM RESIDUAL INTO HAUX'S EIGENBASIS, GATE, TRANSFORM BACK.    ==
!     ==  OFF-DIAGONAL: GEOMETRIC MEAN SQRT(GATE_I*GATE_J) -- SHOWN NOT TO    ==
!     ==  MATTER EMPIRICALLY VS. MAX(GATE_I,GATE_J) FOR THE RAW-WEIGHT       ==
!     ==  VERSION, KEPT FOR SIMPLICITY/CONSISTENCY WITH EDFT$HAUXSTEP_MVP.   ==
!     ==========================================================================
      RESEIG=MATMUL(TRANSPOSE(CONJG(U)),MATMUL(RESIDUAL,U))
      DO I=1,NB
        DO J=1,NB
          RESEIG(I,J)=RESEIG(I,J)*SQRT(MAX(GATE(I)*GATE(J),0.D0))
        ENDDO
      ENDDO
      RESEIG=MATMUL(U,MATMUL(RESEIG,TRANSPOSE(CONJG(U))))
!
!     ==========================================================================
!     ==  DAMPED STEP, RE-HERMITIANIZE (SAME AS EDFT$HAUXSTEP)                ==
!     ==========================================================================
      HAUX(:,:,IKPT_,ISPIN_)=HAUX(:,:,IKPT_,ISPIN_)+CMPLX(LAMBDA,0.D0,KIND=8)*RESEIG(:,:)
      HAUX(:,:,IKPT_,ISPIN_)=0.5D0*(HAUX(:,:,IKPT_,ISPIN_) &
     &                             +TRANSPOSE(CONJG(HAUX(:,:,IKPT_,ISPIN_))))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXSTEP_MVPEXACT(LAMBDA,IKPT_,ISPIN_,CHMPOT,TEMP,RESIDUAL)
!     ******************************************************************
!     **  THE ACTUAL, VERIFIED dF/dHAUX -- NOT A HEURISTIC WEIGHTED       **
!     **  RESIDUAL LIKE EVERY EDFT$HAUXSTEP* ROUTINE ABOVE.               **
!     **                                                                  **
!     **  FOUND 2026-09-05 (SEE THE DEV NOTES AND TESTS/UNITTESTS/     **
!     **  EDFT_GRADIENT_CHECK): AN INDEPENDENT FINITE-DIFFERENCE CHECK OF  **
!     **  EDFT$HAUXSTEP_ADAPTIVE'S RESIDUAL AGAINST THE TRUE GRADIENT OF   **
!     **  THE MERMIN FREE ENERGY F(HAUX)=TR[RHO(HAUX)*HPHYS]-TEMP*S(F)     **
!     **  SHOWED A REAL, dx-INDEPENDENT 28-107% MISMATCH -- NOT TRUNCATION **
!     **  NOISE -- BOTH NEAR AND FAR FROM SELF-CONSISTENCY. ROOT CAUSE,    **
!     **  FOUND BY RE-DERIVING dF/dHAUX FROM FULL FIRST-ORDER              **
!     **  PERTURBATION THEORY OF A HERMITIAN EIGENPROBLEM: THE TRUE        **
!     **  GRADIENT HAS TWO PIECES NO EARLIER EDFT$HAUXSTEP* ROUTINE HAD.   **
!     **                                                                  **
!     **  WITH U,LAMBDA_I = HAUX'S OWN EIGENVECTORS/EIGENVALUES, M=U^H*    **
!     **  HPHYS*U (HPHYS IN HAUX'S EIGENBASIS -- COMPUTED HERE WITHOUT A   **
!     **  SEPARATE HPHYS ARGUMENT, SEE THE "NO SEPARATE HPHYS NEEDED"      **
!     **  NOTE BELOW), F_I THE FERMI OCCUPATIONS AND F'_I=-F_I*(1-F_I)/    **
!     **  TEMP THEIR LAMBDA_I-DERIVATIVE:                                 **
!     **                                                                  **
!     **   DIAGONAL: dF/dHAUX_II = F'_I*[(M_II-LAMBDA_I) - <M-LAMBDA>_F'] **
!     **     WHERE <X>_F' = SUM_K F'_K*X_K / SUM_K F'_K IS A FIXED-        **
!     **     ELECTRON-COUNT LAGRANGE-MULTIPLIER CORRECTION (IT ENFORCES   **
!     **     SUM_I DELTA_F_I=0 UNDER ANY DELTA_HAUX, I.E. THAT THE STEP    **
!     **     DIRECTION ITSELF NEVER TRIES TO CHANGE THE TOTAL ELECTRON     **
!     **     COUNT) -- EVERY EARLIER ROUTINE'S RAW (M_II-LAMBDA_I)         **
!     **     RESIDUAL HAS NO ANALOGUE OF THIS TERM AT ALL.                 **
!     **                                                                  **
!     **   OFF-DIAGONAL (K.NE.L): dF/dHAUX_KL = (F_K-F_L)/(LAMBDA_K-       **
!     **     LAMBDA_L)*M_KL -- MVP'S OWN DEGENERATE-PERTURBATION-THEORY    **
!     **     FORMULA FOR THE EIGENVECTOR-ROTATION CONTRIBUTION (SMOOTHLY   **
!     **     -> F'_AVG*M_KL AS LAMBDA_K->LAMBDA_L, SEE THE DEGENERACY       **
!     **     FALLBACK BELOW). EDFT$HAUXSTEP_MVP/_ADAPTIVE INSTEAD USE       **
!     **     SQRT(WEIGHT_K*WEIGHT_L) -- A GEOMETRIC-MEAN SMOOTHING THEIR   **
!     **     OWN HEADER COMMENTS ALREADY FLAGGED AS "NOT A LITERAL          **
!     **     REPRODUCTION OF MVP'S OWN OFF-DIAGONAL FORMULA" -- A KNOWN,    **
!     **     DELIBERATE SIMPLIFICATION THAT TURNS OUT TO MATTER             **
!     **     QUANTITATIVELY, NOT JUST COSMETICALLY.                        **
!     **                                                                  **
!     **  NO SEPARATE HPHYS ARGUMENT NEEDED: SINCE RESIDUAL=HPHYS-HAUX     **
!     **  (EDFT$HAUXFORCE) AND U^H*HAUX*U=DIAG(LAMBDA) EXACTLY (U IS       **
!     **  HAUX'S OWN EIGENBASIS), U^H*HPHYS*U-DIAG(LAMBDA)=U^H*RESIDUAL*U  **
!     **  EXACTLY -- SO ROTATING RESIDUAL INTO HAUX'S EIGENBASIS ALREADY   **
!     **  GIVES BOTH (M_II-LAMBDA_I) (THE DIAGONAL) AND M_KL, K.NE.L (THE  **
!     **  OFF-DIAGONAL) DIRECTLY, WITH THE SAME SIGNATURE AS EVERY OTHER   **
!     **  EDFT$HAUXSTEP* ROUTINE.                                          **
!     **                                                                  **
!     **  SIGN CONVENTION: THIS COMPUTES +dF/dHAUX DIRECTLY (VALIDATED     **
!     **  AGAINST +FINITE-DIFFERENCE, NO SIGN ASSUMED -- SEE THE TEST      **
!     **  FILE), SO THE STEP APPLIED IS HAUX -= LAMBDA*GRADLAB (DESCENT),   **
!     **  MATCHING THE SIGN CONVENTION OF RESEIG IN EVERY EARLIER STEP      **
!     **  ROUTINE IN THIS FILE (HAUX += LAMBDA*RESEIG, WHERE RESEIG POINTS  **
!     **  ALONG -dF/dHAUX BY VALIDATED CONSTRUCTION).                      **
!     **                                                                  **
!     **  VALIDATED (TESTS/UNITTESTS/EDFT_GRADIENT_CHECK, 2026-09-05):      **
!     **  MATCHES THE FINITE-DIFFERENCE GRADIENT TO 0.00-0.22% AT THE       **
!     **  SMALLEST DX, FOR EVERY DIAGONAL INDEX, IN BOTH A NEAR- AND A      **
!     **  FAR-FROM-SELF-CONSISTENCY TEST CONFIGURATION -- THE dx-SHRINKING  **
!     **  CONVERGENCE SIGNATURE OF A GENUINELY CORRECT ANALYTIC FORMULA,    **
!     **  NOT A FIT TO ONE REGIME.                                         **
!     **                                                                  **
!     **  NOT YET WIRED IN AS A DEFAULT ANYWHERE: ITS STABILITY AS A LIVE,  **
!     **  FED-BACK-INTO-PSI STEP (VS. EDFT$HAUXSTEP_ADAPTIVE'S ALREADY-     **
!     **  CHARACTERIZED LAMBDA_HAUX=0.1/0.3 BEHAVIOR, SEE THE DEV NOTES    **
!     **  .MD) IS NOT YET TESTED -- A CORRECT GRADIENT IS NOT AUTOMATICALLY**
!     **  A STABLE FIXED-LAMBDA STEP, AND THAT CHARACTERIZATION NEEDS ITS   **
!     **  OWN VALIDATION PASS BEFORE THIS REPLACES ANYTHING BY DEFAULT.     **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: LAMBDA
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      REAL(8)   ,INTENT(IN) :: CHMPOT
      REAL(8)   ,INTENT(IN) :: TEMP
      COMPLEX(8),INTENT(IN) :: RESIDUAL(NB,NB)
      REAL(8)                :: EIGVAL(NB)
      COMPLEX(8)             :: U(NB,NB)
      COMPLEX(8)             :: RESROT(NB,NB)   ! U^DAGGER*RESIDUAL*U
      COMPLEX(8)             :: G(NB,NB)        ! GRADIENT, HAUX'S OWN EIGENBASIS
      COMPLEX(8)             :: GRADLAB(NB,NB)  ! ROTATED BACK TO THE LAB FRAME
      REAL(8)                :: FOCC(NB)
      REAL(8)                :: FPRIME(NB)
      REAL(8)                :: SVAR,EX,SUMFP,AVG,DELTA
      INTEGER(4)             :: I,K,L
      REAL(8)   ,PARAMETER   :: DEGENTOL=1.D-8
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXSTEP_MVPEXACT CALLED BEFORE EDFT$INI/INIBLOCK')
        CALL ERROR$STOP('EDFT$HAUXSTEP_MVPEXACT')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXSTEP_MVPEXACT')
      END IF
      IF(TEMP.LE.0.D0) THEN
        CALL ERROR$MSG('TEMP MUST BE POSITIVE')
        CALL ERROR$R8VAL('TEMP',TEMP)
        CALL ERROR$STOP('EDFT$HAUXSTEP_MVPEXACT')
      END IF
!
!     ==========================================================================
!     ==  DIAGONALIZE THE CURRENT HAUX TO GET ITS OWN EIGENBASIS              ==
!     ==========================================================================
      CALL LIB$DIAGC8(NB,HAUX(:,:,IKPT_,ISPIN_),EIGVAL,U)
!
!     ==========================================================================
!     ==  OCCUPATIONS F_I AND THEIR LAMBDA-DERIVATIVE F'_I=-F_I*(1-F_I)/TEMP  ==
!     ==========================================================================
      DO I=1,NB
        SVAR=(EIGVAL(I)-CHMPOT)/TEMP
        IF(SVAR.GT.50.D0) THEN
          FOCC(I)=0.D0
          FPRIME(I)=0.D0
        ELSE IF(SVAR.LT.-50.D0) THEN
          FOCC(I)=1.D0
          FPRIME(I)=0.D0
        ELSE
          EX=EXP(SVAR)
          FOCC(I)=1.D0/(1.D0+EX)
          FPRIME(I)=-FOCC(I)*(1.D0-FOCC(I))/TEMP
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  ROTATE RESIDUAL (=HPHYS-HAUX) INTO HAUX'S EIGENBASIS -- SEE THE     ==
!     ==  "NO SEPARATE HPHYS ARGUMENT NEEDED" NOTE ABOVE.                     ==
!     ==========================================================================
      RESROT=MATMUL(TRANSPOSE(CONJG(U)),MATMUL(RESIDUAL,U))
!
!     ==========================================================================
!     ==  DIAGONAL: FIXED-ELECTRON-COUNT LAGRANGE CORRECTION, THEN THE        ==
!     ==  GRADIENT ITSELF. IF NO BAND IN THIS BLOCK IS FERMI-ACTIVE           ==
!     ==  (SUMFP~0, E.G. A FULLY-ISOLATED CORE-LIKE BLOCK), THE CORRECTION    ==
!     ==  IS MOOT SINCE EVERY F'_I IS ALSO ~0 AND G(I,I) VANISHES REGARDLESS. ==
!     ==========================================================================
      SUMFP=SUM(FPRIME)
      AVG=0.D0
      IF(ABS(SUMFP).GT.1.D-12) THEN
        DO I=1,NB
          AVG=AVG+FPRIME(I)*REAL(RESROT(I,I),KIND=8)
        ENDDO
        AVG=AVG/SUMFP
      END IF
      G=(0.D0,0.D0)
      DO I=1,NB
        G(I,I)=CMPLX(FPRIME(I)*(REAL(RESROT(I,I),KIND=8)-AVG),0.D0,KIND=8)
      ENDDO
!
!     ==========================================================================
!     ==  OFF-DIAGONAL: MVP'S DEGENERATE-PERTURBATION-THEORY FORMULA, WITH A  ==
!     ==  SMOOTH DEGENERACY FALLBACK (LAMBDA_K->LAMBDA_L LIMIT).              ==
!     ==========================================================================
      DO K=1,NB
        DO L=1,NB
          IF(K.EQ.L) CYCLE
          DELTA=EIGVAL(K)-EIGVAL(L)
          IF(ABS(DELTA).GT.DEGENTOL) THEN
            G(K,L)=((FOCC(K)-FOCC(L))/DELTA)*RESROT(K,L)
          ELSE
            G(K,L)=0.5D0*(FPRIME(K)+FPRIME(L))*RESROT(K,L)
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ROTATE BACK TO THE LAB FRAME (G IS HERMITIAN BY CONSTRUCTION: THE   ==
!     ==  OFF-DIAGONAL PREFACTOR IS REAL AND SYMMETRIC UNDER K<->L, AND       ==
!     ==  RESROT(L,K)=CONJG(RESROT(K,L)) SINCE RESIDUAL IS HERMITIAN) -- SO   ==
!     ==  GRADLAB IS A PROPER OPERATOR, NOT JUST A DIAGONAL-ONLY RESULT.      ==
!     ==========================================================================
      GRADLAB=MATMUL(U,MATMUL(G,TRANSPOSE(CONJG(U))))
!
!     ==========================================================================
!     ==  DESCENT STEP (HAUX -= LAMBDA*GRADLAB, SEE THE SIGN-CONVENTION NOTE  ==
!     ==  ABOVE), RE-HERMITIANIZE (SAME AS EVERY OTHER EDFT$HAUXSTEP* STEP).  ==
!     ==========================================================================
      HAUX(:,:,IKPT_,ISPIN_)=HAUX(:,:,IKPT_,ISPIN_)-CMPLX(LAMBDA,0.D0,KIND=8)*GRADLAB(:,:)
      HAUX(:,:,IKPT_,ISPIN_)=0.5D0*(HAUX(:,:,IKPT_,ISPIN_) &
     &                             +TRANSPOSE(CONJG(HAUX(:,:,IKPT_,ISPIN_))))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT_HFLATTEN(NB,HMAT,VEC)
!     ******************************************************************
!     **  FLATTENS A HERMITIAN NBxNB MATRIX INTO NB*NB INDEPENDENT      **
!     **  REAL DEGREES OF FREEDOM (NB DIAGONAL REALS, THEN REAL/IMAG    **
!     **  OF EACH UPPER-TRIANGLE OFF-DIAGONAL ELEMENT) -- NO WASTED     **
!     **  DOF, USED TO FEED EDFT$HAUXSTEP_ANDERSON'S GENERIC REAL-      **
!     **  VECTOR MIXING MACHINERY. INVERSE OF EDFT_HUNFLATTEN.          **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB
      COMPLEX(8),INTENT(IN) :: HMAT(NB,NB)
      REAL(8)   ,INTENT(OUT):: VEC(NB*NB)
      INTEGER(4)             :: I,J,K
!     ******************************************************************
      K=0
      DO I=1,NB
        K=K+1
        VEC(K)=REAL(HMAT(I,I),KIND=8)
      ENDDO
      DO I=1,NB
        DO J=I+1,NB
          K=K+1
          VEC(K)=REAL(HMAT(I,J),KIND=8)
          K=K+1
          VEC(K)=AIMAG(HMAT(I,J))
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT_HUNFLATTEN(NB,VEC,HMAT)
!     ******************************************************************
!     **  INVERSE OF EDFT_HFLATTEN: REBUILDS A HERMITIAN NBxNB MATRIX   **
!     **  FROM ITS NB*NB FLATTENED REAL DEGREES OF FREEDOM. HERMITIAN   **
!     **  BY CONSTRUCTION -- NO SEPARATE RE-HERMITIANIZATION NEEDED.    **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB
      REAL(8)   ,INTENT(IN) :: VEC(NB*NB)
      COMPLEX(8),INTENT(OUT):: HMAT(NB,NB)
      INTEGER(4)             :: I,J,K
      REAL(8)                :: RE,IM
!     ******************************************************************
      K=0
      DO I=1,NB
        K=K+1
        HMAT(I,I)=CMPLX(VEC(K),0.D0,KIND=8)
      ENDDO
      DO I=1,NB
        DO J=I+1,NB
          K=K+1
          RE=VEC(K)
          K=K+1
          IM=VEC(K)
          HMAT(I,J)=CMPLX(RE,IM,KIND=8)
          HMAT(J,I)=CMPLX(RE,-IM,KIND=8)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXSTEP_ANDERSON(IKPT_,ISPIN_,RESIDUAL)
!     ******************************************************************
!     **  ANDERSON-MIXING STEP FOR HAUX -- REPLACES EVERY EDFT$HAUXSTEP* **
!     **  ROUTINE ABOVE'S HAND-TUNED, FIXED LAMBDA_HAUX DAMPING          **
!     **  CONSTANT WITH THE SAME ADAPTIVE, RESTARTED-DIIS ANDERSON       **
!     **  ACCELERATION ALREADY USED (AND ALREADY VALIDATED) FOR THE      **
!     **  ANALOGOUS DENSITY SELF-CONSISTENCY LOOP (WAVES_MIXRHO_         **
!     **  ANDERSON, VESPAW_MIXER.F90, SELECTED VIA !PSIDYN MIXER=        **
!     **  'ANDERSON'). SAME UNDERLYING ALGORITHM (STANDARD TYPE-I        **
!     **  ANDERSON ACCELERATION, EYERT, J. COMPUT. PHYS. 124, 271        **
!     **  (1996)), JUST APPLIED TO HAUX'S OWN FIXED-POINT RESIDUAL       **
!     **  (HAUX=HPHYS[RHO(HAUX)], SEE THE FILE HEADER) INSTEAD OF THE    **
!     **  DENSITY'S.                                                    **
!     **                                                                **
!     **  NO LAMBDA ARGUMENT (UNLIKE EVERY OTHER EDFT$HAUXSTEP*          **
!     **  ROUTINE): THE DAMPING (BETA, PER (K-POINT,SPIN) BLOCK, SEE     **
!     **  EDFTMIXANDERSON_MODULE) IS ENTIRELY SELF-TUNING -- IT SHRINKS  **
!     **  WHENEVER THE RESIDUAL FAILS TO IMPROVE ON ITS BEST VALUE SINCE **
!     **  THE LAST RESTART ("RESTARTED DIIS": THE HISTORY IS ALSO       **
!     **  DISCARDED WHEN THIS HAPPENS), AND GROWS BACK TOWARD BETAMAX    **
!     **  WHILE CONVERGENCE IS STEADY -- DIRECTLY ADDRESSING THE FINDING **
!     **  (SEE THE DEV NOTES) THAT NO SINGLE LAMBDA_HAUX VALUE       **
!     **  WORKED ACROSS SYSTEMS FOR EDFT$HAUXSTEP_MVPEXACT (AL METAL     **
!     **  WANTED ~0.05-0.1, SI WANTED ~0.1 BUT WAS INSENSITIVE TO         **
!     **  OSCILLATION AT ANY TESTED VALUE -- A FIXED CONSTANT CANNOT     **
!     **  SERVE BOTH CASES).                                             **
!     **                                                                **
!     **  ONE INDEPENDENT ANDERSON-MIXING HISTORY PER (K-POINT,SPIN)     **
!     **  BLOCK (UNLIKE DENSITY MIXING'S SINGLE GLOBAL HISTORY): HAUX'S  **
!     **  FIXED POINT IS LOCAL TO EACH BLOCK, SO ITS RESIDUAL HISTORY    **
!     **  SHOULD BE TOO -- POOLING HISTORY ACROSS UNRELATED K-POINTS     **
!     **  WOULD MAKE THE LEAST-SQUARES FIT MEANINGLESS.                  **
!     **                                                                **
!     **  RESIDUAL (=HPHYS-HAUX, FROM EDFT$HAUXFORCE) IS EXACTLY THE     **
!     **  "OUTPUT MINUS INPUT" FIXED-POINT RESIDUAL THE ANDERSON RECIPE  **
!     **  EXPECTS -- NO INTERNAL RECOMPUTATION NEEDED (UNLIKE            **
!     **  WAVES_MIXRHO_ANDERSON, WHICH COMPUTES ITS OWN RES=RHO-XOLD     **
!     **  SINCE ITS CALLER DOESN'T ALREADY HAVE ONE ON HAND).            **
!     **                                                                **
!     **  HAUX IS HERMITIAN (NB COMPLEX NUMBERS' WORTH OF REDUNDANCY),   **
!     **  SO IT AND ITS RESIDUAL ARE FLATTENED TO NB*NB INDEPENDENT      **
!     **  REAL DEGREES OF FREEDOM (EDFT_HFLATTEN/EDFT_HUNFLATTEN,        **
!     **  ABOVE) -- NO WASTED DOF, NO NEED TO MIX A REDUNDANT FULL       **
!     **  NBxNB REAL/IMAGINARY REPRESENTATION. PLAIN FORTRAN ARRAY       **
!     **  INTRINSICS (DOT_PRODUCT), NOT THE BLAS-BACKED LIB$             **
!     **  SCALARPRODUCTR8 THE DENSITY MIXERS USE, SINCE NB*NB IS TINY    **
!     **  (A HANDFUL OF BANDS) COMPARED TO A REAL-SPACE GRID -- NO       **
!     **  PERFORMANCE REASON TO BOTHER WITH THE BLAS WRAPPER HERE.       **
!     ******************************************************************
      USE EDFT_MODULE
      USE EDFTMIXANDERSON_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      COMPLEX(8),INTENT(IN) :: RESIDUAL(NB,NB)
      INTEGER(4)             :: NTOT
      REAL(8)                :: RES(NB*NB)
      REAL(8)                :: XPROP(NB*NB)
      REAL(8)                :: DF(NB*NB,NHISTORY)
      REAL(8)                :: A(NHISTORY,NHISTORY)
      REAL(8)                :: AINV(NHISTORY,NHISTORY)
      REAL(8)                :: B(NHISTORY)
      REAL(8)                :: GAMMA(NHISTORY)
      REAL(8)                :: RESVAR,SVAR,STEPVAR
      INTEGER(4)             :: IH1,IH2,NBACK
      LOGICAL(4)             :: TRESTART
      LOGICAL(4)             :: TDEBUG
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXSTEP_ANDERSON CALLED BEFORE EDFT$INI/INIBLOCK')
        CALL ERROR$STOP('EDFT$HAUXSTEP_ANDERSON')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXSTEP_ANDERSON')
      END IF
      NTOT=NB*NB
!
!     ==========================================================================
!     ==  LAZY ALLOCATION OF THE PER-BLOCK ANDERSON-MIXING HISTORY, FIRST     ==
!     ==  CALL EVER (ANY BLOCK). SIZED BY NB/NKPT/NSPIN, ALREADY FIXED BY     ==
!     ==  EDFT$INI/EDFT$INIBLOCK BY THE TIME THIS IS FIRST CALLED.           ==
!     ==========================================================================
      IF(.NOT.TALLOCATED) THEN
        ALLOCATE(TFIRST(NKPT,NSPIN))
        ALLOCATE(BETA(NKPT,NSPIN))
        ALLOCATE(RESBEST(NKPT,NSPIN))
        ALLOCATE(NVALID(NKPT,NSPIN))
        ALLOCATE(XOLD(NTOT,NHISTORY+1,NKPT,NSPIN))
        ALLOCATE(FOLD(NTOT,NHISTORY+1,NKPT,NSPIN))
        TFIRST=.TRUE.
        BETA=BETAMAX
        RESBEST=HUGE(1.D0)
        NVALID=0
        XOLD=0.D0
        FOLD=0.D0
        TALLOCATED=.TRUE.
      END IF
!
!     ==========================================================================
!     ==  FIRST-EVER CALL FOR *THIS* BLOCK: JUST SEED THE HISTORY WITH THE    ==
!     ==  CURRENT HAUX (NO RESIDUAL/INPUT PAIR TO MIX FROM YET) -- SAME       ==
!     ==  TFIRST-RETURN-EARLY PATTERN AS WAVES_MIXRHO_ANDERSON.               ==
!     ==========================================================================
      IF(TFIRST(IKPT_,ISPIN_)) THEN
        CALL EDFT_HFLATTEN(NB,HAUX(:,:,IKPT_,ISPIN_),XOLD(:,1,IKPT_,ISPIN_))
        TFIRST(IKPT_,ISPIN_)=.FALSE.
        RETURN
      END IF
!
      CALL VERBOSITY$GETL4('DEBUG',TDEBUG)
      CALL EDFT_HFLATTEN(NB,RESIDUAL,RES)
      RESVAR=DOT_PRODUCT(RES,RES)
      IF(TDEBUG) PRINT*,'EDFT ANDERSON-MIXER IKPT=',IKPT_,' ISPIN=',ISPIN_ &
     &                 ,' |RESIDUAL|**2=',RESVAR
!
!     ==========================================================================
!     ==  RESTARTED DIIS: SAME RATIONALE AS WAVES_MIXRHO_ANDERSON             ==
!     ==========================================================================
      TRESTART=(RESVAR.GT.RESBEST(IKPT_,ISPIN_))
      IF(TRESTART) THEN
        BETA(IKPT_,ISPIN_)=MAX(BETA(IKPT_,ISPIN_)*0.5D0,BETAMIN)
        NVALID(IKPT_,ISPIN_)=0
        RESBEST(IKPT_,ISPIN_)=RESVAR
        IF(TDEBUG) PRINT*,'EDFT ANDERSON-MIXER: RESTARTING HISTORY, NEW BETA=' &
     &                    ,BETA(IKPT_,ISPIN_)
      ELSE
        RESBEST(IKPT_,ISPIN_)=RESVAR
        BETA(IKPT_,ISPIN_)=MIN(BETA(IKPT_,ISPIN_)*1.1D0,BETAMAX)
      END IF
!
      FOLD(:,1,IKPT_,ISPIN_)=RES(:)
      NBACK=MIN(NHISTORY,NVALID(IKPT_,ISPIN_))
!
!     ==========================================================================
!     ==  LEAST-SQUARES FIT: MINIMIZE ||RES-SUM_I GAMMA_I*DF_I||**2,          ==
!     ==  DF_I=RES-FOLD(:,I+1) -- SAME RECIPE AS WAVES_MIXRHO_ANDERSON        ==
!     ==========================================================================
      IF(NBACK.EQ.0) THEN
        XPROP=XOLD(:,1,IKPT_,ISPIN_)+BETA(IKPT_,ISPIN_)*RES(:)
      ELSE
        DO IH1=1,NBACK
          DF(:,IH1)=FOLD(:,1,IKPT_,ISPIN_)-FOLD(:,IH1+1,IKPT_,ISPIN_)
        ENDDO
        DO IH1=1,NBACK
          DO IH2=1,NBACK
            A(IH1,IH2)=DOT_PRODUCT(DF(:,IH1),DF(:,IH2))
          ENDDO
        ENDDO
        DO IH1=1,NBACK
          B(IH1)=DOT_PRODUCT(DF(:,IH1),RES(:))
        ENDDO
!       == TIKHONOV REGULARIZATION, SAME RATIONALE AS THE DENSITY MIXERS ==
        SVAR=0.D0
        DO IH1=1,NBACK
          SVAR=SVAR+A(IH1,IH1)
        ENDDO
        SVAR=MAX(SVAR/REAL(NBACK,KIND=8)*1.D-6,1.D-12)
        DO IH1=1,NBACK
          A(IH1,IH1)=A(IH1,IH1)+SVAR
        ENDDO
        CALL LIB$INVERTR8(NBACK,A(1:NBACK,1:NBACK),AINV(1:NBACK,1:NBACK))
        GAMMA=0.D0
        DO IH1=1,NBACK
          DO IH2=1,NBACK
            GAMMA(IH1)=GAMMA(IH1)+AINV(IH1,IH2)*B(IH2)
          ENDDO
        ENDDO
        XPROP=XOLD(:,1,IKPT_,ISPIN_)+BETA(IKPT_,ISPIN_)*RES(:)
        DO IH1=1,NBACK
          XPROP=XPROP-GAMMA(IH1)*((XOLD(:,1,IKPT_,ISPIN_)-XOLD(:,IH1+1,IKPT_,ISPIN_)) &
     &                            +BETA(IKPT_,ISPIN_)*DF(:,IH1))
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  OVERSHOOT SAFEGUARD, SAME IDEA AS THE DENSITY MIXERS                ==
!     ==========================================================================
      STEPVAR=DOT_PRODUCT(XPROP-XOLD(:,1,IKPT_,ISPIN_),XPROP-XOLD(:,1,IKPT_,ISPIN_))
      IF(STEPVAR.GT.100.D0*MAX(RESVAR,1.D-12)) THEN
        XPROP=XOLD(:,1,IKPT_,ISPIN_)+BETA(IKPT_,ISPIN_)*RES(:)
      END IF
!
      CALL EDFT_HUNFLATTEN(NB,XPROP,HAUX(:,:,IKPT_,ISPIN_))
!
!     ==========================================================================
!     ==  PROPAGATE HISTORY (THIS BLOCK ONLY)                                 ==
!     ==========================================================================
      DO IH1=NHISTORY+1,2,-1
        XOLD(:,IH1,IKPT_,ISPIN_)=XOLD(:,IH1-1,IKPT_,ISPIN_)
        FOLD(:,IH1,IKPT_,ISPIN_)=FOLD(:,IH1-1,IKPT_,ISPIN_)
      ENDDO
      XOLD(:,1,IKPT_,ISPIN_)=XPROP
      NVALID(IKPT_,ISPIN_)=NVALID(IKPT_,ISPIN_)+1
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE EDFT$HAUXSTEP_RESGATE(LAMBDA,IKPT_,ISPIN_,CHMPOT,TEMP,RESIDUAL)
!     ******************************************************************
!     **  SAME AS EDFT$HAUXSTEP_ADAPTIVE, EXCEPT THE GATE THAT DECIDES     **
!     **  WHICH BANDS GET THE FULL LAMBDA-STEP IS THE *OR* OF TWO          **
!     **  INDEPENDENT SIGNALS INSTEAD OF JUST ONE. FOUND 2026-09-06, BY     **
!     **  DIRECT INSTRUMENTATION OF HAUX'S OWN EIGENVALUES ITERATION BY     **
!     **  ITERATION ON THE SI TEST CASE (SEE THE DEV NOTES): A SINGLE   **
!     **  BAND'S GATE SAT AT ~6.D-5 (ESSENTIALLY SHUT) FOR OVER 150         **
!     **  ITERATIONS WHILE ITS RESIDUAL STAYED PINNED AT ~0.17 HARTREE --   **
!     **  NEITHER MOVING, A GENUINE SELF-REINFORCING LOCK, NOT SLOW         **
!     **  CONVERGENCE. NUMERICALLY CONFIRMED: LAMBDA*GATE*RESIDUAL          **
!     **  REPRODUCED THE OBSERVED (TINY) PER-ITERATION EIGENVALUE CHANGE    **
!     **  EXACTLY.                                                         **
!     **                                                                   **
!     **  ROOT CAUSE: EDFT$HAUXSTEP_ADAPTIVE'S GATE_I=MIN(1,W_I/WREF) IS A  **
!     **  FUNCTION OF THE BAND'S *CURRENT* EIGENVALUE'S DISTANCE FROM THE   **
!     **  CHEMICAL POTENTIAL ONLY -- IT HAS NO WAY TO KNOW WHETHER THE      **
!     **  RESIDUAL ITSELF IS ACTUALLY LARGE. WHEN A NATURAL ORBITAL'S       **
!     **  IDENTITY HAS DRIFTED FAR FROM WHERE HPHYS SAYS IT SHOULD BE, BUT   **
!     **  ITS *CURRENT* (WRONG) EIGENVALUE STILL SITS DEEP IN THE OCCUPIED   **
!     **  OR EMPTY MANIFOLD, THE GATE READS "SETTLED, LEAVE ALONE" AND       **
!     **  SUPPRESSES EXACTLY THE CORRECTION THAT WOULD FIX IT -- A LOCK      **
!     **  THAT CAN ONLY EVER OPEN AS FAST AS THE (ALREADY-SUPPRESSED)        **
!     **  EIGENVALUE ITSELF CREEPS TOWARD THE FERMI WINDOW.                  **
!     **                                                                   **
!     **  FIX: GATE_I=MAX(GATE_FERMI_I,GATE_RESIDUAL_I), WHERE               **
!     **  GATE_RESIDUAL_I=MIN(1,|RESEIG_II|/RESREF) OPENS DIRECTLY ON THE    **
!     **  SIZE OF THE (ALREADY-COMPUTED) DIAGONAL RESIDUAL ITSELF, NO NEW    **
!     **  QUANTITY NEEDED. EITHER SIGNAL CAN OPEN THE GATE, SO A BAND THAT   **
!     **  IS BOTH SETTLED (FERMI-WISE) *AND* GENUINELY SMALL-RESIDUAL --     **
!     **  THE ORIGINAL "LEAVE IT ALONE" CASE THE GATE WAS DESIGNED FOR --    **
!     **  IS UNAFFECTED; ONLY THE NEWLY-FOUND LOCK CASE (FERMI-GATE SHUT,    **
!     **  RESIDUAL GENUINELY LARGE) NOW ALSO OPENS.                          **
!     **                                                                   **
!     **  RESREF IS RUNTIME-CONFIGURABLE VIA VESPAW_EDFT_RESGATE_RESREF     **
!     **  (DEFAULT 1.D-3 HARTREE -- COMFORTABLY ABOVE A CONVERGED           **
!     **  RESIDUAL'S SCALE (~1.D-6-1.D-4, SEE THE EDIFF/MIXTOL SCALES USED  **
!     **  THROUGHOUT THIS PROJECT) AND COMFORTABLY BELOW THE ~0.17 HARTREE  **
!     **  LOCK FOUND ON THE SI TEST -- NOT YET TUNED AGAINST A REAL RUN THE  **
!     **  WAY WREF WAS FOR EDFT$HAUXSTEP_ADAPTIVE ITSELF).                  **
!     **                                                                   **
!     **  NOT YET THE DEFAULT ANYWHERE -- OPT-IN VIA VESPAW_EDFT_RESGATE,   **
!     **  SAME PATTERN AS EVERY OTHER EXPERIMENTAL STEP IN THIS FILE.       **
!     ******************************************************************
      USE EDFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: LAMBDA
      INTEGER(4),INTENT(IN) :: IKPT_
      INTEGER(4),INTENT(IN) :: ISPIN_
      REAL(8)   ,INTENT(IN) :: CHMPOT
      REAL(8)   ,INTENT(IN) :: TEMP
      COMPLEX(8),INTENT(IN) :: RESIDUAL(NB,NB)
      REAL(8)                :: EIGVAL(NB)
      COMPLEX(8)             :: U(NB,NB)
      COMPLEX(8)             :: RESEIG(NB,NB)
      REAL(8)                :: WEIGHT(NB)
      REAL(8)                :: GATEFERMI(NB),GATERESIDUAL(NB),GATE(NB)
      REAL(8)                :: SVAR,EX
      INTEGER(4)             :: I,J
      LOGICAL(4),SAVE        :: TWREFCHECKED=.FALSE.
      REAL(8)   ,SAVE        :: WREF=1.D-6
      CHARACTER(32)          :: CENVVALWREF
      LOGICAL(4),SAVE        :: TRESREFCHECKED=.FALSE.
      REAL(8)   ,SAVE        :: RESREF=1.D-3
      CHARACTER(32)          :: CENVVALRESREF
!     ******************************************************************
      IF(.NOT.TWREFCHECKED) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_EDFT_ADAPTIVE_WREF',CENVVALWREF)
        IF(LEN_TRIM(CENVVALWREF).GT.0) THEN
          READ(CENVVALWREF,*) WREF
          PRINT*,'VESPAW_EDFT_ADAPTIVE_WREF SET -- USING WREF=',WREF
        END IF
        TWREFCHECKED=.TRUE.
      END IF
      IF(.NOT.TRESREFCHECKED) THEN
        CALL GET_ENVIRONMENT_VARIABLE('VESPAW_EDFT_RESGATE_RESREF',CENVVALRESREF)
        IF(LEN_TRIM(CENVVALRESREF).GT.0) THEN
          READ(CENVVALRESREF,*) RESREF
          PRINT*,'VESPAW_EDFT_RESGATE_RESREF SET -- USING RESREF=',RESREF
        END IF
        TRESREFCHECKED=.TRUE.
      END IF
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('EDFT$HAUXSTEP_RESGATE CALLED BEFORE EDFT$INI/INIBLOCK')
        CALL ERROR$STOP('EDFT$HAUXSTEP_RESGATE')
      END IF
      IF(IKPT_.LT.1.OR.IKPT_.GT.NKPT.OR.ISPIN_.LT.1.OR.ISPIN_.GT.NSPIN) THEN
        CALL ERROR$MSG('IKPT_ OR ISPIN_ OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT_',IKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$I4VAL('ISPIN_',ISPIN_)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('EDFT$HAUXSTEP_RESGATE')
      END IF
      IF(TEMP.LE.0.D0) THEN
        CALL ERROR$MSG('TEMP MUST BE POSITIVE')
        CALL ERROR$R8VAL('TEMP',TEMP)
        CALL ERROR$STOP('EDFT$HAUXSTEP_RESGATE')
      END IF
!
!     ==========================================================================
!     ==  DIAGONALIZE THE CURRENT HAUX TO GET ITS OWN EIGENBASIS             ==
!     ==========================================================================
      CALL LIB$DIAGC8(NB,HAUX(:,:,IKPT_,ISPIN_),EIGVAL,U)
!
!     ==========================================================================
!     ==  FERMI-WINDOW GATE, IDENTICAL TO EDFT$HAUXSTEP_ADAPTIVE'S OWN.       ==
!     ==========================================================================
      DO I=1,NB
        SVAR=(EIGVAL(I)-CHMPOT)/TEMP
        IF(ABS(SVAR).GT.50.D0) THEN
          WEIGHT(I)=0.D0
        ELSE
          EX=EXP(SVAR)
          WEIGHT(I)=4.D0*EX/((1.D0+EX)**2)
        END IF
        GATEFERMI(I)=MIN(1.D0,WEIGHT(I)/WREF)
      ENDDO
!
!     ==========================================================================
!     ==  ROTATE THE RESIDUAL INTO HAUX'S EIGENBASIS FIRST (NEEDED FOR THE    ==
!     ==  NEW RESIDUAL-SIZE GATE TOO, SO COMPUTED ONCE, BEFORE GATING).       ==
!     ==========================================================================
      RESEIG=MATMUL(TRANSPOSE(CONJG(U)),MATMUL(RESIDUAL,U))
!
!     ==========================================================================
!     ==  NEW: RESIDUAL-SIZE GATE, OPENS DIRECTLY ON |RESEIG_II REGARDLESS   ==
!     ==  OF WHERE THIS BAND'S EIGENVALUE CURRENTLY SITS. COMBINED WITH THE  ==
!     ==  FERMI GATE VIA MAX (EITHER SIGNAL CAN OPEN THE GATE) -- SEE THE    ==
!     ==  HEADER COMMENT FOR WHY THIS DOESN'T REOPEN THE OLD MVP REGRESSION. ==
!     ==========================================================================
      DO I=1,NB
        GATERESIDUAL(I)=MIN(1.D0,ABS(RESEIG(I,I))/RESREF)
        GATE(I)=MAX(GATEFERMI(I),GATERESIDUAL(I))
      ENDDO
!
!     ==========================================================================
!     ==  GATE (DIAGONAL DIRECTLY, OFF-DIAGONAL VIA GEOMETRIC MEAN -- SAME    ==
!     ==  CONVENTION AS EDFT$HAUXSTEP_ADAPTIVE), TRANSFORM BACK.              ==
!     ==========================================================================
      DO I=1,NB
        DO J=1,NB
          RESEIG(I,J)=RESEIG(I,J)*SQRT(MAX(GATE(I)*GATE(J),0.D0))
        ENDDO
      ENDDO
      RESEIG=MATMUL(U,MATMUL(RESEIG,TRANSPOSE(CONJG(U))))
!
!     ==========================================================================
!     ==  DAMPED STEP, RE-HERMITIANIZE (SAME AS EDFT$HAUXSTEP_ADAPTIVE)       ==
!     ==========================================================================
      HAUX(:,:,IKPT_,ISPIN_)=HAUX(:,:,IKPT_,ISPIN_)+CMPLX(LAMBDA,0.D0,KIND=8)*RESEIG(:,:)
      HAUX(:,:,IKPT_,ISPIN_)=0.5D0*(HAUX(:,:,IKPT_,ISPIN_) &
     &                             +TRANSPOSE(CONJG(HAUX(:,:,IKPT_,ISPIN_))))
      RETURN
      END
