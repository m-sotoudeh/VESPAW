!      MODULE VERSION_MODULE
!        CHARACTER(256):: VERTYP=_VERTYP
!        CHARACTER(256):: VERINF=_VERINF
!        CHARACTER(256):: VERREV=_VERREV
!        CHARACTER(256):: VERAUT=_VERAUT
!        CHARACTER(256):: VERDAT=_VERDAT
!      END MODULE VERSION_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAW_VERSION()
!     **************************************************************************
!     ** CHECKS CALLING SEQUENCE AND RESOLVES VERSION INFORMATION
!     **************************************************************************
      USE VERSION_MODULE
      USE STRINGS_MODULE
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT
      IMPLICIT NONE
      INTEGER(4)      :: ISVAR
      CHARACTER(512)  :: CNTLNAME
      CHARACTER(32)   :: CH32SVAR
      CHARACTER(256)  :: VERSIONTEXT
      CHARACTER(256)  :: VERSIONINFO
      INTEGER(4)      :: NTASKS,THISTASK
      COMMON/VERSION/VERSIONTEXT
!     **************************************************************************
!     ==========================================================================
!     == THESE LINES CONTAIN INFORMATION THAT CAN BE GREPPED OUT OF THE       ==
!     == EXECUTABLE. THEY HAVE NO FUNCTION IN THE EXECUTION OF THE CODE.      ==
!     ==========================================================================
      VERSIONINFO = '@(#) PAW-VERSION %R% CREATED %U% %E%'
      VERSIONTEXT = 'PROGRAM VERSION %R% CREATED %U% %E%'
!
!     ==========================================================================
!     == STARTUP BANNER (PRINTED ON EVERY INVOCATION, INCLUDING --HELP)       ==
!     ==========================================================================
!     == MPE$INIT RUNS BEFORE THIS ROUTINE: PRINT THE BANNER ON THE FIRST  ==
!     == TASK ONLY, NOT ONCE PER MPI TASK                                   ==
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) CALL VESPAW_PRINTBANNER(OUTPUT_UNIT)
!
!     ==========================================================================
!     == CHECK THAT AT LEAST ONE ARGUMENT IS GIVEN                            ==
!     ==========================================================================
      CALL LIB$NARGS(ISVAR)
      IF(ISVAR.LT.1) THEN
        CALL ERROR$MSG('THE NAME OF THE CONTROLFILE')
        CALL ERROR$MSG('MUST BE GIVEN AS ARGUMENT')
        CALL ERROR$STOP('READIN')
      END IF
!
!     ==========================================================================
!     == CHECK THAT AT LEAST ONE ARGUMENT IS GIVEN                            ==
!     ==========================================================================
      CALL LIB$GETARG(1,CNTLNAME)
      IF(+CNTLNAME.EQ.'--HELP'.OR.+CNTLNAME.EQ.'-H'.OR.+CNTLNAME.EQ.'?') THEN
        WRITE(*,FMT='("USAGE:"/"=====")')
        WRITE(*,FMT='(T5,"1ST ARGUMENT",T25,"ACTION"/T5,32("-"))')
        WRITE(*,FMT='(T5,"CONTROL FILENAME",T25,"EXECUTE ")')
        CH32SVAR=-'--HELP'
        WRITE(*,FMT='(T5,A6,T25,"PRINT HELP INFORMATION")')TRIM(CH32SVAR)
        CH32SVAR=-'-H'
        WRITE(*,FMT='(T5,A2,T25,"PRINT HELP INFORMATION")')TRIM(CH32SVAR)
        WRITE(*,FMT='(T5,"?",T25,"PRINT HELP INFORMATION")')
        CH32SVAR=-'--VERSION'
        WRITE(*,FMT='(T5,A9,T25,"PRINT VERSION INFORMATION")')TRIM(CH32SVAR)
        CH32SVAR=-'--PARMFILE'
        WRITE(*,FMT='(T5,A10,T25,"PRINT PARMFILE")')TRIM(CH32SVAR)
        WRITE(*,FMT='(T5,A10,T25,"ERRORS AFTER THIS LINE ARE irrelevant")')
        CALL ERROR$NORMALSTOP
!
      ELSE IF(+CNTLNAME.EQ.'--VERSION') THEN
        CALL VESPAW_PRINTVERSIONBLOCK(OUTPUT_UNIT)
        WRITE(*,FMT='(T5,A10,T25,"ERRORS AFTER THIS LINE ARE irRELEVANT")')
        CALL ERROR$NORMALSTOP
      ELSE IF(+CNTLNAME.EQ.'--PARMFILE') THEN
        CALL VERSION$WRITEPARMFILE()
        WRITE(*,FMT='(T5,A10,T25,"ERRORS AFTER THIS LINE ARE irRELEVANT")')
        CALL ERROR$NORMALSTOP
      END IF
      return
      end
!
!     ..................................................................
      MODULE PROVENANCE_MODULE
!     **************************************************************************
!     ** STASHES THE SAME VERSION/BACKEND STRINGS VERSION_PRINTDEV PRINTS IN   **
!     ** THE STARTUP BANNER/.PROT HEADER, SO OTHER CODE (E.G. THE HDF5        **
!     ** RESULTS WRITER, src/vespaw_ioroutines.f90) CAN RECORD THE SAME        **
!     ** PROVENANCE FACTS WITHOUT RE-PARSING getversion.sh'S RAW SENTENCES.   **
!     **************************************************************************
      CHARACTER(32) :: PVERSION='UNKNOWN'
      CHARACTER(16) :: PREVISION='UNKNOWN'
      CHARACTER(32) :: PGITSTATE='UNKNOWN'
      CHARACTER(64) :: PBRANCH='UNKNOWN'
      CHARACTER(16) :: PCPPAWBASE='UNKNOWN'
      CHARACTER(64) :: PCOMPILER='UNKNOWN'
      CHARACTER(32) :: PMPI='UNKNOWN'
      CHARACTER(32) :: PBLAS='UNKNOWN'
      CHARACTER(16) :: PFFT='UNKNOWN'
      CHARACTER(16) :: PGPU='UNKNOWN'
      CHARACTER(64) :: PBUILT='UNKNOWN'
      LOGICAL(4)    :: PSET=.FALSE.
      END MODULE PROVENANCE_MODULE
!
!     ..................................................................
      SUBROUTINE PROVENANCE$GET(ID_,VAL_)
!     **************************************************************************
!     ** RETURNS ONE PROVENANCE STRING BY NAME. VALID BEFORE THE FIRST CALL   **
!     ** TO VERSION_PRINTDEV (I.E. BEFORE PUTHEADER/--VERSION HAS RUN ONCE)   **
!     ** ONLY THE DEFAULT 'UNKNOWN' PLACEHOLDERS ARE RETURNED.                **
!     **************************************************************************
      USE PROVENANCE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID_
      CHARACTER(*),INTENT(OUT) :: VAL_
!     **************************************************************************
      IF(ID_.EQ.'VERSION') THEN
        VAL_=PVERSION
      ELSE IF(ID_.EQ.'REVISION') THEN
        VAL_=PREVISION
      ELSE IF(ID_.EQ.'GITSTATE') THEN
        VAL_=PGITSTATE
      ELSE IF(ID_.EQ.'BRANCH') THEN
        VAL_=PBRANCH
      ELSE IF(ID_.EQ.'CPPAWBASE') THEN
        VAL_=PCPPAWBASE
      ELSE IF(ID_.EQ.'COMPILER') THEN
        VAL_=PCOMPILER
      ELSE IF(ID_.EQ.'MPI') THEN
        VAL_=PMPI
      ELSE IF(ID_.EQ.'BLAS') THEN
        VAL_=PBLAS
      ELSE IF(ID_.EQ.'FFT') THEN
        VAL_=PFFT
      ELSE IF(ID_.EQ.'GPU') THEN
        VAL_=PGPU
      ELSE IF(ID_.EQ.'BUILT') THEN
        VAL_=PBUILT
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('PROVENANCE$GET')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE VERSION_PRINTDEV(NFIL,VERINF,VERREV,VERAUT,VERDAT)
!     **************************************************************************
!     ** PARSES THE RAW GIT-DERIVED SENTENCES BUILT BY                        **
!     ** src/buildtools/Version/getversion.sh, ADDS COMPILE-TIME BACKEND      **
!     ** INFO (COMPILER/MPI/BLAS/FFT/GPU), AND PRINTS THE                    **
!     ** "VESPAW VERSION INFORMATION" BLOCK. ONLY VALID WHEN                  **
!     ** VERTYP='DEVELOPMENT VERSION' -- THE 'RELEASE VERSION'/               **
!     ** 'UNKNOWN VERSION' SENTENCES HAVE A DIFFERENT SHAPE AND ARE PRINTED   **
!     ** RAW BY THE CALLER (VESPAW_PRINTVERSIONBLOCK) INSTEAD.                **
!     **                                                                      **
!     ** EXPECTED INPUT SHAPES (SEE getversion.sh):                          **
!     **   VERINF='branch: <BRANCH>; remotes: <REMOTE>;'                     **
!     **   VERREV='revision <N>; <HASH>; <M> changes since last commit'      **
!     **   VERAUT='last commit by <AUTHOR> at <DATE>'                        **
!     **   VERDAT='compiled by <USER> on <HOST> at <DATE>'                   **
!     **                                                                      **
!     ** NFIL IS THE OUTPUT UNIT -- OUTPUT_UNIT FOR THE SCREEN (--VERSION),   **
!     ** OR THE 'PROT' FILE UNIT WHEN CALLED FROM PUTHEADER (paw_ioroutines). **
!     **                                                                      **
!     ** THE 'CP-PAW BASE' COMMIT AND THE 'VERSION' STRING ARE HAND-         **
!     ** MAINTAINED CONSTANTS, NOT DERIVED FROM GIT: 'ABFB75E' IS THE         **
!     ** COMMIT THIS FORK RE-BRANDED FROM CP-PAW AT (SEE THE DEV NOTES),  **
!     ** A FIXED HISTORICAL FACT THAT DOES NOT CHANGE WITH FUTURE COMMITS.    **
!     **************************************************************************
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : COMPILER_VERSION
      USE PROVENANCE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: VERINF,VERREV,VERAUT,VERDAT
      CHARACTER(*),PARAMETER  :: VESPAWVERSION='1.0.0'
      CHARACTER(*),PARAMETER  :: CPPAWBASE='abfb75e'
      CHARACTER(256)          :: SVAR,REST
      CHARACTER(256)          :: BRANCH
      CHARACTER(256)          :: REVHASH,REVCHANGES
      CHARACTER(256)          :: BDATE
      CHARACTER(64)           :: COMPILERSTR
      CHARACTER(32)           :: MPISTR
      CHARACTER(32)           :: BLASSTR
      CHARACTER(16)           :: FFTSTR
      CHARACTER(16)           :: GPUSTR
      CHARACTER(12)           :: LBL
      INTEGER(4)              :: I1,I2
!     **************************************************************************
!     == PARSE VERINF: 'branch: X; remotes: Y;' -- ONLY BRANCH IS KEPT ========
      SVAR=TRIM(VERINF)
      I1=INDEX(SVAR,'branch: ')
      I2=INDEX(SVAR,'; remotes:')
      IF(I1.GT.0.AND.I2.GT.I1) THEN
        BRANCH=SVAR(I1+8:I2-1)
      ELSE
        BRANCH='?'
      END IF
!
!     == PARSE VERREV: 'revision N; HASH; M changes since last commit' ========
      SVAR=TRIM(VERREV)
      I2=INDEX(SVAR,';')
      REVHASH='unknown '
      REVCHANGES='0'
      IF(I2.GT.0) THEN
        REST=ADJUSTL(SVAR(I2+1:))
        I1=INDEX(REST,';')
        IF(I1.GT.0) THEN
          REVHASH=ADJUSTL(REST(1:I1-1))
          REST=ADJUSTL(REST(I1+1:))
          I2=INDEX(REST,' changes since last commit')
          IF(I2.GT.0) REVCHANGES=ADJUSTL(REST(1:I2-1))
        END IF
      END IF
!
!     == PARSE VERDAT: 'compiled by X on Y at Z' -- ONLY BUILD DATE KEPT ======
      SVAR=TRIM(VERDAT)
      I1=INDEX(SVAR,' at ',BACK=.TRUE.)
      IF(I1.GT.0) THEN
        BDATE=ADJUSTL(SVAR(I1+4:))
      ELSE
        BDATE='?'
      END IF
!
!     == COMPILE-TIME BACKEND INFO =============================================
      COMPILERSTR=COMPILER_VERSION()
#IFDEF CPPVARIABLE_PARALLEL
      MPISTR='enabled (parallel build)'
#ELSE
      MPISTR='none (serial build)'
#ENDIF
#IF DEFINED(CPPVAR_BLAS_ESSL)
      BLASSTR='IBM ESSL'
#ELIF DEFINED(CPPVAR_BLAS_MKL)
      BLASSTR='Intel MKL'
#ELIF DEFINED(CPPVAR_BLAS_OPENBLAS)
      BLASSTR='OpenBLAS'
#ELSE
      BLASSTR='generic BLAS/LAPACK'
#ENDIF
#IF DEFINED(CPPVAR_FFT_FFTW3)
      FFTSTR='FFTW3'
#ELIF DEFINED(CPPVAR_FFT_FFTW)
      FFTSTR='FFTW2'
#ELIF DEFINED(CPPVAR_FFT_ESSL)
      FFTSTR='IBM ESSL'
#ELIF DEFINED(CPPVAR_FFT_ACML)
      FFTSTR='AMD ACML'
#ELIF DEFINED(CPPVAR_FFT_CUFFT)
      FFTSTR='cuFFT'
#ELIF DEFINED(CPPVAR_FFT_PACK)
      FFTSTR='FFTPACK'
#ELSE
      FFTSTR='unknown'
#ENDIF
#IF DEFINED(CPPVAR_FFT_CUFFT)
      GPUSTR='CUDA (cuFFT)'
#ELSE
      GPUSTR='none'
#ENDIF
!
!     == PRINT NICELY ==========================================================
      LBL='Version'      ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,VESPAWVERSION
      LBL='Revision'     ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,REVHASH(1:8)
      LBL='Git state'
      IF(TRIM(REVCHANGES).EQ.'0') THEN
        WRITE(NFIL,FMT='(A12,": ","clean")')LBL
      ELSE
        WRITE(NFIL,FMT='(A12,": dirty (",A," changes)")')LBL,TRIM(REVCHANGES)
      END IF
      LBL='Branch'       ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,TRIM(BRANCH)
      LBL='CP-PAW base'  ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,CPPAWBASE
      WRITE(NFIL,FMT='()')
      LBL='Compiler'     ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,TRIM(COMPILERSTR)
      LBL='MPI'          ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,TRIM(MPISTR)
      LBL='BLAS'         ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,TRIM(BLASSTR)
      LBL='FFT'          ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,TRIM(FFTSTR)
      LBL='GPU'          ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,TRIM(GPUSTR)
      WRITE(NFIL,FMT='()')
      LBL='Built'        ;WRITE(NFIL,FMT='(A12,": ",A)')LBL,TRIM(BDATE)
!
!     == STASH FOR PROVENANCE$GET (E.G. THE HDF5 RESULTS WRITER) ===============
      PVERSION=VESPAWVERSION
      PREVISION=REVHASH(1:8)
      IF(TRIM(REVCHANGES).EQ.'0') THEN
        PGITSTATE='clean'
      ELSE
        PGITSTATE='dirty ('//TRIM(REVCHANGES)//' changes)'
      END IF
      PBRANCH=TRIM(BRANCH)
      PCPPAWBASE=CPPAWBASE
      PCOMPILER=TRIM(COMPILERSTR)
      PMPI=TRIM(MPISTR)
      PBLAS=TRIM(BLASSTR)
      PFFT=TRIM(FFTSTR)
      PGPU=TRIM(GPUSTR)
      PBUILT=TRIM(BDATE)
      PSET=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE VESPAW_PRINTBANNER(NFIL)
!     **************************************************************************
!     ** PRINTS THE VESPAW BEE-LOGO + WORDMARK BANNER TO UNIT NFIL.           **
!     ** SHARED BY THE STARTUP SCREEN BANNER (PAW_VERSION, EVERY INVOCATION)  **
!     ** AND THE .PROT FILE HEADER (PUTHEADER, src/paw_ioroutines.f90).       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
!     **************************************************************************
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='(T39,". .")')
      WRITE(NFIL,FMT='(T39,"( )")')
      WRITE(NFIL,FMT='(T37,"<< = >>")')
      WRITE(NFIL,FMT='(T40,"V")')
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='(T20,"__     _______ ____  ____   ___        __")')
      WRITE(NFIL,FMT='(T20,"\ \   / / ____/ ___||  _ \ / \ \      / /")')
      WRITE(NFIL,FMT='(T20," \ \ / /|  _| \___ \| |_) / _ \ \ /\ / / ")')
      WRITE(NFIL,FMT='(T20,"  \ V / | |___ ___) |  __/ ___ \ V  V /  ")')
      WRITE(NFIL,FMT='(T20,"   \_/  |_____|____/|_| /_/   \_\_/\_/   ")')
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='(T9,"Variational electronic structure with projector augmented waves")')
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='(T21,"VESPAW -- developed by Mohsen Sotoudeh")')
      WRITE(NFIL,FMT='(T4,"derived from CP-PAW, (C) Peter E. Bloechl, Clausthal Univ. of Technology")')
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='()')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE VESPAW_PRINTVERSIONBLOCK(NFIL)
!     **************************************************************************
!     ** PRINTS THE BORDERED "VESPAW VERSION INFORMATION" BLOCK TO UNIT NFIL. **
!     ** SHARED BY --VERSION (SCREEN) AND THE .PROT FILE HEADER (PUTHEADER,   **
!     ** src/paw_ioroutines.f90).                                             **
!     **************************************************************************
      USE VERSION_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
!     **************************************************************************
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='("VESPAW VERSION INFORMATION")')
      WRITE(NFIL,FMT='(79("-"))')
      IF(TRIM(VERTYP).EQ.'DEVELOPMENT VERSION') THEN
        CALL VERSION_PRINTDEV(NFIL,VERINF,VERREV,VERAUT,VERDAT)
      ELSE
!       == 'RELEASE VERSION'/'UNKNOWN VERSION' HAVE A DIFFERENT SHAPE  ==
!       == THAN THE GIT-DERIVED SENTENCES VERSION_PRINTDEV PARSES --   ==
!       == NOT EXERCISED BY THIS PROJECT YET, SO JUST PRINT THEM RAW.  ==
        WRITE(NFIL,FMT='(A)')TRIM(VERTYP)
        WRITE(NFIL,FMT='(A)')TRIM(VERINF)
        WRITE(NFIL,FMT='(A)')TRIM(VERREV)
        WRITE(NFIL,FMT='(A)')TRIM(VERAUT)
        WRITE(NFIL,FMT='(A)')TRIM(VERDAT)
      END IF
      WRITE(NFIL,FMT='(79("-"))')
      RETURN
      END
