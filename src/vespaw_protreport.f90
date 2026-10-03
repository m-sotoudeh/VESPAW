!
!.......................................................................
MODULE PROTREPORT_MODULE
!***********************************************************************
!**                                                                   **
!**  STRUCTURED RUN REPORT, WRITTEN TO THE .PROT FILE (FILE ID        **
!**  'REPORT'). LAYOUT AS THE VESPAW.OUT TEMPLATE OF                  **
!**  THE DEV NOTES: INPUT ECHO, 01 SYSTEM ... 09 PROVENANCE,       **
!**  PLUS THE K-POINT/BAND/OCCUPATION TABLE.                           **
!**  USER-FACING UNITS: ENERGIES IN EV, LENGTHS IN ANGSTROM, FORCES IN **
!**  EV/ANGSTROM, STRESS IN GPA. INTERNALLY ATOMIC UNITS ARE KEPT.     **
!**  THE DETAILED LEGACY PROTOCOL IS THE .LOG FILE (FILE ID 'PROT').   **
!**                                                                   **
!**  CALLED FROM THE DRIVER (VESPAW_DRIVER.F90):                       **
!**    PROTREPORT$HEADER          AFTER INITIALIZATION                 **
!**    PROTREPORT$ITERATION(NFI)  AFTER EVERY ELECTRONIC ITERATION     **
!**    PROTREPORT$CAPTUREEIG      IN PRINFO (MPI-COLLECTIVE, ALL TASKS)**
!**    PROTREPORT$FINISH(NFI,...) AT THE END OF THE RUN                **
!**                                                                   **
!***********************************************************************
IMPLICIT NONE
LOGICAL(4)            ,SAVE :: THEADER=.FALSE.  ! HEADER WRITTEN
LOGICAL(4)            ,SAVE :: TFIRSTROW=.TRUE.
INTEGER(4)            ,SAVE :: NITER=0
REAL(8)               ,SAVE :: EOLD=0.D0
REAL(8)               ,SAVE :: DELAST=0.D0
REAL(8)               ,SAVE :: RESLAST=0.D0
REAL(8)               ,SAVE :: TSTART=0.D0      ! WALL CLOCK AT HEADER [S]
REAL(8)               ,SAVE :: TROW=0.D0        ! WALL CLOCK AT LAST ROW [S]
REAL(8)               ,SAVE :: CPUSTART=0.D0
CHARACTER(64)         ,SAVE :: STARTSTAMP=' '
!== EIGENVALUES AND OCCUPATIONS OF THE LAST PRINTED STEP ===============
LOGICAL(4)            ,SAVE :: TEIG=.FALSE.
INTEGER(4)            ,SAVE :: NB=0,NKPT=0,NSPIN=0
REAL(8)   ,ALLOCATABLE,SAVE :: EIG(:,:,:)   !(NB,NKPT,NSPIN) [H]
REAL(8)   ,ALLOCATABLE,SAVE :: OCC(:,:,:)   !(NB,NKPT,NSPIN)
REAL(8)   ,ALLOCATABLE,SAVE :: XK(:,:)      !(3,NKPT)
REAL(8)   ,ALLOCATABLE,SAVE :: WKPT(:)      !(NKPT)
REAL(8)               ,SAVE :: EHOMO=0.D0,ELUMO=0.D0
LOGICAL(4)            ,SAVE :: TMETAL=.FALSE.
!== SPHERE-PROJECTED CHARGES AND MOMENTS PER ATOM AND ANGULAR MOMENTUM ==
LOGICAL(4)            ,SAVE :: TQ=.FALSE.
INTEGER(4)            ,SAVE :: NATQ=0,NDIMDQ=0
REAL(8)   ,ALLOCATABLE,SAVE :: QL(:,:,:)    !(0:3,NAT,NDIMD) S,P,D,F
CHARACTER(*),PARAMETER      :: RULE='(80("="))'
CHARACTER(*),PARAMETER      :: THIN='(80("-"))'
END MODULE PROTREPORT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT_UNIT(NFIL,TWRITE)
!     **************************************************************************
!     ** UNIT OF THE REPORT FILE; ONLY TASK 1 WRITES                          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NFIL
      LOGICAL(4),INTENT(OUT) :: TWRITE
      INTEGER(4)             :: NTASKS,THISTASK
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TWRITE=(THISTASK.EQ.1)
      NFIL=-1
      IF(TWRITE) CALL FILEHANDLER$UNIT('REPORT',NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT_SECTION(NFIL,NUM,TITLE)
      USE PROTREPORT_MODULE, ONLY : RULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NUM
      CHARACTER(*),INTENT(IN) :: TITLE
!     **************************************************************************
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT=RULE)
      WRITE(NFIL,FMT='(A,2X,A)') NUM,TITLE
      WRITE(NFIL,FMT=RULE)
      WRITE(NFIL,FMT='()')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT_KV(NFIL,KEY,VAL)
!     ** ONE "  KEY : VALUE" LINE WITH ALIGNED COLONS                         **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: KEY
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
      WRITE(NFIL,FMT='(2X,A,T34,": ",A)') KEY,TRIM(VAL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT_ECHOFILE(NFIL,ID)
!     **************************************************************************
!     ** COPIES THE FILE WITH FILEHANDLER ID 'ID' VERBATIM INTO THE REPORT    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(512)          :: NAME
      CHARACTER(512)          :: LINE
      INTEGER(4)              :: NFILIN,IOS
!     **************************************************************************
      CALL FILEHANDLER$FILENAME(ID,NAME)
      OPEN(NEWUNIT=NFILIN,FILE=TRIM(NAME),STATUS='OLD',ACTION='READ' &
     &    ,FORM='FORMATTED',IOSTAT=IOS)
      WRITE(NFIL,FMT='("--- ",A," ",54("-"))') TRIM(NAME)
      IF(IOS.NE.0) THEN
        WRITE(NFIL,FMT='("    (FILE COULD NOT BE READ)")')
        RETURN
      END IF
      DO
        READ(NFILIN,FMT='(A)',IOSTAT=IOS) LINE
        IF(IOS.NE.0) EXIT
        WRITE(NFIL,FMT='(4X,A)') TRIM(LINE)
      ENDDO
      CLOSE(NFILIN)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT_FUNCTIONALNAME(NAME)
!     **************************************************************************
!     ** NAME OF THE EXCHANGE-CORRELATION FUNCTIONAL OF THIS RUN, INCLUDING  **
!     ** THE ONSITE HYBRIDS (PBE0R, SCAN0R) AND ORBITAL-TAU SCAN             **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(OUT) :: NAME
      INTEGER(4)               :: ITYPE,IDC
      LOGICAL(4)               :: TLMTO,TTSCAN
      CHARACTER(32)            :: MODUS
!     **************************************************************************
      CALL DFT$GETI4('TYPE',ITYPE)
      SELECT CASE(ITYPE)
        CASE(1);     NAME='LDA (PERDEW-ZUNGER)'
        CASE(2);     NAME='LDA (PERDEW-WANG 92)'
        CASE(3);     NAME='LOCAL HARTREE-FOCK EXCHANGE'
        CASE(7);     NAME='B88 + PW92'
        CASE(8);     NAME='B88 + P86 (PW92 LOCAL)'
        CASE(81);    NAME='BP86'
        CASE(9);     NAME='B88 + PBE CORRELATION'
        CASE(10);    NAME='PBE'
        CASE(-10);   NAME='PBE (PRE-2024 REDUCED GRADIENT)'
        CASE(11);    NAME='RPBE'
        CASE(12);    NAME='PW91'
        CASE(-12);   NAME='PBE EXCHANGE + PW91 CORRELATION'
        CASE(13);    NAME='REVPBE'
        CASE(14);    NAME='PBESOL'
        CASE(20);    NAME='SCAN-L'
        CASE(21);    NAME='SCAN'
        CASE(5001);  NAME='BLYP'
        CASE DEFAULT
          WRITE(NAME,FMT='("DFT TYPE ",I0)') ITYPE
      END SELECT
      CALL POTENTIAL$GETL4('TSCANDELOC',TTSCAN)
      IF(TTSCAN) NAME='SCAN (orbital tau)'
!     == ONSITE HYBRIDS ON TOP OF THE GLOBAL FUNCTIONAL (NAMES OF 2026-10-03): ==
!     ==   SCAN0R   : ORBITAL-TAU SCAN EVERYWHERE + ONSITE EXACT EXCHANGE      ==
!     ==   PBE-SL0R : PBE EVERYWHERE, ONSITE EXCHANGE SCAN-L MIXED WITH EXACT  ==
!     ==              EXCHANGE (DCTYPE='SCANL'; CALLED SCAN0R BEFORE)          ==
!     ==   PBE0R    : PBE EVERYWHERE + ONSITE EXACT EXCHANGE                   ==
      CALL LMTO$GETL4('ON',TLMTO)
      IF(TLMTO) THEN
        CALL LMTO$GETCH('MODUS',MODUS)
        IF(INDEX(MODUS,'HYBRID').GT.0) THEN
          CALL LMTO$GETI4('DCTYPE',IDC)
          IF(TTSCAN.AND.(IDC.EQ.-1.OR.IDC.EQ.21)) THEN
            NAME='SCAN0r (orbital-tau SCAN + onsite exact exchange)'
          ELSE IF(ITYPE.EQ.10.AND.IDC.EQ.20) THEN
            NAME='PBE-SL0r (PBE + onsite SCAN-L/exact exchange)'
          ELSE IF(ITYPE.EQ.10.AND.(IDC.EQ.-1.OR.IDC.EQ.10)) THEN
            NAME='PBE0r (PBE + onsite exact exchange)'
          ELSE
            WRITE(MODUS,FMT='(I0)') IDC
            NAME=TRIM(NAME)//' + onsite hybrid (double counting type '//TRIM(MODUS)//')'
          END IF
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT$HEADER()
!     **************************************************************************
!     ** HEADER, INPUT ECHO, 01 SYSTEM, 02 MODEL, 03 NUMERICAL SETTINGS      **
!     **************************************************************************
      USE PROTREPORT_MODULE
      USE PERIODICTABLE_MODULE, ONLY : PERIODICTABLE$GET
      IMPLICIT NONE
      INTEGER(4)            :: NFIL
      LOGICAL(4)            :: TWRITE,TDYN,TCHK
      CHARACTER(256)        :: SVAL,CNTLNAME,STRCNAME
      CHARACTER(512)        :: CWD
      CHARACTER(128)        :: LINE
      REAL(8)               :: EV,ANGSTROM,RBAS(3,3),VOL,ALEN(3),ANG(3)
      REAL(8)               :: Z,SVAR,EPWPSI,EPWRHO,TEMP,TOTCHA,SPINVAL,EDIFF
      REAL(8)   ,ALLOCATABLE:: POS(:,:),FRAC(:,:)
      INTEGER(4)            :: NAT,IAT,I,J,NSP,ISP,NTASKS,THISTASK,NDIM,NSPINL
      INTEGER(4)            :: NKDIV(3),ISHIFT(3),NKPTL,IZ,NEL
      INTEGER(4)            :: NELEM
      CHARACTER(8)          :: ELEM(200),SYM
      INTEGER(4)            :: NOFELEM(200)
      CHARACTER(32)         :: NAME,OPT,MIX
      CHARACTER(96)         :: FUNC
      CHARACTER(8)          :: DATE
      CHARACTER(10)         :: TIME
      CHARACTER(5)          :: ZONE
      INTEGER(4)            :: VALUES(8)
      REAL(8)               :: RINV(3,3)
!     **************************************************************************
      IF(THEADER) RETURN
      THEADER=.TRUE.
      CALL TIMING_CLOCK(TSTART)
      TROW=TSTART
      CALL CPU_TIME(CPUSTART)
      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
      WRITE(STARTSTAMP,FMT='(I4.4,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2,1X,A)') &
     &      VALUES(1:3),VALUES(5:7),ZONE
      CALL PROTREPORT_UNIT(NFIL,TWRITE)
      IF(.NOT.TWRITE) RETURN
      CALL CONSTANTS('EV',EV)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
!
!     ==========================================================================
!     == BANNER AND RUN IDENTIFICATION                                        ==
!     ==========================================================================
      CALL VESPAW_PRINTBANNER(NFIL)
      CALL PROVENANCE$GET('VERSION',SVAL)
      WRITE(NFIL,FMT='(A)') 'VESPAW '//TRIM(SVAL)
      CALL PROVENANCE$GET('REVISION',SVAL);  CALL PROTREPORT_KV(NFIL,'Revision',SVAL)
      CALL PROVENANCE$GET('GITSTATE',SVAL);  CALL PROTREPORT_KV(NFIL,'Git state',SVAL)
      CALL PROVENANCE$GET('CPPAWBASE',SVAL); CALL PROTREPORT_KV(NFIL,'CP-PAW base',SVAL)
      CALL PROVENANCE$GET('COMPILER',SVAL);  CALL PROTREPORT_KV(NFIL,'Compiler',SVAL)
      CALL PROVENANCE$GET('MPI',SVAL);       CALL PROTREPORT_KV(NFIL,'MPI',SVAL)
      CALL PROVENANCE$GET('BLAS',SVAL);      CALL PROTREPORT_KV(NFIL,'BLAS',SVAL)
      CALL PROVENANCE$GET('FFT',SVAL);       CALL PROTREPORT_KV(NFIL,'FFT',SVAL)
      WRITE(NFIL,FMT='()')
      CALL PROTREPORT_KV(NFIL,'Started',STARTSTAMP)
      CALL FILEHANDLER$FILENAME('CNTL',CNTLNAME)
      CALL FILEHANDLER$FILENAME('STRC',STRCNAME)
      CALL PROTREPORT_KV(NFIL,'Control input',CNTLNAME)
      CALL PROTREPORT_KV(NFIL,'Structure input',STRCNAME)
      CALL GETCWD(CWD)
      CALL PROTREPORT_KV(NFIL,'Working directory',CWD)
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      SVAL=' '
      CALL GET_ENVIRONMENT_VARIABLE('OMP_NUM_THREADS',SVAL)
      IF(SVAL.EQ.' ') SVAL='(default)'
      WRITE(LINE,FMT='(I0," MPI task(s), OMP_NUM_THREADS=",A)') NTASKS,TRIM(SVAL)
      CALL PROTREPORT_KV(NFIL,'Parallelism',LINE)
      CALL PROTREPORT_KV(NFIL,'Units','energies eV, lengths A, forces eV/A, stress GPa')
!
!     ==========================================================================
!     == 00 INPUT (VERBATIM)                                                  ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'00','INPUT')
      CALL PROTREPORT_ECHOFILE(NFIL,'CNTL')
      WRITE(NFIL,FMT='()')
      CALL PROTREPORT_ECHOFILE(NFIL,'STRC')
!
!     ==========================================================================
!     == 01 SYSTEM                                                            ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'01','SYSTEM')
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(POS(3,NAT),FRAC(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,POS)
      CALL CELL$GETR8A('T(0)',9,RBAS)
!     == COMPOSITION ===========================================================
      NELEM=0
      NOFELEM(:)=0
      DO IAT=1,NAT
        CALL ATOMLIST$GETR8('Z',IAT,Z)
        IZ=NINT(Z)
        SYM='X'
        IF(IZ.GT.0) CALL PERIODICTABLE$GET(IZ,'SYMBOL',SYM)
        DO I=2,LEN_TRIM(SYM)
          IF(SYM(I:I).GE.'A'.AND.SYM(I:I).LE.'Z') SYM(I:I)=ACHAR(IACHAR(SYM(I:I))+32)
        ENDDO
        J=0
        DO I=1,NELEM
          IF(ELEM(I).EQ.SYM) J=I
        ENDDO
        IF(J.EQ.0) THEN
          NELEM=NELEM+1
          J=NELEM
          ELEM(J)=SYM
        END IF
        NOFELEM(J)=NOFELEM(J)+1
      ENDDO
      LINE=' '
      DO I=1,NELEM
        IF(NOFELEM(I).EQ.1) THEN
          LINE=TRIM(LINE)//TRIM(ELEM(I))
        ELSE
          WRITE(SVAL,FMT='(I0)') NOFELEM(I)
          LINE=TRIM(LINE)//TRIM(ELEM(I))//TRIM(SVAL)
        END IF
      ENDDO
      WRITE(NFIL,FMT='("Structure")')
      CALL PROTREPORT_KV(NFIL,'Formula (cell)',LINE)
      WRITE(LINE,FMT='(I0)') NAT
      CALL PROTREPORT_KV(NFIL,'Atoms',LINE)
!     == CELL ==================================================================
      VOL=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3))
      DO I=1,3
        ALEN(I)=SQRT(SUM(RBAS(:,I)**2))
      ENDDO
      ANG(1)=ACOS(DOT_PRODUCT(RBAS(:,2),RBAS(:,3))/(ALEN(2)*ALEN(3)))
      ANG(2)=ACOS(DOT_PRODUCT(RBAS(:,3),RBAS(:,1))/(ALEN(3)*ALEN(1)))
      ANG(3)=ACOS(DOT_PRODUCT(RBAS(:,1),RBAS(:,2))/(ALEN(1)*ALEN(2)))
      ANG=ANG*45.D0/ATAN(1.D0)
      WRITE(LINE,FMT='(3F12.6)') ALEN/ANGSTROM
      CALL PROTREPORT_KV(NFIL,'Lattice lengths a,b,c [A]',LINE)
      WRITE(LINE,FMT='(3F12.4)') ANG
      CALL PROTREPORT_KV(NFIL,'Angles alpha,beta,gamma',LINE)
      WRITE(LINE,FMT='(F14.6)') ABS(VOL)/ANGSTROM**3
      CALL PROTREPORT_KV(NFIL,'Volume [A^3]',LINE)
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='(2X,"Lattice vectors [A]")')
      DO I=1,3
        WRITE(NFIL,FMT='(4X,"a",I1,3F16.8)') I,RBAS(:,I)/ANGSTROM
      ENDDO
!     == POSITIONS =============================================================
      CALL LIB$INVERTR8(3,RBAS,RINV)
      FRAC=MATMUL(RINV,POS)
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='(2X,"Atom",T16,"Cartesian [A]",T56,"Fractional")')
      WRITE(NFIL,FMT='(2X,76("-"))')
      DO IAT=1,NAT
        CALL ATOMLIST$GETCH('NAME',IAT,NAME)
        WRITE(NFIL,FMT='(2X,A10,3F11.6,3X,3F9.5)') NAME,POS(:,IAT)/ANGSTROM,FRAC(:,IAT)
      ENDDO
      DEALLOCATE(POS,FRAC)
!     == ELECTRONIC STATE ======================================================
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='("Electronic state")')
      CALL WAVES$GETI4('NDIM',NDIM)
      CALL DYNOCC$GETI4('NSPIN',NSPINL)
      IF(NDIM.EQ.2) THEN
        SVAL='noncollinear (spinors)'
      ELSE IF(NSPINL.EQ.2) THEN
        SVAL='collinear spin polarized'
      ELSE
        SVAL='non-spin-polarized'
      END IF
      CALL PROTREPORT_KV(NFIL,'Spin',SVAL)
      CALL DYNOCC$GETR8('NEL',SVAR)
      WRITE(LINE,FMT='(F12.4)') SVAR
      CALL PROTREPORT_KV(NFIL,'Valence electrons',ADJUSTL(LINE))
      CALL DYNOCC$GETR8('TOTCHA',TOTCHA)
      WRITE(LINE,FMT='(F12.4," e")') TOTCHA
      CALL PROTREPORT_KV(NFIL,'Total charge',ADJUSTL(LINE))
      CALL DYNOCC$GETR8('SPIN',SPINVAL)
      WRITE(LINE,FMT='(F12.4," hbar")') SPINVAL
      CALL PROTREPORT_KV(NFIL,'Total spin (input)',ADJUSTL(LINE))
!
!     ==========================================================================
!     == 02 MODEL                                                             ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'02','MODEL')
      WRITE(NFIL,FMT='("Exchange-correlation")')
      CALL PROTREPORT_FUNCTIONALNAME(FUNC)
      CALL PROTREPORT_KV(NFIL,'Functional',FUNC)
      CALL LMTO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL LMTO$GETR8('HFWEIGHT',SVAR)
        WRITE(LINE,FMT='(F8.4)') SVAR
        CALL PROTREPORT_KV(NFIL,'Onsite exact-exchange fraction',ADJUSTL(LINE))
      END IF
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='("PAW setups")')
      CALL SETUP$GETI4('NSP',NSP)
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETCH('ID',NAME)
        CALL SETUP$GETCH('SETUPID',SVAL)
        CALL PROTREPORT_KV(NFIL,TRIM(NAME),SVAL)
        CALL SETUP$UNSELECT()
      ENDDO
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='("Basis")')
      CALL WAVES$GETR8('EPWPSI',EPWPSI)
      CALL WAVES$GETR8('EPWRHO',EPWRHO)
      WRITE(LINE,FMT='(F10.2," eV  (",F8.2," Ry)")') EPWPSI/EV,2.D0*EPWPSI
      CALL PROTREPORT_KV(NFIL,'Wave-function cutoff',ADJUSTL(LINE))
      WRITE(LINE,FMT='(F10.2," eV  (",F8.2," Ry)")') EPWRHO/EV,2.D0*EPWRHO
      CALL PROTREPORT_KV(NFIL,'Density cutoff',ADJUSTL(LINE))
      CALL KPOINTS$GETI4A('DIV',3,NKDIV)
      CALL KPOINTS$GETI4A('SHIFT',3,ISHIFT)
      CALL DYNOCC$GETI4('NKPT',NKPTL)
      WRITE(LINE,FMT='(I0," x ",I0," x ",I0,"   shift ",3I2,"   (",I0," k-points)")') &
     &      NKDIV,ISHIFT,NKPTL
      CALL PROTREPORT_KV(NFIL,'k-point mesh',LINE)
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='("Occupations")')
      CALL DYNOCC$GETL4('DYN',TDYN)
      IF(TDYN) THEN
        CALL DYNOCC$GETR8('TEMP',TEMP)
        WRITE(LINE,FMT='("Fermi-Dirac (Mermin), kT = ",F10.6," eV")') TEMP/EV
        CALL PROTREPORT_KV(NFIL,'Method',LINE)
      ELSE
        CALL PROTREPORT_KV(NFIL,'Method','fixed')
      END IF
!
!     ==========================================================================
!     == 03 NUMERICAL SETTINGS                                                ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'03','NUMERICAL SETTINGS')
      WRITE(NFIL,FMT='("Self-consistency")')
      CALL WAVES$GETCH('OPTIMIZER',OPT)
      CALL WAVES$GETCH('MIXER',MIX)
      IF(OPT.EQ.' ') OPT='DAMPED WAVE-FUNCTION DYNAMICS'
      CALL PROTREPORT_KV(NFIL,'Wave-function optimizer',OPT)
      IF(MIX.EQ.' ') MIX='none'
      CALL PROTREPORT_KV(NFIL,'Mixer',MIX)
      CALL WAVES$GETR8('EDIFF',EDIFF)
      IF(EDIFF.GT.0.D0) THEN
        WRITE(LINE,FMT='(ES10.2," eV")') EDIFF/EV
        CALL PROTREPORT_KV(NFIL,'Energy tolerance (EDIFF)',ADJUSTL(LINE))
      END IF
      WRITE(NFIL,FMT='()')
      CALL LIB$FLUSHFILE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT$ITERATION(NFI,ETOT,RESIDUAL)
!     **************************************************************************
!     ** ONE ROW OF THE SCF TABLE (SECTION 04)                                **
!     **************************************************************************
      USE PROTREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFI
      REAL(8)   ,INTENT(IN) :: ETOT      ! [H]
      REAL(8)   ,INTENT(IN) :: RESIDUAL
      INTEGER(4)            :: NFIL
      LOGICAL(4)            :: TWRITE
      REAL(8)               :: EV,TNOW
      CHARACTER(14)         :: DESTR
!     **************************************************************************
      CALL PROTREPORT_UNIT(NFIL,TWRITE)
      IF(.NOT.TWRITE) RETURN
      IF(.NOT.THEADER) RETURN
      CALL CONSTANTS('EV',EV)
      CALL TIMING_CLOCK(TNOW)
      IF(TFIRSTROW) THEN
        CALL PROTREPORT_SECTION(NFIL,'04','SCF CONVERGENCE')
        WRITE(NFIL,FMT='(A6,A22,A16,A13,A11)') 'Iter','Energy (eV)','dE (eV)' &
     &                                       ,'Residual','Time (s)'
        WRITE(NFIL,FMT=THIN)
        DESTR='---'
        DELAST=0.D0
      ELSE
        DELAST=ETOT-EOLD
        WRITE(DESTR,FMT='(ES14.3)') DELAST/EV
      END IF
      WRITE(NFIL,FMT='(I6,F22.8,A16,ES13.3,F11.2)') NFI,ETOT/EV,TRIM(DESTR) &
     &                                             ,RESIDUAL,TNOW-TROW
      CALL LIB$FLUSHFILE(NFIL)
      EOLD=ETOT
      RESLAST=RESIDUAL
      TROW=TNOW
      NITER=NITER+1
      TFIRSTROW=.FALSE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT$CAPTUREEIG()
!     **************************************************************************
!     ** KEEPS EIGENVALUES AND OCCUPATIONS OF THE CURRENT STEP FOR THE FINAL **
!     ** K-POINT TABLE. MPI-COLLECTIVE (WAVES$GETEIGDATA GATHERS): CALL ON   **
!     ** ALL TASKS, FROM PRINFO, WHILE THE EIGENVALUES ARE STILL VALID.      **
!     **************************************************************************
      USE PROTREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: NBL,NKPTL,NSPINL
!     **************************************************************************
      CALL DYNOCC$GETI4('NB',NBL)
      CALL DYNOCC$GETI4('NKPT',NKPTL)
      CALL DYNOCC$GETI4('NSPIN',NSPINL)
      IF(ALLOCATED(EIG)) DEALLOCATE(EIG,OCC,XK,WKPT)
      NB=NBL
      NKPT=NKPTL
      NSPIN=NSPINL
      ALLOCATE(EIG(NB,NKPT,NSPIN),OCC(NB,NKPT,NSPIN),XK(3,NKPT),WKPT(NKPT))
      CALL WAVES$GETEIGDATA(NB,NKPT,NSPIN,EIG,OCC,EHOMO,ELUMO,TMETAL)
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT)
!     == BAND EDGES FROM THE OCCUPATION OF EACH STATE (OCC/WEIGHT RELATIVE ===
!     == TO FULL OCCUPATION). WITH SMEARING, ANY PARTIALLY OCCUPIED STATE    ==
!     == (1-99 %) MEANS A METAL; OTHERWISE THE EDGES ARE AT HALF OCCUPATION. ==
      CALL PROTREPORT_BANDEDGES()
      TEIG=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT_BANDEDGES()
      USE PROTREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4) :: IB,IK,IS
      REAL(8)    :: FMAX,F
!     **************************************************************************
      FMAX=2.D0/REAL(NSPIN,KIND=8)
      EHOMO=-1.D+10
      ELUMO=+1.D+10
      TMETAL=.FALSE.
      DO IS=1,NSPIN
        DO IK=1,NKPT
          IF(WKPT(IK).LE.0.D0) CYCLE
          DO IB=1,NB
            F=OCC(IB,IK,IS)/(WKPT(IK)*FMAX)
            IF(F.GT.0.01D0.AND.F.LT.0.99D0) TMETAL=.TRUE.
            IF(F.GE.0.5D0) EHOMO=MAX(EHOMO,EIG(IB,IK,IS))
            IF(F.LT.0.5D0) ELUMO=MIN(ELUMO,EIG(IB,IK,IS))
          ENDDO
        ENDDO
      ENDDO
      IF(ELUMO.LE.EHOMO) TMETAL=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT$ATOMCHARGES(NAT,LMNXX,NDIMD,DENMAT)
!     **************************************************************************
!     ** CHARGE AND MAGNETIZATION OF EACH ATOM, RESOLVED IN ANGULAR MOMENTUM, **
!     ** FROM THE ONE-CENTER DENSITY MATRIX AND THE ALL-ELECTRON PARTIAL     **
!     ** WAVES INSIDE THE ATOMIC SPHERE (RADIUS RAD OF THE SETUP, AS FOR THE **
!     ** PDOS):  Q_L = SUM_{I,J IN L} D_IJ <PHI_I|PHI_J>_SPHERE.             **
!     ** COMPONENT 1 IS THE CHARGE, 2 THE SPIN (NDIMD=2) OR 2..4 THE X,Y,Z   **
!     ** SPIN DENSITY (NDIMD=4). CALLED FROM WAVES$ETOT ON EVERY ITERATION;  **
!     ** THE LAST CALL IS REPORTED.                                          **
!     **************************************************************************
      USE PROTREPORT_MODULE, ONLY : TQ,NATQ,NDIMDQ,QL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT,LMNXX,NDIMD
      COMPLEX(8),INTENT(IN) :: DENMAT(LMNXX,LMNXX,NDIMD,NAT)
      INTEGER(4)            :: NTASKS,THISTASK
      INTEGER(4)            :: IAT,ISP,GID,NR,LNX,LN1,LN2,L,M,LMN1,LMN2,IDIM,IR
      INTEGER(4),ALLOCATABLE:: LOX(:),LMN0(:)
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:),R(:),WORK(:),CUMUL(:)
      REAL(8)               :: RAD,OV
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      IF(.NOT.ALLOCATED(QL)) ALLOCATE(QL(0:3,NAT,NDIMD))
      IF(SIZE(QL,2).NE.NAT.OR.SIZE(QL,3).NE.NDIMD) THEN
        DEALLOCATE(QL)
        ALLOCATE(QL(0:3,NAT,NDIMD))
      END IF
      NATQ=NAT
      NDIMDQ=NDIMD
      QL(:,:,:)=0.D0
      DO IAT=1,NAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        CALL SETUP$GETI4('LNX',LNX)
        CALL SETUP$GETR8('RAD',RAD)
        ALLOCATE(LOX(LNX),LMN0(LNX),AEPHI(NR,LNX),R(NR),WORK(NR),CUMUL(NR))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
        CALL SETUP$UNSELECT()
        CALL RADIAL$R(GID,NR,R)
        LMN0(1)=0
        DO LN1=2,LNX
          LMN0(LN1)=LMN0(LN1-1)+2*LOX(LN1-1)+1
        ENDDO
        DO LN1=1,LNX
          L=LOX(LN1)
          IF(L.GT.3) CYCLE
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            WORK(:)=AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRATE(GID,NR,WORK,CUMUL)
            CALL RADIAL$VALUE(GID,NR,CUMUL,RAD,OV)
            DO M=1,2*L+1
              LMN1=LMN0(LN1)+M
              LMN2=LMN0(LN2)+M
              DO IDIM=1,NDIMD
                QL(L,IAT,IDIM)=QL(L,IAT,IDIM)+REAL(DENMAT(LMN1,LMN2,IDIM,IAT),KIND=8)*OV
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(LOX,LMN0,AEPHI,R,WORK,CUMUL)
      ENDDO
      TQ=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT_QTABLE(NFIL,TITLE,IDIM)
!     ** TABLE OF THE SPHERE-PROJECTED CHARGE OR MAGNETIZATION COMPONENT IDIM **
      USE PROTREPORT_MODULE, ONLY : NATQ,QL
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: TITLE
      INTEGER(4)  ,INTENT(IN) :: IDIM
      INTEGER(4)              :: IAT
      CHARACTER(32)           :: NAME
      REAL(8)                 :: SUML(0:3)
!     **************************************************************************
      WRITE(NFIL,FMT='(A)') TITLE
      WRITE(NFIL,FMT='(2X,A10,5A11)') 'Atom','s','p','d','f','total'
      WRITE(NFIL,FMT='(2X,65("-"))')
      SUML(:)=0.D0
      DO IAT=1,NATQ
        CALL ATOMLIST$GETCH('NAME',IAT,NAME)
        WRITE(NFIL,FMT='(2X,A10,5F11.4)') NAME,QL(:,IAT,IDIM),SUM(QL(:,IAT,IDIM))
        SUML(:)=SUML(:)+QL(:,IAT,IDIM)
      ENDDO
      WRITE(NFIL,FMT='(2X,65("-"))')
      WRITE(NFIL,FMT='(2X,A10,5F11.4)') 'Sum',SUML,SUM(SUML)
      WRITE(NFIL,FMT='()')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PROTREPORT$FINISH(NFI,TSTOPPED,NATF,FORCE,TFORCE,STRESS,TSTRESS,RBAS)
!     **************************************************************************
!     ** SECTIONS 04 (STATUS) TO 09 AND THE SUMMARY                          **
!     **************************************************************************
      USE PROTREPORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFI
      LOGICAL(4),INTENT(IN) :: TSTOPPED
      INTEGER(4),INTENT(IN) :: NATF
      REAL(8)   ,INTENT(IN) :: FORCE(3,NATF) ! [H/BOHR]
      LOGICAL(4),INTENT(IN) :: TFORCE
      REAL(8)   ,INTENT(IN) :: STRESS(3,3)  ! -DE/DEPSILON [H]
      LOGICAL(4),INTENT(IN) :: TSTRESS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4)            :: NFIL,I,IKPT,IB,IS,IAT,NAT,NACC,NREJ,NREA
      LOGICAL(4)            :: TWRITE,TDYN
      REAL(8)               :: EV,ANGSTROM,GPA,ETOT,TNOW,CPUNOW,SVAR,EF
      REAL(8)               :: VOL,SIGMA(3,3),FMAX,FRMS,FABS,TOTWALL
      REAL(8)               :: QUP,QDN
      CHARACTER(64)         :: REASON,STATUS
      CHARACTER(128)        :: LINE
      CHARACTER(32)         :: NAME
      CHARACTER(512)        :: FNAME
      CHARACTER(40)         :: EID
      CHARACTER(32)         :: CNAME
      REAL(8)               :: EVAL,WALL
      INTEGER(4)            :: NE,NENTRY
      REAL(8)   ,ALLOCATABLE:: POSF(:,:)
      REAL(8)               :: DRIFT(3),ETS
!     **************************************************************************
      CALL PROTREPORT_UNIT(NFIL,TWRITE)
      IF(.NOT.TWRITE) RETURN
      IF(.NOT.THEADER) RETURN
      CALL CONSTANTS('EV',EV)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL CONSTANTS('PASCAL',GPA)
      GPA=1.D+9*GPA
      CALL TIMING_CLOCK(TNOW)
      CALL CPU_TIME(CPUNOW)
      CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)
!
!     == STATUS OF THE SCF ======================================================
      CALL STOPIT$GETCH('REASON',REASON)
      IF(INDEX(REASON,'CONVERGED').GT.0) THEN
        STATUS='CONVERGED'
      ELSE IF(TSTOPPED) THEN
        STATUS='STOPPED'
        IF(REASON.EQ.' ') REASON='STOP SIGNAL'
      ELSE
        STATUS='NOT CONVERGED'
        REASON='MAXIMUM NUMBER OF ITERATIONS (NSTEP) REACHED'
      END IF
      WRITE(NFIL,FMT=THIN)
      WRITE(NFIL,FMT='("SCF STATUS: ",A)') TRIM(STATUS)
      CALL PROTREPORT_KV(NFIL,'Reason',REASON)
      WRITE(LINE,FMT='(I0)') NITER
      CALL PROTREPORT_KV(NFIL,'Iterations',LINE)
      WRITE(LINE,FMT='(ES12.3," eV")') DELAST/EV
      CALL PROTREPORT_KV(NFIL,'Last energy change',ADJUSTL(LINE))
      CALL MIXER$GETI4('NACCEPT',NACC)
      CALL MIXER$GETI4('NREJECT',NREJ)
      CALL MIXER$GETI4('NREANCHOR',NREA)
      IF(NACC+NREJ.GT.0) THEN
        WRITE(LINE,FMT='(I0," accepted, ",I0," rejected, ",I0," re-anchored")') &
     &        NACC,NREJ,NREA
        CALL PROTREPORT_KV(NFIL,'Variational steps',LINE)
      END IF
!
!     ==========================================================================
!     == 05 FINAL RESULTS                                                     ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'05','FINAL RESULTS')
      WRITE(NFIL,FMT='("Total energy")')
      CALL DYNOCC$GETL4('DYN',TDYN)
      CALL ENERGYLIST$RETURN('OCCUPATIONAL ENTROPY TERM (-TS)',ETS)
      IF(TDYN) THEN
        WRITE(LINE,FMT='(F22.8," eV")') ETOT/EV
        CALL PROTREPORT_KV(NFIL,'Free energy E-TS',ADJUSTL(LINE))
        WRITE(LINE,FMT='(F22.8," eV")') (ETOT-ETS)/EV
        CALL PROTREPORT_KV(NFIL,'Energy without entropy E',ADJUSTL(LINE))
        WRITE(LINE,FMT='(F22.8," eV")') (ETOT-0.5D0*ETS)/EV
        CALL PROTREPORT_KV(NFIL,'Extrapolated E(T->0)',ADJUSTL(LINE))
      ELSE
        WRITE(LINE,FMT='(F22.8," eV")') ETOT/EV
        CALL PROTREPORT_KV(NFIL,'Total energy',ADJUSTL(LINE))
      END IF
      WRITE(LINE,FMT='(F22.10," Hartree")') ETOT
      CALL PROTREPORT_KV(NFIL,'  (in atomic units)',ADJUSTL(LINE))
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='("Energy components [eV]")')
      CALL ENERGYLIST$GETNUMBER(NE)
      DO I=1,NE
        CALL ENERGYLIST$GETENTRY(I,EID,EVAL)
        IF(ABS(EVAL).LT.1.D-12) CYCLE
        IF(EID.EQ.'TOTAL ENERGY') CYCLE
        WRITE(NFIL,FMT='(2X,A40,F20.8)') EID,EVAL/EV
      ENDDO
!     == BAND EDGES AND FERMI LEVEL ============================================
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT='("Electronic properties")')
      IF(TDYN) THEN
        CALL DYNOCC$GETR8('EFERMI',EF)
        WRITE(LINE,FMT='(F14.6," eV")') EF/EV
        CALL PROTREPORT_KV(NFIL,'Fermi level',ADJUSTL(LINE))
      END IF
      IF(TEIG) THEN
        WRITE(LINE,FMT='(F14.6," eV")') EHOMO/EV
        CALL PROTREPORT_KV(NFIL,'Highest occupied state',ADJUSTL(LINE))
        WRITE(LINE,FMT='(F14.6," eV")') ELUMO/EV
        CALL PROTREPORT_KV(NFIL,'Lowest unoccupied state',ADJUSTL(LINE))
        IF(TMETAL) THEN
          CALL PROTREPORT_KV(NFIL,'Band gap','none (metallic: partially occupied states)')
        ELSE
          WRITE(LINE,FMT='(F14.6," eV")') (ELUMO-EHOMO)/EV
          CALL PROTREPORT_KV(NFIL,'Band gap',ADJUSTL(LINE))
        END IF
      ELSE
        CALL PROTREPORT_KV(NFIL,'Band edges','not available for this run')
      END IF
!
!     ==========================================================================
!     == K-POINTS, BAND ENERGIES AND OCCUPATIONS                              ==
!     ==========================================================================
      IF(TEIG) THEN
        CALL PROTREPORT_SECTION(NFIL,'05b','K-POINTS, BAND ENERGIES [eV] AND OCCUPATIONS')
        DO IKPT=1,NKPT
          WRITE(NFIL,FMT='("k-point ",I0,":  k =",3F10.5,"  (reciprocal lattice)'// &
     &          '  weight =",F10.6)') IKPT,XK(:,IKPT),WKPT(IKPT)
          IF(NSPIN.EQ.2) THEN
            WRITE(NFIL,FMT='(2X,A6,2A16,4X,2A16)') 'Band','E up','occ up','E down','occ down'
          ELSE
            WRITE(NFIL,FMT='(2X,A6,2A16)') 'Band','E','occupation'
          END IF
          DO IB=1,NB
!           == OCC CONTAINS THE K-POINT WEIGHT; PRINT THE OCCUPATION OF THE STATE
            IF(NSPIN.EQ.2) THEN
              WRITE(NFIL,FMT='(2X,I6,F16.6,F16.6,4X,F16.6,F16.6)') IB &
     &              ,EIG(IB,IKPT,1)/EV,OCC(IB,IKPT,1)/WKPT(IKPT) &
     &              ,EIG(IB,IKPT,2)/EV,OCC(IB,IKPT,2)/WKPT(IKPT)
            ELSE
              WRITE(NFIL,FMT='(2X,I6,F16.6,F16.6)') IB,EIG(IB,IKPT,1)/EV &
     &              ,OCC(IB,IKPT,1)/WKPT(IKPT)
            END IF
          ENDDO
          WRITE(NFIL,FMT='()')
        ENDDO
      END IF
!
!     ==========================================================================
!     == 06 SPIN / MAGNETISM                                                  ==
!     ==========================================================================
      IF(TQ.OR.(TEIG.AND.NSPIN.EQ.2)) THEN
        CALL PROTREPORT_SECTION(NFIL,'06','ATOMIC CHARGES AND MAGNETIC MOMENTS')
        IF(TQ) THEN
          WRITE(NFIL,FMT='("Projected onto the all-electron partial waves inside the atomic"' &
     &                //',/"spheres (radius RAD of each setup, as for the PDOS). The sum"' &
     &                //',/"misses the interstitial region.")')
          WRITE(NFIL,FMT='()')
          CALL PROTREPORT_QTABLE(NFIL,'Charge [e]',1)
          IF(NDIMDQ.EQ.2) THEN
            CALL PROTREPORT_QTABLE(NFIL,'Magnetic moment [mu_B]',2)
          ELSE IF(NDIMDQ.EQ.4) THEN
            CALL PROTREPORT_QTABLE(NFIL,'Magnetic moment, x component [mu_B]',2)
            CALL PROTREPORT_QTABLE(NFIL,'Magnetic moment, y component [mu_B]',3)
            CALL PROTREPORT_QTABLE(NFIL,'Magnetic moment, z component [mu_B]',4)
          END IF
        END IF
        IF(TEIG.AND.NSPIN.EQ.2) THEN
          QUP=SUM(OCC(:,:,1))
          QDN=SUM(OCC(:,:,2))
          WRITE(LINE,FMT='(F12.6)') QUP
          CALL PROTREPORT_KV(NFIL,'Spin-up electrons (cell)',ADJUSTL(LINE))
          WRITE(LINE,FMT='(F12.6)') QDN
          CALL PROTREPORT_KV(NFIL,'Spin-down electrons (cell)',ADJUSTL(LINE))
          WRITE(LINE,FMT='(F12.6," mu_B")') QUP-QDN
          CALL PROTREPORT_KV(NFIL,'Total magnetic moment (cell)',ADJUSTL(LINE))
        END IF
      END IF
!
!     ==========================================================================
!     == 07 FORCES AND STRESS                                                 ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'07','FORCES AND STRESS')
      IF(TFORCE) THEN
        NAT=NATF
        ALLOCATE(POSF(3,NAT))
        CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,POSF)
        WRITE(NFIL,FMT='(2X,A8,A33,A39)') 'Atom','Position (A)','Force (eV/A)'
        WRITE(NFIL,FMT='(2X,A8,6A12)') ' ','x','y','z','Fx','Fy','Fz'
        WRITE(NFIL,FMT='(2X,80("-"))')
        FMAX=0.D0
        FRMS=0.D0
        DRIFT(:)=0.D0
        DO IAT=1,NAT
          CALL ATOMLIST$GETCH('NAME',IAT,NAME)
          WRITE(NFIL,FMT='(2X,A8,3F12.6,3F12.6)') NAME,POSF(:,IAT)/ANGSTROM &
     &                                         ,FORCE(:,IAT)/(EV/ANGSTROM)
          DRIFT(:)=DRIFT(:)+FORCE(:,IAT)
          FABS=SQRT(SUM(FORCE(:,IAT)**2))
          FMAX=MAX(FMAX,FABS)
          FRMS=FRMS+FABS**2
        ENDDO
        FRMS=SQRT(FRMS/REAL(MAX(NAT,1),KIND=8))
        WRITE(NFIL,FMT='(2X,80("-"))')
        WRITE(NFIL,FMT='(2X,A8,36X,3F12.6)') 'Sum',DRIFT(:)/(EV/ANGSTROM)
        DEALLOCATE(POSF)
        WRITE(NFIL,FMT='()')
        WRITE(LINE,FMT='(F14.8," eV/A")') FMAX/(EV/ANGSTROM)
        CALL PROTREPORT_KV(NFIL,'Maximum force',ADJUSTL(LINE))
        WRITE(LINE,FMT='(F14.8," eV/A")') FRMS/(EV/ANGSTROM)
        CALL PROTREPORT_KV(NFIL,'RMS force',ADJUSTL(LINE))
      ELSE
        CALL PROTREPORT_KV(NFIL,'Forces','not available for this run')
      END IF
      IF(TSTRESS) THEN
        VOL=ABS(RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &         +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &         +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)))
        SIGMA=STRESS/VOL      ! STRESS=-DE/DEPSILON, SO SIGMA=-(1/V)DE/DEPSILON
        WRITE(NFIL,FMT='()')
        WRITE(NFIL,FMT='(2X,"Stress tensor [GPa] (positive: cell wants to expand)")')
        DO I=1,3
          WRITE(NFIL,FMT='(4X,3F16.6)') SIGMA(I,:)/GPA
        ENDDO
        SVAR=(SIGMA(1,1)+SIGMA(2,2)+SIGMA(3,3))/3.D0
        WRITE(LINE,FMT='(F14.6," GPa")') SVAR/GPA
        CALL PROTREPORT_KV(NFIL,'Pressure (tr sigma/3)',ADJUSTL(LINE))
      ELSE
        CALL PROTREPORT_KV(NFIL,'Stress','not computed (STRESS=T in !GENERIC)')
      END IF
!
!     ==========================================================================
!     == 08 PERFORMANCE                                                       ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'08','PERFORMANCE')
      WRITE(LINE,FMT='(F12.2," s")') TNOW-TSTART
      CALL PROTREPORT_KV(NFIL,'Wall time',ADJUSTL(LINE))
      WRITE(LINE,FMT='(F12.2," s (task 1)")') CPUNOW-CPUSTART
      CALL PROTREPORT_KV(NFIL,'CPU time',ADJUSTL(LINE))
      IF(NITER.GT.0) THEN
        WRITE(LINE,FMT='(F12.3," s")') (TNOW-TSTART)/REAL(NITER,KIND=8)
        CALL PROTREPORT_KV(NFIL,'Time per iteration',ADJUSTL(LINE))
      END IF
!     == MAIN ROUTINES (INCLUSIVE WALL TIME, TASK 1, >= 2 % OF THE TOTAL) =====
      TOTWALL=0.D0
      CALL TIMING$GETNUMBER(NENTRY)
      DO I=1,NENTRY
        CALL TIMING$GETENTRY(I,CNAME,WALL)
        IF(CNAME.EQ.'TOTALCLOCK') TOTWALL=WALL
      ENDDO
      IF(TOTWALL.GT.0.D0) THEN
        WRITE(NFIL,FMT='()')
        WRITE(NFIL,FMT='(2X,"Main routines (inclusive wall time, task 1)")')
        DO I=1,NENTRY
          CALL TIMING$GETENTRY(I,CNAME,WALL)
          IF(CNAME.EQ.'TOTALCLOCK') CYCLE
          IF(WALL.LT.0.02D0*TOTWALL) CYCLE
          WRITE(NFIL,FMT='(4X,A32,F12.2," s",F8.1," %")') CNAME,WALL,100.D0*WALL/TOTWALL
        ENDDO
      END IF
!
!     ==========================================================================
!     == 09 FILES AND PROVENANCE                                              ==
!     ==========================================================================
      CALL PROTREPORT_SECTION(NFIL,'09','FILES AND PROVENANCE')
      CALL FILEHANDLER$FILENAME('REPORT',FNAME)
      CALL PROTREPORT_KV(NFIL,'Report',FNAME)
      CALL FILEHANDLER$FILENAME('PROT',FNAME)
      CALL PROTREPORT_KV(NFIL,'Detailed log',FNAME)
      CALL FILEHANDLER$FILENAME('STRC_OUT',FNAME)
      CALL PROTREPORT_KV(NFIL,'Final structure',FNAME)
      CALL FILEHANDLER$FILENAME('RESTART_OUT',FNAME)
      CALL PROTREPORT_KV(NFIL,'Restart file',FNAME)
      CALL PROVENANCE$GET('REVISION',LINE)
      CALL PROTREPORT_KV(NFIL,'VESPAW revision',LINE)
      CALL PROVENANCE$GET('BUILT',LINE)
      CALL PROTREPORT_KV(NFIL,'Built',LINE)
      LINE=' '
      CALL GET_ENVIRONMENT_VARIABLE('USER',LINE)
      CALL PROTREPORT_KV(NFIL,'User',LINE)
      LINE=' '
      CALL HOSTNM(LINE)
      CALL PROTREPORT_KV(NFIL,'Host',LINE)
!
!     ==========================================================================
!     == SUMMARY                                                              ==
!     ==========================================================================
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT=RULE)
      WRITE(NFIL,FMT='(T31,"CALCULATION SUMMARY")')
      WRITE(NFIL,FMT=RULE)
      WRITE(NFIL,FMT='()')
      CALL PROTREPORT_KV(NFIL,'SCF',STATUS)
      WRITE(LINE,FMT='(I0)') NITER
      CALL PROTREPORT_KV(NFIL,'Iterations',LINE)
      WRITE(LINE,FMT='(F22.8," eV")') ETOT/EV
      CALL PROTREPORT_KV(NFIL,'Total energy',ADJUSTL(LINE))
      IF(TEIG.AND..NOT.TMETAL) THEN
        WRITE(LINE,FMT='(F14.6," eV")') (ELUMO-EHOMO)/EV
        CALL PROTREPORT_KV(NFIL,'Band gap',ADJUSTL(LINE))
      END IF
      IF(TFORCE) THEN
        WRITE(LINE,FMT='(F14.8," eV/A")') FMAX/(EV/ANGSTROM)
        CALL PROTREPORT_KV(NFIL,'Maximum force',ADJUSTL(LINE))
      END IF
      WRITE(LINE,FMT='(F12.2," s")') TNOW-TSTART
      CALL PROTREPORT_KV(NFIL,'Wall time',ADJUSTL(LINE))
      WRITE(NFIL,FMT='()')
      WRITE(NFIL,FMT=RULE)
      WRITE(NFIL,FMT='(T33,"VESPAW FINISHED")')
      WRITE(NFIL,FMT=RULE)
      CALL LIB$FLUSHFILE(NFIL)
      RETURN
      END
