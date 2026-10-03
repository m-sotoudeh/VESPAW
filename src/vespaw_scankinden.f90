!.........1.........2.........3.........4.........5.........6.........7.........8
!***********************************************************************
!**                                                                    **
!**  INTERSTITIAL (PLANE-WAVE) KINETIC-ENERGY-DENSITY CONSTRUCTION FOR **
!**  TRUE SCAN-IN-PAW. MIRRORS WAVES_DENSITY (paw_waves1.f90:4611,     **
!**  read before writing this) AS CLOSELY AS POSSIBLE: SAME SIGNATURE  **
!**  SHAPE, SAME PLANEWAVE$FFT('GTOR',NBH*NDIM,NGL,...) BATCHED-FFT    **
!**  CONVENTION, SAME NDIM=1/NDIM=2 BRANCH STRUCTURE -- BUT BUILDS     **
!**  D(PSI_N)/D(X_D) VIA AN i*G_D MULTIPLICATION IN G-SPACE BEFORE THE **
!**  TRANSFORM, ONCE PER CARTESIAN DIRECTION, SUMMED AFTERWARD (SAME   **
!**  STRUCTURE AS VASP's TAU_PW, ALREADY READ IN THIS PROJECT'S        **
!**  SCAN0r notes, AND AS THE ALREADY-VALIDATED STANDALONE  **
!**  PROTOTYPE SCAN0r notes                   **
!**  verify_interstitial_tau_realpsi.f90, WHICH USED REAL PSI(G) FROM  **
!**  AN INSTRUMENTED RERUN OF THIS EXACT WAVES_VPSI CALL SITE.         **
!**                                                                    **
!**  TAU(R) = (1/2) SUM_N F_N SUM_D |D(PSI_N)/D(X_D)(R)|**2 -- STANDARD**
!**  DEFINITION, SAME CONVENTION AS scan_x_deriv_module.f90's          **
!**  TAU_UNIF=CF*N**(5/3) (10_hamiltonian_tau_spin_and_realpotential/).**
!**                                                                    **
!**  SCOPE CUT, EXPLICIT (MATCHING 10_'s COLLINEAR-ONLY CHOICE):       **
!**  NDIM=1 (SPIN-RESTRICTED OR COLLINEAR-SPIN-POLARIZED, THE ONLY     **
!**  CASE ANY ACTUAL TARGET CALCULATION IN THIS PROJECT USES) IS FULLY **
!**  IMPLEMENTED. TINV (GAMMA-POINT "SUPER WAVE FUNCTION" REAL-VALUED  **
!**  TRICK, WAVES_DENSITY's OWN TINV BRANCH) AND NDIM=2 (NONCOLLINEAR) **
!**  ARE NOT IMPLEMENTED -- BOTH RAISE ERROR$STOP RATHER THAN SILENTLY **
!**  COMPUTE SOMETHING WRONG, MATCHING HOW WAVES_DENSITY ITSELF        **
!**  REFUSES INCOMPATIBLE COMBINATIONS (E.G. TINV+NDIM=2).             **
!**                                                                    **
!**  STATUS: NOT YET CALLED FROM ANYWHERE IN THE PRODUCTION CODE (SAME **
!**  "PORT, NOT YET WIRE INTO THE SCF LOOP" SCOPE AS paw_scantau.f90). **
!***********************************************************************
!
!     ..................................................................
      SUBROUTINE WAVES_KINETICDENSITY(NGL,NRL,NDIM,NB,NBH,F,PSIOFG,TAU)
!     ******************************************************************
!     **  INTERSTITIAL KINETIC-ENERGY DENSITY TAU IN REAL SPACE, SAME  **
!     **  SIGNATURE SHAPE AS WAVES_DENSITY (RHO -> TAU).               **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGL         ! MAX # PLANE WAVES
      INTEGER(4),INTENT(IN)  :: NRL         ! # R-SPACE POINTS
      INTEGER(4),INTENT(IN)  :: NDIM        ! DIMENSION OF THE WAVE FUNCTION
      INTEGER(4),INTENT(IN)  :: NB          ! # BANDS
      INTEGER(4),INTENT(IN)  :: NBH
      REAL(8)   ,INTENT(IN)  :: F(NB)       ! OCCUPATIONS
      COMPLEX(8),INTENT(IN)  :: PSIOFG(NGL,NDIM,NBH)
      REAL(8)   ,INTENT(OUT) :: TAU(NRL,NDIM**2) ! KINETIC-ENERGY DENSITY IN R-SPACE
      COMPLEX(8),ALLOCATABLE :: PSIOFG_D(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSIOFR_D(:,:,:)
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)
      LOGICAL(4)              :: TINV
      INTEGER(4)              :: IBH,IR,IDIR,IG,IDIM
      REAL(8)                 :: F1
      REAL(8)                 :: RE,IM
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
!     ******************************************************************
                               CALL TRACE$PUSH('WAVES_KINETICDENSITY')
!
!     ==================================================================
!     ==  SCOPE CHECKS -- SEE FILE HEADER                             ==
!     ==================================================================
      CALL PLANEWAVE$GETL4('TINV',TINV)
      IF(TINV) THEN
        CALL ERROR$MSG('TINV (SUPER WAVE FUNCTIONS) NOT YET SUPPORTED')
        CALL ERROR$STOP('WAVES_KINETICDENSITY')
      END IF
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('ONLY NDIM=1 (COLLINEAR) IS SUPPORTED SO FAR')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$STOP('WAVES_KINETICDENSITY')
      END IF
      IF(NBH.NE.NB) THEN
        CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
        CALL ERROR$STOP('WAVES_KINETICDENSITY')
      END IF
!
!     ==================================================================
!     ==  GET G-VECTORS FOR THE CURRENTLY-SELECTED PLANE-WAVE BASIS   ==
!     ==================================================================
      ALLOCATE(GVEC(3,NGL))
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
!
!     ==================================================================
!     ==  FOR EACH CARTESIAN DIRECTION: D(PSI)/D(X_D) VIA i*G_D IN    ==
!     ==  G-SPACE, THEN GTOR, THEN ACCUMULATE |D(PSI)/D(X_D)|**2      ==
!     ==================================================================
      ALLOCATE(PSIOFG_D(NGL,NDIM,NBH))
      ALLOCATE(PSIOFR_D(NRL,NDIM,NBH))
      TAU(:,:)=0.D0
      DO IDIR=1,3
        DO IBH=1,NBH
          DO IDIM=1,NDIM
            DO IG=1,NGL
              PSIOFG_D(IG,IDIM,IBH)=CI*GVEC(IDIR,IG)*PSIOFG(IG,IDIM,IBH)
            ENDDO
          ENDDO
        ENDDO
        CALL PLANEWAVE$FFT('GTOR',NBH*NDIM,NGL,PSIOFG_D,NRL,PSIOFR_D)
        DO IBH=1,NB
          F1=F(IBH)
          IF(F1.EQ.0.D0) CYCLE
          DO IR=1,NRL
            RE= REAL(PSIOFR_D(IR,1,IBH),KIND=8)
            IM=AIMAG(PSIOFR_D(IR,1,IBH))
            TAU(IR,1)=TAU(IR,1)+F1*(RE**2+IM**2)
          ENDDO
        ENDDO
      ENDDO
      TAU(:,1)=0.5D0*TAU(:,1)
      DEALLOCATE(GVEC,PSIOFG_D,PSIOFR_D)
                               CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_TAUHPSI(NGL,NRL,NDIM,NBH,VTAU,PSIOFG,HTAUPSI)
!     ******************************************************************
!     **  ADDITIONAL GENERALIZED-KOHN-SHAM HAMILTONIAN TERM REQUIRED   **
!     **  BY A REAL (ORBITAL-TAU-BASED) META-GGA IN THE INTERSTITIAL   **
!     **  REGION:                                                      **
!     **                                                                **
!     **     H_TAU|PSI> = -(1/2) DIV( VTAU(R) GRAD(PSI(R)) )           **
!     **                                                                **
!     **  DERIVATION: E_TAU=INTEGRAL(VTAU*TAU)DR, TAU=(1/2)|GRAD(PSI)|^2**
!     **  (PER BAND, OCCUPATION APPLIED BY THE CALLER LIKE WAVES_VPSI), **
!     **  SO D(E_TAU)/D(PSI*(R)) = -(1/2) DIV(VTAU(R) GRAD(PSI(R))) BY  **
!     **  THE STANDARD FUNCTIONAL-DERIVATIVE/INTEGRATION-BY-PARTS       **
!     **  ARGUMENT (BOUNDARY TERM VANISHES, PERIODIC CELL). VALIDATED   **
!     **  STANDALONE IN PYTHON (FFT-BASED, INDEPENDENT OF THIS FORTRAN  **
!     **  AND OF THE PLANEWAVE OBJECT) AGAINST A DIRECT FINITE          **
!     **  DIFFERENCE OF E_TAU UNDER A RANDOM WAVEFUNCTION PERTURBATION, **
!     **  9.5E-9 RELATIVE AGREEMENT (SCAN0r notes                            **
!     **  SCAN0r notes, a test script) **
!     **  -- CONFIRMS THE SIGN AND FACTOR OF -1/2 BEFORE ANY FORTRAN     **
!     **  EXISTED. THIS ROUTINE IS THE FORTRAN/PLANEWAVE-OBJECT         **
!     **  REALIZATION OF THAT SAME VALIDATED ALGORITHM.                 **
!     **                                                                **
!     **  SAME i*G_D-MULTIPLICATION/FFT TECHNIQUE AS                   **
!     **  WAVES_KINETICDENSITY ABOVE (READ FIRST), APPLIED TWICE: ONCE  **
!     **  FORWARD (GRADIENT, GTOR) PER CARTESIAN DIRECTION, ONCE        **
!     **  BACKWARD (DIVERGENCE, RTOG) AFTER MULTIPLYING BY VTAU(R) IN   **
!     **  REAL SPACE, ACCUMULATING THE THREE DIRECTIONS' CONTRIBUTIONS  **
!     **  IN G-SPACE BEFORE THE FINAL -1/2 SCALING.                    **
!     **                                                                **
!     **  SAME SCOPE CUT AS WAVES_KINETICDENSITY: NDIM=1 ONLY, NO TINV. **
!     **  CALLER (WAVES_VPSI/WAVES_HPSI) MUST ADD HTAUPSI TO ITS OWN    **
!     **  OUTPUT -- THIS ROUTINE DOES NOT TOUCH THE ORDINARY LOCAL-     **
!     **  POTENTIAL OR KINETIC-ENERGY TERMS AT ALL.                    **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGL          ! MAX # PLANE WAVES
      INTEGER(4),INTENT(IN)  :: NRL          ! # R-SPACE POINTS
      INTEGER(4),INTENT(IN)  :: NDIM         ! DIMENSION OF THE WAVE FUNCTION
      INTEGER(4),INTENT(IN)  :: NBH          ! # BANDS (HALF-FILLED/BUNCHED)
      REAL(8)   ,INTENT(IN)  :: VTAU(NRL)    ! D(EXC)/D(TAU) IN REAL SPACE
      COMPLEX(8),INTENT(IN)  :: PSIOFG(NGL,NDIM,NBH)
      COMPLEX(8),INTENT(OUT) :: HTAUPSI(NGL,NDIM,NBH) ! -1/2 DIV(VTAU GRAD(PSI))
      COMPLEX(8),ALLOCATABLE :: PSIOFG_D(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSIOFR_D(:,:,:)
      COMPLEX(8),ALLOCATABLE :: WOFG_D(:,:,:)
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)
      LOGICAL(4)              :: TINV
      INTEGER(4)              :: IBH,IR,IDIR,IG
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
!     ******************************************************************
                               CALL TRACE$PUSH('WAVES_TAUHPSI')
!
!     ==================================================================
!     ==  SCOPE CHECKS -- SEE FILE HEADER / WAVES_KINETICDENSITY       ==
!     ==================================================================
      CALL PLANEWAVE$GETL4('TINV',TINV)
      IF(TINV) THEN
        CALL ERROR$MSG('TINV (SUPER WAVE FUNCTIONS) NOT YET SUPPORTED')
        CALL ERROR$STOP('WAVES_TAUHPSI')
      END IF
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('ONLY NDIM=1 (COLLINEAR) IS SUPPORTED SO FAR')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$STOP('WAVES_TAUHPSI')
      END IF
!
      ALLOCATE(GVEC(3,NGL))
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
!
      ALLOCATE(PSIOFG_D(NGL,NDIM,NBH))
      ALLOCATE(PSIOFR_D(NRL,NDIM,NBH))
      ALLOCATE(WOFG_D(NGL,NDIM,NBH))
      HTAUPSI(:,:,:)=(0.D0,0.D0)
      DO IDIR=1,3
!       == FORWARD: D(PSI)/D(X_IDIR) IN G-SPACE, THEN TO REAL SPACE ====
        DO IBH=1,NBH
          DO IG=1,NGL
            PSIOFG_D(IG,1,IBH)=CI*GVEC(IDIR,IG)*PSIOFG(IG,1,IBH)
          ENDDO
        ENDDO
        CALL PLANEWAVE$FFT('GTOR',NBH*NDIM,NGL,PSIOFG_D,NRL,PSIOFR_D)
!       == MULTIPLY BY VTAU(R) IN REAL SPACE ===========================
        DO IBH=1,NBH
          DO IR=1,NRL
            PSIOFR_D(IR,1,IBH)=VTAU(IR)*PSIOFR_D(IR,1,IBH)
          ENDDO
        ENDDO
!       == BACK TO G-SPACE, ACCUMULATE i*G_IDIR * (THIS DIRECTION) =====
        CALL PLANEWAVE$FFT('RTOG',NBH*NDIM,NGL,WOFG_D,NRL,PSIOFR_D)
        DO IBH=1,NBH
          DO IG=1,NGL
            HTAUPSI(IG,1,IBH)=HTAUPSI(IG,1,IBH)+CI*GVEC(IDIR,IG)*WOFG_D(IG,1,IBH)
          ENDDO
        ENDDO
      ENDDO
      HTAUPSI(:,:,:)=-0.5D0*HTAUPSI(:,:,:)
      DEALLOCATE(GVEC,PSIOFG_D,PSIOFR_D,WOFG_D)
                               CALL TRACE$POP
      RETURN
      END
