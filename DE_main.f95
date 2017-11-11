!------------------------------
!..   Fortran .f95 MAIN     ..
!------------------------------
! Uses a differential evolution (DE) algorthim to optimize the Griewank Function.
! This DE can be classified as DE/rand/1/bin with dither.
! License: https://github.com/ian-mmm/differential-evolution_f95/blob/master/LICENSE
!  -- required modules: DE_mod.f95
!  -- outputs:

PROGRAM main !% % % % % % % % %

! Load modules
USE DE_mod

IMPLICIT NONE

! Precision is set in DE_mod: dp = SELECTED_REAL_KIND(12, 100), sp = SELECTED_REAL_KIND(6, 35)

! INPUTS: *++*++*++*++*++*++*++*++*++*
    INTEGER, PARAMETER, DIMENSION(12)       :: seedA =(/ 11981, 31601, 90971, 854099, 295, 4593, &
                                                                369, 12437, 75, 6213, 532, 1989 /)
    INTEGER, PARAMETER                      :: nop = 4 !number of parameters
! *++*++*++*++*++*++*++*++*++*++*++*

! Variables:
    REAL(dp)                                :: start, finish
    INTEGER                                 :: NP, T
    REAL(dp)                                :: Fl, Fh, Cr, fval
    REAL(dp), DIMENSION( nop )              :: BetaUB, BetaLB, BetaHAT

!======================================================
! Seed random number generator
    CALL RANDOM_SEED(put = seedA )

! Start the stopwatch
    CALL CPU_TIME( start )

! MORE INPUTS: *++*++*++*++*++*++*++*++*++*
    ! Initial bounds for BETA:
        BetaLB(1:nop) = -1000.0D0
        BetaUB(1:nop) = 1000.0D0

    ! Optimization Control Variables:
        ! --Storn and Price (1997), J. of Global Opt., found these rules of thumbs: F \in [0.5, 1.0], Cr \in [0.8, 1.0], Np = 10*D

            ! Number of population vectors:
                NP = 10 * nop + 10

            ! Number of generations
                T = 10000

            ! Crossover variable: Cr \in [0,1]
                Cr = 0.85D0 ! manually set

            ! Scale factor determined by dither:
                Fl = SQRT(1.0D0 - 0.50D0* Cr )
                Fh = 0.950D0 ! L,H order should be fine as long as Cr > 0.2
! *++*++*++*++*++*++*++*++*++*++*++*

    CALL Diff_Evol( nop, NP, betaLB, betaUB, T, Fl, Fh, Cr, 1, BetaHAT, fval )
        ! -- SUBROUTINE Diff_Evol( nop, NP, LB, UB, T, F_lo, F_hi, Cr, PD, BH_best, F_best)

        PRINT*, "------"
    PRINT*, fval, BetaHat
        PRINT*, "------"

    ! Calculate time:
        CALL CPU_TIME( finish )
        PRINT '("Time = ",f15.2," seconds.")', finish - start

PRINT*, 'DONE ! ! ! ! ! ! ! '
END PROGRAM main !% % % % % % % % %
!======================================================
