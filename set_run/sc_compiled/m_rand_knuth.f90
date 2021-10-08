MODULE m_rand_knuth


  ! Random number generator


  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2000-09-10  Time: 16:37:48
  ! Latest revision - 16 January 2003

  ! FORTRAN 77 version of "ran_array"
  ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
  ! including the MODIFICATIONS made in the 9th printing (2002)
  ! ********* see the book for explanations and caveats! *********
  ! Author: Steve Kifowit
  ! http://ourworld.compuserve.com/homepages/steve_kifowit
  ! with modifications by Alan Miller to rnarry and rnstrt based upon
  ! Knuth's code.

  ! For Donald Knuth's Fortran 77 versions, go to:
  ! http://www-cs-faculty.stanford.edu/~knuth/programs
  ! Look for frng.f and frngdb.f

  IMPLICIT NONE
  INTEGER, PARAMETER  :: kk=100, ll=37, mm=2**30, tt=70, kkk=kk+kk-1
  INTEGER, SAVE       :: ranx(kk)

CONTAINS


  SUBROUTINE rand_knuth(u, n)

    ! Generate an array of n real values between 0 and 1.

    REAL   , INTENT(OUT)    :: u(:)
    INTEGER, INTENT(IN)  :: n

    ! Local array
    INTEGER  :: aa(n)

    CALL rnarry(aa, n)
    u(1:n) = SCALE( REAL(aa), -30)

    RETURN
  END SUBROUTINE rand_knuth



  SUBROUTINE rnarry(aa, n)

    ! Generate an array of n integers between 0 and 2^30-1.

    INTEGER, INTENT(OUT)  :: aa(:)
    INTEGER, INTENT(IN)   :: n

    ! Local variables
    INTEGER  :: j

    aa(1:kk) = ranx(1:kk)
    DO  j = kk + 1, n
       aa(j) = aa(j-kk) - aa(j-ll)
       IF (aa(j) < 0) aa(j) = aa(j) + mm
    END DO
    DO  j=1,ll
       ranx(j) = aa(n+j-kk) - aa(n+j-ll)
       IF (ranx(j) < 0) ranx(j) = ranx(j) + mm
    END DO
    DO  j=ll+1,kk
       ranx(j) = aa(n+j-kk) - ranx(j-ll)
       IF (ranx(j) < 0) ranx(j) = ranx(j) + mm
    END DO

    RETURN
  END SUBROUTINE rnarry



  SUBROUTINE rand_knuth_start(seed)

    ! Initialize integer array ranx using the input seed.

    INTEGER*8, INTENT(IN)  :: seed

    ! Local variables
    INTEGER  :: x(kkk), j, ss, sseed, t

    IF (seed < 0) THEN
       sseed = mm - 1 - MOD(-1-seed, mm)
    ELSE
       sseed = MOD(seed, mm)
    END IF
    ss = sseed - MOD(sseed,2) + 2
    DO  j=1, kk
       x(j) = ss
       ss = ISHFT(ss, 1)
       IF (ss >= mm) ss = ss - mm + 2
    END DO
    x(kk+1:kkk) = 0
    x(2) = x(2)+1
    ss = sseed
    t = tt - 1
10  DO  j=kk, 2, -1
       x(j+j-1) = x(j)
    END DO
    DO  j = kkk, kk + 1, -1
       x(j-(kk-ll)) = x(j-(kk-ll)) - x(j)
       IF (x(j-(kk-ll)) < 0) x(j-(kk-ll)) = x(j-(kk-ll)) + mm
       x(j-kk) = x(j-kk) - x(j)
       IF (x(j-kk) < 0) x(j-kk) = x(j-kk) + mm
    END DO
    IF (MOD(ss,2) == 1) THEN
       DO  j=kk, 1, -1
          x(j+1) = x(j)
       END DO
       x(1) = x(kk+1)
       x(ll+1) = x(ll+1) - x(kk+1)
       IF (x(ll+1) < 0) x(ll+1) = x(ll+1) + mm
    END IF
    IF (ss /= 0) THEN
       ss = ISHFT(ss, -1)
    ELSE
       t = t - 1
    END IF
    IF (t > 0) GO TO 10

    DO  j=1, ll
       ranx(j+kk-ll) = x(j)
    END DO
    DO  j=ll+1,kk
       ranx(j-ll) = x(j)
    END DO

    DO  j = 1, 10
       CALL rnarry(x,kkk)
    END DO

    RETURN
  END SUBROUTINE rand_knuth_start


  ! Initialization subroutine

  subroutine rand_knuth_init(seed1)

    integer*8 :: seed1
    
    call rand_knuth_start(seed1)


    return
  end subroutine rand_knuth_init


END MODULE m_rand_knuth

