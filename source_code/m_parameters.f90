!================================================================================
! M_PARAMETERS - module for all parameters in the calculation:
!                such as array dimensions, reynolds numbers, switches/flags etc.
!
! Time-stamp: <2009-08-19 11:50:11 (chumakov)>
! Time-stamp: <2008-11-20 17:27:59 MST (vladimirova)>
!================================================================================

module m_parameters

  use m_openmpi
  use m_io
  implicit none

  ! --- problem related

  character*10 :: run_name

  ! --- input filep arameters
  integer :: nx,ny,nz, nz_all ! Dimensions of the problem
  integer :: nxyz, nxyz_all
  integer :: n_scalars        ! # of scalars
  real*8  :: time             ! time of simulation
  real*8  :: dx, dy, dz
  integer :: kmax


  integer :: ITIME, ITMIN, ITMAX, IPRINT1, IPRINT2, IWRITE4

  real*8  :: TMAX, TRESCALE, TSCALAR, RE, nu, dt

  ! now many times to rescale teh velocities
  integer :: NRESCALE

  integer :: flow_type

  logical :: variable_dt

  integer :: isp_type, ir_exp, force_type

  real*8  :: peak_wavenum

  real*8 :: famp

  integer :: kfmax ! Maximum wavenumber for forcing (integer)

  real*8 :: courant


  integer  :: dealias

  integer :: det_rand
  real*8  :: RN1, RN2, RN3

  ! particle-related

  ! indicator that says which particle tracking scheme to use:
  ! 0 = trilinear
  ! 1 = spectral
  ! 2 = tricubic
  ! trilinear by default
  integer :: particles_tracking_scheme = 0
  real*8  :: starttime_particles

  ! sometimes we want to advect particles by locally averaged field
  ! the following variables address that concern
  real*8  :: particles_filter_size

  ! number of particles assigned to the processor
  ! and the total number of particles
  integer(kind=MPI_INTEGER_KIND) :: np, np1, nptot

  ! If using Large Eddy Simulation (LES), the LES model ID is here
  integer :: les_model

  integer, allocatable :: scalar_type(:)
  real*8, allocatable  :: pe(:), sc(:), ir_exp_sc(:), peak_wavenum_sc(:), reac_sc(:)


  ! constants
  real*8  :: zip=0.0d0, half=0.5d0
  real*8  :: one=1.0d0,two=2.0d0,three=3.d0,four=4.d0,five=5.d0, six=6.d0


  integer :: last_dump

  ! --- supporting stuff
  logical      :: there
  logical      :: fos, fov
  integer      :: ierr
  real*8       :: PI, TWO_PI

  logical      :: int_scalars, int_particles


  ! --- number of LES variables in the arrays (initialized to zero)
  integer :: n_les = 0

  ! benchmarking tools
  logical :: benchmarking=.false.
  integer (kind=8) :: i81, i82, bm(12)

!================================================================================
contains
!================================================================================

  subroutine m_parameters_init

    implicit none

    call get_run_name

    ! constants
    PI     = four * atan(one)
    TWO_PI = two * PI

    ! switches
    int_scalars = .false.

    call read_input_file

    ! maximum resolved wavenumber
    if (dealias.eq.0) then
       kmax = nx/3
    elseif (dealias.eq.1) then
       kmax = floor(real(nx,8) / three * sqrt(two))
    else
       write(out,*) "*** M_PARAMETERS_INIT: wrong dealias flag: ",dealias
       call flush(out)
       call my_exit(-1)
    end if


    write(out,*) "kmax = ",kmax
    call flush(out)

  end subroutine m_parameters_init

!================================================================================

  subroutine get_run_name
    implicit none

    character*80 :: tmp_str
    integer      :: iargc

    ! reading the run_name from the command line
    if(iargc().eq.0) then
       call getarg(0,tmp_str)
       write(out,*) 'Format: ',trim(tmp_str),' <run name>'
       write(*,*)      'Format: ',trim(tmp_str),' <run name>'
       call flush(out)
       call MPI_FINALIZE(ierr)
       stop
    end if
    call getarg(1,run_name)
    if(len_trim(run_name).ne.10) then
       write(out,*) 'Run name: "',run_name,'"'
       write(out,*) '          "1234567890"'
       write(out,*) 'Length of run name is less than 10, sorry.'
       call MPI_FINALIZE(ierr)
       stop
    end if
    write(out,*) 'Run name: "',run_name,'"'
    call flush(out)

  end subroutine get_run_name


!================================================================================

  subroutine read_input_file

    implicit none

    logical :: there
    integer :: n
    integer*4  :: passed, passed_all
    character*80 :: str_tmp

    ! making sure the input file is there
    inquire(file='input.in', exist=there)
    if(.not.there) then
       write(out,*) '*** cannot find the input file'
       call flush(out)
       call my_exit(-1)
    end if

    ! now the variable "passed" will show if the parameters make sense
    passed = 1


    ! -------------------------------------------------
    ! reading parameters from the input file
    ! and checking them for consistency
    ! -------------------------------------------------
    open(in,file='input.in',form='formatted')
    read(in,*)
    read(in,*)
    read(in,*)

    read(in,*,ERR=9000) nx,ny,nz_all
    read(in,*)


    nz = nz_all/numprocs
    if (nz*numprocs.ne.nz_all) then
       write(out,*) '*** wrong nz_all:', nz_all, &
            '*** should be divisible by numprocs:',numprocs
       call flush(out)
       passed = 0
    end if
    write(out,'(70(''=''))')
    write(out,"('NX,NY,NZ_ALL', 3i4)") nx,ny,nz_all
    write(out,"('NX,NY,NZ    ', 3i4)") nx,ny,nz
    call flush(out)

    dx = 2.0d0 * PI / dble(nx)
    dy = 2.0d0 * PI / dble(ny)
    dz = 2.0d0 * PI / dble(nz_all)

    ! -------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) ITMIN
    write(out,*) 'ITMIN =   ',ITMIN
    last_dump = ITMIN

    read(in,*,ERR=9000,END=9000) ITMAX
    write(out,*) 'ITMAX =   ',ITMAX

    read(in,*,ERR=9000,END=9000) IPRINT1
    write(out,*) 'IPRINT1=   ',IPRINT1

    read(in,*,ERR=9000,END=9000) IPRINT2
    write(out,*) 'IPRINT2=   ',IPRINT2

    read(in,*,ERR=9000,END=9000) IWRITE4
    write(out,*) 'IWRITE4=  ',IWRITE4
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! ------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) TMAX
    write(out,*) 'TMAX     =',TMAX

    read(in,*,ERR=8000,END=9000) TRESCALE, NRESCALE
100 write(out,*) 'TRESCALE, NRESCALE =',TRESCALE, NRESCALE

    read(in,*,ERR=9000,END=9000) TSCALAR
    write(out,*) 'TSCALAR  =',TSCALAR
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

!      if(TSCALAR.le.TRESCALE) then
!        TSCALAR = TRESCALE
!        write(out,*) '*** RESET: TSCALAR = ',TSCALAR
!      end if

    ! ------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) flow_type
    write(out,*) 'flow_type    ', flow_type
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! ------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) RE
    write(out,*) 'RE    =   ',RE

    nu = 1.0d0/RE

    read(in,*,ERR=9000,END=9000) DT
    write(out,*) 'DT    =   ',DT
    if (dt.lt.0.0d0) then
       variable_dt = .false.
       dt = -dt
    else
       variable_dt = .true.
    end if

    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! ------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) isp_type
    write(out,*) 'isp_type=   ', isp_type

    read(in,*,ERR=9000,END=9000) ir_exp
    write(out,*) 'ir_exp   =   ', ir_exp

    read(in,*,ERR=9000,END=9000) peak_wavenum
    write(out,*) 'peak_wavenum =   ',peak_wavenum
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! ------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) force_type
    write(out,*) 'force_type', force_type

    read(in,*,ERR=9000,END=9000) kfmax
    write(out,*) 'kfmax =   ',kfmax

    read(in,*,ERR=9000,END=9000) FAMP
    write(out,*) 'FAMP  =   ',FAMP

    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

!!$    c------------------------------------------------------------
!!$
!!$    read(in,*,ERR=9000,END=9000) IRESET
!!$    write(out,*) 'IRESET=   ',IRESET
!!$
!!$    read(in,*,ERR=9000,END=9000) INEWSC
!!$    write(out,*) 'INEWSC=   ',INEWSC
!!$    read(in,*)
!!$
!!$    c------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) dealias
    write(out,*) 'dealias = ',dealias
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! -------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) det_rand
    write(out,*) 'det_rand =',det_rand

    read(in,*,ERR=9000,END=9000) RN1
    write(out,*) 'RN1      =',RN1

    read(in,*,ERR=9000,END=9000) RN2
    write(out,*) 'RN2      =',RN2

    read(in,*,ERR=9000,END=9000) RN3
    write(out,*) 'RN3      =',RN3
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! -------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) nptot
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    if (.not.task_split .and. nptot > 0) then
       write(out,*) "tasks are not split, making nptot=0"
       nptot = 0
    end if
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG

    write(out,*) 'nptot    =',nptot


    read(in,*,ERR=9000,END=9000) particles_tracking_scheme
    write(out,*) 'particles_tracking_scheme', particles_tracking_scheme

    select case (particles_tracking_scheme)
    case (0)
       write(out,*) '--- Trilinear tracking'
    case (1)
       write(out,*) '--- CINT (cubic interpolation on integer nodes)'
    case (2)
       write(out,*) '--- Spectral tracking (CAUTION: SLOW!)'
    case default
       write(out,*) 'don''t recognize particle tracking:', &
                  particles_tracking_scheme
       write(out,*) 'reset to zero'
       particles_tracking_scheme = 0
    end select
    call flush(out)


! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    if (particles_tracking_scheme .gt. 1) stop 'Cannot do this particle tracking'
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG



    read(in,*,ERR=9000,END=9000) starttime_particles
    write(out,*) 'starttime_particles: ',starttime_particles

    read(in,*,ERR=9000,END=9000) particles_filter_size
    write(out,*) 'particles_filter_size:',particles_filter_size

    if (particles_filter_size .gt. zip .and. particles_filter_size .lt. three*dx) then
       write(out,*) "particles_filter_size is too small (less than 3*dx)"
       write(out,*) particles_filter_size, three*dx
       call flush(out)
       call my_exit(-1)
    end if

    read(in,*)


    ! -------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) les_model
    write(out,*) 'les_model    =',les_model
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! making sure that if the LES mode is on, the dealiasing is 3/2-rule
    if (les_model .gt. 0 .and. dealias .ne. 0) then
       dealias = 0
       write(out,*) "*** LES mode, changing dealias to 0."
       call flush(out)
    end if

    ! -------------------------------------------------------------

    read(in,*,ERR=9000,END=9000) n_scalars
    write(out,*) '# of scalars:', n_scalars
    read(in,*)
    write(out,"(70('-'))")
    call flush(out)

    ! ------------------------------------------------------------

    ! if there are scalars, then read them one by one
    if (n_scalars>0) then
       read(in,'(A)',ERR=9000,END=9000) str_tmp
       write(out,*) str_tmp
       call flush(out)

       ! reading parameters of each scalar
       allocate(scalar_type(n_scalars), pe(n_scalars), sc(n_scalars), &
            ir_exp_sc(n_scalars), peak_wavenum_sc(n_scalars), &
            reac_sc(n_scalars), stat=ierr)
       if (ierr.ne.0) passed = 0

       do n = 1,n_scalars
          read(in,*,ERR=9000,END=9000) scalar_type(n), sc(n), ir_exp_sc(n), &
               peak_wavenum_sc(n), reac_sc(n)
          write(out,'(9x,i4,1x,4(f8.3,1x))') scalar_type(n), sc(n), ir_exp_sc(n), &
               peak_wavenum_sc(n), reac_sc(n)
               call flush(out)

          PE(n) = nu/SC(n)       ! INVERSE Peclet number

       end do
    end if

    ! -------------------------------------------------------------

    ! closing the input file
    close(in)
    write(out,'(70(''=''))')
    call flush(out)

    ! defining the rest of the parameters

    nxyz = nx * ny * nz
    nxyz_all = nx * ny * nz_all

    ! ------------------------------------------------------------


!--------------------------------------------------------------------------------
!  Checking if the task splitting conflicts with particle advection.  Currently
!  we canot have split=never and have particles.  This is to be resolved later,
!  now my head is spinning already.
!--------------------------------------------------------------------------------
    if (.not.task_split .and. nptot.gt.0) then
       write(out,*) "*** READ_INPUT_FILE: Cannot have .not.task_split and nptot > 0.  Stopping"
       call flush(out)
       passed = 0
    end if
!--------------------------------------------------------------------------------


    count = 1
    call MPI_REDUCE(passed,passed_all,count,MPI_INTEGER4,MPI_MIN,0,MPI_COMM_WORLD,mpi_err)
    count = 1
    call MPI_BCAST(passed_all,count,MPI_INTEGER4,0,MPI_COMM_WORLD,mpi_err)

    if (passed.lt.one) then
       write(out,*) "not passed the check, stopping"
       call flush(out)
       stop
    end if

    return

!--------------------------------------------------------------------------------
!  ERROR PROCESSING
!--------------------------------------------------------------------------------

8000 continue
    NRESCALE = 0
    if (TRESCALE.gt.zip) NRESCALE = 1
    write(out,*) "*** NRESCALE IS AUTOMATICALLY ASSIGNED to be ONE"
    call flush(out)
    goto 100

9000 continue
    write(out,*)'An error was encountered while reading input file'
    call flush(out)
    stop
  end subroutine read_input_file

!================================================================================
end module m_parameters
