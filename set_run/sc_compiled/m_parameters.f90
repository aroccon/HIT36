module m_parameters
use m_openmpi
implicit none
character*10 :: run_name
! --- input filep parameters
integer :: nx,ny,nz, nz_all ! Dimensions of the problem
integer :: nxyz, nxyz_all
integer :: n_scalars        ! # of scalars
real*8  :: time             ! time of simulation
real*8  :: dx, dy, dz
integer :: kmax
integer :: ITIME, ITMIN, ITMAX, IPRINT1, IPRINT2, IWRITE4
real*8  :: TMAX, TRESCALE, TSCALAR, RE, nu, dt
! now many times to rescale teh velocities
integer :: NRESCALE,variable_dt,flow_type
integer :: isp_type, ir_exp, force_type
real*8  :: peak_wavenum,famp,courant
integer :: kfmax ! Maximum wavenumber for forcing (integer)
integer  :: dealias
integer :: det_rand
real*8  :: RN1, RN2, RN3
character*6  :: file_ext
character*80 :: fname
! particle-related
! indicator that says which particle tracking scheme to use:
! 0 = trilinear
! 1 = spectral
! 2 = tricubic
! trilinear by default
integer :: particles_tracking_scheme = 0
real*8  :: starttime_particles
! sometimes we want to advect particles by locally averaged field
real*8  :: particles_filter_size
! number of particles assigned to the processor and total number of particles
integer(kind=MPI_INTEGER_KIND) :: np, np1, nptot
! If using Large Eddy Simulation (LES), the LES model ID is here
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
real*8       :: pi=3.14159265358979, twopi=6.28318530717958
logical      :: int_scalars, int_particles
! --- number of LES variables in the arrays (initialized to zero)
integer :: n_les = 0
! benchmarking tools
logical :: benchmarking=.false.
integer (kind=8) :: i81, i82, bm(12)




contains



subroutine m_parameters_init

implicit none

int_scalars = .false.

call read_input_file

! maximum resolved wavenumber
if (dealias.eq.0) then
  kmax = nx/3
elseif (dealias.eq.1) then
  kmax = floor(real(nx,8) / three * sqrt(two))
else
  write(*,*) "*** M_PARAMETERS_INIT: wrong dealias flag: ",dealias
  call my_exit(-1)
end if

write(*,*) "kmax = ",kmax

end subroutine m_parameters_init










subroutine read_input_file
implicit none
logical :: there
integer :: n
integer*4  :: passed, passed_all
character*80 :: str_tmp

inquire(file='input.f90', exist=there)
if(.not.there) then
  write(*,*) '*** cannot find the input file'
  call my_exit(-1)
end if

! now the variable "passed" will show if the parameters make sense
passed = 1

open(10,file='input.f90',form='formatted')
read(10,*)
read(10,*)
read(10,*)
read(10,*,ERR=9000) nx,ny,nz_all
read(10,*)
nz = nz_all/numprocs
if (nz*numprocs.ne.nz_all) then
  write(*,*) '*** wrong nz_all:', nz_all, &
  '*** should be divisible by numprocs:',numprocs
  passed = 0
end if
write(11,'(70(''=''))')
write(11,"('NX,NY,NZ_ALL', 3i4)") nx,ny,nz_all
write(11,"('NX,NY,NZ    ', 3i4)") nx,ny,nz

dx = twopi/dble(nx)
dy = twopi/dble(ny)
dz = twopi/dble(nz_all)

read(10,*,ERR=9000,END=9000) ITMIN
last_dump = ITMIN
read(10,*,ERR=9000,END=9000) ITMAX
read(10,*,ERR=9000,END=9000) IPRINT1
read(10,*,ERR=9000,END=9000) IPRINT2
read(10,*,ERR=9000,END=9000) IWRITE4
read(10,*)
read(10,*,ERR=9000,END=9000) TMAX


write(*,*) 'ITMAX =    ',ITMAX
write(*,*) 'IPRINT1=   ',IPRINT1
write(*,*) 'IPRINT2=   ',IPRINT2
write(*,*) 'IWRITE4=   ',IWRITE4
write(*,*) 'TMAX     =',TMAX

    read(10,*,ERR=8000,END=9000) TRESCALE, NRESCALE
100 write(*,*) 'TRESCALE, NRESCALE =',TRESCALE, NRESCALE

    read(10,*,ERR=9000,END=9000) TSCALAR
    write(*,*) 'TSCALAR  =',TSCALAR
    read(10,*)
    write(*,"(70('-'))")
    read(10,*,ERR=9000,END=9000) flow_type
    write(*,*) 'flow_type    ', flow_type
    read(10,*)
    write(*,"(70('-'))")

    ! ------------------------------------------------------------

    read(10,*,ERR=9000,END=9000) RE
    write(*,*) 'RE    =   ',RE
    nu = 1.0d0/RE
    read(10,*,ERR=9000,END=9000) DT
    write(*,*) 'DT    =   ',DT
    if (dt.lt.0.0d0) then
       variable_dt = .false.
       dt = -dt
    else
       variable_dt = .true.
    end if

    read(10,*)
    write(*,"(70('-'))")

    ! ------------------------------------------------------------

    read(10,*,ERR=9000,END=9000) isp_type
    write(11,*) 'isp_type=   ', isp_type

    read(10,*,ERR=9000,END=9000) ir_exp
    write(*,*) 'ir_exp   =   ', ir_exp

    read(10,*,ERR=9000,END=9000) peak_wavenum
    write(*,*) 'peak_wavenum =   ',peak_wavenum
    read(10,*)
    write(*,"(70('-'))")

    ! ------------------------------------------------------------

    read(10,*,ERR=9000,END=9000) force_type
    write(*,*) 'force_type', force_type

    read(10,*,ERR=9000,END=9000) kfmax
    write(*,*) 'kfmax =   ',kfmax

    read(10,*,ERR=9000,END=9000) FAMP
    write(*,*) 'FAMP  =   ',FAMP

    read(10,*)
    write(*,"(70('-'))")

!!$    c------------------------------------------------------------
!!$
!!$    read(10,*,ERR=9000,END=9000) IRESET
!!$    write(out,*) 'IRESET=   ',IRESET
!!$
!!$    read(10,*,ERR=9000,END=9000) INEWSC
!!$    write(out,*) 'INEWSC=   ',INEWSC
!!$    read(10,*)
!!$
!!$    c------------------------------------------------------------

    read(10,*,ERR=9000,END=9000) dealias
    write(*,*) 'dealias = ',dealias
    read(10,*)
    write(*,"(70('-'))")

    ! -------------------------------------------------------------

    read(10,*,ERR=9000,END=9000) det_rand
    write(*,*) 'det_rand =',det_rand

    read(10,*,ERR=9000,END=9000) RN1
    write(*,*) 'RN1      =',RN1

    read(10,*,ERR=9000,END=9000) RN2
    write(*,*) 'RN2      =',RN2

    read(10,*,ERR=9000,END=9000) RN3
    write(*,*) 'RN3      =',RN3
    read(10,*)
    write(*,"(70('-'))")

    ! -------------------------------------------------------------

    read(10,*,ERR=9000,END=9000) nptot
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    if (.not.task_split .and. nptot > 0) then
       write(*,*) "tasks are not split, making nptot=0"
       nptot = 0
    end if
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG

    write(*,*) 'nptot    =',nptot


    read(10,*,ERR=9000,END=9000) particles_tracking_scheme
    write(*,*) 'particles_tracking_scheme', particles_tracking_scheme

    select case (particles_tracking_scheme)
    case (0)
       write(*,*) '--- Trilinear tracking'
    case (1)
       write(*,*) '--- CINT (cubic interpolation on integer nodes)'
    case (2)
       write(*,*) '--- Spectral tracking (CAUTION: SLOW!)'
    case default
       write(*,*) 'don''t recognize particle tracking:', &
                  particles_tracking_scheme
       write(*,*) 'reset to zero'
       particles_tracking_scheme = 0
    end select


! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    if (particles_tracking_scheme .gt. 1) stop 'Cannot do this particle tracking'
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG



    read(10,*,ERR=9000,END=9000) starttime_particles
    write(*,*) 'starttime_particles: ',starttime_particles

    read(10,*,ERR=9000,END=9000) particles_filter_size
    write(*,*) 'particles_filter_size:',particles_filter_size

    if (particles_filter_size .gt. zip .and. particles_filter_size .lt. three*dx) then
       write(*,*) "particles_filter_size is too small (less than 3*dx)"
       write(*,*) particles_filter_size, three*dx
       call my_exit(-1)
    end if

    read(10,*)


    ! -------------------------------------------------------------

    read(10,*,ERR=9000,END=9000)
    read(10,*)
    write(*,"(70('-'))")

    read(10,*,ERR=9000,END=9000) n_scalars
    write(*,*) '# of scalars:', n_scalars
    read(10,*)
    write(*,"(70('-'))")

    ! ------------------------------------------------------------

    ! if there are scalars, then read them one by one
    if (n_scalars>0) then
       read(10,'(A)',ERR=9000,END=9000) str_tmp
       write(*,*) str_tmp

       ! reading parameters of each scalar
       allocate(scalar_type(n_scalars), pe(n_scalars), sc(n_scalars), &
            ir_exp_sc(n_scalars), peak_wavenum_sc(n_scalars), &
            reac_sc(n_scalars), stat=ierr)
       if (ierr.ne.0) passed = 0

       do n = 1,n_scalars
          read(10,*,ERR=9000,END=9000) scalar_type(n), sc(n), ir_exp_sc(n), &
               peak_wavenum_sc(n), reac_sc(n)
          write(*,'(9x,i4,1x,4(f8.3,1x))') scalar_type(n), sc(n), ir_exp_sc(n), &
               peak_wavenum_sc(n), reac_sc(n)

          PE(n) = nu/SC(n)       ! INVERSE Peclet number

       end do
    end if

    ! -------------------------------------------------------------

    ! closing the input file
    close(10)
    write(*,'(70(''=''))')

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
       write(*,*) "*** READ_INPUT_FILE: Cannot have .not.task_split and nptot > 0.  Stopping"
       passed = 0
    end if
!--------------------------------------------------------------------------------


    count = 1
    call MPI_REDUCE(passed,passed_all,count,MPI_INTEGER4,MPI_MIN,0,MPI_COMM_WORLD,mpi_err)
    count = 1
    call MPI_BCAST(passed_all,count,MPI_INTEGER4,0,MPI_COMM_WORLD,mpi_err)

    if (passed.lt.one) then
       write(*,*) "not passed the check, stopping"
       stop
    end if

    return

!--------------------------------------------------------------------------------
!  ERROR PROCESSING
!--------------------------------------------------------------------------------

8000 continue
    NRESCALE = 0
    if (TRESCALE.gt.zip) NRESCALE = 1
    write(*,*) "*** NRESCALE IS AUTOMATICALLY ASSIGNED to be ONE"
    goto 100

9000 continue
    write(*,*)'An error was encountered while reading input file'
    stop
end subroutine read_input_file

end module m_parameters
