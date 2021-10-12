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
real*8  :: TMAX, TRESCALE, TSCALAR, re, nu, dt
! now many times to rescale teh velocities
integer :: nrescale,flow_type
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
  kmax = floor(real(nx,8)/3.d0*sqrt(2.d0))
else
  if (myid_world .eq. 0) write(*,*) "Wrong dealias flag: ", dealias
  call my_exit(-1)
end if
if (myid_world .eq. 0) write(*,*) "Maximum resolved wavenumber - kmax = ", kmax

end subroutine m_parameters_init










subroutine read_input_file
implicit none
logical :: there
integer :: n
integer*4  :: passed, passed_all
character*80 :: str_tmp

!check if file input.f90 is there.
inquire(file='input.f90', exist=there)
if(.not.there) then
  if (myid_world .eq. 0) write(*,*) 'Cannot find input.f90 file'
  call my_exit(-1)
end if

! now the variable "passed" will show if the parameters make sense
passed = 1

open(10,file='input.f90',form='formatted')
read(10,*) nx
read(10,*)
read(10,*) itmin
read(10,*) itmax
read(10,*) iprint1
read(10,*) iprint2
read(10,*) iwrite4
read(10,*)
read(10,*) tmax
read(10,*) trescale, nrescale
read(10,*) tscalar
read(10,*)
read(10,*) flow_type
read(10,*)
read(10,*) re
read(10,*) dt
read(10,*)
read(10,*) isp_type
read(10,*) ir_exp
read(10,*) peak_wavenum
read(10,*)
read(10,*) force_type
read(10,*) kfmax
read(10,*) famp
read(10,*)
read(10,*) dealias
read(10,*)
read(10,*) det_rand
read(10,*) RN1
read(10,*) RN2
read(10,*) RN3
read(10,*)
read(10,*) nptot
read(10,*) particles_tracking_scheme
read(10,*) starttime_particles
read(10,*) particles_filter_size
read(10,*)
read(10,*) n_scalars
read(10,*)
if (n_scalars .gt. 0) then
  allocate(scalar_type(n_scalars), pe(n_scalars), sc(n_scalars), &
  ir_exp_sc(n_scalars), peak_wavenum_sc(n_scalars), &
  reac_sc(n_scalars))
  do n = 1,n_scalars
     read(10,*) scalar_type(n), sc(n), ir_exp_sc(n), &
          peak_wavenum_sc(n), reac_sc(n)
     write(*,'(9x,i4,1x,4(f8.3,1x))') scalar_type(n), sc(n), ir_exp_sc(n), &
          peak_wavenum_sc(n), reac_sc(n)
     pe(n) = nu/sc(n)
  end do
end if
close(10)


! Performing some checks on the parameters
! number of nodes
ny=nx
nz_all=nx
nz = nz_all/numprocs
if (nz*numprocs.ne.nz_all) then
  if (myid_world.eq.0)  write(*,*) 'Wrong combination of grid noes and number of MPI tasks'
  passed = 0
end if

! task_split
if (.not.task_split .and. nptot .gt. 0) then
   if (myid_world.eq.0) write(*,*) "Tasks are not split, making nptot=0"
   nptot = 0
end if

! particle_filter size
if (particles_filter_size .gt. 0.0d0 .and. particles_filter_size .lt. 3.d0*dx) then
   if (myid_world.eq.0) write(*,*) "particles_filter_size is too small (less than 3*dx)"
   call my_exit(-1)
end if

! defining the rest of the parameters
nxyz = nx * ny * nz
nxyz_all = nx * ny * nz_all
dx = twopi/dble(nx)
dy = twopi/dble(ny)
dz = twopi/dble(nz_all)
last_dump = itmin
nu = 1.0d0/re

count = 1
call MPI_REDUCE(passed,passed_all,count,MPI_INTEGER4,MPI_MIN,0,MPI_COMM_WORLD,mpi_err)
count = 1
call MPI_BCAST(passed_all,count,MPI_INTEGER4,0,MPI_COMM_WORLD,mpi_err)

if (passed.lt.one) then
  write(*,*) "Not passed the check, stopping"
  stop
end if

!print on screen some of the parameters
if (myid_world.eq.0)  write(*,*) 'Number of nodes NX x NY x NZ:', nx, ny, nz_all
if (myid_world.eq.0)  write(*,*) 'First timestep              :', itmin
if (myid_world.eq.0)  write(*,*) 'Last timestep               :', itmax
if (myid_world.eq.0)  write(*,*) 'Statistics computed every   :', iprint1
if (myid_world.eq.0)  write(*,*) 'Dump restart files          :', iprint2
if (myid_world.eq.0)  write(*,*) 'Dump output files           :', iwrite4
if (myid_world.eq.0)  write(*,*) 'Flow type (0=decay,1=forced):', flow_type
if (myid_world.eq.0)  write(*,*) 'Force type                  :', isp_type
if (myid_world.eq.0)  write(*,*) 'Reynolds number             :', re
if (myid_world.eq.0)  write(*,*) 'Time step                   :', dt
if (myid_world.eq.0)  write(*,*) 'Number of particles         :', nptot
if (myid_world.eq.0)  write(*,*) 'Number of passive sclars    :', n_scalars


return
end subroutine read_input_file






end module m_parameters
