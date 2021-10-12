module m_work
use m_parameters
implicit none
real*4,allocatable :: tmp4(:,:,:)
real*8, allocatable :: wrk(:,:,:,:)
real*8, allocatable :: wrkp(:,:)
real*8, allocatable :: rhs_old(:,:,:,:)


contains

subroutine m_work_init
use m_parameters
implicit none

if (task.eq.'hydro')  then
  call m_work_allocate(max(6,3+n_scalars+2))
elseif (task.eq.'stats')  then
  call m_work_allocate(6)
elseif (task.eq.'parts') then
  ! required recources are different for the "parts" part of the code
  ! we need several layers of velocities for velocity interpolation
  select case (particles_tracking_scheme)
  case (0)
    allocate(wrk(1:nx+2,1:ny,1:3,0:0),stat=ierr)
    !write(*,*) "Allocated wrk(1:nx+2,1:ny,1:3,0:0)", ierr
    wrk = 0.0d0
  case (1)
    allocate(wrk(1:nx+2,1:ny,1:3,1:3),stat=ierr)
    !write(*,*) "Allocated wrk(1:nx+2,1:ny,1:3,1:3)", ierr
    wrk = 0.0d0
  case default
    stop 'Wrong particles_tracking_scheme'
  end select
       allocate(tmp4(nx,ny,nz), stat=ierr)
       !write(*,*) "Allocated tmp4."
  else
    !write(*,*) 'TASK variable is set in such a way that I dont know how to allocate work arrays'
    !write(*,*) 'task = ',task
    call my_exit(-1)
end if

return
end subroutine m_work_init






!  allocating and defining the prescribed number of arrays
subroutine m_work_allocate(number)
implicit none
integer :: number,i

! array that is needed for output (nx,ny,nz)
allocate(tmp4(nx,ny,nz))

! main working array, needed for FFT etc, so (nx+2,ny,nz)
allocate(wrk(nx+2,ny,nz,0:number))

! array for the spare RHS for Adams-Bashforth time-stepping scheme methods
if (task.eq.'hydro') then
  allocate(rhs_old(nx+2,ny,nz,3+n_scalars))
end if

tmp4 = 0.0
wrk = zip
if (allocated(rhs_old)) rhs_old = 0.d0

return
end subroutine m_work_allocate




subroutine m_work_exit
implicit none

if (allocated(tmp4)) deallocate(tmp4)
if (allocated(wrk)) deallocate(wrk)
if (allocated(wrkp)) deallocate(wrkp)
write(*,*) 'deallocated wrk arrays'

return
end subroutine m_work_exit


end module m_work
