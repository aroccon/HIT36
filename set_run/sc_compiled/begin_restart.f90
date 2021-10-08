subroutine begin_restart

  use m_parameters
  use m_particles
  use m_fields
  use m_timing
  implicit none

  real*8, allocatable :: zhopa(:,:,:,:)
  integer :: i


  write(out,"(70('='))") 
  write(out,*) '                    RESTART '
  write(out,"(70('='))") 
  call flush(out)

  ITIME = ITMIN
  call get_file_ext

  ! reading the restart file
  if (task.eq.'hydro') call restart_read_parallel

  ! broadcasting the current simulation time from the restart file
  call MPI_BCAST(time,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)

  ! redefining the timestep.
  ! the code uses Adams-Bashforth time-stepping scheme,
  ! which is 2nd order accurate in time.  After the restart, thefirst
  ! timestep is done using a simple Euler scheme (1st order in time).
  ! To help maintain any sensible accuracy, we need to start with
  ! the timestep which is smaller than the last.
  if (variable_dt) then
     dt = half * dt
     call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
     write(out,*) "Making the timestep smaller: ",dt
     call flush(out)
  end if


!!$!================================================================================
!!$!  Checking the parallel read speed (reading 100 times)
!!$!  Currently the speedup factor from using parallel read is about 2.5
!!$!================================================================================
!!$  if(task.eq.'hydro') then
!!$     call m_timing_check
!!$     write(out,*) 'Start! ',cpu_min,cpu_sec
!!$     do i = 1,200
!!$        call restart_read! _parallel
!!$     end do
!!$     call m_timing_check
!!$     write(out,*) 'Finish!',cpu_min,cpu_sec
!!$
!!$  end if
!!$  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!!$  stop 'checking the restart read speed'
!!$
!!$!================================================================================
!!$!  Checking the accuracy of the parallel read 
!!$!  (reading parallel first, then serial and comparing)
!!$!================================================================================
!  if (task.eq.'hydro') then
!
!     call restart_read_parallel
!
!     allocate(zhopa(nx+2,ny,nz,3+n_scalars), stat=ierr)
!     if (ierr.ne.0) stop 'cannot allocate zhopa'
!     zhopa = zip
!     zhopa = fields
!
!     call restart_read
!
!     print "(10e15.6)",maxval(abs(zhopa-fields))
!  end if
!  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!  stop 'checking the restart read'
!!$!================================================================================

  ! deciding whether we advance scalars or not
  if (n_scalars.gt.0  .and. time.gt.TSCALAR) then
     int_scalars = .true.
     write(out,"('Advancing ',i3,' scalars.')") n_scalars
  end if

  if (task.eq.'parts') call particles_init

  return
end subroutine begin_restart
