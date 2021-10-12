program hit36
use m_openmpi
use m_parameters
use m_fields
use m_work
use x_fftw
use m_stats
use m_timing
use m_force
use m_rand_knuth
use m_particles
use m_filter_xfftw
implicit none
integer :: n
character :: sym

! Setting CPU time to zero
call m_timing_init
! Splitting MPI (Hydro, Stats, parts)
call m_openmpi_init
! Read input file input.f90
call m_parameters_init
call m_fields_init
call m_work_init
! allocating and initializing FFTW arrays
call x_fftw_allocate(1)
call x_fftw_init

call m_stats_init
!init forcing (if activated)
call m_force_init

! allocating and initializing particles
if (task.eq.'parts') then
   call particles_init
   write(file_ext,"(i6.6)") itime
   if (itime.eq.0) call particles_restart_write_binary
end if

!Starting from the beginning or from the saved flowfield
if(itmin.eq.0) then
   call begin_new
else
   call begin_restart
endif

! checking divergence
if (task.eq.'hydro') call divergence

! indicators whether to use first-order in time
! for velocities and scalars
fov = .true.
fos = .true.

! need to dealias the fields at the beginning
if (task.eq.'hydro') call dealias_all

!Temporal loop
do itime=itmin+1,itmax

  !FLOW PART
  !define file_extension
  write(file_ext,"(i6.6)") itime
  if (task.eq.'hydro') then
    ! taking care of rescaling when running decaying turbulence
    ! if the time just was divisible by TRESCALE
    if (flow_type.eq.0 .and. floor((time-dt)/TRESCALE) .lt. floor(time/TRESCALE)) then
      ! ...and if we haven't rescaled NRESCALE times
      if (floor(time/TRESCALE) .le. NRESCALE .and. itime.ne.1) then
        call velocity_rescale
        ! after rescaling, the time-sceping needs to be first order
        fov = .true.; fos = .true.
        if (.not. task_split .and. mod(itime,iprint1).eq.0) call stat_main
      end if
    end if

    ! RHS for scalars
    call rhs_scalars

    ! Velocities in x-space, send them to parts for LPT
    if (task_split) call fields_to_parts

    ! advance scalars - either Euler or Adams-Bashforth
    if (int_scalars) then
      n = 3 + n_scalars
      if (fos) then
        rhs_old(:,:,:,4:n) = wrk(:,:,:,4:n)
        fields(:,:,:,4:n) = fields(:,:,:,4:n) + dt * rhs_old(:,:,:,4:n)
        fos = .false.
      else
        fields(:,:,:,4:n) = fields(:,:,:,4:n) + &
        dt * ( 1.5d0 * wrk(:,:,:,4:n)  - 0.5d0 * rhs_old(:,:,:,4:n) )
        rhs_old(:,:,:,4:n) = wrk(:,:,:,4:n)
      end if
    end if

    ! RHS for velocities
    call rhs_velocity

    ! Computing forcing (if forced flow)
    if (flow_type.eq.1) call force_velocity

    ! advance velocity - either Euler or Adams-Bashforth
    if (fov) then
      rhs_old(:,:,:,1:3) = wrk(:,:,:,1:3)
      fields(:,:,:,1:3) = fields(:,:,:,1:3) + dt * rhs_old(:,:,:,1:3)
      fov = .false.
    else
      fields(:,:,:,1:3) = fields(:,:,:,1:3) + &
      dt * ( 1.5d0 * wrk(:,:,:,1:3)  - 0.5d0 * rhs_old(:,:,:,1:3) )
      rhs_old(:,:,:,1:3) = wrk(:,:,:,1:3)
    end if

    ! solve for pressure and update velocities so they are incompressible
    call pressure

    ! advance time
    time = time + dt

    ! write the restart file if it's the time
    if (mod(itime,iprint2).eq.0) call restart_write_parallel

    ! CPU usage statistics
    if (mod(itime,iprint1).eq.0) then
      call m_timing_check
      if (mod(itime,iwrite4).eq.0) then
        sym = "*"
      else
        sym = " "
      end if
      if (myid_world.eq.0) write(*,9000) itime,time,courant,cpu_hrs,cpu_min,cpu_sec,sym
    end if

    if (mod(itime,iprint1).eq.0 .or. mod(itime,iwrite4).eq.0) then
      ! send the velocities to the "stats" part of the code for statistics
      if (task_split) call fields_to_stats
      ! checking if we need to stop the calculations due to simulation time
      if (time.gt.TMAX) call my_exit(1)
      ! checking if we need to start advancing scalars
      if (n_scalars.gt.0 .and. .not.int_scalars .and. time.gt.TSCALAR) then
        int_scalars = .true.
        call init_scalars
        !write(*,*) "Starting to move the scalars."
      end if
    end if
  end if

  !STATS PART
  if (task.eq.'stats' .or. .not.task_split) then
    if (mod(itime,iprint1).eq.0 .or. mod(itime,iwrite4).eq.0) then
      ! if this is a separate set of processors, then...
      if (task_split) then
      ! checking if we need to stop the calculations due to simulation time
        if (time.gt.tmax) call my_exit(1)
      end if

      ! these are executed regardless of the processor configuration
      if (task_split) call fields_to_stats
      if (mod(itime,iprint1).eq.0) call stat_main
      if (mod(itime,iwrite4).eq.0) call write_output

    end if
  end if


  !PARTICLE PARTS
  !  This is not enabled to work when task_split=.false.
  !  Currently the particles can be calculated only if we split the tasks due to
  !  requirements on the wrk array sizes in the particle interpolation routines.

  if (task.eq.'parts') then
    call fields_to_parts
    if (int_particles) then
      call particles_move
      if (mod(itime,iwrite4).eq.0) call particles_restart_write_binary
    end if

    if (mod(itime,iprint1).eq.0 .or. mod(itime,iwrite4).eq.0) then
      if (TIME.gt.TMAX) call my_exit(1)
    end if
  end if

  ! every iteration checking
  if (mod(ITIME,10).eq.0) then
    if (myid_world.eq.0) call m_timing_check
    count = 1
    call MPI_BCAST(cpu_min_total,count,MPI_INTEGER4,0,MPI_COMM_WORLD,mpi_err)
  end if

end do
! end of temporal loop

! In a case when we've gone to ITMAX, write the restart file
itime = itime-1
if (task.eq.'hydro') call restart_write_parallel
call my_exit(0)
call m_openmpi_exit
stop

9000 format('Iteration=',i6,3x,'t=',f8.4,4x,3x,'Courant= ',f6.4, &
          2x,'CPU:(',i4.4,':',i2.2,':',i2.2,')',x,a1,x,a3)

end program hit36
