module m_openmpi
implicit none
include 'mpif.h'
! Uncomment this for the systems that do not have OpenMPI
! In OpenMPI, the parameter MPI_INTEGER_KIND is defined in 'mpif.h'
! With other MPI implementations, this parameter has to be defined manually.
!  integer MPI_INTEGER_KIND
!  parameter (MPI_INTEGER_KIND = 4)
! --- MPI variables
logical :: iammaster
integer(kind=MPI_INTEGER_KIND) :: myid_world, numprocs_world
integer(kind=MPI_INTEGER_KIND) :: numprocs_hydro, numprocs_stats, numprocs_parts
integer(kind=MPI_INTEGER_KIND) :: myid, numprocs, master, mpi_err, mpi_info
integer(kind=MPI_INTEGER_KIND) :: id_to, id_from, tag, count
integer(kind=MPI_INTEGER_KIND) :: id_root_hydro, id_root_stats, id_root_parts
! communicator for separate tasks
integer(kind=MPI_INTEGER_KIND) :: MPI_COMM_TASK
! exclusive communicator for root processes of tasks
integer(kind=MPI_INTEGER_KIND) :: MPI_COMM_ROOTS
integer (kind=MPI_INTEGER_KIND) :: sendtag, recvtag
integer (kind=MPI_INTEGER_KIND) :: request, request1, request2, request3, mpi_request
integer (kind=MPI_INTEGER_KIND) :: id_l, id_r
integer (kind=mpi_INTEGER_KIND) :: mpi_status(MPI_STATUS_SIZE)
integer(kind=MPI_INTEGER_KIND) :: color, key
character*5 :: task, split="never"
!character*10 :: run_name_local
logical :: task_split

contains





subroutine m_openmpi_init
implicit none
integer (kind=mpi_INTEGER_KIND) :: n
integer*4 :: np_local
integer :: i

! initializing MPI
call mpi_init(mpi_err)
call mpi_comm_size(MPI_COMM_WORLD,numprocs_world,mpi_err)
call mpi_comm_rank(MPI_COMM_WORLD,myid_world,mpi_err)

! split="split" means that hydro, statistics and particles are assigned three
! separate process groups (they differ by the char*5 parameter "task").
! split="never" (default if the parameter is missing) means that all
! processes do all tasks (not working for particles)
if (split == "never") task_split = .false.
if (split == "split") task_split = .true.

! check task_split value
! if task_split = .false., we perform a ficticious split with uniform color
if (.not. task_split) then
  if (myid_world.eq.0) print*,'Parallelization strategy: No split'
  task = 'hydro'
  color = 0
  myid = myid_world
!  goto 1000
end if

! if task_split = .true., we perform a ficticious splitting
if (task_split) then

! check if we need to split hydro-stats-parts or only hydro and parts.
! We read the number of particles
  if (myid_world.eq.0) then
    open(99,file='input.f90')
    do i = 1,35
    read(99,*)
    end do
    read(99,*) np_local
    close(99)
  end if
  count = 1
  call MPI_BCAST(np_local,count,MPI_INTEGER4,0,MPI_COMM_WORLD,mpi_err)

  ! No particles, the split between the hydro and the stats part is 2/3 and 1/3.
  ! Total number of processor sneeds to be 3*2^n
  if (np_local.eq.0) then
    if (myid_world.eq.0) print*,'Parallelization strategy: Split - 2/3 Hydro and 1/3 Stats'
    if (int(numprocs_world/3)*3 .eq. numprocs_world) then
      numprocs_hydro = numprocs_world * 2/3
      numprocs_stats = numprocs_world - numprocs_hydro
      numprocs_parts = 0
    else if (2**floor(log(real(numprocs_world))/log(2.d0)) .eq. numprocs_world) then
      print*, 'numprocs_world is 2^n, allocating half for hydro: ',numprocs_world
      numprocs_hydro = numprocs_world / 2
      numprocs_stats = numprocs_world - numprocs_hydro
      numprocs_parts = 0
    else
      numprocs_hydro = 2**floor(log(real(numprocs_world))/log(2.d0))
      numprocs_stats = numprocs_world - numprocs_hydro
      numprocs_parts = 0
    end if

    id_root_hydro = 0
    id_root_stats = numprocs_hydro
    id_root_parts = 0
  end if

  ! There are particles, the split between hydro, stats and parts is 1/2, 1/4, 1/4.
  if (np_local.gt.0) then
    if (myid_world.eq.0) print*,'Parallelization strategy: Split - 1/2 Hydro and 1/4 Stats and 1/4 Particles'
    numprocs_hydro = numprocs_world / 2
    numprocs_stats = numprocs_world / 4
    numprocs_parts = numprocs_world / 4

    id_root_hydro = 0
    id_root_stats = numprocs_hydro
    id_root_parts = numprocs_hydro + numprocs_stats
  end if

   ! splitting the communicator into several parts
  if (myid_world.lt.numprocs_hydro) then
    task = 'hydro'
    color = 0
    myid = myid_world
  elseif (myid_world.ge.numprocs_hydro .and. myid_world .lt. numprocs_hydro+numprocs_stats) then
    task = 'stats'
    color = 1
    myid = myid_world - numprocs_hydro
  else
    task = 'parts'
    color = 2
    myid = myid_world - numprocs_hydro - numprocs_stats
  end if
end if

! Fictiuos (task_split .false.) or real splitting (task_split .true.)
!1000 continue
! Call MPI to do the splitting
call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,myid,MPI_COMM_TASK,mpi_err)
call MPI_COMM_SIZE(MPI_COMM_TASK,numprocs,mpi_err)
call MPI_COMM_RANK(MPI_COMM_TASK,myid,mpi_err)

! each task will have its master process
master = 0
iammaster = .false.
if (myid.eq.master) iammaster=.true.

return
end subroutine m_openmpi_init









subroutine m_openmpi_exit

call MPI_COMM_FREE(MPI_COMM_TASK,mpi_err)
call MPI_FINALIZE(mpi_err)

return
end subroutine m_openmpi_exit









end module m_openmpi
