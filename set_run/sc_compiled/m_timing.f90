module m_timing

  integer*8 :: cpu0, cpu1, cpu2, dcpu
  integer*4 :: cpu_sec, cpu_min, cpu_hrs, cpu_min_total
  integer*4 :: job_runlimit

!================================================================================

contains
!================================================================================
  subroutine get_job_runlimit

    use m_openmpi
    use m_io
    implicit none

    logical :: there

    ! make the default job runlimit to be 6 months
    job_runlimit = 6 * 30 * 24 * 60

    ! make the default job runlimit to be 12 hours
    job_runlimit = 12 * 60

    write(out,*) 'job_runlimit (default):',job_runlimit

    if(myid_world.eq.0) then

       inquire(file='job_parameters.txt',exist=there)
       if (there) then
          open(98,file='job_parameters.txt')
          read(98,*,err=5) job_runlimit
5         close(98)
       end if
    end if
    call MPI_BCAST(job_runlimit,1,MPI_INTEGER4,0,MPI_COMM_WORLD,mpi_err)

    write(out,*) 'job_runlimit: ',job_runlimit
    call flush(out)

    return

  end subroutine get_job_runlimit

!================================================================================
!================================================================================

  subroutine m_timing_init

    call system_clock(cpu0,dcpu)
    return
  end subroutine m_timing_init

!================================================================================
!================================================================================
  subroutine m_timing_check


    use m_openmpi
    implicit none

    call system_clock(cpu1,dcpu)
    cpu_sec = (cpu1-cpu0)/dcpu
    cpu_min = cpu_sec/60;   cpu_min_total = cpu_min
    cpu_hrs = cpu_min/60
    cpu_min = mod(cpu_min,60)
    cpu_sec = mod(cpu_sec,60)

!!$    call MPI_BCAST(cpu_hrs,1,MPI_INTEGER4,master,MPI_COMM_TASK,mpi_err)
!!$    call MPI_BCAST(cpu_min,1,MPI_INTEGER4,master,MPI_COMM_TASK,mpi_err)

    return
  end subroutine m_timing_check

end module m_timing
