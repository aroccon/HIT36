subroutine my_exit(reason)
use m_openmpi
use m_fields
use m_work
use m_particles
implicit none
integer :: reason


!write(*,"('my_exit with reason ',i4)") reason


if (reason.ge.0) then
  if (task.eq.'hydro') call restart_write_parallel
  if (task.eq.'parts') call particles_restart_write_binary
end if


!  write(out,*) "Done."
!  call flush(out)
!  close(out)
stop

return
end subroutine my_exit
