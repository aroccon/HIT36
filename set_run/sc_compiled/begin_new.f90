subroutine begin_new

  use m_openmpi
  use m_parameters
  use m_fields
  use m_stats
  implicit none


  ! defining time
  TIME = zip

  ! deciding if we advance scalars or not
  if (TSCALAR.le.zip .and. n_scalars.gt.0) int_scalars = .true.

  ! defining the iteration number
  ITIME = 0
  file_ext = '000000'



  if (task.eq.'hydro') then
    call init_velocity
    if (n_scalars.gt.0) call init_scalars
    call io_write_4
  end if

  if (task_split) then
     call fields_to_stats
     if (task.eq.'stats') call stat_main
  end if

  return
end subroutine begin_new
