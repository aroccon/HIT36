!================================================================================
! MY_EXT - de-initializes everything, called at the end the main program
!          or whenever we need to stop the program anywhere
!
!          reason - the passed variable that describes the particular reason why
!                   we want to stop the program.
!================================================================================
subroutine my_exit(reason)

  use m_openmpi
  use m_io
  use m_fields
  use m_work
!  use x_fftw
  use m_particles

  implicit none
  integer :: reason


  write(out,"('my_exit with reason ',i4)") reason

  select case (reason)

  case(0)

     write(out,*) '---------------------------------------------'
     write(out,*) '          NORMAL TERMIATION'
     write(out,*) '---------------------------------------------'

  case(1)

     write(out,*) '---------------------------------------------'
     write(out,*) '           TIME TERMINATION'
     write(out,*) '---------------------------------------------'
     call flush(out)

  case(2)

     write(out,*) '---------------------------------------------'
     write(out,*) '          RUN-TIME TERMINATION'
     write(out,*) '---------------------------------------------'

  case(3)

     write(out,*) '---------------------------------------------'
     write(out,*) '          USER TERMINATION'
     write(out,*) '---------------------------------------------'

  case default 

     write(out,*) '---------------------------------------------'
     write(out,*) '      TERMINATION FOR NO APPARENT REASON'
     write(out,*) '---------------------------------------------'

  end select


  if (reason.ge.0) then
     if (task.eq.'hydro') call restart_write_parallel
     if (task.eq.'parts') call particles_restart_write_binary
  end if


  write(out,*) "Done."
  call flush(out)
  close(out)
  stop


!  call m_fields_exit
!  call m_work_exit
!  call x_fftw_allocate(-1)
!  call m_io_exit
!  call m_openmpi_exit


  return
end subroutine my_exit

!!$!================================================================================
!!$! MY_INIT - initializes everything, called at the beginning of the main program
!!$!================================================================================
!!$subroutine my_init
!!$
!!$!  --- modules used
!!$  use m_openmpi
!!$  use m_io
!!$  use m_parameters
!!$  use m_fields
!!$  use m_work
!!$  use x_fftw
!!$  use m_filter_xfftw
!!$  implicit none
!!$  integer :: iargc
!!$
!!$!  --- initializing
!!$  call openmpi_init
!!$  call io_init
!!$  call parameters_init
!!$  call fields_init
!!$  call work_init(15)
!!$  call x_fftw_allocate(1)
!!$  call x_fftw_init
!!$
!!$!  --- getting the filter size
!!$  if(iargc().eq.0) then
!!$     call getarg(0,tmp_str)
!!$     write(out,*) 'Format: ',trim(tmp_str),' <filter_size>'
!!$     write(*,*)      'Format: ',trim(tmp_str),' <filter_size>'
!!$     call my_exit(-1)
!!$  end if
!!$
!!$  ftype = 2
!!$  call getarg(1,txt7)
!!$  read(txt7,*) filter_size
!!$  write(out,*) 'Filter size = ',filter_size
!!$  call flush(out)
!!$
!!$  call filter_xfftw_init
!!$
!!$return
!!$end subroutine my_init
