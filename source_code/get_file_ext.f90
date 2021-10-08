!================================================================================
!================================================================================
  subroutine get_file_ext
    use m_parameters, only : ITIME
    use m_io,         only : file_ext
    implicit none

    write(file_ext,"(i6.6)") itime
    return
  end subroutine get_file_ext

