subroutine write_tmp4

  use m_parameters
  use m_io
  use m_work
  implicit none

  integer :: my_out=13, i,j,k
  integer*4 :: sizes(3)

!======================================================================

  ! --- defining the size of whole array
  sizes(1)=nx;   sizes(2)=ny;   sizes(3)=nz*numprocs;
  count = nx*ny*nz

  if (myid.ne.master) then
     id_to = master
     tag = myid
     call MPI_SEND(tmp4,count,MPI_REAL4,master,tag,MPI_COMM_TASK,mpi_err)
  else
!!     open(my_out,file=fname,form='binary')
     open(my_out,file=fname,form='unformatted', access='stream')
     write(my_out) sizes(1:3)
     write(my_out) (((tmp4(i,j,k),i=1,nx),j=1,ny),k=1,nz)

     do id_from=1,numprocs-1
        tag = id_from
        call MPI_RECV(tmp4,count,MPI_REAL4,id_from,tag,MPI_COMM_TASK,mpi_status,mpi_err)
        write(my_out) (((tmp4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
     end do
     close(my_out)
  end if

!======================================================================  

  return
end subroutine write_tmp4

!======================================================================
subroutine write_tmp4_all

  use m_parameters
  use m_io
  use m_work
  implicit none

  integer :: my_out=13, i,j,k
  integer(kind=MPI_INTEGER_KIND) :: sizes(3), fh
  integer(kind=MPI_OFFSET_KIND)  :: offset

!======================================================================

  ! --- defining the size of whole array
  sizes(1)=nx;   sizes(2)=ny;   sizes(3)=nz*numprocs;

  ! --- writing into the file with appropriate offset
  call MPI_INFO_CREATE(mpi_info,mpi_err)
  if(ierr.ne.0) stop '*** WRITE_TMP4_ALL: cannot create mpi_info'

  call MPI_FILE_OPEN(MPI_COMM_TASK,fname,MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,fh,mpi_err)
  if (myid.eq.0) call MPI_FILE_WRITE_AT(fh,0,sizes,3,MPI_INTEGER4,mpi_status,mpi_err)
  offset = 12 + myid*nx*ny*nz * 4
  count = nx * ny * nz
  call MPI_FILE_WRITE_AT_ALL(fh,offset,tmp4,count,MPI_REAL4,mpi_status,mpi_err)
  call MPI_FILE_CLOSE(fh,mpi_err)
  call MPI_INFO_FREE(mpi_info,mpi_err)

!======================================================================  

  return
end subroutine write_tmp4_all

