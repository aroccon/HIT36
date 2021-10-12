subroutine restart_write_parallel
use m_openmpi
use m_parameters
use m_fields
use m_work
use x_fftw
implicit none
integer :: n, nums_out, i, j, k
integer :: MST
integer(kind=MPI_INTEGER_KIND) :: fh
integer(kind=MPI_OFFSET_KIND)  :: offset
real*8, allocatable :: sctmp8(:,:,:)
integer*4 :: nx1, ny1, nz1, nums1, nles1
real*8 :: ST

if (itime.eq.last_dump) return
! how many scalars to write
nums_out = 0
if (int_scalars) nums_out =  n_scalars

! using wrk array
wrk(:,:,:,1:3+nums_out) = fields(:,:,:,1:3+nums_out)

fname = run_name//'.64.'//file_ext

  ! allocating the temporary array sctmp8
  allocate(sctmp8(nx,ny,nz),stat=ierr)
  if (ierr.ne.0) stop '*** RESTART_READ_PARALLEL: cannot allocate sctmp8'
  sctmp8 = zip

  ! opening the file
  call MPI_INFO_CREATE(mpi_info, mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_TASK,fname,MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,fh,mpi_err)

  ! the master node writes the header with parameters
  if (myid.eq.0) then
     nx1 = nx;  ny1 = ny;  nz1 = nz_all;  nums1 = nums_out;
     count = 1
     ST = zip
     call MPI_FILE_WRITE(fh,   nx1, count, MPI_INTEGER4, mpi_status, mpi_err)
     call MPI_FILE_WRITE(fh,   ny1, count, MPI_INTEGER4, mpi_status, mpi_err)
     call MPI_FILE_WRITE(fh,   nz1, count, MPI_INTEGER4, mpi_status, mpi_err)
     call MPI_FILE_WRITE(fh, nums1, count, MPI_INTEGER4, mpi_status, mpi_err)
!     call MPI_FILE_WRITE(fh, nles1, count, MPI_INTEGER4, mpi_status, mpi_err)
     call MPI_FILE_WRITE(fh,  TIME, count, MPI_REAL8, mpi_status, mpi_err)
     call MPI_FILE_WRITE(fh,    DT, count, MPI_REAL8, mpi_status, mpi_err)
  end if

  ! all nodes write their stuff into the file
  ! first writing the velocities and passive scalars
  writing_fields: do n = 1, 3 + nums_out

     offset = 36 + (n-1)*nx*ny*nz_all*8 + myid*nx*ny*nz*8
     count = nx * ny * nz
     ! note that we want the data from the restart to have dimensions (nx,ny,nz),
     ! while the fields array has fimensions (nx+2,ny,nz).
     ! this is an artefact of the times when the code used to write the variables in real space.
     ! that is why we need to duplicate each field in the sctmp array first, and then
     ! write sctmp8 into the file with appropriate offset

     sctmp8(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,n)

     call MPI_FILE_WRITE_AT(fh, offset, sctmp8, count, MPI_REAL8, mpi_status, mpi_err)

  end do writing_fields


  call MPI_FILE_CLOSE(fh, mpi_err)
  call MPI_INFO_FREE(mpi_info, mpi_err)
  deallocate(sctmp8)


  last_dump = ITIME

  return
end subroutine restart_write_parallel




subroutine restart_read_parallel
use m_openmpi
use m_parameters
use m_fields
use m_work
use x_fftw
use m_particles
implicit none
integer*4    :: nx1,ny1,nz1, nums1, nles1, nums_read
integer      :: i, j, k, n, n_skip

  integer(kind=MPI_INTEGER_KIND) :: fh
  integer(kind=MPI_OFFSET_KIND)  :: offset
  real*8, allocatable :: sctmp8(:,:,:)

  ! checking if the restart file exists
  fname = run_name//'.64.'//file_ext
  inquire(file=fname,exist=there)
  if(.not.there) then
!     write(out,*) '*** error: Cannot find file : '//trim(fname)
     stop
  end if

!  write(out,*) 'Reading from the file (par): ',trim(fname)
!  call flush(out)

  ! ----------------------------------------------------------------------
  ! first reading the parameters from the restart file.
  ! the root process opens it and reads the parameters, then broadcasts
  ! the parameters.  After that it's decided if the parameters make sense,
  ! how many scalars to read etc.
  ! ----------------------------------------------------------------------

  if (myid.eq.0) then
     open(91,file=fname,form='unformatted',access='stream')
     read(91) nx1, ny1, nz1, nums1, nles1, TIME, DT
     close(91)
  end if

  call MPI_BCAST(nx1,  1,MPI_INTEGER4,0,MPI_COMM_TASK,mpi_err)
  call MPI_BCAST(ny1,  1,MPI_INTEGER4,0,MPI_COMM_TASK,mpi_err)
  call MPI_BCAST(nz1,  1,MPI_INTEGER4,0,MPI_COMM_TASK,mpi_err)
  call MPI_BCAST(nums1,1,MPI_INTEGER4,0,MPI_COMM_TASK,mpi_err)
  call MPI_BCAST(nles1,1,MPI_INTEGER4,0,MPI_COMM_TASK,mpi_err)

  call MPI_BCAST(TIME,1,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
  call MPI_BCAST(  DT,1,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)

  ! checking if the array sizes are the same in .in file and restart file
  if (nx.ne.nx1 .or. ny.ne.ny1 .or. nz_all.ne.nz1) then
!     write(out,*) '*** error: Dimensions are different'
!     write(out,*) '***     .in file: ',nx,ny,nz_all
!     write(out,*) '*** restart file: ',nx1,ny1,nz1
!     call flush(out)
     call my_exit(-1)
  end if



!-----------------------------------------------------------------------
!     dealing with scalars.
!
!     The number of scalars can be varied throughout the simulation.
!     If the restart file has fewer scalars than the
!     .in file, the scalars are added and initialized according to
!     their description in the .in-file.
!     If the restart file has more scalars than the .in file, the
!     extra scalars are dropped.
!
!     in short, whatever is specfied in the .in file, prevails.
!-----------------------------------------------------------------------

  if (n_scalars.lt.nums1) then

  !   write(out,*) ' WARNING: nums in restart file:',nums1
  !   write(out,*) '          nums in .in file    :',n_scalars
  !   write(out,'(''Losing '',i3,'' scalars.'')') nums1-n_scalars
  !   call flush(out)
     nums_read = n_scalars

  else if (n_scalars.gt.nums1) then

  !   write(out,*) ' WARNING: nums in restart file:',nums1
  !   write(out,*) '          nums in .in file    :',n_scalars
  !   write(out,'(''Adding '',i3,'' scalars.'')') n_scalars-nums1
  !   call flush(out)
     nums_read = nums1

     ! initializing the added scalars
     do n=nums1+1,n_scalars
        call init_scalar(n)
     end do

  else

     nums_read = n_scalars

  end if

!----------------------------------------------------------------------
!  ------------ reading the stuff from the restart file ---------------
!----------------------------------------------------------------------

  ! allocating the temporary array sctmp8
  allocate(sctmp8(nx,ny,nz),stat=ierr)
  if (ierr.ne.0) stop '*** RESTART_READ_PARALLEL: cannot allocate sctmp8'
  sctmp8 = zip

  ! opening the file
  call MPI_INFO_CREATE(mpi_info, mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_TASK,fname,MPI_MODE_RDONLY,mpi_info,fh,mpi_err)

  ! note that the data from the restart file has dimensions (nx,ny,nz),
  ! while the fields array has fimensions (nx+2,ny,nz).
  ! that is why we need to read each field in the sctmp array first, and then
  ! rearrange it and put into the fields array.

  reading_fields: do n = 1, 3 + nums_read

  !   write(out,"('Reading variable # ',i3)") n
  !   call flush(out)

     offset = 36 + (n-1)*nx*ny*nz_all*8 + myid*nx*ny*nz*8
     count = nx * ny * nz
!     call MPI_FILE_READ_AT(fh, offset, sctmp8, count, MPI_REAL8, mpi_status, mpi_err)
     call MPI_FILE_READ_AT_ALL(fh, offset, sctmp8, count, MPI_REAL8, mpi_status, mpi_err)
     fields(1:nx,1:ny,1:nz,n) = sctmp8(1:nx,1:ny,1:nz)
  end do reading_fields


  call MPI_FILE_CLOSE(fh, mpi_err)
  call MPI_INFO_FREE(mpi_info, mpi_err)
  deallocate(sctmp8)

!  write(out,*) "Restart file successfully read."
!  call flush(out)

  return
end subroutine restart_read_parallel
