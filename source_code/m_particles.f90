!================================================================================
module m_particles
!================================================================================
!  Incompressible pseudo-spectral code with scalars
!
!  The module contains subroutines for Largangian particle tracking.
!
!
!================================================================================
!                                    VARIABLES
!================================================================================

  use m_parameters
  use m_io
  implicit none

  ! logical variable that indicates if we integrate the particles
!  logical :: int_particles

  ! when to release the particles
  real*8 :: time_p

  ! array which contains the particles' positions in the grid units
  real*8, allocatable :: xyzp(:,:)

  ! array with second-layer velocities for particles
  real*8, allocatable :: uvwp(:,:)

  ! array with particles' tags (1 to nptot)
  integer*4, allocatable :: ipart(:)

  ! array with particles' addresses (which processors hold the particle)
  integer*4, allocatable :: myid_part(:), itmp_part(:)

  ! min and max z-coordinate for a slab
  real*8 :: zmin,zmax

  ! temp arrays for particle interpolation
  real*8, allocatable :: wrk1p(:,:), wrk2p(:,:), wrk3p(:,:)

  ! tmp arrays for CINT interpolation procedure (see particles_move_cint.f90)
!  real*8, allocatable :: 

  ! needed to avoid deadlock in send/receive particles
  real*8, allocatable :: xyzp1(:,:), uvwp1(:,:)
  integer*4, allocatable :: ipart1(:)

  ! number of particles to be sent/received to/from neighbours when moving
  ! the particles
  integer*4  :: np_send_u, np_send_d, np_get_u, np_get_d

  ! the last iteration at which the particles were written out
  integer :: particles_last_dump = -1

!================================================================================
!================================================================================

contains

!================================================================================
!                                   SUBROUTINES
!================================================================================

  subroutine particles_allocate

    use m_io
    implicit none

    write(out,*) 'Particle allocation'
    call flush(out)

    ! allocate the arrays
    if (.not.allocated(xyzp)) then
       allocate(xyzp(3,nptot),uvwp(3,nptot),ipart(nptot),itmp_part(nptot),&
            myid_part(nptot),&
            xyzp1(3,nptot),uvwp1(3,nptot),ipart1(nptot),stat=ierr)
       if (ierr.ne.0) stop 'particle allocation'
    else
       write(out,*) 'arrays xyzp etc are already allocated.'
       call flush(out)
    end if


    ! allocating temp arrays for CINT interpolation
    if (particles_tracking_scheme.eq.2) then
       write(out,*) 'Allocating temp arrays for CINT interpolation'
       call flush(out)
       if (.not.allocated(wrk1p)) then
          allocate(wrk1p(nx,ny), wrk2p(nx,ny), wrk3p(nx,ny),stat=ierr)
          wrk1p = zip; wrk2p = zip; wrk3p = zip;
          if (ierr.ne.0) stop 'allocation of particle temp arrays'
       else
          write(out,*) 'arrays for CINT are already allocated.'
          call flush(out)
       end if
    end if

    return
  end subroutine particles_allocate

!================================================================================
!================================================================================
!================================================================================

  subroutine particles_deallocate
    implicit none

    write(out,*) 'Particle deallocation.'
    call flush(out)

    ! allocate the arrays
    if (allocated(xyzp)) deallocate(xyzp,uvwp,ipart,itmp_part,myid_part)

    if (allocated(wrk1p)) deallocate(wrk1p,wrk2p,wrk3p,xyzp1,uvwp1,ipart1)

    return
  end subroutine particles_deallocate


!================================================================================
!================================================================================
!================================================================================
!   PARTICLES_GET_MYID - fills the array myid_part, which contains the IDs of
!                       the processes that house the particles.
!================================================================================
  subroutine particles_get_myid


    implicit none

    integer :: n

    ! filling up the array which contains id's of processes that are 
    ! responsible for each particle

    itmp_part = 0
    myid_part = 0

    if (np.gt.0) then
       do n = 1,np
          itmp_part( ipart(n) ) =  myid
       end do
    end if

    call MPI_REDUCE(itmp_part,myid_part,nptot,MPI_INTEGER4,MPI_SUM,master,MPI_COMM_TASK,ierr)

    call MPI_BCAST(myid_part,nptot,MPI_INTEGER4,master,MPI_COMM_TASK,ierr)

    return
  end subroutine particles_get_myid

!================================================================================
!================================================================================
!================================================================================
!  PARTICLES_INIT  - initializes the particles
!================================================================================
  subroutine particles_init
!--------------------------------------------------------------------------------
!  particle initialization: 
!   - aray allocation
!   - reading coordinates of particles from the file <run_name>.pt
!   - shifing the particles so hey are in [0:2*pi)
!   - distributing the particles among CPUs
!--------------------------------------------------------------------------------

    use m_filter_xfftw

    implicit none

    logical :: init_start = .false., init_internally = .false.
    integer :: i,j,k,npthird,n
    real*8  :: sctmp

!--------------------------------------------------------------------------------

    ! first, some foolproofing
    int_particles = .false.

    ! if the total # of particles is zero, return without initialization
    if (nptot.eq.0) then
       write(out,*) ' nptot = 0, no particles initialized.'
       call flush(out)
       return
    end if

!--------------------------------------------------------------------------------
!  choose which way to initialize particles:
!  1. read from file <run_name>.pt (manually given coordinates)
!  2. generate uniformly spaced particles in the computational domain
!  (in this case the number of particles is taken to be nptot=2^n, where
!  n is such that nptot is close to what is specified in the input file)
!--------------------------------------------------------------------------------


    ! if the file pt.<file_Ext> is there, then read particles from it and return.
    ! whatever is specified in the pt.<file_ext>, supercedes anything else.

    inquire(file='pt.'//file_ext,exist=there)
    if (there) then
       write(out,*) 'Detected file  pt.'//file_ext
       write(out,*) 'Particles will be read from the file.'
       call particles_restart_read_binary
       return
    end if

    ! if the restart file is not there, we assume that it's a
    ! brand new particle simulation so we read the particle
    ! coordinates from the "run_name".pt file
    ! ---------  OR  ----------
    ! define particles as uiformly distributed over the domain


    ! check if the <run_name>.pt is in the directory 
    inquire(file=run_name//'.pt',exist=there)
    if (there) then
       write(out,*) 'Reading particles from the file '//run_name//'.pt'
       call flush(out)
       init_start = .true.
    else
       ! if not reading from the file, redefine nptot
       init_internally = .true.
       npthird = int(dble(nptot)**(1.d0/3.d0))
       if (nptot-npthird**3.ge.(npthird+1)**3-nptot) npthird = npthird + 1
       nptot = npthird**3
       nptot = min(nptot,2**30)
       write(out,*) 'Redefined nptot to be a perfect cube:',nptot
       if (nptot.eq.2**30) write(out,*) "REACHED MAXIMUM: 2**30"
       call flush(out)
    end if


    write(out,*) 'Initializing',nptot,' particles'
    call flush(out)

    call particles_allocate


    if (myid.eq.0) then
       ! read the whole array of particles
       if (init_start) then
          open(99,file=run_name//'.pt')
          read(99,*,ERR=9000,END=9000) ((xyzp(i,j),i=1,3),j=1,nptot)
          close(99)
       end if

       ! ... or define it manually

       sctmp = two*PI / real(npthird,8)

       if (init_internally) then
          n = 0
          do k = 1,npthird
             do j = 1,npthird
                do i = 1,npthird
                   n = n + 1
                   xyzp(1,n) = (dble(i-1)+half) * sctmp
                   xyzp(2,n) = (dble(j-1)+half) * sctmp
                   xyzp(3,n) = (dble(k-1)+half) * sctmp
                end do
             end do
          end do
       end if

    end if


    ! Broadcast the array of particles
    call MPI_BCAST(xyzp,nptot*3,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)



    ! definition of particles' identifiers
    do i = 1,nptot
       ipart(i) = i
    end do

!!$  write(out,*) 'Particles read in:'
!!$  write(out,'(3f20.10)') ((xyzp(i,j),i=1,3),j=1,nptot)
!!$  call flush(out)

    ! REARRANGING THE PARTICLES

    ! makign sure the xyz are in [0,2*pi)
    do j=1,nptot
       do i = 1,3
          if (xyzp(i,j) .ge. 2.0d0*PI .or. xyzp(i,j) .lt. 0.0d0) then
             xyzp(i,j) = xyzp(i,j) - 2.0d0*PI * dble(floor(xyzp(i,j)/(2.0d0*PI)))
          end if
       end do
    end do

!!$  write(out,*) 'Particles normalized:',nptot
!!$  write(out,'(i5,3f20.10)') ((ipart(j),(xyzp(i,j),i=1,3)),j=1,nptot)

    ! removing particles that do not correspond to this slab
    np = nptot
    j = 1
    zmin = dble(myid*nz)     * 2.0d0*PI / dble(nz_all)
    zmax = dble((myid+1)*nz) * 2.0d0*PI / dble(nz_all)

    do while (j.le.np)
       if (xyzp(3,j).ge.zmin .and. xyzp(3,j).lt.zmax) then
          j = j + 1
       else
          xyzp(:,j) = xyzp(:,np)
          ipart(j) = ipart(np)
          np = np-1
       end if
    end do

!!$  write(out,*) 'Particles left in:',np
!!$  write(out,'(i5,3f20.10)') ((ipart(j),(xyzp(i,j),i=1,3)),j=1,np)

    ! Recalculating the particles' coordinates in grid cells
    do j = 1,np
       xyzp(1,j) = xyzp(1,j) / (2.0d0 * PI / dble(nx)) + 1.0d0
       xyzp(2,j) = xyzp(2,j) / (2.0d0 * PI / dble(ny)) + 1.0d0
       xyzp(3,j) = xyzp(3,j) / (2.0d0 * PI / dble(nx)) + 1.0d0 - dble(myid*nz)
    end do
    write(out,*) 'Particles left in, rescaled:',np
!!$  write(out,'(i5,3f20.10)') ((ipart(j),(xyzp(i,j),i=1,3)),j=1,np)


    ! initializing uvwp with -999.00
    uvwp = -999.0d0

!!$    ! making a subdirectory "particles"
!!$    if (iammaster) call system('mkdir -p particles')

    ! in case if the particles are advected by locally averaged field
    ! (particles_filter_size > 0), define the filter
    if (particles_filter_size > 0.d0) call filter_xfftw_init

    return
9000 write(out,*) '*** PARTICLES_INIT: ERROR in reading file <run_title>.pt !!!'
    stop
  end subroutine particles_init

!===========================================================================
!===========================================================================
!  Subroutine that writes out the particles coordinates in a file
!  during the restart write
!===========================================================================
  subroutine particles_restart_write


    implicit none

    integer(kind=MPI_INTEGER_KIND) :: id, count, my_status(MPI_STATUS_SIZE)
    integer*4 :: np_rec
    integer   :: np_cur, ntmp


    if (.not.int_particles) return

    write(out,*) 'writing the particle restart file'
    call flush(out)

    ! ------------------------------------------------------------------------------
    ! the slave nodes send their particles' IDs and coordinates to the master node
    ! the master node assembles the particles IDs and coordinates in one array and
    ! writes them out
    ! ------------------------------------------------------------------------------
    if (myid.ne.0) then

       ! First send the number of particles that the process holds to the master process
       np_rec = np
       call MPI_SEND(np_rec,1,MPI_INTEGER4,0,3*myid,MPI_COMM_TASK,mpi_err)

!    write(out,*) 'Sent np=',np_rec,' to master',mpi_err
!    call flush(out)

       ! if the number of particles is positive, send their IDs and coordinates
       if (np.gt.0) then
          count = np
          call MPI_SEND(ipart(1),count,MPI_INTEGER4,0,3*myid+1,MPI_COMM_TASK,mpi_err)

!      write(out,*) 'Sent ipart to master',mpi_err
!      call flush(out)

          count = 3*np
          call MPI_SEND(xyzp(1,1),count,MPI_REAL8,0,3*myid+2,MPI_COMM_TASK,mpi_err)

!      write(out,*) 'Sent xyzp to master',mpi_err
!      call flush(out)

       end if

    else

       np_cur = np+1
       do id = 1,numprocs-1
          ! master node receives the number of particles from slave node
          call MPI_RECV(np_rec,1,MPI_INTEGER4,id,3*id,MPI_COMM_TASK,my_status,mpi_err)

!      write(out,*) 'Received ',np_rec,' from ',id,':',mpi_err
!      call flush(out)

          ! if there are any particles to receive, receive the IDs and coordinates
          if (np_rec.gt.0) then
             count = np_rec
             call MPI_RECV(ipart(np_cur),count,MPI_INTEGER4,id,3*id+1,MPI_COMM_TASK,my_status,mpi_err)

!        write(out,*) 'Received ipart from ',id,':',mpi_err
!        call flush(out)

             count = 3 * np_rec
             call MPI_RECV(xyzp(1,np_cur),count,MPI_REAL8,id,3*id+2,MPI_COMM_TASK,my_status,mpi_err)

!        write(out,*) 'Received xyzp from ',id,':',mpi_err
!        call flush(out)

             ! since the xyzp contains particle coordinates in terms of the cells, the third
             ! coordinate should be augmented by id*nz.  this is done so later we can just
             ! multiply the third coordinate by dz to get the real particle coordinate
             do ntmp = np_cur,np_cur+np_rec
                xyzp(3,ntmp) = xyzp(3,ntmp) + dble(nz*id)
             end do
          end if
          np_cur = np_cur + np_rec

       end do
    end if

    ! ------------------------------------------------------------------------------
    ! Now the master node opens the particle restart file <run_name>.pt.<file_ext>
    ! and writes the particles with their IDs in it. 
    ! ------------------------------------------------------------------------------

    if (myid.eq.0) then
       open(89,file=run_name//'.pt.'//file_ext)
       do ntmp = 1,nptot
          write(89,'(i6,x,10d16.8)') ipart(ntmp),&
               (xyzp(1,ntmp)-1.d0)*dx,(xyzp(2,ntmp)-1.d0)*dy,(xyzp(3,ntmp)-1.d0)*dz
       end do
       close(89,status='keep')
       write(out,*) 'particle restart file written.'
       call flush(out)
    end if


    return
  end subroutine particles_restart_write


!===========================================================================
!===========================================================================
!  Subroutine that reads in the particles from the particle restart file
!===========================================================================
  subroutine particles_restart_read


    implicit none

    logical :: there
    character*80 :: fname
    integer :: n,i,j


    if (nptot.eq.0) return

    ! allocate the arrays
    if (.not.allocated(xyzp)) then
       allocate(xyzp(3,nptot),uvwp(3,nptot),ipart(nptot),itmp_part(nptot),&
            myid_part(nptot),stat=ierr)
       if (ierr.ne.0) stop 'particle allocation'
    else
       write(out,*) 'PARTICLES_RESTART_READ: xyzp allocated already!'
       call flush(out)
    end if

    fname = run_name//'.pt.'//file_ext
    write(out,*) 'Reading the particles from file '//fname
    call flush(out)

!-------------------------------------------------------------------------
!  Reading and broadcasting the particles' IDs and coordinates
!-------------------------------------------------------------------------

    inquire(file=fname,exist=there)
    if(.not.there) then
       write(out,*) 'could not find the particle restart file ',fname
       write(out,*) 'exiting.'
       call flush(out)
       call my_exit(-1)
    end if

    if (myid.eq.0) then
       open(81,file=fname)
       do n = 1,nptot
          read(81,*,ERR=9000,END=9000) ipart(n),xyzp(1:3,n)
       end do
       close(81)
    end if

    call MPI_BCAST(ipart,nptot,MPI_INTEGER4,0,MPI_COMM_TASK,mpi_err)
    call MPI_BCAST(xyzp,nptot*3,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)


    write(out,*) 'Particles normalized:',nptot
    do j = 1,nptot
       write(out,'(i5,3f20.10)') ipart(j),(xyzp(i,j),i=1,3)
    end do

!-------------------------------------------------------------------------
!  removing particles that do not correspond to this slab
!-------------------------------------------------------------------------

    write(out,*) 'removing particles that do not correspond to this slab'
    call flush(out)

    np = nptot
    j = 1
    zmin = dble(myid*nz)     * 2.0d0*PI / dble(nz_all)
    zmax = dble((myid+1)*nz) * 2.0d0*PI / dble(nz_all)

    do while (j.le.np)
       if (xyzp(3,j).ge.zmin .and. xyzp(3,j).lt.zmax) then
          j = j + 1
       else
          xyzp(:,j) = xyzp(:,np)
          ipart(j) = ipart(np)
          np = np-1
       end if
    end do

    write(out,*) 'Particles left in:',np
    write(out,*) 'Particles normalized:',nptot
    do j = 1,nptot
       write(out,'(i5,3f20.10)') ipart(j),(xyzp(i,j),i=1,3)
    end do

    ! Recalculating the particles' coordinates in grid cells
    do j = 1,np
       xyzp(1,j) = xyzp(1,j) / (2.0d0 * PI / dble(nx)) + 1.0d0
       xyzp(2,j) = xyzp(2,j) / (2.0d0 * PI / dble(ny)) + 1.0d0
       xyzp(3,j) = xyzp(3,j) / (2.0d0 * PI / dble(nx)) + 1.0d0 - dble(myid*nz)
    end do
    write(out,*) 'Particles left in, rescaled:',np
    write(out,*) 'Particles normalized:',nptot
    do j = 1,nptot
       write(out,'(i5,3f20.10)') ipart(j),(xyzp(i,j),i=1,3)
    end do


    ! initializing uvwp with -999.00
    uvwp = -999.0d0

    return
9000 write(out,*) 'Error reading particle restart file ',fname
    stop
  end subroutine particles_restart_read
!================================================================================
!================================================================================
!================================================================================
!  Subroutine that writes out the particles coordinates in a BINARY restart file
!================================================================================
  subroutine particles_restart_write_binary


    implicit none

    integer(kind=MPI_INTEGER_KIND) :: id, count, my_status(MPI_STATUS_SIZE)
    integer*4 :: np_rec
    integer   :: np_cur, ntmp

    ! if the last iteration at which the particles were written out is the
    ! current iteration, then skip the writing
    if (particles_last_dump .eq. itime) return

    ! otherwise proceed and redefine the particles_last_dump
    particles_last_dump = itime

!!$    write(out,*) 'writing the BINARY particle restart file'
!!$    call flush(out)

    ! ------------------------------------------------------------------------------
    ! the slave nodes send their particles' IDs and coordinates to the master node
    ! the master node assembles the particles IDs and coordinates in one array and
    ! writes them out
    ! ------------------------------------------------------------------------------
    if (myid.ne.0) then

       ! First send the number of particles that the process holds to the master process
       np_rec = np
       call MPI_SEND(np_rec,1,MPI_INTEGER4,0,3*myid,MPI_COMM_TASK,mpi_err)

       ! if the number of particles is positive, send their IDs and coordinates
       if (np.gt.0) then
          count = np
          call MPI_SEND(ipart(1),count,MPI_INTEGER4,0,3*myid+1,MPI_COMM_TASK,mpi_err)

          count = 3*np
          call MPI_SEND(xyzp(1,1),count,MPI_REAL8,0,3*myid+2,MPI_COMM_TASK,mpi_err)

       end if

    else

       np_cur = np+1
       do id = 1,numprocs-1
          ! master node receives the number of particles from slave node
          call MPI_RECV(np_rec,1,MPI_INTEGER4,id,3*id,MPI_COMM_TASK,my_status,mpi_err)

          ! if there are any particles to receive, receive the IDs and coordinates
          if (np_rec.gt.0) then
             count = np_rec
             call MPI_RECV(ipart(np_cur),count,MPI_INTEGER4,id,3*id+1,MPI_COMM_TASK,my_status,mpi_err)

             count = 3 * np_rec
             call MPI_RECV(xyzp(1,np_cur),count,MPI_REAL8,id,3*id+2,MPI_COMM_TASK,my_status,mpi_err)

             ! since the xyzp contains particle coordinates in terms of the cells, the third
             ! coordinate should be augmented by id*nz.  this is done so later we can just
             ! multiply the third coordinate by dz to get the real particle coordinate
             do ntmp = np_cur,np_cur+np_rec
                xyzp(3,ntmp) = xyzp(3,ntmp) + dble(nz*id)
             end do
          end if
          np_cur = np_cur + np_rec

       end do
    end if

    ! ------------------------------------------------------------------------------
    ! Now the master node opens the particle BINARY file pt.<file_ext>
    ! and writes the particles with their IDs in it. 
    ! ------------------------------------------------------------------------------

    if (myid.eq.0) then
!!       open(89,file='pt.'//file_ext,form='binary')
       open(89,file='pt.'//file_ext,form='unformatted',access='stream')
       write(89) nptot
       do ntmp = 1,nptot
          write(89) ipart(ntmp),&
               (xyzp(1,ntmp)-1.d0)*dx,(xyzp(2,ntmp)-1.d0)*dy,(xyzp(3,ntmp)-1.d0)*dz
       end do
       ! writing time stamp in the particle output file
       write(89) TIME

       close(89,status='keep')
       write(out,'("particle BINARY file written:",i7,e15.6)') itime, time
       call flush(out)
    end if


    return
  end subroutine particles_restart_write_binary
!================================================================================


!===========================================================================
!  Subroutine that reads in the particles from the particle restart file
!===========================================================================
  subroutine particles_restart_read_binary

    implicit none

    integer(kind=MPI_INTEGER_KIND) :: nptot1

    logical :: there
    character*80 :: fname
    integer :: n,i,j


    fname = 'pt.'//file_ext
    write(out,*) 'Reading the particles from binary file '//fname
    call flush(out)

!-------------------------------------------------------------------------
!  Reading and broadcasting the particles' IDs and coordinates
!-------------------------------------------------------------------------

!  Checking that the number of particles in the binary file is correct  
    if (myid.eq.0) then
       open(81,file=fname,form='unformatted',access='stream')
       read(81) nptot
    end if
    call MPI_BCAST(nptot,1,MPI_INTEGER,0,MPI_COMM_TASK,mpi_err)

    write(out,*) 'nptot = ',nptot
    call flush(out)

    ! allocating arrays
    call particles_allocate


    ! Readign and broadcasting particles  
    if (myid.eq.0) then
       do n = 1,nptot
          read(81,ERR=9000,END=9000) ipart(n),xyzp(1:3,n)
       end do
       close(81)
    end if

    call MPI_BCAST(ipart,nptot,MPI_INTEGER4,0,MPI_COMM_TASK,mpi_err)
    call MPI_BCAST(xyzp,nptot*3,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)


    write(out,*) 'Particles normalized:',nptot
!!$    write(out,'(i5,3f20.10)') ((ipart(j),(xyzp(i,j),i=1,3)),j=1,nptot)

!-------------------------------------------------------------------------
!  removing particles that do not correspond to this slab
!-------------------------------------------------------------------------

    write(out,*) 'removing particles that do not correspond to this slab'
    call flush(out)

    np = nptot
    j = 1
    zmin = dble(myid*nz)     * 2.0d0*PI / dble(nz_all)
    zmax = dble((myid+1)*nz) * 2.0d0*PI / dble(nz_all)

    do while (j.le.np)
       if (xyzp(3,j).ge.zmin .and. xyzp(3,j).lt.zmax) then
          j = j + 1
       else
          xyzp(:,j) = xyzp(:,np)
          ipart(j) = ipart(np)
          np = np-1
       end if
    end do

    write(out,*) 'Particles left in:',np
!!$    write(out,'(i5,3f20.10)') ((ipart(j),(xyzp(i,j),i=1,3)),j=1,np)

    ! Recalculating the particles' coordinates in grid cells
    do j = 1,np
       xyzp(1,j) = xyzp(1,j) / (2.0d0 * PI / dble(nx)) + 1.0d0
       xyzp(2,j) = xyzp(2,j) / (2.0d0 * PI / dble(ny)) + 1.0d0
       xyzp(3,j) = xyzp(3,j) / (2.0d0 * PI / dble(nx)) + 1.0d0 - dble(myid*nz)
    end do
!!$    write(out,*) 'Particles left in, rescaled:',np
!!$    write(out,'(i5,3f20.10)') ((ipart(j),(xyzp(i,j),i=1,3)),j=1,np)


    ! initializing uvwp with -999.00
    uvwp = -999.0d0

    return
9000 write(out,*) 'Error reading particle restart file ',fname
    stop
  end subroutine particles_restart_read_binary


!================================================================================
!================================================================================
!================================================================================
  subroutine particles_move
    use m_fields
    use m_work
    use x_fftw
    use m_filter_xfftw
    implicit none
    integer :: m, n

    if (task.ne.'parts') return


    ! if the particles are advected by locally averaged velocities
    ! instead of fully resolved velocities, then we need to filter
    ! (locally average) these velocities.  The indicator of this case is
    ! particles_filter_size > 0.d0
    ! in this case the velocities that are sent to the "parts" part of the code
    ! are given in the Fourier space.  So we need to 
    ! 1) filter them
    ! 2) transform them back to real space

    locally_averaged_velocity: if (particles_filter_size .gt. 0.d0) then
       ! filter each field and transform to real space
       do n = 1,3
          call filter_xfftw_fields(n)
          call xFFT3d_fields(-1,n)
       end do
    end if locally_averaged_velocity

    ! now shifting particles according to the interpolation scheme and advection scheme
    select case (particles_tracking_scheme)
    case (0)
       call particles_interpolate_trilinear
    case (1)
       call particles_interpolate_cint
    case default
       stop 'wrong particles_tracking_scheme'
    end select

    call particles_update_slabs

    return
  end subroutine particles_move


!================================================================================
!  PARTICLE VELOCITY INTERPOLATION AND MOVING
!
!  This is set up in the following way:
!  
!  The interpolation subroutines do the following:
!  1) find particles' velocities using some interpolation methof
!  2) move the particles according to the numerical scheme 
!     (Euler or Adams-Bashforth)
!  3) Calculate the number of particles that went out of the slab:
!     np_send_d  - # of particles that needs to be sent "down", to slab # myid-1
!     np_send_u  - # of particles that needs to be sent "up", to slab # myid+1
!
!  Then the subroutine particles_update_slabs should be called, which
!  1) rearranges particles according to their destination (here, down, up)
!  2) inform the neighbor slabs about the number of particles to send
!  3) transfer particles between the slabs
!
!  May the Force be with us.
!================================================================================
!================================================================================
  subroutine particles_interpolate_trilinear

    use m_fields
    use m_work

    implicit none

    integer :: i,j,k,i1,j1,k1,n
    real*8  :: xp,yp,zp,up,vp,wp
    real*8  :: c1,c2,c3,c4,c5,c6,c7,c8

    ! note that we assume that fields1...3 contain REAL velocities (in X-space)

    ! to interpolate velocities of the particles that are on the upper z-edge of 
    ! the slice, we need the first layer of velocities from the neibour process.  
    ! This is done by sending them to the process with number (myid-1) and storing
    ! in the first three slices (k=1,3) of wrk0

    id_to   = mod(myid-1+numprocs,numprocs)
    id_from = mod(myid+1,numprocs)

    ! sending u
    sendtag = 3 * myid
    recvtag = 3 * id_from

    count = (nx+2)*ny
    call MPI_IRECV(  wrk(1,1,1,0), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,1,1), count, MPI_REAL8, id_to,   sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    ! sending v
    sendtag = 3 * myid    + 1
    recvtag = 3 * id_from + 1

    count = (nx+2)*ny
    call MPI_IRECV(  wrk(1,1,2,0), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,1,2), count, MPI_REAL8, id_to,   sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    ! sending w
    sendtag = 3 * myid    + 2
    recvtag = 3 * id_from + 2

    count = (nx+2)*ny
    call MPI_IRECV(  wrk(1,1,3,0), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,1,3), count, MPI_REAL8, id_to,   sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    ! # of particles to be sent/received to/from left and right neibors
    np_send_d = 0
    np_send_u = 0
    np_get_d = 0
    np_get_u = 0


    ! interpolating the velocities and moving
    do n = 1,np

       ! getting the cell where the particle is at
       i = floor(xyzp(1,n))
       j = floor(xyzp(2,n))
       k = floor(xyzp(3,n))

       ! getting fractional coordinates of the particle
       xp = xyzp(1,n) - dble(i)
       yp = xyzp(2,n) - dble(j)
       zp = xyzp(3,n) - dble(k)

       ! cpefficients for trilinear interpolation
       c1 = (one-xp) * (one-yp) * (one-zp)
       c2 = (xp)     * (one-yp) * (one-zp)
       c3 = (one-xp) * (yp)     * (one-zp)
       c4 = (xp)     * (yp)     * (one-zp)
       c5 = (one-xp) * (one-yp) * (zp)
       c6 = (xp)     * (one-yp) * (zp)
       c7 = (one-xp) * (yp)     * (zp)
       c8 = (xp)     * (yp)     * (zp)

       ! accounting for periodicity in x and y directions
       if (i.eq.nx) then
          i1 = 1
       else
          i1 = i + 1
       end if

       if (j.eq.ny) then
          j1 = 1
       else
          j1 = j + 1
       end if

       ! interpolating velocities to get the velocity of the particle

       up = c1 * fields(i,j ,k,1) + c2 * fields(i1,j ,k,1) + &
            c3 * fields(i,j1,k,1) + c4 * fields(i1,j1,k,1)
       vp = c1 * fields(i,j ,k,2) + c2 * fields(i1,j ,k,2) + &
            c3 * fields(i,j1,k,2) + c4 * fields(i1,j1,k,2)
       wp = c1 * fields(i,j ,k,3) + c2 * fields(i1,j ,k,3) + &
            c3 * fields(i,j1,k,3) + c4 * fields(i1,j1,k,3)

       if (k.eq.nz) then
          ! if the particle is in the last layer, add velocities from wrk0
          up = up + &
               c5 * wrk(i,j ,1,0) + c6 * wrk(i1,j ,1,0) + &
               c7 * wrk(i,j1,1,0) + c8 * wrk(i1,j1,1,0)
          vp = vp + &
               c5 * wrk(i,j ,2,0) + c6 * wrk(i1,j ,2,0) + &
               c7 * wrk(i,j1,2,0) + c8 * wrk(i1,j1,2,0)
          wp = wp + &
               c5 * wrk(i,j ,3,0) + c6 * wrk(i1,j ,3,0) + &
               c7 * wrk(i,j1,3,0) + c8 * wrk(i1,j1,3,0)
       else
          ! if the particle is inside the slice, straightforward interpolation
          up = up + &
               c5 * fields(i,j ,k+1,1) + c6 * fields(i1,j ,k+1,1) + &
               c7 * fields(i,j1,k+1,1) + c8 * fields(i1,j1,k+1,1)
          vp = vp + &
               c5 * fields(i,j ,k+1,2) + c6 * fields(i1,j ,k+1,2) + &
               c7 * fields(i,j1,k+1,2) + c8 * fields(i1,j1,k+1,2)
          wp = wp + &
               c5 * fields(i,j ,k+1,3) + c6 * fields(i1,j ,k+1,3) + &
               c7 * fields(i,j1,k+1,3) + c8 * fields(i1,j1,k+1,3)
       end if

       ! now move the particle

       if (uvwp(1,n).eq.-999.0d0) then
          ! if the second slice is unavailable, use Euler
          xp = xp + dt * up / dx
          yp = yp + dt * vp / dy
          zp = zp + dt * wp / dz
       else
          ! if there is the second slice, use Adams-Bashforth
          xp = xp + (1.50d0*up - 0.50d0*uvwp(1,n)) * dt/dx
          yp = yp + (1.50d0*vp - 0.50d0*uvwp(2,n)) * dt/dy
          zp = zp + (1.50d0*wp - 0.50d0*uvwp(3,n)) * dt/dz
       end if

       ! update the velocities
       uvwp(1,n) = up
       uvwp(2,n) = vp
       uvwp(3,n) = wp

       ! apply periodicity in x- and y-directions
       ! also count the # of particles to be sent to neibour slices
       ! note that we assume that the particle cannot move further than one cell
       ! from its curent position

       ! x direction
       if (xp.gt.1.0d0) then
          xp = xp - 1.0d0
          i = i+1
          if (i.gt.nx) i = i - nx
       end if
       if (xp.lt.0.0d0) then
          xp = xp + 1.0d0
          i = i-1
          if (i.lt.1) i = i + nx
       end if

       ! y direction
       if (yp.gt.1.0d0) then
          yp = yp - 1.0d0
          j = j+1
          if (j.gt.ny) j = j - ny
       end if
       if (yp.lt.0.0d0) then
          yp = yp + 1.0d0
          j = j-1
          if (j.lt.1) j = j + ny
       end if

       ! z direction
       if (zp.gt.1.0d0) then
          zp = zp - 1.0d0
          k = k+1
          if (k.gt.nz) np_send_u = np_send_u + 1
       end if
       if (zp.lt.0.0d0) then
          zp = zp + 1.0d0
          k = k-1
          if (k.lt.1) np_send_d = np_send_d + 1
       end if


       ! updating the particle coordinates.  
       ! note that z-coordinate can be outside of the slab
       xyzp(1,n) = xp + dble(i)
       xyzp(2,n) = yp + dble(j)
       xyzp(3,n) = zp + dble(k)

    end do

    return
  end subroutine particles_interpolate_trilinear


!================================================================================
!================================================================================
!  Moving Lagrangian particles.  The velocity is interpolated using tricubic
!  interpolation method, found in Wikipedia  :)
!================================================================================
  subroutine particles_interpolate_cint

    use m_fields
    use m_work
    implicit none

    ! local variables
    integer :: ip,jp,kp,i,j,k,i1,j1,k1,n
    real*8  :: xp,yp,zp,up,vp,wp
    real*8  :: sctmp
    real*8 :: uloc(4,4,4), vloc(4,4,4), wloc(4,4,4)
!    real*8 :: cint_3d

    real*8 :: r

    ! making sure that the thickness of the slabs is mroe than 1 slice
    if (nz.lt.2) then
       write(out,*) "PARTICLES_MOVE_CINT: nz is less than 2, stopping"
       call flush(out)
       stop 'particles_move_cint *** nz less than 2'
    end if

    ! assuming that array fields1...3 has the velocities in X-space

    ! because this method involves the adjacent 26 cells (a 4x4 cube of points
    ! in which the particle occupies the central cell), we need to pass the
    ! first two slices "down" and the last single slice (nz) "up"

    ! sending/receiving the single slice 
    ! sending "up"
    ! receiving from "down"

    count = (nx+2) * ny
    id_from = mod(myid - 1 + numprocs, numprocs);  
    id_to   = mod(myid + 1           , numprocs);  

    sendtag = 6 * myid    + 0
    recvtag = 6 * id_from + 0
    call MPI_IRECV(  wrk(1,1,1, 1), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,nz,1), count, MPI_REAL8, id_to  , sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    sendtag = 6 * myid    + 1
    recvtag = 6 * id_from + 1
    call MPI_IRECV(  wrk(1,1,1, 2), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,nz,2), count, MPI_REAL8, id_to  , sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    sendtag = 6 * myid    + 2
    recvtag = 6 * id_from + 2
    call MPI_IRECV(  wrk(1,1,1, 3), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,nz,3), count, MPI_REAL8, id_to  , sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    ! sending/receiving two slicee 
    ! sending "down"
    ! receiving from "up"

    count = 2 * (nx+2) * ny
    id_from = mod(myid + 1           , numprocs);  
    id_to   = mod(myid - 1 + numprocs, numprocs);  

    sendtag = 6 * myid    + 3
    recvtag = 6 * id_from + 3
    call MPI_IRECV(  wrk(1,1,2,1), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,1,1), count, MPI_REAL8, id_to  , sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    sendtag = 6 * myid    + 4
    recvtag = 6 * id_from + 4
    call MPI_IRECV(  wrk(1,1,2,2), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,1,2), count, MPI_REAL8, id_to  , sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    sendtag = 6 * myid    + 5
    recvtag = 6 * id_from + 5
    call MPI_IRECV(  wrk(1,1,2,3), count, MPI_REAL8, id_from, recvtag, MPI_COMM_TASK, request, mpi_err)
    call MPI_SEND(fields(1,1,1,3), count, MPI_REAL8, id_to  , sendtag, MPI_COMM_TASK, mpi_err)
    call MPI_WAIT(request, mpi_status, mpi_err)

    ! # of particles to be sent/received to/from "down" and "up"
    np_send_u = 0
    np_send_d = 0
    np_get_u = 0
    np_get_d = 0

    do n = 1,np

       ! getting the cell where the particle is at
       ip = floor(xyzp(1,n))
       jp = floor(xyzp(2,n))
       kp = floor(xyzp(3,n))

       ! getting fractional coordinates of the particle
       xp = xyzp(1,n) - dble(ip)
       yp = xyzp(2,n) - dble(jp)
       zp = xyzp(3,n) - dble(kp)

       ! -------------------------------------------------------------------
       ! Getting the velocities of the particles by tricubic interpolation
       ! (for particularities, see notes from 12/13/07)

       ! fillin out arrays uloc, vloc, and wloc(4,4,4) for tricubic interpolation
       ! taking into account periodicity in x and y.  

       do k = 1,4
          k1 = kp - 2 + k

          do j = 1,4
             j1 = jp - 2 + j
             if (j1.lt.1 ) j1 = j1 + ny
             if (j1.gt.ny) j1 = j1 - ny

             do i = 1,4
                i1 = ip - 2 + i
                if (i1.lt.1 ) i1 = i1 + nx
                if (i1.gt.nx) i1 = i1 - nx

                if (k1.lt.1) then
                   uloc(i,j,k) = wrk(i1,j1,1,1)
                   vloc(i,j,k) = wrk(i1,j1,1,2)
                   wloc(i,j,k) = wrk(i1,j1,1,3)
                else if (k1.gt.nz) then
                   uloc(i,j,k) = wrk(i1,j1,k1-nz+1,1)
                   vloc(i,j,k) = wrk(i1,j1,k1-nz+1,2)
                   wloc(i,j,k) = wrk(i1,j1,k1-nz+1,3)
                else
                   uloc(i,j,k) = fields(i1,j1,k1,1)
                   vloc(i,j,k) = fields(i1,j1,k1,2)
                   wloc(i,j,k) = fields(i1,j1,k1,3)
                end if

             end do
          end do
       end do

       ! now we have arrays uloc, vloc, wloc and can interpolate

       up = cint_3d(uloc,xp,yp,zp)
       vp = cint_3d(vloc,xp,yp,zp)
       wp = cint_3d(wloc,xp,yp,zp)

       ! now move the particle

       if (uvwp(1,n).eq.-999.0d0) then
          ! if the second slice is unavailable, use Euler
          xp = xp + dt * up / dx
          yp = yp + dt * vp / dy
          zp = zp + dt * wp / dz
       else
          ! if there is the second slice, use Adams-Bashforth
          xp = xp + (1.50d0*up - 0.50d0*uvwp(1,n)) * dt/dx
          yp = yp + (1.50d0*vp - 0.50d0*uvwp(2,n)) * dt/dy
          zp = zp + (1.50d0*wp - 0.50d0*uvwp(3,n)) * dt/dz
       end if

       uvwp(1,n) = up
       uvwp(2,n) = vp
       uvwp(3,n) = wp


       ! apply periodicity in x- and y-directions
       ! also count the # of particles to be sent to neibour slices
       ! note that we assume that the particle cannot move further than one cell
       ! from its curent position

       ! x direction
       if (xp.gt.1.0d0) then
          xp = xp - 1.0d0
          ip = ip+1
          if (ip.gt.nx) ip = ip - nx
       end if
       if (xp.lt.0.0d0) then
          xp = xp + 1.0d0
          ip = ip-1
          if (ip.lt.1) ip = ip + nx
       end if

       ! y direction
       if (yp.gt.1.0d0) then
          yp = yp - 1.0d0
          jp = jp+1
          if (jp.gt.ny) jp = jp - ny
       end if
       if (yp.lt.0.0d0) then
          yp = yp + 1.0d0
          jp = jp-1
          if (jp.lt.1) jp = jp + ny
       end if

       ! z direction
       if (zp.gt.1.0d0) then
          zp = zp - 1.0d0
          kp = kp+1
          if (kp.gt.nz) np_send_u = np_send_u + 1
       end if
       if (zp.lt.0.0d0) then
          zp = zp + 1.0d0
          kp = kp-1
          if (kp.lt.1) np_send_d = np_send_d + 1
       end if


       ! updating the particle coordinates.  
       ! note that z-coordinate can be outside of the slice
       xyzp(1,n) = xp + dble(ip)
       xyzp(2,n) = yp + dble(jp)
       xyzp(3,n) = zp + dble(kp)

    end do


    return
  end subroutine particles_interpolate_cint



!================================================================================
!================================================================================
  subroutine particles_update_slabs

    use m_openmpi
    use m_work
    implicit none

    integer (kind=MPI_INTEGER_KIND) :: id_down, id_up

    integer :: i, j, n, iii


    id_down = mod(myid-1+numprocs,numprocs)
    id_up   = mod(myid+1,numprocs)

!!$    write(out,*) "np_sends:", np_send_d, np_send_u
!!$    call flush(out)

    ! let the neighbors know how many particles went to them
    call MPI_SEND(np_send_d,1,MPI_INTEGER4,id_down,2*myid  ,MPI_COMM_TASK,mpi_err)
    call MPI_SEND(np_send_u,1,MPI_INTEGER4,id_up  ,2*myid+1,MPI_COMM_TASK,mpi_err)

    ! learn how many particles will be transferred to this slice from neighbors
    call MPI_RECV(np_get_d,1,MPI_INTEGER4,id_down,2*id_down+1,MPI_COMM_TASK,mpi_status,mpi_err)
    call MPI_RECV(np_get_u,1,MPI_INTEGER4,id_up  ,2*id_up    ,MPI_COMM_TASK,mpi_status,mpi_err)


!!$    write(out,*) "np_gets:", np_get_d, np_get_u
!!$    call flush(out)

    ! rearranging particles in the following order:
    ! 1. particles that stay
    ! 2. particles that leave down
    ! 3. particles that leave up

    if (np_send_u.gt.0) then
       i = 1
       n = np
       do while (i.le.n)
          if (xyzp(3,i).ge.dble(nz+1)) then
             wrk(1:3,1,1,0) = xyzp(1:3,i)
             wrk(4:6,1,1,0) = uvwp(1:3,i)
             j = ipart(i)

             xyzp(:,i) = xyzp(:,n)
             uvwp(:,i) = uvwp(:,n)
             ipart(i) = ipart(n)

             xyzp(1:3,n) = wrk(1:3,1,1,0)
             uvwp(1:3,n) = wrk(4:6,1,1,0)
             xyzp(3,n) = xyzp(3,n) - dble(nz)
             ipart(n) = j

             n = n - 1

          else
             i = i + 1
          end if
       end do
    end if

    if (np_send_d.gt.0) then
       i = 1
       n = np - np_send_u
       do while (i.le.n)
          if (xyzp(3,i).lt.1.0d0) then
             wrk(1:3,1,1,0) = xyzp(1:3,i)
             wrk(4:6,1,1,0) = uvwp(1:3,i)
             j = ipart(i)

             xyzp(:,i) = xyzp(:,n)
             uvwp(:,i) = uvwp(:,n)
             ipart(i) = ipart(n)

             xyzp(1:3,n) = wrk(1:3,1,1,0)
             uvwp(1:3,n) = wrk(4:6,1,1,0)
             xyzp(3,n) = xyzp(3,n) + dble(nz)
             ipart(n) = j

             n = n - 1
          else
             i = i + 1
          end if
       end do
    end if


    ! sending and receiving particiles from adjacent slabs

    ! extra arrays (a dumb way to avoid deadlocking in send/receive when
    ! number of particles is large and the send/receive is not completed in the
    ! single handshake between the two processes)

    xyzp1  = xyzp
    uvwp1  = uvwp
    ipart1 = ipart

    np1 = np
    np = np - np_send_d - np_send_u

    ! posting receive from down
    if (np_get_d.gt.0) then
       call MPI_IRECV(xyzp(1,np+1),3*np_get_d,MPI_REAL8 ,id_down,6*id_down+0,MPI_COMM_TASK,request1,mpi_err)
       call MPI_IRECV(uvwp(1,np+1),3*np_get_d,MPI_REAL8 ,id_down,6*id_down+1,MPI_COMM_TASK,request2,mpi_err)
       call MPI_IRECV( ipart(np+1),np_get_d,MPI_INTEGER8,id_down,6*id_down+2,MPI_COMM_TASK,request3,mpi_err)

!!$      write(out,*) 'posted receive from ',id_down
!!$      call flush(out)

    end if

    ! sending particles to up
    if (np_send_u.gt.0) then

!!$      write(out,*) 'xyzp1'
!!$      write(out,"(3e15.6)") (xyzp1(:,iii),iii=1,np)
!!$      write(out,*) 'uvwp1'
!!$      write(out,"(3e15.6)") (uvwp1(:,iii),iii=1,np)
!!$      write(out,*) 'ipart1'
!!$      write(out,"(i5)") (ipart(iii),iii=1,np)
!!$      call flush(out)

       call MPI_SEND(xyzp1(1,np1-np_send_u+1),3*np_send_u,MPI_REAL8 ,id_up,6*myid+0,MPI_COMM_TASK,mpi_err)
       call MPI_SEND(uvwp1(1,np1-np_send_u+1),3*np_send_u,MPI_REAL8 ,id_up,6*myid+1,MPI_COMM_TASK,mpi_err)
       call MPI_SEND( ipart1(np1-np_send_u+1),np_send_u,MPI_INTEGER4,id_up,6*myid+2,MPI_COMM_TASK,mpi_err)

!!$      write(out,*) 'sent',np_send_u,' to the right:', ipart1(np1-np_send_u+1:np1)
!!$      call flush(out)

       np1 = np1 - np_send_u

    end if

    ! completing the receive from the down
    if (np_get_d.gt.0) then
       call MPI_WAIT(request1,mpi_status,mpi_err)
       call MPI_WAIT(request2,mpi_status,mpi_err)
       call MPI_WAIT(request3,mpi_status,mpi_err)

!!$      write(out,*) 'received',np_get_d,' from the left:', ipart(np+1:np+np_get_d)
!!$      call flush(out)

       np  = np  + np_get_d

    end if



    ! posting receive from the up
    if (np_get_u.gt.0) then
       call MPI_IRECV(xyzp(1,np+1),3*np_get_u,MPI_REAL8   ,id_up,6*id_up+3,MPI_COMM_TASK,request1,mpi_err)
       call MPI_IRECV(uvwp(1,np+1),3*np_get_u,MPI_REAL8   ,id_up,6*id_up+4,MPI_COMM_TASK,request2,mpi_err)
       call MPI_IRECV( ipart(np+1),np_get_u  ,MPI_INTEGER8,id_up,6*id_up+5,MPI_COMM_TASK,request3,mpi_err)

!!$      write(out,*) 'posted receive from ',id_up
!!$      call flush(out)

    end if


    ! sending particles to down
    if (np_send_d.gt.0) then
       call MPI_SEND(xyzp1(1,np1-np_send_d+1),3*np_send_d,MPI_REAL8   ,id_down,6*myid+3,MPI_COMM_TASK,mpi_err)
       call MPI_SEND(uvwp1(1,np1-np_send_d+1),3*np_send_d,MPI_REAL8   ,id_down,6*myid+4,MPI_COMM_TASK,mpi_err)
       call MPI_SEND( ipart1(np1-np_send_d+1),np_send_d  ,MPI_INTEGER4,id_down,6*myid+5,MPI_COMM_TASK,mpi_err)

!!$      write(out,*) 'sent',np_send_d,' to the left:', ipart(np1-np_send_d+1:np1)
!!$      call flush(out)

    end if

    ! completing the receive from up
    if (np_get_u.gt.0) then
       call MPI_WAIT(request1,mpi_status,mpi_err)
       call MPI_WAIT(request2,mpi_status,mpi_err)
       call MPI_WAIT(request3,mpi_status,mpi_err)

!!$      write(out,*) 'received',np_get_u,' from the right:', ipart(np+1:np+np_get_u)
!!$      call flush(out)

       np = np + np_get_u

    end if



    ! writing out particles' locations to separate files for each particle
!!$    if (mod(ITIME,IWRITE4).eq.0) then

    if (nptot.lt.1000) then
       do i=1,np

          write(fname,"('p.',i4.4)") ipart(i)
          open(99,file=fname,position='append')
          write(99,'(i6.6,x, i4, 12e16.8)') &
               itime,myid,time,(xyzp(1:2,i)-1.0d0)*dx,(xyzp(3,i)-1.0d0+dble(myid*nz))*dz,&
               uvwp(:,i)
          close(99)
       end do
    end if

!!$    end if




    return
  end subroutine particles_update_slabs


!================================================================================
!================================================================================
!================================================================================
!================================================================================

!================================================================================
  function cint_3d(a,x,y,z) result (blah)

    implicit none
    real*8 :: blah
    real*8 :: a(4,4,4), x, y, z

    logical :: compute_r
    integer :: j, k
    real*8 :: t2(4,4), t1(4), r(4)

    blah = 0.0d0

    compute_r = .true.
    do k = 1,4
       do j = 1,4
          t2(j,k) = cint(a(1,j,k),a(2,j,k),a(3,j,k),a(4,j,k),x,r,compute_r)
          if (compute_r) compute_r = .false.
       end do
    end do

    compute_r = .true.
    do k = 1,4
       t1(k) = cint(t2(1,k),t2(2,k),t2(3,k),t2(4,k),y,r,compute_r)
       if (compute_r) compute_r = .false.
    end do

    blah = cint(t1(1),t1(2),t1(3),t1(4),z,r,.true.)

    return
  end function cint_3d


!================================================================================

  function cint(p1,p2,p3,p4,x,r,compute_r) result (blah)
    implicit none
    logical :: compute_r
    real*8 :: p1,p2,p3,p4,x,r(4)

    real*8 :: sctmp

    real*8 :: blah

    ! If flag=.true., compute vector R(x)
    ! If flag=.false., consider vector R(x) defined and do not compute it
    if (compute_r) then
       sctmp = (x - 1.d0) * (x - 2.d0)
       r(1) = -x * sctmp
       r(2) =  3.d0 * (x + 1.d0) * sctmp
       sctmp = (x + 1) * x 
       r(3) = -3.d0 * sctmp * (x - 2.d0)
       r(4) = sctmp * (x - 1.d0)
    end if

    blah = (r(1)*p1 + r(2)*p2 + r(3)*p3 + r(4)*p4) / 6.0d0

    return
  end function cint



end module m_particles
