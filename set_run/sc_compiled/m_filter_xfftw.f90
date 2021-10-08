!======================================================================
!  Module that contains filtering rocedure for the part of the code that
!  addvects lagrangian particles
!
!  Time-stamp: <2009-05-20 11:50:52 (chumakov)>
!======================================================================
module m_filter_xfftw

  use m_parameters
  use x_fftw
  use m_work
  implicit none

  ! filter type (Gaussian by default)
  integer :: filter_type = 2

  ! filter size
  real*8 :: filter_size

  ! filter function in Fourier space
  real*8, allocatable :: filter_g(:,:,:)

!======================================================================
contains
!======================================================================

!======================================================================
!  filering subroutine for 2*pi^3-periodic domain
!  IN: ss(n,:,:,:)
!  OUT: ss(n:,;,:)
!======================================================================
  subroutine filter_xfftw(n)

    use m_io
    implicit none

    integer n

    integer :: ix, iy, iz

    real*8 :: a,b,c,d

!----------------------------------------------------------------------

!!$    call xFFT3d(1,n)

    ! multiplying wrk(:,:,:,n) by filter_g
    ! remember both are in the complex form
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx + 1, 2
             a = wrk(ix,iz,iy,n)
             b = wrk(ix+1,iz,iy,n)
             c = filter_g(ix,iz,iy)
             d = filter_g(ix+1,iz,iy)

             wrk(ix    , iz, iy, n) = a*c - b*d
             wrk(ix + 1, iz, iy, n) = b*c + a*d
          end do
       end do
    end do

!!$    ! perform inverse FFT
!!$    call xFFT3d(-1,n)

    return
  end subroutine filter_xfftw

!======================================================================
!  filering subroutine for 2*pi^3-periodic domain
!  IN: ss(n,:,:,:)
!  OUT: ss(n:,;,:)
!======================================================================
  subroutine filter_xfftw_fields(n)

    use m_io
    use m_fields
    implicit none

    integer n

    integer :: ix, iy, iz

    real*8 :: a,b,c,d

!----------------------------------------------------------------------

!!$    call xFFT3d(1,n)

    ! multiplying wrk(:,:,:,n) by filter_g
    ! remember both are in the complex form
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx + 1, 2
             a = fields(ix,iz,iy,n)
             b = fields(ix+1,iz,iy,n)
             c = filter_g(ix,iz,iy)
             d = filter_g(ix+1,iz,iy)

             fields(ix    , iz, iy, n) = a*c - b*d
             fields(ix + 1, iz, iy, n) = b*c + a*d
          end do
       end do
    end do

!!$    ! perform inverse FFT
!!$    call xFFT3d(-1,n)

    return
  end subroutine filter_xfftw_fields

!======================================================================
!======================================================================
!======================================================================
!  subroutine that initializes the filter function filter_g
!  IN: filter_type - type of the filter
!      delta - filter's characteristic width
!
!  OUT: filter_g - normalized FFT of the filter function
!======================================================================
  subroutine filter_xfftw_init

    use m_parameters
    use m_io
    use m_work
    implicit none

    real*8  delta

    integer :: di,dj,dk,i,j,k, m
    integer :: idelta,sum1
    integer :: ii(8)
    real*8  :: rx,ry,rz
    real*8  :: const1,const2

    filter_size = particles_filter_size
    delta = filter_size

    if (task.eq.'hydro') delta = dx*four

    if (delta .lt. dx*three) then
       write(out,*) "filter_xfftw_init: delta is too small:",delta,dx
       write(out,*) "Must be at least 3*dx = ",three*dx
       call flush(out)
       call my_exit(-1)
    end if

    write(out,"('initializing the filter, filter_size = ',2e15.6)") delta, dx
    call flush(out)

    ! the main idea is as follows:
    ! 1) create the filter kernel in one of the fields array
    ! 2) transform it into the Fourier space
    ! 3) put it into the array filter_g 

    ! number of the slice in the fields array that will be used
    m = LBOUND(fields,4)


!---------------------------------------------------------------
!   define the filtering function
!---------------------------------------------------------------
    case_filter_type: select case (filter_type)

!-------------------------------------------------------------------
!  tophat filter
!-------------------------------------------------------------------
    case (1)
       write(out,*) '-- tophat filter, delta =',delta

       write(out,*) "Tophat filter is not working at the moment."
       write(out,*) "We're sorry for the inconvenience, stopping."
       call my_exit(-1)

       fields(:,:,:,m) = zip

       idelta = delta / dx / 2
       ! normalization constant
       const1 = real((2*idelta+1)**3,8)
       const2 = 1.0d0 / const1
       do k = 1,nz
          dk = min(myid*nz+k-1,nz*numprocs-(myid*nz+k)+1)
          do j = 1,ny
             dj = min(j-1,ny-j+1)
             do i = 1,nx
                di = min(i-1,nx-i+1)
                fields(i,j,k,m) = zip
                if (di.le.idelta .and. dj.le.idelta .and. dk.le.idelta) then
                   fields(i,j,k,m) = const2
                end if
             end do
          end do
       end do
!-------------------------------------------------------------------
!  Gaussian filter
!-------------------------------------------------------------------
    case (2)
       write(out,*) '-- Gaussian filter, delta =',delta
       call flush(out)

       fields(:,:,:,m) = zip

       const1 = 6.0d0 / delta**2
       const2 = sqrt(const1/PI)**3 *dx**3

       do k = 1,nz
          dk = min(myid*nz+k-1,nz*numprocs-(myid*nz+k)+1)
          rz = dx * real(dk,8)
          do j = 1,ny
             dj = min(j-1,ny-j+1)
             ry = dx * real(dj,8)
             do i = 1,nx
                di = min(i-1,nx-i+1)
                rx = dx * real(di,8)
                rx = rx*rx+ry*ry+rz*rz
                fields(i,j,k,m) = const2*exp(-const1*rx)
                if (fields(i,j,k,m) < 1.e-20) fields(i,j,k,m) = zip
             end do
          end do
       end do

!-------------------------------------------------------------------
!  linear filter
!-------------------------------------------------------------------
    case(3)
       write(out,*) '-- linear filter, delta =',delta
!       if(myid.eq.0) print*,'-- linear filter, delta =',delta
       const1 = 24.0d0 / (PI*delta**3)  *dx**3
       const2 = delta**2 / 4.0d0
       do k = 1,nz
          dk = min(myid*nz+k-1,nz*numprocs-(myid*nz+k)+1)
          rz = dx * real(dk,8)
          do j = 1,ny
             dj = min(j-1,ny-j+1)
             ry = dx * real(dj,8)
             do i = 1,nx
                di = min(i-1,nx-i+1)
                rx = dx * real(di,8)

                rx = rx*rx+ry*ry+rz*rz
                fields(i,j,k,m) = zip
                if (rx.le.const2) fields(i,j,k,m) = &
                     const1*(1.0d0-2.0d0*sqrt(rx)/delta)
             end do
          end do
       end do

    case default
       print *,'FILTER_FFT_INIT: wrong filter_type:',filter_type
       stop
    end select case_filter_type

    ! outputting the sum of all elements of G.
    ! It should equal 1.0.
    const1 = sum(fields(1:nx,1:ny,1:nz,m))
    call MPI_REDUCE(const1,const2,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)
    if(myid.eq.0) write(out,*) '-- NORM OF G: ',const2

    ! compute FFT of g
    call xFFT3d_fields(1,m)

    ! putting the FFT of the filter into filter_g.
    ! allocating the filter array
    if (.not.allocated(filter_g)) then
       allocate(filter_g(nx+2,ny,nz), stat=ierr)
       if (ierr /= 0) stop 'cannot allocate filter_g'
       filter_g = zip
       write(out,*) 'allocated filter_g'
       call flush(out)
    end if
    filter_g(:,:,:) = fields(:,:,:,m)

!!$!  no need to normalize filter_g because our implementation of FFT
!!$!  normalizes the result anyways
!!$    const1 = one/real(nx**3,8)
!!$    do k=1,nz
!!$       do j=1,ny
!!$          do i=1,nx+2
!!$             filter_g(i,j,k) = wrk(i,j,k,m)*const1
!!$          end do
!!$       end do
!!$    end do


    return
  end subroutine filter_xfftw_init


end module m_filter_xfftw
