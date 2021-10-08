subroutine pressure

  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i, j, k, n
  real*8  :: div1, div2, lapl1, lapl2, p1, p2

  do k = 1,nz
    do j = 1,ny
      do i = 1,nx+1,2

        ! getting the divergence (i*k*\hat(u))
        ! remember that in the Fourier space indicies go as (ix,iz,iy)

        ! getting divergence
        div2 =     akx(i  )*fields(i  ,j,k,1) + aky(k)*fields(i  ,j,k,2) + akz(j)*fields(i  ,j,k,3)
        div1 = - ( akx(i+1)*fields(i+1,j,k,1) + aky(k)*fields(i+1,j,k,2) + akz(j)*fields(i+1,j,k,3) )

        ! inverce laplace operator
        lapl1 =  akx(i  )**2 + aky(k)**2 + akz(j)**2
        lapl2 =  akx(i+1)**2 + aky(k)**2 + akz(j)**2

        if (lapl1.eq.0.d0) lapl1 = 9e20
        if (lapl2.eq.0.d0) lapl2 = 9e20

        ! calculating pressure
        p1 = - div1 / lapl1
        p2 = - div2 / lapl2

        ! Taking derivatives of the pressure and subtracting from the corresponding velocities
        fields(i  ,j,k,1) = fields(i  ,j,k,1) + p2 * akx(i+1)
        fields(i+1,j,k,1) = fields(i+1,j,k,1) - p1 * akx(i  )

        fields(i  ,j,k,2) = fields(i  ,j,k,2) + p2 * aky(k)
        fields(i+1,j,k,2) = fields(i+1,j,k,2) - p1 * aky(k)

        fields(i  ,j,k,3) = fields(i  ,j,k,3) + p2 * akz(j)
        fields(i+1,j,k,3) = fields(i+1,j,k,3) - p1 * akz(j)

      end do
    end do
  end do


  return
end subroutine pressure



subroutine divergence

  use m_openmpi
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i,j,k
  real*8 :: dmin, dmax, d1


  wrk(:,:,:,1:3) = fields(:,:,:,1:3)

  call x_derivative(1,'x',4)
  call x_derivative(2,'y',5)
  call x_derivative(3,'z',6)

  call xFFT3d(-1,4)
  call xFFT3d(-1,5)
  call xFFT3d(-1,6)

  wrk(:,:,:,0) = wrk(:,:,:,4) + wrk(:,:,:,5) + wrk(:,:,:,6)

  d1 = minval(wrk(1:nx,:,:,0))
  call MPI_REDUCE(d1,dmin,1,MPI_REAL8,MPI_MIN,0,MPI_COMM_TASK,mpi_err)
  d1 = maxval(wrk(1:nx,:,:,0))
  call MPI_REDUCE(d1,dmax,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_TASK,mpi_err)

  if (myid.eq.0) then
    write(out,*) 'divergence:',dmin,dmax
!!$    print *, 'divergence:',dmin,dmax
    call flush(out)
  end if

  return
end subroutine divergence



