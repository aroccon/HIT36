subroutine rhs_velocity
use m_openmpi
use m_parameters
use m_fields
use m_work
use x_fftw
implicit none
integer :: i, j, k, n, nx3
real*8  :: t1(0:6), rtmp, wnum2, rnx3




! getting the Courant number (on the master process only)
wrk(:,:,:,4) = abs(wrk(:,:,:,1)) + abs(wrk(:,:,:,2)) + abs(wrk(:,:,:,3))
rtmp = maxval(wrk(1:nx,:,:,4))
call MPI_REDUCE(rtmp,courant,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_TASK,mpi_err)
!if (variable_dt) then
!  count = 1
!  call MPI_BCAST(courant,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
!end if
courant = courant * dt / dx



!  Calculating the right-hand side for the velocities
!  Two options available: the standard 2/3 rule (dealias=0) or phase shift and truncation (dealias=1).
two_thirds_rule: if (dealias.eq.0) then
! getting all 6 products of velocities
do k = 1,nz
  do j = 1,ny
    do i = 1,nx
      t1(1) = wrk(i,j,k,1) * wrk(i,j,k,1) !u*u
      t1(2) = wrk(i,j,k,1) * wrk(i,j,k,2) !u*v
      t1(3) = wrk(i,j,k,1) * wrk(i,j,k,3) !u*w
      t1(4) = wrk(i,j,k,2) * wrk(i,j,k,2) !v*v
      t1(5) = wrk(i,j,k,2) * wrk(i,j,k,3) !v*W
      t1(6) = wrk(i,j,k,3) * wrk(i,j,k,3) !W*w
      do n = 1,6
        wrk(i,j,k,n) = t1(n)
      end do
    end do
  end do
end do


! converting the products to the Fourier space
do n = 1,6
  call xFFT3d(1,n)
end do

! Building the RHS.
! 1) put into wrk arrays the convectove terms (that will be multiplyed by "i"
! 2) and the factor that corresponds to the diffusion
! Do not forget that in Fourier space the indicies are (ix, iz, iy)
do k = 1,nz
  do j = 1,ny
    do i = 1,nx+2
      t1(1) = - ( akx(i) * wrk(i,j,k,1) + aky(k) * wrk(i,j,k,2) + akz(j) * wrk(i,j,k,3))
      t1(2) = - ( akx(i) * wrk(i,j,k,2) + aky(k) * wrk(i,j,k,4) + akz(j) * wrk(i,j,k,5))
      t1(3) = - ( akx(i) * wrk(i,j,k,3) + aky(k) * wrk(i,j,k,5) + akz(j) * wrk(i,j,k,6))
      t1(4) = - nu * ( akx(i)**2 + aky(k)**2 + akz(j)**2 )
      do n = 1,4
        wrk(i,j,k,n) = t1(n)
      end do
    end do
  end do
end do

! now take the actual fields from fields(:,:,:,:) and calculate the RHSs
! at this moment the contains of wrk(:,:,:,1:3) are the convective terms in the RHS
! which are not yet multiplied by "i"
! wrk(:,:,:,4) contains the Laplace operator in Fourier space.  To get the diffusion term
 ! we need to take wrk(:,:,:,4) and multiply it by the velocity

t1(6) = real(kmax,8)
do k = 1,nz
  do j = 1,ny
    do i = 1,nx+1,2
      ! If the dealiasing option is 2/3-rule (dealias=0) then we retain the modes
      ! inside the cube described by $| k_i | \leq  k_{max}$, $i=1,2,3$.
      ! The rest of the modes is purged
      if (ialias(i,j,k) .gt. 0) then
        ! setting the Fourier components to zero
        wrk(i  ,j,k,1:3) = zip
        wrk(i+1,j,k,1:3) = zip
      else
        ! RHS for u, v and w
        do n = 1,3
          ! taking the convective term, multiply it by "i"
          ! (see how it's done in x_fftw.f90)
          ! and adding the diffusion term
          rtmp =           - wrk(i+1,j,k,n) + wrk(i  ,j,k,4) * fields(i  ,j,k,n)
          wrk(i+1,j,k,n) =   wrk(i  ,j,k,n) + wrk(i+1,j,k,4) * fields(i+1,j,k,n)
          wrk(i  ,j,k,n) = rtmp
        end do
      end if
    end do
  end do
end do

end if two_thirds_rule

!--------------------------------------------------------------------------------
!  The second option (dealias=1).  All pairwise products of velocities are
!  dealiased using one phase shift of (dx/2,dy/2,dz/2).
!--------------------------------------------------------------------------------
  phase_shifting: if (dealias.eq.1) then

     ! work parameters
     wrk(:,:,:,0) = zip

     ! getting all 6 products of velocities
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx

              t1(1) = wrk(i,j,k,1) * wrk(i,j,k,1)
              t1(2) = wrk(i,j,k,1) * wrk(i,j,k,2)
              t1(3) = wrk(i,j,k,1) * wrk(i,j,k,3)
              t1(4) = wrk(i,j,k,2) * wrk(i,j,k,2)
              t1(5) = wrk(i,j,k,2) * wrk(i,j,k,3)
              t1(6) = wrk(i,j,k,3) * wrk(i,j,k,3)
              do n = 1,6
                 wrk(i,j,k,n) = t1(n)
              end do
           end do
        end do
     end do

     ! converting the products to the Fourier space
     do n = 1,6
        call xFFT3d(1,n)
     end do

     ! Building the RHS.
     ! First, put into wrk arrays the convectove terms (that will be multiplyed by "i"
     ! later) and the factor that corresponds to the diffusion

     ! Do not forget that in Fourier space the indicies are (ix, iz, iy)
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+2

              t1(1) = - ( akx(i) * wrk(i,j,k,1) + aky(k) * wrk(i,j,k,2) + akz(j) * wrk(i,j,k,3) )
              t1(2) = - ( akx(i) * wrk(i,j,k,2) + aky(k) * wrk(i,j,k,4) + akz(j) * wrk(i,j,k,5) )
              t1(3) = - ( akx(i) * wrk(i,j,k,3) + aky(k) * wrk(i,j,k,5) + akz(j) * wrk(i,j,k,6) )
              ! putting a factor from the diffusion term into t1(4) (and later in wrk4)
              t1(4) = - nu * ( akx(i)**2 + aky(k)**2 + akz(j)**2 )
              do n = 1,4
                 wrk(i,j,k,n) = t1(n)
              end do
           end do
        end do
     end do

     ! now use the actual fields from fields(:,:,:,:) to calculate the RHSs

     ! at this moment the contains of wrk(:,:,:,1:3) are the convective terms in the RHS
     ! which are not yet multiplied by "i"

     ! wrk(:,:,:,4) contains the Laplace operator in Fourier space.  To get the diffusion term
     ! we need to take wrk(:,:,:,4) and multiply it by the velocity

     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2

              ! If the dealiasing option is (dealias=1) then we retain the modes
              ! for which no more than one component of the k-vector is larger than nx/3.
              ! The rest of the modes is purged.

              if  (ialias(i,j,k) .gt. 1) then
                 ! setting the Fourier components to zero
                 wrk(i  ,j,k,1:3) = zip
                 wrk(i+1,j,k,1:3) = zip
              else
                 ! RHS for u, v and w
                 do n = 1,3
                    ! taking the HALF of the convective term, multiply it by "i"
                    ! and adding the diffusion term
                    rtmp =           - 0.5d0 * wrk(i+1,j,k,n) + wrk(i  ,j,k,4) * fields(i  ,j,k,n)
                    wrk(i+1,j,k,n) =   0.5d0 * wrk(i  ,j,k,n) + wrk(i+1,j,k,4) * fields(i+1,j,k,n)
                    wrk(i  ,j,k,n) = rtmp
                 end do
              end if

           end do
        end do
     end do


!--------------------------------------------------------------------------------
!  Second part of the phase shifting technique
!--------------------------------------------------------------------------------

     ! since wrk1...3 are taken by parts of RHS constructed earlier, we can use
     ! only wrk0 and wrk4...6.

     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2

              ! computing sines and cosines for the phase shift of dx/2,dy/2,dz/2
              ! and putting them into wrk0
              wrk(i  ,j,k,0) = cos(half*(akx(i  )+aky(k)+akz(j))*dx)
              wrk(i+1,j,k,0) = sin(half*(akx(i+1)+aky(k)+akz(j))*dx)

              ! wrk4 will have phase-shifted u
              wrk(i  ,j,k,4) = fields(i  ,j,k,1) * wrk(i,j,k,0) - fields(i+1,j,k,1) * wrk(i+1,j,k,0)
              wrk(i+1,j,k,4) = fields(i+1,j,k,1) * wrk(i,j,k,0) + fields(i  ,j,k,1) * wrk(i+1,j,k,0)

              ! wrk5 will have phase-shifted v
              wrk(i  ,j,k,5) = fields(i  ,j,k,2) * wrk(i,j,k,0) - fields(i+1,j,k,2) * wrk(i+1,j,k,0)
              wrk(i+1,j,k,5) = fields(i+1,j,k,2) * wrk(i,j,k,0) + fields(i  ,j,k,2) * wrk(i+1,j,k,0)

           end do
        end do
     end do

     ! transforming u+ and v+ into X-space
     call xFFT3d(-1,4)
     call xFFT3d(-1,5)

     ! now wrk4 and wrk5 contain u+ and v+

     ! getting (u+)*(u+) in real space, converting it to Fourier space,
     ! phase shifting back and adding -0.5*(the results)  to the RHS for u
     wrk(:,:,:,6) = wrk(:,:,:,4)**2
     call xFFT3d(1,6)
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              rtmp           = wrk(i  ,j,k,6) * wrk(i,j,k,0) + wrk(i+1,j,k,6) * wrk(i+1,j,k,0)
              wrk(i+1,j,k,6) = wrk(i+1,j,k,6) * wrk(i,j,k,0) - wrk(i  ,j,k,6) * wrk(i+1,j,k,0)
              wrk(i  ,j,k,6) = rtmp
           end do
        end do
     end do
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              wrk(i  ,j,k,1) = wrk(i  ,j,k,1) + 0.5d0 * akx(i+1) * wrk(i+1,j,k,6)
              wrk(i+1,j,k,1) = wrk(i+1,j,k,1) - 0.5d0 * akx(i  ) * wrk(i  ,j,k,6)
           end do
        end do
     end do

     ! getting (u+)*(v+) in real space, converting it to Fourier space,
     ! phase shifting back and adding -0.5*(the results)  to the RHSs for u and v
     wrk(:,:,:,6) = wrk(:,:,:,4)*wrk(:,:,:,5)
     call xFFT3d(1,6)
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              rtmp           = wrk(i  ,j,k,6) * wrk(i,j,k,0) + wrk(i+1,j,k,6) * wrk(i+1,j,k,0)
              wrk(i+1,j,k,6) = wrk(i+1,j,k,6) * wrk(i,j,k,0) - wrk(i  ,j,k,6) * wrk(i+1,j,k,0)
              wrk(i  ,j,k,6) = rtmp
           end do
        end do
     end do
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              wrk(i  ,j,k,1) = wrk(i  ,j,k,1) + 0.5d0 * aky(k) * wrk(i+1,j,k,6)
              wrk(i+1,j,k,1) = wrk(i+1,j,k,1) - 0.5d0 * aky(k) * wrk(i  ,j,k,6)

              wrk(i  ,j,k,2) = wrk(i  ,j,k,2) + 0.5d0 * akx(i+1) * wrk(i+1,j,k,6)
              wrk(i+1,j,k,2) = wrk(i+1,j,k,2) - 0.5d0 * akx(i  ) * wrk(i  ,j,k,6)
           end do
        end do
     end do

     ! getting (v+)*(v+) in real space, converting it to Fourier space,
     ! phase shifting back and adding -0.5*(the results)  to the RHS for v
     wrk(:,:,:,6) = wrk(:,:,:,5)**2
     call xFFT3d(1,6)
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              rtmp           = wrk(i  ,j,k,6) * wrk(i,j,k,0) + wrk(i+1,j,k,6) * wrk(i+1,j,k,0)
              wrk(i+1,j,k,6) = wrk(i+1,j,k,6) * wrk(i,j,k,0) - wrk(i  ,j,k,6) * wrk(i+1,j,k,0)
              wrk(i  ,j,k,6) = rtmp
           end do
        end do
     end do
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              wrk(i  ,j,k,2) = wrk(i  ,j,k,2) + 0.5d0 * aky(k) * wrk(i+1,j,k,6)
              wrk(i+1,j,k,2) = wrk(i+1,j,k,2) - 0.5d0 * aky(k) * wrk(i  ,j,k,6)
           end do
        end do
     end do

     ! now get the (w+) in wrk6
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              ! wrk6 will have phase-shifted w
              wrk(i  ,j,k,6) = fields(i  ,j,k,3) * wrk(i,j,k,0) - fields(i+1,j,k,3) * wrk(i+1,j,k,0)
              wrk(i+1,j,k,6) = fields(i+1,j,k,3) * wrk(i,j,k,0) + fields(i  ,j,k,3) * wrk(i+1,j,k,0)
           end do
        end do
     end do
     ! transforming w+ into X-space
     call xFFT3d(-1,6)

     ! at this point wrk4..6 contain (u+), (v+) and (w+) in real space.
     ! the combinations that we have not dealt with are: uw, vw and ww.
     ! we'll deal with all three of them at once.

     ! first get all three of these in wrk4...6 and
     wrk(:,:,:,4) = wrk(:,:,:,4) * wrk(:,:,:,6)
     wrk(:,:,:,5) = wrk(:,:,:,5) * wrk(:,:,:,6)
     wrk(:,:,:,6) = wrk(:,:,:,6)**2

     ! transform them into Fourier space
     call xFFT3d(1,4)
     call xFFT3d(1,5)
     call xFFT3d(1,6)

     ! phase shift back to origianl grid and add to corresponding RHSs
     do n = 4,6
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2
                 rtmp           = wrk(i  ,j,k,n) * wrk(i,j,k,0) + wrk(i+1,j,k,n) * wrk(i+1,j,k,0)
                 wrk(i+1,j,k,n) = wrk(i+1,j,k,n) * wrk(i,j,k,0) - wrk(i  ,j,k,n) * wrk(i+1,j,k,0)
                 wrk(i  ,j,k,n) = rtmp
              end do
           end do
        end do
     end do

     ! adding to corresponding RHSs
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2

              ! If the dealiasing option is (dealias=1) then we retain the modes
              ! for which no more than one component of the k-vector is larger than nx/3.
              ! The rest of the modes is purged.

              if  (ialias(i,j,k) .lt. 2) then

                 wrk(i  ,j,k,1) = wrk(i  ,j,k,1) + 0.5d0 * akz(j) * wrk(i+1,j,k,4)
                 wrk(i+1,j,k,1) = wrk(i+1,j,k,1) - 0.5d0 * akz(j) * wrk(i  ,j,k,4)

                 wrk(i  ,j,k,2) = wrk(i  ,j,k,2) + 0.5d0 * akz(j) * wrk(i+1,j,k,5)
                 wrk(i+1,j,k,2) = wrk(i+1,j,k,2) - 0.5d0 * akz(j) * wrk(i  ,j,k,5)

                 wrk(i  ,j,k,3) = wrk(i  ,j,k,3) + 0.5d0 * &
                      (akx(i+1)*wrk(i+1,j,k,4) + aky(k)*wrk(i+1,j,k,5) + akz(j)*wrk(i+1,j,k,6))
                 wrk(i+1,j,k,3) = wrk(i+1,j,k,3) - 0.5d0 * &
                      (akx(i  )*wrk(i  ,j,k,4) + aky(k)*wrk(i  ,j,k,5) + akz(j)*wrk(i  ,j,k,6))

              else
                 wrk(i:i+1,j,k,1) = zip
                 wrk(i:i+1,j,k,2) = zip
                 wrk(i:i+1,j,k,3) = zip
              end if

           end do
        end do
     end do

  end if phase_shifting

  ! if performing large eddy simulations, call LES subroutine to augment
  ! the right hand side for velocioties
!  les_active: if (les) then
!     call les_rhs_velocity
!  end if les_active

return
end subroutine rhs_velocity
