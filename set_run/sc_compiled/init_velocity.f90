subroutine init_velocity

  use m_openmpi
  use m_parameters
  use m_io
  use m_fields
  use m_work
  use x_fftw
  use m_rand_knuth
  use RANDu

  implicit none

  integer :: i, j, k, n
  integer*8 :: seed1, seed2, i8, j8, k8

  integer :: time_array(8)


  real,       allocatable   :: rr(:)
  real*8,     allocatable :: e_spec(:), e_spec1(:)
  integer *8, allocatable :: hits(:), hits1(:)

  integer   :: n_shell
  real*8    :: sc_rad1, sc_rad2

  real*8 :: wmag, wmag2, ratio, fac, fac2


!--------------------------------------------------------------------------------
!  First, if it's a Taylor-Green vortex, then initialize and quit
!--------------------------------------------------------------------------------
  if (isp_type .eq. -1) then
     call init_velocity_taylor_green
     return
  end if


!================================================================================
  allocate( e_spec(kmax), e_spec1(kmax), rr(nx+2), hits(kmax), hits1(kmax), stat=ierr)
  if (ierr.ne.0) stop "cannot allocate the init_velocity arrays"



  write(out,*) 'generating random velocities'
  call flush(out)


!-------------------------------------------------------------------------------
!  Generate the velocities
!-------------------------------------------------------------------------------

  ! initialize the random number sequence by the seed from the first processor
  if (myid.eq.0) call system_clock(seed1,seed2)
  if (myid.eq.0) then
     call date_and_time(values=time_array)
     seed1 = time_array(8)
  end if
  count = 1
  call MPI_BCAST(seed1,count,MPI_INTEGER8,0,MPI_COMM_TASK,mpi_err)

!!$  seed1 = 23498675
!!$  call rand_knuth_init(seed1)
!!$  write(out,*) 'seed1 = ',seed1
!!$  call flush(out)

  seed1 = RN1
  write(out,*) "RANDOM SEED FOR VELOCITIES = ", seed1
  call flush(out)
  rseed = real(seed1,8)
  fac = random(-rseed)




  ! bringing the processors to their own places in the random sequence
  ! ("6" is there because we're generating six fields
  ! using seed1 because it's int*8

!!$  write(out,*) "Will scroll down to my initial position", myid*ny*nz*6
!!$  call flush(out)

!!$  do i = 1,myid*ny*nz*6
!!$     call rand_knuth(rr,nx+2)
!!$  end do

  do i8 = 1,myid*(nx+2)*ny*nz*6
     fac = random(rseed)
  end do

!!$  write(out,*) "Scrolled."
!!$  call flush(out)

  ! now filling the arrays wrk1...wrk6 
  do n = 1,6
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+2
              wrk(i,j,k,n) = random(rseed)
           end do
        end do
     end do
  end do

  ! just in case, bringing the random numbers to the same 
  ! point in the sequence again
  do seed1 = 1,int((numprocs-myid-1)*(nx+2)*ny*nz*6,8)
     fac = random(rseed)
  end do

  ! making three random arrays with Gaussian PDF 
  ! out of the six arrays that we generated
  wrk(:,:,:,1:3) = sqrt(-two*log(wrk(:,:,:,1:3))) * sin(TWO_PI*wrk(:,:,:,4:6))

  ! --- Making three arrays that have Gaussian PDF and the incompressibility property

  ! go to Fourier space
  do n = 1,3
     call xFFT3d(1,n)
  end do

  ! assemble the arrays in wrk4..6, only the wavenumbers below kmax
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+1,2

           n_shell = nint(sqrt(akx(i)**2 + aky(k)**2 + akz(j)**2))


           if (n_shell .gt. 0 .and. n_shell .le. kmax) then

              wrk(i  ,j,k,4) = - (aky(k)*wrk(i+1,j,k,3) - akz(j)*wrk(i+1,j,k,2))
              wrk(i+1,j,k,4) =    aky(k)*wrk(i  ,j,k,3) - akz(j)*wrk(i  ,j,k,2)

              wrk(i  ,j,k,5) = - (akz(j)*wrk(i+1,j,k,1) - akx(i+1)*wrk(i+1,j,k,3))
              wrk(i+1,j,k,5) =    akz(j)*wrk(i  ,j,k,1) - akx(i  )*wrk(i  ,j,k,3)

              wrk(i  ,j,k,6) = - (akx(i+1)*wrk(i+1,j,k,2) - aky(k)*wrk(i+1,j,k,1))
              wrk(i+1,j,k,6) =    akx(i  )*wrk(i  ,j,k,2) - aky(k)*wrk(i  ,j,k,1)

           else
              wrk(i:i+1,j,k,4:6) = zip
           end if

        end do
     end do
  end do

  fields(:,:,:,1:3) = wrk(:,:,:,4:6)

!-------------------------------------------------------------------------------
!     Making the spectrum to be what it should 
!-------------------------------------------------------------------------------


  ! --- first get the energy spectrum (copied from m_stat.f90)

  ! need this normalization factor because the FFT is unnormalized
  fac = one / real(nx*ny*nz_all)**2

  e_spec1 = zip
  e_spec = zip
  hits = 0
  hits1 = 0

  ! assembling the total energy in each shell and number of hits in each shell
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx

           n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
           if (n_shell .gt. 0 .and. n_shell .le. kmax) then
              fac2 = fac * (fields(i,j,k,1)**2 + fields(i,j,k,2)**2 + fields(i,j,k,3)**2)
              if (akx(i).eq.0.d0) fac2 = fac2 * 0.5d0
              e_spec1(n_shell) = e_spec1(n_shell) + fac2
           end if

        end do
     end do
  end do

  ! reducing the number of hits and energy to two arrays on master node
  count = kmax
  call MPI_REDUCE(e_spec1,e_spec,count,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)
  count = kmax
  call MPI_BCAST(e_spec,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)

!-------------------------------------------------------------------------------
!  Now make the spectrum to be as desired
!-------------------------------------------------------------------------------

  ! first, define the desired spectrum
  do k = 1,kmax

     wmag = real(k, 8)
     ratio = wmag / peak_wavenum

     if (isp_type.eq.0) then
        ! Plain Kolmogorov spectrum
        e_spec1(k) = wmag**(-5.d0/3.d0)

     else if (isp_type.eq.1) then
        ! Exponential spectrum
        e_spec1(k) =  ratio**3 / peak_wavenum * exp(-3.0D0*ratio)

     else if (isp_type.eq.3) then
        ! Von Karman spectrum
        fac = two * PI * ratio
        e_spec1(k) = fac**4 / (one + fac**2)**3

     else
        write(out,*) "ERROR: WRONG INITIAL SPECTRUM TYPE: ",isp_type
        call flush(out)
        stop

     end if
  end do

!  normalize it so it has the unit total energy
  e_spec1 = e_spec1 / sum(e_spec1(1:kmax))

  ! now go over all Fourier shells and multiply the velocities in a shell by
  ! the sqrt of ratio of the resired to the current spectrum
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+2

           n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
           if (n_shell .gt. 0 .and. n_shell .le. kmax .and. e_spec(n_shell) .gt. zip) then
              fields(i,j,k,1:3) = fields(i,j,k,1:3) * sqrt(e_spec1(n_shell)/e_spec(n_shell))
           else
              fields(i,j,k,1:3) = zip
           end if

        end do
     end do
  end do

  write(out,*) "Generated the velocities."
  call flush(out)


  ! deallocate work arrays
  deallocate(e_spec, e_spec1, rr, hits, hits1, stat=ierr)
  return
end subroutine init_velocity


!================================================================================
!  Initialize the velocities with Taylor-Green vortex
!================================================================================
subroutine init_velocity_taylor_green

  use m_openmpi
  use m_parameters
  use m_io
  use m_fields
  use m_work
  use x_fftw

  implicit none

  logical :: verbose = .true.

  integer :: i, j, k, n
  real*8 :: xx, yy, zz

  if (verbose) write(out,*) " --- Initial velocity field is Taylor-Green vortex"
  if (verbose) call flush(out)

  do k = 1, nz
     zz = real(nz*myid + k - 1, 8) * dz
     do j = 1, ny
        yy = real(j-1, 8) * dy
        do i = 1, nx
           xx = real(i-1,8) *dx
           wrk(i,j,k,1) =   sin(xx) * cos(yy) * cos(zz)
           wrk(i,j,k,2) = - cos(xx) * sin(yy) * cos(zz)
           wrk(i,j,k,3) = zip
        end do
     end do
  end do

  call xFFT3D(1, 1)
  call xFFT3D(1, 2)
  call xFFT3D(1, 3)
  fields(:,:,:,1:3) = wrk(:,:,:,1:3)

  if (verbose) write(out,*) " --- initialized."
  if (verbose) call flush(out)

  return
end subroutine init_velocity_taylor_green


