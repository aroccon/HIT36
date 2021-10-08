!================================================================================
!================================================================================
!  Initialization of passive scalars
!================================================================================

subroutine init_scalars

  use m_parameters
  use m_io
  use m_fields
  use x_fftw, only : ialias

  implicit none

  integer :: n_scalar, i, j, k

  write(out,*) 'Generating scalars'
  call flush(out)

  do n_scalar = 1,n_scalars
     call init_scalar(n_scalar)
  end do

  write(out,*) "Generated the scalars."
  call flush(out)

  ! now making sure that the scalars do not have any high
  ! Fourier harmonics by zeroing out everything that has a wavenumber
  ! that potentially can produce aliasing
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx+2
           if (ialias(i,j,k).gt.0) fields(i,j,k,4:3+n_scalars) = zip
        end do
     end do
  end do

  return

end subroutine init_scalars

!================================================================================
!================================================================================

subroutine init_scalar(n_scalar)

  use m_parameters
  use m_io
  use m_fields

  implicit none

  integer :: n_scalar, ic_type, sc_type

  write(out,*) 'Generating scalar #', n_scalar
  call flush(out)

  sc_type = scalar_type(n_scalar)
  ic_type = sc_type - (sc_type/100)*100

  if (ic_type.eq.0) then   
     ! gradient source - no need for initial conditions
     ! thus making the initial scalar field zero
     write(out,*) "The scalar type = 0, nothing to generate"
     call flush(out)
     fields(:,:,:,n_scalar+3) = zip

  elseif (ic_type.lt.10) then
     ! if the last two digits of the scalar type are less than 10,
     ! the scalar initial conditions are generated based on the
     ! particular spectrum of the scalar
     call init_scalar_spectrum(n_scalar)

  else
     ! if the last two digits are bigger than 10, the scalar is generated
     ! in physical space and then transformed in the Fourier space
     call init_scalar_space(n_scalar)

  end if
  return
end subroutine init_scalar
!================================================================================
!================================================================================
subroutine init_scalar_spectrum(n_scalar)
!================================================================================

  use m_openmpi
  use m_parameters
  use m_io
  use m_fields
  use m_work
  use x_fftw
  use m_rand_knuth
  use RANDu

  implicit none

  integer :: i, j, k, n, n_scalar
  integer *8 :: i8

  real*8, allocatable :: e_spec(:), e_spec1(:), rr(:)
  integer *8, allocatable :: hits(:), hits1(:)

  integer   :: n_shell
  real*8    :: sc_rad1, sc_rad2

  real*8 :: wmag, wmag2, ratio, fac, fac2

!--------------------------------------------------------------------------------
  write(out,*) " Generating scalar # ",n_scalar
  call flush(out)

  ! Initializing the random sequence with the seed RN2
  fac = random(-RN2)

  ! allocate work arrays
  allocate( e_spec(kmax), e_spec1(kmax), hits(kmax), hits1(kmax), &
       rr(nx+2), stat=ierr)

  ! bringing the processors to their own places in the random sequence
  ! ("2" is there because we're generating two random number fields
  ! for each scalar field
  ! using i8 because it's int*8
  do i8 = 1,myid*(nx+2)*ny*nz*2
     fac = random(RN2)
  end do

  ! now filling the arrays wrk1, wrk2
  do n = 1,2
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+2
              wrk(i,j,k,n) = random(RN2)
           end do
        end do
     end do
  end do

  ! bringing the random numbers to the same 
  ! point in the sequence again
  do i8 = 1,(numprocs-myid-1)*(nx+2)*ny*nz*2
     fac = random(RN2)
  end do

  ! making  random array with Gaussian PDF 
  ! out of the two arrays that we generated
  wrk(:,:,:,3) = sqrt(-two*log(wrk(:,:,:,1))) * sin(TWO_PI*wrk(:,:,:,2))

  ! go to Fourier space
  call xFFT3d(1,3)

!-------------------------------------------------------------------------------
!     Calculating the scalar spectrum
!-------------------------------------------------------------------------------

  ! need this normalization factor because the FFT is unnormalized
  fac = one / real(nx*ny*nz_all)**2

  e_spec1 = zip
  e_spec = zip
  hits = 0
  hits1 = 0

  ! assembling the scalar energy in each shell and number of hits in each shell
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx

           n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
           if (n_shell .gt. 0 .and. n_shell .le. kmax) then
              fac2 = fac * wrk(i,j,k,3)**2
              if (akx(i).eq.0.d0) fac2 = fac2 * 0.5d0
              e_spec1(n_shell) = e_spec1(n_shell) + fac2
           end if
        end do
     end do
  end do

  ! reducing the number of hits and energy to two arrays on master node
  count = kmax
  call MPI_REDUCE(e_spec1,e_spec,count,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)
  ! broadcasting the spectrum
  count = kmax
  call MPI_BCAST(e_spec,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)

!-------------------------------------------------------------------------------
!  Now make the spectrum to be as desired
!-------------------------------------------------------------------------------
  ! first, define the desired spectrum
  do k = 1,kmax

     wmag = real(k, 8)
     ratio = wmag / peak_wavenum_sc(n_scalar)

     if (scalar_type(n_scalar).eq.0) then
        ! Plain Kolmogorov spectrum
        e_spec1(k) = wmag**(-5.d0/3.d0)

     else if (scalar_type(n_scalar).eq.1 .or. scalar_type(n_scalar).eq.3) then
        ! Exponential spectrum
        e_spec1(k) =  ratio**3 / peak_wavenum_sc(n_scalar) * exp(-3.0D0*ratio)

     else if (scalar_type(n_scalar).eq.2) then
        ! Von Karman spectrum
        fac = two * PI * ratio
        e_spec1(k) = fac**4 / (one + fac**2)**3

     else
        write(out,*) "INIT_SCALARS: WRONG INITIAL SPECTRUM TYPE: ",scalar_type(n_scalar)
        call flush(out)
        stop

     end if
  end do

  ! normalize it so it has the unit total energy
  e_spec1 = e_spec1 / sum(e_spec1(1:kmax))

  ! now go over all Fourier shells and multiply the velocities in a shell by
  ! the sqrt of ratio of the resired to the current spectrum
  fields(:,:,:,3+n_scalar) = zip

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+2

           n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
           if (n_shell .gt. 0 .and. n_shell .le. kmax .and. e_spec(n_shell) .gt. zip) then
              fields(i,j,k,3+n_scalar) = wrk(i,j,k,3) * sqrt(e_spec1(n_shell)/e_spec(n_shell))
           else
              fields(i,j,k,3+n_scalar) = zip
           end if

        end do
     end do
  end do

!-------------------------------------------------------------------------------
!   Creating scalars with double-delta PDF
!-------------------------------------------------------------------------------

  if (scalar_type(n_scalar).eq.3) then
     wrk(:,:,:,0) = fields(:,:,:,3+n_scalar)
     call xFFT3d(-1,0)

     ! making it double-delta (0.9 and -0.9)
     wrk(:,:,:,0) = sign(one,wrk(:,:,:,0)) * 0.9d0
     call xFFT3d(1,0)

     ! smoothing it by zeroing out high harmonics
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+2
              n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
              if (n_shell .eq. 0 .or. n_shell .ge. kmax*2/3) then
                 wrk(i,j,k,0) = zip
              end if
           end do
        end do
     end do
     fields(:,:,:,3+n_scalar) = wrk(:,:,:,0)

  end if

  ! deallocate work arrays
  deallocate(e_spec, e_spec1, rr, hits, hits1, stat=ierr)

  return
end subroutine init_scalar_spectrum

!================================================================================
!================================================================================
!================================================================================
!================================================================================
subroutine init_scalar_space(n_scalar)
!================================================================================

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i, k, n_scalar, sc_type, ic_type, nfi
  real*8  :: zloc, sctmp, h, xx

  nfi = 3 + n_scalar

  write(out,*) " Generating scalar # ", n_scalar
  call flush(out)

  sc_type = scalar_type(n_scalar)
  ic_type = sc_type - (sc_type/100)*100

  select case (ic_type)

!---------------------------------------------------
! single slab of the scalar
!---------------------------------------------------
  case(11) 

     ! how much to smear out the interface
     ! the minimum interface length is set to 2*dz
     ! the maximum is set to 2*PI/(peak wavenumber for the scalar, taken from
     ! the .in file)
     h = max(2.*dz, two*PI/peak_wavenum_sc(n_scalar))

     ! creating array of scalar
     do k = 1,nz
        zloc = dble(myid*nz + k-1) * dz
        sctmp = tanh((zloc-PI*0.5)/h) - tanh((zloc-PI*1.5)/h) - one
        wrk(:,:,k,0) = sctmp
     end do

     ! FFT of the scalar
     call xFFT3d(1,0)

     ! putting it into the scalar array
     fields(:,:,:, 3+n_scalar) = wrk(:,:,:,0)
     if (iammaster) fields(1,1,1,3+n_scalar) = zip

!---------------------------------------------------
! two slabs of the scalar
!---------------------------------------------------
  case(12)

     write(out,*) "-- Double-slab scalar in real space"
     call flush(out)

     ! how much to smear out the interface
     h = max(8.*dz, PI/8.d0)

     ! creating array of scalar
     do i = 1,nx
        xx = dble(i-1) * dx
        sctmp = tanh((xx-PI*0.25)/h) - tanh((xx-PI*0.75)/h) + tanh((xx-PI*1.25)/h) - tanh((xx-PI*1.75)/h)
        wrk(i,:,:,0) = sctmp - one
        ! now it is between -1 and 1
     end do

     ! FFT of the scalar
     call xFFT3d(1,0)

     ! putting it into the scalar array
     fields(:,:,:, 3+n_scalar) = wrk(:,:,:,0)

     ! making sure that the mean is ero
     if (iammaster) fields(1,1,1,3+n_scalar) = zip

!---------------------------------------------------
!  N slabs of the scalar
!  actually more like N sinusoidal waves
!---------------------------------------------------
  case(13)
     write(out,*) "-- Multi-slab scalar in real space"
     call flush(out)

     ! making sure that we can support the desired numberof waves with our FFT
     peak_wavenum_sc(n_scalar) = min(peak_wavenum_sc(n_scalar),real(nx/8)) 

     ! definition of initial scalar
     fields(:,:,:,nfi) = zip
     do i = 1, nx
        fields(i,:,:,nfi) = sin(peak_wavenum_sc(n_scalar) * dx * real(i-1))
     end do
     call xFFT3d_fields(1,nfi)


  case default
     write(out,*) "INIT_SCALARS: UNEXPECTED SCALAR TYPE: ", scalar_type(n_scalar)
     call flush(out)
     stop
  end select

  write(out,*) "Initialized the scalars."
  call flush(out)

  return

end subroutine init_scalar_space
