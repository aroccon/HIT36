!================================================================================
!  Subroutine that rescales the current velocities.  The rescales velocity
!  will have the spectrum that is defined in the input file via parameters
!  isp_type (spectrum type) and peak_wavenum (peak wavenumber).
!
!  The program is copies from a part of init_velocity.f90.  Some time in the
!  future we should make it one routine that is called in init_velocity.
!   
!  Time-stamp: <2010-01-25 17:01:55 (chumakov)>
!================================================================================

subroutine velocity_rescale

  use m_openmpi
  use m_parameters
  use m_io
  use m_fields
  use x_fftw

  implicit none

  integer :: i, j, k, n

  real,       allocatable   :: rr(:)
  real*8,     allocatable :: e_spec(:), e_spec1(:)
  integer *8, allocatable :: hits(:), hits1(:)

  integer   :: n_shell
  real*8    :: sc_rad1, sc_rad2

  real*8 :: wmag, wmag2, ratio, fac


  ! if Taylor-Green, return
  if (isp_type.eq.-1) return


!================================================================================
  allocate( e_spec(kmax), e_spec1(kmax), rr(nx+2), hits(kmax), hits1(kmax), stat=ierr)
  if (ierr.ne.0) stop "cannot allocate the init_velocity arrays"

  write(out,*) 'Rescaling the velocities'
  call flush(out)

!-------------------------------------------------------------------------------
!     Making the spectrum to be what is prescribed in the input file <...>.in 
!-------------------------------------------------------------------------------

  ! --- first get the energy spectrum (copied from m_stat.f90)

  ! need this normalization factor because the FFT is unnormalized
  fac = one / real(nx*ny*nz_all)**2

  ! zeroing out the arrays
  e_spec1 = zip
  e_spec = zip
  hits = 0
  hits1 = 0

  ! finding the total energy in each shell and number of hits in each shell
  ! storing them in arrays hits1 and e_spec1
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx

           n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
           if (n_shell .gt. 0 .and. n_shell .le. kmax) then
              hits1(n_shell) = hits1(n_shell) + 1
              e_spec1(n_shell) = e_spec1(n_shell) + &
                   fac * (fields(i,j,k,1)**2 + fields(i,j,k,2)**2 + fields(i,j,k,3)**2)
           end if
        end do
     end do
  end do

  ! reducing the number of hits and energy to two arrays on master node
  count = kmax
  call MPI_REDUCE(hits1,hits,count,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)
  count = kmax
  call MPI_REDUCE(e_spec1,e_spec,count,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)

  ! now the master node counts the energy density in each shell
  if (myid.eq.0) then
     fac = four/three * PI / two
     do k = 1,kmax
        sc_rad1 = real(k,8) + half
        sc_rad2 = real(k,8) - half
        if (k.eq.1) sc_rad2 = 0.d0
        if (hits(k).gt.0) then
           e_spec(k) = e_spec(k) / hits(k) * fac * (sc_rad1**3 - sc_rad2**3)
        else
           e_spec(k) = zip
        end if
     end do
  end if

  ! broadcasting the spectrum
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

  write(out,*) "Rescaled the velocities."
  call flush(out)

  ! deallocate work arrays
  deallocate(e_spec, e_spec1, rr, hits, hits1, stat=ierr)

  return
end subroutine velocity_rescale
