subroutine io_write_4

  ! Writing out the velocities and scalars in X-space
  ! to the real*4 file

  use m_parameters
  use m_fields
  use m_work
  use x_fftw
  use m_les

  implicit none

  integer :: n_out, n, i, j, k
  real*8  :: wmag2, rkmax2

  ! every variable will undergo a mode truncation for all modes
  ! that are higher than kmax.  This will ensure that the written
  ! variables are isotropic

  rkmax2 = real(kmax,8)**2

  ! number of variables to write out
  n_out = 3
  if (int_scalars) n_out = n_out + n_scalars
  if (les .and. n_les>0) n_out = n_out + n_les

  ! putting all variables in wrk array
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+2
           wmag2 = akx(i)**2 + aky(k)**2 + akz(j)**2

           if (wmag2 .gt. rkmax2) then
              wrk(i,j,k,1:n_out) = zip
           else
              wrk(i,j,k,1:n_out) = fields(i,j,k,1:n_out)
           end if

        end do
     end do
  end do

  ! velocities
  call xFFT3d(-1,1)
  fname = 'u.'//file_ext
  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,1)
  call write_tmp4

  call xFFT3d(-1,2)
  fname = 'v.'//file_ext
  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,2)
  call write_tmp4

  call xFFT3d(-1,3)
  fname = 'w.'//file_ext
  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,3)
  call write_tmp4

  ! scalars
  if (int_scalars) then
     do n = 1,n_scalars
        call xFFT3d(-1,3+n)
        write(fname,"('sc',i2.2,'.',a6)") n,file_ext
        tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,3+n)
        call write_tmp4

     end do
  end if

  ! LES quantities
  if (les) then
     ! turbulent viscosity
     if (allocated(turb_visc)) then
        write(fname,"('nu_t.',a6)") file_ext
        tmp4 = turb_visc
        call write_tmp4
     end if

     if (n_les > 0) then
        do n = 1, n_les
           call xFFT3d(-1,3+n_scalars+n)
           write(fname,"('les',i1,'.',a6)") n,file_ext
           tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,3+n_scalars+n)
           call write_tmp4
        end do
     end if

  end if

  return
end subroutine io_write_4
