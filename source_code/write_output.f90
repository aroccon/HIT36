subroutine write_output
use m_parameters
use m_fields
use m_work
use x_fftw
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
fname = '../results/u_'//file_ext//'.dat'
tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,1)
call write_tmp4

call xFFT3d(-1,2)
fname = '../results/v_'//file_ext//'.dat'
tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,2)
call write_tmp4

call xFFT3d(-1,3)
fname = '../results/w_'//file_ext//'.dat'
tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,3)
call write_tmp4

! scalars
if (int_scalars) then
  do n = 1,n_scalars
    call xFFT3d(-1,3+n)
    write(fname,"('../results/sc',i2.2,'.',a6)") n,file_ext
    tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,3+n)
    call write_tmp4
  end do
end if

return
end subroutine write_output
