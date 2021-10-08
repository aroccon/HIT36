subroutine dealias_all

  use m_parameters
  use m_fields
  use x_fftw
  implicit none

  integer :: i,j,k
  real*8  :: akmax

  ! "2/3-rule".  For good explanation, see the following paper:
  ! R.S.Rogallo, "Numerical Exoeriments in Homogeneous Turbuelnce"
  ! NASA Technical Memorandum, NASA 1981.  
  ! Or contact authors of the code, they should have a PDF copy somewhere.
  !
  ! truncating all the modes in which at least one component of the k-vector
  ! has magnitude that is larger than N/3

  akmax = real(kmax,8)

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+2
           if  (abs(akx(i)) .gt. akmax .or. &
                abs(aky(k)) .gt. akmax .or. &
                abs(akz(j)) .gt. akmax) then
              fields(i,j,k,1:3+n_scalars) = zip
           end if

        end do
     end do
  end do

  return
end subroutine dealias_all
