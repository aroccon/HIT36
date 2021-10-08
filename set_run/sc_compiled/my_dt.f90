subroutine my_dt

  use m_openmpi
  use m_parameters
  implicit none

  integer :: n
  real*8  :: courant_max=0.1

  if (flow_type.le.1) then

     ! Courant number is calculated every iteration in rhsv.f

     if (courant.lt.courant_max) then
        DT = DT * 1.01d0
        DT = min(DT,0.01)
     else
        DT = DT/1.05d0
     end if

!----------------------------------------------------------------------
!  make sure DT is appropriate for scalars
!  DT < 0.09 * dx^2*Pe
!----------------------------------------------------------------------
     if (int_scalars) then
        do n=1,n_scalars
!!$           DT = min(DT,0.09*dx**2/PE(n))
           DT = min(DT,DT*courant_max/courant/sc(n))
        end do
     end if
  end if

  return
end subroutine my_dt
