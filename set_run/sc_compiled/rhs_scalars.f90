subroutine rhs_scalars

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw
  use m_timing
  use m_les

  implicit none

  integer :: i, j, k, n, n1, n2, nv, ns_lo, ns_hi
  real*8  :: rtmp1, rtmp2, wnum2, r11, r12, r21, r22, r31, r32

  ! calculate turbulent viscosity, if any
  if (les) call les_get_turb_visc

  ! If we're not advancing scalars, put the real-space velocities 
  ! in wrk1...3 and return
  ! same is done if dealias=0, that is, we use 2/3 rule.
  ! if dealias=1 (phase shifts) then this is done later in the subroutine

  ! also if doing LES and n_les (the # of les-related auxilary scalars) is
  ! greater than 0, then we need to transport these LES-related scalars,
  ! even if n_scalars=0.

  ! thus we do the obligatory part (tranfer velocities to X-space) and quit
  ! only when we are not transporting any scalars and if there are no
  ! LES-related scalars.
  if (.not.int_scalars .and. n_les .eq. 0) then
     ! converting velocities to the real space and returning
     wrk(:,:,:,1:3) = fields(:,:,:,1:3)
     do n = 1,3
        call xFFT3d(-1,n)
     end do
     return
  end if

  ! making the RHS for all scalars zero
  wrk(:,:,:,4:3+n_scalars+n_les) = zip

!--------------------------------------------------------------------------------
!  If dealias=0, performing the 2/3 rule dealiasing on scalars
!--------------------------------------------------------------------------------

  if (dealias.eq.0) then

     ! converting velocities to the real space
     wrk(:,:,:,1:3) = fields(:,:,:,1:3)
     do n = 1,3
        call xFFT3d(-1,n)
     end do

     ! Do each scalar one at a time.  Keep the velocities in wrk1:3 intact 
     ! because they are needed later.

     ! first we need to know which scalars do we want to transport.
     ! There are three cases:
     ! (0) No LES extra scalars, and passive scalars have not been initialized yet
     ! Thus this subroutine calculates the IFFT of velocities and exits.
     ! This is taken care of earlier.
     ! (1) Both passive scalars and LES extra scalars are transported
     ! This is possible when n_les > 0 and int_scalars=.true.
     ! (2) Only LES extra scalars are transported
     ! This is possible if and only if (.not.int_scalars .and. n_les>0)
     ! 
     ! The last two cases are taken care of by prescribing ns_lo and ns_hi, the 
     ! smallest and largest number of the scalar that needs to be transported.

     ns_lo = 1; 
     ns_hi = n_scalars + n_les
     if (.not.int_scalars) ns_lo = n_scalars + 1

     do n = ns_lo, ns_hi

        wrk(:,:,:,0) = fields(:,:,:,3+n)
        call xFFT3d(-1,0)

        ! Products of the scalar and velocities
        do i = 1,3
           wrk(:,:,:,n+2+i) = wrk(:,:,:,0) * wrk(:,:,:,i)
           call xFFT3d(1,n+2+i)
        end do

        ! Assembling the RHS in wrk(:,:,:,3+n)
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2

                 ! If the dealiasing option is 2/3-rule (dealias=0) then we retain the modes
                 ! inside the cube described by $| k_i | \leq  k_{max}$, $i=1,2,3$.
                 ! The rest of the modes is purged

                 if (ialias(i,j,k) .gt. 0) then
                    ! all the wavenumbers that are greater than kmax get zeroed out
                    wrk(i  ,j,k,3+n) = zip
                    wrk(i+1,j,k,3+n) = zip

                 else
                    ! taking the convective term, multiply it by "i" 
                    ! (see how it's done in x_fftw.f90)
                    ! and adding the diffusion term

                    ! also using the fact that the waveunmbers for (i,j,k) are the same
                    ! as wavenumbers for (i+1,j,k)

                    ! i * (a + ib) + d = -b + ia + d
                    rtmp1 =   akx(i+1)*wrk(i+1,j,k,3+n) + aky(k)*wrk(i+1,j,k,4+n) + akz(j)*wrk(i+1,j,k,5+n)
                    rtmp2 =   akx(i  )*wrk(i  ,j,k,3+n) + aky(k)*wrk(i  ,j,k,4+n) + akz(j)*wrk(i  ,j,k,5+n)

                    wnum2 = akx(i)**2 + aky(k)**2 + akz(j)**2

                    wrk(i  ,j,k,3+n) =   rtmp1 - pe(n) * wnum2*fields(i  ,j,k,3+n)
                    wrk(i+1,j,k,3+n) = - rtmp2 - pe(n) * wnum2*fields(i+1,j,k,3+n)

                 end if
              end do
           end do
        end do

        ! Now adding the reaction part (for scalars only, not for LES-related quantities
        ! that are formally scalars with indicies n_scalars+1...n_scalars+n_les )
        if (n .le. n_scalars) then
           if (scalar_type(n).ge.100) then
              call add_reaction(n)
              call dealias_rhs(3+n)
           end if
        end if

     end do

  end if

!--------------------------------------------------------------------------------
!  If dealias=1, performing the phase shift and truncation on scalars
!  The main ideology is as follows.  First we evaluate the phase-shifted
!  quantities, then not phase shifted.  This is done in order to have more
!  working space: at the beginning we have all work arrays available, but at the
!  end we must put sin/cos factors in wrk0 and u.v.w in real space in wrk1...3.
!  They are going to be used again in rhs_velocity later.
!--------------------------------------------------------------------------------
  phase_shifting_dealiasing: if (dealias.eq.1) then

     ! define the sin/cos factors that are used in phase shifting.
     ! computing sines and cosines for the phase shift of dx/2,dy/2,dz/2 
     ! and putting them into wrk0
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              wrk(i  ,j,k,0) = cos(half*(akx(i  )+aky(k)+akz(j))*dx)
              wrk(i+1,j,k,0) = sin(half*(akx(i+1)+aky(k)+akz(j))*dx)
           end do
        end do
     end do

     ! Doing all scalars at the same time.  We can do this because the
     ! number of work arrays that are available to us is 3+n+2.  First
     ! three later will be taken by  velocities, and the last two are
     ! primary work arrays here.

     ! First, get phase shifted velocities and scalars
     do n = 1, 3 + n_scalars + n_les
        ! phase-shifting the quantity
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2
                 wrk(i  ,j,k,n) = fields(i  ,j,k,n) * wrk(i,j,k,0) - fields(i+1,j,k,n) * wrk(i+1,j,k,0)
                 wrk(i+1,j,k,n) = fields(i+1,j,k,n) * wrk(i,j,k,0) + fields(i  ,j,k,n) * wrk(i+1,j,k,0)
              end do
           end do
        end do
        ! transforming it to real space
        call xFFT3d(-1,n)
     end do

     ! now we have two vacant arrays: n_scalars+n_les+4 and n_scalars+n_les+5.  Work in them
     n1 = n_scalars + n_les + 4 
     n2 = n_scalars + n_les + 5

     ! do one scalar at a time
     phase_shifted_rhs: do n = 4, 3 + n_scalars + n_les

        ! getting all three products of phase-shifted scalar and phase-shifted velocities
        ! using three work arrays: n, n1 and n2
        wrk(:,:,:,n1) = wrk(:,:,:,n) * wrk(:,:,:,1)
        wrk(:,:,:,n2) = wrk(:,:,:,n) * wrk(:,:,:,2)
        wrk(:,:,:,n ) = wrk(:,:,:,n) * wrk(:,:,:,3)
        ! transforming them to Fourier space
        call xFFT3d(1,n1)  
        call xFFT3d(1,n2)  
        call xFFT3d(1,n )

        ! phase shifting them back using wrk0 and adding -0.5*(ik) to the RHS for the scalar
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2

                 if (ialias(i,j,k) .gt. 1) then
                    wrk(i:i+1,j,k,n) = zip
                 else
                    ! (u+)*(phi+) phase shifted back
                    r11 = wrk(i  ,j,k,n1) * wrk(i,j,k,0) + wrk(i+1,j,k,n1) * wrk(i+1,j,k,0)
                    r12 = wrk(i+1,j,k,n1) * wrk(i,j,k,0) - wrk(i  ,j,k,n1) * wrk(i+1,j,k,0)
                    ! (v+)*(phi+) phase shifted back
                    r21 = wrk(i  ,j,k,n2) * wrk(i,j,k,0) + wrk(i+1,j,k,n2) * wrk(i+1,j,k,0)
                    r22 = wrk(i+1,j,k,n2) * wrk(i,j,k,0) - wrk(i  ,j,k,n2) * wrk(i+1,j,k,0)
                    ! (w+)*(phi+) phase shifted back
                    r31 = wrk(i  ,j,k,n ) * wrk(i,j,k,0) + wrk(i+1,j,k,n ) * wrk(i+1,j,k,0)
                    r32 = wrk(i+1,j,k,n ) * wrk(i,j,k,0) - wrk(i  ,j,k,n ) * wrk(i+1,j,k,0)
                    ! adding -0.5*(ik)*(the result) to the RHSs for the scalar
                    wrk(i  ,j,k,n) = + 0.5d0 * ( akx(i+1)*r12 + aky(k)*r22 + akz(j)*r32 )
                    wrk(i+1,j,k,n) = - 0.5d0 * ( akx(i  )*r11 + aky(k)*r21 + akz(j)*r31 )
                 end if

              end do
           end do
        end do

     end do phase_shifted_rhs

     ! at this moment wrk4...3+n_scalars+n_les contain the half of the convective term, which was
     ! obtained from the phase shifted quantities.  Now we need to add the other half of the
     ! convective term (the one that is obtained by multiplication of no-phase-shifted stuff)

     ! first get the velocities into the real space and put them in wrk1...3
     ! this should remain in there untouched, to be used in rhs_velocity later
     do n = 1,3
        wrk(:,:,:,n) = fields(:,:,:,n)
        call xFFT3d(-1,n)
     end do

     ! now do scalars one at a time since we don't have enough storage to do them 
     ! all at once

     not_phase_shifted_rhs: do n = 4, 3 + n_scalars + n_les

        ! get the scalar into the real space and put it in the wrk(0)
        wrk(:,:,:,0) = fields(:,:,:,n)
        call xFFT3d(-1,0)

        ! calculate the product of the scalar with velocitiy components
        wrk(:,:,:,n1) = wrk(:,:,:,0) * wrk(:,:,:,1)
        wrk(:,:,:,n2) = wrk(:,:,:,0) * wrk(:,:,:,2)
        wrk(:,:,:, 0) = wrk(:,:,:,0) * wrk(:,:,:,3)
        ! transform it to the Fourier space
        call xFFT3d(1,n1)
        call xFFT3d(1,n2)
        call xFFT3d(1, 0)

        ! add the -0.5*(ik)*results to the RHS, along with the diffusion term
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2

                 if (ialias(i,j,k) .lt. 2) then

                    rtmp1 =   akx(i+1)*wrk(i+1,j,k,n1) + aky(k)*wrk(i+1,j,k,n2) + akz(j)*wrk(i+1,j,k,0)
                    rtmp2 =   akx(i  )*wrk(i  ,j,k,n1) + aky(k)*wrk(i  ,j,k,n2) + akz(j)*wrk(i  ,j,k,0)

                    wnum2 = akx(i)**2 + aky(k)**2 + akz(j)**2
                    wrk(i  ,j,k,n) = wrk(i  ,j,k,n) + 0.5d0 * rtmp1 - pe(n-3) * wnum2*fields(i  ,j,k,n)
                    wrk(i+1,j,k,n) = wrk(i+1,j,k,n) - 0.5d0 * rtmp2 - pe(n-3) * wnum2*fields(i+1,j,k,n)

                 end if
              end do
           end do
        end do

     end do not_phase_shifted_rhs

     ! ------------------------------------------------------------
     ! Add the reaction rates to the RHS     
     ! 
     ! The LES-related scalars are not affected because the 
     ! upper bound for n is 3+n_scalars, not 3+n_scalars+n_les
     ! -----------------------------------------------------------

     reaction_rates: do n = 4, 3+n_scalars

        if (scalar_type(n-3) .gt. 300) then

!!$           ! putting phase shifted scalar in wrk(n1)
!!$           do k = 1,nz
!!$              do j = 1,ny
!!$                 do i = 1,nx+1,2
!!$                    wrk(i  ,j,k,n1) = fields(i  ,j,k,n) * wrk(i,j,k,0) - fields(i+1,j,k,n) * wrk(i+1,j,k,0)
!!$                    wrk(i+1,j,k,n1) = fields(i+1,j,k,n) * wrk(i,j,k,0) + fields(i  ,j,k,n) * wrk(i+1,j,k,0)
!!$                 end do
!!$              end do
!!$           end do
!!$           ! transforming it to real space
!!$           call xFFT3d(-1,n1)
!!$           ! putting the phase-shifted reaction rate in wrk(n2)
!!$           call scalar_reaction_rate(n1,n2)
!!$           ! transforming to Fourier space
!!$           call xFFT3d(1,n2)
!!$           ! phase shifting back and adding a half of it to the RHS
!!$           do k = 1,nz
!!$              do j = 1,ny
!!$                 do i = 1,nx+1,2
!!$                    if (ialias(i,j,k) .le. 1) then
!!$                       ! phase shifting back
!!$                       r11 = wrk(i  ,j,k,n2) * wrk(i,j,k,0) + wrk(i+1,j,k,n2) * wrk(i+1,j,k,0)
!!$                       r12 = wrk(i+1,j,k,n2) * wrk(i,j,k,0) - wrk(i  ,j,k,n2) * wrk(i+1,j,k,0)
!!$                       ! adding 0.5*(the result) to the RHSs for the scalar
!!$                       wrk(i  ,j,k,n) = wrk(i  ,j,k,n) + 0.5d0 * r11
!!$                       wrk(i+1,j,k,n) = wrk(i+1,j,k,n) + 0.5d0 * r12
!!$                    end if
!!$                 end do
!!$              end do
!!$           end do
!!$
!!$           ! scond part: doing the same  thing with not-phase-shifted scalar
!!$           ! putting it in wrk(n1)
!!$           wrk(:,:,:,n1) = fields(:,:,:,n)
!!$           ! transforming it to real space
!!$           call xFFT3d(-1,n1)
!!$           ! putting the phase-shifted reaction rate in wrk(n2)
!!$           call scalar_reaction_rate(n1,n2)
!!$           ! transforming to Fourier space
!!$           call xFFT3d(1,n2)
!!$           ! phase shifting back and adding a half of it to the RHS
!!$           do k = 1,nz
!!$              do j = 1,ny
!!$                 do i = 1,nx+1,2
!!$                    if (ialias(i,j,k) .le. 1) then
!!$                       ! phase shifting back
!!$                       r11 = wrk(i  ,j,k,n2) * wrk(i,j,k,0) + wrk(i+1,j,k,n2) * wrk(i+1,j,k,0)
!!$                       r12 = wrk(i+1,j,k,n2) * wrk(i,j,k,0) - wrk(i  ,j,k,n2) * wrk(i+1,j,k,0)
!!$                       ! adding 0.5*(the result) to the RHSs for the scalar
!!$                       wrk(i  ,j,k,n) = wrk(i  ,j,k,n) + 0.5d0 * r11
!!$                       wrk(i+1,j,k,n) = wrk(i+1,j,k,n) + 0.5d0 * r12
!!$                    end if
!!$                 end do
!!$              end do
!!$           end do

           if (scalar_type(n-3) .eq. 311) then

              ! since the reaction rate is cubic, we need to apply some severe truncation.

              ! putting the scalar in wrk(n1)
              wrk(:,:,:,n1) = fields(:,:,:,n)
              ! truncating it so only  the modes with |k_i| < nx/4 remain
              do k = 1,nz
                 do j = 1,ny
                    do i = 1,nx+1,2
                       if (abs(akx(i)).gt.nx/4 .or. abs(aky(k)).gt.nx/4 .or. abs(akz(j)).gt.nx/4) then
                          wrk(i  ,j,k,n1) = zip
                          wrk(i+1,j,k,n1) = zip
                       end if
                    end do
                 end do
              end do

              ! transforming it to real space
              call xFFT3d(-1,n1)

              ! self-adjusting bistable reaction needs the mean value of the scalar
              rtmp1 = fields(1,1,1,n)/nxyz_all
              call MPI_BCAST(rtmp1, 1, MPI_REAL8, 0, MPI_COMM_TASK, mpi_err)

              ! getting the reaction rate
              wrk(:,:,:,n2) = reac_sc(n-3) * (one - wrk(:,:,:,n1)**2) * (wrk(:,:,:,n1) - rtmp1)

              ! transforming it to Fourier space
              call xFFT3d(1,n2)

              ! adding to the RHS
              wrk(:,:,:,n) = wrk(:,:,:,n) + wrk(:,:,:,n2)

           else
              write(out,*) "The scalar type has a reaction rate that is not supported yet:", scalar_type(n-3)
              call flush(out)
              call my_exit(-1)
           end if
        end if
     end do reaction_rates

  end if phase_shifting_dealiasing

     ! special case - passive scalar with the uniform gradient as a source
     ! adding the source term - the first component of velocity, because we assume
     ! that the uniform gradient has slope 1 and direction in the x-direction
     gradient_source: do n = 1, n_scalars
        if (scalar_type(n) .eq. 0) wrk(:,:,:,n+3) = wrk(:,:,:,n+3) - fields(:,:,:,1)
     end do gradient_source



!--------------------------------------------------------------------------------
!  Add LES to the RHS of all the scalars
!--------------------------------------------------------------------------------
  les_active: if (les) then
     call les_rhs_scalars
  end if les_active


  return
end subroutine rhs_scalars

!================================================================================
!================================================================================
subroutine add_reaction(n)

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: n, rtype
  real*8  :: scmean, rrate

  ! reaction type
  rtype =  scalar_type(n)/100

  ! raction rate 
  rrate = reac_sc(n)

  select case (rtype)
  case (1)

     ! KPP reaction rate
     wrk(:,:,:,0) = rrate * (1.d0 - wrk(:,:,:,0)**2)

  case (2)

     ! symmetric bistable       
     wrk(:,:,:,0) = rrate * (1.d0 - wrk(:,:,:,0)**2) * wrk(:,:,:,0)

  case (3)

     ! self-adjusting bistable
     scmean = fields(1,1,1,3+n)/nxyz_all
     call MPI_BCAST(scmean, 1, MPI_REAL8, 0, MPI_COMM_TASK, mpi_err)

     wrk(:,:,:,0) = rrate * (1.d0 - wrk(:,:,:,0)**2) * &
          (wrk(:,:,:,0) - scmean)
  case default

     write(out,*) "Unknown reaction rate"
     call flush(out)

     stop

  end select


  ! FFT the reaction into the Fourier space
  call xFFT3d(1,0)

  ! Adding reaction to the RHS in wrk(:,:,:,3+n)

  wrk(:,:,:,3+n) = wrk(:,:,:,3+n) + wrk(:,:,:,0)

end subroutine add_reaction

!================================================================================

subroutine dealias_rhs(n)

  use m_io
  use m_parameters
  use m_work
  use x_fftw

  implicit none

  integer :: i, j, k, n
  real*8  :: wnum2, akmax

  akmax = real(kmax,8)

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+1,2

           if ( abs(akx(i)).gt.akmax .or. &
                abs(aky(k)).gt.akmax .or. &
                abs(akz(j)).gt.akmax ) then

              wrk(i  ,j,k,n) = zip
              wrk(i+1,j,k,n) = zip
           end if

        end do
     end do
  end do

  return

end subroutine dealias_rhs

!================================================================================
!================================================================================
!================================================================================
!================================================================================
subroutine test_rhs_scalars

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i,j,k, n
  real*8 :: a,b,c, x,y,z
  if (task.eq.'hydro') then

     a = 1.d0
     b = 5.d0
     c = 17.d0

     do k = 1,nz
        do j = 1,ny
           do i = 1,nx

              x = dx*real(i-1)
              y = dx*real(j-1)
              z = dx*real(myid*nz + k-1)

              wrk(i,j,k,1) = sin(a * x)
              wrk(i,j,k,2) = sin(b * y)
              wrk(i,j,k,3) = sin(c * z)
              wrk(i,j,k,4) = cos(a * x)
           end do
        end do
     end do


     do n = 1,4
        call xFFT3d(1,n)
        fields(:,:,:,n) = wrk(:,:,:,n)
     end do


     nu = .5d0

     call rhs_scalars

     print *,'got rhs'

     call xFFT3d(-1,4)


     do k = 1,nz
        do j = 1,ny
           do i = 1,nx

              x = dx*real(i-1)
              y = dx*real(j-1)
              z = dx*real(myid*nz + k-1)


              ! checking 
              wrk(i,j,k,0) = -a*cos(2.*a*x) - cos(a*x)*(b*cos(b*y) + c*cos(c*z) + nu*a**2)

           end do
        end do
     end do

!!$     tmp4(:,:,:) = wrk(1:nx,:,:,4)
!!$     fname = 'r1.arr'
!!$     call write_tmp4
!!$
!!$     tmp4(:,:,:) = wrk(1:nx,:,:,0)
!!$     fname = 'r0.arr'
!!$     call write_tmp4




     wrk(:,:,:,0) = abs(wrk(1:nx,:,:,0) - wrk(1:nx,:,:,4))

     print *,'Maximum error is ',maxval(wrk(1:nx,:,:,0))

!!$     tmp4(:,:,:) = wrk(1:nx,:,:,3) - wrk(1:nx,:,:,6)
!!$     fname = 'e3.arr'
!!$     call write_tmp4
!!$


  end if
  return
end subroutine test_rhs_scalars
