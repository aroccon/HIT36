!================================================================================
!  m_les - the module that contains all subroutines and variables that are
!  needed to introduce Large Eddy Simulation (LES) into the code.
!
!  The behaviour of the module is governed by the variable "les_mode" from the
!  module m_parameters.f90
!
!  Time-stamp: <2010-03-18 13:52:33 (chumakov)>
!================================================================================
module m_les

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw
  implicit none


!================================================================================
!  Arrays and constants
!================================================================================

  ! indicator, whether we use LES at all
  logical :: les = .false.

  ! model switch "les_mode" is contained by m_parameters
  ! integer(kind=4) :: les_model

  ! array for turbulent viscosity
  real(kind=8), allocatable :: turb_visc(:,:,:)

  ! LES sources for velocities
  real(kind=8), allocatable :: vel_source_les(:,:,:,:)

  ! separate array for the subgrid stress
  real(kind=8), allocatable :: tauij(:,:,:,:)

  ! Smagorinsky constant
  real(kind=8) :: c_smag = 0.18

  ! Scaling constant for the lag-model
  real(kind=8) :: C_T = 1.d0

  ! Scaling constant for the mixed model
  real(kind=8) :: C_mixed

  ! test filter width
  real(kind=8) :: les_delta

  ! array with the test filter in it
  real(kind=8), allocatable :: filter_g(:,:,:)

  ! model indicator for the output
  character*3 :: les_model_name = '   '

  ! TEMP variables:
  ! - production
  real(kind=8) :: energy, production, B, dissipation

!================================================================================
!                            SUBROUTINES
!================================================================================
contains
!================================================================================
!  Allocation of LES arrays
!================================================================================
  subroutine m_les_init

    implicit none

    integer :: n
    real*8, allocatable :: sctmp(:)


    ! if les_model=0, do not initialize anything and return
    if (les_model==0) return

    write(out,*) 'Initializing LES...'
    call flush(out)

    ! depending on the value of les_model, initialize different things
    les = .true.
    ! initialize the filter width to be equal to the grid spaxing
    les_delta = dx
    write(out,*) 'LES_DELTA = ',les_delta
    ! initializeing stuff based on the model switch "les_model"
    select case (les_model)

!--------------------------------------------------------------------------------
!  Smagorinsky model
!--------------------------------------------------------------------------------
    case(1)
       write(out,*) ' - Smagorinsky model: initializing the eddy viscosity'
       call flush(out)
       allocate(turb_visc(nx,ny,nz),stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 0
       les_model_name = " SM"

!--------------------------------------------------------------------------------
!  Dynamic localization model
!--------------------------------------------------------------------------------
    case(2)
       write(out,*) ' - DL model: initializing the eddy viscosity and adding extra transport equation'
       call flush(out)
       allocate(turb_visc(nx,ny,nz),stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 1
       les_model_name = "DLM"

!--------------------------------------------------------------------------------
!  Dynamic localization model + lag model for the subgrid energy dissipation
!--------------------------------------------------------------------------------
    case(3)
       ! Dynamic Localization model + lag model for the dissipation
       write(out,*) ' - DL model + lag-model for dissipation'
       call flush(out)
       allocate(turb_visc(nx,ny,nz),stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 3
       les_model_name = "DLL"

!--------------------------------------------------------------------------------
!  Dynamic Structure model + algebraic model for dissipation
!--------------------------------------------------------------------------------

    case(4)

       write(out,*) ' - DSt model + algebraic model for dissipation'
       call flush(out)
       allocate(turb_visc(nx,ny,nz), vel_source_les(nx+2,ny,nz,3), stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 1
       les_model_name = "DST"

       ! Fot this model we need a filter.  The filter cannot be initialized
       ! without previous initialization of fields array.  So initialization of
       ! the filter is done in m_les_begin

       ! Note that for this model we need a bigger wrk array
       ! this is taken care of in m_work.f90

!--------------------------------------------------------------------------------
!  Dynamic Structure model + lag model for dissipation
!--------------------------------------------------------------------------------
    case(5)

       write(out,*) ' - DSt model + lag model for dissipation'
       call flush(out)
       allocate(turb_visc(nx,ny,nz), vel_source_les(nx+2,ny,nz,3), stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 3
       les_model_name = "STL"

       ! Fot this model we need a filter.  The filter cannot be initialized
       ! without previous initialization of fields array.  So initialization of
       ! the filter is done in m_les_begin

       ! Note that for this model we need a bigger wrk array
       ! this is taken care of in m_work.f90

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
       inquire(file = 'c_t.in', exist = there)
       if (.not.there) then
          write(out,*) "Cannot find the file 'c_t.in', exiting"
          call my_exit(-1)
       end if
       if (iammaster) then
          open(900,file='c_t.in')
          read(900,*) C_T
          close(900)
       end if
       count = 1
       call MPI_BCAST(C_T,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
       write(out,*) "DSTM + lag model WITH C_T = ",C_T
       call flush(out)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-

!--------------------------------------------------------------------------------
!  Mixed model: Dynamic Structure + C-mixed * Dynamic Localization model
!  Dissipation model is algebraic
!--------------------------------------------------------------------------------
    case(6)

       write(out,*) ' - MIXED MODEL: DSTM + DLM'
       write(out,*) ' -              Dissipation is algebraic'
       call flush(out)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
       inquire(file = 'c_mixed.in', exist = there)
       if (.not.there) then
          write(out,*) "Cannot find the file 'c_mixed.in', exiting"
          call my_exit(-1)
       end if
       if (iammaster) then
          open(900,file='c_mixed.in')
          read(900,*) C_mixed
          close(900)
       end if
       count = 1
       call MPI_BCAST(C_mixed,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
       write(out,*) "MIXED MODEL WITH C_MIXED = ",C_mixed
       call flush(out)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-


       allocate(turb_visc(nx,ny,nz), vel_source_les(nx+2,ny,nz,3), stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 1
       les_model_name = "MMA"

       ! Fot this model we need a filter.  The filter cannot be initialized
       ! without previous initialization of fields array.  So initialization of
       ! the filter is done in m_les_begin

       ! Note that for this model we need a bigger wrk array
       ! this is taken care of in m_work.f90


!--------------------------------------------------------------------------------
!  Mixed model: Dynamic Structure + some viscosity (about 15% of the usual)
!  Dissipation model is lag-model
!--------------------------------------------------------------------------------
    case(7)

       write(out,*) ' - MIXED MODEL: DSTM + eddy viscosity'
       write(out,*) ' -              Lag-model for Dissipation'
       call flush(out)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
       inquire(file = 'c_mixed.in', exist = there)
       if (.not.there) then
          write(out,*) "Cannot find the file 'c_mixed.in', exiting"
          call my_exit(-1)
       end if
       if (iammaster) then
          open(900,file='c_mixed.in')
          read(900,*) C_mixed
          close(900)
       end if
       count = 1
       call MPI_BCAST(C_mixed,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
       write(out,*) "MIXED MODEL WITH C_MIXED = ",C_mixed
       call flush(out)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-


       allocate(turb_visc(nx,ny,nz), vel_source_les(nx+2,ny,nz,3), stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 3
       les_model_name = "MML"

       ! Fot this model we need a filter.  The filter cannot be initialized
       ! without previous initialization of fields array.  So initialization of
       ! the filter is done in m_les_begin

       ! Note that for this model we need a bigger wrk array
       ! this is taken care of in m_work.f90

!--------------------------------------------------------------------------------
!  LES MODEL # 10
!
!  One-equation mixed model for the momentum closure (DSTM + C_mixed * SM)
!  Lag model for the closure of the k-equation
!  Harlow model for closure of the scalar tranport equation
!--------------------------------------------------------------------------------
    case(10)
       write(out,*) '   LES MODEL # 10'
       write(out,*) ' - MIXED MODEL: DSTM + eddy viscosity'
       write(out,*) '                Lag-model for Dissipation'
       write(out,*) ' - SCALARS: Harlow model'
       write(out,*) ' ============> file c_mixed.in is still required'
       call flush(out)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
       inquire(file = 'c_mixed.in', exist = there)
       if (.not.there) then
          write(out,*) "Cannot find the file 'c_mixed.in', making default 0.5"
          c_mixed = 0.5d0
       else
          if (iammaster) then
             open(900,file='c_mixed.in')
             read(900,*) C_mixed
             close(900)
          end if
          count = 1
          call MPI_BCAST(C_mixed,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
       end if
       write(out,*) "MIXED MODEL WITH C_MIXED = ",C_mixed
       call flush(out)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-


       allocate(turb_visc(nx,ny,nz), tauij(nx+2,ny,nz,6), stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the auxiliary arrays'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = zip
       tauij = zip

       ! number of extra LES quantities is 2: (BT), (eps * T)
       n_les = 2
       les_model_name = "#10"

       ! For this model we need a filter.  The filter cannot be initialized
       ! without previous initialization of fields array.  So initialization of
       ! the filter is done in m_les_begin

       ! Note that for this model we need a bigger wrk array
       ! this is taken care of in m_work.f90

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------


    case default
       write(out,*) 'M_LES_INIT: invalid value of les_model:',les_model
       call flush(out)
       call my_exit(-1)
    end select

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

    write(out,*) "n_les = ", n_les

    ! If the number of LES quantities is non-zero, then for the sake of modularity
    ! we need to have SC and PE numbers defined for those quantities.  For this we
    ! need to re-allocate the arrays SC and PE
    additional_scalars: if (n_les .gt. 0) then
       write(out,*) "Adding elements to arrays PE and SC for the LES-related scalars..."
       call flush(out)
       allocate(sctmp(1:n_scalars+n_les), stat=ierr)
       sctmp(1:n_scalars) = sc(1:n_scalars)
       sctmp(n_scalars+1:n_scalars+n_les) = one
       if (allocated(sc)) deallocate(sc)
       allocate(sc(1:n_scalars+n_les))
       sc = sctmp
       if (allocated(pe)) deallocate(pe)
       allocate(pe(1:n_scalars+n_les))
       pe = nu / sc
       deallocate(sctmp)
       write(out,*) " ...done."
       call flush(out)
    end if additional_scalars

    write(out,*) 'initialized LES.'
    call flush(out)

    return
  end subroutine m_les_init

!================================================================================
!================================================================================
!  initialization of the LES arrays - part 2
!  definition of the arrays
!
!  called after the restart
!================================================================================
  subroutine m_les_begin

    use m_filter_xfftw
    implicit none

    write(out,*) "M_LES_BEGIN..."
    call flush(out)

    select case(les_model)

!--------------------------------------------------------------------------------
!  if les_model=0, do not initialize anything and return
!  also don't do anything if the model is Smagorinsky model - it does not
!  need any additional initialization beside the array allocation which has
!  been already done.
!--------------------------------------------------------------------------------
    case(0) 
       return
    case(1)
       return

!--------------------------------------------------------------------------------
!  Dynamic localization model
!--------------------------------------------------------------------------------
    case(2)
       write(out,*) "-- DLM model"
       call flush(out)

       if (itime.eq.0) then
          write(out,*) "-- Initializing k_sgs"
          call flush(out)
          call m_les_dlm_k_init
       end if

!--------------------------------------------------------------------------------
    case(3)
       write(out,*) "-- DLM model with lag model for dissipation"
       call flush(out)

       if (itime.gt.0) return

       ! call m_les_dlm_k_init

       ! COMMENT OUT THE FOLLOWING IN ORDER TO INITIALIZE B AND EPSILON AS ZERO

       ! making initial epsilon = k^(3/2)/Delta everywhere
       ! since k=const=0.5, just change one entry in epsilon array
       ! Note that the array itself contains not epsilon but (epsilon * T_epsilon)
       ! so the contents of the array is not presicely k^(3/2)/Delta. Some math is involved.
       ! (comment out to start from zero dissipation)
       ! if (iammaster) fields(1,1,1,3+n_scalars+3) = C_T * 0.5 * real(nxyz_all)

       ! Now initial conditions for B.  We want B to be same as epsilon
       ! (kind of "starting from steady state"), but again the array contains B*T_B, not just B.
       ! Current implementation is T_B = 1/|S|.  To get |S|, we call m_les_src_k_dlm, which
       ! gives us |S|^2 in wrk0.
       ! call m_les_k_src_dlm
       ! now wrk0 contains |S|^2 in x-space, and we can use it to get B
       ! fields(:,:,:,3+n_scalars+2) = 0.5**1.5d0 / les_delta / sqrt(wrk(:,:,:,0))
       ! call xFFT3d_fields(1,3+n_scalars+2)

       ! INITIALIZING K_SGS, B*T_B and eps*T_eps
       ! Initializing them so that k_sgs = B*T_b + eps*T_eps
       ! in fact this makes the equation for k_sgs unnecessary but we're keeping it for
       ! debug purposes and such

       write(out,*) "   Initializing k_sgs = 0.1, B*T_B = eps*T_eps = 0.05"
       call flush(out)
       ! definition of k_sgs = 0.1 everywhere
       if (iammaster) fields(1,1,1,3+n_scalars+1) = 0.1d0 * real(nxyz_all)
       ! definition of eps*T_eps = B*T_B = 0.5 k_sgs
       if (iammaster) fields(1,1,1,3+n_scalars+2) = 0.5d0 * fields(1,1,1,3+n_scalars+1)
       if (iammaster) fields(1,1,1,3+n_scalars+3) = 0.5d0 * fields(1,1,1,3+n_scalars+1)

!--------------------------------------------------------------------------------

    case(4)

       write(out,*) "-- Dynamic Structure model with algebraic model for dissipation"
       call flush(out)

       ! initializing filtering arrays
       ! because filter_xfftw_init uses fields(1) as a temporary array, we need
       ! to store it before we initialize the filter, and the restore it to
       ! what it was.
       wrk(:,:,:,LBOUND(wrk,4)) = fields(:,:,:,LBOUND(fields,4))
       call filter_xfftw_init
       fields(:,:,:,LBOUND(fields,4)) = wrk(:,:,:,LBOUND(wrk,4))

       if (itime.eq.0) then
          write(out,*) "-- Initializing k_sgs = 0.2"
          call flush(out)

          if (iammaster) fields(1,1,1,3+n_scalars+1) = 0.2d0 * real(nxyz_all)
       end if

!--------------------------------------------------------------------------------

    case(5)

       write(out,*) "-- Dynamic Structure model with lag-model for dissipation"
       call flush(out)

       ! initializing filtering arrays
       ! because filter_xfftw_init uses fields(1) as a temporary array, we need
       ! to store it before we initialize the filter, and the restore it to
       ! what it was.
       wrk(:,:,:,LBOUND(wrk,4)) = fields(:,:,:,LBOUND(fields,4))
       call filter_xfftw_init
       fields(:,:,:,LBOUND(fields,4)) = wrk(:,:,:,LBOUND(wrk,4))


       if (itime.eq.0) then
          write(out,*) "-- Initializing k_sgs = 0.5"
          call flush(out)

          if (iammaster) fields(1,1,1,3+n_scalars+1) = 0.5d0 * real(nxyz_all)
          ! definition of eps*T_eps = 0, B*T_B = k_sgs
          if (iammaster) fields(1,1,1,3+n_scalars+2) = fields(1,1,1,3+n_scalars+1)
          if (iammaster) fields(1,1,1,3+n_scalars+3) = 0.d0*fields(1,1,1,3+n_scalars+1)
       end if

!--------------------------------------------------------------------------------

    case(6)

       write(out,*) "-- MIXED MODEL (Dynamic Structure model + Smagorinsky)"
       write(out,*) "               Algebraic model for dissipation"
       call flush(out)

       ! initializing filtering arrays
       ! because filter_xfftw_init uses fields(1) as a temporary array, we need
       ! to store it before we initialize the filter, and the restore it to
       ! what it was.
       write(out,*) "Initializing filter"
       call flush(out)
       wrk(:,:,:,LBOUND(wrk,4)) = fields(:,:,:,LBOUND(fields,4))
       call filter_xfftw_init
       fields(:,:,:,LBOUND(fields,4)) = wrk(:,:,:,LBOUND(wrk,4))

       if (itime.eq.0) then

          write(out,*) "-- Initializing k_sgs = 0.5"
          call flush(out)

          if (iammaster) fields(1,1,1,3+n_scalars+1) = 0.5d0 * real(nxyz_all)

       end if

!--------------------------------------------------------------------------------

    case(7)

       write(out,*) "-- MIXED MODEL (Dynamic Structure model + Smagorinsky)"
       write(out,*) "               Lag-model for dissipation"
       call flush(out)

       ! initializing filtering arrays
       ! because filter_xfftw_init uses fields(1) as a temporary array, we need
       ! to store it before we initialize the filter, and the restore it to
       ! what it was.
       write(out,*) "   Initializing filter"
       call flush(out)
       wrk(:,:,:,LBOUND(wrk,4)) = fields(:,:,:,LBOUND(fields,4))
       call filter_xfftw_init
       fields(:,:,:,LBOUND(fields,4)) = wrk(:,:,:,LBOUND(wrk,4))


       if (itime.eq.0) then

          write(out,*) "-- Initializing k_sgs = 0.5, (BT)=k_s, (Eps T) = 0"
          call flush(out)

          ! Initializing k
          if (iammaster) fields(1,1,1,3+n_scalars+1) = 0.5d0 * real(nxyz_all)
          ! initializing (BT)
          if (iammaster) fields(1,1,1,3+n_scalars+2) = fields(1,1,1,3+n_scalars+1)
          ! initializing (eps T)
          if (iammaster) fields(1,1,1,3+n_scalars+3) = 0.d0*fields(1,1,1,3+n_scalars+1)

       end if

!--------------------------------------------------------------------------------

    case(10)

       write(out,*) "-- MIXED MODEL (Dynamic Structure model + Smagorinsky)"
       write(out,*) "               Lag-model for dissipation"
       call flush(out)

       ! initializing filtering arrays
       ! because filter_xfftw_init uses fields(1) as a temporary array, we need
       ! to store it before we initialize the filter, and the restore it to
       ! what it was.
       write(out,*) "   Initializing filter"
       call flush(out)
       wrk(:,:,:,LBOUND(wrk,4)) = fields(:,:,:,LBOUND(fields,4))
       call filter_xfftw_init
       fields(:,:,:,LBOUND(fields,4)) = wrk(:,:,:,LBOUND(wrk,4))

       if (itime.eq.0) then
          write(out,*) "-- (BT) = 0.1, (Eps T) = 0.1"
          call flush(out)
          ! No need to initialize k_s, since k_s = (BT)+(eps T)
          ! if (iammaster) fields(1,1,1,3+n_scalars+1) = 0.5d0 * real(nxyz_all)
          ! initializing (BT)
          if (iammaster) fields(1,1,1,3+n_scalars+1) = 0.1d0 * real(nxyz_all)
          ! initializing (eps T)
          if (iammaster) fields(1,1,1,3+n_scalars+2) = 0.1d0 * real(nxyz_all)
       end if

!--------------------------------------------------------------------------------

    case default
       write(out,*) "M_LES_BEGIN: invalid value of les_model: ",les_model
       call flush(out)
       call my_exit(-1)
    end select


    return
  end subroutine m_les_begin

!================================================================================
!================================================================================
!  Adding LES sources to the RHS of velocities
!================================================================================
  subroutine les_rhs_velocity

    use x_fftw
    implicit none

    integer :: i, j, k

    select case (les_model)
    case(1:3)
       call les_rhsv_turb_visc
    case(4:5)
       ! Dynamic Structure model
       ! add the velocity sources to RHS for velocitieis
       wrk(:,:,:,1:3) = wrk(:,:,:,1:3) + vel_source_les(:,:,:,1:3)

    case(6)
       ! Mixed model (DSTM + DLM)
       ! add the velocity sources to RHS for velocitieis
       wrk(:,:,:,1:3) = wrk(:,:,:,1:3) + vel_source_les(:,:,:,1:3)
       ! also apply turbulent viscosity
       call les_rhsv_turb_visc

    case(7)
       ! Mixed model (DSTM + DLM) + lag-model for dissipation of k_sgs
       ! add the velocity sources to RHS for velocitieis
       wrk(:,:,:,1:3) = wrk(:,:,:,1:3) + vel_source_les(:,:,:,1:3)
       ! also apply turbulent viscosity to velocities
       call les_rhsv_turb_visc

    case(10)
       ! Mixed model (DSTM + DLM) + lag-model for dissipation of k_sgs
       ! add the velocity sources to RHS for velocitieis
       ! The sources are directly computed from the array tauij
       ! so we do not need a separate array for them
       call les_add_vel_source_from_tauij

    case default
       write(out,*) 'LES_RHS_VELOCITY: invalid value of les_model:',les_model
       call flush(out)
       call my_exit(-1)
    end select

!--------------------------------------------------------------------------------
!  Making sure that we are not getting any aliasing errors. That is done by 
!  making RHS zero for wavenumbers that are aliased.
!--------------------------------------------------------------------------------
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 2
             if (ialias(i,j,k) .ne. 0) then
                wrk(i,j,k,1:3) = zip
             end if
          end do
       end do
    end do
!--------------------------------------------------------------------------------

    return
  end subroutine les_rhs_velocity

!================================================================================
!================================================================================
!  Adding LES sources to the RHS of scalars
!================================================================================
  subroutine les_rhs_scalars
    use x_fftw
    implicit none

    integer :: n, i, j, k

    ! note that the turbulent viscosity itself is computed in rhs_scalars.f90
    ! here we only modify the RHSs for scalars in case we're running LES

    select case (les_model)
!--------------------------------------------------------------------------------
    case(1)
       call m_les_rhss_turb_visc
!--------------------------------------------------------------------------------
    case(2)
       ! -- Dynamic Localization model with algebraic model for dissipation
       call m_les_k_src_dlm
       call m_les_k_diss_algebraic
       call m_les_rhss_turb_visc
!--------------------------------------------------------------------------------
    case(3)
       ! -- Dynamic Localization model with lag-model for dissipation
       call m_les_k_src_dlm
       call m_les_lag_model_sources
       call m_les_rhss_turb_visc
!--------------------------------------------------------------------------------
    case(4) 
       ! -- Dynamic Structure Model with algebraic model for dissipation
       ! First taking care of the passive scalars (don't have them for now)
       if (n_scalars .gt. 0) then
          write(out,*) "*** Current version of the code cannot transport scalars"
          write(out,*) "*** with the les_model=4 (Dynamic Structure Model)"
          write(out,*) "Please specify a different LES model."
          call flush(out)
          call my_exit(-1)
       end if
       ! now taking care of the LES-related scalars
       ! for the Dynamic Structure model, k-equation has turbulent viscosity
       ! so need to add that term.

       call m_les_dstm_vel_k_sources
       call m_les_k_diss_algebraic
       call m_les_rhss_turb_visc
!--------------------------------------------------------------------------------
    case(5) 
       ! -- Dynamic Structure Model with lag-model for dissipation
       ! First taking care of the passive scalars (don't have them for now)
       if (n_scalars .gt. 0) then
          write(out,*) "*** Current version of the code cannot transport scalars"
          write(out,*) "*** with the les_model=5 (Dynamic Structure Model)"
          write(out,*) "Please specify a different LES model."
          call flush(out)
          call my_exit(-1)
       end if
       ! now taking care of the LES-related scalars
       ! for the Dynamic Structure model, k-equation has turbulent viscosity
       ! so need to add that term.

       ! getting sources for velocities and production for k_s and B
       call m_les_dstm_vel_k_sources  
       ! getting |S|^2 and putting it into wrk0 (for the timescale for B)
       call m_les_k_src_dlm
       ! getting the sources and sinks for (BT) and (eps T) and a sink for k_s
       call m_les_lag_model_sources
       ! diffusing k_s, (BT) and (epsilon*T) with turbulent viscosity
       call m_les_rhss_turb_visc

!--------------------------------------------------------------------------------
    case(6) 
       ! Mixed model (Dynamic Structure model + a fraction of Dynamic Localization model)
       ! The fraction is given by the constant C_mixed, which is read from the file
       ! This is taken care of in les_get_turb_visc

       ! First taking care of the passive scalars (don't have them for now)
       if (n_scalars .gt. 0) then
          write(out,*) "*** Current version of the code cannot transport scalars"
          write(out,*) "*** with the les_model=6 (Mixed Model)"
          write(out,*) "Please specify a different LES model."
          call flush(out)
          call my_exit(-1)
       end if
       ! now taking care of the LES-related scalars
       ! for the mixed model, scalars and velocities have turbulent viscosity
       ! so need to add that term.

       ! getting DSTM sources for velocities and production for k_sgs
       call m_les_dstm_vel_k_sources  
       ! getting the DLM transfer term turb_visc*|S|^2 and adding to the RHS for k_sgs
       call m_les_k_src_dlm
       ! diffusing k_sgs with turbulent viscosity
       call m_les_rhss_turb_visc
       ! algebraic model for dissipation of k_sgs: k^{3/2}/Delta
       call m_les_k_diss_algebraic

!--------------------------------------------------------------------------------
    case(7) 
       ! Mixed model (Dynamic Structure model + a fraction of Dynamic Localization model)
       ! The fraction is given by the constant C_mixed, which is read from the file
       ! This is taken care of in les_get_turb_visc

       ! Dissipation is via lag model

       ! First taking care of the passive scalars (don't have them for now)
       if (n_scalars .gt. 0) then
          write(out,*) "*** Current version of the code cannot transport scalars"
          write(out,*) "*** with the les_model=6 (Mixed Model)"
          write(out,*) "Please specify a different LES model."
          call flush(out)
          call my_exit(-1)
       end if

       ! now taking care of the LES-related scalars
       ! for the Dynamic Structure model, k-equation has turbulent viscosity
       ! so need to add that term.

       ! getting sources for velocities and production for k_s and B
       call m_les_dstm_vel_k_sources  
       ! getting |S|^2 and putting it into wrk0 (for the timescale for B)
       call m_les_k_src_dlm
       ! getting the sources and sinks for (BT) and (eps T) and a sink for k_s
       call m_les_lag_model_sources
       ! diffusing k_s, (BT) and (epsilon*T) with turbulent viscosity
       call m_les_rhss_turb_visc

!--------------------------------------------------------------------------------
    case(10) 
       ! Mixed model (Dynamic Structure model + a fraction of Smagorinsky model)
       ! The fraction is given by the constant C_mixed, which is read from the file
       ! c_mixed.in

       ! Dissipation term in the k-equation is closed via lag-model

       ! first we add the turbulent diffusion to (BT) and (eps T)
       ! (no turbulent diffusion for the passive scalars, since they are taken
       ! care of by the Harlow model later)
       call m_les_rhss_turb_visc1(3+n_scalars+1)
       call m_les_rhss_turb_visc1(3+n_scalars+2)

       ! now calculate the SGS stress tauij and calculate the source term Pi 
       ! for (BT), the scalar with the number 3+n_scalars+1
       ! Add the source term for (BT) to the RHS for (BT)
       ! Store the SGS stress in the array tauij.
       call m_les_get_tauij_dstm_source_BT

       ! get the SGS fluxes for all passive scalars, using Harlow model
       if (n_scalars.gt.0) call m_les_rhss_sgs_flux_harlow

       ! get the lag-model sources for (BT) and (eps T) and add those to RHSs
       ! the "no_k" means that we calculate k_sgs as (BT)+(eps T), effectively
       ! reducing the number of extra equations.
       call m_les_lag_model_sources_no_k

!--------------------------------------------------------------------------------
    case default
       write(out,*) 'LES_RHS_SCALARS: invalid value of les_model:',les_model
       call flush(out)
       call my_exit(-1)
    end select

!--------------------------------------------------------------------------------
!  Making sure that we are not getting any aliasing errors. That is done by 
!  making RHS zero for wavenumbers that produce aliasing.
!--------------------------------------------------------------------------------
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 2
             if (ialias(i,j,k) .ne. 0) then
                wrk(i,j,k,4:3+n_scalars+n_les) = zip
             end if
          end do
       end do
    end do

!--------------------------------------------------------------------------------
!  outputting the LES quantities in the file les.gp
!--------------------------------------------------------------------------------
    if(iammaster) then
       if  (mod(itime,iprint1)==0) then
          open(999,file='les.gp', position='append')
          if (n_les>0) then
             if (les_model.lt.10) then
                energy = fields(1,1,1,3+n_scalars+1) / real(nxyz_all)
             else
                energy = (fields(1,1,1,3+n_scalars+1) + fields(1,1,1,3+n_scalars+2))/ real(nxyz_all)
             end if
             write(999,"(i6,x,10e15.6)") itime, time, energy, production, B, dissipation
             close(999)
          end if
       end if
       B = zip
       production = zip
       dissipation = zip
    end if
!--------------------------------------------------------------------------------
    return
  end subroutine les_rhs_scalars

!================================================================================
!================================================================================
!  calculating LES sources for velocity field and adding them to the 
!  RHS's (wrk1...3)
!================================================================================
  subroutine les_rhsv_turb_visc

    use x_fftw
    implicit none

    integer :: i, j, k, n
    real(kind=8) :: rtmp

    ! due to memory constraints we have only three work arrays wrk4..6,
    ! because the first three wrk1..3 contain already comptued velocity RHS's.

    ! Calculating S_11, S_12, S_13
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             ! S_11, du/dx
             wrk(i  ,j,k,4) = - akx(i+1) * fields(i+1,j,k,1)
             wrk(i+1,j,k,4) =   akx(i  ) * fields(i  ,j,k,1)

             ! S_12, 0.5 (du/dy + dv/dx)
             wrk(i  ,j,k,5) = -half * ( aky(k) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,2) )
             wrk(i+1,j,k,5) =  half * ( aky(k) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,2) )

             ! S_13, 0.5 (du/dz + dw/dx)
             wrk(i  ,j,k,6) = -half * ( akz(j) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,6) =  half * ( akz(j) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,3) )
          end do
       end do
    end do

    ! Multiplying them by  -2 * turbulent viscosity to get tau_11, tau_12, tau_13
    do n = 4,6
       call xFFT3d(-1,n)
       wrk(1:nx,1:ny,1:nz,n) = - two * turb_visc(1:nx,1:ny,1:nz) * wrk(1:nx,1:ny,1:nz,n) 
       call xFFT3d(1,n)
    end do

    ! Taking d/dx tau_11,  d/dy tau_12, d/dz tau_13 and subtracting from the current RHS
    ! note the sign reversal (-/+) because we subtract this from RHS
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! Cutting off any wave modes that can introduce aliasing into the velocities
             ! This "dealiasing" is done in the most restrictive manner: only Fourier modes
             ! that have ialias=0 are added
             if (ialias(i,j,k) .eq. 0) then

                wrk(i  ,j,k,1) = wrk(i  ,j,k,1) + akx(i+1)*wrk(i+1,j,k,4) + aky(k)*wrk(i+1,j,k,5) + akz(j)*wrk(i+1,j,k,6)
                wrk(i+1,j,k,1) = wrk(i+1,j,k,1) - akx(i  )*wrk(i  ,j,k,4) - aky(k)*wrk(i  ,j,k,5) - akz(j)*wrk(i  ,j,k,6)

                wrk(i  ,j,k,2) = wrk(i  ,j,k,2) + akx(i+1)*wrk(i+1,j,k,5)
                wrk(i+1,j,k,2) = wrk(i+1,j,k,2) - akx(i  )*wrk(i  ,j,k,5)

                wrk(i  ,j,k,3) = wrk(i  ,j,k,3) + akx(i+1)*wrk(i+1,j,k,6)
                wrk(i+1,j,k,3) = wrk(i+1,j,k,3) - akx(i  )*wrk(i  ,j,k,6)

             end if

          end do
       end do
    end do


    ! Calculating S_22, S_23, S_33
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! S_22, dv/dy
             wrk(i  ,j,k,4) = - aky(k) * fields(i+1,j,k,2)
             wrk(i+1,j,k,4) =   aky(k) * fields(i  ,j,k,2)

             ! S_23, 0.5 (dv/dz + dw/dy)
             wrk(i  ,j,k,5) = - half * ( akz(j) * fields(i+1,j,k,2) + aky(k) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,5) =   half * ( akz(j) * fields(i  ,j,k,2) + aky(k) * fields(i  ,j,k,3) )

             ! S_33, de/dz
             wrk(i  ,j,k,6) = - akz(j) * fields(i+1,j,k,3)
             wrk(i+1,j,k,6) =   akz(j) * fields(i  ,j,k,3)
          end do
       end do
    end do

    ! Multiplying them by -2 * turbulent viscosity to get tau_22, tau_23, tau_33
    do n = 4,6
       call xFFT3d(-1,n)
       wrk(1:nx,1:ny,1:nz,n) = - two * turb_visc(1:nx,1:ny,1:nz) * wrk(1:nx,1:ny,1:nz,n) 
       call xFFT3d(1,n)
    end do

    ! Taking
    ! d/dy tau_22, d/dz tau_23
    ! d/dy tau_23, d/dz tau_33 
    ! and subtracting from the current RHS for v and w
    ! note the sign reversal (-/+) because we subtract this from RHS
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! Cutting off any wave modes that can introduce aliasing into the velocities
             ! This "dealiasing" is done in the most restrictive manner: only Fourier modes
             ! that have ialias=0 are added
             if (ialias(i,j,k) .eq. 0) then
                wrk(i  ,j,k,2) = wrk(i  ,j,k,2) + aky(k)*wrk(i+1,j,k,4) + akz(j)*wrk(i+1,j,k,5)
                wrk(i+1,j,k,2) = wrk(i+1,j,k,2) - aky(k)*wrk(i  ,j,k,4) - akz(j)*wrk(i  ,j,k,5)

                wrk(i  ,j,k,3) = wrk(i  ,j,k,3) + aky(k)*wrk(i+1,j,k,5) + akz(j)*wrk(i+1,j,k,6)
                wrk(i+1,j,k,3) = wrk(i+1,j,k,3) - aky(k)*wrk(i  ,j,k,5) - akz(j)*wrk(i  ,j,k,6)
             end if

          end do
       end do
    end do

    return
  end subroutine les_rhsv_turb_visc

!================================================================================
!================================================================================
!    calculating LES sources for scalars and adding them to the RHS's 
!    wrk4...3+n_scalars+n_les
!
!    case when SGS stress tau_ij is modeled using turbulent viscosity
!    the turb. viscosity is supposed to be in the array turb_visc(nx,ny,nz) 
!================================================================================
  subroutine m_les_rhss_turb_visc
    implicit none
    integer :: i
    do i = 4, 3+n_scalars+n_les
       call m_les_rhss_turb_visc1(i)

    end do
    return
  end subroutine m_les_rhss_turb_visc

!================================================================================
!  Turbulent viscosity for one single field #n (or, scalar number n-3)
!================================================================================
  subroutine m_les_rhss_turb_visc1(n)

    use x_fftw
    implicit none

    integer :: i, j, k, n, tmp1, tmp2
    character :: dir

    ! have two arrays available as work arrays: wrk 3+n_scalars+n_les+1 and +2
    tmp1 = 3 + n_scalars + n_les + 1
    tmp2 = 3 + n_scalars + n_les + 2

    ! computing the second derivative, multiplying it by turb_visc
    ! and adding to the RHS (that is contained in wrk(n))

    wrk(:,:,:,tmp1) = fields(:,:,:,n)
    wrk(:,:,:,0) = zip

    directions: do i = 1,3
       if (i.eq.1) dir = 'x'
       if (i.eq.2) dir = 'y'
       if (i.eq.3) dir = 'z'

       ! taking the first derivatite, multiplying by the turb_visc
       ! then taking another derivative and adding the result to wrk0
       call x_derivative(tmp1, dir, tmp2)
       call xFFT3d(-1,tmp2)
       wrk(1:nx, 1:ny, 1:nz, tmp2) = wrk(1:nx, 1:ny, 1:nz, tmp2) * turb_visc(1:nx, 1:ny, 1:nz)
       call xFFT3d(1,tmp2)
       call x_derivative(tmp2, dir, tmp2)

       ! following Yoshizawa and Horiuti (1985), the viscosity in the scalar equation
       ! is twice the viscosity used in the production of k_sgs
       ! (Journal of Phys. Soc. Japan, V.54 N.8, pp.2834-2839)
       ! wrk(1:nx, 1:ny, 1:nz, tmp2) = two * wrk(1:nx, 1:ny, 1:nz, tmp2)

       wrk(:,:,:,0) = wrk(:,:,:,0) + wrk(:,:,:,tmp2)
    end do directions

    ! adding the d/dx ( nu_t d phi/dx) to the RHS for the scalar (field #n)
    ! only adding the Fourier modes that are not producing any aliasing
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+2
             if (ialias(i,j,k) .eq. 0) wrk(i,j,k,n) = wrk(i,j,k,n) + wrk(i,j,k,0)
          end do
       end do
    end do

    return
  end subroutine m_les_rhss_turb_visc1

!================================================================================
!================================================================================
!  Calculation of turbulent viscosity turb_visc(:,:,:)
!================================================================================
  subroutine les_get_turb_visc

    implicit none

    select case (les_model)
    case(1,4:5)
       if (mod(itime,iprint1).eq.0) write(out,*) "Smagorinsky viscosity"
       call flush(out)
       call les_get_turb_visc_smag

    case(2:3)
       call les_get_turb_visc_dlm

    case(6)
       call les_get_turb_visc_smag
       ! making turbulent viscosity a fraction of what it is since this is a
       ! mixed model
       turb_visc = C_mixed * turb_visc

    case(7)
       call les_get_turb_visc_smag
       ! making turbulent viscosity a fraction of what it is since this is a
       ! mixed model
       turb_visc = C_mixed * turb_visc

    case(10)
       ! calculating the turbulent viscosity and saving the strain in the
       ! array for the SGS stress
       call les_get_turb_visc_smag_save_sij

    case default
       write(out,*) "LES_GET_TURB_VISC: les_model: ", les_model
       write(out,*) "LES_GET_TURB_VISC: Not calculating turb_visc"
       call flush(out)
       call my_exit(-1)
    end select

    return
  end subroutine les_get_turb_visc


!================================================================================
!  Calculation of turbulent viscosity turb_visc(:,:,:) - Smagorinsky model
!================================================================================
!================================================================================
  subroutine les_get_turb_visc_smag

    use x_fftw
    implicit none

    integer :: i, j, k, n
    real*8 :: c_smag = 0.18_8

    ! due to memory constraints we have only three work arrays wrk4..6,
    ! because the first three wrk1..3 contain already comptued velocity RHS's.

    ! Calculating S_11, S_12, S_13
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             ! S_11, du/dx
             wrk(i  ,j,k,4) = - akx(i+1) * fields(i+1,j,k,1)
             wrk(i+1,j,k,4) =   akx(i  ) * fields(i  ,j,k,1)

             ! S_12, 0.5 (du/dy + dv/dx)
             wrk(i  ,j,k,5) = -half * ( aky(k) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,2) )
             wrk(i+1,j,k,5) =  half * ( aky(k) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,2) )

             ! S_13, 0.5 (du/dz + dw/dx)
             wrk(i  ,j,k,6) = -half * ( akz(j) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,6) =  half * ( akz(j) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,3) )
          end do
       end do
    end do

    ! Converting them to real space and adding to turb_visc(:,:,:)
    do n = 4,6
       call xFFT3d(-1,n)
    end do
    turb_visc(1:nx,1:ny,1:nz) =                                   wrk(1:nx,1:ny,1:nz,4)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,5)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,6)**2

    ! Calculating S_22, S_23, S_33
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! S_22, dv/dy
             wrk(i  ,j,k,4) = - aky(k) * fields(i+1,j,k,2)
             wrk(i+1,j,k,4) =   aky(k) * fields(i  ,j,k,2)

             ! S_23, 0.5 (dv/dz + dw/dy)
             wrk(i  ,j,k,5) = - half * ( akz(j) * fields(i+1,j,k,2) + aky(k) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,5) =   half * ( akz(j) * fields(i  ,j,k,2) + aky(k) * fields(i  ,j,k,3) )

             ! S_33, dw/dz
             wrk(i  ,j,k,6) = - akz(j) * fields(i+1,j,k,3)
             wrk(i+1,j,k,6) =   akz(j) * fields(i  ,j,k,3)
          end do
       end do
    end do

    ! Converting them to real space and adding to turb_visc(:,:,:)
    do n = 4,6
       call xFFT3d(-1,n)
    end do
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) +       wrk(1:nx,1:ny,1:nz,4)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,5)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) +       wrk(1:nx,1:ny,1:nz,6)**2
    ! now turb_visc contains S_ij S_ij

    ! Finishing up the turbulent viscosiy
    ! making it (C_s Delta)^2 |S|, where |S| = sqrt(2 S_{ij} S_{ij})
    turb_visc = sqrt( two * turb_visc)
    turb_visc  = turb_visc * (c_smag * les_delta)**2

    return
  end subroutine les_get_turb_visc_smag

!================================================================================
!  Calculation of turbulent viscosity turb_visc(:,:,:) - Smagorinsky model
!  Saving S_ij in the array tauij(:,:,:,:) for further needs.  This is the only
!  difference between this program and les_get_turb_visc_smag
!================================================================================
!================================================================================
  subroutine les_get_turb_visc_smag_save_sij

    use x_fftw
    implicit none

    integer :: i, j, k, n
    real*8 :: c_smag = 0.18_8

    ! due to memory constraints we have only three work arrays wrk4..6,
    ! because the first three wrk1..3 contain already comptued velocity RHS's.

    ! Calculating S_11, S_12, S_13
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             ! S_11, du/dx
             wrk(i  ,j,k,4) = - akx(i+1) * fields(i+1,j,k,1)
             wrk(i+1,j,k,4) =   akx(i  ) * fields(i  ,j,k,1)

             ! S_12, 0.5 (du/dy + dv/dx)
             wrk(i  ,j,k,5) = -half * ( aky(k) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,2) )
             wrk(i+1,j,k,5) =  half * ( aky(k) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,2) )

             ! S_13, 0.5 (du/dz + dw/dx)
             wrk(i  ,j,k,6) = -half * ( akz(j) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,6) =  half * ( akz(j) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,3) )
          end do
       end do
    end do

    ! Converting them to real space and adding to turb_visc(:,:,:)
    do n = 4,6
       call xFFT3d(-1,n)
    end do
    turb_visc(1:nx,1:ny,1:nz) =                                   wrk(1:nx,1:ny,1:nz,4)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,5)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,6)**2

    ! Saving S_11, S_12, S_13
    tauij(:,:,:,1) = wrk(:,:,:,4)
    tauij(:,:,:,2) = wrk(:,:,:,5)
    tauij(:,:,:,3) = wrk(:,:,:,6)

    ! Calculating S_22, S_23, S_33
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! S_22, dv/dy
             wrk(i  ,j,k,4) = - aky(k) * fields(i+1,j,k,2)
             wrk(i+1,j,k,4) =   aky(k) * fields(i  ,j,k,2)

             ! S_23, 0.5 (dv/dz + dw/dy)
             wrk(i  ,j,k,5) = - half * ( akz(j) * fields(i+1,j,k,2) + aky(k) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,5) =   half * ( akz(j) * fields(i  ,j,k,2) + aky(k) * fields(i  ,j,k,3) )

             ! S_33, dw/dz
             wrk(i  ,j,k,6) = - akz(j) * fields(i+1,j,k,3)
             wrk(i+1,j,k,6) =   akz(j) * fields(i  ,j,k,3)
          end do
       end do
    end do

    ! Converting them to real space and adding to turb_visc(:,:,:)
    do n = 4,6
       call xFFT3d(-1,n)
    end do
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) +       wrk(1:nx,1:ny,1:nz,4)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,5)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) +       wrk(1:nx,1:ny,1:nz,6)**2
    ! now turb_visc contains S_ij S_ij

    ! Saving S_22, S_23, S_33
    tauij(:,:,:,4) = wrk(:,:,:,4)
    tauij(:,:,:,5) = wrk(:,:,:,5)
    tauij(:,:,:,6) = wrk(:,:,:,6)
    ! now tauij(:,:,:,1:6) contains S_ij

    ! Finishing up the turbulent viscosiy
    ! making it (C_s Delta)^2 |S|, where |S| = sqrt(2 S_{ij} S_{ij})
    turb_visc = sqrt( two * turb_visc)
    turb_visc  = turb_visc * (c_smag * les_delta)**2

    return
  end subroutine les_get_turb_visc_smag_save_sij

!================================================================================
!================================================================================
!  Calculation of turbulent viscosity turb_visc(:,:,:) - DLM model
!================================================================================
  subroutine les_get_turb_visc_dlm

    use x_fftw
    implicit none
    integer :: i,j,k
    real*8  :: rkmax2, wmag2

!!$    real*8 :: C_k = 0.05d0 ! This is taken from Yoshizawa and Horiuti (1985)
    real*8 :: C_k = 0.1d0 ! This is what works for this code.  Dunno why...
    real*8 :: sctmp, sctmp1

!!$    write(out,*) "Calculating turbulent viscosity using DLM model"
!!$    call flush(out)

    wrk(:,:,:,0) = fields(:,:,:,3+n_scalars+1)

![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    rkmax2 = real(kmax**2,8)
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx+2
             wmag2 = akx(i)**2 + aky(k)**2 + akz(j)**2
             if (wmag2 .gt. rkmax2) wrk(i,j,k,0) = zip
          end do
       end do
    end do
!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

    call xFFT3d(-1,0)

    ! C alculating the minimum value of k_sgs.  If it's less that zero, clip it
    ! at zero.
    sctmp1 = minval(wrk(1:nx,1:ny,1:nz,0))
    count = 1
    call MPI_REDUCE(sctmp1,sctmp,count,MPI_REAL8,MPI_MIN,0,MPI_COMM_TASK,mpi_err)
    call MPI_BCAST(sctmp,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
    if (iammaster) write(697,"(i6,x,e15.6)") itime, sctmp
    if (iammaster .and. mod(itime,10).eq.0) call flush(697)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-

    if (sctmp.lt.zip) then
       ! write(out,*) 'LES_GET_TURB_VISC_DLM: minval of k is less than 0:',sctmp
       ! call flush(out)
       ! call my_exit(-1)
       wrk(:,:,:,0) = max(wrk(:,:,:,0), zip)       
    end if

    turb_visc(1:nx,1:ny,1:nz) = C_k * les_delta * sqrt(wrk(1:nx,1:ny,1:nz,0))

!!$![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
!!$    if (mod(itime,iwrite4).eq.0) then
!!$       tmp4 = turb_visc
!!$       write(fname,"('nut1.',i6.6)") itime
!!$       call write_tmp4
!!$
!!$       tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,0)
!!$       write(fname,"('k1.',i6.6)") itime
!!$       call write_tmp4
!!$
!!$    end if
!!$!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

    ! If the minimum of k_sgs is less than zero, we want to clip it
    ! mercilessly (already done in wrk0) then put it back in Fourier space
    if (sctmp.lt.zip) then
       call xFFT3d(1,0)
       fields(:,:,:,3+n_scalars+1) = wrk(:,:,:,0)
    end if

!!$    write(out,*) "      turbulent viscosity calculated"
!!$    call flush(out)

    return
  end subroutine les_get_turb_visc_dlm


!================================================================================
!  Subroutine that initializes k_sgs for the Dynamic Localization Model (DLM)
!================================================================================
  subroutine m_les_dlm_k_init

    use x_fftw
    implicit none

    integer :: i, j, k, n_k

    write(out,*) "m_les_dlm_k_init: initializing k_sgs"
    call flush(out)

    n_k = 3 + n_scalars + 1
    fields(:,:,:,n_k) = zip

    ! initializing it as a constant
    write(out,*) "m_les_dlm_k_init: initialized k=0.1"
    call flush(out)
    fields(:,:,:,n_k) = 0.1d0
    call xFFT3d_fields(1,n_k)
    return
  end subroutine m_les_dlm_k_init

!================================================================================
!================================================================================
!  Subroutine that calculates the source for k_sgs:   - tau_{ij} S_{ij} 
!  for the case of Dynamic Localization Model (DLM)
!
!  If called while les_model=5, then just calculate |S|^2 and put it into wrk0
!================================================================================
  subroutine m_les_k_src_dlm

    use x_fftw
    implicit none
    integer :: n1, n2, k_n, i, j, k

    ! The source itself is nothing but 2 * nu_t * S_{ij} * S_{ij}
    ! So the main hassle is to calculate S_{ij}, since nu_t is available from
    ! the array turb_visc.

    ! have two arrays available as work arrays: wrk 3+n_scalars+n_les+1 and +2
    ! also have wrk0 in which we will assemble the source at the end
    n1 = 3 + n_scalars + n_les + 1
    n2 = 3 + n_scalars + n_les + 2

    wrk(:,:,:,0) = zip
    wrk(:,:,:,n1) = zip
    wrk(:,:,:,n2) = zip

    ! calculating the S_{ij} S_{ij}.  Note that when calculating derivatives,
    ! we only process those Fourier modes that won't introduce aliasing when
    ! the quantity is squared.  These modes are given by ialias(i,j,k)=0

    ! calculating S_11, S_12
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             if (ialias(i,j,k).eq.0) then
                ! S_11, du/dx
                wrk(i  ,j,k,n1) = - akx(i+1) * fields(i+1,j,k,1)
                wrk(i+1,j,k,n1) =   akx(i  ) * fields(i  ,j,k,1)
                ! S_12, 0.5 (du/dy + dv/dx)
                wrk(i  ,j,k,n2) = -half * ( aky(k) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,2) )
                wrk(i+1,j,k,n2) =  half * ( aky(k) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,2) )
             end if
          end do
       end do
    end do
    ! converting to real space, squaring and adding to wrk0
    call xFFT3d(-1,n1);  call xFFT3d(-1,n2); 
    wrk(:,:,:,0) = wrk(:,:,:,0) + wrk(:,:,:,n1)**2 + two*wrk(:,:,:,n2)**2

    ! calculating S_13, S_22
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             if(ialias(i,j,k).eq.0) then
                ! S_13, 0.5 (du/dz + dw/dx)
                wrk(i  ,j,k,n1) = -half * ( akz(j) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,3) )
                wrk(i+1,j,k,n1) =  half * ( akz(j) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,3) )
                ! S_22, dv/dy
                wrk(i  ,j,k,n2) = - aky(k) * fields(i+1,j,k,2)
                wrk(i+1,j,k,n2) =   aky(k) * fields(i  ,j,k,2)
             end if
          end do
       end do
    end do
    ! converting to real space, squaring and adding to wrk0
    call xFFT3d(-1,n1);  call xFFT3d(-1,n2); 
    wrk(:,:,:,0) = wrk(:,:,:,0) + two*wrk(:,:,:,n1)**2 + wrk(:,:,:,n2)**2

    ! calculating S_23, S_33
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             if(ialias(i,j,k).eq.0) then
                ! S_23, 0.5 (dv/dz + dw/dy)
                wrk(i  ,j,k,n1) = - half * ( akz(j) * fields(i+1,j,k,2) + aky(k) * fields(i+1,j,k,3) )
                wrk(i+1,j,k,n1) =   half * ( akz(j) * fields(i  ,j,k,2) + aky(k) * fields(i  ,j,k,3) )
                ! S_33, dw/dz
                wrk(i  ,j,k,n2) = - akz(j) * fields(i+1,j,k,3)
                wrk(i+1,j,k,n2) =   akz(j) * fields(i  ,j,k,3)
             end if
          end do
       end do
    end do
    ! converting to real space, squaring and adding to wrk0
    call xFFT3d(-1,n1);  call xFFT3d(-1,n2); 
    wrk(:,:,:,0) = wrk(:,:,:,0) + two*wrk(:,:,:,n1)**2 + wrk(:,:,:,n2)**2

    ! at this point wrk0 contains S_{ij} S_{ij} in real space.
    ! need to multiply by two and multiply by turb_visc to get the source
    ! (the energy transfer term).  Assemble the transfer term in wrk(n1).
    ! NOTE: We do not touch wrk0 because we want to preserve S_{ij}S_{ij}
    ! for other routines.
    wrk(:,:,:,0) = two * wrk(:,:,:,0)

    ! if the subroutine is called with les_model=5, then the only part needed is the
    ! calculation of the |S|^2 in wrk0.  So now we check if les_model=5 and exit of it is
    if (les_model.eq.5) return

    ! continue calculating the transfer term
    wrk(:,:,:,n1) = wrk(:,:,:,0)
    wrk(1:nx,1:ny,1:nz,n1) = wrk(1:nx,1:ny,1:nz,n1) * turb_visc(1:nx,1:ny,1:nz)

!!$!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>DEBUG+
!!$    ! writing out the energy transfer due to the viscosity
!!$    if (mod(itime,iwrite4).eq.0) then
!!$       tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,n1)
!!$       write(fname,"('pi_nu.',i6.6) ") itime
!!$       call write_tmp4
!!$    end if
!!$!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>DEBUG-



    ! convert the transfer term to Fourier space
    call xFFT3d(1,n1)

    ! adding this energy transfer term to the RHS for k_sgs
    ! the RHS for k_sgs is supposed to be in wrk(3+n_scalars+1)
    k_n = 3 + n_scalars + 1

    ! saving the energy transfer to be output later
    if (iammaster) production = production + wrk(1,1,1,n1) / real(nxyz_all)

    ! adding the source for k_sgs
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 2
             if (ialias(i,j,k).eq.0) wrk(i,j,k,k_n) = wrk(i,j,k,k_n) + wrk(i,j,k,n1)
          end do
       end do
    end do

    ! if les_model=3 (DLM model + lag model for epsilon)
    ! if les_model=7 (DSTM+DLM model + lag model for epsilon)
    ! then add the source term to the RHS for B
    if (les_model.eq.3 .or. les_model.eq.7) then
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx + 2
                if (ialias(i,j,k).eq.0) wrk(i,j,k,k_n+1) = wrk(i,j,k,k_n+1) + wrk(i,j,k,n1)
             end do
          end do
       end do
    end if

    return
  end subroutine m_les_k_src_dlm

!================================================================================
!================================================================================
!  Dissipation term in k-equation: simple algebraic model k^(3/2)/Delta
!================================================================================
  subroutine m_les_k_diss_algebraic

    use x_fftw
    implicit none
    integer :: n_k, i, j, k
    real*8 :: sctmp, sctmp1

    ! the "field number" for k_sgs
    n_k = 3 + n_scalars + 1

    ! get the SGS kinetic energy in x-space
    wrk(:,:,:,0) = fields(:,:,:,n_k)

    ! first zero the modes that can produce aliasing 
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+2
             if (ialias(i,j,k).gt.0) wrk(i,j,k,0) = zip
          end do
       end do
    end do
    ! then convert to x-space
    call xFFT3d(-1,0)

    ! check the minimum value of k
    sctmp1 = minval(wrk(1:nx,1:ny,1:nz,0))
    count = 1
    call MPI_REDUCE(sctmp1,sctmp,count,MPI_REAL8,MPI_MIN,0,MPI_COMM_TASK,mpi_err)
    call MPI_BCAST(sctmp,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
!!$    if (sctmp.lt.zip) then
!!$       write(out,*) itime,"Minimum value of k is ", sctmp
!!$       call flush(out)
!!$    end if

    ! calculating the dissipation rate
    wrk(:,:,:,0) = max(zip, wrk(:,:,:,0))
    wrk(:,:,:,0) = wrk(:,:,:,0)**1.5d0 / les_delta
    call xFFT3d(1,0)

    ! subtracting the dissipation rate from the RHS for k
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (ialias(i,j,k).eq.0) wrk(i,j,k,n_k) = wrk(i,j,k,n_k) - wrk(i,j,k,0)
          end do
       end do
    end do
!!$    wrk(:,:,:,n_k) = wrk(:,:,:,n_k) - wrk(:,:,:,0)

    ! saving dissipation for output
    if (iammaster) dissipation = dissipation + wrk(1,1,1,0) / real(nxyz_all)


    return
  end subroutine m_les_k_diss_algebraic

!================================================================================
!================================================================================
!  Model No. 3: 
!  - Dynamic Localization model for tau_{ij} with constant coefficients
!  - Lag-model for the dissipation term (extra two equations)
!  - turbulent viscosity for all scalars and velocities = sqr(k) * Delta
!================================================================================
!================================================================================
  subroutine m_les_lag_model_sources

    use x_fftw

    implicit none
    integer :: nk, n1, n2, i, j, k


    ! the "field number" for k_sgs
    nk = 3 + n_scalars + 1
    ! the numbers for two work arrays
    n1 = 3 + n_scalars + n_les + 1
    n2 = 3 + n_scalars + n_les + 2

    ! Getting B from (B T_B).
    ! Currently T_B = 1/|S|, and |S|^2 is contained in wrk0 from m_les_k_src_dlm.
    ! - getting (B T_B) to real space
    wrk(:,:,:,n1) = fields(:,:,:,nk+1)
    call xFFT3d(-1,n1)
    ! - Dividing by T_B (multiplying by |S|)
    wrk(:,:,:,n1) = wrk(:,:,:,n1) * sqrt(wrk(:,:,:,0))
    ! - converting back to Fourier space
    call xFFT3d(1,n1)

    ! Getting epsilon from (epsilon T_epsilon)
    wrk(:,:,:,n2) = fields(:,:,:,nk+2)
    call xFFT3d(-1,n2)
    ! Currently T_epsilon = C_T Delta^(2/3) / epsilon^(1/3).
    ! Solving for epsilon:
    wrk(:,:,:,n2) = max(wrk(:,:,:,n2), zip)
    wrk(:,:,:,n2) = wrk(:,:,:,n2)**1.5D0 / (les_delta * C_T**1.5d0)
    call xFFT3d(1,n2)

    ! saving B and dissipation for output later
    if (iammaster) B = B + wrk(1,1,1,n1) / real(nxyz_all)
    if (iammaster) dissipation = dissipation + wrk(1,1,1,n2) / real(nxyz_all)


    ! Now we have B and epsilon, so we can update the RHS for k, B and epsilon
    ! with the sources.  The energy transfer term (Pi) was added to RHSs for
    ! k and B in the m_les_k_src_dlm.  Now adding the rest of the terms

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 2
             if (ialias(i,j,k) .eq. 0) then

                ! updating the RHS for k_sgs (subtracting epsilon)
                wrk(i,j,k,nk) = wrk(i,j,k,nk) - wrk(i,j,k,n2) 

                ! updating the RHS for B (adding Pi and subtracting B)
                ! note that Pi is already added in subroutines 
                ! m_les_dstm_vel_k_sources and m_les_k_src_dlm
                wrk(i,j,k,nk+1) = wrk(i,j,k,nk+1) - wrk(i,j,k,n1)

                ! updating the RHS for epsilon (adding B and subtracting epsilon)
                wrk(i,j,k,nk+2) = wrk(i,j,k,nk+2) + wrk(i,j,k,n1) - wrk(i,j,k,n2)

             end if
          end do
       end do
    end do

    return
  end subroutine m_les_lag_model_sources



!================================================================================
!================================================================================
!  Dynamic Structure Model for tau_{ij}
!================================================================================
!================================================================================

!================================================================================
!  Subroutine that calculates the LES sources for the velocitieis
!================================================================================

  subroutine m_les_dstm_vel_k_sources

    use x_fftw
    use m_filter_xfftw

    implicit none
    real*8 :: fac, rtmp1, rtmp2
    real*8 :: fs, fs1, bs, bs1
    integer :: n1, n2, n3, n4, n5
    integer :: nn, i, j, k, ii, jj, kk, nk
    character :: dir_i, dir_j
    logical :: diagonal


    ! temporary array to store the complete k-source and write it out
    real*4, allocatable :: k_source(:,:,:)
    allocate(k_source(1:nx,1:ny,1:nz))

    ! there are FIVE working arrays that we can use: wrk0 and
    ! wrk(3+n_scalars+n_les+1....+4).  The array n(:) will contain the indicies.
    ! in comments we'll refer to the arrays as wrk1...5
    n1 = 0;
    n2 = 3 + n_scalars + n_les + 1
    n3 = 3 + n_scalars + n_les + 2
    n4 = 3 + n_scalars + n_les + 3
    n5 = 3 + n_scalars + n_les + 4
    nk = 3 + n_scalars + 1

    ! converting k_sgs to x-space and placing it in wrk1 
    wrk(:,:,:,n1) = fields(:,:,:,nk)

    call xFFT3d(-1,n1)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
    ! compute the minimum value of k
    rtmp1 = minval(wrk(1:nx,:,:,n1))
    count = 1
    call MPI_REDUCE(rtmp1,rtmp2,count,MPI_REAL8,MPI_MIN,0,MPI_COMM_TASK,mpi_err)
    if (myid.eq.0 .and. mod(itime,iprint1).eq.0) then
       write(699,"(i6,x,e15.6)") itime, rtmp2
       call flush(699)
    end if
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-

    ! Assembling first part of L_ii in wrk2
    wrk(:,:,:,n3) = fields(:,:,:,1)
    wrk(:,:,:,n4) = fields(:,:,:,2)
    wrk(:,:,:,n5) = fields(:,:,:,3)
    call xFFT3d(-1,n3)
    call xFFT3d(-1,n4)
    call xFFT3d(-1,n5)
    wrk(:,:,:,n2) = wrk(:,:,:,n3)**2 + wrk(:,:,:,n4)**2 + wrk(:,:,:,n5)**2
    call xFFT3d(1,n2)
    call filter_xfftw(n2)
    call xFFT3d(-1,n2)

    ! Putting u, v, w in wrk3..5 and filtering them
    wrk(:,:,:,n3) = fields(:,:,:,1)
    wrk(:,:,:,n4) = fields(:,:,:,2)
    wrk(:,:,:,n5) = fields(:,:,:,3)
    call filter_xfftw(n3)
    call filter_xfftw(n4)
    call filter_xfftw(n5)
    call xFFT3d(-1,n3)
    call xFFT3d(-1,n4)
    call xFFT3d(-1,n5)

    ! Now subtracting the second part of L_ii into wrk2.
    wrk(:,:,:,n2) = wrk(:,:,:,n2) &
         - wrk(:,:,:,n3)**2 - wrk(:,:,:,n4)**2 - wrk(:,:,:,n5)**2


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
    if (mod(itime,iwrite4).eq.0) then
       tmp4(1:nx,:,:) = wrk(1:nx,:,:,n1)
       write(fname,"('k.',i6.6)") itime
       call write_tmp4
       tmp4(1:nx,:,:) = wrk(1:nx,:,:,n2)
       write(fname,"('L.',i6.6)") itime
       call write_tmp4
    end if
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-


    ! now put the scaling factor of the DStM in wrk1
    ! the scaling factor is 2*k/L_ii
!wrk(:,:,:,n1) = two * wrk(:,:,:,n1)  / max(wrk(:,:,:,n2),1.d-15)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG+
    ! Now we want to clip the scaling factor to be non-negative.
    ! This is done in order to make the energy transfer from the DSTM 
    ! non-negative at the places where k_sgs = 0.    
    ! Basically we shut down the DSTM energy transfer where k<=0 
    ! and let the diffusion work.  
    ! 
    ! Formixed models, the viscous part of the model should provide
    ! enough positive energy transfer to over come the nagetiveness.
    ! Although that does not happen as fast as we want.
    wrk(:,:,:,n1) = two * max(wrk(:,:,:,n1),zip)  / max(wrk(:,:,:,n2),1.d-15)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG-



    if (mod(itime,iwrite4).eq.0) then
       tmp4(1:nx,:,:) = wrk(1:nx,:,:,n1)
       write(fname,"('F.',i6.6)") itime
       call write_tmp4
    end if

    ! now going over tau_{ij} and S_{ij}, term by term.
    ! cycling over i and j.  I know this is inefficient but this results in
    ! the minimal amount of code to debud and I don't have much time to
    ! optimize the performance right now.

    ! when calculated, the LES source/sink terms for velocities are placed in
    ! the special array vel_source_les(:,:,:,:)
    vel_source_les = zip

    k_source = zip

    ! forward/backward scatter
    fs = zip
    fs1 = zip
    bs = zip
    bs1 = zip

    ! cycling over i and j
    direction_i: do i = 1,3

       if (i.eq.1) dir_i = 'x'
       if (i.eq.2) dir_i = 'y'
       if (i.eq.3) dir_i = 'z'

       direction_j: do j = 1,3

          if (j.eq.1) dir_j = 'x'
          if (j.eq.2) dir_j = 'y'
          if (j.eq.3) dir_j = 'z'

          wrk(:,:,:,n2) = fields(:,:,:,i)
          wrk(:,:,:,n3) = fields(:,:,:,j)
          wrk(:,:,:,n4) = fields(:,:,:,i)
          wrk(:,:,:,n5) = fields(:,:,:,j)
          call xFFT3d(-1,n2)
          call xFFT3d(-1,n3)
          call filter_xfftw(n4)
          call filter_xfftw(n5)
          call xFFT3d(-1,n4)
          call xFFT3d(-1,n5)

          ! hat(u_i u_j) -> wrk2
          wrk(:,:,:,n2) = wrk(:,:,:,n2) * wrk(:,:,:,n3)           
          call xFFT3d(1,n2)
          call filter_xfftw(n2)
          call xFFT3d(-1,n2)

          ! L_{ij} -> wrk2
          wrk(:,:,:,n2) = wrk(:,:,:,n2) - wrk(:,:,:,n4) * wrk(:,:,:,n5)

          ! tau_{ij} -> wrk2
          wrk(:,:,:,n2) = wrk(:,:,:,n2) * wrk(:,:,:,n1)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG +
          ! writing out tau_{ij}
          if (mod(itime,iwrite4).eq.0) then
             tmp4(1:nx,:,:) = wrk(1:nx,:,:,n2)
             write(fname,"('tau',i1,i1,'.',i6.6)") i,j,itime
             call write_tmp4
          end if
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG -


          ! The source in k-equation is really - du_i/dx_j * tau_{ij} 
          wrk(:,:,:,n3) = fields(:,:,:,i)
          call x_derivative(n3,dir_j,n3)
          call xFFT3d(-1,n3)
          wrk(:,:,:,n5) = - wrk(:,:,:,n3) * wrk(:,:,:,n2)


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG +
          k_source(:,:,:) = k_source(:,:,:) + wrk(1:nx,:,:,n5)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG -

          ! converting both tau_{ij} and souce for k_sgs to Fourier space
          call xFFT3d(1,n2)
          call xFFT3d(1,n5)

          ! calculating source for velocities u_i and u_j:
          ! for u_i : - d/dx_j tau_{ij} -> wrk3
!!$          ! for u_j : - d/dx_i tau_{ij} -> wrk4
          wrk(:,:,:,n3) = - wrk(:,:,:,n2)
!!$          wrk(:,:,:,n4) = - wrk(:,:,:,n2)
          call x_derivative(n3,dir_j,n3)
!!$          call x_derivative(n4,dir_i,n4)

          ! now adding the sources to the RHSs
          ! adding only the wavenumbers that do not produce aliasing
          adding_sources: do kk = 1,nz
             do jj = 1,ny
                do ii = 1,nx+2
                   if (ialias(ii,jj,kk).eq.0) then

                      ! velocity sources
                      vel_source_les(ii,jj,kk,i) = vel_source_les(ii,jj,kk,i) + wrk(ii,jj,kk,n3)
                      ! k_source
                      wrk(ii,jj,kk,nk) = wrk(ii,jj,kk,nk) + wrk(ii,jj,kk,n5)
                      ! B-source due to the DSTM (for models #5 and #7)
                      if (les_model.eq.5 .or. les_model.eq.7) then
                         wrk(ii,jj,kk,nk+1) = wrk(ii,jj,kk,nk+1) + wrk(ii,jj,kk,n5)
                      end if

!!$                      ! if i.ne.j that is, we need to add some more stuff
!!$                      if (j > i) then
!!$                         ! velocity sources
!!$                         vel_source_les(ii,jj,kk,j) = vel_source_les(ii,jj,kk,j) + wrk(ii,jj,kk,n4)
!!$                         ! k_source
!!$                         wrk(ii,jj,kk,nk) = wrk(ii,jj,kk,nk) + wrk(ii,jj,kk,n5)
!!$                         ! B-source (for model #5)
!!$                         if (les_model.eq.5) then
!!$                            wrk(ii,jj,kk,nk+1) = wrk(ii,jj,kk,nk+1) + wrk(ii,jj,kk,n5)
!!$                         end if
!!$                      end if

                   end if
                end do
             end do
          end do adding_sources

!!$          if (iammaster) then
!!$             write(720+3*(j-1)+i,"(i6,x,10e15.6)") itime, vel_source_les(1,1,1,:)
!!$             call flush(720+3*(j-1)+i)
!!$          end if


          ! saving production for output
          if (iammaster) production = production + wrk(1,1,1,n5) / real(nxyz_all)
!!$          if (iammaster .and. j>i) production = production + wrk(1,1,1,n5) / real(nxyz_all)

          ! doing the budget: counting the positive and negative production 
          ! (i.e., forward and backward scatter)
          call xFFT3d(-1,n5)
          do kk=1,nz
             do jj = 1,ny
                do ii = 1,nx
                   if (wrk(ii,jj,kk,n5) .gt. zip) then
                      fs1 = fs1 + wrk(ii,jj,kk,n5)
                   else
                      bs1 = bs1 + wrk(ii,jj,kk,n5)
                   end if
                end do
             end do
          end do


!!$! --------------------------------------------------
!!$        wrk(:,:,:,n5) = wrk(:,:,:,3+n_scalars+1)
!!$        call xFFT3d(-1,n5)
!!$        tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,n5)
!!$        write(fname,"('source',i1,i1)") i,j
!!$        call write_tmp4
!!$! --------------------------------------------------



       end do direction_j
    end do direction_i

    ! writing out int he file fort.698 forward scatter, back scatter produced by the DSTM
    count = 1
    call MPI_REDUCE(fs1,fs,count,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)
    call MPI_REDUCE(bs1,bs,count,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)
    if (myid.eq.0 .and. mod(itime,iprint1).eq.0) then
       fs = fs / real(nxyz_all,8)
       bs = bs / real(nxyz_all,8)
       write(698,"(i6,x,3e15.6)") itime, fs, bs
       call flush(698)
    end if

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG +
    ! writing out k_source
    if (mod(itime,iwrite4).eq.0) then
       tmp4 = k_source
       write(fname,"('pi.',i6.6)") itime
       call write_tmp4
    end if
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><> DEBUG -


    if (allocated(k_source)) deallocate(k_source)

  end subroutine m_les_dstm_vel_k_sources


!================================================================================
!================================================================================
!  Calculating the subgrid stress tauij(:,:,:,:) using mixed model: DSTM+C_m*SM.  
!  Assume that the SGS dissipation is closed using the lag-model.
!  Using the fact that S_ij is calculated prior to calling this subroutine
!  and is stored in the array tauij, we calculate the source term
!  for (BT), the scalar number 3+n_scalars+1, and add it to the RHS for (BT).
!================================================================================
  subroutine m_les_get_tauij_dstm_source_BT

    use m_filter_xfftw

    implicit none
    integer :: n1, n2, n3, n4, n5
    integer :: n, i, j, k
    integer :: n_bt, n_et

    ! there are FIVE working arrays that we can use: wrk0 and
    ! wrk(3+n_scalars+n_les+1....+4).  We use indicies n1,...,n5
    ! in comments we'll refer to the arrays as wrk1...5
    n1 = 0;
    n2 = 3 + n_scalars + n_les + 1
    n3 = 3 + n_scalars + n_les + 2
    n4 = 3 + n_scalars + n_les + 3
    n5 = 3 + n_scalars + n_les + 4

    ! the numbers for the field (BT)
    n_bt = 3 + n_scalars + 1
    n_et = 3 + n_scalars + 2

    ! L_11, L_22, L_33 -> wrk1, wrk2, wrk3
    do i = 1,3
       wrk(:,:,:,n4) = fields(:,:,:,i)
       call xFFT3d(-1,n4)
       wrk(:,:,:,n4) = wrk(:,:,:,n4)**2
       call xFFT3d(1,n4)
       call filter_xfftw(n4)
       call xFFT3d(-1,n4)

       wrk(:,:,:,n5) = fields(:,:,:,i)
       call filter_xfftw(n5)
       call xFFT3d(-1,n5)

       ! putting L_ii (no summation implied) in wrk(i)
       if (i.eq.1) n = n1
       if (i.eq.2) n = n2
       if (i.eq.3) n = n3
       wrk(:,:,:,n) = wrk(:,:,:,n4) - wrk(:,:,:,n5)**2
    end do

    ! putting the trace of L_ij, L_ii (summation implied) into wrk4
    wrk(:,:,:,n4) = wrk(:,:,:,n1) + wrk(:,:,:,n2) + wrk(:,:,:,n3)

    ! putting k_sgs into wrk5 and converting it to x-space
    ! note that k_sgs in this model is the sum (BT) + (eps T)
    wrk(:,:,:,n5) = fields(:,:,:,n_bt) + fields(:,:,:,n_et)
    call xFFT3d(-1,n5)

    ! compute the scaling factor 2*k_sgs/L_kk, and put it into wrk4
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if (wrk(i,j,k,n4) .lt. 1e-10) then
                wrk(i,j,k,n4) = zip
             else
                wrk(i,j,k,n4) = two * wrk(i,j,k,n5) / wrk(i,j,k,n4)
             end if
          end do
       end do
    end do

    ! now we want to do two things simultaneously:
    ! 1) compute tau_ij = 2k_sgs/L_kk L_ij - C_mixed * 2 * turb_visc * S_ij
    ! 2) compute the energy transfer term tau_ij S_ij
    ! the way to do this:
    ! - for each ij, compute -L_ij S_ij and add it to wrk5
    ! - after that compute the true tau_ij and put it in tauij(...)
    ! - when all indicies (ij) are processed, multiply the wrk5 by wrk4 (scaling factor)
    ! - after that add C_mixed turb_visc |S|^2 to wrk5
    ! thus the array tauij will contain tau_ij (full version, both parts of the mixed model)
    ! and wrk5 will contain -(tau_ij S_ij) = Pi, the energy transfer term

    ! computing the diagonal elements of tauij
    ! note that right now tauij contains the elements of S_ij
    ! and S_ij should be a part of tauij.  Thus the memory acrobatics.

    ! putting -(L_ij S_ij) into wrk5 (only diagonal elements so far)
    wrk(:,:,:,n5) = &
         - wrk(:,:,:,n1) * tauij(:,:,:,1) &
         - wrk(:,:,:,n2) * tauij(:,:,:,4) &
         - wrk(:,:,:,n3) * tauij(:,:,:,6) 

    ! computing the diagonal elements of tau_ij
    tauij(1:nx,:,:,1) = &
         wrk(1:nx,:,:,n4) * wrk(1:nx,:,:,n1) - &
         C_mixed * two * turb_visc(1:nx,:,:) * tauij(1:nx,:,:,1)
    tauij(1:nx,:,:,4) = &
         wrk(1:nx,:,:,n4) * wrk(1:nx,:,:,n2) - &
         C_mixed * two * turb_visc(1:nx,:,:) * tauij(1:nx,:,:,4)
    tauij(1:nx,:,:,6) = &
         wrk(1:nx,:,:,n4) * wrk(1:nx,:,:,n3) - &
         C_mixed * two * turb_visc(1:nx,:,:) * tauij(1:nx,:,:,6)

    ! currently wrk4 has the scaling factor from DSTM model: (2k_sgs/L_kk)
    ! and wrk5 has the part of the energy transfer term given by DSTM model 
    ! (only diagonal terms so far),
    ! not the whole mixed model
    ! this means that the free arrays are wrk1...3

    ! computing the off-diagonal parts of tau_ij and their contribution to the energy tansfer
    ! the index n will be the corresponding index in the tauij array, since 
    ! the elements of tauij are : tau_11, tau_12, 13, 22, 23, 33.

    ! NOTE: This is done in an inefficient way, with a good possibility of a speedup
    ! because filtering of all velocities is done twice.

    n = 0
    do i = 1,3
       do j = i,3
          n = n + 1
          if (i.lt.j) then

             ! computing L_ij
             wrk(:,:,:,n1) = fields(:,:,:,i)
             wrk(:,:,:,n2) = fields(:,:,:,j)
             call xFFT3d(-1,n1)
             call xFFT3d(-1,n2)
             wrk(:,:,:,n1) = wrk(:,:,:,n1) * wrk(:,:,:,n2)
             call xFFT3d(1,n1)
             call filter_xfftw(n1)
             call xFFT3d(-1,n1)

             wrk(:,:,:,n2) = fields(:,:,:,i)
             wrk(:,:,:,n3) = fields(:,:,:,j)
             call filter_xfftw(n2)
             call filter_xfftw(n3)
             call xFFT3d(-1,n2)
             call xFFT3d(-1,n3)

             ! calculating L_ij
             wrk(:,:,:,n1) = wrk(:,:,:,n1) - wrk(:,:,:,n2) * wrk(:,:,:,n3)

             ! adding the (-L_ij S_ij) to the transfer term in wrk5
             ! note that since these are off-diagonal terms, they enter the
             ! summation twice, thus we need to multiply them by two
             wrk(:,:,:,n5) = wrk(:,:,:,n5) - two * wrk(:,:,:,n1) * tauij(:,:,:,n)

             ! computing the full tau_ij and putting it into tauij(n)
             tauij(1:nx,:,:,n) = &
                  wrk(1:nx,:,:,n4) * wrk(1:nx,:,:,n1) - &
                  C_mixed * two * turb_visc(1:nx,:,:) * tauij(1:nx,:,:,n)

          end if
       end do
    end do

    ! now tauij(1..6) contain the full mixed model for tau_ij
    ! and wrk5 contains the sum:  -L_ij S_ij
    ! to compute the full energy transfer term we need to multiply wrk5 by
    ! the scalaing factor first and then add
    ! C_mixed * turb_visc * |S|^2.  Since the expression for turb_visc
    ! contains |S|, we can algebraically solve for |S| and add the missing
    ! part without re-calculating S_ij.

    ! multiplying wrk5 by the scaling factor
    wrk(:,:,:,n5) = wrk(:,:,:,n5) * wrk(:,:,:,n4)

    ! adding the part of Pi that comes from the viscous part of the model for tau_ij
    wrk(1:nx,:,:,n5) = wrk(1:nx,:,:,n5) + C_mixed * turb_visc(1:nx,:,:)**3 / (c_smag*les_delta)**4

    ! converting the energy transfer to F-space and adding to the RHS for (BT)
    call xFFT3d(1,n5)
    wrk(:,:,:,n_bt) = wrk(:,:,:,n_bt) + wrk(:,:,:,n5)

    ! saving the mean energy transfer <Pi> for statistics file
    if (iammaster) production = wrk(1,1,1,n5) / real(nxyz_all)

    return
  end subroutine m_les_get_tauij_dstm_source_BT


!================================================================================
!================================================================================
!  Calculate the SGS fluxes for all scalars using Harlow model.
!  The timescale is 1/|S|.  |S| is taken from the turbulent viscosity that was
!  calculated earlier and is in the array turb_visc.  The formula for turb_visc:
!  turb_visc = (c_smag * Delta)^2 |S| 
!  The formula for the SGS flux of a scalar using Harlow model is
!  flux = d/dx ( 1/|S| tau_ij d(phi)/dx_i ) 
!================================================================================
  subroutine m_les_rhss_sgs_flux_harlow

    implicit none

    integer :: n, n1, n2, n3, n4, n5
    integer :: i, j, k

    ! these five work arrays are free, so we're using them
    n1 = 0;
    n2 = 3 + n_scalars + n_les + 1
    n3 = 3 + n_scalars + n_les + 2
    n4 = 3 + n_scalars + n_les + 3
    n5 = 3 + n_scalars + n_les + 4

    ! first, get the |S| from the turb_visc array and put it in wrk5
    wrk(1:nx,:,:,n5) = turb_visc(1:nx,:,:) / (c_smag * les_delta)**2 
    ! inverting it so we have the timescale 1/|S|
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if (abs(wrk(i,j,k,n5)) .lt. 1e-10) then
                wrk(i,j,k,n5) = zip
             else
                wrk(i,j,k,n5) = 1.d0 / wrk(i,j,k,n5)
             end if
          end do
       end do
    end do
    ! now wrk5 contains the timescale for the Harlow model

    ! Calculating the SGS fluxes for scalars.
    ! The active scalars (BT) and (eps T) do not need SGS transport, because the SGS transport
    ! does not enter the equation for k_sgs.

    ! doing one scalar at a time
    sgs_fluxes_scalars: do n = 1, n_scalars

       ! get all derivatives of the scalar
       wrk(:,:,:,n1) = fields(:,:,:,3+n)
       call x_derivative(n1,'z',n3)
       call x_derivative(n1,'y',n2)
       call x_derivative(n1,'x',n1)
       ! convert them to X-space
       call xFFT3d(-1,n1)
       call xFFT3d(-1,n2)
       call xFFT3d(-1,n3)

       ! component by component, calculate the vector tau_ij d(phi)/dx_j
       ! multiply it by the timescale
       ! take the derivative d/dx_i
       ! add to the RHS of the scalar transport equation

       ! i = 1 
       wrk(:,:,:,n4) = tauij(:,:,:,1)*wrk(:,:,:,n1) + tauij(:,:,:,2)*wrk(:,:,:,n2) + tauij(:,:,:,3)*wrk(:,:,:,n3)
       wrk(:,:,:,n4) = wrk(:,:,:,n4) * wrk(:,:,:,n5)
       call xFFT3d(1,n4)
       call x_derivative(n4,'x',n4)
       ! adding to the RHS
       wrk(:,:,:,3+n) = wrk(:,:,:,3+n) + wrk(:,:,:,n4)

       ! i = 2 
       wrk(:,:,:,n4) = tauij(:,:,:,2)*wrk(:,:,:,n1) + tauij(:,:,:,4)*wrk(:,:,:,n2) + tauij(:,:,:,5)*wrk(:,:,:,n3)
       wrk(:,:,:,n4) = wrk(:,:,:,n4) * wrk(:,:,:,n5)
       call xFFT3d(1,n4)
       call x_derivative(n4,'y',n4)
       ! adding to the RHS
       wrk(:,:,:,3+n) = wrk(:,:,:,3+n) + wrk(:,:,:,n4)

       ! i = 3 
       wrk(:,:,:,n4) = tauij(:,:,:,3)*wrk(:,:,:,n1) + tauij(:,:,:,5)*wrk(:,:,:,n2) + tauij(:,:,:,6)*wrk(:,:,:,n3)
       wrk(:,:,:,n4) = wrk(:,:,:,n4) * wrk(:,:,:,n5)
       call xFFT3d(1,n4)
       call x_derivative(n4,'z',n4)
       ! adding to the RHS
       wrk(:,:,:,3+n) = wrk(:,:,:,3+n) + wrk(:,:,:,n4)

    end do sgs_fluxes_scalars

    return
  end subroutine m_les_rhss_sgs_flux_harlow
!================================================================================
!================================================================================
!  Sources for the lag-model quantities (BT) and (eps*T) that arise from the 
!  lag-model:
!  *  for (BT) the source is (-B)
!  * for (eps T) the sources are (+ B) and (- eps)
!  
!  This version caters to the model that does NOT have a separate transport
!  equation for k_sgs.  Thus the number are 
!  3+n_scalars+1 for (BT), and 3+n_scalars+2 for (eps*T) 
!================================================================================
  subroutine m_les_lag_model_sources_no_k

    use x_fftw

    implicit none
    integer :: n, n1, n2, n3, n4, n5, n_bt, n_et
    integer :: i, j, k

    ! these five work arrays are free, so we're using them
    n1 = 0;
    n2 = 3 + n_scalars + n_les + 1
    n3 = 3 + n_scalars + n_les + 2
    n4 = 3 + n_scalars + n_les + 3
    n5 = 3 + n_scalars + n_les + 4

    ! the number for BT and (eps T)
    n_bt = 3 + n_scalars + 1
    n_et = 3 + n_scalars + 2


    ! Currently T_B = 1/|S|.  We need the inverse, which is |S|.
    ! compute it from turb_visc = (C_smag * les_delta)**2 * |S|
    wrk(1:nx,:,:,n5) = turb_visc(1:nx,:,:) / (C_smag * les_delta)**2

    ! Getting B from (B T_B).
    ! - getting (B T_B) to real space
    wrk(:,:,:,n1) = fields(:,:,:,n_bt)
    call xFFT3d(-1,n1)
    ! - Dividing by T_B (multiplying by |S|)
    wrk(:,:,:,n1) = wrk(:,:,:,n1) * wrk(:,:,:,n5)
    ! - converting back to Fourier space
    call xFFT3d(1,n1)

    ! Getting epsilon from (epsilon T_epsilon)
    wrk(:,:,:,n2) = fields(:,:,:,n_et)
    call xFFT3d(-1,n2)
    ! Currently T_epsilon = C_T Delta^(2/3) / epsilon^(1/3).
    ! Solving for epsilon:
    wrk(:,:,:,n2) = max(wrk(:,:,:,n2), zip)
    wrk(:,:,:,n2) = wrk(:,:,:,n2)**1.5D0 / (les_delta * C_T**1.5d0)
    call xFFT3d(1,n2)

    ! Now we have B and epsilon, so we can update the RHS for (BT) and (eps * T)
    ! with the sources.  

    ! updating the RHS for B (subtracting B)
    wrk(:,:,:,n_bt) = wrk(:,:,:,n_bt) - wrk(:,:,:,n1)

    ! updating the RHS for (epsilon T) (adding B and subtracting epsilon)
    wrk(:,:,:,n_et) = wrk(:,:,:,n_et) + wrk(:,:,:,n1) - wrk(:,:,:,n2)

    ! saving B and dissipation for output later
    if (iammaster) B = wrk(1,1,1,n1) / real(nxyz_all)
    if (iammaster) dissipation = wrk(1,1,1,n2) / real(nxyz_all)

    return
  end subroutine m_les_lag_model_sources_no_k

!================================================================================
!================================================================================
!  Calculate the LES velocity sources from the array tauij(:,:,:,:)
!================================================================================
  subroutine les_add_vel_source_from_tauij

    implicit none
    integer :: n, n1, n2, n3, n4, n5
    integer :: i, j, k

    ! these five work arrays are free, so we're using them
    n1 = 0;
    n2 = 3 + n_scalars + n_les + 1
    n3 = 3 + n_scalars + n_les + 2
    n4 = 3 + n_scalars + n_les + 3
    n5 = 3 + n_scalars + n_les + 4

    ! First get tauij in Fourier space, all of them
    do n = 1, 6
       wrk(:,:,:,n1) = tauij(:,:,:,n)
       call xFFT3d(1,n1)
       tauij(:,:,:,n) = wrk(:,:,:,n1)
    end do

    ! Source for u
    wrk(:,:,:,n1) = tauij(:,:,:,1)
    wrk(:,:,:,n2) = tauij(:,:,:,2)
    wrk(:,:,:,n3) = tauij(:,:,:,3)
    call x_derivative(n1,'x',n1)
    call x_derivative(n2,'y',n2)
    call x_derivative(n3,'z',n3)
    wrk(:,:,:,1) = wrk(:,:,:,1) - wrk(:,:,:,n1) - wrk(:,:,:,n2) - wrk(:,:,:,n3)

    ! Source for v
    wrk(:,:,:,n1) = tauij(:,:,:,2)
    wrk(:,:,:,n2) = tauij(:,:,:,4)
    wrk(:,:,:,n3) = tauij(:,:,:,5)
    call x_derivative(n1,'x',n1)
    call x_derivative(n2,'y',n2)
    call x_derivative(n3,'z',n3)
    wrk(:,:,:,2) = wrk(:,:,:,2) - wrk(:,:,:,n1) - wrk(:,:,:,n2) - wrk(:,:,:,n3)

    ! Source for w
    wrk(:,:,:,n1) = tauij(:,:,:,3)
    wrk(:,:,:,n2) = tauij(:,:,:,5)
    wrk(:,:,:,n3) = tauij(:,:,:,6)
    call x_derivative(n1,'x',n1)
    call x_derivative(n2,'y',n2)
    call x_derivative(n3,'z',n3)
    wrk(:,:,:,3) = wrk(:,:,:,3) - wrk(:,:,:,n1) - wrk(:,:,:,n2) - wrk(:,:,:,n3)

    return
  end subroutine les_add_vel_source_from_tauij

!================================================================================
!================================================================================
end module m_les
