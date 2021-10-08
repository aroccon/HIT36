module m_force

  use m_parameters, only : kfmax, famp, force_type


  integer*4 :: n_forced_nodes, n_forced_nodes_total

  ! coordinated of the forced nodes
  integer, allocatable :: ifn(:), jfn(:), kfn(:), k_shell(:)

contains

  subroutine m_force_init

    use m_parameters
    use x_fftw

    implicit none

    integer :: i, j, k, n, n_shell

    ! if flow is not forced, return
    if (flow_type .ne. 1) return

    if (task.ne.'hydro') return

    select case (force_type)

    case (1)
       ! Machiels forcing (see article in PRL #79(18) p.3411)
       write(out,*) "Forcing #1: Machiels forcing - setting up"
       call flush(out)

       ! find out how many nodes are we forcing and book them
       n_forced_nodes = 0
       n_forced_nodes_total = 0
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
                if (n_shell .gt. 0 .and. n_shell .le. kfmax) then
                   n_forced_nodes = n_forced_nodes + 1
                end if
             end do
          end do
       end do

       ! reducing to the master process to find out the total number of forced nodes
!!$       write(out,*) 'before reducing'; call flush(out)
       count = 1
       call MPI_REDUCE(n_forced_nodes, n_forced_nodes_total, count,  &
            MPI_INTEGER4,MPI_SUM,0,MPI_COMM_TASK,mpi_err)

       ! writing out the # of forced nodes
       write(out,*) 'Number of forced nodes for this process:',n_forced_nodes
       if (myid.eq.0) write(out,*) ' total number:',n_forced_nodes_total
       call flush(out)

       ! allocating arrays for the coordinates of the forced nodes
       allocate(ifn(n_forced_nodes), jfn(n_forced_nodes), kfn(n_forced_nodes), &
            k_shell(n_forced_nodes))

       ! filling up the arrays
       n = 0
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
                if (n_shell .gt. 0 .and. n_shell .le. kfmax) then
                   n = n + 1
                   ifn(n) = i
                   jfn(n) = j
                   kfn(n) = k
                   k_shell(n) = n_shell
                end if
             end do
          end do
       end do

!!$       ! writing out the nodes
!!$       do n = 1,n_forced_nodes
!!$          write(out,"(3i4)") ifn(n),jfn(n),kfn(n)
!!$       end do
!!$       call flush(out)

    case default

       write(out,*) 'WRONG FORCE TYPE:',force_type
       write(out,*) 'STOPPING'
       call flush(out)
       stop

    end select


    return
  end subroutine m_force_init

!================================================================================

  subroutine force_velocity

    ! adding forcing to the arrays wrk(:,:,:,1:3) that already contain the RHS for velocities


    use m_openmpi
    use m_io
    use m_parameters
    use m_fields
    use m_work
    use x_fftw
    use m_stats

    implicit none

    integer :: i, j, k, n_shell, n

    real*8 :: fac, fac2

    select case (force_type)


    case (1)
       ! Machiels forcing (see article in PRL #79(18) p.3411)
       ! write(out,*) "Machiels forcing"; call flush(out)

       e_spec = zip
       e_spec1 = zip
       hits = 0
       hits1 = 0

       ! need this normalization factor because the FFT is unnormalized
       fac = one / real(nx*ny*nz_all)**2

       ! assembling the total energy in each shell
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
                if (n_shell .gt. 0 .and. n_shell .le. kfmax) then
                   fac2 = fac * (fields(i,j,k,1)**2 + fields(i,j,k,2)**2 + fields(i,j,k,3)**2)
                   if (akx(i).eq.0.d0) fac2 = fac2 * 0.5d0
                   e_spec1(n_shell) = e_spec1(n_shell) + fac2
                end if
             end do
          end do
       end do
       ! reducing the number of hits and energy to two arrays on master node
       count = kfmax
       call MPI_REDUCE(e_spec1,e_spec,count,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)

       ! getting the total energy in the region [0:kfmax] by integrating the spectrum
       if (myid.eq.0) energy = sum(e_spec(1:kfmax))

       ! broadcasting the current energy in the forcing range of wavenumbers
       count = 1
       call MPI_BCAST(energy,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)

       ! now applying the forcing to the RHS for velocities (wrk(:,:,:,1:3))

       fac = FAMP / energy

       do n = 1,n_forced_nodes
          n_shell = k_shell(n)
          i = ifn(n)
          j = jfn(n)
          k = kfn(n)

          wrk(i,j,k,1) = wrk(i,j,k,1) + fac * fields(i,j,k,1)
          wrk(i,j,k,2) = wrk(i,j,k,2) + fac * fields(i,j,k,2)
          wrk(i,j,k,3) = wrk(i,j,k,3) + fac * fields(i,j,k,3)
          
       end do
       

    case default
       write(out,*) "WRONG FORCE_TYPE:",force_type
       stop
    end select

    return
  end subroutine force_velocity


end module m_force
