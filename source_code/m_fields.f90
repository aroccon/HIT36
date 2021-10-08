module m_fields

  implicit none

  real*8, allocatable :: fields(:,:,:,:)

!================================================================================
contains
!================================================================================

  subroutine m_fields_init

    use m_io
    use m_parameters
    implicit none

    integer :: n

    n = 3 + n_scalars + n_les

    allocate(fields(nx+2,ny,nz,n), stat=ierr)
    if (ierr.ne.0) then
       write(out,*) "Cannot allocate fields, stopping."
       call my_exit(-1)
    end if
    fields = zip

    write(out,"('Allocated ',i3,' fields.')") n
    call flush(out)


    return
  end subroutine m_fields_init

!================================================================================

  subroutine m_fields_exit
    use m_io
    implicit none

    if (allocated(fields)) deallocate(fields)

    write(out,*) 'fields deallocated.'
    call flush(out)

    return
  end subroutine m_fields_exit

!================================================================================


!--------------------------------------------------------------------------------
!  Subroutine that broadcasts the array "fields" from the hydro part to the 
!  "stats" part of the code
!--------------------------------------------------------------------------------
  subroutine fields_to_stats

    use m_openmpi
    use m_parameters
    use m_io
    implicit none

    integer :: n_field, n_proc, k, ratio, n_scalars_bcast

    ! broadcasting time from the hydro root process to the whole world
!!$    write(out,*) "Broadcasting fields to stats part."
!!$    call flush(out)


    ! first send it to the root process of the stats part
    count = 1
    tag = 0
    if (iammaster) then
       if (task.eq.'hydro') call MPI_SEND(TIME,count,MPI_REAL8,id_root_stats,tag,MPI_COMM_WORLD,mpi_err)
       if (task.eq.'stats') call MPI_RECV(TIME,count,MPI_REAL8,id_root_hydro,tag,MPI_COMM_WORLD,mpi_status,mpi_err)
!!$       write(out,*) "Exchanged information between master processors."
!!$       call flush(out)
    end if
    ! then broadcast it over the "stats" communicator
    if (task.eq.'stats') then
       call MPI_BCAST(TIME,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
       ! checking if we need to start advancing scalars
       if (.not. int_scalars .and. TIME .gt. TSCALAR) then
          int_scalars = .true.
          write(out,*) "Starting to move the scalars."
          call flush(out)
       end if
    end if
!!$    write(out,*) "Broadcasted to slave processors."
!!$    call flush(out)

    ! figuring out how many scalars to broadcast:
    ! 0     if we do not move scalars
    ! all   if we move scalars
    n_scalars_bcast = 0
    if (int_scalars) n_scalars_bcast = n_scalars
!!$    write(out,*) "Number of scalars to broadcast:", n_scalars_bcast
!!$    call flush(out)


    if (numprocs_hydro .ge. numprocs_stats) then

       ! using the code structure:
       ! since the size of the arrays is always 2^n, there is always 2^k slabs
       ! of a "hydro" array that correspond to q slab of the "stats" array

       ratio = numprocs_hydro / numprocs_stats

       select case (task)

       case ('hydro')

          ! sending the fields array to the corresponding process in "stats"
          do n_field = 1,3+n_scalars_bcast
             count = (nx+2)*ny*nz
             id_to = numprocs_hydro + floor(real(myid_world) / real(ratio))
             tag = myid_world*(3+n_scalars_bcast) + n_field-1

!!$             write(out,*) "Sending :", n_field, count, id_to, tag
!!$             call flush(out)

             call MPI_ISEND(fields(1,1,1,n_field),count,MPI_REAL8,id_to,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
             call MPI_WAIT(mpi_request,mpi_status,mpi_err)

!!$             write(out,*) "Sent."
!!$             call flush(out)

          end do
       case ('stats')

          ! receiving fields from hydro processors
          do n_proc = 0,ratio-1
             id_from = myid*ratio + n_proc
             count = (nx+2)*ny*nz/ratio
             k = (nz/ratio) * n_proc + 1
             do n_field = 1,3+n_scalars_bcast
                tag = id_from*(3+n_scalars_bcast) + n_field-1

!!$                write(out,*) "Receiving :", n_field, count, id_from, tag
!!$                call flush(out)

                call MPI_IRECV(fields(1,1,k,n_field),count,MPI_REAL8,&
                     id_from,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
                call MPI_WAIT(mpi_request,mpi_status,mpi_err)

!!$                write(out,*) "Received."
!!$                call flush(out)

             end do
          end do

       end select

    else

       ! now doing the same for the case when more processors are involved in 
       ! the "stat" part than in "hydro" part

       ratio = numprocs_stats / numprocs_hydro

       select case (task)

       case ('hydro')

          do n_proc = 0,ratio-1
             count = (nx+2)*ny*nz/ratio
             id_to = numprocs_hydro + myid_world*ratio + n_proc
             k = (nz/ratio) * n_proc + 1
             do n_field = 1,3+n_scalars_bcast
                tag = id_to*(3+n_scalars_bcast) + n_field-1
                call MPI_ISEND(fields(1,1,k,n_field),count,MPI_REAL8,id_to,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
                call MPI_WAIT(mpi_request,mpi_status,mpi_err)
             end do
          end do

       case ('stats')

          do n_field = 1,3+n_scalars_bcast
             count = (nx+2)*ny*nz
             id_from = floor(real(myid) / real(ratio))
             tag = myid_world*(3+n_scalars_bcast) + n_field-1
             call MPI_RECV(fields(1,1,1,n_field),count,MPI_REAL8,&
                  id_from,tag,MPI_COMM_WORLD,mpi_status,mpi_err)
             call MPI_WAIT(mpi_request,mpi_status,mpi_err)
          end do

       end select

    end if

!!$    write(out,*) "broadcasted fields to stats"
!!$    call flush(out)

    return
  end subroutine fields_to_stats



!================================================================================
!================================================================================
!--------------------------------------------------------------------------------
!  Subroutine that broadcasts the arrays that contain velocities in the x-space
!  to the "parts" part of the code, for tracking particles
!  The velocities are contained in the arrays wrk1...3
!  The whole ideology remains similar to the subroutine fields_to_stats, except
!  for the fact that the parts part of the code receives the velocities into
!  the "fields" array.
!--------------------------------------------------------------------------------
  subroutine fields_to_parts

    use m_openmpi
    use m_parameters
    use m_io
    use m_work
    implicit none

    integer :: n_field, n_proc, k, ratio, n_scalars_bcast

    ! if there are zero particles, return
    if (nptot.eq.0) return


    ! broadcasting time and timestep (dt) from the hydro root process
    ! first send it to the root process of the "parts" part
    count = 1
    if (iammaster) then
       tag = 0
       if (task.eq.'hydro') call MPI_SEND(TIME,count,MPI_REAL8,id_root_parts,tag,MPI_COMM_WORLD,mpi_err)
       if (task.eq.'parts') call MPI_RECV(TIME,count,MPI_REAL8,id_root_hydro,tag,MPI_COMM_WORLD,mpi_status,mpi_err)
       tag = 1
       if (task.eq.'hydro') call MPI_SEND(dt,count,MPI_REAL8,id_root_parts,tag,MPI_COMM_WORLD,mpi_err)
       if (task.eq.'parts') call MPI_RECV(dt,count,MPI_REAL8,id_root_hydro,tag,MPI_COMM_WORLD,mpi_status,mpi_err)
    end if
    ! then broadcast them over the "parts" communicator
    if (task.eq.'parts') call MPI_BCAST(TIME,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
    if (task.eq.'parts') call MPI_BCAST(dt  ,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)

    ! if have not yet started moving particles, return
    if (TIME.lt.starttime_particles) return

    ! if this is the first timestep when we need to start moving particles, 
    ! change the int_particle variable
    if (TIME.ge.starttime_particles .and. .not. int_particles) then
       write(out,*) 'Starting to move particles'
       call flush(out)
       int_particles = .true.
    end if

    ! figure out if we broadcast scalars (currently not)
    n_scalars_bcast = 0
    ! if (int_scalars) n_scalars_bcast = n_scalars


    if (numprocs_hydro .ge. numprocs_parts) then

       ! using the code structure:
       ! since the size of the arrays is always 2^n, there is always 2^k slabs
       ! of a "hydro" array that correspond to q slab of the "parts" array

       ratio = numprocs_hydro / numprocs_parts

       select case (task)

       case ('hydro')

          ! sending the fields array to the corresponding process in "stats"
          do n_field = 1,3+n_scalars_bcast
             count = (nx+2)*ny*nz
             id_to = id_root_parts + floor(real(myid_world) / real(ratio))
             tag = myid_world*(3+n_scalars_bcast) + n_field-1

             ! if the paricles are advected by fully resolved velocity
             ! ( that is, particles_filter_size=0) then send the fully resolved
             ! velocity to the "parts" task
!   Else, if the particles are advected by locally averaged velofity, send
!   the velocities in the Fourier form.  They will be locally averaged and
!   processed by the "parts" part of the code
             if (particles_filter_size .le. 0.d0) then
                call MPI_ISEND(wrk(1,1,1,n_field),count,MPI_REAL8,id_to,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
             else
                call MPI_ISEND(fields(1,1,1,n_field),count,MPI_REAL8,id_to,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
!!$
!!$                write(out,*) id_to, n_field, fields(:,1,1,n_field)
!!$                call flush(out)
!!$
             end if
             call MPI_WAIT(mpi_request,mpi_status,mpi_err)
          end do
       case ('parts')

          ! receiving fields from hydro processors
          do n_proc = 0,ratio-1
             id_from = myid*ratio + n_proc
             count = (nx+2)*ny*nz/ratio
             k = (nz/ratio) * n_proc + 1
             do n_field = 1,3+n_scalars_bcast
                tag = id_from*(3+n_scalars_bcast) + n_field-1

                call MPI_IRECV(fields(1,1,k,n_field),count,MPI_REAL8,&
                     id_from,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
                call MPI_WAIT(mpi_request,mpi_status,mpi_err)
!!$
!!$                write(out,*) 'rec',id_from, n_field, fields(:,1,k,n_field)
!!$                call flush(out)
!!$
             end do
          end do

       end select

    else

       ! now doing the same for the case when more processors are involved in 
       ! the "parts" part than in "hydro" part

       ratio = numprocs_parts / numprocs_hydro

       select case (task)

       case ('hydro')

          do n_proc = 0,ratio-1
             count = (nx+2)*ny*nz/ratio
             id_to = id_root_parts + myid_world*ratio + n_proc
             k = (nz/ratio) * n_proc + 1
             do n_field = 1,3+n_scalars_bcast
                tag = id_to*(3+n_scalars_bcast) + n_field-1
                if (particles_filter_size .le. 0.d0) then
                   call MPI_ISEND(wrk(1,1,k,n_field),count,MPI_REAL8,id_to,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
                else
                   call MPI_ISEND(fields(1,1,k,n_field),count,MPI_REAL8,id_to,tag,MPI_COMM_WORLD,mpi_request,mpi_err)
                end if
                call MPI_WAIT(mpi_request,mpi_status,mpi_err)
             end do
          end do

       case ('parts')

          do n_field = 1,3+n_scalars_bcast
             count = (nx+2)*ny*nz
             id_from = floor(real(myid) / real(ratio))
             tag = myid_world*(3+n_scalars_bcast) + n_field-1
             call MPI_RECV(fields(1,1,1,n_field),count,MPI_REAL8,&
                  id_from,tag,MPI_COMM_WORLD,mpi_status,mpi_err)
             call MPI_WAIT(mpi_request,mpi_status,mpi_err)
          end do

       end select

    end if

    return
  end subroutine fields_to_parts



!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================

!!$
!!$  subroutine check_bcast
!!$
!!$    use m_parameters
!!$    use m_io
!!$    use m_work
!!$    implicit none
!!$
!!$    integer :: i,j,k,n
!!$    real*8 :: zyu
!!$
!!$    ! defining the fields to be cosines (testing)
!!$    zyu = 1.0
!!$    if (task.eq.'hydro') then
!!$       do n = 1,3
!!$          do k = 1,nz
!!$             do j = 1,ny
!!$                do i = 1,nx
!!$
!!$
!!$                   fields(i,j,k,n) = sin(dble(n*i)*dx)
!!$
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end if
!!$
!!$
!!$    call m_fields_bcast
!!$
!!$
!!$    do n = 1,3
!!$
!!$       if (task.eq.'hydro') write(fname,"('hydro',i1,'.arr')") n
!!$       if (task.eq.'stats') write(fname,"('stats',i1,'.arr')") n
!!$
!!$       tmp4(:,:,:) = fields(1:nx,:,:,n)
!!$       call write_tmp4
!!$    end do
!!$
!!$    return
!!$  end subroutine check_bcast
!!$



end module m_fields


