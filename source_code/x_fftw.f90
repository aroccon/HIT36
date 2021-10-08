!==============================================================================!
!
!  Fast Fourier Transform that uses FFTW3 library,
!  pseudospectral DNS code
!  Copyright (C) 2006 Sergei Chumakov, Natalia Vladimirova, Misha Stepanov
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the
!    Free Software Foundation, Inc.,
!    51 Franklin Street, Fifth Floor,
!    Boston, MA 02110-1301, USA
!
!==============================================================================!

MODULE x_fftw

!==============================================================================!
!  VARIABLES
!==============================================================================!
  use m_parameters
  use m_io
  use m_fields
  use m_work
  use m_timing
  implicit none

  ! FFTW parameters that do not change
  integer(kind=8), parameter :: FFTW_ESTIMATE = 0
  integer(kind=8), parameter :: FFTW_FORWARD = -1
  integer(kind=8), parameter :: FFTW_BACKWARD = 1

  ! dedicated arrays for parallel FFT
  real(kind=8), allocatable :: xy_sheet(:, :), buff(:, :, :), z_stick(:)

  real*8, allocatable :: buff2(:,:,:)

  ! order of message passing between processors
  integer(kind=4), allocatable :: order(:), order_matrix(:, :)

  ! the arrays to store FFTW plans for the 2D r2c or c2r steps
  ! plan_r2c(1..nz, 1..n_scalars)
  integer(kind=8), allocatable :: plan_r2c(:, :), plan_c2r(:, :)
  integer(kind=8), allocatable :: plan_r2c_f(:, :), plan_c2r_f(:, :)

  ! FFTW plans for the 1D c2c forward/backward steps
  integer(kind=8) :: plan_f_c2c, plan_b_c2c

  ! k-vectors ("a" added as arrays are real
  real(kind=8), allocatable :: akx(:), aky(:), akz(:)
  real(kind=8), allocatable :: coskx2(:), cosky2(:), coskz2(:)
  real(kind=8), allocatable :: sinkx2(:), sinky2(:), sinkz2(:)
  integer(kind=4), allocatable :: rezkax(:), rezkay(:), rezkaz(:)

  ! auxiliary parameters
  integer(kind=4) :: nx21
  real(kind=8)    :: norm

  ! array that contains indicator of aliasing when products are taken
  integer(kind=1), allocatable :: ialias(:,:,:)

!==============================================================================!
!==============================================================================!
CONTAINS
!==============================================================================!
!==============================================================================!
!  SUBROUTINES
!==============================================================================!
!==============================================================================!
!  Subroutine that allocates/deallocates the FFTW arrays
!==============================================================================!
  subroutine X_FFTW_ALLOCATE(flag)

    implicit none
    integer :: flag

    if (flag == 1) then

!!$       print *,'FFTW_ALLOCATE: size(wrk,4) = ',size(wrk,4)
!!$       print *,'FFTW_ALLOCATE: bounds(wrk,4) = ',LBOUND(wrk,4),UBOUND(wrk,4)
!!$       write(out,*) 'FFTW_ALLOCATE: size(wrk,4) = ',size(wrk,4)
!!$       write(out,*) 'FFTW_ALLOCATE: bounds(wrk,4) = ',LBOUND(wrk,4),UBOUND(wrk,4)
!!$       call flush(out)

       allocate(&
            plan_r2c(nz, LBOUND(wrk,4):UBOUND(wrk,4)), &
            plan_c2r(nz, LBOUND(wrk,4):UBOUND(wrk,4)), &
            plan_r2c_f(nz, LBOUND(fields,4):UBOUND(fields,4)), &
            plan_c2r_f(nz, LBOUND(fields,4):UBOUND(fields,4)), &
            xy_sheet(nx, ny), buff(nx + 2, nz, nz), &
            z_stick(2 * nz_all), akx(nx + 2), aky(nz), akz(nz_all), & 
            rezkax(nx + 2), rezkay(nz), rezkaz(nx), &
            coskx2(nx + 2), cosky2(nz), coskz2(nx), &
            sinkx2(nx + 2), sinky2(nz), sinkz2(nx), &
            order(numprocs - 1), &
            buff2(nx+2, nz, nz), &
            ialias(nx+2, ny, nz), stat = ierr) 



!!$       write(out,*) 'Size of plan_r2c:',SIZE(plan_r2c,1),SIZE(plan_r2c,2)
!!$       write(out,*) 'Size of plan_c2r:',SIZE(plan_c2r,1),SIZE(plan_c2r,2)
!!$       call flush(out)

       if (ierr /= 0) then
          write (out, *) '*** X_FFTW_ALLOCATE: cannot allocate'
          call my_exit(-1)
       end if

       write(out,*) "x_fftw_allocated."
       call flush(out)

       ! assigning temporary values to allocated arrays
       plan_r2c = 0
       plan_c2r = 0
       xy_sheet = zip
       buff = zip
       buff2 = zip
       z_stick = zip
       order = 0
       ialias = 0



    elseif (flag == -1) then
       if (allocated(plan_r2c)) then
          deallocate(plan_r2c, plan_c2r, plan_r2c_f, plan_c2r_f, &
               xy_sheet, buff, z_stick, order, &
               akx, aky, akz, rezkax, rezkay, rezkaz, coskx2, cosky2, coskz2,&
               sinkx2, sinky2, sinkz2, buff2, ialias)
       end if
       write(out,*) "x_fftw_deallocated."
       call flush(out)
    else
       write (out, *) '*** X_FFTW_ALLOCATE: Wrong value of flag:', flag
       call my_exit(-1)
    end if

    return
  end subroutine X_FFTW_ALLOCATE

!==============================================================================!
!  Program that initializes the auxilary arrays for FFT
!==============================================================================!

  subroutine x_fftw_init

    implicit none

    integer :: itmp, ix, iy, iz, n, i, j, k
    real *8 :: rnx3

    ! write(out, *) 'Initializing FFT arrays.'
    ! call flush(out)

!------------------------------------------------------------------------------!
!  filling up the array order_matrix
!------------------------------------------------------------------------------!
    allocate(order_matrix(numprocs, numprocs), stat = ierr)
    if (ierr /= 0) then
       write(out,*) '***  X_FFTW_INIT: Cannot allocate order_matrix'
       call my_exit(-1)
    end if

    order_matrix(1, 1) = 0
    itmp = 1

    do while (itmp < numprocs)

       do ix = 1, itmp
          do iy = 1, itmp

             order_matrix(ix + itmp, iy) = order_matrix(ix, iy) + itmp
             order_matrix(ix, iy + itmp) = order_matrix(ix, iy) + itmp
             order_matrix(ix + itmp, iy + itmp) = order_matrix(ix, iy)

          end do
       end do

       itmp = 2 * itmp

    end do

!------------------------------------------------------------------------------!
!  filling the array order and deallocating order_matrix
!------------------------------------------------------------------------------!
    do ix = 1, numprocs
       if (order_matrix(ix, myid + 1) /= 0) then
          order(order_matrix(ix, myid + 1)) = ix - 1
       end if
    end do
    deallocate(order_matrix)

!------------------------------------------------------------------------------!
!  initializing FFTW plans for the 1st step in FFT --- 2D r2c
!------------------------------------------------------------------------------!
!!$    print *,'FFTW_INIT: size(wrk,4) = ',size(wrk,4)
!!$    print *,'FFTW_INIT: bounds(wrk,4) = ',LBOUND(wrk,4),UBOUND(wrk,4)

    do n = LBOUND(wrk,4),UBOUND(wrk,4)
       do iz = 1, nz
          call DFFTW_PLAN_DFT_R2C_2D(plan_r2c(iz, n), nx, ny, &
               xy_sheet, wrk(1, 1, iz, n), &
               FFTW_ESTIMATE)
       end do
    end do

    ! separately initializing the plans for FFT of the array "fields"
    do n = LBOUND(fields,4),UBOUND(fields,4)
       do iz = 1, nz
          call DFFTW_PLAN_DFT_R2C_2D(plan_r2c_f(iz, n), nx, ny, &
               xy_sheet, fields(1, 1, iz, n), &
               FFTW_ESTIMATE)
       end do
    end do

!------------------------------------------------------------------------------!
!  initializing FFTW plan for the 2nd step in FFT --- 1D c2c
!------------------------------------------------------------------------------!
    call DFFTW_PLAN_DFT_1D(plan_f_c2c, nz_all, z_stick, z_stick, &
         FFTW_FORWARD, FFTW_ESTIMATE)

!------------------------------------------------------------------------------!
!  initializing FFTW plan for the 1st step in IFFT --- 1D c2c
!------------------------------------------------------------------------------!
    call DFFTW_PLAN_DFT_1D(plan_b_c2c, nz_all, z_stick, z_stick, &
         FFTW_BACKWARD, FFTW_ESTIMATE)

!------------------------------------------------------------------------------!
!  initializing FFTW plans for the 2nd step in IFFT --- 2D c2r
!------------------------------------------------------------------------------!
    do n = LBOUND(wrk,4),UBOUND(wrk,4)
       do iz = 1, nz
          call DFFTW_PLAN_DFT_C2R_2D(plan_c2r(iz, n), nx, ny, &
               wrk(1, 1, iz, n), xy_sheet, &
               FFTW_ESTIMATE)
       end do
    end do

    ! separately intializing the plans for FFT of the array "fields"
    do n = LBOUND(fields,4),UBOUND(fields,4)
       do iz = 1, nz
          call DFFTW_PLAN_DFT_C2R_2D(plan_c2r_f(iz, n), nx, ny, &
               fields(1, 1, iz, n), xy_sheet, &
               FFTW_ESTIMATE)
       end do
    end do

!------------------------------------------------------------------------------!
!  initializing some useful constants
!------------------------------------------------------------------------------!
    nx21 = nx / 2 + 1
    norm = one / real(nx * ny * nz_all, 8)

!------------------------------------------------------------------------------!
!  filling up the wavenumber arrays akx, aky, akz
!  filling up the wavenumber arrays rezkax, rezkay, rezkaz
!------------------------------------------------------------------------------!
    ! in Fourier space it is (nx / 2 + 1) complex numbers along kx-axis
    do ix = 1, nx + 1, 2
       akx(ix) = real((ix - 1) / 2, 8)
       akx(ix + 1) = akx(ix)
       coskx2(ix) = dcos(half * akx(ix))
       sinkx2(ix) = dsin(half * akx(ix))
       coskx2(ix + 1) = coskx2(ix)
       sinkx2(ix + 1) = sinkx2(ix)
       rezkax(ix) = 0
       if (dabs(akx(ix)) > (real(nz_all, 8)) / 3.0D0) rezkax(ix) = 1
    end do
    ! in Fourier space ky-axis is distributed among the processors
    do iy = 1, nz
       aky(iy) = real(myid * nz + iy - 1, 8)
       if (aky(iy) > (0.5D0 * real(ny, 8))) aky(iy) = aky(iy) - real(ny, 8)
       cosky2(iy) = dcos(half * aky(iy))
       sinky2(iy) = dsin(half * aky(iy))
       rezkay(iy) = 0
       if (dabs(aky(iy)) > (real(ny, 8)) / 3.0D0) rezkay(iy) = 1
    end do
    ! in Fourier space the z wavenumbers are aligned along the second index
    do iz = 1, ny
       akz(iz) = real(iz - 1, 8)
       if (akz(iz) > (0.5D0 * real(nz_all, 8))) akz(iz) = akz(iz) - real(nz_all, 8)
       coskz2(iz) = dcos(half * akz(iz))
       sinkz2(iz) = dsin(half * akz(iz))
       rezkaz(iz) = 0
       if (dabs(akz(iz)) > (real(nz_all, 8)) / 3.0D0) rezkaz(iz) = 1
    end do

    ! Definition of the array ialias.
    ! The array ialias is just the number of wavenumbers at (i,j,k) that have
    ! their magnitude higher than nx/3.  This is needed in dealiasing procedures.
    rnx3 = real(nx/3, 8)
    do k = 1,nz
       if (abs(aky(k)) .gt. rnx3) ialias(:,:,k) = 1
       do j = 1,ny
          if (abs(akz(j)) .gt. rnx3) ialias(:,j,k) = ialias(:,j,k) + 1
          do i = 1,nx+2
             if (abs(akx(i)) .gt. rnx3) ialias(i,j,k) = ialias(i,j,k) + 1
          end do
       end do
    end do


    write(out,*) "x_fftw arrays are intiialized."
    call flush(out)

    return
  end subroutine x_fftw_init


!==============================================================================!
!  Subroutine that performs the FFT of a 3-D variable.  The variable is
!  contained within the array "wrk(:, :, :, n)".  Note that the
!  result of FFT has different coordinate arrangement: in physical
!  space it is (x, y, z), and in Fourier space it is (kx, kz, ky).
!  Details can be extracted from very graphic comments in the body of
!  the subroutine.
!==============================================================================!

  subroutine xxxFFT3d(flag, n)

    use m_openmpi
    implicit none

    integer :: flag, n
    integer :: ix, iy, iz

    real(kind=8) :: rtmp

    integer(kind=MPI_INTEGER_KIND) :: iproc

    if (flag == 1) then

!------------------------------------------------------------------------------!
!  Direct FFT, step 1: 2-D real-to-complex transform
!------------------------------------------------------------------------------!
!  
!   R2C           (# = ny) A y          A y                 (# = ny) A k_y
!                          |            |                            |
!          +---+---+---+---+            +            +---+---+---+---+
!         /   /   /   /   /|           /|           /   /   /   /   /|
!        /   /   /   / wrk      xy_sheet        /   /   /   / wrk
!       /   /   /   /   /  |         /  |         /   /   /   /   /  |
!      +---+---+---+---+   |        +   |        +---+---+---+---+   |
!      |   |   |   |   |   |        |   |  R2C   |   |   |   |   |   |
!      |   |   |   |   |   | -----> |   | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |   |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |   +    <---|   |   |   |   |   +
!  z   |   |   |   |   |  /         |  /     z   |   |   |   |   |  /
!      |   |   |   |   | /          | /          |   |   |   |   | /
!      |   |   |   |   |/           |/           |   |   |   |   |/
!      +---+---+---+---+            +            +---+---+---+---+
!                     /            /                            /
!           (# = nx) V x          V x         (# = nx / 2 + 1) V k_x
!
!        3   2   1   0  --- myid                   3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!

       do iz = 1, nz
          do iy = 1, ny
             do ix = 1, nx
                xy_sheet(ix, iy) = wrk(ix, iy, iz, n)
             end do
          end do

!!$          write(out,*) 'before r2c',iz,n,plan_r2c(iz, n),size(plan_r2c,1),size(plan_r2c,2)
!!$          call flush(out)

          call DFFTW_EXECUTE(plan_r2c(iz, n))
       end do
!------------------------------------------------------------------------------!
!  Direct FFT, step 2:  transposing the variable via MPI messaging
!------------------------------------------------------------------------------!
!
!   MPI           (# = ny) A k_y   +---+          +---+          (# = nz_all) A z
!                          |      /buff|         /buff|                   |
!          +---+---+---+---+     /   / +        /   / +   +---+---+---+---+
!         /   /+++/   /   /|    /   / /  MPI   /   / /   /   /   /   /   /|
!        /   /+++/   / wrk  +---+ /  ----> +---+ /   / ../.  /   / wrk
!       /   /+++/   /   /  |   |   |/         |   |/   / ../.. /   /   /  |
!      +---+---+---+---+   |   +---+          +---+   +---+---+---+---+   |
!      |   |+++|   |   |   |    A                \    |...|.  |   |   |   |
!      |   |   \____   |   ____/                  `--->+++|   |   |   |   |
!      |   |   |   |\_____/|                          |+++|   |   |   |   |
!  <---|   |   |   |   |   +                      <---|   |   |   |   |   +
!  z   |   |   |   |   |  /                       k_y |   |   |   |   |  /
!  (nz_all / numprocs)|   | /                        (# = ny / numprocs)|   | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_x                                          V k_x
!
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!
       ! - "diagonal" messages, no need to use MPI

!!$          write(out,*) 'before diagonal message'
!!$          call flush(out)


       do iz = 1, nz - 1
          do iy = iz + 1, nz
             do ix = 1, nx + 2
                rtmp = wrk(ix, myid * nz + iy, iz, n)
                wrk(ix, myid * nz + iy , iz, n) = &
                     & wrk(ix, myid * nz + iz, iy, n)
                wrk(ix, myid * nz + iz, iy, n) = rtmp
             end do
          end do
       end do
       ! - sending and receiving MPI messages

!!$          write(out,*) 'before MPI message'
!!$          call flush(out)

       count = (nx+2) * ny * nz

       do iproc = 1, numprocs - 1
          do iy = 1, nz
             do iz = 1, nz
                buff(:, iz, iy) = wrk(:, order(iproc) * nz + iy, iz, n)
             end do
          end do
          call MPI_SENDRECV_REPLACE(buff, &
               & (nx + 2) * nz * nz, MPI_REAL8, &
               & order(iproc), myid * numprocs + order(iproc), &
               & order(iproc), order(iproc) * numprocs + myid, &
               & MPI_COMM_TASK, mpi_status, mpi_err) 

          do iy = 1, nz
             do iz = 1, nz
                wrk(:, order(iproc) * nz + iz, iy, n) = buff(:, iz, iy)
             end do
          end do
       end do
!------------------------------------------------------------------------------!
!  Direct FFT, step 3: one-dimensional complex-to-complex FFT
!------------------------------------------------------------------------------!
!
!   C2C->                  A z                                            A k_z
!                          |        A z      A k_z                        |
!          +---+---+---+---+        |        |            +---+---+---+---+
!         /   /   /   /   /|        +        +           /   /   /   /   /|
!        /   /   /   / wrk       |        |          /   /   /   / wrk
!       /   /   /   /   /  |     z_stick  z_stick      /   /   /   /   /  |
!      +---+---+---+---+   |        |        |        +---+---+---+---+   |
!      |   |   |   |   |   |        |  C2C   |        |   |   |   |   |   |
!      |   |   |   |   |   | -----> | -----> | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |        |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |        |    <---|   |   |   |   |   +
!  k_y |   |   |   |   |  /         +        +    k_y |   |   |   |   |  /
!      |   |   |   |   | /                            |   |   |   |   | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_x                                          V k_x
!  
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!  
!------------------------------------------------------------------------------!

!!$          write(out,*) 'before 1d fft'
!!$          call flush(out)

       do iy = 1, nz
          do ix = 1, nx21
             do iz = 1, nz_all
                z_stick(2 * iz - 1) = wrk(2 * ix - 1, iz, iy, n)
                z_stick(2 * iz) = wrk(2 * ix, iz, iy, n)
             end do

!!$          write(out,*) 'before 1d fft',iy
!!$          call flush(out)

             call DFFTW_EXECUTE(plan_f_c2c)
             do iz = 1, nz_all
                wrk(2 * ix - 1, iz, iy, n) = z_stick(2 * iz - 1)
                wrk(2 * ix, iz, iy, n) = z_stick(2 * iz)
             end do
          end do
       end do
!------------------------------------------------------------------------------!

    elseif (flag == -1) then

!------------------------------------------------------------------------------!
!  Inverse FFT, step 1: one-dimensionsal complex-to-complex transform
!------------------------------------------------------------------------------!
!   
!   <-C2C                  A k_z                                          A z
!                          |        A k_z    A z                          |
!          +---+---+---+---+        |        |            +---+---+---+---+
!         /   /   /   /   /|        +        +           /   /   /   /   /|
!        /   /   /   / wrk       |        |          /   /   /   / wrk
!       /   /   /   /   /  |     z_stick  z_stick      /   /   /   /   /  |
!      +---+---+---+---+   |        |        |        +---+---+---+---+   |
!      |   |   |   |   |   |        |  C2C   |        |   |   |   |   |   |
!      |   |   |   |   |   | -----> | -----> | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |        |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |        |    <---|   |   |   |   |   +
!  k_y |   |   |   |   |  /         +        +    k_y |   |   |   |   |  /
!      |   |   |   |   | /                            |   |   |   |   | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_x                                          V k_x
!  
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!
       do iy = 1, nz
          do ix = 1, nx21
             do iz = 1, nz_all
                z_stick(2 * iz - 1) = wrk(2 * ix - 1, iz, iy, n)
                z_stick(2 * iz) = wrk(2 * ix, iz, iy, n)
             end do
             call DFFTW_EXECUTE(plan_b_c2c)
             do iz = 1, nz_all
                wrk(2 * ix - 1, iz, iy, n) = z_stick(2 * iz - 1)
                wrk(2 * ix, iz, iy, n) = z_stick(2 * iz)
             end do
          end do
       end do
!------------------------------------------------------------------------------!
!  Inverse FFT, step 2: transposing the variable via MPI messaging
!------------------------------------------------------------------------------!
!
!   MPI       (# = nz_all) A z     +---+          +---+          (# = ny) A k_y
!                          |      /buff|         /buff|                   |
!          +---+---+---+---+     /   / +        /   / +   +---+---+---+---+
!         /   /+++/   /   /|    /   / /  MPI   /   / /   /   /   /   /   /|
!        /   /+++/   / wrk     +---+ /  ----> +---+ /   / ../.  /   / wrk
!       /   /+++/   /   /  |   |   |/         |   |/   / ../.. /   /   /  |
!      +---+---+---+---+   |   +---+          +---+   +---+---+---+---+   |
!      |   |+++|   |   |   |    A                \    |...|.  |   |   |   |
!      |   |   \____   |   ____/                  `--->+++|   |   |   |   |
!      |   |   |   |\_____/|                          |+++|   |   |   |   |
!  <---|   |   |   |   |   +                      <---|   |   |   |   |   +
!  k_y |   |   |   |   |  /                       z   |   |   |   |   |  /
!  (# = ny / numprocs)|   | /                        (# = nz_all / numprocs)|   | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_x                                          V k_x
!
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!
       ! - "diagonal" messages, no need to use MPI
       do iy = 1, nz - 1
          do iz = iy + 1, nz
             do ix = 1, nx + 2
                rtmp = wrk(ix, myid * nz + iz, iy, n)
                wrk(ix, myid * nz + iz , iy, n) = &
                     & wrk(ix, myid * nz + iy, iz, n)
                wrk(ix, myid * nz + iy, iz, n) = rtmp
             end do
          end do
       end do
       ! - sending and receiving MPI messages
       do iproc = 1, numprocs - 1
          do iz = 1, nz
             do iy = 1, nz
                buff(:, iy, iz) = wrk(:, order(iproc) * nz + iz, iy, n)
             end do
          end do
          call MPI_SENDRECV_REPLACE(buff, &
               & (nx + 2) * nz * nz, MPI_REAL8, &
               & order(iproc), myid * numprocs + order(iproc), &
               & order(iproc), order(iproc) * numprocs + myid, &
               & MPI_COMM_TASK, mpi_status, mpi_err) 
          do iz = 1, nz
             do iy = 1, nz
                wrk(:, order(iproc) * nz + iy, iz, n) = buff(:, iy, iz)
             end do
          end do
       end do
!------------------------------------------------------------------------------!
!  Inverse FFT, step 3: 2-D complex-to-real transform
!------------------------------------------------------------------------------!
!
!   C2R           (# = ny) A k_y        A y                 (# = ny) A y
!                          |            |                            |
!          +---+---+---+---+            +            +---+---+---+---+
!         /   /   /   /   /|           /|           /   /   /   /   /|
!        /   /   /   /    wrk      xy_sheet        /   /   /   / wrk
!       /   /   /   /   /  |         /  |         /   /   /   /   /  |
!      +---+---+---+---+   |        +   |        +---+---+---+---+   |
!      |   |   |   |   |   |  C2R   |   |        |   |   |   |   |   |
!      |   |   |   |   |   | -----> |   | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |   |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |   +    <---|   |   |   |   |   +
!  z   |   |   |   |   |  /         |  /     z   |   |   |   |   |  /
!      |   |   |   |   | /          | /          |   |   |   |   | /
!      |   |   |   |   |/           |/           |   |   |   |   |/
!      +---+---+---+---+            +            +---+---+---+---+
!                     /            /                            /
!   (# = nx / 2 + 1) V k_x        V x                 (# = nx) V x
!
!        3   2   1   0  --- myid                   3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!
       do iz = 1, nz
          call DFFTW_EXECUTE(plan_c2r(iz, n))
          do iy = 1, ny
             do ix = 1, nx
                wrk(ix, iy, iz, n) = norm * xy_sheet(ix, iy)
             end do
          end do
       end do
!------------------------------------------------------------------------------!

    end if

    return
  end subroutine xxxFFT3d


!==============================================================================!
!  Subroutine that calculates the derivative.
!
!  Takes variable from wrk(:,:,:,n), differentiates it and put into
!  wrk(:,:,:,nto).  All happens in Fourier space
!==============================================================================!

  subroutine x_derivative(n,axis,nto)

    implicit none
    integer :: ix, iy, iz, n, nto
    real(kind=8) :: rtmp
    character :: axis

    ! in Fourier space it is (kx, kz, ky)
    ! here we multiply by k-vector, will multiply by i later
    select case (axis)
    case ('x')
       do iy = 1, nz
          do iz = 1, nz_all
             do ix = 1, nx + 2
                wrk(ix, iz, iy, nto) = akx(ix) * wrk(ix, iz, iy, n)
             end do
          end do
       end do
    case ('y')
       do iy = 1, nz
          wrk(:, :, iy, nto) = aky(iy) * wrk(:, :, iy, n)
       end do
    case ('z')
       do iy = 1, nz
          do iz = 1, nz_all
             wrk(:, iz, iy, nto) = akz(iz) * wrk(:, iz, iy, n)
          end do
       end do
    case default
       write (out, *) '*** x_derivative: wrong value of axis: ', axis
       call my_exit(-1)
    end select

    ! multiplying wrk(ix) + i wrk(ix + 1) by i
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx + 1, 2
             rtmp = -wrk(ix + 1, iz, iy, nto)
             wrk(ix + 1, iz, iy, nto) = wrk(ix, iz, iy, nto)
             wrk(ix, iz, iy, nto) = rtmp
          end do
       end do
    end do


  end subroutine x_derivative


!==============================================================================!
!  Subroutine --- dealiasing.
!==============================================================================!

  subroutine x_dealiasing(n1, n2)

    implicit none

    integer :: n1, n2, tmp1, tmp2, tmp
    integer :: ix, iy, iz

    tmp = 0
    tmp1 = 1
    tmp2 = 2

    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx + 2
             if (rezkax(ix) + rezkay(iy) + rezkaz(iz) > 1) then
                wrk(ix, iz, iy, n1) = zip
                wrk(ix, iz, iy, n2) = zip
             end if
             wrk(ix, iz, iy, tmp1) = wrk(ix, iz, iy, n1)
             wrk(ix, iz, iy, tmp2) = wrk(ix, iz, iy, n2)
          end do
       end do
    end do

    call xFFT3d(-1, tmp1)
    call xFFT3d(-1, tmp2)

    do iz = 1, nz
       do iy = 1, ny
          do ix = 1, nx
             wrk(ix, iy, iz, tmp) = wrk(ix, iy, iz, tmp1) * &
                  & wrk(ix, iy, iz, tmp2)
          end do
       end do
    end do

!   phase-shift, x direction
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx21, 2
             wrk(ix, iz, iy, tmp1) = coskx2(ix) * &
                  wrk(ix, iz, iy, n1) - &
                  & sinkx2(ix) * &
                  wrk(ix + 1, iz, iy, n1)
             wrk(ix + 1, iz, iy, tmp1) = coskx2(ix) * &
                  wrk(ix + 1, iz, iy, n1) + &
                  & sinkx2(ix) * &
                  wrk(ix, iz, iy, n1)
             wrk(ix, iz, iy, tmp2) = coskx2(ix) * &
                  wrk(ix, iz, iy, n2) - &
                  & sinkx2(ix) * &
                  wrk(ix + 1, iz, iy, n2)
             wrk(ix + 1, iz, iy, tmp2) = coskx2(ix) * &
                  wrk(ix + 1, iz, iy, n2) + &
                  & sinkx2(ix) * &
                  wrk(ix, iz, iy, n2)
          end do
       end do
    end do

    call xFFT3d(-1, tmp1)
    call xFFT3d(-1, tmp2)

    do iz = 1, nz
       do iy = 1, ny
          do ix = 1, nx
             wrk(ix, iy, iz, tmp) = wrk(ix, iy, iz, tmp) + &
                  wrk(ix, iy, iz, tmp1) * &
                  & wrk(ix, iy, iz, tmp2)
          end do
       end do
    end do

!   phase-shift, y direction
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx21, 2
             wrk(ix, iz, iy, tmp1) = cosky2(iy) * &
                  wrk(ix, iz, iy, n1) - &
                  & sinky2(iy) * &
                  wrk(ix + 1, iz, iy, n1)
             wrk(ix + 1, iz, iy, tmp1) = cosky2(iy) * &
                  wrk(ix + 1, iz, iy, n1) + &
                  & sinky2(iy) * &
                  wrk(ix, iz, iy, n1)
             wrk(ix, iz, iy, tmp2) = cosky2(iy) * &
                  wrk(ix, iz, iy, n2) - &
                  & sinky2(iy) * &
                  wrk(ix + 1, iz, iy, n2)
             wrk(ix + 1, iz, iy, tmp2) = cosky2(iy) * &
                  wrk(ix + 1, iz, iy, n2) + &
                  & sinky2(iy) * &
                  wrk(ix, iz, iy, n2)
          end do
       end do
    end do

    call xFFT3d(-1, tmp1)
    call xFFT3d(-1, tmp2)

    do iz = 1, nz
       do iy = 1, ny
          do ix = 1, nx
             wrk(ix, iy, iz, tmp) = wrk(ix, iy, iz, tmp) + &
                  wrk(ix, iy, iz, tmp1) * &
                  & wrk(ix, iy, iz, tmp2)
          end do
       end do
    end do

!   phase-shift, z direction
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx21, 2
             wrk(ix, iz, iy, tmp1) = coskz2(iz) * &
                  wrk(ix, iz, iy, n1) - &
                  & sinkz2(iz) * &
                  wrk(ix + 1, iz, iy, n1)
             wrk(ix + 1, iz, iy, tmp1) = coskz2(iz) * &
                  wrk(ix + 1, iz, iy, n1) + &
                  & sinkz2(iz) * &
                  wrk(ix, iz, iy, n1)
             wrk(ix, iz, iy, tmp2) = coskz2(iz) * &
                  wrk(ix, iz, iy, n2) - &
                  & sinkz2(iz) * &
                  wrk(ix + 1, iz, iy, n2)
             wrk(ix + 1, iz, iy, tmp2) = coskz2(iz) * &
                  wrk(ix + 1, iz, iy, n2) + &
                  & sinkz2(iz) * &
                  wrk(ix, iz, iy, n2)
          end do
       end do
    end do

    call xFFT3d(-1, tmp1)
    call xFFT3d(-1, tmp2)

    do iz = 1, nz
       do iy = 1, ny
          do ix = 1, nx
             wrk(ix, iy, iz, tmp) = 0.25D0 * (wrk(ix, iy, iz, tmp) + &
                  wrk(ix, iy, iz, tmp1) * &
                  & wrk(ix, iy, iz, tmp2))
          end do
       end do
    end do

    call xFFT3d(1, tmp)

  end subroutine x_dealiasing

!==============================================================================!
!==============================================================================!
!==============================================================================!
!==============================================================================!
!  Subroutine that performs the FFT of a 3-D variable.  The variable is
!  contained within the array "wrk(:, :, :, n)".  Note that the
!  result of FFT has different coordinate arrangement: in physical
!  space it is (x, y, z), and in Fourier space it is (kx, kz, ky).
!  Details can be extracted from very graphic comments in the body of
!  the subroutine.
!==============================================================================!

  subroutine xFFT3d(flag, n)

    use m_openmpi
    use m_io
    implicit none

    integer :: flag, n
    integer :: ix, iy, iz, i, j, k

    INTEGER (KIND=MPI_INTEGER_KIND) :: tag1, tag2, iproc

    real(kind=8) :: rtmp

    if (flag == 1) then

!------------------------------------------------------------------------------!
!  Direct FFT, step 1: 2-D real-to-complex transform
!------------------------------------------------------------------------------!

       if(benchmarking) then
          call system_clock(i81,dcpu)
          bm(11) = bm(11) - i81
       end if

       do iz = 1, nz
          do iy = 1, ny
             do ix = 1, nx
                xy_sheet(ix, iy) = wrk(ix, iy, iz, n)
             end do
          end do
          call DFFTW_EXECUTE(plan_r2c(iz, n))
       end do

       if (benchmarking) then
          call system_clock(i82,dcpu)
          bm(1) = bm(1) + i82 - i81
          i81 = i82
       end if

!------------------------------------------------------------------------------!
!  Direct FFT, step 2:  transposing the variable via MPI messaging
!------------------------------------------------------------------------------!
       ! - "diagonal" messages, no need to use MPI

       do iz = 1, nz - 1
          do iy = iz + 1, nz
             do ix = 1, nx + 2
                rtmp = wrk(ix, myid * nz + iy, iz, n)
                wrk(ix, myid * nz + iy , iz, n) = &
                     & wrk(ix, myid * nz + iz, iy, n)
                wrk(ix, myid * nz + iz, iy, n) = rtmp
             end do
          end do
       end do

       ! - sending and receiving MPI messages

       count = (nx+2) * nz * nz

       do iproc = 1, numprocs - 1

!!$          ! The following order  should be implemented in the future, 
!!$          ! because it shaves off about 4% of time
!!$          id_to   = mod(myid+iproc,numprocs)
!!$          id_from = mod(myid-iproc+numprocs,numprocs)

          id_to   = order(iproc)
          id_from = id_to

          tag1 = myid*numprocs + iproc
          tag2 = id_from*numprocs + iproc

          do iy = 1, nz
             do iz = 1, nz
                buff(:, iz, iy) = wrk(:, id_to * nz + iy, iz, n)
             end do
          end do

          call MPI_SENDRECV(&
               buff,  count, MPI_REAL8, id_to,   tag1, &
               buff2, count, MPI_REAL8, id_from, tag2, &
               MPI_COMM_TASK, mpi_status, mpi_err)

          do iy = 1, nz
             do iz = 1, nz
                wrk(:, id_from * nz + iz, iy, n) = buff2(:, iz, iy)
             end do
          end do

       end do

       if (benchmarking) then
          call system_clock(i82,dcpu)
          bm(2) = bm(2) + i82 - i81
          i81 = i82
       end if

!------------------------------------------------------------------------------!
!  Direct FFT, step 3: one-dimensional complex-to-complex FFT
!------------------------------------------------------------------------------!

       do iy = 1, nz
          do ix = 1, nx21
             do iz = 1, nz_all
                z_stick(2 * iz - 1) = wrk(2 * ix - 1, iz, iy, n)
                z_stick(2 * iz) = wrk(2 * ix, iz, iy, n)
             end do

             call DFFTW_EXECUTE(plan_f_c2c)
             do iz = 1, nz_all
                wrk(2 * ix - 1, iz, iy, n) = z_stick(2 * iz - 1)
                wrk(2 * ix, iz, iy, n) = z_stick(2 * iz)
             end do
          end do
       end do

       if (benchmarking) then
          call system_clock(i82,dcpu)
          bm(3) = bm(3) + i82 - i81
          bm(11) = bm(11) + i82
       end if

!------------------------------------------------------------------------------!

    elseif (flag == -1) then

!------------------------------------------------------------------------------!
!  Inverse FFT, step 1: one-dimensionsal complex-to-complex transform
!------------------------------------------------------------------------------!
       if(benchmarking) then
          call system_clock(i81,dcpu)
          bm(12) = bm(12) - i81
       end if


       do iy = 1, nz
          do ix = 1, nx21
             do iz = 1, nz_all
                z_stick(2 * iz - 1) = wrk(2 * ix - 1, iz, iy, n)
                z_stick(2 * iz) = wrk(2 * ix, iz, iy, n)
             end do
             call DFFTW_EXECUTE(plan_b_c2c)
             do iz = 1, nz_all
                wrk(2 * ix - 1, iz, iy, n) = z_stick(2 * iz - 1)
                wrk(2 * ix, iz, iy, n) = z_stick(2 * iz)
             end do
          end do
       end do

       if (benchmarking) then
          call system_clock(i82,dcpu)
          bm(4) = bm(4) + i82 - i81
          i81 = i82
       end if
!------------------------------------------------------------------------------!
!  Inverse FFT, step 2: transposing the variable via MPI messaging
!------------------------------------------------------------------------------!
       ! - "diagonal" messages, no need to use MPI
       do iy = 1, nz - 1
          do iz = iy + 1, nz
             do ix = 1, nx + 2
                rtmp = wrk(ix, myid * nz + iz, iy, n)
                wrk(ix, myid * nz + iz , iy, n) = &
                     & wrk(ix, myid * nz + iy, iz, n)
                wrk(ix, myid * nz + iy, iz, n) = rtmp
             end do
          end do
       end do


       ! - sending and receiving MPI messages

       count = (nx+2) * nz * nz

       do iproc = 1, numprocs - 1

!!$          ! This order should be implemented in the future because this saves
!!$          ! about 4% of wallclock time for FFT
!!$          id_to   = mod(myid+iproc,numprocs)
!!$          id_from = mod(myid-iproc+numprocs,numprocs)

          id_to   = order(iproc)
          id_from = order(iproc)

          tag1 = myid*numprocs + iproc
          tag2 = id_from*numprocs + iproc

          do iz = 1, nz
             do iy = 1, nz
                buff(:, iy, iz) = wrk(:, id_to * nz + iz, iy, n)
             end do
          end do

          call MPI_SENDRECV(&
               buff,  count, MPI_REAL8, id_to,   tag1, &
               buff2, count, MPI_REAL8, id_from, tag2, &
               MPI_COMM_TASK, mpi_status, mpi_err)

          do iz = 1, nz
             do iy = 1, nz
                wrk(:, id_from * nz + iy, iz, n) = buff2(:, iy, iz)
             end do
          end do
       end do

       if (benchmarking) then
          call system_clock(i82,dcpu)
          bm(5) = bm(5) + i82 - i81
          i81 = i82
       end if

!------------------------------------------------------------------------------!
!  Inverse FFT, step 3: 2-D complex-to-real transform
!------------------------------------------------------------------------------!
       do iz = 1, nz
          call DFFTW_EXECUTE(plan_c2r(iz, n))
          do iy = 1, ny
             do ix = 1, nx
                wrk(ix, iy, iz, n) = norm * xy_sheet(ix, iy)
             end do
          end do
       end do

       if (benchmarking) then
          call system_clock(i82,dcpu)
          bm(6) = bm(6) + i82 - i81
          bm(12) = bm(12) + i82
       end if
!------------------------------------------------------------------------------!

    end if

    return
  end subroutine xFFT3d

!================================================================================



!==============================================================================!
!  Subroutine that performs the FFT of a 3-D variable.  The variable is
!  contained within the array "fields(:, :, :, n)".  Note that the
!  result of FFT has different coordinate arrangement: in physical
!  space it is (x, y, z), and in Fourier space it is (kx, kz, ky).
!==============================================================================!

  subroutine xFFT3d_fields(flag, n)

    use m_openmpi
    use m_io
    use m_fields
    implicit none

    integer :: flag, n
    integer :: ix, iy, iz, i, j, k

    INTEGER (KIND=MPI_INTEGER_KIND) :: tag1, tag2, iproc

    real(kind=8) :: rtmp

    if (flag == 1) then

!------------------------------------------------------------------------------!
!  Direct FFT, step 1: 2-D real-to-complex transform
!------------------------------------------------------------------------------!
       do iz = 1, nz
          do iy = 1, ny
             do ix = 1, nx
                xy_sheet(ix, iy) = fields(ix, iy, iz, n)
             end do
          end do
          call DFFTW_EXECUTE(plan_r2c_f(iz, n))
       end do

!------------------------------------------------------------------------------!
!  Direct FFT, step 2:  transposing the variable via MPI messaging
!------------------------------------------------------------------------------!
       ! - "diagonal" messages, no need to use MPI

       do iz = 1, nz - 1
          do iy = iz + 1, nz
             do ix = 1, nx + 2
                rtmp = fields(ix, myid * nz + iy, iz, n)
                fields(ix, myid * nz + iy , iz, n) = &
                     & fields(ix, myid * nz + iz, iy, n)
                fields(ix, myid * nz + iz, iy, n) = rtmp
             end do
          end do
       end do

       ! - sending and receiving MPI messages

       count = (nx+2) * nz * nz

       do iproc = 1, numprocs - 1

!!$          ! The following order  should be implemented in the future, 
!!$          ! because it shaves off about 4% of time
!!$          id_to   = mod(myid+iproc,numprocs)
!!$          id_from = mod(myid-iproc+numprocs,numprocs)

          id_to   = order(iproc)
          id_from = id_to

          tag1 = myid*numprocs + iproc
          tag2 = id_from*numprocs + iproc

          do iy = 1, nz
             do iz = 1, nz
                buff(:, iz, iy) = fields(:, id_to * nz + iy, iz, n)
             end do
          end do

          call MPI_SENDRECV(&
               buff,  count, MPI_REAL8, id_to,   tag1, &
               buff2, count, MPI_REAL8, id_from, tag2, &
               MPI_COMM_TASK, mpi_status, mpi_err)

          do iy = 1, nz
             do iz = 1, nz
                fields(:, id_from * nz + iz, iy, n) = buff2(:, iz, iy)
             end do
          end do

       end do

!------------------------------------------------------------------------------!
!  Direct FFT, step 3: one-dimensional complex-to-complex FFT
!------------------------------------------------------------------------------!

       do iy = 1, nz
          do ix = 1, nx21
             do iz = 1, nz_all
                z_stick(2 * iz - 1) = fields(2 * ix - 1, iz, iy, n)
                z_stick(2 * iz    ) = fields(2 * ix    , iz, iy, n)
             end do

             call DFFTW_EXECUTE(plan_f_c2c)
             do iz = 1, nz_all
                fields(2 * ix - 1, iz, iy, n) = z_stick(2 * iz - 1)
                fields(2 * ix    , iz, iy, n) = z_stick(2 * iz)
             end do
          end do
       end do
!------------------------------------------------------------------------------!

    elseif (flag == -1) then

!------------------------------------------------------------------------------!
!  Inverse FFT, step 1: one-dimensionsal complex-to-complex transform
!------------------------------------------------------------------------------!
       do iy = 1, nz
          do ix = 1, nx21
             do iz = 1, nz_all
                z_stick(2 * iz - 1) = fields(2 * ix - 1, iz, iy, n)
                z_stick(2 * iz    ) = fields(2 * ix    , iz, iy, n)
             end do
             call DFFTW_EXECUTE(plan_b_c2c)
             do iz = 1, nz_all
                fields(2 * ix - 1, iz, iy, n) = z_stick(2 * iz - 1)
                fields(2 * ix    , iz, iy, n) = z_stick(2 * iz)
             end do
          end do
       end do
!------------------------------------------------------------------------------!
!  Inverse FFT, step 2: transposing the variable via MPI messaging
!------------------------------------------------------------------------------!
       ! - "diagonal" messages, no need to use MPI
       do iy = 1, nz - 1
          do iz = iy + 1, nz
             do ix = 1, nx + 2
                rtmp = fields(ix, myid * nz + iz, iy, n)
                fields(ix, myid * nz + iz , iy, n) = &
                     & fields(ix, myid * nz + iy, iz, n)
                fields(ix, myid * nz + iy, iz, n) = rtmp
             end do
          end do
       end do


       ! - sending and receiving MPI messages

       count = (nx+2) * nz * nz

       do iproc = 1, numprocs - 1

!!$          ! This order should be implemented in the future because this saves
!!$          ! about 4% of wallclock time for FFT
!!$          id_to   = mod(myid+iproc,numprocs)
!!$          id_from = mod(myid-iproc+numprocs,numprocs)

          id_to   = order(iproc)
          id_from = order(iproc)

          tag1 = myid*numprocs + iproc
          tag2 = id_from*numprocs + iproc

          do iz = 1, nz
             do iy = 1, nz
                buff(:, iy, iz) = fields(:, id_to * nz + iz, iy, n)
             end do
          end do

          call MPI_SENDRECV(&
               buff,  count, MPI_REAL8, id_to,   tag1, &
               buff2, count, MPI_REAL8, id_from, tag2, &
               MPI_COMM_TASK, mpi_status, mpi_err)

          do iz = 1, nz
             do iy = 1, nz
                fields(:, id_from * nz + iy, iz, n) = buff2(:, iy, iz)
             end do
          end do
       end do
!------------------------------------------------------------------------------!
!  Inverse FFT, step 3: 2-D complex-to-real transform
!------------------------------------------------------------------------------!
       do iz = 1, nz
          call DFFTW_EXECUTE(plan_c2r_f(iz, n))
          do iy = 1, ny
             do ix = 1, nx
                fields(ix, iy, iz, n) = norm * xy_sheet(ix, iy)
             end do
          end do
       end do
!------------------------------------------------------------------------------!

    end if

    return
  end subroutine xFFT3d_fields

!================================================================================

end MODULE X_FFTW
