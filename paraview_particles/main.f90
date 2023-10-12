program read_to_paraview
use commondata
implicit none

integer :: i,j,k
logical :: check
character(len=40) :: namefile

!! read input
open(10,file='../sc_compiled/input.f90',form='formatted')
read(10,*) nx
read(10,*)
read(10,*) nstart
read(10,*) nend
read(10,*)
read(10,*)
read(10,*) dump
do i=1,25
  read(10,*) !ignore 22 lines
end do
read(10,*) nptot
close(10)

!! overide (if there are problems)
!nstart=100
!dump=200
!nend=200
!nx=64

ny=nx
nz=nx

! read particles data
do i=nstart,nend,dump
 call read_particles(i)
enddo


end program read_to_paraview
