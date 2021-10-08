program read_to_paraview

use commondata

implicit none

integer :: nstart,nstep,nend
integer :: i,j,k
logical :: check
character(len=40) :: namefile

nstart=0
nstep=200
nend=200

nx=64
ny=64
nz=64

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

dx=2*pi/(nx-1)
dy=2*pi/(ny-1)
dz=2*pi/(nz-1)

x=0.0d0
y=0.0d0
z=0.0d0

do i=1,nx-1
  x(i+1)=x(1)+dx
enddo
do j=1,ny-1
  y(j+1)=y(j)+dy
enddo
do k=1,nz-1
  z(k+1)=z(k)+dz
enddo

! read fluid data
do i=nstart,nstep,nend
 call read_fields(nstep)
enddo

deallocate(x,y,z)

end program read_to_paraview
