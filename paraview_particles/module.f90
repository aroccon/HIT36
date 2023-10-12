module commondata
 integer :: nx, ny, nz
 integer :: nstart,nend,dump
 integer*8 :: nptot
 double precision, parameter :: pi=3.14159265358979
 double precision, allocatable, dimension(:) :: pp(:,:)
end module commondata
