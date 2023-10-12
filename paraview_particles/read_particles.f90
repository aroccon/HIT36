subroutine read_particles(nstep)

use commondata

integer :: nstep, idp
integer :: dummy
character(len=40) :: namedir,namefile
character(len=6) :: numfile
character(len=3) :: setnum
logical :: check



namedir='../results/'
write(numfile,'(i6.6)') nstep


namefile=trim(namedir)//'p_'//numfile//'.dat'
write(*,*) namefile
inquire(file=trim(namefile),exist=check)

write(*,*) "Reading"

allocate(pp(3,nptot))

if(check.eqv..true.)then
write(*,*) 'Reading step ',nstep,' out of ',nend,' , particles'
  !reading particles position
  namefile=trim(namedir)//'p_'//numfile//'.dat'
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) dummy
  !write(*,*) "dummy", dummy !debug only
  do i=1,nptot
     !read id particle + position
     read(666) idp, pp(1:3,i)
  enddo
  close(666,status='keep')

  ! generate paraview output file
  call generate_output(nstep)
endif

deallocate(pp)

return
end
