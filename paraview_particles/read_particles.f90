subroutine read_particles(nstep)

use commondata

integer :: nstep
character(len=40) :: namedir,namefile
character(len=6) :: numfile
character(len=3) :: setnum
logical :: check



namedir='../results/'
write(numfile,'(i6.6)') nstep


! check if u file exists; if u exists we assume that also v and w exist
namefile=trim(namedir)//'p_'//numfile//'.dat'
!write(*,*) namefile
inquire(file=trim(namefile),exist=check)


allocate(pp(3,nptot))

if(check.eqv..true.)then
write(*,*) 'Reading step ',nstep,' out of ',nend,' , particles'
  !reading particles position
  namefile=trim(namedir)//'p_'//numfile//'.dat'
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666,*)
  read(666,*) pp
  read(666,*)
  close(666,status='keep')
  ! generate paraview output file
  ! This should be updated
  call generate_output(nstep)
endif

deallocate(pp)

return
end
