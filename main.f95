!### MolecularDynamic November 2019
!### include
include 'MyTools.f95'
!###

!####
!#Program start
!####

program main
!### use of include
use WConst
use MyVector
use Data
use LenModule
use AlgRevEiler
use AlgPbcRevEiler
use AlgVerlet

!###
implicit none

!### Global
integer, parameter :: dimen = 8
real , parameter :: startDistance = 1.5
type(atomsData) :: mainAtomsData
type(atomsData) :: pdcAtomsData
!###

integer :: i,j,k,t


mainAtomsData%aCount = dimen * dimen * dimen
pdcAtomsData%aCount = dimen * dimen * dimen
allocate(mainAtomsData%posAtoms(mainAtomsData%aCount))
allocate(mainAtomsData%velAtoms(mainAtomsData%aCount))
allocate(mainAtomsData%forceAtoms(mainAtomsData%aCount))

allocate(pdcAtomsData%posAtoms(9*pdcAtomsData%aCount))
allocate(pdcAtomsData%velAtoms(pdcAtomsData%aCount))
allocate(pdcAtomsData%forceAtoms(pdcAtomsData%aCount))


do i = 1, dimen
  do j = 1, dimen
    do k = 1, dimen
    mainAtomsData%velAtoms( (i-1)*dimen**2 + (j-1)*dimen + k ) = mainAtomsData%velAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )*0.0
    mainAtomsData%posAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )%x = (i-1)*startDistance*LSigma
    mainAtomsData%posAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )%y = (j-1)*startDistance*LSigma
    mainAtomsData%posAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )%z = (k-1)*startDistance*LSigma

    pdcAtomsData%velAtoms( (i-1)*dimen**2 + (j-1)*dimen + k ) = pdcAtomsData%velAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )*0.0
    pdcAtomsData%posAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )%x = (i-1)*startDistance*LSigma
    pdcAtomsData%posAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )%y = (j-1)*startDistance*LSigma
    pdcAtomsData%posAtoms( (i-1)*dimen**2 + (j-1)*dimen + k )%z = (k-1)*startDistance*LSigma
    end do
  end do
end do

!call startRevEiler(5000,0.01,mainAtomsData)
call initAtomsData(dimen,mainAtomsData,startDistance,10000,0.01)

write(*,*) (calcPEnergy(mainAtomsData) + calcKEnergy(mainAtomsData))

open(12,file = '../test1.xyz')
write(12,*) mainAtomsData%aCount
write(12,*) 'somethingnew'
do i  = 1 , mainAtomsData%aCount
  write(12,*) 'Cr', mainAtomsData%posAtoms(i)%x, mainAtomsData%posAtoms(i)%y, mainAtomsData%posAtoms(i)%z
end do
close(12)


!### Global tick

open(12,file = '../newtest/test2.xyz')
write(12,*) mainAtomsData%aCount
write(12,*) 'somethingnew'
do i  = 1 , mainAtomsData%aCount
  write(12,*) 'Cr', mainAtomsData%posAtoms(i)%x, mainAtomsData%posAtoms(i)%y, mainAtomsData%posAtoms(i)%z
end do
close(12)


write(*,*) "that's fine it's working"

end program main

!###

subroutine writeData(lMainAtomsData,i)
  use Data
  implicit none
  type(atomsData) :: lMainAtomsData
  integer :: i
  character(len = 10) :: name

  open (15,file ='../newtest/tmp')
  write(15,'(I0.3,A4)') i,'.xyz'
  close(15)
  open (15,file ='../newtest/tmp')
  read(15,*) name
  close(15)
  write(*,*) name

  open(14,file = '../newtest/'//name)

  write(14,*) lMainAtomsData%aCount
  write(14,*) 'somethingnew'
  do i  = 1 , lMainAtomsData%aCount
    write(14,*) 'Cr', lMainAtomsData%posAtoms(i)%x, lMainAtomsData%posAtoms(i)%y, lMainAtomsData%posAtoms(i)%z
  end do
  close(14)
end subroutine writeData
