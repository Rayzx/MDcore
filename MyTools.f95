module WConst
  implicit none

  real, parameter :: LSigma = 1
  real, parameter :: LTerminate = 3
  real, parameter :: LEpsilon = 1
  real, parameter :: m = 1
end module WConst


module MyVector
  implicit none

  type vector3
    real  :: x
    real  :: y
    real  :: z
  end type vector3

  interface operator (+)
    procedure vAdd
  end interface operator (+)
  interface operator (-)
    procedure vSub
  end interface operator (-)
  interface operator (*)
    procedure vDot
  end interface operator (*)

  contains
  !###
  real  function vectorAbs(cords)
    implicit none
    type(vector3) :: cords

    vectorAbs = sqrt(cords%x * cords%x + cords%y * cords%y + cords%z * cords%z)
  end function vectorAbs
  !###
  function vAdd(v1,v2)
    implicit none
    type(vector3), intent(in) :: v1,v2
    type(vector3) :: vAdd

    vAdd%x = v1%x + v2%x
    vAdd%y = v1%y + v2%y
    vAdd%z = v1%z + v2%z
  end function vAdd
  !###
  function vSub(v1,v2)
    implicit none
    type(vector3), intent(in) :: v1,v2
    type(vector3) :: vSub

    vSub%x = v1%x - v2%x
    vSub%y = v1%y - v2%y
    vSub%z = v1%z - v2%z
  end function vSub
  !###
  function vDot(v1,alpha)
    implicit none
    type(vector3), intent(in) :: v1
    real , intent(in) :: alpha
    type(vector3) :: vDot

    vDot%x = v1%x * alpha
    vDot%y = v1%y * alpha
    vDot%z = v1%z * alpha
  end function vDot
  !###
end module MyVector


module Data
  use MyVector
  implicit none

  type atomsData
    integer :: aCount
    type(vector3), dimension(:) , allocatable :: posAtoms
    type(vector3), dimension(:) , allocatable :: velAtoms
    type(vector3), dimension(:) , allocatable :: forceAtoms
  end type atomsData
end module Data


module LenModule
  use WConst
  use MyVector
  use Data
  implicit none

  contains
  !###
  real  function calcLenD(lDistance)
    implicit none
    real  :: lDistance
    if (lDistance.le.LTerminate) then
      if (lDistance.le.0.01) then
        calcLenD = 0
        write(*,*) 'too close'
      else
        calcLenD = 4*LEpsilon*( (LSigma/lDistance)**12 - (LSigma/lDistance)**6 )
      end if
    else
      calcLenD = 0
    end if
  end function
  !###
  real  function calcLenV(vAtom1,vAtom2)
    implicit none
    type(vector3) :: vAtom1,vAtom2
    real  :: lDistance

    lDistance = vectorAbs(vAtom1-vAtom2)
    calcLenV = calcLenD(lDistance)
  end function
  !###
  real function calcKEnergy(lMainAtomsData)
    implicit none
    type(atomsData) :: lMainAtomsData
    integer :: i

    calcKEnergy = 0.0
    do i = 1, lMainAtomsData%aCount
      calcKEnergy = calcKEnergy + m*vectorAbs( lMainAtomsData%velAtoms(i) )**2/2
    end do
  end function
  !###
  real function calcPEnergy(lMainAtomsData)
    implicit none
    type(atomsData) :: lMainAtomsData
    integer :: i,j

    calcPEnergy = 0
    do i = 1, lMainAtomsData%aCount
      do j = 1, lMainAtomsData%aCount
        if (i.ne.j) then
        calcPEnergy = calcPEnergy + calcLenV(lMainAtomsData%posAtoms(i),lMainAtomsData%posAtoms(j))
        endif
      end do
    end do
    calcPEnergy = calcPEnergy / 2
  end function
  !###
  real  function calcLenForceX(vAtomTo,vAtomFrom)
    implicit none
    type(vector3) :: vAtomTo, vAtomFrom
    real  :: lDistance, lX

    lDistance = vectorAbs(vAtomTo-vAtomFrom)
    if (lDistance.le.LTerminate) then
      lX = - vAtomTo%x + vAtomFrom%x
      calcLenForceX = -24*LEpsilon/LSigma*( 2*(LSigma/lDistance)**13 - (LSigma/lDistance)**7 )*lX/lDistance
    else
      calcLenForceX = 0
    endif
  end function
  !###
  real  function calcLenForceY(vAtomTo,vAtomFrom)
    implicit none
    type(vector3) :: vAtomTo, vAtomFrom
    real  :: lDistance, lY

    lDistance = vectorAbs(vAtomTo-vAtomFrom)
    if (lDistance.le.LTerminate) then
      lY = - vAtomTo%y + vAtomFrom%y
      calcLenForceY = -24*LEpsilon/LSigma*( 2*(LSigma/lDistance)**13 - (LSigma/lDistance)**7 )*lY/lDistance
    else
      calcLenForceY = 0
    end if
  end function
  !###
  real  function calcLenForceZ(vAtomTo,vAtomFrom)
    implicit none
    type(vector3) :: vAtomTo, vAtomFrom
    real  :: lDistance, lZ

    lDistance = vectorAbs(vAtomTo-vAtomFrom)
    if (lDistance.le.LTerminate) then
      lZ = - vAtomTo%z + vAtomFrom%z
      calcLenForceZ = -24*LEpsilon/LSigma*( 2*(LSigma/lDistance)**13 - (LSigma/lDistance)**7 )*lZ/lDistance
    else
      calcLenForceZ = 0
    endif
  end function
  !###
  real function calcPbcEnergy(lMainAtomsData)
    implicit none
    type(atomsData) :: lMainAtomsData
    integer :: i,j
    calcPbcEnergy = 0
    do i = 1, lMainAtomsData%aCount
      do j = 1, lMainAtomsData%aCount*27
        if (i.ne.j) then
        calcPbcEnergy = calcPbcEnergy + calcLenV(lMainAtomsData%posAtoms(i),lMainAtomsData%posAtoms(j))
        endif
      end do
    end do
    calcPbcEnergy = calcPbcEnergy / 2
  end function
  !###
end module LenModule


module AlgRevEiler
  use MyVector
  use data
  use WConst
  use LenModule
  implicit none

  contains
  !###
  subroutine startRevEiler(n,dt,mainAtomsData)
    implicit none
    integer :: i,n
    real :: dt
    type(atomsData) :: mainAtomsData

    write(*,*) n
    open(11, file = '../outE')
    do i=1, n
      call updateForces(mainAtomsData) !re
      call updateVel(mainAtomsData,dt)
      write(11,*) i*dt, (calcPEnergy(mainAtomsData) + calcKEnergy(mainAtomsData))
      call updatePos(mainAtomsData,dt)
      !write(*,*) 'hello11s'
      call writeData1(mainAtomsData,i)
    end do
    write(*,*) i
    close(11)
  end subroutine startRevEiler
  !###
  subroutine updateForces(lMainAtomsData)
    implicit none
    type(atomsData) lMainAtomsData
    integer :: i,j

    do i = 1, lMainAtomsData%aCount
      lMainAtomsData%forceAtoms(i)%x = 0.0
      lMainAtomsData%forceAtoms(i)%y = 0.0
      lMainAtomsData%forceAtoms(i)%z = 0.0
      do j = 1, lMainAtomsData%aCount
        if (i.ne.j) then
          lMainAtomsData%forceAtoms(i)%x = lMainAtomsData%forceAtoms(i)%x &
            + calcLenForceX( lMainAtomsData%posAtoms(i), lMainAtomsData%posAtoms(j) )
          lMainAtomsData%forceAtoms(i)%y = lMainAtomsData%forceAtoms(i)%y &
            + calcLenForceY( lMainAtomsData%posAtoms(i), lMainAtomsData%posAtoms(j) )
          lMainAtomsData%forceAtoms(i)%z = lMainAtomsData%forceAtoms(i)%z &
            + calcLenForceZ( lMainAtomsData%posAtoms(i), lMainAtomsData%posAtoms(j) )
        endif
      end do
    end do
  end subroutine updateForces
  !###
  subroutine updateVel(lMainAtomsData, dt)
    implicit none
    type(atomsData) :: lMainAtomsData
    real  :: dt
    integer :: i,j

    do i = 1, lMainAtomsData%aCount
      lMainAtomsData%velAtoms(i)%x = lMainAtomsData%velAtoms(i)%x + lMainAtomsData%forceAtoms(i)%x/m*dt
      lMainAtomsData%velAtoms(i)%y = lMainAtomsData%velAtoms(i)%y + lMainAtomsData%forceAtoms(i)%y/m*dt
      lMainAtomsData%velAtoms(i)%z = lMainAtomsData%velAtoms(i)%z + lMainAtomsData%forceAtoms(i)%z/m*dt
    end do
  end subroutine updateVel
  !###
  subroutine updatePos(lMainAtomsData, dt)
    implicit none
    type(atomsData) :: lMainAtomsData
    real :: dt
    integer :: i,j

    do i = 1, lMainAtomsData%aCount
      lMainAtomsData%posAtoms(i)%x = lMainAtomsData%posAtoms(i)%x + lMainAtomsData%velAtoms(i)%x*dt
      lMainAtomsData%posAtoms(i)%y = lMainAtomsData%posAtoms(i)%y + lMainAtomsData%velAtoms(i)%y*dt
      lMainAtomsData%posAtoms(i)%z = lMainAtomsData%posAtoms(i)%z + lMainAtomsData%velAtoms(i)%z*dt
    end do
  end subroutine updatePos
  !###

  subroutine writeData1(lMainAtomsData,i)
    implicit none
    type(atomsData) :: lMainAtomsData
    integer :: i,j
    character(len = 10) :: name

    open (15,file ='../newtest/tmp')
    write(15,'(I0.3,A4)') i,'.xyz'
    close(15)
    open (15,file ='../newtest/tmp')
    read(15,*) name
    close(15)
    write(*,*) name

    open(14,file = '../newtest/'//name)

    write(14,*) lmainAtomsData%aCount
    write(14,*) 'somethingnew'
    do j  = 1 , lmainAtomsData%aCount
      write(14,*) 'Cr', lmainAtomsData%posAtoms(j)%x, lmainAtomsData%posAtoms(j)%y, lmainAtomsData%posAtoms(j)%z
    end do
    close(14)
  end subroutine writeData1
end module AlgRevEiler


module AlgPbcRevEiler
  use AlgRevEiler
  use MyVector
  use data
  use WConst
  use LenModule
  implicit none
  type(atomsData) :: pbcAtomsData
  real :: lstartDistance

  contains
  !###
  subroutine initAtomsData(dimen, lMainAtomsData, inStartDistance,n,dt)
    implicit none
    integer :: dimen,i,n
    type(atomsData) :: lMainAtomsData
    real :: inStartDistance,dt

    lstartDistance = inStartDistance
    pbcAtomsData%aCount = dimen * dimen * dimen
    allocate(pbcAtomsData%posAtoms(27*pbcAtomsData%aCount))
    allocate(pbcAtomsData%velAtoms(pbcAtomsData%aCount))
    allocate(pbcAtomsData%forceAtoms(pbcAtomsData%aCount))
    do i = 1, pbcAtomsData%aCount
      pbcAtomsData%posAtoms(i) = lMainAtomsData%posAtoms(i)
      pbcAtomsData%velAtoms(i) = lMainAtomsData%velAtoms(i)
      pbcAtomsData%forceAtoms(i) = lMainAtomsData%forceAtoms(i)
    end do
    call makeImages(dimen)
    call printData
    call startPbcRevEiler(n,dt,dimen)
  end subroutine initAtomsData
  !###
  subroutine startPbcRevEiler(n,dt,dimen)
    implicit none
    integer :: i,n, dimen
    real :: dt

    write(*,*) n
    open(17, file = '../outE3')
    do i=1, n
      call pbcUpdateForce
      call updateVel(pbcAtomsData,dt)
      write(17,*) i*dt, (calcPbcEnergy(pbcAtomsData) + calcKEnergy(pbcAtomsData))
      call updatePos(pbcAtomsData,dt)
      call resolvePos(dimen)
      call makeImages(dimen)
      call writeData2(pbcAtomsData,i,dimen)
    end do
    write(*,*) i
    close(17)
  end subroutine startPbcRevEiler
  !###
  subroutine makeImages(dimen)
    integer :: i,k,j,l,dimen
    do l = 1, 3
    do j = 1, 3
    do k = 1, 3
      if ((k+3*j-3+9*(l-1)).lt.14) then
        do i = 1, pbcAtomsData%aCount
          pbcAtomsData%posAtoms(i+pbcAtomsData%aCount*(k+3*j-3+9*(l-1)))%x &
              = pbcAtomsData%posAtoms(i)%x + dimen*LSigma*lstartDistance*(k-2)
          pbcAtomsData%posAtoms(i+pbcAtomsData%aCount*(k+3*j-3+9*(l-1)))%y &
              = pbcAtomsData%posAtoms(i)%y + dimen*LSigma*lstartDistance*(j-2)
          pbcAtomsData%posAtoms(i+pbcAtomsData%aCount*(k+3*j-3+9*(l-1)))%z &
              = pbcAtomsData%posAtoms(i)%z + dimen*LSigma*lstartDistance*(l-2)
        end do
      end if
      if ((k+3*j-3+9*(l-1)).gt.14) then
        do i = 1, pbcAtomsData%aCount
          pbcAtomsData%posAtoms(i+pbcAtomsData%aCount*(k+3*j-3-1+9*(l-1)))%x &
              = pbcAtomsData%posAtoms(i)%x + dimen*LSigma*lstartDistance* (k-2)
          pbcAtomsData%posAtoms(i+pbcAtomsData%aCount*(k+3*j-3-1+9*(l-1)))%y &
              = pbcAtomsData%posAtoms(i)%y + dimen*LSigma*lstartDistance* (j-2)
          pbcAtomsData%posAtoms(i+pbcAtomsData%aCount*(k+3*j-3-1+9*(l-1)))%z &
              = pbcAtomsData%posAtoms(i)%z + dimen*LSigma*lstartDistance* (l-2)
        end do
      end if
    end do
    end do
    end do
  end subroutine makeImages
  !###
  subroutine printData
    implicit none
    integer :: i
    open(15,file = '../test3.xyz')
    write(15,*) pbcAtomsData%aCount*9
    write(15,*) 'somethingnew'
    do i  = 1 , pbcAtomsData%aCount*9
      write(15,*) 'Cr', pbcAtomsData%posAtoms(i)%x, pbcAtomsData%posAtoms(i)%y, pbcAtomsData%posAtoms(i)%z
    end do
    close(15)
  end subroutine printData
  !###
  subroutine resolvePos(dimen)
    implicit none
    integer :: i,dimen
    do i = 1, pbcAtomsData%aCount
      if ((pbcAtomsData%posAtoms(i)%x).lt.(-1.0*lstartDistance*LSigma/2)) then
        pbcAtomsData%posAtoms(i)%x = pbcAtomsData%posAtoms(i)%x + lstartDistance*LSigma*dimen
      end if
      if ((pbcAtomsData%posAtoms(i)%x).gt.(lstartDistance*LSigma*dimen-lstartDistance*LSigma/2)) then
        pbcAtomsData%posAtoms(i)%x = pbcAtomsData%posAtoms(i)%x - lstartDistance*LSigma*dimen
      end if
      if ((pbcAtomsData%posAtoms(i)%y).lt.(-1.0*lstartDistance*LSigma/2)) then
        pbcAtomsData%posAtoms(i)%y = pbcAtomsData%posAtoms(i)%y + lstartDistance*LSigma*dimen
      end if
      if ((pbcAtomsData%posAtoms(i)%y).gt.(lstartDistance*LSigma*dimen-lstartDistance*LSigma/2)) then
        pbcAtomsData%posAtoms(i)%y = pbcAtomsData%posAtoms(i)%y - lstartDistance*LSigma*dimen
      end if
      if ((pbcAtomsData%posAtoms(i)%z).lt.(-1.0*lstartDistance*LSigma/2)) then
        pbcAtomsData%posAtoms(i)%z = pbcAtomsData%posAtoms(i)%z + lstartDistance*LSigma*dimen
      end if
      if ((pbcAtomsData%posAtoms(i)%z).gt.(lstartDistance*LSigma*dimen-lstartDistance*LSigma/2)) then
        pbcAtomsData%posAtoms(i)%z = pbcAtomsData%posAtoms(i)%z - lstartDistance*LSigma*dimen
      end if
    end do
  end subroutine resolvePos
  !###
  subroutine pbcUpdateForce
    implicit none
    integer :: i,j

    do i = 1, pbcAtomsData%aCount
      pbcAtomsData%forceAtoms(i)%x = 0.0
      pbcAtomsData%forceAtoms(i)%y = 0.0
      pbcAtomsData%forceAtoms(i)%z = 0.0
      do j = 1, pbcAtomsData%aCount*27
        if (i.ne.j) then
          pbcAtomsData%forceAtoms(i)%x = pbcAtomsData%forceAtoms(i)%x &
            + calcLenForceX( pbcAtomsData%posAtoms(i), pbcAtomsData%posAtoms(j) )
          pbcAtomsData%forceAtoms(i)%y = pbcAtomsData%forceAtoms(i)%y &
            + calcLenForceY( pbcAtomsData%posAtoms(i), pbcAtomsData%posAtoms(j) )
          pbcAtomsData%forceAtoms(i)%z = pbcAtomsData%forceAtoms(i)%z &
            + calcLenForceZ( pbcAtomsData%posAtoms(i), pbcAtomsData%posAtoms(j) )
        endif
      end do
    end do
  end subroutine pbcUpdateForce
  !###
  subroutine writeData2(lMainAtomsData,i,dimen)
    implicit none
    type(atomsData) :: lMainAtomsData
    integer :: i,j,dimen
    character(len = 10) :: name

    open (19,file ='../newtest/tmp')
    write(19,'(I0.3,A4)') i,'.xyz'
    close(19)
    open (19,file ='../newtest/tmp')
    read(19,*) name
    close(19)
    write(*,*) name

    open(19,file = '../newtest/'//name)

    write(19,*) lmainAtomsData%aCount*27 + 4
    write(19,*) 'somethingnew'
    do j  = 1 , lmainAtomsData%aCount*27
      write(19,*) 'Cr', lmainAtomsData%posAtoms(j)%x, lmainAtomsData%posAtoms(j)%y, lmainAtomsData%posAtoms(j)%z
    end do
    write(19,*) 'w', -1.0*lstartDistance*LSigma/2, -1.0*lstartDistance*LSigma/2, -0.1
    write(19,*) 'w', lstartDistance*LSigma*dimen-lstartDistance*LSigma/2, -1.0*lstartDistance*LSigma/2, -0.1
    write(19,*) 'w', lstartDistance*LSigma*dimen-lstartDistance*LSigma/2, &
    lstartDistance*LSigma*dimen-lstartDistance*LSigma/2, -0.1
    write(19,*) 'w', -1.0*lstartDistance*LSigma/2, lstartDistance*LSigma*dimen-lstartDistance*LSigma/2, -0.1
    close(19)
  end subroutine writeData2
  !###
end module

module AlgVerlet
  use WConst
  use MyVector
  use Data
  use LenModule
  implicit none

  contains
  subroutine startAlgVerlet(n,dt,mainAtomsData)
    implicit none
    integer :: i,n,j
    real :: dt
    type(atomsData) :: befAtomsData !# -1
    type(atomsData) :: mainAtomsData !# 0
    type(atomsData) :: afterAtomsData  !# +1

    befAtomsData = mainAtomsData
    afterAtomsData = mainAtomsData
    open(13, file = '../outV')
    do j=1, n
      do i=1, MainAtomsData%aCount
        befAtomsData%posAtoms(i) = mainAtomsData%posAtoms(i)
        mainAtomsData%posAtoms(i) = afterAtomsData%posAtoms(i)
        call updateForces(mainAtomsData)
        afterAtomsData%posAtoms(i) = mainAtomsData%posAtoms(i)*(2.0) - befAtomsData%posAtoms(i) + mainAtomsData%forceAtoms(i)*dt*dt
        mainAtomsData%velAtoms(i) = (afterAtomsData%posAtoms(i) - befAtomsData%posAtoms(i))*(1/(2*dt))
      end do
      write(11,*) j*dt, (calcPEnergy(mainAtomsData) + calcKEnergy(mainAtomsData))
    end do
    close(13)
  end subroutine startAlgVerlet
  !###
  subroutine updateForces(lMainAtomsData)
    implicit none
    type(atomsData) lMainAtomsData
    integer :: i,j

    do i = 1, lMainAtomsData%aCount
      lMainAtomsData%forceAtoms(i)%x = 0.0
      lMainAtomsData%forceAtoms(i)%y = 0.0
      lMainAtomsData%forceAtoms(i)%z = 0.0
      do j = 1, lMainAtomsData%aCount
        if (i.ne.j) then
          lMainAtomsData%forceAtoms(i)%x = lMainAtomsData%forceAtoms(i)%x &
            + calcLenForceX( lMainAtomsData%posAtoms(i), lMainAtomsData%posAtoms(j) )
          lMainAtomsData%forceAtoms(i)%y = lMainAtomsData%forceAtoms(i)%y &
            + calcLenForceY( lMainAtomsData%posAtoms(i), lMainAtomsData%posAtoms(j) )
          lMainAtomsData%forceAtoms(i)%z = lMainAtomsData%forceAtoms(i)%z &
            + calcLenForceZ( lMainAtomsData%posAtoms(i), lMainAtomsData%posAtoms(j) )
        endif
      end do
    end do
  end subroutine updateForces
  !###
end module AlgVerlet
