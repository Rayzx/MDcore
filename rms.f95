


program rms
  implicit none
  real, dimension(100) :: res
  real :: calcRMS
  integer :: i, n
  character(len = 10) :: name

  !char(len = 8) :: fileName
  !char(len = 5), parameter :: xyz='.xyz'

  write(*,*) calcRMS(10)



  open (15,file ='./newtest/tmp')
  write(15,'(I0.3,A4)') i,'.xyz'
  close(15)
  open (15,file ='./newtest/tmp')
  read(15,*) name
  close(15)
  write(*,*) name


end program rms



real function calcRMS(n)
  real, dimension(:), allocatable :: a1, a2
  real :: c,s
  integer :: i,n

  open(11,file="outE", status='old', action='read')

  allocate(a1(n))
  allocate(a2(n))

  do i = 1, n
    read(11,*) a1(i), a2(i)
  end do
    close(11)

  do i = 1, n
    s = s + a2(i)
  end do
  s = s/n
  c = 0
  do i = 1, n
  c = (s-a2(i))*(s-a2(i))
  end do
  c = sqrt(c/n)
  calcRMS = log10(c)
end function calcRMS
