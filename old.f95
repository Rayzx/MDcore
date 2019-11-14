program main
implicit none
!#########FLAGS###########!
integer, parameter :: isWriteCordsToFile = 0
integer, parameter :: isConsoleMode = 1
integer, parameter :: isMinMode = 0
integer, parameter :: isAuMode = 0
integer, parameter :: isInterMode = 0
!##########################!
real, parameter :: cul = 3.60
integer, parameter :: dimen = 11

!####t###
integer :: flag1
real :: minPotent , this , previous , resol , delta , di, en, ren
!####t###


real :: calcDist
real :: calcMors
real :: calcLen
real :: calcFPotent
real :: calcCpotent
real :: calcNpotent

real :: sumPotent , fl , EOA, EOV

integer :: i,j,k,c, N1, Cn, RelaxCount

real :: X1,Y1,Z1

real, dimension(:) , allocatable :: Xcu,Ycu,Zcu,Enu
integer, dimension(:) , allocatable :: Relax, RelaxEnd

allocate(Xcu(dimen ** 3 + 3 * dimen * (dimen - 1) ** 2))
allocate(Ycu(dimen ** 3 + 3 * dimen * (dimen - 1) ** 2))
allocate(Zcu(dimen ** 3 + 3 * dimen * (dimen - 1) ** 2))
allocate(Enu(dimen ** 3 + 3 * dimen * (dimen - 1) ** 2))
allocate(Relax(dimen ** 3 + 3 * dimen * (dimen - 1) ** 2))

call GenerateCubic(cul , dimen , 1 , Xcu , Ycu , Zcu)


        !###calcPotent###!
sumPotent = 0
if (isConsoleMode.eq.1) then

        N1 = dimen**3 / 2 +1
        X1 = Xcu(N1)
        Y1 = Ycu(N1)
        Z1 = Zcu(N1)

        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                sumPotent = sumPotent + calcMors ( calcDist( X1,Y1,Z1,Xcu(i),Ycu(i),Zcu(i) ) )
        end do

        write(*,*) 'Central Atom potential = ', sumPotent , N1

        sumPotent = 0
        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                sumPotent = sumPotent + calcMors ( calcDist( Xcu(1),Ycu(1),Zcu(1),Xcu(i),Ycu(i),Zcu(i) ) )
        end do

        write(*,*) 'Corner Atom potential = ', sumPotent, 1

        sumPotent = 0
        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                sumPotent = sumPotent + calcMors &
                ( calcDist( Xcu(dimen ** 3),Ycu(dimen ** 3),Zcu(dimen ** 3),Xcu(i),Ycu(i),Zcu(i) ) )
        end do

        write(*,*) 'Corner Atom potential test= ', sumPotent, dimen ** 3

        sumPotent = 0
        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                sumPotent = sumPotent + calcMors ( calcDist( Xcu(dimen),Ycu(dimen),Zcu(dimen),Xcu(i),Ycu(i),Zcu(i) ) )
        end do

        write(*,*) 'Corner Atom potential test2= ', sumPotent , dimen


        sumPotent = 0
        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        sumPotent = sumPotent + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                end do
        end do

        write(*,*) 'Full Atom potential = ', sumPotent
        fl = sumPotent

        sumPotent = 0
        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        if ( (i.ne.N1).and.(j.ne.N1) ) then
                                sumPotent = sumPotent + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                        end if
                end do
        end do

        write(*,*) 'Full Atom potential without central test= ', sumPotent
        write(*,*) 'Delta =', sumPotent-fl
        write(*,*) 'sublim =', fl / (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2)


        Xcu(N1) = 1000
        Ycu(N1) = 1000
        Zcu(N1) = 1000

        sumPotent = 0
        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        sumPotent = sumPotent + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                end do
        end do

        write(*,*) 'Full Atom potential  Magic = ', sumPotent

        write(*,*) 'Energu of vacansy = ', sumPotent - fl  + fl / (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2)

endif


        !###min mode###!
if (isMinMode.eq.1) then


       open (12,file = 'min3.txt')
       do c = 1, 400
                call GenerateCubic(3.0+7.0/400*c , dimen , 1 , Xcu , Ycu , Zcu)
                sumPotent = 0
                do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                                sumPotent = sumPotent + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                        end do
                end do
                sumPotent = sumPotent

                write(*,*) 'Full Atom potential = ', sumPotent
                write(12,*) 3.0+7.0/400*c , sumPotent
        end do
        close(12)
end if


         !###Inter mode###!
if (isInterMode.eq.1) then
        N1 = dimen**3 / 2 +1
                X1 = Xcu(N1)
                Y1 = Ycu(N1)
                Z1 = Zcu(N1)

        sumPotent = 0
                do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                                sumPotent = sumPotent + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                        end do
                end do

                write(*,*) 'Full Atom potential = ', sumPotent

        fl = 0
                do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                                if ( (i.ne.N1).and.(j.ne.N1) ) then
                                        fl = fl + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                                end if
                        end do
                end do

                write(*,*) 'Full Atom potential without central = ', fl

                write(*,*) 'sublim =', fl / (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2)
                EOV = fl - sumPotent  + sumPotent / (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2)
                write(*,*) 'Energu of vacansy = ', EOV

                Xcu(N1) = Xcu(N1) + cul/2
                Ycu(N1) = Ycu(N1) + cul/2
                Zcu(N1) = Zcu(N1) + cul/2

        fl = 0
                do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                                        fl = fl + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                        end do
                end do

                EOA = fl - sumPotent
                write(*,*) 'Energy of atom = ', EOA

                write(*,*) 'equal = ', EOA / EOV
endif


        !###out crystal###!
if (isWriteCordsToFile.eq.1) then
        write(*,*) 'writing to file'
        open(11 ,file = 'out.xyz')

        write(11,*) dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
        write(11,*) 'O is also Cu but was exchanged for more contrast '

        do i = 1, dimen ** 3
                write(11,*) 'Cu' , Xcu( i ) , Ycu( i ) , Zcu( i )
        end do
        do i = dimen ** 3 + 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                write(11,*) 'O' , Xcu( i ) , Ycu( i ) , Zcu( i )
        end do

        close(11)
endif


        !###MinAuMode###!
if (isAuMode.eq.1) then
        call GenerateCubic(2.0, dimen , 1 , Xcu , Ycu , Zcu)
        previous =  calcFPotent(cul, dimen, Xcu, Ycu, Zcu)
        call GenerateCubic(2.1, dimen , 1 , Xcu , Ycu , Zcu)
        this =  calcFPotent(cul, dimen, Xcu, Ycu, Zcu)
        di = 2.1
        flag1 = 1
        resol = 0.4
        en = previous
        do i  = 1, 5
                do while ((this-previous).le.0)
                        previous = this
                        di = di + resol
                        call GenerateCubic(di, dimen , 1 , Xcu , Ycu , Zcu)
                        this = calcCPotent(dimen, Xcu, Ycu, Zcu)
                end do

                if (previous.le.en) then
                        en = previous
                        ren = di
                end if

                write(*,*) en
                resol = resol/3.0
                previous = this
                di = di - resol
                call GenerateCubic(di, dimen , 1 , Xcu , Ycu , Zcu)
                this = calcCPotent(dimen, Xcu, Ycu, Zcu)

                do while ((this-previous).le.0)
                        previous = this
                        di = di - resol
                        call GenerateCubic(di, dimen , 1 , Xcu , Ycu , Zcu)
                        this = calcCPotent(dimen, Xcu, Ycu, Zcu)
                end do

                if (previous.le.en) then
                        en = previous
                        ren = di
                end if
                resol = resol/3
                previous = this
                di = di + resol
                call GenerateCubic(di, dimen , 1 , Xcu , Ycu , Zcu)
                this = calcCPotent( dimen, Xcu, Ycu, Zcu)
        end do

        write(*,*) 'full potent = ',  en, ren
endif


        !###TestMode###!
N1 = dimen**3 / 2 +1
        X1 = Xcu(N1)
        Y1 = Ycu(N1)
        Z1 = Zcu(N1)
c = 1
do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
        if ( ( Xcu(i).le.(X1+2) ).and.( Xcu(i).ge.(X1-2) ).and.( Ycu(i).le.(Y1+2) ).and.(&
        Ycu(i).ge.(Y1-2) ).and.( Zcu(i).le.(Z1+2) ).and.( Zcu(i).ge.(Z1-2) ).and.(i.ne.N1) ) then
        !write(*,*) i, c
        Relax(c) = i
        c = c+1
        end if
end do
do i = 1, c-1
        write(*,*) Relax(i)
enddo
allocate(RelaxEnd(c))    !###allocate(RelaxEnd(c-1))
do i = 1, c-1
        RelaxEnd(i) = Relax(i)
enddo
do i = 1, c-1
        write(*,*) RelaxEnd(i)
enddo
RelaxCount = c
RelaxEnd(c) = N1
Relax(c) = N1

!### hard test
!Xcu(N1) = 10000
!Ycu(N1) = 10000
!Zcu(N1) = 10000

Xcu(N1) = X1*1.05
Ycu(N1) = Y1*1.05
Zcu(N1) = Z1*1.05


!RelaxCount = c
!RelaxEnd(c) = N1
!####
!Xcu(N1) = Xcu(N1) + 0.05*cul
!Ycu(N1) = Ycu(N1) + 0.05*cul
!Zcu(N1) = Zcu(N1) + 0.05*cul

!call GenerateRelax (cul,dimen,Xcu,Ycu,Zcu,RelaxCount,RelaxEnd)

!Xcu(N1) = X1
!Ycu(N1) = Y1
!Zcu(N1) = Z1


write(*,*) c
write(*,*) 'writing to file'
        open(11 ,file = 'out2.xyz')
        write(11,*) dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
        write(11,*) 'O is also Cu but was exchanged for more contrast '
        c = 1
        do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                if ((i).eq.(Relax(c))) then
                        write(11,*) 'O' , Xcu( i ) , Ycu( i ) , Zcu( i )
                        c = c + 1
                else
                        if (i.le.(dimen**3)) then
                                write(11,*) 'Cu' , Xcu( i ) , Ycu( i ) , Zcu( i )
                        else
                                write(11,*) 'W' , Xcu( i ) , Ycu( i ) , Zcu( i )
                        end if
                endif

        end do
        close(11)


end

real function calcDist(Xc1,Yc1,Zc1,Xc2,Yc2,Zc2)
implicit none
        real :: Xc1,Yc1,Zc1,Xc2,Yc2,Zc2

calcDist = sqrt( (Xc1-Xc2)**2 + (Yc1-Yc2)**2 + (Zc1-Zc2)**2  )
end function


real function calcMors(dist)
implicit none
        real :: dist
        real, parameter :: Dkoef = 0.3429
        real, parameter :: Akoef = 1.3588
        real, parameter :: r0 = 2.866
if (dist.le.(0.1)) then
        calcMors = 0
else
        calcMors = Dkoef * ( exp( -2 * Akoef * ( dist - r0 ) )  - 2 * exp(  -Akoef * ( dist - r0 ) ) )
endif

end function


real function calcLen(dist)
implicit none
        real :: dist
        real, parameter :: Skoef = 2.27
        real, parameter :: Ekoef = 2.332
if (dist.le.(0.1)) then
        calcLen = 0
else
        calcLen = Ekoef*( (Skoef/dist)**12 - (Skoef/dist)**6)
endif

end function


real function calcFPotent(cul, dimen, Xcu, Ycu, Zcu)
implicit none
        integer :: dimen
        real :: calcMors
        real :: calcDist
        real :: cul
        real , dimension (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2) :: Xcu , Ycu , Zcu
        integer :: i, j


        calcFPotent = 0
                do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        do j = i, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                                calcFPotent = calcFPotent + calcMors ( calcDist( Xcu(i),Ycu(i),Zcu(i),Xcu(j),Ycu(j),Zcu(j) ) )
                        end do
                end do


end function


real function calcCPotent(dimen, Xcu, Ycu , Zcu)
implicit none
        integer :: dimen
        real :: calcMors
        real :: calcDist
        real , dimension (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2) :: Xcu , Ycu , Zcu
        integer :: i
        real :: N1, X1, Y1, Z1

        N1 = dimen**3 / 2 +1
                X1 = Xcu(N1)
                Y1 = Ycu(N1)
                Z1 = Zcu(N1)
        calcCPotent = 0
                do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        calcCPotent = calcCPotent + calcMors ( calcDist( X1,Y1,Z1,Xcu(i),Ycu(i),Zcu(i) ) )
                end do
end function


real function calcNPotent(dimen, Xcu, Ycu , Zcu, N)
implicit none
        integer :: dimen
        real :: calcMors
        real :: calcDist
        real , dimension (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2) :: Xcu , Ycu , Zcu
        integer :: i , N
        real ::  X1, Y1, Z1
                X1 = Xcu(N)
                Y1 = Ycu(N)
                Z1 = Zcu(N)
        calcNPotent = 0
                do i = 1, dimen ** 3 + 3 * dimen * (dimen - 1) ** 2
                        calcNPotent = calcNPotent + calcMors ( calcDist( X1,Y1,Z1,Xcu(i),Ycu(i),Zcu(i) ) )
                end do
end function

subroutine GenerateCubic ( cul, dimen , doesGCR , Xcu , Ycu , Zcu)
implicit none
        integer :: doesGCR
        real  :: cul
        real , dimension (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2) :: Xcu , Ycu , Zcu
        integer :: dimen
        integer :: i , j , k, c
        if (doesGCR.eq.1) then
                c=0
                do i = 1,  dimen
                        do j = 1, dimen
                                do k=1, dimen
                                        c = c + 1
                                        Xcu(c) =  cul*(i-1)
                                        Ycu(c) =  cul*(j-1)
                                        Zcu(c) =  cul*(k-1)
                                end do
                        end do
                end do


                do i = 1,  dimen
                        do j = 1, dimen -1
                                do k=1, dimen -1
                                        c = c + 1
                                        Xcu(c) =  cul * (i - 1)
                                        Ycu(c) =  cul * (j - 1) + 0.5 * cul
                                        Zcu(c) =  cul * (k - 1) + 0.5 * cul
                                end do
                        end do
                end do

                do i = 1,  dimen -1
                        do j = 1, dimen
                                do k=1, dimen -1
                                        c = c + 1
                                        Xcu(c) =  cul * (i - 1) + 0.5 * cul
                                        Ycu(c) =  cul * (j - 1)
                                        Zcu(c) =  cul * (k - 1) + 0.5 * cul
                                end do
                        end do
                end do

                do i = 1,  dimen - 1
                        do j = 1, dimen - 1
                                do k= 1, dimen
                                        c = c + 1
                                        Xcu(c) =  cul * (i - 1) + 0.5 * cul
                                        Ycu(c) =  cul * (j - 1) + 0.5 * cul
                                        Zcu(c) =  cul * (k - 1)
                                end do
                        end do
                end do
        end if


end subroutine GenerateCubic

subroutine GenerateRelax ( cul , dimen , Xcu , Ycu , Zcu , RelaxCount , Relax)
implicit none
real :: calcNPotent
        real  :: cul
        integer :: dimen, RelaxCount
        real , dimension (dimen ** 3 + 3 * dimen * (dimen - 1) ** 2) :: Xcu , Ycu , Zcu
        integer , dimension (RelaxCount) :: Relax
        real , dimension (RelaxCount) :: Xb, Yb, Zb
        integer :: i , j , k, c, N, f
        real :: dh, E0, dEx, dEy, dEz, dx,dy,dz, D
        dh = 0.5

        do i  = 1, RelaxCount
                Xb(i)  = Xcu(Relax(i))
                Yb(i)  = Ycu(Relax(i))
                Zb(i)  = Zcu(Relax(i))
        end do

        do c = 1, 50
        write(*,*) ' one '
        do i = 1 , RelaxCount
                f = 0
                dh = 0.1
                N = Relax(i)
                E0 = calcNPotent(dimen,Xcu,Ycu,Zcu,N)

                write(*,*) N


                do while(f.eq.0)
                        Xcu(N)=Xcu(N)+dh
                        dEx = (calcNPotent(dimen,Xcu,Ycu,Zcu,N) - E0)/dh
                        Xcu(N)=Xcu(N)-dh

                        write(*,*) ' ', dEx

                        Ycu(N)=Ycu(N)+dh
                        dEy = (calcNPotent(dimen,Xcu,Ycu,Zcu,N) - E0)/dh
                        Ycu(N)=Ycu(N)-dh

                        Zcu(N)=Zcu(N)+dh
                        dEz = (calcNPotent(dimen,Xcu,Ycu,Zcu,N) - E0)/dh
                        Zcu(N)=Zcu(N)-dh

                        D = sqrt(dEx*dEx + dEy*dEy + dEz*dEz) + 0.00001

                        dx = -dh*dEx/D
                        dy = -dh*dEy/D
                        dz = -dh*dEz/D
                        write(*,*) '     ', dx, dy, dz

                        Xcu(N)=Xcu(N)+dx
                        Ycu(N)=Ycu(N)+dy
                        Zcu(N)=Zcu(N)+dz

                        if ( (calcNPotent(dimen, Xcu, Ycu, Zcu, N)).le.E0) then
                                f = 1
                        else
                                Xcu(N)=Xcu(N)-dx
                                Ycu(N)=Ycu(N)-dy
                                Zcu(N)=Zcu(N)-dz
                                dh = dh / 2
                                if  (dh.le.0.001) then
                                        f =1
                                end if
                        end if
                end do
        end do
        end do


end subroutine GenerateRelax

