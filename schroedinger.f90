subroutine mappedschroed(npnts,xint,mass,potential,coordderiv,energy,wavefunction)
!=======================================================================
!
!  Name:        schroedinger_1d
!  Description: Solves one-dimensional Schroedinger equation
!               using Fourier Grid Hamiltonian method
!
!  Input parameters:
!
!     npnts - number of grid points (should be power of two)
!     nwf   - number of eigenvectors to calculate/store
!     xint  - integration interval in bohr
!     mass  - mass of the particle in electron mass units
!     potential(1:npnts)  - potential energy on the grid (au)
!     coordderiv(i=1:3,j=1:npnts) - derivative of Cartesian coords wrt scaled coords
!                                   i = derivative order;
!                                   i = grid point.
!
!  Output parameters:
!
!     energy(1:npnts) - energy levels (kcal/mol)
!     wavefunction(i=1:npnts,j=1:nwf) - wavefunctions;
!                                       i - grid point;
!                                       j - quantum number.
!
!  Accuracy: integer*4 and real*8
!
!  Date: 11/26/2003
!
!=======================================================================
implicit none

! input/output variables
integer, intent(in)                          :: npnts
real*8,  intent(in)                          :: xint
real*8,  intent(in)                          :: mass
real*8,  intent(in),  dimension(npnts)       :: potential
real*8,  intent(in),  dimension(3,npnts)     :: coordderiv
real*8,  intent(out), dimension(npnts)       :: energy
real*8,  intent(out), dimension(npnts,npnts) :: wavefunction

! parameters
real*8, parameter :: au2kcal  = 627.5095d0
real*8, parameter :: kcal2au  = 1.d0/au2kcal
real*8, parameter :: bohr2ang = 0.529177249d0
real*8, parameter :: ang2bohr = 1.d0/bohr2ang
real*8, parameter :: pi       = 3.14159265358979d0
real*8, parameter :: hbar     = 1.d0
real*8, parameter :: tol      = 1.d-6

! local variables
integer :: i, j, ierr
real*8  :: xintke
real*8  :: scal
real*8  :: hamiltonian(npnts*(npnts+1)/2)
real*8  :: ekinetic(npnts*(npnts+1)/2)
logical :: power2

!integer :: elemcount
!real*8  :: testKE(npnts,npnts)
!real*8  :: testH(npnts,npnts)
!real*8  :: aux(npnts,npnts)
!real*8  :: junk(npnts)
!real*8  :: workq(1)

real*8, allocatable :: work(:)

hamiltonian = 0.d0
ekinetic = 0.d0
energy = 0.d0
wavefunction = 0.d0

!aux = 0.d0
!junk = 0.d0
!workq = 0.d0
!testKE = 0.0d0
!testH = 0.0d0

! check whether npnts is an integer power of two
if ( npnts.le.0 .or. .not.power2(npnts) ) then
   write(*,*) 'number of grid points (',npnts,') is invalid...'
   stop
endif

! Calculate kinetic energy matrix (Fast Fourier transform)
xintke=npnts*xint/(npnts-1)
!xintke=xint
call mappedgridke(npnts,xintke,mass,coordderiv(1,:),ekinetic)

!! ARS( testing
!   elemcount=1
!   do i=1,npnts
!   do j=i,npnts
!     testKE(i,j)=ekinetic(elemcount)
!     testKE(j,i)=ekinetic(elemcount)
!     elemcount=elemcount+1
!   end do
!   end do
!   open(11,file="mat_T.dat")
!   call primat(testKE,npnts,npnts,npnts,6,9,11,"Kinetic Energy Matrix")
!   close(11)
!! )

do i=1,npnts
   hamiltonian(i+(i-1)*(2*npnts-i)/2) = potential(i) + hbar/(2.0d+00*mass) * (7.0d+00*coordderiv(2,i)*coordderiv(2,i)/(4.0d+00*coordderiv(1,i)*coordderiv(1,i)*coordderiv(1,i)*coordderiv(1,i)) - coordderiv(3,i)/(2.0d+00*coordderiv(1,i)*coordderiv(1,i)*coordderiv(1,i)))
enddo

!! ARS( testing
!   elemcount=1
!   do i=1,npnts
!   do j=i,npnts
!     testH(i,j)=hamiltonian(elemcount)
!     testH(j,i)=hamiltonian(elemcount)
!     elemcount=elemcount+1
!   end do
!   end do
!   open(11,file="mat_V.dat")
!   call primat(testH,npnts,npnts,npnts,6,9,11,"Potential Energy Matrix")
!   close(11)
!! )

hamiltonian = hamiltonian + ekinetic

!! ARS( symmetrize
!   do i=1,npnts-1
!   do j=i+1,npnts
!     junk(1)=(hamiltonian(i,j)+hamiltonian(j,i))/2.0d+00
!     hamiltonian(i,j)=junk(1)
!     hamiltonian(j,i)=junk(1)
!   end do
!   end do
!! )

!! ARS( testing
!   elemcount=1
!   do i=1,npnts
!   do j=i,npnts
!     testH(i,j)=hamiltonian(elemcount)
!     testH(j,i)=hamiltonian(elemcount)
!     elemcount=elemcount+1
!   end do
!   end do
!   open(11,file="mat_H.dat")
!   call primat(testH,npnts,npnts,npnts,6,9,11,"Hamiltonian Matrix")
!   close(11)
!! )

!do i=1,npnts-1
!do j=i+1,npnts
!  if (abs(testH(i,j)-testH(j,i)).gt.tol) then
!   write(*,*) "asymmetry in hamiltonian"
!   write(*,*) "   i,j,h(i,j),h(j,i):",i,j,testH(i,j),testH(j,i)
!  end if
!end do
!end do

allocate(work(3*npnts))
ierr=0
call dspev("v","l",npnts,hamiltonian,energy,wavefunction,npnts,work,ierr)
deallocate(work)
if(ierr.ne.0) then
 write(*,*) "Error in dspev"
 return
end if

!! query work space
!aux(:,:)=testH(:,:)
!ierr=0
!call dsyev("v","l",npnts,aux,npnts,energy,workq,-1,ierr)
!if(ierr.ne.0) then
! write(*,*) "Error in dsyev query"
! return
!end if
!allocate(work(int(workq(1))))
!
!! diagonalize
!aux(:,:)=testH(:,:)
!ierr=0
!call dsyev("v","l",npnts,aux,npnts,energy,work,int(workq(1)),ierr)
!deallocate(work)
!if(ierr.ne.0) then
! write(*,*) "Error in dsyev"
! return
!end if

do i=1,npnts ! grid points
  scal=1.0d+00/dsqrt(coordderiv(1,i))
  call dscal(npnts,scal,wavefunction(i,:),1)
end do

!! query work space
!aux(:,:)=hamiltonian(:,:)
!ierr=0
!call dgeev("n","v",npnts,aux,npnts,energy,junk,junk,npnts,wavefunction,npnts,workq,-1,ierr)
!if(ierr.ne.0) then
! write(*,*) "Error in dgeev query"
! return
!end if
!allocate(work(int(workq(1))))
!
!! diagonalize
!aux(:,:)=hamiltonian(:,:)
!ierr=0
!call dgeev("n","v",npnts,aux,npnts,energy,junk,junk,npnts,wavefunction,npnts,work,int(workq(1)),ierr)
!if(ierr.ne.0) then
! write(*,*) "Error in dgeev"
! return
!end if
!deallocate(work)

energy = energy*au2kcal

return
end subroutine mappedschroed


logical function power2(n)
! Checks whether the given integer (N) is a power of two

implicit none

integer, intent(in) :: n

integer, parameter :: itwo = 2
integer :: nw

power2 = .false.

nw = n
if (nw.lt.2) then
   power2 = .false.
   return
else
   power2 = .true.
endif

do while (nw.gt.1)
   if (mod(nw,itwo).ne.0) then
      power2 = .false.
      return
   else
      nw = nw/itwo
   endif
enddo

return
end function power2


subroutine mappedgridke(n,width,mass,deriv,KEmat)
implicit none

integer, intent(in)  :: n
real*8,  intent(in)  :: width ! integration interval
real*8,  intent(in)  :: mass
real*8,  intent(in)  :: deriv(n)
real*8,  intent(out) :: KEmat(n*(n+1)/2)

real*8, parameter :: zero=0.0d+00, one=1.0d+00, two=2.0d+00
real*8, parameter :: hbar=one

integer :: nn, i, j, currelem
real*8  :: pi
real*8  :: delk, tempre, tempim
real*8  :: pmom(n), fourwf1(2*n), fourwf2(2*n)

pi=two*two*datan(one)
delk = two*pi/width
nn = n/2
do i=1,nn
   pmom(i)=(i-1)*delk*hbar           ! p*|k_i> = k_i*hbar*|k_i> where k_i=(i-1)*delk
   pmom(i+nn)=(i-1-nn)*delk*hbar
enddo

currelem=1 ! Index for storing in KEmat (lower triangle of kinetic energy matrix)

do i=1,n

! Set column vector e_i
! Second contribution scaled by 1/J^2 prior to FFT
  do j=1,n
    if (j == i) then
     fourwf1(2*j-1)=one          ! Real Part
     fourwf2(2*j-1)=one/(deriv(j)*deriv(j))
    else
     fourwf1(2*j-1)=zero
     fourwf2(2*j-1)=zero
    end if
    fourwf1(2*j)=zero            ! Imaginary Part
    fourwf2(2*j)=zero
  end do

! FFT to momentum DVR
  call four1(fourwf1,n,1)
  call four1(fourwf2,n,1)

! Apply square of momentum operator = pmom(k) * pmom(k)
  do j=1,n
    fourwf1(2*j-1)=pmom(j)*pmom(j)*fourwf1(2*j-1)        ! Real Part
    fourwf2(2*j-1)=pmom(j)*pmom(j)*fourwf2(2*j-1)
    fourwf1(2*j)=pmom(j)*pmom(j)*fourwf1(2*j)            ! Imaginary Part
    fourwf2(2*j)=pmom(j)*pmom(j)*fourwf2(2*j)
  enddo

! FFT to coordinate DVR
  call four1(fourwf1,n,-1)
  call four1(fourwf2,n,-1)
  fourwf1=fourwf1/n
  fourwf2=fourwf2/n

! First contribution scaled by 1/J^2 after FFT
  do j=1,n
    fourwf1(2*j-1)=fourwf1(2*j-1)/(deriv(j)*deriv(j))
    fourwf1(2*j)=fourwf1(2*j)/(deriv(j)*deriv(j))
  end do

  do j=i,n
    KEmat(currelem)=(fourwf1(2*j-1)+fourwf2(2*j-1))/(two*two*mass)
    currelem=currelem+1
  end do

end do

return
end subroutine mappedgridke


subroutine FOUR1(DATA,NN,ISIGN)
!---1D FFT routine (Numerical Recipes)

implicit none

integer, intent(in) :: nn, isign
real*8,  intent(inout), dimension(2*nn) :: data

integer :: n, j, i, m, mmax, istep
real*8 :: WR, WI, WPR, WPI, WTEMP, THETA, TEMPR, TEMPI

N = 2*NN
J = 1

DO I=1,N,2
   IF (J.GT.I) THEN
      TEMPR = DATA(J)
      TEMPI = DATA(J+1)
      DATA(J) = DATA(I)
      DATA(J+1) = DATA(I+1)
      DATA(I) = TEMPR
      DATA(I+1) = TEMPI
   ENDIF
   M = N/2
   do while ((M.GE.2).AND.(J.GT.M))
      J = J - M
      M = M/2
   enddo
   J = J + M
enddo

MMAX = 2

do while (N.GT.MMAX)
   ISTEP = 2*MMAX
   THETA = 6.28318530717959D0/(ISIGN*MMAX)
   WPR = -2.D0*DSIN(0.5D0*THETA)**2
   WPI = DSIN(THETA)
   WR = 1.D0
   WI = 0.D0
   DO M=1,MMAX,2
      DO I=M,N,ISTEP
         J = I + MMAX
         TEMPR = WR*DATA(J) - WI*DATA(J+1)
         TEMPI = WR*DATA(J+1) + WI*DATA(J)
         DATA(J) = DATA(I) - TEMPR
         DATA(J+1) = DATA(I+1) - TEMPI
         DATA(I) = DATA(I) + TEMPR
         DATA(I+1) = DATA(I+1) + TEMPI
      enddo
      WTEMP = WR
      WR = WR*WPR - WI*WPI + WR
      WI = WI*WPR + WTEMP*WPI + WI
   enddo
   MMAX = ISTEP
enddo

return
end subroutine four1

