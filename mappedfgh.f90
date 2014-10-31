program gridnacme

! Inputs:
!    potential.in - (see schroedinger.f90) adiabatic electronic ground state potential in au on scaled grid
!
! Outputs:
!    DA-energies.out - double adiabatic vibronic energies (assumes no nonadiabatic effects)
!    hamiltonian_nog.dat - vibronic Hamiltonian in double adiabatic product basis neglecting second-order terms
!    hamiltonian.dat - vibronic Hamiltonian in double adiabatic product basis
!    protwf1.out - four lowest-energy proton vibrational states moving in adiabatic electronic potential 1
!    protwf2.out - four lowest-energy proton vibrational states moving in adiabatic electronic potential 2
!    energies.out - 'exact' vibronic state energies in kcal/mol and cm^-1
!    wf1.out - expansion coefficients for ground vibronic state
!    wf2.out - expansion coefficients for first excited vibronic state
!
! Features:
!    Assumes first order nonadiabatic coupling d_ik for i\in{1,2} and k>2 are negligible
!    Normalizes and phases proton vibrational wavefunctions

use cubicspline
implicit none

interface
  subroutine mappedschroed(npnts,xint,mass,potential,coordderiv,energy,wavefunction)
  integer, intent(in)                          :: npnts
  real*8,  intent(in)                          :: xint
  real*8,  intent(in)                          :: mass
  real*8,  intent(in),  dimension(npnts)       :: potential
  real*8,  intent(in),  dimension(3,npnts)     :: coordderiv
  real*8,  intent(out), dimension(npnts)       :: energy
  real*8,  intent(out), dimension(npnts,npnts) :: wavefunction
  end subroutine mappedschroed
end interface

integer             :: ngrid
integer             :: nwf
integer             :: i, j, k, l, m, n, ierr, repcount
real*8              :: totwidth, width, stotwidth, swidth
real*8              :: mass, wfnrm
real*8              :: dnrm2
character*3         :: istring
real*8, parameter   :: zero=0.0d+00, one=1.0d+00, two=2.0d+00
real*8, parameter   :: bohr2ang=0.529177249d+00, ang2bohr=one/bohr2ang, tolerance=1.0d-05
real*8, parameter   :: au2kcal=627.5095d+00, kcal2au=one/au2kcal, au2cm=219474.6, cm2au=one/au2cm
real*8, allocatable :: coords(:), scoords(:), sderiv(:,:)
real*8, allocatable :: pot(:), en(:)
real*8, allocatable :: wf(:,:)

! Read in electronic potential in au
open(unit=301,file='potential.in')
read(301,*) totwidth                               ! Size of integration interval[ang]
read(301,*) mass                                   ! Mass of quantum particle[emu]
read(301,*) ngrid                                  ! Number of grid points
if (mod(ngrid,2) /= 0) then
 write(*,*) "Number of grid points,",ngrid," is not even"
 STOP
end if
allocate (coords(ngrid), pot(ngrid))       ! Coords(grid)[bohr], Potential(grid)[au]
coords=zero
pot=zero
do k=1,ngrid
  read(301,*) coords(k), pot(k)
  coords(k)=coords(k)*ang2bohr
  pot(k)=pot(k)
end do
close(301)
width=totwidth/dble(ngrid-1)*ang2bohr                    ! grid spacing[bohr]

! Define the scaled coords to run uniformly along a grid [-1,1]au
allocate (scoords(ngrid))
stotwidth=two
swidth=stotwidth/dble(ngrid-1)
do i=1,ngrid/2
  scoords(i)=-stotwidth/two+swidth*(i-1)
  scoords(ngrid-i+1)=stotwidth/two-swidth*(i-1)
end do

! Spline coords mapping and calculate necessary derivative (unitless since both coords given in au)
call splineinit(ngrid,scoords,coords)
call writespline(16*ngrid)

allocate(sderiv(3,ngrid)) ! (order of deriv,grid pt)
do i=1,ngrid
  sderiv(1,i)=splinefuncderiv(1,scoords(i))
  sderiv(2,i)=splinefuncderiv(2,scoords(i))
  sderiv(3,i)=splinefuncderiv(3,scoords(i))
end do

call cleanspline

! Calculate wavefunctions
allocate (en(ngrid),wf(ngrid,ngrid))  ! Energies(quantum_num)[kcal/mol], Wavefunctions(grid,quantum_num)

call mappedschroed(ngrid,stotwidth,mass,pot,sderiv,en,wf)

! Write energies
open(unit=303,file="energies.out")
write(303,'(A12,2(A20))') "State", "En(kcal/mol)", "En(cm^(-1))"
do k=1,ngrid
  write(303,'(I12,2(G20.12))') k, en(k), en(k)*kcal2au*au2cm
end do
close(303)

! Process and print out a maximum of 200 wavefunctions
nwf=min(200,ngrid)

! Normalize wavefunctions
do j=1,nwf
  wfnrm=dnrm2(ngrid,wf(:,j),1)
  call dscal(ngrid,one/wfnrm,wf(:,j),1)
end do 

! Write first nwf wavefunctions
do i=1,nwf
  write(istring,'(I3.3)') i
  open(unit=303,file='wf'//istring//'.out')
  write(303,'(A1,2(A20))') "#","Coords (Ang)","Wavefunction"
  do k=1,ngrid
    write(303,'(1X,2(G20.12))') coords(k)*bohr2ang,wf(k,i)
  end do
  close(303)
end do

deallocate(sderiv,scoords,coords)
deallocate(pot,wf,en)

end program gridnacme

