program gridnacme

! Inputs:
!    potential.in - adiabatic electronic ground state potential in eV on scaled grid in bohr
!
! Outputs:
!    energies.out - state energies in eV
!    wf????.out   - min(ngrid,2000) number of wavefunctions printed out on grid (bohr)
!
! Features:
!    Normalizes wavefunctions
!    Requires LAPACK for digonalization and up to BLAS level 2

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
character*4         :: istring
real*8, parameter   :: zero=0.0d+00, one=1.0d+00, two=2.0d+00
real*8, parameter   :: bohr2ang=0.529177249d+00, ang2bohr=one/bohr2ang
real*8, parameter   :: au2ev=27.21138505d+00, ev2au=one/au2ev
real*8, parameter   :: au2kcal=627.5095d+00, kcal2au=one/au2kcal
real*8, parameter   :: au2cm=219474.6, cm2au=one/au2cm
real*8, parameter   :: normtol=1.0d-10
real*8, allocatable :: coords(:), scoords(:), sderiv(:,:)
real*8, allocatable :: pot(:), en(:)
real*8, allocatable :: wf(:,:)

! Read in electronic potential in au
open(unit=301,file='potential.in')
read(301,*) totwidth                               ! Size of integration interval[bohr]
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
  coords(k)=coords(k)
  pot(k)=pot(k)*ev2au    ! schroedinger.f90 requires au
end do
close(301)
width=totwidth/dble(ngrid-1)                    ! grid spacing[bohr]

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
allocate (en(ngrid),wf(ngrid,ngrid))  ! Energies(quantum_num)[au], Wavefunctions(grid,quantum_num)

call mappedschroed(ngrid,stotwidth,mass,pot,sderiv,en,wf)

! Write energies
open(unit=303,file="energies.out")
write(303,'(A1,1X,A10,A20)') "#", "State", "En(eV)"
do k=1,ngrid
  write(303,'(I12,G20.12)') k, en(k)*au2ev
end do
close(303)

! Process and print out a maximum of 2000 wavefunctions
nwf=min(2000,ngrid)

! Normalize wavefunctions
do j=1,nwf
  wfnrm=dnrm2(ngrid,wf(:,j),1)
  if(wfnrm > normtol) then
   call dscal(ngrid,one/wfnrm,wf(:,j),1)
  else
   write(*,*) "Wavefunction norm is too small to scale (i,<wf_i|wf_i>):",j,wfnrm
  end if
end do 

! Write first nwf wavefunctions
do i=1,nwf
  write(istring,'(I4.4)') i
  open(unit=303,file='wf'//istring//'.out')
  write(303,'(A1,2(A20))') "#","Coords (bohr)","Wavefunction"
  do k=1,ngrid
    write(303,'(1X,2(G20.12))') coords(k),wf(k,i)
  end do
  close(303)
end do

deallocate(sderiv,scoords,coords)
deallocate(pot,wf,en)

end program gridnacme

