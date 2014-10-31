
module cubicspline

! Module global variables
  implicit none
  private

  integer             :: n
  real*8, allocatable :: x(:),f(:),s(:,:) ! s(i,j) gives the a_j'th coefficient of polynomial S_i

! Module functions/subroutines
  public :: splineinit,writespline,cleanspline
  public :: splinefunc,splinefuncderiv,splinefuncint

contains

  subroutine splineinit(ngrid,coords,func)
    implicit none

    integer, intent(in) :: ngrid
    real*8, dimension(:), intent(in) :: coords(1:ngrid),func(1:ngrid)
    integer :: i

    n=ngrid
    allocate(x(1:n),f(1:n),s(1:n-1,1:4))
    do i=1,n
      x(i)=coords(i)
      f(i)=func(i)
    end do

    call setspline

    return
  end subroutine splineinit


  subroutine setspline
    implicit none

    integer :: i,j
    integer :: istat
    real*8, dimension(:) :: matdiag(1:n), matsubdiag(1:n-1), D(1:n)
    real*8, parameter :: zero=0.0d+00, one=1.0d+00, two=2.0d+00, three=3.0d+00, four=4.0d+00

    ! Set a_1=f(x_i) for all S_i
    do i=1,n-1
      s(i,1)=f(i)
    end do

    ! Set up coefficient matrix for D=[D_i]^t
    matsubdiag=one
    matdiag(1)=two
    do i=2,n-1
      matdiag(i)=four
    end do
    matdiag(n)=two

    ! Set up solution vector
    D(1)=three*(f(2)-f(1))
    do i=2,n-1
      D(i)=three*(f(i+1)-f(i-1))
    end do
    D(n)=three*(f(n)-f(n-1))

    ! Solve AD=b (D initially given as b, overwritten with solution)
    istat=0
    call dptsv(n,1,matdiag,matsubdiag,D,n,istat)
    if (istat.ne.0) then
     write(*,*) "Error in dptsv"
    end if

    ! Set a_2,a_3,a_4 for all S_i
    do i=1,n-1
      s(i,2)=D(i)
      s(i,3)=three*(f(i+1)-f(i))-two*D(i)-D(i+1)
      s(i,4)=two*(f(i)-f(i+1))+D(i)+D(i+1)
    end do

    return
  end subroutine setspline


  subroutine writespline(ngrid)
    implicit none

    integer, intent(in) :: ngrid

    integer             :: i, splinengrid
    real*8              :: splinewidth
    real*8, allocatable :: splinecoords(:)

    allocate(splinecoords(1:ngrid))

    splinewidth=(x(n)-x(1))/dble(ngrid-1)
    do i=1,ngrid
      splinecoords(i)=x(1)+(i-1)*splinewidth
    end do

    open(unit=303,file="spline.dat")
    do i=1,ngrid
      write(303,*) splinecoords(i), splinefunc(splinecoords(i))
    end do
    close(303)

    deallocate(splinecoords)

    return
  end subroutine writespline


  double precision function splinefunc(arg) result(val)
    implicit none

    real*8, intent(in) :: arg
    integer :: i
    real*8  :: coord

    if ((arg.lt.x(1)).or.(arg.gt.x(n))) then
     write(*,*) "x-value requested is out of bounds: ", arg
    else
     do i=1,n-1
       if ((arg.ge.x(i)).and.(arg.le.x(i+1))) then
        coord=(arg-x(i))/(x(i+1)-x(i))
        val=s(i,1)+s(i,2)*coord+s(i,3)*coord*coord+s(i,4)*coord*coord*coord
       end if
     end do
    end if

    return
  end function splinefunc


  double precision function splinefuncderiv(ord,arg) result(val)
    implicit none

    integer, intent(in) :: ord
    real*8,  intent(in) :: arg
    integer :: i
    real*8  :: coord,jacob
    real*8, parameter :: zero=0.0d+00, one=1.0d+00, two=2.0d+00, three=3.0d+00

    val=zero

    if ((arg < x(1)).or.(arg > x(n))) then
     write(*,*) "x-value requested is out of bounds: ", arg
     return
    end if

    if ((ord < 1).or.(ord > 3)) then
     write(*,*) "derivative order requested is out of bounds: ",ord
     return
    end if

    do i=1,n-1
      if ((arg.ge.x(i)).and.(arg.le.x(i+1))) then
       coord=(arg-x(i))/(x(i+1)-x(i))
       jacob=one/(x(i+1)-x(i))
       if (ord == 1) then
        val=jacob*(s(i,2)+two*s(i,3)*coord+three*s(i,4)*coord*coord)
       else if (ord == 2) then
        val=jacob*jacob*(two*s(i,3)+two*three*s(i,4)*coord)
       else if (ord == 3) then
        val=jacob*jacob*jacob*(two*three*s(i,4))
       end if
      end if
    end do

    return
  end function splinefuncderiv


  double precision function splinefuncint(a,b) result(val)
    implicit none

    real*8, intent(in) :: a,b
    integer :: i
    integer :: lowerknotmin,lowerknotplus
    integer :: upperknotmin,upperknotplus
    integer :: currlowerknot,currupperknot
    real*8  :: lower,upper,dummylower,dummyupper,tlower,tupper
    logical :: inverted
    real*8, parameter :: zero=0.0d+00, two=2.0d+00, three=3.0d+00, four=4.0d+00

    if (a.lt.b) then
     inverted=.false.
     lower=a
     upper=b
    else
     inverted=.true.
     lower=b
     upper=a
    end if

    if ((lower.lt.x(1)).or.(upper.gt.x(n))) then
     write(*,*) "One or more integration limits is out of bounds: ", lower, upper
    else

     do i=1,n-1
       if ((lower.ge.x(i)).and.(lower.le.x(i+1))) then
        lowerknotmin=i
        lowerknotplus=i+1
       end if
       if ((upper.ge.x(i)).and.(upper.le.x(i+1))) then
        upperknotmin=i
        upperknotplus=i+1
       end if
     end do

     val=zero
     do i=lowerknotmin,upperknotmin
       currlowerknot=i
       currupperknot=i+1
       if (currlowerknot.eq.lowerknotmin) then
        dummylower=lower
       else
        dummylower=x(currlowerknot)
       end if
       if (currupperknot.eq.upperknotplus) then
        dummyupper=upper
       else
        dummyupper=x(currupperknot)
       end if
       tlower=(dummylower-x(currlowerknot))/(x(currupperknot)-x(currlowerknot))
       tupper=(dummyupper-x(currlowerknot))/(x(currupperknot)-x(currlowerknot))
       val=val+(x(currupperknot)-x(currlowerknot))*(s(i,1)*(tupper-tlower)+s(i,2)/two*(tupper**2-tlower**2)+s(i,3)/three*(tupper**3-tlower**3)+s(i,4)/four*(tupper**4-tlower**4))
     end do

     if (inverted) then
      val=-val
     end if

    end if

    return
  end function splinefuncint


  subroutine cleanspline
    implicit none

    deallocate(s,f,x)

    return
  end subroutine cleanspline

end module cubicspline

