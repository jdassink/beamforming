!!
!!  Auxiliary Numerical recipes routines 
!!

  module sa_math
  implicit none

  double precision, parameter :: pi = 3.1415926536
  double complex, parameter :: im = (0.,1.)

  contains

  subroutine jacobi(A,N,D,V,NROT)
  ! Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986
   integer ip, iq, i, j
   integer N,NROT
   real*8  A(1:N,1:N),D(1:N),V(1:N,1:N)
   real*8, pointer :: B(:), Z(:)
   real*8  c,g,h,s,sm,t,tau,theta,tresh

   allocate(B(100))
   allocate(Z(100))

  do ip=1, N    !initialize V to identity matrix
    do iq=1, N
      V(ip,iq)=0.d0
    end do
      V(ip,ip)=1.d0
  end do
  do ip=1, N
    B(ip)=A(ip,ip)
    D(ip)=B(ip)
    Z(ip)=0.d0
  end do
  NROT=0
  do i=1, 50
    sm=0.d0
    do ip=1, N-1     !sum off-diagonal elements
      do iq=ip+1, N
        sm=sm+DABS(A(ip,iq))
      end do
    end do
    if(sm==0.d0) return  !normal return
    if(i.lt.4) then
      tresh=0.2d0*sm**2
    else
      tresh=0.d0
    end if
    do ip=1, N-1
      do iq=ip+1, N
        g=100.d0*DABS(A(ip,iq))
! after 4 sweeps, skip the rotation if the off-diagonal element is small
        if((i.gt.4).and.(DABS(D(ip))+g.eq.DABS(D(ip))) &
                .and.(DABS(D(iq))+g.eq.DABS(D(iq)))) then
                  A(ip,iq)=0.d0
        else if(DABS(A(ip,iq)).gt.tresh) then
          h=D(iq)-D(ip)
          if(DABS(h)+g.eq.DABS(h)) then
            t=A(ip,iq)/h
          else
            theta=0.5d0*h/A(ip,iq)
            t=1.d0/(DABS(theta)+DSQRT(1.d0+theta**2))
            if(theta.lt.0.d0) t=-t
          end if
          c=1.d0/DSQRT(1.d0+t**2)
          s=t*c
          tau=s/(1.d0+c)
          h=t*A(ip,iq)
          Z(ip)=Z(ip)-h
          Z(iq)=Z(iq)+h
          D(ip)=D(ip)-h
          D(iq)=D(iq)+h
          A(ip,iq)=0.d0
          do j=1, ip-1
            g=A(j,ip)
            h=A(j,iq)
            A(j,ip)=g-s*(h+g*tau)
            A(j,iq)=h+s*(g-h*tau)
          end do
          do j=ip+1, iq-1
            g=A(ip,j)
            h=A(j,iq)
            A(ip,j)=g-s*(h+g*tau)
            A(j,iq)=h+s*(g-h*tau)
          end do
          do j=iq+1, N
            g=A(ip,j)
            h=A(iq,j)
            A(ip,j)=g-s*(h+g*tau)
            A(iq,j)=h+s*(g-h*tau)
          end do
          do j=1, N
            g=V(j,ip)
            h=V(j,iq)
            V(j,ip)=g-s*(h+g*tau)
            V(j,iq)=h+s*(g-h*tau)
          end do
          NROT=NROT+1
        end if !if ((i.gt.4)...
      end do !main iq loop
    end do !main ip loop
    do ip=1, N
      B(ip)=B(ip)+Z(ip)
      D(ip)=B(ip)
      Z(ip)=0.d0
    end do
  end do !main i loop
  write(6,*) ' 50 iterations !'
  return
  end subroutine


  subroutine spline (x,y,n,yp1,ypn,y2)
    !-----------------------------------------------------------------------
    ! given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
    ! y_i = f(x_i), with x_1 < x_2 < ... < x_n, and given values yp1 ypn for
    ! the first derivative of the interpolating function at points 1 and n,
    ! this routine returns an array y2(1:n) of length n which contains the
    ! second derivatives of the interpolating function at the tabulated
    ! points x_i. If yp1 and ypn are equal 1.d30 or larger, the routine
    ! is signalled to set the corresponding boundary conditions for a
    ! natural spline, with zero second derivative on that boundary.
    ! (c) modified from Numerical recipes
    !-----------------------------------------------------------------------
    integer                                    :: n,i,k
    integer, parameter                         :: idouble = kind(1.0d0)
    integer, parameter                         :: isingle = kind(1.0)
    real(idouble), dimension(n), intent(in)    :: x,y
    real(idouble), dimension(n), intent(out)   :: y2
    real(idouble), dimension(n)                :: u
    real(idouble), intent(in)                  :: yp1,ypn
    real(idouble)                              :: p,qn,sig,un

    if (yp1.gt..99e30) then
      y2(1) = 0.
      u(1)  = 0.
    else
      y2(1) = -0.5
      u(1)  = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i = 2,n-1
      sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p     = sig*y2(i-1)+2.
      y2(i) = (sig-1.)/p
      u(i)  = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))             &
&          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (ypn.gt..99e30) then
      qn = 0.
      un = 0.
    else
      qn = 0.5
      un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do k = n-1,1,-1
      y2(k) = y2(k)*y2(k+1)+u(k)
    enddo
    return
  end subroutine spline


  function splint (xa,ya,y2a,n,x)
!-----------------------------------------------------------------------
! given the arrays xa(1:n) and ya(1:n) of length n, which tabulate
! a function (with the xa's in ascending order), and given the array
! y2(1:n), which is the output  from SPLINE, and given a value of x,
! the routine returns a cubic-spline interpolated value of y.
! (c) modified from Numerical recipes
!-----------------------------------------------------------------------
    integer                                    :: n,k,khi,klo
    integer, parameter                         :: idouble = kind(0.0d0)
    integer, parameter                         :: isingle = kind(0.0)
    real(idouble), intent(in)                  :: x
    real(idouble), dimension(n), intent(in)    :: xa,ya,y2a
    real(idouble)                              :: a,b,h,splint

    klo = 1
    khi = n
1   if (khi-klo .gt. 1) then
      k=int(dfloat(khi+klo)/2.d0)
      if (xa(k) .gt. x) then
        khi = k
      else
        klo = k
      endif
      goto 1
    endif
    h = xa(khi) - xa(klo)

    if (h.eq.0.) stop 'bad xa input in splint'
    
    a      = (xa(khi)-x)/h
    b      = (x-xa(klo))/h
    splint = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)                     &
 &    + (b**3-b)*y2a(khi))*(h**2)/6.d0
    
    return
  end function splint

end module
