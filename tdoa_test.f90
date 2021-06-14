  ! Time-difference of Arrival (TDOA) analysis
  ! based on "Atmospheric infrasound parameter estimation", by Szuberla and Olson, JASA 2003

  ! The slowness vector is estimated using a least-squares estimater.

  ! Furthermore, the Fisher ratio and the Q-factor is computed, discriminating between
  ! near-field and far-field events
  ! Also, uncertainty values in bearing and velocity are computed from the error ellipse

  ! Author     : Jelle D. Assink
  ! Date       : November 2010
 
  module mod_tdoa_test
  use sa_math
  implicit none
  
  contains 

  subroutine convert_slowness(px,py,bearing,trcvel)
    double precision px, py, bearing, trcvel

    bearing = modulo((180./pi)*Atan2(px, py),360.)
    trcvel  = 1./sqrt(px**2+py**2)
  end subroutine


  subroutine get_co_array(r,n_instr,coarray,n_pairs)
   integer i, j, n_pairs, n_instr, dummy
   double precision, allocatable, dimension(:,:) :: r, coarray

   n_pairs = 0
   dummy = 0
   do i=1,n_instr
    do j=i+1,n_instr
     n_pairs = n_pairs + 1
     coarray(n_pairs,1) = r(i,1)-r(j,1)
     coarray(n_pairs,2) = r(i,2)-r(j,2)
     if (i == 3 .OR. j == 3) then
      coarray(n_pairs,1) = 0.
      coarray(n_pairs,2) = 0.
      dummy = dummy + 1
     endif
    end do
   end do
   n_pairs = n_pairs - dummy
  end subroutine

  subroutine get_XtX_matrices(coarray,n_pairs,XtX,XtXinv)
   integer n_pairs, i
   double precision, allocatable, dimension(:,:) :: coarray, XtX, XtXinv
   double precision determinant

   XtX(1,1) = 0.
   XtX(2,2) = 0.
   XtX(1,2) = 0.
   XtX(2,1) = 0.
   do i=1,n_pairs
    XtX(1,1) = XtX(1,1) + coarray(i,1)*coarray(i,1)
    XtX(1,2) = XtX(1,2) + coarray(i,1)*coarray(i,2)
    XtX(2,1) = XtX(2,1) + coarray(i,2)*coarray(i,1)
    XtX(2,2) = XtX(2,2) + coarray(i,2)*coarray(i,2)
   end do

   determinant = XtX(1,1)*XtX(2,2)-XtX(2,1)*XtX(1,2)
   XtXinv(1,1) = (1./determinant)*XtX(2,2)
   XtXinv(2,2) = (1./determinant)*XtX(1,1)
   XtXinv(1,2) = (-1./determinant)*XtX(2,1)
   XtXinv(2,1) = (-1./determinant)*XtX(1,2)
  end subroutine

  subroutine estimate_slowness(XtXinv,coarray,n_pairs,dt,px,py,error)
   integer n_pairs, j
   double precision px, py, error, Xf
   double precision, allocatable, dimension(:,:) :: XtXinv, coarray
   double precision, allocatable, dimension(:)   :: dt, XtTau

   allocate(XtTau(2))
   j = 1
   XtTau(1) = 0.
   XtTau(2) = 0.
   do while (j < n_pairs+1)
    XtTau(1) = XtTau(1) + coarray(j,1)*dt(j)
    XtTau(2) = XtTau(2) + coarray(j,2)*dt(j)
    j = j + 1
   end do 
   px = XtXinv(1,1)*XtTau(1) + XtXinv(1,2)*XtTau(2)
   py = XtXinv(2,1)*XtTau(1) + XtXinv(2,2)*XtTau(2)
   if (px .eq. 0) then
    px = 1.0E-5
   endif
   if (py .eq. 0) then
    py = 1.0E-5
   endif

   error = 0
   j = 1
   do while (j < n_pairs+1)
    Xf    = coarray(j,1)*px + coarray(j,2)*py
    error = error + (dt(j)-Xf)**2
    j     = j + 1
   end do
   deallocate(XtTau) 
  end subroutine
  
  subroutine estimate_sigma_Q(n_pairs,coarray,XtXinv,error,sigma,Q)
   integer i, j, n_pairs
   double precision error, sigma, Q, trc_R
   double precision, allocatable, dimension(:,:) :: coarray, XtXinv, Rmat

   allocate(Rmat(n_pairs,n_pairs))
   trc_R = 0
   do i=1,n_pairs
    do j=1,n_pairs
     Rmat(i,j) = (coarray(i,1)*XtXinv(1,1)+coarray(i,2)*XtXinv(2,1))*coarray(j,1) + &
  &   (coarray(i,1)*XtXinv(1,2)+coarray(i,2)*XtXinv(2,2))*coarray(j,2)
     if (i == j) then
      trc_R = trc_R + Rmat(i,j) 
     endif
    end do
   end do

   sigma = sqrt(error/(n_pairs-trc_R))
   Q     = sigma*20.0
   deallocate(Rmat)
  end subroutine

  subroutine compute_uncertainty_params(XtXinv,lambda_x,lambda_y,angle)
   integer nrot
   double precision lambda_x,lambda_y, angle
   double precision, allocatable, dimension(:,:) :: a, XtXinv, v
   double precision, allocatable, dimension(:) :: d
 
   allocate(d(2)) 
   allocate(v(2,2))
   allocate(a(2,2))
   a = XtXinv
   call jacobi(a,2,d,v,nrot)
   lambda_y = minval(d,1)
   lambda_x = maxval(d,1)
   angle    = (180./pi)*Atan(v(2,maxloc(d,1))/v(1,maxloc(d,1)))
   deallocate(v)
   deallocate(d)
   deallocate(a)
  end subroutine

  subroutine compute_uncertainties(sigma,lambda_x,lambda_y,angle,px,py,vel_dev,bear_dev)
   integer i, j, n, imax, jmax
   double precision lambda_x, lambda_y, angle, px, py, vel_dev, bear_dev, sigma
   double precision dchi2, term1, term2, px_min, px_max, tmp, dx, dy
   double precision dmax, phmax, d, ph
   double precision, allocatable, dimension (:) :: px_, py_

   dchi2 = 5.9915
   term1 = sigma**2/lambda_y
   term2 = lambda_x/lambda_y
   n     = 20000;

   px_min  = px - px;
   px_max  = px + px;

   allocate(px_(2*n))
   allocate(py_(2*n))
   j = 0
   do i=1,2*n
    if (i<=n) then
     tmp = px_min + (px_max-px_min)*i/n
     if (term1*dchi2-term2*(tmp-px)**2 >= 0.) then
      j      = j+1
      px_(j) = tmp
      py_(j) = py + sqrt(term1*dchi2-term2*(tmp-px)**2)
     endif
    else
     tmp = px_min + (px_max-px_min)*(i-n)/n
     if (term1*dchi2-term2*(tmp-px)**2 >= 0.) then
      j      = j+1
      px_(j) = tmp
      py_(j) = py -sqrt(term1*dchi2-term2*(tmp-px)**2)
     endif
    endif
   end do

   do i=1,j
    dx = px + cos((pi/180.)*angle)*(px_(i)-px)-sin((pi/180.)*angle)*(py_(i)-py)
    dy = py + sin((pi/180.)*angle)*(px_(i)-px)+cos((pi/180.)*angle)*(py_(i)-py)
    px_(i) = dx
    py_(i) = dy
   end do

   dmax  = 0;
   phmax = 0;
   imax  = 1;
   jmax  = 1;
   do i=1,j
    d  = sqrt((py_(i)-py)**2+(px_(i)-px)**2)
    ph = abs((180./pi)*atan2(px_(i),py_(i))-(180./pi)*atan2(px,py))
    if (d > dmax) then
     imax = i
     dmax = d
    endif
    if (ph > phmax) then
        jmax = i
        phmax = ph
    endif
   enddo

   vel_dev  = abs(1./sqrt(px_(imax)**2+py_(imax)**2)-1./sqrt(px**2+py**2))
   bear_dev = abs((180./pi)*atan2(px_(jmax),py_(jmax))-(180./pi)*atan2(px,py))
  end subroutine


  subroutine printUsage()
    write(6,*) ' '
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    write(6,*) 'Example time-difference of arrival (TDOA) analysis'
    write(6,*) ''
    write(6,*) ' Runs the Olson - Inframatics, 2004 case study'
    write(6,*) ''
    write(6,*) 'Usage: tdoa_test'
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    call exit(0)
  end subroutine

  end module mod_tdoa_test

!--------------------------------------------------------------------------------------
!-------------------------------START OF MAIN PROGRAM----------------------------------
!--------------------------------------------------------------------------------------
  program tdoa_test
    use mod_tdoa_test
    implicit none

    integer   i, j, n_instr, n_pairs, n_pairs_rdx
    double precision px, py, bearing, velocity
    double precision error, sigma, Q, lambda_x, lambda_y, angle
    double precision vel_dev, bear_dev

    double precision, allocatable, dimension(:,:) :: r
    double precision, allocatable, dimension(:,:) :: coarray, XtX, XtXinv
    double precision, allocatable, dimension(:)   :: dt

    n_instr = 8
    allocate(r(n_instr,2))

!   Test the case study mentioned in the Inframatics paper

    r(1:n_instr,1) = (/ 0.0000, 0.9546, 0.5671, -0.5718, -0.9565, 0.0003, 0.0882, -0.0895 /)
    r(1:n_instr,2) = (/ 0.0000,-0.6857,-1.8072, -1.8152, -0.6757,-0.8947,-1.0487, -1.0500 /)
    r = r*1000

    n_pairs = 0
    do i=1,n_instr
     do j=i+1,n_instr
      n_pairs = n_pairs + 1
     end do
    end do

    allocate(dt(n_pairs))
    allocate(coarray(n_pairs,2))
    allocate(XtX(2,2))
    allocate(XtXinv(2,2))

    write(6,*) ''
    write(6,*) ''
    write(6,*) ' Time-delay of arrival (TDOA) example analysis ...'
    write(6,*) '-------------------------------------------------'
    write(6,*) ''

    call get_co_array(r,n_instr,coarray,n_pairs_rdx)
    dt = (/ 3.30, 0.0, 4.50, 0.70, 2.55, 3.15, 2.90, 0.0, 1.20, -2.55, &
&          -0.75, -0.15, -0.40, 0.0, 0.0, 0.0, 0.0, 0.0, -3.75, -1.90, &
&          -1.35, -1.60,  1.85, 2.40, 2.15, 0.55, 0.35, -0.25 /)


    call get_XtX_matrices(coarray,n_pairs,XtX,XtXinv)
    call compute_uncertainty_params(XtX,lambda_x,lambda_y,angle)

    call estimate_slowness(XtXinv,coarray,n_pairs,dt,px,py,error)   
    call estimate_sigma_Q(n_pairs_rdx,coarray,XtXinv,error,sigma,Q)
    call compute_uncertainties(sigma,lambda_x,lambda_y,angle,px,py,vel_dev,bear_dev)
    call convert_slowness(px,py,bearing,velocity)

    write(6,'(a24,f8.2,a4)') ' Back-azimuth         : ', bearing,  ' deg'
    write(6,'(a24,f8.2,a4)') ' Azimuth error        : ', bear_dev, ' deg'
    write(6,'(a24,f8.2,a4)') ' Trace velocity       : ', velocity, ' m/s'
    write(6,'(a24,f8.2,a4)') ' Trace velocity error : ', vel_dev,  ' m/s'
    write(6,'(a42)'        ) ' ---------------------------------------- '
    write(6,'(a24,e9.2,a4)') ' Least-squares error  : ', error   , ''
    write(6,'(a24,f8.5,a4)') ' Sigma tau            : ', sigma   , '   s'
    write(6,'(a24,f8.2,a4)') ' Q factor             : ', Q       , ''

    deallocate(r)
    deallocate(dt)
    deallocate(XtX)
    deallocate(XtXinv)
    deallocate(coarray)
  end program tdoa_test
