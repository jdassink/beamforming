  ! Time-difference of Arrival (TDOA) analysis
  ! based on "Atmospheric infrasound parameter estimation", by Szuberla and Olson, JASA 2003

  ! The slowness vector is estimated using a least-squares estimater.

  ! Furthermore, the Fisher ratio and the Q-factor is computed, discriminating between
  ! near-field and far-field events
  ! Also, uncertainty values in bearing and velocity are computed from the error ellipse

  ! Author     : Jelle D. Assink
  ! Date       : November 2010
 
  module mod_tdoa
  use sa_array

  contains 

  subroutine get_timelags(counts_bin,binsize,fbinsize,srate,n_instr,dt,correlation)
   implicit none
   integer i, j, n_instr, n_pairs, binsize, fbinsize, ndim
   double precision c_max, srate, t_max
   double precision, allocatable, dimension (:)   :: dt, correlation
   double precision, allocatable, dimension (:)   :: tau, c_wrap
   double precision, allocatable, dimension (:,:) :: counts_bin, s1, s2

   ndim = 1
   allocate(s1(binsize,ndim))
   allocate(s2(binsize,ndim))
   allocate(tau(fbinsize))
   allocate(c_wrap(fbinsize))

   n_pairs = 0
   do i = 1,n_instr
    do j = i+1,n_instr
     s1(1:binsize,ndim) = counts_bin(1:binsize,i)
     s2(1:binsize,ndim) = counts_bin(1:binsize,j)

     call xcorr(s1,s2,ndim,binsize,fbinsize,srate,tau,c_wrap,t_max,c_max)     
     
     n_pairs = n_pairs + 1
     dt(n_pairs) = -1.0*t_max
     correlation(n_pairs)= c_max
    end do 
   end do

   deallocate(s1)
   deallocate(s2)
   deallocate(tau)
   deallocate(c_wrap)
  end subroutine

  subroutine estimate_slowness(XtXinv,coarray,n_pairs,dt,px,py,error)
   implicit none
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
  
  subroutine estimate_sigma_Q(n_pairs,coarray,XtXinv,srate,error,sigma,Q)
   implicit none
   integer i, j, n_pairs
   double precision srate, error, sigma, Q, trc_R
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
   Q     = sigma*srate
   deallocate(Rmat)
  end subroutine

  subroutine compute_uncertainty_params(XtX,pi,lambda_x,lambda_y,angle)
   implicit none
   integer nrot
   double precision lambda_x,lambda_y, angle, pi
   double precision, allocatable, dimension(:,:) :: a, XtX, v
   double precision, allocatable, dimension(:) :: d
 
   allocate(d(2)) 
   allocate(v(2,2))
   allocate(a(2,2))
   a = XtX
   call jacobi(a,2,d,v,nrot)
   lambda_y = minval(d,1)
   lambda_x = maxval(d,1)
   angle    = (180./pi)*Atan(v(2,maxloc(d,1))/v(1,maxloc(d,1)))
   deallocate(v)
   deallocate(d)
   deallocate(a)
  end subroutine

  subroutine compute_uncertainties(pi,sigma,lambda_x,lambda_y,angle,px,py,vel_dev,bear_dev)
   implicit none
   integer i, j, n, imax, jmax
   double precision lambda_x, lambda_y, angle, px, py, vel_dev, bear_dev, sigma
   double precision dchi2, term1, term2, px_min, px_max, tmp, pi, dx, dy
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

   deallocate(px_)
   deallocate(py_)
  end subroutine


  subroutine get_cmd_parameters(n_instr,tbinsize,overlap,oversampling,r,file_format,timeseries)
    implicit none
    integer i, n_instr, eof, overlap, iargc, oversampling
    double precision tbinsize
    character(40), allocatable, dimension (:) :: timeseries
    character(40) coordinates, file_format, dummy
    double precision, allocatable, dimension(:,:) :: r

    call getarg(1, dummy)
    read(dummy,*,IOSTAT=eof) coordinates
    if (eof < 0) then
      call printUsage()
    endif
    call get_coordinates(coordinates,n_instr,r)

    if (iargc() .ne. 5 + n_instr) then
      deallocate(r)
      deallocate(timeseries)
      call printUsage()
    endif              
    call getarg(2, dummy)
    read(dummy,*) tbinsize
    call getarg(3, dummy)
    read(dummy,*) overlap
    call getarg(4, dummy)
    read(dummy,*) oversampling
    call getarg(5, dummy)
    read(dummy,*) file_format
    i = 1
    do while ( i < n_instr+1 )
      call getarg(i+5, timeseries(i) )
      i = i+1
    end do

    if (overlap > 99) then
      write(6,*)  ' '
      write(6,*)  'ERROR: Overlap size cannot exceed 99%'
      write(6,*)  ' - Exiting ...'
      write(6,*)  ' '
      call exit(0)
    endif
  end subroutine


  subroutine printUsage()
    write(6,*) ' '
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    write(6,*) 'Time-difference of arrival (TDOA) analysis for ASCII/SAC formatted data'
    write(6,*) '________________________________________________________________________'
    write(6,*) ''
    write(6,*) 'Usage:'
    write(6,*) '------'
    write(6,*) 'tdoa <coordinate file> <binsize> <overlap> <oversampling> <ascii|sac> <files>'
    write(6,*) ' '
    write(6,*) ' - Order of coordinates must match the order of input files'
    write(6,*) ' - Sample rate must be provided in the input file'
    write(6,*) ' '
    write(6,*) ' '
    write(6,*) 'Example: '
    write(6,*) '-------- '
    write(6,*) 'tdoa stationtable 10.0 50 10 sac DBN*sac'
    write(6,*) ' '
    write(6,*) ' - Time-difference of arrival processing (max. 24 hours of data)'
    write(6,*) ' - Timebins of 10.0 s, with 50% overlap = 5.0 s'
    write(6,*) ' - Cross correlation function is oversampled using factor of 10'
    write(6,*) ' - Binary SAC input is expected, files DBN*.sac'
    write(6,*) ''
    write(6,*) ' The output file <tdoa.dat> format is:'
    write(6,*) ' '
    write(6,*) '   1. Timebin'
    write(6,*) '   2. Fisher-ratio'
    write(6,*) '   3. Back azimuth   [deg]'
    write(6,*) '   4. Trace velocity [m/s]'
    write(6,*) '   5. Slowness X     [s/m]'
    write(6,*) '   6. Slowness Y     [s/m]'
    write(6,*) '   7. RMS amplitude'
    write(6,*) '   8. # of instruments used in analysis'
    write(6,*) '   9. Center frequency  [Hz]'
    write(6,*) '  10. ----------------------------------------'
    write(6,*) '  11. Mean of cross-correlation maximum (MCCM)'
    write(6,*) '  12. Least squares error'
    write(6,*) '  13. Sigma_t        [s]'
    write(6,*) '  14. Q-value'
    write(6,*) '  15. Std. error azimuth        [deg]'
    write(6,*) '  16. Std. error trace velocity [m/s]'
    write(6,*) ' '
    write(6,*) ' High Q values are indicative of near field arrivals '
    write(6,*) ' '
    write(6,*) ' Bestbeam output is written out as SAC file <bestbeam.sac>'
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    call exit(0)
  end subroutine

  end module mod_tdoa

!--------------------------------------------------------------------------------------
!-------------------------------START OF MAIN PROGRAM----------------------------------
!--------------------------------------------------------------------------------------
  program tdoa
    use mod_tdoa
    implicit none

    integer   i, n_instr, alloc_instr, alloc_samples, n_pairs
    integer   binsize, fbinsize, bin, n_bins, overlap, start_sample, end_sample, oversampling
    double precision tbinsize, fisher, px, py, bearing, velocity, pi, mccm
    double precision error, sigma, Q, lambda_x, lambda_y, angle
    double precision vel_dev, bear_dev, prms, time, fsrate
    character(40) file_format, fid_beam

    double precision, allocatable, dimension(:,:) :: r, counts, counts_bin, bestbeam
    double precision, allocatable, dimension(:,:) :: coarray, XtX, XtXinv
    double precision, allocatable, dimension(:)   :: dt, c
    character(40), allocatable, dimension (:)     :: timeseries

    type(Frequency) :: freq
    type(WFHeader) :: header

    alloc_instr = 50
    alloc_samples = 2*24*3600*250
    pi = Acos(-1.)

    allocate(timeseries(alloc_instr))
    allocate(r(alloc_instr,3))
    call get_cmd_parameters(n_instr,tbinsize,overlap,oversampling,r,file_format,timeseries)

    allocate(counts(alloc_samples,n_instr))
    allocate(bestbeam(alloc_samples,1))
    
    write(6,*) ' '
    write(6,*) ' '
    write(6,*) 'Reading files ...'
    write(6,*) '----------------------------------'
    call get_timeseries(n_instr,file_format,timeseries,alloc_samples,header,counts)

    binsize = nint(tbinsize*header%srate)
    !fbinsize = 2**(nint(log10(1.0*oversampling*binsize)/log10(2.)))
    fbinsize = oversampling*binsize

    if (fbinsize .lt. binsize) then
     fbinsize = binsize
    endif

    fsrate = (1.0*fbinsize/binsize)*header%srate
    n_pairs = 0
    do i=n_instr-1,1,-1
     n_pairs = n_pairs + i
    end do

    allocate(counts_bin(binsize,n_instr))
    allocate(dt(n_pairs))
    allocate(c(n_pairs))
    allocate(coarray(n_pairs,2))
    allocate(XtX(2,2))
    allocate(XtXinv(2,2))

    write(6,*) ''
    write(6,*) ''
    write(6,*) 'Time-delay of arrival (TDOA) analysis ...'
    write(6,*) '-----------------------------------------'
    write(6,*) ''

    write(6,'(a32,f8.2,a4)') ' Time-domain sampling : ', header%srate, ' Hz'
    write(6,'(a32,f8.2,a4)') ' Time bin length      : ', tbinsize, ' s'
    write(6,'(a32,f8.2,a4)') ' Time bin overlap     : ', overlap/100.*tbinsize, ' s'
    write(6,'(a32,f8.2,a4)') ' Interpolation factor : ', 1.0*fbinsize/binsize, ' x'
    write(6,'(a32,f8.2,a4)') ' Processing sampling  : ', fsrate, ' Hz'

    open (unit=30, file='tdoa.dat')
    fid_beam = 'bestbeam.sac'

    call get_co_array(r,n_instr,coarray)
    call get_XtX_matrices(coarray,n_pairs,XtX,XtXinv)
    call compute_uncertainty_params(XtX,pi,lambda_x,lambda_y,angle)

    overlap = binsize - binsize*overlap/100.
    n_bins = (binsize/overlap)*(header%n_samples/binsize)

    do bin = 2, n_bins - 1
      start_sample = (bin-1)*overlap+1
      end_sample   = start_sample + (binsize-1)
      time         = header%tzero + ((start_sample + end_sample-1) / 2.0 ) / header%srate

      counts_bin(1:binsize,1:n_instr) = counts(start_sample:end_sample,1:n_instr)
      call window_data(n_instr,binsize,counts_bin)
      call get_timelags(counts_bin,binsize,fbinsize,fsrate,n_instr,dt,c)
      mccm = sum(c)/max(1,size(c))
   
      call estimate_slowness(XtXinv,coarray,n_pairs,dt,px,py,error)   
      call estimate_sigma_Q(n_pairs,coarray,XtXinv,header%srate,error,sigma,Q)
      call compute_uncertainties(pi,sigma,lambda_x,lambda_y,angle,px,py,vel_dev,bear_dev)

      call compute_f_ratio(start_sample,binsize,n_instr,header%srate,counts,r,px,py,fisher)
      call get_bestbeam(bin,start_sample,header,counts,binsize,n_instr,overlap,r,px,py,bestbeam,prms,freq)
      call convert_slowness(px,py,bearing,velocity)

      write (30,99) time, fisher, bearing, velocity, px, py, prms, n_instr, &
&             freq%center, ' ---', mccm, error, sigma, Q, bear_dev, vel_dev
99    format (f10.2, 1x, f10.2, 1x, f8.2, 1x, f10.2, 1x, e11.3, 1x, e11.3, 1x, e11.3, 1x, i2, &
&                1x, f8.3 , a4, 1x, f4.2, 1x, e11.3, 1x, e11.3, 1x, f10.2, 1x, f10.2, 1x, f10.2)
    end do

    call write_sac(fid_beam,header,bestbeam)

    close(30)

    deallocate(r)
    deallocate(dt)
    deallocate(c)
    deallocate(counts)
    deallocate(XtX)
    deallocate(XtXinv)
    deallocate(coarray)
    deallocate(bestbeam)
    deallocate(counts_bin)
    deallocate(timeseries)

  end program tdoa
