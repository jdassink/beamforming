  ! Basic beamforming routines, shared by timefisher and freqfisher
 
  module sa_array
  use sa_io
  use sa_fourier
  implicit none
  
  type :: Slowness
   double precision dc, dth
   double precision bearing_min, bearing_max, trcvel_min, trcvel_max
   double precision, allocatable, dimension (:) :: bearing, trcvel, px, py
   integer, allocatable, dimension (:) :: bi, ci
   integer n_beams, n_c, n_b
   character(40) grid_type, tele_local, max_all
  end type

  contains

  subroutine get_coordinates(fid,n_instr,r)
    character(40) fid, dummy
    integer n_instr
    double precision, allocatable, dimension(:,:) :: r
    double precision elev, edepth
    LOGICAL :: file_exists

    inquire(file=fid, EXIST=file_exists)

    if (file_exists) then
      open (unit=11, file=fid)
      n_instr = 0
      do while (.true.)
        n_instr = n_instr + 1
        read (11,*, end=991) dummy, r(n_instr,1), r(n_instr,2), elev, edepth, r(n_instr,4), r(n_instr,5)
        r(n_instr,3) = elev + edepth
      end do
      991 continue
      n_instr = n_instr - 1
      close(unit=11)
    else
      write(6,*) 'ERROR: station coordinate file not found. Exiting ...'
      write(6,*) ''
      call exit(1)
    endif
    ! Take first element of array as beam reference
    r(1:n_instr,1) = r(1:n_instr,1) - r(1,1)
    r(1:n_instr,2) = r(1:n_instr,2) - r(1,2)    
    ! Take middle point of array as beam reference
    !r(1:n_instr,1) = r(1:n_instr,1) - sum(r(1:n_instr,1))/n_instr
    !r(1:n_instr,2) = r(1:n_instr,2) - sum(r(1:n_instr,2))/n_instr
  end subroutine


  subroutine get_co_array(r,n_instr,coarray)
   integer i, j, n_pairs, n_instr
   double precision, allocatable, dimension(:,:) :: r, coarray

   n_pairs = 0
   do i=1,n_instr
    do j=i+1,n_instr
     n_pairs = n_pairs + 1 
     coarray(n_pairs,1) = r(i,1)-r(j,1)
     coarray(n_pairs,2) = r(i,2)-r(j,2)
    end do
   end do
  end subroutine


  subroutine convert_slowness(px,py,bearing,trcvel)
    double precision px, py, bearing, trcvel

    bearing = modulo((180./pi)*Atan2(px, py),360.)
    trcvel  = 1./sqrt(px**2+py**2)
  end subroutine


  subroutine convert_to_slowness(bearing,trcvel,px,py)
    double precision px, py, bearing, trcvel

    px = sin((pi/180.0)*bearing) / trcvel
    py = cos((pi/180.0)*bearing) / trcvel
  end subroutine


  subroutine span_slowness_grid(s_grid)
    integer i, trcvel_grid, linear_grid, n_c_local
    double precision alpha, delta, linear_max, c_tr, px, py
    type(Slowness) :: s_grid

    s_grid%grid_type = 'cylindrical'

    if (s_grid%grid_type == 'cylindrical') then
      linear_max  = 450.0

      if (s_grid%tele_local == 'tele') then
        trcvel_grid = nint((s_grid%trcvel_max  - s_grid%trcvel_min  ) / s_grid%dc ) + 1

      elseif (s_grid%tele_local == 'local') then
        linear_grid = nint((linear_max  - s_grid%trcvel_min  ) / s_grid%dc ) + 1
        delta       = s_grid%trcvel_max - linear_max
        alpha       = s_grid%dc / linear_max * 5
        n_c_local   = nint(log((delta/linear_max) + 1.) / alpha + 1)
        trcvel_grid = linear_grid + n_c_local

      else
          write (6,'(a43,a10)') 'ERROR: No support for trace velocity grid: ', s_grid%tele_local
          call exit(-1)
      endif

      s_grid%n_b     = nint((s_grid%bearing_max - s_grid%bearing_min ) / s_grid%dth) + 1
      s_grid%n_c     = trcvel_grid
      s_grid%n_beams = s_grid%n_b * s_grid%n_c

      ! generate bearing and trace velocity grids
      allocate(s_grid%bearing( s_grid%n_b ))
      allocate(s_grid%trcvel ( s_grid%n_c ))

      do i=1,trcvel_grid
        if ( (s_grid%tele_local == 'local') .and. (i .gt. linear_grid) ) then
          c_tr  = linear_max * exp(alpha*(i-linear_grid))

        else
          c_tr  = s_grid%trcvel_min + s_grid%dc*(i-1)

        endif
        s_grid%trcvel(i) = c_tr
      end do
      
      do i=1,s_grid%n_b
        s_grid%bearing(i) = s_grid%bearing_min + s_grid%dth*(i-1)
      end do


      ! generate indices to find the bearing and trace velocity from the beam number
      allocate(s_grid%bi( s_grid%n_beams ))
      allocate(s_grid%ci( s_grid%n_beams ))
      ! generate slowness values just once
      allocate(s_grid%px( s_grid%n_beams ))
      allocate(s_grid%py( s_grid%n_beams ))

      do i=1,s_grid%n_beams
        call s_grid_1d_to_2d_index(s_grid%n_b,i,s_grid%bi(i),s_grid%ci(i))
        call convert_to_slowness( s_grid%bearing(s_grid%bi(i)),   &
  &                               s_grid%trcvel( s_grid%ci(i)),   &
  &                               s_grid%px(i),  s_grid%py(i) )
      end do

    endif
  end subroutine

  subroutine s_grid_1d_to_2d_index(width_grid,k,bi,ci)
   integer k, bi, ci, width_grid

   bi = modulo((k-1),width_grid) + 1
   ci = (k-1) / width_grid + 1
  end subroutine


  subroutine beamgrid_maximum(results,k_max)
   integer k_max, tmp(1)
   double precision, allocatable, dimension (:) :: results
   
   tmp      = maxloc(results(:))
   k_max    = tmp(1)
  end subroutine


  subroutine beamgrid(statistic,start_sample,srate,binsize,s_grid,n_instr,r,counts,results)
    integer start_sample, binsize, k, n_instr
    double precision  det_stat, srate, px, py
    double precision, allocatable, dimension(:) :: results
    double precision, allocatable, dimension(:,:) :: r, counts
    character(40) statistic
    type(Slowness) :: s_grid
   
    !$omp parallel &
    !$omp& default (private) &
    !$omp& shared  (statistic,results,s_grid,start_sample,binsize,n_instr,srate,counts,r)
    !$omp do
    do k=1,s_grid%n_beams
      px = s_grid%px(k)
      py = s_grid%py(k)
      
      if (statistic .eq. 'fisher') then
        call compute_f_ratio(start_sample,binsize,n_instr,srate,counts,r,px,py,det_stat)
      
      else if (statistic .eq. 'correlation') then
        call compute_mccm(start_sample,binsize,n_instr,srate,counts,r,px,py,det_stat)
      endif

      results(k) = det_stat
    end do
    !$omp end do
    !$omp end parallel

  end subroutine


  subroutine compute_f_ratio(start_sample,binsize,n_instr,srate,counts,r,px,py,f_ratio)
    integer start_sample, binsize, c, ss_c, es_c, n_instr
    double precision f_dof, term_1, term_2, term_3, term_4, f_ratio, srate, px, py

    double precision, allocatable, dimension(:,:) :: r, counts
    double precision :: sum_c(binsize), sum_c_sq(binsize)

    f_dof = (1.*binsize*(n_instr-1))/(1.*n_instr*(binsize-1))
    sum_c(:) = 0.
    sum_c_sq(:) = 0.

    do c=1,n_instr
      ss_c     = start_sample - nint(srate*(px*r(c,1) + py*r(c,2)))
      es_c     = ss_c + (binsize-1)
      sum_c    = sum_c    + counts(ss_c:es_c,c)
      sum_c_sq = sum_c_sq + counts(ss_c:es_c,c)**2
    end do

    term_1 = sum(sum_c**2)
    term_2 = sum(sum_c)**2 / binsize
    term_3 = sum(sum_c_sq)
    term_4 = term_1 / n_instr
    f_ratio = f_dof*(term_1 - term_2)/(term_3 - term_4)

  end subroutine


  subroutine compute_mccm(start_sample,binsize,n_instr,srate,counts,r,px,py,mccm)
    integer start_sample, binsize, c, ss_c, es_c, n_instr
    double precision px, py, srate, mccm

    double precision, allocatable, dimension(:,:) :: r, counts
    double precision :: sum_c(binsize)

    sum_c(:) = 0.

    do c=1,n_instr
      ss_c     = start_sample - nint(srate*(px*r(c,1) + py*r(c,2)))
      es_c     = ss_c + (binsize-1)
      sum_c    = sum_c + counts(ss_c:es_c,c)
    end do
    sum_c = sum_c / n_instr

    mccm = maxval( sum_c )

  end subroutine


  subroutine get_bestbeam(bin,start_sample,header,counts,binsize,n_instr,overlap,r,px,py,bestbeam,prms,p2p,freq)
    integer binsize, n_instr, bin, start_sample, end_sample, overlap, c, ss_c, es_c
    integer, parameter :: ndim = 1
    double precision prms, p2p, px, py
    double precision, allocatable, dimension (:,:) :: r, counts, bb, bestbeam
    type(Frequency) :: freq
    type(WFHeader) :: header

    allocate(bb(binsize,ndim))
    
    bb(1:binsize,1:ndim) = 0.

    do c=1,n_instr
      ss_c = start_sample - nint(header%srate*(px*r(c,1) + py*r(c,2)))
      es_c = ss_c + (binsize-1)
      bb(1:binsize,1) = bb(1:binsize,1) + counts(ss_c:es_c,c)
    end do

    bb = bb / n_instr

    prms = sqrt(sum(bb(:,1)**2)/size(bb) )
    p2p = maxval(bb(:,1)) - minval(bb(:,1))

    ! Measures from Barnes et al., 1993, Geophysics
    call get_freq_parameters(binsize,header%srate,bb,freq)

    if (mod(bin,int(binsize/overlap)) .eq. 0) then
      end_sample = start_sample + (binsize-1)
      bestbeam(start_sample:end_sample,1) = bb(1:binsize,1)
    endif    

    deallocate(bb)
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


  subroutine get_cc_traces(counts_bin,binsize,fbinsize,zpd,srate,n_instr,cc_traces)
   integer i, j, n_instr, n_pairs, binsize, fbinsize, zpd
   integer, parameter :: ndim = 1
   double precision c_max, srate, t_max
   double precision, allocatable, dimension (:)   :: tau, c_wrap
   double precision, allocatable, dimension (:,:) :: counts_bin, s1, s2
   double precision, allocatable, dimension (:,:) :: cc_traces

   allocate(s1(binsize,ndim))
   allocate(s2(binsize,ndim))
   allocate(tau(fbinsize))
   allocate(c_wrap(fbinsize))

   n_pairs = 0
   do i = 1,n_instr
    do j = i+1,n_instr
     n_pairs = n_pairs + 1
     s1(1:binsize,ndim) = counts_bin(1:binsize,i)
     s2(1:binsize,ndim) = counts_bin(1:binsize,j)

     call xcorr(s1,s2,ndim,binsize,fbinsize,srate,tau,c_wrap,t_max,c_max)
     cc_traces(1:fbinsize+2*zpd,n_pairs) = 0.0
     cc_traces(zpd+1:zpd+fbinsize,n_pairs) = c_wrap(:)
    end do 
   end do

   deallocate(s1)
   deallocate(s2)
   deallocate(tau)
   deallocate(c_wrap)

  end subroutine

  end module
