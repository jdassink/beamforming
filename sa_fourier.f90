  ! Basic seismo-acoustic Fourier functions
 
  module sa_fourier
  use sa_math
  implicit none

  include 'fftw3.f'
  
  type :: Frequency
  double precision dominant, center, bandwith, rms
  double precision f_min, f_max, f_step, average
  integer smoother
  end type
  
  contains

  subroutine compute_dft_1d(m,n,input,output)
    ! Transforms M 1-D double precision arrays with N samples each
    ! to M 1-D double complex Fourier transformed arrays
    integer n, m, i
    integer*8 plan
    double precision, allocatable, dimension(:,:) :: input
    double complex, allocatable, dimension(:,:) :: output

    do i=1,m
      call dfftw_plan_dft_r2c_1d(plan,n,input(:,i),output(:,i),FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan,input(:,i),output(:,i))
      call dfftw_destroy_plan(plan)
    end do
  end subroutine

  subroutine compute_idft_1d(m,n,input,output)
    ! Inverse transform M 1-D double precision arrays with N samples each
    ! to M 1-D double complex Fourier transformed arrays
    integer n, m, i
    integer*8 iplan
    double precision, allocatable, dimension(:,:) :: output
    double complex, allocatable, dimension(:,:) :: input

    do i=1,m
      call dfftw_plan_dft_c2r_1d(iplan,n,input(:,i),output(:,i),FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(iplan,input(:,i),output(:,i))
      call dfftw_destroy_plan(iplan)
    end do
  end subroutine

  subroutine get_freq_parameters(n_samples,srate,input,freq)
    ! From Arthur E. Barnes, Geophysics, 1993
    integer n_samples, fr
    double precision srate, df, max_ampl, ampl, ampl_, Power, WPower, f
    double precision, allocatable, dimension (:,:) :: input
    double complex, allocatable, dimension(:,:) :: ft_input
    type(Frequency) :: freq

    allocate(ft_input(n_samples,1))

    call compute_dft_1d(1,n_samples,input,ft_input)

    df = 1.0*srate/n_samples
    
    max_ampl = 0.0
    f        = (2-1)*df
    ampl     = abs(ft_input(2,1))
    ampl_    = abs(ft_input(n_samples/2+1,1))
    WPower   = (0.5 * df) * f * (ampl**2 + ampl_**2)
    Power    = (0.5 * df) *     (ampl**2 + ampl_**2)

    do fr=3,n_samples/2
      ampl = abs(ft_input(fr,1))
      f    = (fr-1)*df
      if (ampl .ge. max_ampl) then
        freq%dominant = (fr-1)*df
        max_ampl = ampl
      endif
      WPower = WPower + df * ( f * ampl**2 )
      Power  = Power  + df * (     ampl**2 )
    end do
    
    freq%center = WPower / Power

    ! Start again for the bandiwth, as it is dependent on central frequency
    f      = (2-1)*df
    ampl   = abs(ft_input(2,1))
    ampl_  = abs(ft_input(n_samples/2+1,1))

    WPower = (0.5 * df) * (ampl**2 + ampl_**2) * (f - freq%center)**2

    do fr=3,n_samples/2
       ampl = abs(ft_input(fr,1))
       f    = (fr-1)*df
       WPower = WPower + df * ( ampl**2 * (f - freq%center)**2 )
    end do

    freq%bandwith = sqrt(WPower / Power)
    freq%rms = sqrt(freq%center**2 + freq%bandwith**2)
    
    deallocate(ft_input)
  end subroutine


  subroutine xcorr(s1,s2,ndim,size,fft_size,srate,tau,c_wrap,t_max,c_max)
   integer i, j, size, fft_size, fft_real, ndim
   double precision c_max, t_max, srate
   
   double precision, allocatable, dimension(:)   :: tau, c_wrap
   double precision, allocatable, dimension(:,:) :: s1, s2, c_t, b_t, a_t 
   double complex  , allocatable, dimension(:,:) :: a, b, c, s1_fft, s2_fft

   allocate(s1_fft(size,ndim))
   allocate(s2_fft(size,ndim))
   allocate(a(fft_size,ndim))
   allocate(b(fft_size,ndim))
   allocate(c(fft_size,ndim))
   allocate(a_t(fft_size,ndim))
   allocate(b_t(fft_size,ndim))
   allocate(c_t(fft_size,ndim))

   call compute_dft_1d(ndim,size,s1,s1_fft)
   call compute_dft_1d(ndim,size,s2,s2_fft)

   ! The transformed array only has entries up to binsize_2+1
   ! The other half, reserved for negative frequencies, is empty
   ! as a real input signal is Fourier transformed

   fft_real = size/2
   do i=1,fft_real+1
    c(i,ndim) = s1_fft(i,ndim) * conjg(s2_fft(i,ndim))/fft_real
    a(i,ndim) = s1_fft(i,ndim) * conjg(s1_fft(i,ndim))/fft_real
    b(i,ndim) = s2_fft(i,ndim) * conjg(s2_fft(i,ndim))/fft_real
   end do

   ! Optional zeropadding in frequency domain for upsampling
   do i=fft_real+2,fft_size
    c(i,ndim) = dcmplx(0.0,0.0)
    a(i,ndim) = dcmplx(0.0,0.0)
    b(i,ndim) = dcmplx(0.0,0.0)
   end do
  
   c(1,ndim) = dcmplx(real(c(1,ndim)),real(c(fft_real+1,ndim)))
   a(1,ndim) = dcmplx(real(a(1,ndim)),real(a(fft_real+1,ndim)))
   b(1,ndim) = dcmplx(real(b(1,ndim)),real(b(fft_real+1,ndim)))

   call compute_idft_1d(ndim,fft_size,a,a_t)
   call compute_idft_1d(ndim,fft_size,b,b_t)
   call compute_idft_1d(ndim,fft_size,c,c_t)

   ! Re-arrange the tau and cross correlation coefficient matrix
   ! so that the center is where zero is.
   ! The cross-correlation function is normalized so that the
   ! auto-correlation at zero-lag time is unity.
   fft_real = fft_size/2
   do i=1,fft_real
    j         = i+fft_real
    c_wrap(j) = c_t(i,ndim)/sqrt((a_t(1,ndim)*b_t(1,ndim)))
    tau(j)    = (i-1)*1.0/srate
   end do
   do i=fft_real+1,fft_size
    j         = i-fft_real
    c_wrap(j) = c_t(i,ndim)/sqrt((a_t(1,ndim)*b_t(1,ndim)))
    tau(j)    = (i-fft_size-1)*1.0/srate
   end do

   c_max  = 0.
   t_max  = 0.
   
   do j=1,fft_size
    if (c_wrap(j) > c_max) then
      c_max = c_wrap(j)
      t_max = tau(j)
    endif
   end do

   deallocate(a)
   deallocate(b)
   deallocate(c)
   deallocate(a_t)
   deallocate(b_t)
   deallocate(c_t)
   deallocate(s1_fft)
   deallocate(s2_fft)
  end subroutine xcorr 

  subroutine window_data(n_instr,binsize,counts_bin)
    integer   i, n_instr
    integer   binsize
    double precision hann
    double precision, allocatable, dimension(:,:) :: counts_bin

    do i=1,binsize
      hann = 0.5*(1.0-cos((2*pi*(i-1))/(binsize-1)))
      counts_bin(i,1:n_instr) = counts_bin(i,1:n_instr)*hann*2.0
    end do
  end subroutine

  subroutine zeropadding(n_instr,start_sample,tbinsize,fbinsize,counts,counts_bin)
    integer   n_instr
    integer   tbinsize, fbinsize, start_sample, end_sample
    double precision, allocatable, dimension(:,:) :: counts, counts_bin

    end_sample = start_sample + (tbinsize-1)
    counts_bin(1:fbinsize,1:n_instr) = 0.0
    counts_bin(1:tbinsize,1:n_instr) = counts(start_sample:end_sample,1:n_instr)
  end subroutine

  end module
