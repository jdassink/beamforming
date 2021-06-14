  ! Frequency-Domain Fisher Detector with command line parameters
  ! based on "Fast Frequency-Wavenumber Analysis and Fisher Signal Detection
  ! in Real-Time Infrasonic Array Data Processing", by Smart & Flinn, GJI 1971

  ! Author     : Jelle D. Assink (assink@knmi.nl)
  ! Date       : November 2010

  module fk_fisher
  use string_utility
  use sa_array
  use sa_fourier
  implicit none

  contains 

  subroutine compute_freqfisher(freq,fr,df,n_instr,s_grid,r,ffisher_matrix,results)
    integer n_instr, k, fr
    double precision fisher_ratio, df, px, py, e_w, e_wp_abs
    double precision, allocatable, dimension (:) :: results
    double precision, allocatable, dimension (:,:) :: r
    double complex, allocatable, dimension(:,:) :: ffisher_matrix

    type(Slowness) :: s_grid
    type(Frequency) :: freq

    !$omp  parallel &
    !$omp& default(private) &
    !$omp& shared(r,results,s_grid,n_instr,ffisher_matrix, freq, fr, df)
    !$omp do
    do k=1,s_grid%n_beams
      px = s_grid%px(k)
      py = s_grid%py(k)
      call compute_fk_terms(n_instr,px,py,r,freq,fr,df,ffisher_matrix,e_wp_abs,e_w)

      fisher_ratio = ( e_wp_abs / (e_w-e_wp_abs) ) * (n_instr-1)
      results(k) = fisher_ratio

    end do
    !$omp end do 
    !$omp end parallel

  end subroutine

  subroutine compute_fk_terms(n_instr,px,py,r,freq,fri,df,ffm,e_wp_abs,e_w)
    ! Calculates the sums of E(omega) and E(omega,p)
    ! Optionally, estimates of E(omega) and E(omega,p) are averaged over frequency
    integer n_instr, fri, fr, frs, fre, instr, dfr, frcnt
    double precision px, py, omega, df
    double precision e_w, dt, e_wp_abs
    double complex e_wp
    double complex, allocatable, dimension(:,:) :: ffm
    double precision, allocatable, dimension (:,:) :: r
    
    type(Frequency) :: freq

    dfr = nint(freq%f_step/df)     
    frs = fri
    fre = fri + dfr*freq%smoother - 1

    e_w   = 0.                        ! Total power estimate
    e_wp_abs = 0.
    freq%average = 0.
    frcnt = 0

    do fr = frs, fre, dfr
      e_wp  = dcmplx(0.,0.)             ! Signal estimate

      frcnt = frcnt + 1
      omega = 2*pi*(fr-1)*df
      freq%average = freq%average + (fr-1)*df

      do instr=1,n_instr
        dt = px*r(instr,1) + py*r(instr,2)
        e_wp = e_wp + ffm(fr,instr)*exp(-im*omega*dt)
        e_w  = e_w  + abs(ffm(fr,instr))**2
      end do
      e_wp_abs = e_wp_abs + abs(e_wp / n_instr)**2
    end do

    freq%average = freq%average / frcnt
    e_wp_abs = e_wp_abs / frcnt
    e_w = ( e_w / n_instr ) / frcnt
  end subroutine

  subroutine get_cmd_parameters(n_instr,r,tbinsize,overlap,freq,s_grid,file_format,timeseries)
    integer i, n_instr, eof, overlap, iargc
    double precision tbinsize
    character(40), allocatable, dimension (:) :: timeseries
    character(40) coordinates, dummy, file_format
    double precision, allocatable, dimension(:,:) :: r

    type(Slowness) :: s_grid
    type(Frequency) :: freq

    call getarg(1, dummy)
    read(dummy,*,IOSTAT=eof) coordinates
    if (eof < 0) then
    call printUsage()
    endif
    call get_coordinates(coordinates,n_instr,r)

    if (iargc() .ne. 15 + n_instr) then
     deallocate(r)
     deallocate(timeseries)
     call printUsage()
    endif             
    call getarg(2, dummy)
    read(dummy,*) tbinsize
    call getarg(3, dummy)
    read(dummy,*) overlap
    call getarg(4, dummy)
    read(dummy,*) freq%f_min
    call getarg(5, dummy)
    read(dummy,*) freq%f_max
    call getarg(6, dummy)
    read(dummy,*) freq%f_step
    call getarg(7, dummy)
    read(dummy,*) freq%smoother
    call getarg(8, dummy)
    read(dummy,*) s_grid%bearing_min
    call getarg(9, dummy)
    read(dummy,*) s_grid%bearing_max
    call getarg(10, dummy)
    read(dummy,*) s_grid%dth
    call getarg(11, dummy)
    read(dummy,*) s_grid%trcvel_min
    call getarg(12, dummy)
    read(dummy,*) s_grid%trcvel_max
    call getarg(13, dummy)
    read(dummy,*) s_grid%dc
    call getarg(14, dummy)
    dummy = StrLowCase(dummy)
    read(dummy,*) s_grid%tele_local
    call getarg(15, dummy)
    dummy = StrLowCase(dummy)
    read(dummy,*) file_format
    i = 1
    do while ( i < n_instr+1 )
    call getarg(i+15, timeseries(i) )
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
    write(6,*) ' Frequency-domain Fisher detector and beamformer for SAC/ASCII data'
    write(6,*) '________________________________________________________________________'
    write(6,*) ''
    write(6,*) 'Usage:'
    write(6,*) '------'
    write(6,*) 'freqfisher <coordinate file> <binsize> <overlap>'
    write(6,*) '                  <f_min> <f_max> <f_step> <f_averaging>'
    write(6,*) '                  <th_min> <th_max> <dth> <c_min> <c_max> <dc>'
    write(6,*) '                  <tele|local> <ascii|sac> <files>'
    write(6,*) ' '
    write(6,*) ' - Order of coordinates must match the order of input files'
    write(6,*) ' - Sample rate must be provided in the input file'
    write(6,*) ' '
    write(6,*) ' '
    write(6,*) 'Example: '
    write(6,*) '-------- '
    write(6,*) 'freqfisher stationtable 20.0 50'
    write(6,*) '           0.1 10.0 0.01 10'
    write(6,*) '           0.0 360.0 2.0 300.0 450.0 5.0'
    write(6,*) '           tele sac DBN*.sac'
    write(6,*) ' '
    write(6,*) ' - FK Fisher processing (max. 24 hours of data)'
    write(6,*) ' - Timebins of 20.0 s, with 50% overlap = 10.0 s'
    write(6,*) ' - Processing between 0.1 and 10.0 Hz, in steps of 0.01 Hz'
    write(6,*) ' - It will average over 10 frequency bands (output every 0.1 Hz).'
    write(6,*) ' - The beamforming resolution is 5 m/s and 2 deg.'
    write(6,*) ' - Trace velocity range: 300-450 m/s; bearing range: 0-360 deg.'
    write(6,*) '   -  [ tele: only search in specified trace velocity range     ]'
    write(6,*) '   -  [ local: include a search for (near) vertical propagation ]'
    write(6,*) ' - Binary SAC input is expected, files DBN*.sac'
    write(6,*) ''
    write(6,*) 'The output file <freqfisher.dat> format is:'
    write(6,*) ' '
    write(6,*) '   1. Timebin'
    write(6,*) '   2. Frequency      [Hz]'
    write(6,*) '   3. Fisher-ratio'
    write(6,*) '   4. Back azimuth   [deg]'
    write(6,*) '   5. Trace velocity [m/s]'
    write(6,*) '   6. # of instruments used in analysis'
    write(6,*) '   7. Power spectral density [amp^2/Hz]'
    write(6,*) ' '
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    call exit(0)
  end subroutine

  end module

!--------------------------------------------------------------------------------------

  program freqfisher
    use fk_fisher
    implicit none

    integer   fr, n_instr, fr_min, fr_max, dfr, alloc_instr, alloc_samples, k
    integer   binsize, fbinsize, bin, n_bins, overlap, start_sample, end_sample
    double precision tbinsize, fisher, bearing, trcvel
    double precision px, py, time, df, e_w, e_wp_abs

    double precision, allocatable, dimension(:) :: results
    double precision, allocatable, dimension(:,:) :: r, counts, counts_bin
    double complex, allocatable, dimension(:,:) :: ffisher_matrix
    character(40), allocatable, dimension (:) :: timeseries
    character(40) file_format

    type(Slowness)  :: s_grid
    type(Frequency) :: freq
    type(WFHeader) :: header

    alloc_instr = 50
    alloc_samples = 2*24*3600*250

    allocate(r(alloc_instr,5))
    allocate(timeseries(alloc_instr))

    call get_cmd_parameters(n_instr,r,tbinsize,overlap,freq,s_grid,file_format,timeseries)
    call span_slowness_grid(s_grid)

    allocate(counts(alloc_samples,n_instr))
    allocate(results(s_grid%n_beams))
    
    write(6,*) ' '
    write(6,*) ' '
    write(6,*) 'Reading files ...'
    write(6,*) '----------------------------------'
    call get_timeseries(n_instr,file_format,timeseries,alloc_samples,header,counts)

    ! Zeropadding parameter (should be a power of 2)
    binsize = nint(tbinsize*header%srate)
    fbinsize = 2**nint(log10(1.0*header%srate/freq%f_step)/log10(2.0)+1)
    if (fbinsize .lt. binsize) then
     fbinsize = 2**nint(log10(1.0*binsize)/log10(2.0))
     if (fbinsize .lt. binsize) then
      fbinsize = 2**(nint(log10(1.0*binsize)/log10(2.0))+1)
     endif
    endif

    df       = 1.0*header%srate/fbinsize
    fr_min   = nint(freq%f_min/df)+1 
    fr_max   = nint(freq%f_max/df)+1
    dfr      = nint(freq%f_step/df) 

    if (freq%f_min .lt. 1./tbinsize) then
      write(6,*)  'ERROR: Window size is too small to resolve minimum frequency.'
      write(6,*)  ''
      write(6,*)  ' - This can be changed by increasing the binsize'
      write(6,96) ' - f_min must be greater than ', 1./tbinsize, ' Hz (now: ', freq%f_min, ' Hz).'
      write(6,*)  ' - Exiting ...'
      96 format (a31, e10.3, a10, e10.3, a5)
      call exit(0)
    endif

    allocate(counts_bin(fbinsize,n_instr))
    allocate(ffisher_matrix(fbinsize,n_instr))

    write(6,*) ''
    write(6,*) ''
    write(6,*) 'F-K Fisher beamforming ...'
    write(6,*) '--------------------------'
    write(6,*) ''

    open (unit=30, file='freqfisher.dat')

    write(6,'(a32,f8.2,a4)') ' Time-domain sampling : ', header%srate, ' Hz'
    write(6,'(a32,f8.2,a4)') ' Time bin length      : ', tbinsize, ' s'
    write(6,'(a32,f8.2,a4)') ' Time bin overlap     : ', overlap/100.*tbinsize, ' s'
    write(6,*) ''
    write(6,'(a32,f6.3,a4)') ' Frequency band sampling  : ', df, ' Hz'
    write(6,'(a32,f6.3,a4)') ' Frequency band stepping  : ', freq%f_step, ' Hz'
    write(6,'(a32,f6.3,a4)') ' Frequency band averaging  : ', freq%smoother * freq%f_step, ' Hz'
    write(6,'(a32,f5.1,a4)') ' Trace velocity resolution : ', s_grid%dc, ' m/s'
    write(6,'(a32,f5.1,a4)') ' Bearing resolution : ', s_grid%dth, ' deg'
    write(6,*) ''

    write(6,'(a32,f8.3,a3,f8.3,a5)') ' Frequency band : [ ', freq%f_min, ' - ', freq%f_max, ' ] Hz'
    write(6,'(a32,f8.2,a3,f8.2,a6)') ' Trace velocity domain : [ ', s_grid%trcvel_min, ' - ', s_grid%trcvel_max, ' ] m/s'
    write(6,'(a32,f8.2,a3,f8.2,a6)') ' Bearing domain : [ ', s_grid%bearing_min, ' - ', s_grid%bearing_max, ' ] deg'

    overlap = binsize - binsize*overlap/100.
    n_bins = header%n_samples/overlap

    do bin = 1 , n_bins 
      ffisher_matrix(1:fbinsize,1:n_instr) = dcmplx(0.0)

      start_sample = (bin-1)*overlap+1
      end_sample   = start_sample + (binsize-1)
      time         = header%tzero + ((start_sample + end_sample-1) / 2.0 ) / header%srate

      call zeropadding(n_instr,start_sample,binsize,fbinsize,counts,counts_bin)
      call window_data(n_instr,binsize,counts_bin)
      call compute_dft_1d(n_instr,fbinsize,counts_bin,ffisher_matrix)

      ! Scale FFTW result appropriately for positive frequencies and binsize
      ffisher_matrix = 2*ffisher_matrix / binsize

      fr = fr_min
      do while ( ( fr + dfr*freq%smoother - 1 ) < fr_max )
        call compute_freqfisher(freq,fr,df,n_instr,s_grid,r,ffisher_matrix,results)
        call beamgrid_maximum(results,k)
        px = s_grid%px(k)
        py = s_grid%py(k)
        call compute_fk_terms(n_instr,px,py,r,freq,fr,df,ffisher_matrix,e_wp_abs,e_w)

        write(30,99) time, freq%average, results(k), s_grid%bearing(s_grid%bi(k)), s_grid%trcvel(s_grid%ci(k)), n_instr, e_wp_abs
        99 format(f10.2, 1x, f10.5, 1x, f15.2, 1x, f9.2, 1x, f9.2, 1x, i2, 1x, e15.8)

        fr = fr + dfr*freq%smoother
      end do
    write(30,*) ''

    end do
    close (unit=30)

    deallocate(r)
    deallocate(s_grid%bearing)
    deallocate(s_grid%trcvel)
    deallocate(s_grid%bi)
    deallocate(s_grid%ci)
    deallocate(s_grid%px)
    deallocate(s_grid%py)
    deallocate(counts)
    deallocate(counts_bin)
    deallocate(ffisher_matrix)
    deallocate(results)
    deallocate(timeseries)

  end program freqfisher

