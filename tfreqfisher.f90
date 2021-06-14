  ! !!!! BETA VERSION !!!!
  ! To-do: make input from timefisher.dat more sophisticated

  ! Frequency-Domain Fisher Detector with command line parameters
  ! Slowness data is provided by <timefisher.dat>

  ! based on "Fast Frequency-Wavenumber Analysis and Fisher Signal Detection
  ! in Real-Time Infrasonic Array Data Processing", by Smart & Flinn, GJI 1971
  ! Author     : Jelle D. Assink (assink@knmi.nl)
  ! Date       : November 2010

  module tfk_fisher
  use string_utility
  use sa_array
  use sa_fourier
  implicit none

  contains 

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

  subroutine get_cmd_parameters(n_instr,r,tbinsize,overlap,freq,file_format,timeseries)
    integer i, n_instr, eof, overlap, iargc
    double precision tbinsize
    character(40), allocatable, dimension (:) :: timeseries
    character(40) coordinates, dummy, file_format
    double precision, allocatable, dimension(:,:) :: r

    type(Frequency) :: freq

    call getarg(1, dummy)
    read(dummy,*,IOSTAT=eof) coordinates
    if (eof < 0) then
    call printUsage()
    endif
    call get_coordinates(coordinates,n_instr,r)

    if (iargc() .ne. 8 + n_instr) then
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
    dummy = StrLowCase(dummy)
    read(dummy,*) file_format
    i = 1
    do while ( i < n_instr+1 )
    call getarg(i+8, timeseries(i) )
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
    write(6,*) '___________________________________________________________________________'
    write(6,*) ' '
    write(6,*) ' Driven Frequency-domain Fisher detector and beamformer for SAC/ASCII data'
    write(6,*) '___________________________________________________________________________'
    write(6,*) ''
    write(6,*) 'Usage:'
    write(6,*) '------'
    write(6,*) 'tfreqfisher <coordinate file> <binsize> <overlap>'
    write(6,*) '                  <f_min> <f_max> <f_step> <f_averaging>'
    write(6,*) '                  <ascii|sac> <files>'
    write(6,*) ' '
    write(6,*) ' - Order of coordinates must match the order of input files'
    write(6,*) ' - Sample rate must be provided in the input file'
    write(6,*) ' - Slowness is provided by <timefisher.dat>'
    write(6,*) ' '
    write(6,*) 'Example: '
    write(6,*) '-------- '
    write(6,*) 'tfreqfisher stationtable 20.0 50'
    write(6,*) '           0.1 10.0 0.01 10'
    write(6,*) '           sac DBN*.sac'
    write(6,*) ' '
    write(6,*) ' - FK Fisher processing (max. 24 hours of data)'
    write(6,*) ' - Timebins of 20.0 s, with 50% overlap = 10.0 s'
    write(6,*) ' - Processing between 0.1 and 10.0 Hz, in steps of 0.01 Hz'
    write(6,*) ' - It will average over 10 frequency bands (output every 0.1 Hz).'
    write(6,*) ' - Binary SAC input is expected, files DBN*.sac'
    write(6,*) ''
    write(6,*) 'The output file <tfreqfisher.dat> format is:'
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

  program tfreqfisher
    use tfk_fisher
    implicit none

    integer   fr, n_instr, fr_min, fr_max, dfr, alloc_instr, alloc_samples, k, ntf, tmp(1)
    integer   binsize, fbinsize, bin, n_bins, overlap, start_sample, end_sample
    double precision tbinsize, fisher, bearing, trcvel, px, py, time, df, e_w, e_wp_abs

    double precision, allocatable, dimension(:,:) :: r, counts, counts_bin, tff
    double complex, allocatable, dimension(:,:) :: ffisher_matrix
    character(40), allocatable, dimension (:) :: timeseries
    character(40) file_format, dummy

    type(Frequency) :: freq
    type(WFHeader) :: header

    alloc_instr = 50
    alloc_samples = 2*24*3600*250

    allocate(r(alloc_instr,5))
    allocate(timeseries(alloc_instr))

    call get_cmd_parameters(n_instr,r,tbinsize,overlap,freq,file_format,timeseries)

    allocate(counts(alloc_samples,n_instr))
    
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
    write(6,*) 'F-K Fisher timefisher-driven beamforming (beta)...'
    write(6,*) '--------------------------------------------------'
    write(6,*) ''

    open (unit=30, file='tfreqfisher.dat')
    open (unit=60, file='timefisher.dat')

    allocate(tff(8000000,3))
    ntf = 0.
    do while (.true.)
      ntf = ntf + 1
      read (60,*, end=991) tff(ntf,1), dummy, dummy, dummy, tff(ntf,2), tff(ntf,3)
    end do
    991 continue
    ntf = ntf - 1
    close(unit=60)

    write(6,'(a32,f8.2,a4)') ' Time-domain sampling : ', header%srate, ' Hz'
    write(6,'(a32,f8.2,a4)') ' Time bin length      : ', tbinsize, ' s'
    write(6,'(a32,f8.2,a4)') ' Time bin overlap     : ', overlap/100.*tbinsize, ' s'
    write(6,*) ''
    write(6,'(a32,f6.3,a4)') ' Frequency band sampling  : ', df, ' Hz'
    write(6,'(a32,f6.3,a4)') ' Frequency band stepping  : ', freq%f_step, ' Hz'
    write(6,'(a32,f6.3,a4)') ' Frequency band averaging  : ', freq%smoother * freq%f_step, ' Hz'
    write(6,*) ''

    write(6,'(a32,f8.3,a3,f8.3,a5)') ' Frequency band : [ ', freq%f_min, ' - ', freq%f_max, ' ] Hz'

    overlap = binsize - binsize*overlap/100.
    n_bins = header%n_samples/overlap

    do bin = 1 , n_bins 
      ffisher_matrix(1:fbinsize,1:n_instr) = dcmplx(0.0)

      start_sample = (bin-1)*overlap+1
      end_sample   = start_sample + (binsize-1)
      time         = header%tzero + ((start_sample + end_sample-1) / 2.0 ) / header%srate

      ! Load slowness from timefisher.dat file
      tmp = minloc(abs(time-tff(:,1)))
      ntf = tmp(1)
      px = tff(ntf,2)
      py = tff(ntf,3)
      call convert_slowness(px,py,bearing,trcvel)

      call zeropadding(n_instr,start_sample,binsize,fbinsize,counts,counts_bin)
      call window_data(n_instr,binsize,counts_bin)
      call compute_dft_1d(n_instr,fbinsize,counts_bin,ffisher_matrix)

      ! Scale FFTW result appropriately for positive frequencies and binsize
      ffisher_matrix = 2*ffisher_matrix / binsize

      fr = fr_min
      do while ( ( fr + dfr*freq%smoother - 1 ) < fr_max )
        call compute_fk_terms(n_instr,px,py,r,freq,fr,df,ffisher_matrix,e_wp_abs,e_w)
        fisher = ( e_wp_abs / (e_w-e_wp_abs) ) * (n_instr-1)

        write(30,99) time, freq%average, fisher, bearing, trcvel, n_instr, e_wp_abs
        99 format(f10.2, 1x, f10.5, 1x, f15.2, 1x, f9.2, 1x, f9.2, 1x, i2, 1x, e15.8)

        fr = fr + dfr*freq%smoother
      end do
    write(30,*) ''

    end do
    close (unit=30)

    deallocate(r)
    deallocate(counts)
    deallocate(counts_bin)
    deallocate(ffisher_matrix)
    deallocate(timeseries)

  end program tfreqfisher

