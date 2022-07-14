  ! Time-Domain Fisher Detector with command line parameters
  ! based on "Multiple Signal Correlators", by Melton & Bailey, Geophysics 1957

  ! Author     : Jelle D. Assink (assink@knmi.nl)
  ! Date       : November 2010
 
  module td_fisher
  use string_utility
  use sa_array
  implicit none

  contains 

  subroutine get_cmd_parameters(n_instr,tbinsize,overlap,s_grid,r,file_format,timeseries)
    integer i, n_instr, eof, overlap, iargc
    double precision tbinsize
    character(40), allocatable, dimension (:) :: timeseries
    character(40) coordinates, file_format, dummy
    double precision, allocatable, dimension(:,:) :: r

    type(Slowness) :: s_grid

    call getarg(1, dummy)
    read(dummy,*,IOSTAT=eof) coordinates
    if (eof < 0) then
      call printUsage()
    endif
    call get_coordinates(coordinates,n_instr,r)

    if (iargc() .ne. 12 + n_instr) then
      deallocate(r)
      deallocate(timeseries)
      call printUsage()
    endif             
    call getarg(2, dummy)
    read(dummy,*) tbinsize
    call getarg(3, dummy)
    read(dummy,*) overlap
    call getarg(4, dummy)
    read(dummy,*) s_grid%bearing_min
    call getarg(5, dummy)
    read(dummy,*) s_grid%bearing_max
    call getarg(6, dummy)
    read(dummy,*) s_grid%dth
    call getarg(7, dummy)
    read(dummy,*) s_grid%trcvel_min
    call getarg(8, dummy)
    read(dummy,*) s_grid%trcvel_max
    call getarg(9, dummy)
    read(dummy,*) s_grid%dc
    call getarg(10, dummy)
    dummy = StrLowCase(dummy)
    read(dummy,*) s_grid%tele_local
    call getarg(11, dummy)
    dummy = StrLowCase(dummy)
    read(dummy,*) s_grid%max_all
    call getarg(12, dummy)
    dummy = StrLowCase(dummy)
    read(dummy,*) file_format
    i = 1
    do while ( i < n_instr+1 )
      call getarg(i+12, timeseries(i) )
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
    write(6,*) 'Time-domain Fisher detector and beamformer for ASCII/SAC formatted data'
    write(6,*) '________________________________________________________________________'
    write(6,*) ''
    write(6,*) 'Usage:'
    write(6,*) '------'
    write(6,*) 'timefisher <coordinate file> <binsize> <overlap>'
    write(6,*) '                  <th_min> <th_max> <dth> <c_min> <c_max> <dc>'
    write(6,*) '                  <tele|local> <max|all> <ascii|sac> <input ASCII/SAC>'
    write(6,*) ' '
    write(6,*) ' '
    write(6,*) ' - Order of coordinates must match the order of input files'
    write(6,*) ' - Sample rate must be provided in the input file'
    write(6,*) ' - max  | all   : toggle between printing full grid or just Fisher maximum'
    write(6,*) ' - tele | local : toggle between far or near-field grid'
    write(6,*) ' '
    write(6,*) 'Example: '
    write(6,*) '-------- '
    write(6,*) 'timefisher stations 10.0 50 0.0 360.0 2.0 300.0 450.0 5.0 tele max sac DBN*.sac'
    write(6,*) ' '
    write(6,*) ' - Time-domain Fisher processing (max. 24 hours of data)'
    write(6,*) ' - Timebins of 10.0 s, with 50% overlap = 5.0 s'
    write(6,*) ' - The beamforming resolution is 5 m/s and 2 deg.'
    write(6,*) ' - Trace velocity range: 300-450 m/s; bearing range: 0-360 deg.'
    write(6,*) ' - Only one maximum Fisher per timebin is printed.'
    write(6,*) ' - Binary SAC input is expected, files DBN*.sac'
    write(6,*) ''
    write(6,*) ' The output file <timefisher.dat> format is:'
    write(6,*) ' '
    write(6,*) '   1. Timebin'
    write(6,*) '   2. Fisher-ratio'
    write(6,*) '   3. Back azimuth   [deg]'
    write(6,*) '   4. Trace velocity [m/s]'
    write(6,*) '   5. Slowness X     [s/m]'
    write(6,*) '   6. Slowness Y     [s/m]'
    write(6,*) '   7. RMS amplitude'
    write(6,*) '   8. Peak-to-peak amplitude'    
    write(6,*) '   9. # of instruments used in analysis'
    write(6,*) '  10. Center frequency     [Hz]'
    write(6,*) '  11. Frequency bandwidth  [Hz]'
    write(6,*) '  12. RMS frequency        [Hz]'
    write(6,*) ' '
    write(6,*) ' Bestbeam output is written out as SAC file <bestbeam.sac>'
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    call exit(0)
  end subroutine

  end module

!--------------------------------------------------------------------------------------

  program timefisher
    use td_fisher
    implicit none

    integer   k, n_instr, alloc_samples, alloc_instr, ci, bi
    integer   bin, n_bins, binsize, overlap, start_sample, end_sample
    double precision fisher, px, py, bearing, trcvel, prms, p2p, time, tbinsize

    type(Slowness) :: s_grid
    type(Frequency) :: freq
    type(WFHeader) :: header

    double precision, allocatable, dimension(:) :: results
    double precision, allocatable, dimension(:,:) :: r, counts, bestbeam
    character(40), allocatable, dimension (:) :: timeseries
    character(40) file_format, fid_beam, statistic

    statistic = 'fisher'
    alloc_instr = 50
    alloc_samples = 2*24*3600*250

    allocate(timeseries(alloc_instr))
    allocate(r(alloc_instr,5))

    call get_cmd_parameters(n_instr,tbinsize,overlap,s_grid,r,file_format,timeseries)
    call span_slowness_grid(s_grid)

    allocate(counts(alloc_samples,n_instr))
    allocate(bestbeam(alloc_samples,1))
    allocate(results(s_grid%n_beams))

    write(6,*) ' '
    write(6,*) ' '
    write(6,*) 'Reading files ...'
    write(6,*) '-----------------'
    call get_timeseries(n_instr,file_format,timeseries,alloc_samples,header,counts)

    write(6,*) ''
    write(6,*) ''
    write(6,*) 'Time-domain Fisher beamforming ...'
    write(6,*) '----------------------------------'
    write(6,*) ''

    write(6,'(a32,f8.2,a4)') ' Time-domain sampling : ', header%srate, ' Hz'
    write(6,'(a32,f8.2,a4)') ' Time bin length      : ', tbinsize, ' s'
    write(6,'(a32,f8.2,a4)') ' Time bin overlap     : ', overlap/100.*tbinsize, ' s'
    write(6,*) ''
    write(6,'(a32,f5.1,a4)') ' Trace velocity resolution : ', s_grid%dc, ' m/s'
    write(6,'(a32,f5.1,a4)') ' Bearing resolution : ', s_grid%dth, ' deg'
    write(6,*) ''

    write(6,'(a32,f8.2,a3,f8.2,a6)') ' Trace velocity domain : [ ', s_grid%trcvel_min, ' - ', s_grid%trcvel_max, ' ] m/s'
    write(6,'(a32,f8.2,a3,f8.2,a6)') ' Bearing domain : [ ', s_grid%bearing_min, ' - ', s_grid%bearing_max, ' ] deg'
    write(6,*) ''
    write(6,'(a32,a3)') ' Output type : ', s_grid%max_all

    open (unit=30, file='timefisher.dat')
    fid_beam = 'bestbeam.sac'
    
    binsize = nint(tbinsize*header%srate)
    overlap = binsize - binsize*overlap/100.
    n_bins = (binsize/overlap)*(header%n_samples/binsize)

    do bin = 2, n_bins - 1
      start_sample  = (bin-1)*overlap+1
      end_sample    = start_sample + (binsize-1)
      time          = header%tzero + ((start_sample + end_sample-1) / 2.0 ) / header%srate

      call beamgrid(statistic,start_sample,header%srate,binsize,s_grid,n_instr,r,counts,results)

      if (s_grid%max_all .eq. 'max') then
       call beamgrid_maximum(results,k)
       px = s_grid%px(k)
       py = s_grid%py(k)

       call get_bestbeam(bin,start_sample,header,counts,binsize,n_instr,overlap,r,px,py,bestbeam,prms,p2p,freq)

       if ((s_grid%bi(k) /= 1) .and. &
       &   (s_grid%bi(k) /= s_grid%bi(size(s_grid%bi))) .and. &
       &   (s_grid%ci(k) /= 1) .and. &
       &   (s_grid%ci(k) /= s_grid%ci(size(s_grid%ci)))) then
        write (30,98) time, results(k), s_grid%bearing(s_grid%bi(k)), &
  &                  s_grid%trcvel(s_grid%ci(k)), px, py, prms, p2p, &
  &                  n_instr, freq%center, freq%bandwith, freq%rms
        98 format (f10.2, 1x, f10.2, 1x, f8.2, 1x, f8.2, 1x, e15.8, 1x, e15.8, 1x, &
    &               e15.8, 1x, e15.8, 1x, i2, 1x, f8.3, 1x, f8.3, 1x, f8.3)
       endif

      elseif (s_grid%max_all .eq. 'all') then
       do k=1,s_grid%n_beams
        write (30,99) time, results(k), s_grid%bearing(s_grid%bi(k)), s_grid%trcvel(s_grid%ci(k)), n_instr
        99 format (f10.2, 1x, f10.2, 1x, f8.2, 1x, f8.2, 1x, i2)
       end do
       write (30,*) ''
      endif
    
    end do

    if (s_grid%max_all .eq. 'max') then
     call write_sac(fid_beam,header,bestbeam)
    endif

    close(30)

    deallocate(s_grid%bearing)
    deallocate(s_grid%trcvel)
    deallocate(s_grid%bi)
    deallocate(s_grid%ci)
    deallocate(s_grid%px)
    deallocate(s_grid%py)
    
    deallocate(r)
    deallocate(counts)
    deallocate(bestbeam)
    deallocate(results)
    deallocate(timeseries)

  end program timefisher
