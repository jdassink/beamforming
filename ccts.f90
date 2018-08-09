  ! Cross-correlation Trace Stacking (CCTS) analysis
  ! based on "The European Arctic.. A Lab for Seismoacoustic Studies"
  ! by Gibbons et al., SRL. 2015

  ! Author     : Jelle D. Assink
  ! Date       : November 2010
 
  module mod_ccts
  use sa_array

  contains 

  subroutine get_cmd_parameters(n_instr,tbinsize,overlap,oversampling,s_grid,r,tele_local,max_all,file_format,timeseries)
    implicit none
    integer i, n_instr, eof, overlap, iargc, oversampling
    double precision tbinsize
    character(40), allocatable, dimension (:) :: timeseries
    character(40) coordinates, file_format, tele_local, max_all, dummy
    double precision, allocatable, dimension(:,:) :: r

    type(Slowness) :: s_grid

    call getarg(1, dummy)
    read(dummy,*,IOSTAT=eof) coordinates
    if (eof < 0) then
      call printUsage()
    endif
    call get_coordinates(coordinates,n_instr,r)

    if (iargc() .ne. 13 + n_instr) then
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
    read(dummy,*) s_grid%bearing_min
    call getarg(6, dummy)
    read(dummy,*) s_grid%bearing_max
    call getarg(7, dummy)
    read(dummy,*) s_grid%dth
    call getarg(8, dummy)
    read(dummy,*) s_grid%trcvel_min
    call getarg(9, dummy)
    read(dummy,*) s_grid%trcvel_max
    call getarg(10, dummy)
    read(dummy,*) s_grid%dc
    call getarg(11, dummy)
    read(dummy,*) tele_local
    call getarg(12, dummy)
    read(dummy,*) max_all
    call getarg(13, dummy)
    read(dummy,*) file_format
    i = 1
    do while ( i < n_instr+1 )
      call getarg(i+13, timeseries(i) )
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
    write(6,*) 'Cross-correlation Trace Stacking (CCTS) for ASCII/SAC formatted data'
    write(6,*) '________________________________________________________________________'
    write(6,*) ''
    write(6,*) 'Usage:'
    write(6,*) '------'
    write(6,*) 'ccts <coordinate file> <binsize> <overlap> <oversampling>'
    write(6,*) '                  <th_min> <th_max> <dth> <c_min> <c_max> <dc>'
    write(6,*) '                  <tele|local> <max|all> <ascii|sac> <files>'
    write(6,*) ' '
    write(6,*) ' - Order of coordinates must match the order of input files'
    write(6,*) ' - Sample rate must be provided in the input file'
    write(6,*) ' - max  | all   : toggle between printing full grid or just correlation maximum'
    write(6,*) ' - tele | local : toggle between far or near-field grid'
    write(6,*) ' '
    write(6,*) 'Example: '
    write(6,*) '-------- '
    write(6,*) 'ccts stationtable 10.0 50 2'
    write(6,*) '           0.0 360.0 2.0 300.0 450.0 5.0'
    write(6,*) '           tele max sac DBN*.sac'
    write(6,*) ' '
    write(6,*) ' - CCTS processing (max. 24 hours of data)'
    write(6,*) ' - Timebins of 10.0 s, with 50% overlap = 5.0 s'
    write(6,*) ' - Factor 2 frequency domain oversampling'
    write(6,*) ' - The beamforming resolution is 5 m/s and 2 deg.'
    write(6,*) ' - Trace velocity range: 300-450 m/s; bearing range: 0-360 deg.'
    write(6,*) ' - Only one maximum correlation maximum per timebin is printed.'
    write(6,*) ' - Binary SAC input is expected, files DBN*.sac'
    write(6,*) ''
    write(6,*) 'The output file <ccts.dat> format is:'
    write(6,*) ' '    
    write(6,*) '   1. Timebin'
    write(6,*) '   2. Maximum of cross-correlation (CC traces)'
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
    write(6,*) '  13. Fisher-ratio (waveform data)'
    write(6,*) ' '
    write(6,*) ' Bestbeam output is written out as SAC file <bestbeam.sac>'
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    call exit(0)
  end subroutine

  end module mod_ccts

!--------------------------------------------------------------------------------------
!-------------------------------START OF MAIN PROGRAM----------------------------------
!--------------------------------------------------------------------------------------
  program ccts
    use mod_ccts
    implicit none

    integer   k, n_instr, alloc_samples, alloc_instr
    integer   bin, n_bins, overlap, start_sample, end_sample
    integer   binsize, fbinsize, i, n_pairs, oversampling, zpd
    double precision tbinsize, fisher, mccm, px, py, bearing, trcvel, prms, p2p, time, fsrate


    type(Slowness) :: s_grid
    type(Frequency) :: freq
    type(WFHeader) :: header

    double precision, allocatable, dimension(:,:) :: r, counts, counts_bin, bestbeam, results
    double precision, allocatable, dimension(:,:) :: coarray, cc_traces
    character(40), allocatable, dimension (:) :: timeseries
    character(40) file_format, tele_local, max_all, fid_beam, statistic

    statistic = 'fisher'
    alloc_instr = 50
    alloc_samples = 2*24*3600*250

    allocate(timeseries(alloc_instr))
    allocate(r(alloc_instr,3))

    call get_cmd_parameters(n_instr,tbinsize,overlap,oversampling,s_grid,r,tele_local,max_all,file_format,timeseries)
    call span_slowness_grid(s_grid,tele_local)

    allocate(counts(alloc_samples,n_instr))
    allocate(bestbeam(alloc_samples,1))
    allocate(results(s_grid%n_beams,3))

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

    ! Add 10% zeropadding to accomodate beamforming of the cross-correlation traces
    zpd = fbinsize/10

    n_pairs = 0
    do i=n_instr-1,1,-1
     n_pairs = n_pairs + i
    end do

    allocate(counts_bin(binsize,n_instr))
    allocate(cc_traces(fbinsize+2*zpd,n_pairs))
    allocate(coarray(n_pairs,2))

    write(6,*) ''
    write(6,*) ''
    write(6,*) 'Cross-correlation Trace Stacking (CCTS) analysis ...'
    write(6,*) '----------------------------------------------------'
    write(6,*) ''

    write(6,'(a32,f8.2,a4)') ' Time-domain sampling : ', header%srate, ' Hz'
    write(6,'(a32,f8.2,a4)') ' Time bin length      : ', tbinsize, ' s'
    write(6,'(a32,f8.2,a4)') ' Time bin overlap     : ', overlap/100.*tbinsize, ' s'
    write(6,'(a32,f8.2,a4)') ' Interpolation factor : ', 1.0*fbinsize/binsize, ' x'
    write(6,'(a32,f8.2,a4)') ' Processing sampling  : ', fsrate, ' Hz'

    write(6,*) ''
    write(6,'(a32,f5.1,a4)') ' Trace velocity resolution : ', s_grid%dc, ' m/s'
    write(6,'(a32,f5.1,a4)') ' Bearing resolution : ', s_grid%dth, ' deg'
    write(6,*) ''

    write(6,'(a32,f8.2,a3,f8.2,a6)') ' Trace velocity domain : [ ', s_grid%trcvel_min, ' - ', s_grid%trcvel_max, ' ] m/s'
    write(6,'(a32,f8.2,a3,f8.2,a6)') ' Bearing domain : [ ', s_grid%bearing_min, ' - ', s_grid%bearing_max, ' ] deg'
    write(6,*) ''
    write(6,'(a32,a3)') ' Output type : ', max_all

    open (unit=30, file='ccts.dat')
    fid_beam = 'bestbeam.sac'

    overlap = binsize - binsize*overlap/100.
    n_bins = (binsize/overlap)*(header%n_samples/binsize)
    
    call get_co_array(r,n_instr,coarray)

    do bin = 2, n_bins - 1
      start_sample = (bin-1)*overlap+1
      end_sample   = start_sample + (binsize-1)
      time         = header%tzero + ((start_sample + end_sample-1) / 2.0 ) / header%srate

      counts_bin(1:binsize,1:n_instr) = counts(start_sample:end_sample,1:n_instr)
      call window_data(n_instr,binsize,counts_bin)
      call get_cc_traces(counts_bin,binsize,fbinsize,zpd,fsrate,n_instr,cc_traces)

      call beamgrid(statistic,zpd+1,fsrate,fbinsize,s_grid,n_pairs,coarray,cc_traces,results)

      if (max_all == 'max') then
       call beamgrid_maximum(results,px,py,mccm)
       call convert_slowness(px,py,bearing,trcvel)
       call compute_f_ratio(start_sample,binsize,n_instr,header%srate,counts,r,px,py,fisher)
       call get_bestbeam(bin,start_sample,header,counts,binsize,n_instr,overlap,r,px,py,bestbeam,prms,p2p,freq)

       write (30,98) time, mccm, bearing, trcvel, px, py, prms, p2p, n_instr, freq%center, freq%bandwith, freq%rms, fisher
       98 format (f10.2, 1x, f10.2, 1x, f8.2, 1x, f8.2, 1x, e15.8, 1x, e15.8, 1x, e15.8, 1x, e15.8, 1x, i2, &
&                                                                  1x, f8.3, 1x, f8.3, 1x, f8.3, 1x, f10.2 )

      elseif (max_all == 'all') then
       do k=1,s_grid%n_beams
        mccm = results(k,1)
        px = results(k,2)
        py = results(k,3)
        call convert_slowness(px,py,bearing,trcvel)
        !call get_bestbeam(bin,start_sample,header,counts,binsize,n_instr,overlap,r,px,py,bestbeam,prms,freq)

        write (30,99) time, mccm, bearing, trcvel, px, py, n_instr
        99 format (f10.2, 1x, f10.2, 1x, f8.2, 1x, f8.2, 1x, e15.8, 1x, e15.8, 1x, i2)
       end do
       write (30,*) ''
      endif
    end do

    if (max_all == 'max') then
     call write_sac(fid_beam,header,bestbeam)
    endif

    close(30)

    deallocate(r)
    deallocate(timeseries)
    deallocate(counts)
    deallocate(coarray)
    deallocate(s_grid%px)
    deallocate(s_grid%py)
    deallocate(cc_traces)
    deallocate(counts_bin)
    deallocate(bestbeam)
    deallocate(results)

  end program ccts
