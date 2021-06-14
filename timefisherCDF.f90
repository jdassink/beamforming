  ! Time-Domain Fisher Detector with command line parameters
  ! based on "Multiple Signal Correlators", by Melton & Bailey, Geophysics 1957

  ! Author     : Jelle D. Assink (assink@knmi.nl)
  ! Date       : November 2010
 
  module td_fisherCDF
  use string_utility
  use sa_array
  use netcdf
  implicit none

  contains 

  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

  subroutine get_cmd_parameters(n_instr,tbinsize,overlap,t0_epoch,s_grid,r,file_format,timeseries)
    integer i, n_instr, eof, overlap, iargc, t0_epoch
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
    read(dummy,*) t0_epoch
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
    dummy = StrLowCase(dummy)
    read(dummy,*) s_grid%tele_local
    call getarg(12, dummy)
    dummy = StrLowCase(dummy)
    read(dummy,*) s_grid%max_all
    call getarg(13, dummy)
    dummy = StrLowCase(dummy)
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
    write(6,*) 'Time-domain Fisher detector and beamformer for ASCII/SAC formatted data'
    write(6,*) '________________________________________________________________________'
    write(6,*) ''
    write(6,*) 'Usage:'
    write(6,*) '------'
    write(6,*) 'timefisherCDF <coordinate file> <binsize> <overlap> <t0_epoch>'
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
    write(6,*) 'timefisherCDF stations 10.0 50 0 0.0 360.0 2.0 300.0 450.0 5.0 tele max sac DBN*.sac'
    write(6,*) ' '
    write(6,*) ' - Time-domain Fisher processing (max. 24 hours of data)'
    write(6,*) ' - Timebins of 10.0 s, with 50% overlap = 5.0 s'
    write(6,*) ' - Reference time is 0 seconds since 1970-01-01 00:00:00'
    write(6,*) ' - The beamforming resolution is 5 m/s and 2 deg.'
    write(6,*) ' - Trace velocity range: 300-450 m/s; bearing range: 0-360 deg.'
    write(6,*) ' - Only one maximum Fisher per timebin is printed.'
    write(6,*) ' - '
    write(6,*) ' - Binary SAC input is expected, files DBN*.sac'
    write(6,*) ''
    write(6,*) ' The output file is <timefisher.nc>'
    write(6,*) '________________________________________________________________________'
    write(6,*) ' '
    call exit(0)
  end subroutine

  end module

!--------------------------------------------------------------------------------------

  program timefisherCDF
    use td_fisherCDF
    implicit none

    integer   k, n_instr, alloc_samples, alloc_instr, ci, bi
    integer   bin, n_bins, binsize, overlap, start_sample, end_sample, t0_epoch
    double precision fisher, px, py, bearing, trcvel, prms, p2p, tbinsize

    type(Slowness) :: s_grid
    type(Frequency) :: freq
    type(WFHeader) :: header

    double precision, allocatable, dimension(:) :: results, time
    double precision, allocatable, dimension(:,:) :: r, counts, bestbeam
    double precision, allocatable, dimension(:,:,:) :: resultsCDF
    character(40), allocatable, dimension (:) :: timeseries
    character(40) file_format, fid_beam, statistic
    character*200 x_m, y_m, z_m

    ! NetCDF stuff
    character (len = *), parameter :: nc_file = 'timefisher.nc'
    character (len = *), parameter :: UNITS = "units"
    character (len = *), parameter :: F_UNITS = "-"
    character (len = *), parameter :: X_UNITS = "degrees"
    character (len = *), parameter :: Y_UNITS = "m/s"
    character (len = *), parameter :: T_UNITS = "seconds"

    character (len = *), parameter :: F_NAME="Fisher_ratio"
    character (len = *), parameter :: X_NAME="back_azimuth"
    character (len = *), parameter :: Y_NAME="trace_velocity"
    character (len = *), parameter :: T_NAME="time"

    integer, parameter :: NDIMS = 3
    ! When we create netCDF files, variables and dimensions, we get back an ID for each one.
    integer :: ncid, x_dimid, y_dimid, t_dimid, dimids(NDIMS)
    integer :: f_varid, x_varid, y_varid, t_varid
    call check( nf90_create(nc_file, NF90_CLOBBER, ncid) )

    statistic = 'fisher'
    alloc_instr = 50
    alloc_samples = 2*24*3600*250

    allocate(timeseries(alloc_instr))
    allocate(r(alloc_instr,5))

    call get_cmd_parameters(n_instr,tbinsize,overlap,t0_epoch,s_grid,r,file_format,timeseries)
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

    write(6,'(a32,i12,a4)')  ' Reference epoch      : ', t0_epoch, ' s'
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

    !open (unit=30, file='timefisher.dat')
    fid_beam = 'bestbeam.sac'
    
    binsize = nint(tbinsize*header%srate)
    overlap = binsize - binsize*overlap/100.
    n_bins = (binsize/overlap)*(header%n_samples/binsize)

    ! Prepare an array for results in NetCDF format
    allocate(resultsCDF(s_grid%n_b,s_grid%n_c,n_bins))
    allocate(time(n_bins))

    do bin = 2, n_bins - 1
      start_sample  = (bin-1)*overlap+1
      end_sample    = start_sample + (binsize-1)
      time(bin)     = header%tzero + ((start_sample + end_sample-1) / 2.0 ) / header%srate

      call beamgrid(statistic,start_sample,header%srate,binsize,s_grid,n_instr,r,counts,results)

      if (s_grid%max_all .eq. 'max') then
       call beamgrid_maximum(results,k)
       px = s_grid%px(k)
       py = s_grid%py(k)

       call get_bestbeam(bin,start_sample,header,counts,binsize,n_instr,overlap,r,px,py,bestbeam,prms,p2p,freq)

!       write (30,98) time(bin), results(k), s_grid%bearing(s_grid%bi(k)), s_grid%trcvel(s_grid%ci(k)), px, py,  &
!  &                        prms, p2p, n_instr, freq%center, freq%bandwith, freq%rms
!       98 format (f10.2, 1x, f10.2, 1x, f8.2, 1x, f8.2, 1x, e15.8, 1x, e15.8, 1x,   &
!  &                        e15.8, 1x, e15.8, 1x, i2, 1x, f8.3, 1x, f8.3, 1x, f8.3)

      elseif (s_grid%max_all .eq. 'all') then
        resultsCDF(:,:,bin) = reshape(results, (/s_grid%n_b, s_grid%n_c/))
      endif
    
    end do

    call check(nf90_put_att(ncid, NF90_GLOBAL, "algorithm", "Time-domain Fisher beamforming"))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "reference_1", "Melton & Bailey, 1957 - https://doi.org/10.1190/1.1438390"))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "reference_2", "Evers, 2008 - https://tinyurl.com/yaxgeoet"))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "author", "assink@knmi.nl"))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "array_latitude_deg", r(1,4)))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "array_longitude_deg", r(1,5)))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "array_number_of_elements", n_instr))
    write (x_m,'(*(f10.2,A))') (r(k,1), ", ",k=1,n_instr)
    write (y_m,'(*(f10.2,A))') (r(k,2), ", ",k=1,n_instr)
    write (z_m,'(*(f10.2,A))') (r(k,3), ", ",k=1,n_instr)
    call check(nf90_put_att(ncid, NF90_GLOBAL, "x_m", trim(x_m)))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "y_m", trim(y_m)))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "z_m", trim(z_m)))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "sampling_rate_Hz", nint(header%srate)))

    call check(nf90_put_att(ncid, NF90_GLOBAL, "start_time_epoch_s", t0_epoch))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "binsize_s", tbinsize))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "overlap_percentage", 100-nint(100.0*overlap/binsize)))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "trace_velocity_minimum_mps", s_grid%trcvel_min))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "trace_velocity_maximum_mps", s_grid%trcvel_max))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "trace_velocity_resolution_mps", s_grid%dc))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "back_azimuth_minimum_deg", s_grid%bearing_min))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "back_azimuth_maximum_deg", s_grid%bearing_max))
    call check(nf90_put_att(ncid, NF90_GLOBAL, "back_azimuth_resolution_deg", s_grid%dth))

    call check( nf90_def_dim(ncid, "backazi",  s_grid%n_b, x_dimid) )
    call check( nf90_def_dim(ncid, "trcvel" ,  s_grid%n_c, y_dimid) )
    call check( nf90_def_dim(ncid, "time"   ,      n_bins, t_dimid) )
    dimids =  (/ x_dimid, y_dimid, t_dimid /)

    call check( nf90_def_var(ncid, X_NAME, NF90_FLOAT, x_dimid, x_varid) )
    call check( nf90_def_var(ncid, Y_NAME, NF90_FLOAT, y_dimid, y_varid) )
    call check( nf90_def_var(ncid, T_NAME, NF90_FLOAT, t_dimid, t_varid) )
    call check( nf90_def_var(ncid, F_NAME, NF90_FLOAT,  dimids, f_varid) )

    call check( nf90_put_att(ncid, x_varid, UNITS, X_UNITS) )
    call check( nf90_put_att(ncid, y_varid, UNITS, Y_UNITS) )
    call check( nf90_put_att(ncid, t_varid, UNITS, T_UNITS) )
    call check( nf90_put_att(ncid, f_varid, UNITS, F_UNITS) )

    call check( nf90_enddef(ncid) )

    if (s_grid%max_all .eq. 'max') then
      call write_sac(fid_beam,header,bestbeam)
      
    elseif (s_grid%max_all .eq. 'all') then
      call check( nf90_put_var(ncid, t_varid, time) )
      call check( nf90_put_var(ncid, x_varid, s_grid%bearing) )
      call check( nf90_put_var(ncid, y_varid, s_grid%trcvel) )
      call check( nf90_put_var(ncid, f_varid, resultsCDF) )
    endif

    call check( nf90_close(ncid) )
    !close(30)

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
    deallocate(resultsCDF)

  end program timefisherCDF
