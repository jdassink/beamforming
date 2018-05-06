  ! Basic seismo-acoustic IO routines, shared by timefisher and freqfisher
 
  module sa_io
  use string_utility
  use f90sac

  type :: WFHeader
  integer nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec, n_samples
  double precision srate, tzero
  end type

  contains

  subroutine get_timeseries(n_instr,file_format,timeseries,alloc_samples,header,counts)
    implicit none
    integer n_instr, instr, min_samples, alloc_samples
    character(40), allocatable, dimension (:) :: timeseries
    character(40) file_format
    double precision, allocatable, dimension(:,:) :: counts
    type(WFHeader) :: header

    file_format = StrLowCase(file_format)
    min_samples = 0

    do instr=1,n_instr
      if (file_format == "ascii") then
        call read_ascii(timeseries(instr),instr,alloc_samples,header,counts)
      elseif (file_format == "sac") then
        call read_sac(timeseries(instr),instr,alloc_samples,header,counts)
      else
        write (6,'(a33,a10)') 'ERROR: No support for file type: ', file_format
        call exit(-1)
      endif
      
      ! Truncate timeseries to the max number of samples in the file.
      ! Files must start at the same time
      if (instr == 1) then
        min_samples = header%n_samples
      endif
      if (header%n_samples < min_samples) then
        min_samples = header%n_samples
      endif
    end do

    header%n_samples = min_samples
  end subroutine

  subroutine read_ascii(fid,instr,alloc_samples,header,counts)
    implicit none
    character(40) fid
    integer instr, alloc_samples
    double precision    time
    double precision, allocatable, dimension(:,:) :: counts
    LOGICAL :: file_exists
    type(WFHeader) :: header

    write(6,*) 'Reading file : ', fid
    inquire(file=fid, EXIST=file_exists)

    if (file_exists) then
      open (unit=11, file=fid)
      header%n_samples = 0
    
      do while (.true.)
        header%n_samples = header%n_samples + 1
        
        if (header%n_samples-1 > alloc_samples ) then
          write (6,*) ''
          write (6,*) 'ERROR: Too many samples in file. Exiting ...'
          call exit(-1)
        endif
        
        read (11, *, end=999) time, counts(header%n_samples,instr)
        
        if (header%n_samples .EQ. 1) then
          header%tzero = time
        endif
        
        if (header%n_samples .EQ. 2) then
          header%srate = 1.0/(time - header%tzero)
        endif

      end do
      999 continue
      close (unit=11)
      header%n_samples = header%n_samples - 1

    else
      write(6,*) 'ERROR: ASCII timeseries file not found. Exiting ...'
      write(6,*) ''
      call exit(-1)
    endif

  end subroutine

  subroutine read_sac(fid,instr,alloc_samples,header,counts)
    implicit none
    character(40) fid
    double precision, allocatable, dimension(:,:) :: counts
    integer instr, alloc_samples
    type(SACTrace) :: tr
    type(WFHeader) :: header

    write(6,*) 'Reading SAC file : ', fid

    call f90sac_readtrace(fid,tr)
    counts(1:tr%npts,instr) = tr%trace(1:tr%npts)

    header%n_samples = int(tr%npts)
    if (header%n_samples > alloc_samples ) then
      write (6,*) ''
      write (6,*) 'ERROR: Too many samples in file. Exiting ...'
      call exit(-1)
    endif

    header%nzyear = int(tr%nzyear)
    header%nzjday = int(tr%nzjday)
    header%nzhour = int(tr%nzhour)
    header%nzmin  = int(tr%nzmin)
    header%nzsec  = int(tr%nzsec)
    header%nzmsec = int(tr%nzmsec)
    header%tzero  = dble(tr%b)
    header%srate  = 1./dble(tr%delta)
  end subroutine


  subroutine write_sac(fid,header,bestbeam)
    implicit none
    character(40) fid
    double precision, allocatable, dimension(:,:) :: bestbeam
    real delta
    type(SACTrace) :: tr
    type(WFHeader):: header

    delta = real(1./header%srate)
    call f90sac_newtrace(header%n_samples, delta, tr)

    write(6,*) ''
    write(6,*) ''
    write(6,*) 'Writing SAC file : ', trim(fid), ' ...'
    write(6,*) '-----------------------------------'
    write(6,*) ''

    tr%npts   = header%n_samples
    tr%delta  = delta
    tr%nzyear = header%nzyear
    tr%nzjday = header%nzjday
    tr%nzhour = header%nzhour
    tr%nzmin  = header%nzmin
    tr%nzsec  = header%nzsec
    tr%nzmsec = header%nzmsec
    
    !tr%b = 0.
    !tr%e = tr%b + (tr%npts-1)*tr%delta
    tr%kstnm = 'BESTBEAM'
    tr%kcmpnm = 'M'
    tr%trace(1:tr%npts) = real(bestbeam(1:header%n_samples,1))

    call f90sac_writetrace(fid,tr)

  end subroutine

  end module
