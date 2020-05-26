c************************************************************************
c			IO_LYAP_INIT
c************************************************************************
c Get data for displacement vectors of shadow particles and 
c use this data with the tp positions and vels. to return with
c initialized values for their posns, vels, initial phase space distances and
c the quantities lrsum(*).
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c                 ntp           ==>  number of massive bodies (int scalar)
c              xht,yht,zht      ==>  initial position of TP in Helio coord 
c                                    (real arrays)
c              vxht,vyht,vzht   ==>  initial position of TP in Helio coord 
c                                    (real arrays)
c
c             Output:
c              xsh,ysh,zsh      ==>  initial position of shawdow part in Helio coord 
c                                    (real arrays)
c              vxsh,vysh,vzsh   ==>  initial position of SP in Helio coord 
c                                    (real arrays)
c                   dist0       ==>  initial separation between TP and SP(real array)
c                   lrsum       ==>  initial (0.0) run of the log of gamma
c                                    (real array)
c                   dtnorm      ==>  How often the dist is renoramized
c                                    (real scalar)
c                   iul         ==>  Unit that lyap will write to
c
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    5/5/93 
c Last revision:  2/21/94  HFL

        subroutine io_lyap_init(infile,ntp,xht,yht,zht,vxht,vyht,vzht,
     &                xsh,ysh,zsh,vxsh,vysh,vzsh,dist0,lrsum,dtnorm,iul)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	integer ntp,iul
	character*(*) infile
	real*8 xht(ntp),yht(ntp),zht(ntp)
	real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...    Output
	real*8 xsh(ntp),ysh(ntp),zsh(ntp)
	real*8 vxsh(ntp),vysh(ntp),vzsh(ntp)
	real*8 dist0(ntp),lrsum(ntp),dtnorm

c...   Internal
	integer nsh,n,ierr
	real*8 xeps,yeps,zeps,vxeps,vyeps,vzeps

c-----
c...  Executable code

	write(*,*) 'Shadow particle data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)
        write(*,*) ' '

	read(7,*) nsh,dtnorm
	if(nsh.ne.ntp) then
	  write(*,*) '***STOP !!! - nsh is not equal to ntp ***'
	  call util_exit(1)
	endif

	do n=1,ntp
	  read(7,*) xeps,yeps,zeps,vxeps,vyeps,vzeps
	  xsh(n) = xht(n) + xeps
	  ysh(n) = yht(n) + yeps
	  zsh(n) = zht(n) + zeps
	  vxsh(n) = vxht(n) + vxeps
	  vysh(n) = vyht(n) + vyeps
	  vzsh(n) = vzht(n) + vzeps
	  dist0(n)= sqrt((xht(n) - xsh(n))**2 + (yht(n) -ysh(n))**2
     &                  +(zht(n) - zsh(n))**2  +(vxht(n) - vxsh(n))**2
     &                 +(vyht(n) - vysh(n))**2 +(vzht(n) - vzsh(n))**2)
	  lrsum(n) =0.d0
	enddo

	close(unit=7)

        call io_open(iul,'lyap.out','new','formatted',ierr)
        if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in lyap_step: '
           write(*,*) '    Could not open lyap output file',ierr
           call util_exit(1)
        endif
        close(iul)

	return
	end    ! io_lyap_init
c-------------------------------------------------------------------------
