c**********************************************************************
c			IO_INIT_TP.F
c**********************************************************************
c Read in test particle data
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c
c             Output:
c                 ntp           ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c              xht,yht,zht      ==>  initial position in Helio coord 
c                                    (real arrays)
c              vxht,vyht,vzht   ==>  initial position in Helio coord 
c                                    (real arrays)
c               istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0  active
c                                      istat(i,1) = 1 not
c               rstat           ==>  status of the test paricles
c                                      (2d  real array)
c
c
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  12/22/95  HFL

	subroutine io_init_tp(infile,ntp,xht,yht,zht,vxht,vyht,
     &     vzht,istat,rstat)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile

c...    Output
	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
        real*8 rstat(NTPMAX,NSTATR)
	integer istat(NTPMAX,NSTAT)
	integer ntp

c...   Internal
	integer i,j,ierr,ns

c-----
c...  Executable code      

	write(*,*) 'Test particle file called ',infile
        call io_open(7,infile,'old','formatted',ierr)

	read(7,*) ntp

        if(ntp.gt.NTPMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_tp: '
           write(*,*) '   The number of test bodies,',ntp,','
           write(*,*) '   is too large, it must be less than',NTPMAX
           call util_exit(1)
        endif

	write(*,*) ' '
	write(*,*) 'ntp : ',ntp

        if(ntp.eq.0) then
           close(unit = 7)
           write(*,*) ' '
           return
        endif               ! <===== NOTE

c...   Determine the number of istat and rstat variables.  In what follows,
c...   we assume that they are the same.

        call io_getns(7,ns)
        
        if(ns.ne.NSTAT) then
           write(*,*) 'Warning:  The size of istat and rstat arrays is '
           write(*,*) '          not NSTAT=',NSTAT,', but is ',ns
        endif

c Start again:
        rewind(7)
        read(7,*) ntp 

c Read in the x's and v's and istat(*,*)
	  write(*,*) ' '
	  do  i=1,ntp
	    read(7,*) xht(i),yht(i),zht(i)
	    read(7,*) vxht(i),vyht(i),vzht(i)
	    read(7,*) (istat(i,j),j=1,ns)
	    read(7,*) (rstat(i,j),j=1,ns)
            do j=ns+1,NSTAT
               istat(i,j) = 0
            enddo
            do j=ns+1,NSTATR
               rstat(i,j) = 0.0d0
            enddo
	  enddo

	close(unit = 7)
        write(*,*) ' '

	return
	end    ! io_init_tp.f
c-----------------------------------------------------------------

