c**********************************************************************
c			IO_DUMP_TP.F
c**********************************************************************
c Dump test particle data
c
c             Input:
c                 dtpfile       ==>  Name of file to write to (character*80)
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
c Last revision: 2/25/94 HFL

	subroutine io_dump_tp(dtpfile,ntp,xht,yht,zht,vxht,vyht,
     &          vzht,istat,rstat)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
        real*8 rstat(NTPMAX,NSTATR)
	integer istat(NTPMAX,NSTAT)
	integer ntp
	character*(*) dtpfile

c...   Internal
	integer i,j,ierr

c-----
c...  Executable code      

        call io_open(7,dtpfile,'unknown','formatted',ierr)

	write(7,*) ntp
	do  i=1,ntp
	    write(7,123) xht(i),yht(i),zht(i)
	    write(7,123) vxht(i),vyht(i),vzht(i)
	    write(7,*) (istat(i,j),j=1,NSTAT)
	    write(7,123) (rstat(i,j),j=1,NSTATR)
	enddo
123	format(4(1p1e23.16,1x))

	close(unit = 7)

	return
	end    ! io_dump_tp
c-----------------------------------------------------------------

