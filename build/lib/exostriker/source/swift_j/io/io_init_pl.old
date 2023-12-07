c************************************************************************
c                              IO_INIT_PL.F
c************************************************************************
c IO_INIT_PL reads in the data for the Sun and planets 
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 iflgchk        ==>  bit 5 set ==>  include J2 and J4 terms
c
c             Output:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in Helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial position in Helio coord 
c                                    (real arrays)
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision: 3/3/94  HFL

	subroutine io_init_pl(infile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &     vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile
	integer iflgchk
        logical*2 lclose

c...    Output
	real*8 mass(NPLMAX),rplsq(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
	integer nbod

c...   Internal
	integer j,ierr
        real rpl

c-----
c...  Executable code      

	write(*,*) 'Planet data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)

c Read number of planets
	read(7,*) nbod

        if(nbod.gt.NPLMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_pl: '
           write(*,*) '   The number of massive bodies,',nbod,','
           write(*,*) '   is too large, it must be less than',NPLMAX
           call util_exit(1)
        endif

	write(*,23) nbod
23	format(/,'Number of bodies (incl. the Sun) is ',i3,/,
     &   'For each, list mass ',/,
     &   'Followed by x,y,z,vx,vy,vz : '/)

c For each planet read mass, 
c and helioc. position and vel .
        if(btest(iflgchk,5))  then ! bit 5 is set
           read(7,*) mass(1),j2rp2,j4rp4
        else
           read(7,*) mass(1)
           j2rp2 = 0.0d0
           j4rp4 = 0.0d0
        endif
        read(7,*) xh(1),yh(1),zh(1)
        read(7,*) vxh(1),vyh(1),vzh(1)
        write(*,*) mass(1)
        write(*,*) xh(1),yh(1),zh(1)
        write(*,*) vxh(1),vyh(1),vzh(1)
        rplsq(1) = 0.0d0

        if(  (xh(1).ne.0.0d0) .or.
     &       (yh(1).ne.0.0d0) .or.
     &       (zh(1).ne.0.0d0) .or.
     &       (vxh(1).ne.0.0d0) .or.
     &       (vyh(1).ne.0.0d0) .or.
     &       (vzh(1).ne.0.0d0) ) then
           write(*,*) ' SWIFT ERROR: in io_init_pl: '
           write(*,*) '   Input MUST be in heliocentric coordinates '
           write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
           call util_exit(1)
        endif

	do j=2,nbod
           if(lclose) then
              read(7,*) mass(j),rpl
              write(*,*) mass(j),rpl
              rplsq(j) = rpl*rpl
           else
              read(7,*) mass(j)
              write(*,*) mass(j)
           endif
	  read(7,*) xh(j),yh(j),zh(j)
	  read(7,*) vxh(j),vyh(j),vzh(j)
	  write(*,*) xh(j),yh(j),zh(j)
	  write(*,*) vxh(j),vyh(j),vzh(j)
	enddo

	close(unit = 7)
	return
	end     ! io_init_pl.f
c--------------------------------------------------------------------------

