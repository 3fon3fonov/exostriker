c************************************************************************
c                         IO_DUMP_PL.F
c************************************************************************
c Dumps the data for the Sun and planets 

c
c             Input:
c                 dplfile       ==>  Name of file to write to (character*80)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in Helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial position in Helio coord 
c                                    (real arrays)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 iflgchk       ==>  bit 5 set ==>  include J2 and J4 terms
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  2/22/94 HFL

	subroutine io_dump_pl(dplfile,nbod,mass,xh,yh,zh,
     &     vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	real*8 mass(NPLMAX),rplsq(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
	integer nbod, iflgchk
	character*(*) dplfile
        logical*2 lclose

c...   Internal
	integer j,ierr
        real*8 rpl

c-----
c...  Executable code      

        call io_open(7,dplfile,'unknown','formatted',ierr)

	write(7,*) nbod

        if(btest(iflgchk,5))  then ! bit 5 is set
           write(7,123) mass(1),j2rp2,j4rp4
        else
           write(7,123) mass(1)
        endif
        write(7,123) xh(1),yh(1),zh(1)
        write(7,123) vxh(1),vyh(1),vzh(1)

	do j=2,nbod
           if(lclose) then
              rpl = sqrt(rplsq(j))
              write(7,123) mass(j),rpl
           else
              write(7,123) mass(j)
           endif
           write(7,123) xh(j),yh(j),zh(j)
           write(7,123) vxh(j),vyh(j),vzh(j)
	enddo
 123    format(3(1p1e23.16,1x))

	close(unit = 7)
	return
	end    ! io_dump_pl.f
c--------------------------------------------------------------------------

