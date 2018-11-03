c************************************************************************
c                          IO_DUMP_PARAM.F
c************************************************************************
c IO_DUMP_PARAM dumps out the parameters for the integration. 
c
c      Input:
c       dparfile      ==>  Name of file to write to (character*80)
c            t0       ==> Initial time (real scalar)
c            tstop    ==> final time (real scalar)
c            dt       ==> time step  (real scalar)
c            dtout    ==> time between binary outputs (real scalar)
c            dtdump   ==> time between dumps  (real scalar)
c            iflgchk  ==>  =0 don't run diagnostic routines
c                         !=0 run them
c      rmin,rmax      ==>  maximum and min distance from Sun
c                                if <0  then don't check
c                                    (real scalar)
c       rmaxu         ==>  maximum distance from Sun in not bound
c                                 if <0  then don't check
c                                      (real scalar)
c       qmin          ==> Smallest perihelion distance
c                                 if <0  then don't check
c                                      (real scalar)
c       lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c       outfile       ==>  Name of binary output file (character*80)

c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  5/10/94 HFL

	subroutine io_dump_param(dparfile,t,tstop,dt,dtout,dtdump,
     &           iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)	

	include '../swift.inc'
	include 'io.inc'

c...   Inputs
	real*8 t,tstop,dt
	integer iflgchk
	real*8 dtout,dtdump
	real*8 rmin,rmax,rmaxu,qmin
        logical*2 lclose
	character*(*) outfile,dparfile

c...  Internals
        character*1 lflg(0:IO_NBITS-1),cclose
        integer i,ierr

c-----
c...  Executable code 

c Open parameter data file for the dump
        call io_open(7,dparfile,'unknown','formatted',ierr)

	write(7,*) t,tstop,dt
	write(7,*) dtout,dtdump

        do i=0,IO_NBITS-1
           if(btest(iflgchk,i))  then 
              lflg(i) = 'T'
           else
              lflg(i) = 'F'
           endif
        enddo

        write(7,1000) (lflg(i),i=IO_NBITS-1,0,-1)
 1000   format(100(a1,1x))

        if(btest(iflgchk,4))  then ! bit 4 is set
           if(lclose) then
              cclose = 'T'
           else
              cclose = 'F'
           endif
           write(7,*) rmin,rmax,rmaxu,qmin,' ',cclose
        endif

        if(btest(iflgchk,0).or.btest(iflgchk,1))  then 
           write(7,2000) outfile
 2000      format(a)
        endif

        write(7,2000) 'append'

	close(unit = 7)

	return
	end     ! io_dump_param
c____________________________________________________________________________
c
c




