c************************************************************************
c                          IO_INIT_PARAM.F
c************************************************************************
c INIT_PARAM reads in the parameters for the integration. 
c
c      Input:
c            infile   ==> File name to read from (character*80)
c
c      Output:
c            t0       ==> Initial time (real scalar)
c            tstop    ==> final time (real scalar)
c            dt       ==> time step  (real scalar)
c            dtout    ==> time between binary outputs (real scalar)
c            dtdump   ==> time between dumps  (real scalar)
c            iflgchk  ==>  =0 don't run diagnostic routines
c                          bit 0 set ==>  write int*2 binary data file
c                          bit 1 set ==>  write real*4 binary file 
c                          bit 2 set ==>  calc energy of system wrt time
c                          bit 3 set ==>  calc jacobi of the test particles
c                          bit 4 set ==>  check if particles are removed
c                          bit 5 set ==>  include J2 and J4 terms
c      rmin,rmax      ==>  maximum and min distance from Sun
c                                if <0  then don't check
c                                    (real scalar)
c      rmaxu          ==>  maximum distance from Sun in not bound
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
c       fopenstat     ==>  The status flag for the open statements of the
c                          output files.  Must be one of the following:
c                                 new      (die if the file exists)
c                                 append   (add to what is there)
c                                 unknown  (just write over what is there)
c                                 (character*80)
c
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  5/10/94  HFL

        subroutine io_init_param(infile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile

c...  Outputs: 
	integer iflgchk
	real*8 t0,tstop,dt
	real*8 dtout,dtdump
	real*8 rmin,rmax,rmaxu,qmin
        logical*2 lclose
	character*80 outfile,fopenstat

c...  Internals
        logical*1 lflg(0:IO_NBITS-1)
        integer i,ierr

c-----
c...  Executable code 

	write(*,*) 'Parameter data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)

	read(7,*) t0,tstop,dt
	write(*,*) 't0,tstop,dt : ',t0,tstop,dt
	read(7,*) dtout,dtdump
	write(*,*) 'dtout,dtdump : ',dtout,dtdump
        read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

        iflgchk=0
        do i=0,IO_NBITS-1
           if(lflg(i)) then
              iflgchk = ibset(iflgchk,i)
           endif
        enddo

        write(*,*) (lflg(i),i=IO_NBITS-1,0,-1),' = ',iflgchk

        if(btest(iflgchk,0) .and. btest(iflgchk,1))  then 
           write(*,*) ' SWIFT ERROR: in io_init_param:'
           write(*,*) '    Invalid logical flags '
           write(*,*) '    You cannot request that both a real and ',
     &                '       an integer binary file be written '
           call util_exit(1)
        endif

        if(btest(iflgchk,4))  then ! bit 4 is set
           read(7,*) rmin,rmax,rmaxu,qmin,lclose
           write(*,*) 'rmin,rmax,rmaxu,qmin,lclose :',
     &          rmin,rmax,rmaxu,qmin,lclose
        else
           rmin = -1.0
           rmax = -1.0
           rmaxu = -1.0
           qmin = -1.0
           lclose = .false.
        endif

        if(btest(iflgchk,0) .or. btest(iflgchk,1))  then 
           read(7,999) outfile
 999       format(a)
           write(*,*) 'outfile : ', outfile
           write(*,*) ' '
        endif

        read(7,999) fopenstat
        if(  (fopenstat(1:3).ne.'new') .and. 
     &       (fopenstat(1:3).ne.'NEW') .and.
     &       (fopenstat(1:7).ne.'unknown') .and. 
     &       (fopenstat(1:7).ne.'UNKNOWN') .and.
     &       (fopenstat(1:6).ne.'append') .and. 
     &       (fopenstat(1:6).ne.'APPEND') ) then
           write(*,*) ' SWIFT ERROR: in io_init_param:'
           write(*,*) '    Invalid status flag:',fopenstat(1:7),':'
           call util_exit(1)
        endif
        
	close(unit = 7)

	return
	end     ! io_init_param
c____________________________________________________________________________
c
c

