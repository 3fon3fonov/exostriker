c************************************************************************
c			IO_MVITS_INIT
c************************************************************************
c This subroutine reads the TS.IN files
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 dt            ==>  time step
c                 mass          ==>  mass of the massive particles (real array)
c
c             Output:
c                 dtpl               ==> timestep for planet integration.
c                                             (real scalar)
c                 nstep              ==> largest of previous (int scalar)
c                 dtpli              ==> individual timestep for planets.
c                                             (real aray)
c                 massi              ==>  mass of bodies wrt t (2d real array)
c                 nbodi              ==>  number of massive bodies wrt t
c                                           (int array)
c                 idpl               ==>  list of what planets are active wrt t
c                                           (int array)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    12/16/93
c Last revision: 2/21/94

      subroutine io_mvits_init(infile,nbod,dt,mass,dtpl,nstep,dtpli,
     &     massi,nbodi,idpl)

      include '../swift.inc'
      include '../mvits/mvits.inc'
      include 'io.inc'

c...  Input
      integer nbod
      character*(*) infile
      real*8 dt,mass(NPLMAX)

c...  Output
      real*8 dtpl
      real*8 massi(NPLMAX,0:NDTMAX),dtpli(NPLMAX)
      integer idpl(NPLMAX,0:NDTMAX),nbodi(0:NDTMAX),nstep

c...  Internals
      integer i,j,nrat,ierr
      integer nstepi(NPLMAX),nbloc,iflg

c-----
c...  Executable code

      write(*,*) 'Timestep particle data file is ',infile
      call io_open(7,infile,'old','formatted',ierr)
      write(*,*) ' '

      read(7,*) nbloc
      if(nbloc.ne.nbod) then
	 write(*,*) 'ERROR: - wrong number of planets in timestep file'
	 call util_exit(1)
       endif

c...   ignore the first because it is the Sun
       read(7,*) nstepi(1)

       nstep = 1
       do i=2,nbod
          read(7,*) nstepi(i)
          iflg = 0
          do j = 0,NUMDT
            if( nstepi(i) .eq. (2**j)) then
               iflg = 1
            endif
         enddo
         if(iflg.eq.0) then
            write(*,*) 'ERROR: - time step reductions are NOT a 2^n '
            call util_exit(1)
         endif
         dtpli(i) = dt/float(nstepi(i))
         nstep = max(nstep,nstepi(i))
      enddo

      dtpl = dt/float(nstep)

c...  set up 1st time step
      do i=1,nbod
         idpl(i,0) = i
         massi(i,0) = mass(i)
      enddo
      nbodi(0) = nbod

      nstepi(1) = nstep      ! we need this to be true

      do i=1,nstep
         nbodi(i) = 0
         do j=1,nbod
            nrat = nstep/nstepi(j)
            if( mod(i,nrat) .eq. 0 ) then
               nbodi(i) = nbodi(i) + 1
               massi(nbodi(i),i) = mass(j)
               idpl(nbodi(i),i) = j
            endif
         enddo
      enddo

      return
      end      ! io_mvits_init
c-------------------------------------------------------------------------
