c************************************************************************
c                              IO_INIT_PL_SYMBA.F
c************************************************************************
c IO_INIT_PL_SYMBA reads in the data for the Sun and planets for 
c symba routines
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl_symba
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
c                 rpl           ==>  physical size of planet
c                                    (real array)
c                 rhill         ==>  size of planet's hills sphere
c                                    (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c
c Remarks: Based on io_init_pl 
c Authors:  Hal Levison
c Date:    11/21/96
c Last revision: 1/10/97

      subroutine io_init_pl_symba(infile,lclose,iflgchk,nbod,mass,
     &     xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

      include '../swift.inc'
      include 'io.inc'
        
c...  Input
      character*(*) infile
      integer iflgchk
      logical*2 lclose

c...  Output
      real*8 mass(NTPMAX),rpl(NTPMAX),j2rp2,j4rp4
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX),rhill(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)
      integer nbod

c...  Internal
      integer j,ierr,ibad
      real*8 r2hill(NTPMAX),rhrat

c-----
c...  Executable code      

      write(*,*) 'Planet data file is ',infile
      call io_open(7,infile,'old','formatted',ierr)

c Read number of planets
      read(7,*) nbod

      if(nbod.gt.NTPMAX) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   The number of massive bodies,',nbod,','
         write(*,*) '   is too large, it must be less than',NTPMAX
         call util_exit(1)
      endif
      
      write(*,23) nbod
 23   format(/,'Number of bodies (incl. the Sun) is ',i3)
      
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
      rpl(1) = 0.0d0
      rhill(1) = 0.0d0
      
      if(  (xh(1).ne.0.0d0) .or.
     &     (yh(1).ne.0.0d0) .or.
     &     (zh(1).ne.0.0d0) .or.
     &     (vxh(1).ne.0.0d0) .or.
     &     (vyh(1).ne.0.0d0) .or.
     &     (vzh(1).ne.0.0d0) ) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   Input MUST be in heliocentric coordinates '
         write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
         call util_exit(1)
      endif

      do j=2,nbod
         if(lclose) then
            read(7,*) mass(j),rhill(j),rpl(j)
         else
            read(7,*) mass(j),rhill(j)
         endif
         read(7,*) xh(j),yh(j),zh(j)
         read(7,*) vxh(j),vyh(j),vzh(j)
      enddo

      close(unit = 7)

c...  check to see if the hills spheres are ok
      call util_hills(nbod,mass,xh,yh,zh,vxh,vyh,vzh,r2hill) 
      ibad = 0
      do j=2,nbod
         rhrat = rhill(j)/sqrt(r2hill(j))
         if( (rhrat.gt.2.0) .or. (rhrat.lt.0.5) ) then
            ibad = ibad + 1
         endif
      enddo

      if(ibad.ne.0) then
         write(*,*) 'Warning in io_init_pl_symba:'
         write(*,*) '   Hill''s spheres are not consistent on ',
     &        ibad,' objects'
      endif
      
      return
      end                       ! io_init_pl_symba.f
c--------------------------------------------------------------------------

