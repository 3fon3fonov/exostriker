c*************************************************************************
c                            IO_DISCARD_MASS
c*************************************************************************
c Write out information about a discarded massive body.
c
c             Input:
c                 init          ==>  initiize flag if = 0 initialize and return
c                                                     = 1 run through 
c                 id            ==> particle number (int scalar)
c                 time          ==>  current time (real scalar)
c                 m1            ==>  Mass of pl (real scalar)
c                 r1            ==>  Radius of pl 2 (real scalar)
c                 x1,y1,z1      ==>  current position of pl 1 in helio coord 
c                                    (real scalar)
c                 vx1,vy1,vz1   ==>  current velocity of pl 1 in helio coord 
c                                    (real scalar)
c                 iu            ==> IO unit (int scalar)
c                 iwhy          ==> reason for discard (int scalar)
c                 fopenstat     ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c Remarks: 
c Authors:  Hal Levison 
c Date:    12/30/96
c Last revision: 

      subroutine io_discard_mass(init,time,id,m1,r1,x1,y1,z1,vx1,vy1,
     &     vz1,iu,iwhy,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iwhy,iu,init,id
      real*8 time
      real*8 m1,r1
      real*8 x1,y1,z1
      real*8 vx1,vy1,vz1
      character*(*) fopenstat

c...  Internals
      integer ierr

c----
c...  Executable code 

      if(init.eq.0) then

         call io_open(iu,'discard_mass.out',fopenstat,
     &        'FORMATTED',ierr)

c...     if there was an error and fopenstat='append' then
c...     try to open as new
         if(ierr.ne.0) then  
            if( (fopenstat(1:6).eq.'append') .or. 
     &           (fopenstat(1:6).eq.'APPEND') ) then
               call io_open(iu,'discard_mass.out','new','FORMATTED',
     &              ierr)
            endif
         endif

         if(ierr.ne.0) then
            write(*,*) ' SWIFT ERROR: in io_discard_mass: '
            write(*,*) '    Could not open discard output file'
            call util_exit(1)
         endif
         return   ! <=== NOTE!!!!
      else
         call io_open(iu,'discard_mass.out','append','FORMATTED',
     &        ierr)
      endif

      write(iu,1000) time,iwhy
 1000 format(1x,1p1e23.16,1x,i4)

      write(iu,2000) id,m1,r1
 2000 format('-1',1x,i3,1x,2(1p1e23.16,1x))
      write(iu,3000) x1,y1,z1
 3000 format(3(1p1e23.16,1x))
      write(iu,3000) vx1,vy1,vz1

      close(unit = iu)
      return
      end                       ! io_discard_mass.f
c--------------------------------------------------------------------------

