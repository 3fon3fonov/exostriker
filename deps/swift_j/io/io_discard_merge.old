c*************************************************************************
c                            IO_DISCARD_MERGE
c*************************************************************************
c Write out information about a merger.
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 ip1,ip2       ==>  planets to merge (real scalar)
c                 m1            ==>  Mass of pl 1 (real scalar)
c                 r1            ==>  Radius of pl 1 (real scalar)
c                 x1,y1,z1      ==>  current position of pl 1 in helio coord 
c                                    (real arrays)
c                 vx1,vy1,vz1   ==>  current velocity of pl 1 in helio coord 
c                                    (real arrays)
c                 m2            ==>  Mass of pl 2 (real scalar)
c                 r2            ==>  Radius of pl 2 (real scalar)
c                 x2,y2,z2      ==>  current position of pl 2 in helio coord 
c                                    (real arrays)
c                 vx2,vy2,vz2   ==>  current velocity of pl 2 in helio coord 
c                                    (real arrays)
c                 mn            ==>  Mass of new pl  (real scalar)
c                 rn            ==>  Radius of new pl (real scalar)
c                 xn,yn,zn      ==>  current position of new pl in helio coord 
c                                    (real arrays)
c                 vxn,vyn,vzn   ==>  current velocity of new pl in helio coord 
c                                    (real arrays)
c                 nleft           ==>  number of active test bodies(int scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    12/30/96
c Last revision: 

      subroutine io_discard_merge(time,ip1,ip2,m1,r1,x1,y1,z1,vx1,vy1,
     &     vz1,m2,r2,x2,y2,z2,vx2,vy2,vz2,mn,rn,xn,yn,zn,vxn,vyn,vzn)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer ip1,ip2
      real*8 time
      real*8 m1,r1
      real*8 x1,y1,z1
      real*8 vx1,vy1,vz1
      real*8 m2,r2
      real*8 x2,y2,z2
      real*8 vx2,vy2,vz2
      real*8 mn,rn
      real*8 xn,yn,zn
      real*8 vxn,vyn,vzn

c...  Internals
      integer ierr,iu

c----
c...  Executable code 

      iu = 40

      call io_open(iu,'discard_mass.out','append','FORMATTED',ierr)

      write(iu,1000) time
 1000 format(1x,1p1e23.16,'  2')

      write(iu,2000) ip1,m1,r1
 2000 format('-1',1x,i3,1x,2(1p1e23.16,1x))
      write(iu,3000) x1,y1,z1
 3000 format(3(1p1e23.16,1x))
      write(iu,3000) vx1,vy1,vz1

      write(iu,2000) ip2,m2,r2
      write(iu,3000) x2,y2,z2
      write(iu,3000) vx2,vy2,vz2

      write(iu,4000) ip1,mn,rn
 4000 format('+1',1x,i3,1x,2(1p1e23.16,1x))
      write(iu,3000) xn,yn,zn
      write(iu,3000) vxn,vyn,vzn

      close(unit = iu)
      return
      end                       ! io_discard_merge.f
c--------------------------------------------------------------------------

