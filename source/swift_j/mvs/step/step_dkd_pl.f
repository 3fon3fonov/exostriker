c*************************************************************************
c                            STEP_DKD_PL.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a DRIFT than a KICK than a DRIFT.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                                    not used here !!!
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xhm,yhm,zhm   ==>  massive part position in middle of step
c                                       (real arrays)
c
c Remarks: Adopted from martin's nbwh.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/12/93

      subroutine step_dkd_pl(i1st,nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &                 vxh,vyh,vzh,dt,xhm,yhm,zhm)	


      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Output Only
      real*8 xhm(nbod),yhm(nbod),zhm(nbod)

c...  Internals:
      integer i
      real*8 dth 
      real*8 axh(NPLMAX),ayh(NPLMAX),azh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
      real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  Convert to jacobi coords
      call coord_h2j(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &          xj,yj,zj,vxj,vyj,vzj)

c...  Take a drift forward 0.5dt to begin the step (the first kick = 0).
      call drift(nbod,mass,xj,yj,zj,vxj,vyj,vzj,dth)

c...  After drift, compute helio. xh and vh for acceleration calculations
      call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &             xh,yh,zh,vxh,vyh,vzh)

c...  save position
      do i=1,nbod
         xhm(i) = xh(i)
         yhm(i) = yh(i)
         zhm(i) = zh(i)
      enddo

c...  Get the accelerations in helio frame.
      call getacch(nbod,mass,j2rp2,j4rp4,xj,yj,zj,xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a full dt 
      call kickvh(nbod,vxh,vyh,vzh,axh,ayh,azh,dt)

c...  Convert the helio. vels. to Jac. vels. (the Jac. positions are unchanged)
      call coord_vh2vj(nbod,mass,vxh,vyh,vzh,vxj,vyj,vzj)

c..   Drift again in Jacobi coords for the final half-step 0.5dt	  
      call drift(nbod,mass,xj,yj,zj,vxj,vyj,vzj,dth)

c...  Convert back to heliocentric
      call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xh,yh,zh,vxh,vyh,vzh)

      return
      end   ! step_dkd_pl
c---------------------------------------------------------------------

