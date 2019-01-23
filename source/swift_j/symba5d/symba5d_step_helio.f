c*************************************************************************
c                            SYMBA5D_STEP_HELIO.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a LINEAR DRIFT, then a KICK, then a DRAG, then a KEPLER DRIFT,
c then a DRAG, then a KICK, and finally a LINEAR DRIFT.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c                 rpl           ==>  physical size of a planet
c                                    (real array)
c                 kdrag0,nkdrag ==>  Drag constant and exponent
c                                    (real scalars)
c                 eta0,neta     ==>  Gas velocity parameter and exponent
c                                    (real scalars)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c Remarks: Based on helio_step_pl.f but does not pass the intermediate
c          positions and velocities back for the TP to use.
c Authors:  Hal Levison 
c Date:    3/20/97
c Last revision: 12/22/99

      subroutine symba5d_step_helio(i1st,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,rpl,kdrag0,nkdrag,eta0,neta)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,nbodm
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 rpl(nbod),kdrag0,nkdrag,eta0,neta

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Internals:
      real*8 dth 
      real*8 axh(NTPMAX),ayh(NTPMAX),azh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX),msys
      real*8 ptxb,ptyb,ptzb            ! Not used here
      real*8 ptxe,ptye,ptze

      save axh,ayh,azh     ! Note this !!
      save vxb,vyb,vzb     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  Apply gas drag for a half dt with both x and v in helio frame
      call helio_drag(nbod,mass,xh,yh,zh,vxh,vyh,vzh,dth,rpl,
     &     kdrag0,nkdrag,eta0,neta)
      call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxb,ptyb,ptzb)

c...  Get the accelerations in helio frame. if frist time step
      call symba5_helio_getacch(nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c..   Drift in helio coords for the full step 
      call helio_drift(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)

c...  Get the accelerations in helio frame. if frist time step
      call symba5_helio_getacch(nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxe,ptye,ptze)

c...  Apply gas drag for a half dt with both x and v in helio frame
      call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)
      call helio_drag(nbod,mass,xh,yh,zh,vxh,vyh,vzh,dth,rpl,
     &     kdrag0,nkdrag,eta0,neta)

      return
      end   ! symba5d_step_helio
c---------------------------------------------------------------------

