c*************************************************************************
c                            TU4_STEP.F
c*************************************************************************
c This subroutine has same i/o as STEP_KDK but here we use a 4th
c order T + U decomposition of the Hamiltonian like the old symplectic
c method (4th order) used by Gladman and Duncan.
c Steps both massive and test particles in the same subroutine.
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c
c Remarks:  Based on Martin's NB4M routines and Hal's STEP_KDK
c Authors:  Martin Duncan 
c Date:    3/8/93
c Last revision: 2/24/94

      subroutine tu4_step(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer j
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      real*8 axb(NPLMAX),ayb(NPLMAX),azb(NPLMAX)
      real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
      real*8 axbt(NTPMAX),aybt(NTPMAX),azbt(NTPMAX)
      real*8 msys

      real*8 w0,w1,asymp(4),bsymp(4),adt,bdt

      save asymp,bsymp
c----
c...  Executable code 

c If first time, load up the symplectic coeffs asymp(*) and bsymp(*)

      if (i1st .eq.0) then

         w1= 1.d0/(2.d0 -(2.d0**(1.d0/3.d0)))
         w0= 1.d0 - 2.d0*(w1)
         asymp(1) = 0.5*w1
         asymp(4) = asymp(1)
         asymp(2) = 0.5*(w0 + w1)
         asymp(3) = asymp(2)
         bsymp(1) = 0.0d0
         bsymp(2) = w1
         bsymp(3) = w0
         bsymp(4) = w1

         i1st = 1
      endif

c...  Convert to barycentric coords
      call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &     xb,yb,zb,vxb,vyb,vzb,msys)
      call coord_h2b_tp(ntp,xht,yht,zht,vxht,vyht,vzht,
     &     xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &     xbt,ybt,zbt,vxbt,vybt,vzbt)
      
c Take a drift forward to begin the leapfrog (the first kick = 0).
      adt = asymp(1)*dt

      call tu4_ldrift(nbod,xb,yb,zb,vxb,vyb,vzb,adt)
      call tu4_ldrift_tp(ntp,xbt,ybt,zbt,vxbt,vybt,vzbt,
     &     adt,istat)

c Get the accelerations, kick and drift three times
      do j=2,4
         call tu4_getaccb(nbod,mass,j2rp2,j4rp4,xb,yb,zb,axb,ayb,azb)
         call tu4_getaccb_tp(nbod,mass,j2rp2,j4rp4,xb,yb,zb,
     &        ntp,xbt,ybt,zbt,istat,axbt,aybt,azbt)

         bdt = bsymp(j)*dt
         call tu4_vkickb(nbod,vxb,vyb,vzb,axb,ayb,azb,bdt)
         call tu4_vkickb_tp(ntp,vxbt,vybt,vzbt,axbt,aybt,
     &        azbt,bdt,istat)

         adt = asymp(j)*dt
         call tu4_ldrift(nbod,xb,yb,zb,vxb,vyb,vzb,adt)
         call tu4_ldrift_tp(ntp,xbt,ybt,zbt,vxbt,vybt,
     &            vzbt,adt,istat)
      enddo

c...  Convert back to helio. coords at the end of the step
      call coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &     xh,yh,zh,vxh,vyh,vzh)
      call coord_b2h_tp(ntp,xbt,ybt,zbt,vxbt,vybt,vzbt,
     &     xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &     xht,yht,zht,vxht,vyht,vzht)
      
      return
      
      end   ! step_tu4
c------------------------------------------------------------------------
