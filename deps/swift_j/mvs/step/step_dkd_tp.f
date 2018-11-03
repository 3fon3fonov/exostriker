c*************************************************************************
c                            STEP_DKD_TP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a DRIFT than a KICK than a DRIFT.
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                                     not used here !!!
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xhm,yhm,zhm   ==>  massive part position in middle of step
c                                       (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adopted from martin's nbwh.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/12/93

      subroutine step_dkd_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xhm,yhm,zhm,xht,yht,zht,vxht,vyht,vzht,istat,dt)	

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 xhm(nbod),yhm(nbod),zhm(nbod)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
      real*8 dth 
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  Take a drift forward 0.5dt to begin the step (the first kick = 0).
      call drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dth,istat)	

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xhm,yhm,zhm,
     &                xht,yht,zht,istat,axht,ayht,azht)

c...  Apply a heliocentric kick for a full dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dt) 

c..   Drift again in Jacobi coords for the final half-step 0.5dt	  
      call drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dth,istat)	

      return
      end   ! step_dkd_tp
c---------------------------------------------------------------------

