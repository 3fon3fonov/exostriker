c*************************************************************************
c                            HELIO_STEP_TP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xbeg,ybeg,zbeg ==>  massive part position at beginning of dt
c                                       (real arrays)
c                 xend,yend,zend ==>  massive part position at end of dt
c                                       (real arrays)
c                 ptxb,ptyb,ptzb ==> beg momentum of sun
c                                       (real scalars)
c                 ptxe,ptye,ptze ==> end momentum of sun
c                                       (real scalars)
c                vxsb,vxsb,vxsb  ==> Initial vel of the Sun: tp's need this
c                                         (real scalars)
c                vxse,vxse,vxse  ==> final vel of the Sun: tp's need this
c                                       (real scalars)
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
c Remarks: Based on step_kdk_tp
c Authors:  Hal Levison 
c Date:    11/14/96
c Last revision: 

      subroutine helio_step_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,xend,yend,zend,ptxb,ptyb,ptzb,ptxe,ptye,
     &     ptze,vxsb,vysb,vzsb,vxse,vyse,vzse,
     &     xht,yht,zht,vxht,vyht,vzht,istat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 ptxb,ptyb,ptzb
      real*8 ptxe,ptye,ptze
      real*8 vxsb,vysb,vzsb,vxse,vyse,vzse

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
c      integer i
      real*8 dth 
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)

      save axht,ayht,azht     ! Note this !!
      save vxbt,vybt,vzbt

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then 
c...     Convert velocities to bary
         call coord_vh2b_tp(ntp,vxht,vyht,vzht,vxsb,vysb,vzsb,
     &        vxbt,vybt,vzbt)
         i1st = 1               ! turn this off
      endif

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift_tp(ntp,ptxb,ptyb,ptzb,dth,
     &     istat,xht,yht,zht)

c...  Get the accelerations in helio frame.
      call helio_getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,xht,yht,zht,istat,axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxbt,vybt,vzbt,axht,ayht,azht,istat,dth) 

c...  Take a drift forward full step
      call drift_tp(ntp,mass(1),xht,yht,zht,vxbt,vybt,vzbt,dt,istat)	

c...  Get the accelerations in helio frame.
      call helio_getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,
     &     xend,yend,zend,xht,yht,zht,istat,axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxbt,vybt,vzbt,axht,ayht,azht,istat,dth) 

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift_tp(ntp,ptxe,ptye,ptze,dth,
     &     istat,xht,yht,zht)

c...   Put back to helio
      call coord_vb2h_tp(ntp,istat,vxbt,vybt,vzbt,vxse,vyse,vzse,
     &     vxht,vyht,vzht)

      return
      end   ! helio_step_tp
c---------------------------------------------------------------------

