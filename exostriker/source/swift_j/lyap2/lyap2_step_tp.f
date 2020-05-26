c*************************************************************************
c                            LYAP2_STEP_TP
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
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c              dxht,dyht,dzht    ==>  initial separation in position
c                                      (real arrays)
c              dvxht,dvyht,dvzht ==>  initial separation in velocity
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
c              dxht,dyht,dzht    ==>  final separation in position
c                                      (real arrays)
c              dvxht,dvyht,dvzht ==>  final separation in velocity
c                                        (real arrays)
c
c Remarks: Adopted from step_kdk_tp.f
c Authors:  Hal Levison 
c Date:    7/11/95
c Last revision: 

      subroutine lyap2_step_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,
     &              dxht,dyht,dzht,dvxht,dvyht,dvzht,istat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 dxht(ntp),dyht(ntp),dzht(ntp)
      real*8 dvxht(ntp),dvyht(ntp),dvzht(ntp)

c...  Internals:
c      integer i
      real*8 dth 
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)
      real*8 daxht(NTPMAX),dayht(NTPMAX),dazht(NTPMAX)

      save axht,ayht,azht     ! Note this !!
      save daxht,dayht,dazht     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xbeg,ybeg,zbeg,
     &        xht,yht,zht,istat,axht,ayht,azht)
          call lyap2_acc_tp(nbod,ntp,mass,j2rp2,j4rp4,xbeg,ybeg,zbeg,
     &         xht,yht,zht,dxht,dyht,dzht,istat,daxht,dayht,dazht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 
      call kickvh_tp(ntp,dvxht,dvyht,dvzht,daxht,dayht,dazht,istat,dth) 

c...  Take a drift forward full step
      call lyap2_drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dxht,
     &     dyht,dzht,dvxht,dvyht,dvzht,dt,istat)

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xend,yend,zend,
     &                  xht,yht,zht,istat,axht,ayht,azht)
      call lyap2_acc_tp(nbod,ntp,mass,j2rp2,j4rp4,xend,yend,zend,
     &     xht,yht,zht,dxht,dyht,dzht,istat,daxht,dayht,dazht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 
      call kickvh_tp(ntp,dvxht,dvyht,dvzht,daxht,dayht,dazht,istat,dth) 

      return
      end   ! lyap2_step_tp
c---------------------------------------------------------------------

