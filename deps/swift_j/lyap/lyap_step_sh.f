c*************************************************************************
c                            LYAP_STEP_SH.F
c*************************************************************************
c Just like STEP_KDK_TP except we use shadow particles.
c CAN'T JUST USE STEP_KDK_TP BECAUSE WE SAVE ACCELS IN THE SUBROUTINE
c BETWEEN CALLS IF I1ST + 1.
c To avoid errors, tho, I will call the coords xht,...vzht
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
c Remarks: Adopted from martin's nbwhnew.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/25/94

      subroutine lyap_step_sh(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istat,dt)	

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

c...  Internals:
c      integer i
      real*8 dth 
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)

      save axht,ayht,azht     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xbeg,ybeg,zbeg,
     &                  xht,yht,zht,istat,axht,ayht,azht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

c...  Take a drift forward full step
      call drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dt,istat)	

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xend,yend,zend,
     &                  xht,yht,zht,istat,axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

      return
      end   ! lyap_step_sh
c---------------------------------------------------------------------

