c*************************************************************************
c                            RMVS_STEP_IN_TP.F
c*************************************************************************
c This subroutine takes a step for tp in inner region.
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
c                 ipl            ==>  the planet that is currently in the 
c                                      center (integer scalar)
c           aoblxb,aoblyb,aoblzb ==> acceleration of the Sun on the central pl
c                                    at beginning of dt due to J2 and J4
c                                      (real scalars)
c           aoblxe,aoblye,aoblze ==> acceleration of the Sun on the central pl
c                                    at end of dt  due to J2 and J4
c                                         (real scalars)
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Taken from step_kdk_tp.f
c Authors:  Hal Levison 
c Date:    2/24/94
c Last revision: 

      subroutine rmvs_step_in_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istat,dt,ipl,
     &              aoblxb,aoblyb,aoblzb,aoblxe,aoblye,aoblze)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st,ipl
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 aoblxb,aoblyb,aoblzb,aoblxe,aoblye,aoblze

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
      real*8 dth,j2rp2i,j4rp4i
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)

      save axht,ayht,azht     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  disable J2 and J4 terms in the getacch_tp calls
      j2rp2i = 0.0d0
      j4rp4i = 0.0d0

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_tp(nbod,ntp,mass,j2rp2i,j4rp4i,xbeg,ybeg,zbeg,
     &                  xht,yht,zht,istat,axht,ayht,azht)
          call rmvs_obl_acc(nbod,ntp,ipl,mass,j2rp2,j4rp4,xbeg,
     &         ybeg,zbeg,aoblxb,aoblyb,aoblzb,xht,yht,zht,istat,
     &         axht,ayht,azht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

c...  Take a drift forward full step
      call drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dt,istat)	

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,ntp,mass,j2rp2i,j4rp4i,xend,yend,zend,
     &                  xht,yht,zht,istat,axht,ayht,azht)
      call rmvs_obl_acc(nbod,ntp,ipl,mass,j2rp2,j4rp4,xend,yend,
     &     zend,aoblxe,aoblye,aoblze,xht,yht,zht,istat,
     &     axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

      return
      end   ! rmvs_step_in_tp
c---------------------------------------------------------------------

