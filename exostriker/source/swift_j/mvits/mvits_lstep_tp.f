c*************************************************************************
c                            MVITS_LSTEP_TP.F
c*************************************************************************
c This subroutine takes a step in helio coord using independant time steps
c ONLY DOES TEST PARTICLES
c
c  LEFT HAND SIDE
c
c             Input:
c                 i1st               ==>  = 0 if first step; = 1 not 
c                                             (int scalar)
c                 ntp                ==>  number of test bodies (int scalar)
c                 xstor,ystor,zstor  ==>  massive part position wrt time
c                                          (2d real arrays)
c                 cosrot,sinrot      ==>  cos and sin of rotation angle of tp
c                                          (real arrays)
c                 dtpl               ==> timestep for planet integration.
c                                             (real scalar)
c                 nstep              ==> largest of previous (int scalar)
c                 dtpli              ==> individual timestep for planets.
c                                             (real aray)
c                 dtpli2             ==> 2 X above
c                                             (real aray)
c                 massi              ==>  mass of bodies wrt t (2d real array)
c                 j2rp2,j4rp4        ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                          (real scalars)
c                 nbodi              ==>  number of massive bodies wrt t
c                                           (int array)
c                 idpl               ==> list link of planet id wrt t 
c                                          (2d int array)
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
c                 axht,ayht,azht ==>  accel of individual planets on tp 
c                                       (2d real arrays)
c
c Remarks: Adopted from step_kdk
c Authors:  Hal Levison 
c Date:    12/16/93
c Last revision: 2/24/94

      subroutine mvits_lstep_tp(i1st,ntp,xstor,ystor,zstor,cosrot,
     &     sinrot,dtpl,nstep,dtpli,dtpli2,massi,j2rp2,j4rp4,nbodi,
     &     idpl,xht,yht,zht,vxht,vyht,vzht,istat,dt,axht,ayht,azht)

      include '../swift.inc'
      include 'mvits.inc'

c...  Inputs Only: 
      integer ntp,i1st
      integer idpl(NPLMAX,0:NDTMAX)
      real*8 massi(NPLMAX,0:NDTMAX),dt,j2rp2,j4rp4
      real*8 xstor(NPLMAX,0:NDTMAX),ystor(NPLMAX,0:NDTMAX)
      real*8 zstor(NPLMAX,0:NDTMAX),dtpli(NPLMAX),dtpl,dtpli2(NPLMAX)
      real*8 cosrot(NTPMAX,2),sinrot(NTPMAX,2)
      integer nbodi(0:NDTMAX),nstep

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 axht(NTPMAX,NPLMAX),ayht(NTPMAX,NPLMAX),azht(NTPMAX,NPLMAX)

c...  Internals:
      integer i

c----
c...  Executable code 

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call mvits_getacch_tp(nbodi(0),ntp,massi(1,0),j2rp2,j4rp4,
     &        xstor(1,0),ystor(1,0),zstor(1,0),xht,yht,zht,istat,
     &        axht,ayht,azht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick 
      call mvits_kickvh_tp(nbodi(0),ntp,vxht,vyht,vzht,axht,ayht,
     &     azht,istat,dtpli,idpl(1,0)) 

c...  rotate forwards
      do i=1,nstep-1
         call mvits_rotate(1.0d0,ntp,cosrot(1,1),sinrot(1,1),xht,yht,
     &        vxht,vyht)
         call mvits_getacch_tp(nbodi(i),ntp,massi(1,i),j2rp2,j4rp4,
     &        xstor(1,i),ystor(1,i),zstor(1,i),xht,yht,zht,istat,
     &        axht,ayht,azht)
         call mvits_kickvh_tp(nbodi(i),ntp,vxht,vyht,vzht,axht,ayht,
     &        azht,istat,dtpli2,idpl(1,i)) 
      enddo

c...  rotate back the big angle
      if(nstep.ne.1) then
         call mvits_rotate(-1.0d0,ntp,cosrot(1,2),sinrot(1,2),xht,yht,
     &     vxht,vyht)
      endif

c...  Take a drift forward full step
      call drift_tp(ntp,massi(1,1),xht,yht,zht,vxht,vyht,vzht,dt,istat)	

c...  now do outside
      i = nstep
      call mvits_getacch_tp(nbodi(0),ntp,massi(1,i),j2rp2,j4rp4,
     &     xstor(1,i),ystor(1,i),zstor(1,i),xht,yht,zht,istat,
     &     axht,ayht,azht)
      call mvits_kickvh_tp(nbodi(0),ntp,vxht,vyht,vzht,axht,ayht,
     &     azht,istat,dtpli,idpl(1,0)) 

      return
      end   ! mvits_lstep_tp
c---------------------------------------------------------------------


