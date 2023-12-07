c*************************************************************************
c                            MVITS_STEP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  Uses Individual timesteps
c for tp and planets.  Note that the planets are all integrated together
c with the smallest time step.  The test particles all have the 
c timestep dt.
c
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
c                 dt            ==>  time step
c                 rstat         ==>  status of the test paricles
c                                      (2d real array)
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
c Authors:  Hal Levison 
c Date:    12/16/93
c Last revision: 2/24/94

      subroutine mvits_step(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt)

      include '../swift.inc'
      include 'mvits.inc'

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
      integer i1stloc,i,icall
      real*8 xstor(NPLMAX,0:NDTMAX),ystor(NPLMAX,0:NDTMAX)
      real*8 massi(NPLMAX,0:NDTMAX)
      real*8 zstor(NPLMAX,0:NDTMAX),dtpli(NPLMAX),dtpli2(NPLMAX)
      real*8 cosrot(NTPMAX,2),sinrot(NTPMAX,2)
      real*8 dtpl
      integer idpl(NPLMAX,0:NDTMAX),nbodi(0:NDTMAX),nstep,istep,i1stpl
      character*80 indtfile

      real*8 axht(NTPMAX,NPLMAX),ayht(NTPMAX,NPLMAX),azht(NTPMAX,NPLMAX)

      data i1stloc/0/

      save i1stloc,cosrot,sinrot,dtpl,nstep,dtpli,dtpli2
      save icall,massi,nbodi,idpl
      save axht,ayht,azht     ! Note this !!

c----
c...  Executable code 

c...  If this is the first time through read time ratios. also calculate
c...  rotation sin and cos
      if(i1stloc.eq.0) then
        i1stloc = 1
        icall = 1
        write(*,*) 'Enter name of timestep data file : '
        read(*,999) indtfile
999 	format(a)
	call io_mvits_init(indtfile,nbod,dt,mass,
     &       dtpl,nstep,dtpli,massi,nbodi,idpl)
        do i=1,nbod
           dtpli2(i) = 2.*dtpli(i)
        enddo
	call mvits_init_rot(ntp,dtpl,nstep,mass(1),xht,yht,zht,
     &       vxht,vyht,vzht,cosrot,sinrot)
        write(*,*) ' Planet integeration will take ',
     &              nstep,' steps per dt '
        write(*,*) ' For test particle integration: '
        write(*,*) '  Planet      dt '
        do i=2,nbod
           write(*,*) i, '   ',dtpli(i)
        enddo
        write(*,*) ' '
        write(*,*) ' MVITS version 2 '
        write(*,*) ' CONTINUE: '
      endif
      
c...  Now do the planets
      do i=1,nbod
         xstor(i,0) = xh(i)
         ystor(i,0) = yh(i)
         zstor(i,0) = zh(i)
      enddo

      i1stpl = i1st
      do istep=1,nstep
         call step_kdk_pl(i1stpl,nbod,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,dtpl)
         do i=1,nbodi(istep)
            xstor(i,istep) = xh(idpl(i,istep))
            ystor(i,istep) = yh(idpl(i,istep))
            zstor(i,istep) = zh(idpl(i,istep))
         enddo
      enddo

      if(icall.gt.0) then
         call mvits_lstep_tp(i1st,ntp,xstor,ystor,zstor,cosrot,
     &        sinrot,dtpl,nstep,dtpli,dtpli2,massi,j2rp2,j4rp4,nbodi,
     &        idpl,xht,yht,zht,vxht,vyht,vzht,istat,dt,axht,ayht,azht)
      else
         call mvits_rstep_tp(i1st,ntp,xstor,ystor,zstor,cosrot,
     &        sinrot,dtpl,nstep,dtpli,dtpli2,massi,j2rp2,j4rp4,nbodi,
     &        idpl,xht,yht,zht,vxht,vyht,vzht,istat,dt,axht,ayht,azht)
      endif

      icall = -1*icall

      return
      end   ! mvits_step
c------------------------------------------------------------------------

