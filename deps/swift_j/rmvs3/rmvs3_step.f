c*************************************************************************
c                            RMVS3_STEP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c both massive and test particles
c INCLUDES close approuches between test particles and planets
c
c VERSION 3 of RMVS
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial planet position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial planet velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial tp position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final planet position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final planet velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c
c
c Remarks: Based on rmvs2_step.f
c Authors:  Hal Levison 
c Date:    7/10/96
c Last revision: 

      subroutine rmvs3_step(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt)	

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

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
      integer i1sttp,i1stpl,i1sto,icflg,i,j,np
      real*8 xtmp(NPLMAX,NTENC),ytmp(NPLMAX,NTENC)
      real*8 ztmp(NPLMAX,NTENC),rts,peri(NTPMAX)
      real*8 vxtmp(NPLMAX,NTENC),vytmp(NPLMAX,NTENC)
      real*8 vztmp(NPLMAX,NTENC)
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)
      Integer nenc(NPLMAX),itpenc(NTPMAX,NPLMAX),ienc(NTPMAX)
      integer istattmp(NTPMAX,NSTAT),isperi(NTPMAX)

c----
c...  Executable code 

      i1sttp = i1st
      i1sto = i1st

c...  are there any encounters?
      rts = RHSCALE*RHSCALE
      call rmvs3_chk(nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,
     &       zht,vxht,vyht,vzht,istat,dt,rts,icflg,nenc,itpenc,ienc)
c...     nenc and itpenc not used here!

c.... if not just do a normal step and leave
      if(icflg.eq.0) then
         call step_kdk(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,xh,yh,zh,
     &        vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,dt)	

         do i=1,ntp
            if(istat(i,1).eq.0) then
               istat(i,2) = 0
            endif
         enddo

         return       !  NOTE AN EXIT
      endif

c...  ENCOUNTER STUFF FROM HERE ON!!!!!

c...  save initial x and v of planets if there are planocentric enc
      do i=1,nbod
            xbeg(i) = xh(i)
            ybeg(i) = yh(i)
            zbeg(i) = zh(i)
            vxbeg(i) = vxh(i)
            vybeg(i) = vyh(i)
            vzbeg(i) = vzh(i)
      enddo

c... do a full step for the planets
      i1stpl = i1st
      call step_kdk_pl(i1stpl,nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,dt)	

c...  save the final position and velocity of planets
      do i=1,nbod
         xend(i) = xh(i)
         yend(i) = yh(i)
         zend(i) = zh(i)
         vxend(i) = vxh(i)
         vyend(i) = vyh(i)
         vzend(i) = vzh(i)
      enddo

c...  Now do the interpolation for intermediate steps
      call rmvs3_interp(nbod,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,dt,mass(1),NTENC,
     &     xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp)

c...  do the integration
      call rmvs3_step_out(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,xtmp,ytmp,ztmp,vxtmp,
     &     vytmp,vztmp,xht,yht,zht,vxht,vyht,vzht,istat,ienc,dt,
     &     isperi,peri)


c...  As of this point all the test particles that are involved in an
c...  encounter have been moved.  But not the ones that have not.
c...  so move those,  BUT NOT the onces in the encounter


c...  make a temporary istat array so only they are active
      do i=1,ntp
         if(istat(i,1).eq.0) then
            if(ienc(i).ne.0) then
               istattmp(i,1) = 1 ! don't integrate it
            else
               istattmp(i,1) = 0 ! integrate it
            endif
            do j=2,NSTAT
               istattmp(i,j) = 0
            enddo
         else
            istattmp(i,1) = 1   ! don't integrate it
         endif
       enddo

c...  do a full step
      i1sto = 0      ! we need to recalculate accel arrays
      call step_kdk_tp(i1sto,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istattmp,dt)	

c...  fix up the istat array
      do i=1,ntp
         if(istattmp(i,2) .ne. 0) then   ! danby screwed up
            istat(i,1) = 1
            do j=2,NSTAT
               istat(i,j) = istattmp(i,j)
            enddo
          endif
       enddo

c...  put the enc info into istat
      do i=1,ntp
         if(istat(i,1).eq.0) then
            if(ienc(i).lt.0) then
               istat(i,2) = ienc(i)
               istat(i,3) = abs(ienc(i))
               if(isperi(i).eq.0) then
                  rstat(i,3) = peri(i)
                  np = NSTATP + abs(ienc(i)) - 1
                  if(rstat(i,np).eq.0) then
                     rstat(i,np) = peri(i)
                  else
                     rstat(i,np) = min(rstat(i,np),peri(i))
                  endif
               else
                  rstat(i,3) = 0.0d0
               endif
            else if(ienc(i).gt.0) then
               istat(i,2) = ienc(i)
            else
               istat(i,2) = 0
            endif
         endif
      enddo

c...  we MUST make sure that the saved accel arrays are ok
c...  calculate them again
       i1st = 0 

      return
      end   ! step_enc
c------------------------------------------------------------------------






