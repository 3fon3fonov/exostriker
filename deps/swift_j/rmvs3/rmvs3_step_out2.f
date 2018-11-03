c*************************************************************************
c                            RMVS3_STEP_OUT2.F
c*************************************************************************
c This subroutine takes a full dt step in helio coord for test particles
c in the outer region of an encounter.  It also sets up and call
c the inner region inetgration if necessary
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c             j2rp2,j4rp4        ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c             xbeg,ybeg,zbeg     ==>  initial planet position in helio coord 
c                                         (real arrays)
c             vxbeg,vybeg,vzbeg  ==>  initial planet velcoity in helio coord 
c                                         (real arrays)
c             xend,yend,zend     ==>  final planet position in helio coord 
c                                         (real arrays)
c             vxend,vyend,vzend  ==>  final planet velcoity in helio coord 
c                                         (real arrays)
c             xht,yht,zht        ==>  initial tp position in helio coord 
c                                      (real arrays)
c             vxht,vyht,vzht     ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c             istat              ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c             ienc0               ==> ienc0(j) = 0 if tp j not involved in enc 
c                                       in outer region: = planet# if it is. 
c                                         (integer array)
c                                              =  0 No (integer scalar)  
c                 dt            ==>  time step (real sclar)
c             Output:
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c                                      NOTE: only the tp in the outer region
c                                            will have their x and v's changed 
c             istat              ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c             ienc0              ==> Multiplied by -1 if particle entered inner
c                                     region (integer array)
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real array)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    7/10/96
c Last revision: 

      subroutine rmvs3_step_out2(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &           xbeg,ybeg,zbeg,xend,yend,zend,vxbeg,vybeg,
     &           vzbeg,vxend,vyend,vzend,xht,
     &     yht,zht,vxht,vyht,vzht,istat,dt,ienc0,isperi,peri)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      integer ienc0(NTPMAX)

c...  Outputs:
      real*8 peri(NTPMAX)
      integer isperi(NTPMAX)

c...  Internals
      integer i1sttp,i,j,istattmp(NTPMAX,NSTAT),icflg,i1sto
      real*8 rts
      Integer nenc(NPLMAX),itpenc(NTPMAX,NPLMAX),ienc(NTPMAX)
      real*8 xtmp(NPLMAX,NTPHENC),ytmp(NPLMAX,NTPHENC)
      real*8 ztmp(NPLMAX,NTPHENC)
      real*8 vxtmp(NPLMAX,NTPHENC),vytmp(NPLMAX,NTPHENC)
      real*8 vztmp(NPLMAX,NTPHENC)


c----
c...  Executable code 

      i1sttp = i1st

c...  are there any encounters?
      rts = RHPSCALE*RHPSCALE
      call rmvs3_chk(nbod,ntp,mass,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xht,yht,zht,vxht,vyht,vzht,istat,dt,rts,icflg,nenc,
     &     itpenc,ienc)

c...  keep track of encounters in the inner region
      call rmvs3_elog(ntp,icflg,ienc,istat)

c.... if not just do a normal step and leave
      if(icflg.eq.0) then

          call step_kdk_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istat,dt)	

         return       !  NOTE AN EXIT
      endif

c...  INNER ENCOUNTER STUFF FROM HERE ON!!!!!

c...  Now do the interpolation for intermediate steps
      call rmvs3_interp(nbod,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,dt,mass(1),NTPHENC,
     &     xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp)

c...  Do the inner integration
      call rmvs3_step_in(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,xtmp,
     &     ytmp,ztmp,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,xht,yht,zht,
     &     vxht,vyht,vzht,istat,nenc,itpenc,isperi,peri,dt)
      
      do i=1,ntp
         if(istat(i,1).eq.0) then
            istat(i,2) = 0
         endif
      enddo

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
            if(ienc(i).ne.0) then
               ienc0(i) = -1 * abs(ienc0(i))
            endif
         endif
      enddo

c...  we MUST make sure that the saved accel arrays are ok
c...  calculate them again
       i1st = 0 

      return
      end   ! step_enc
c------------------------------------------------------------------------

