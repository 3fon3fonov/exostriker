c*************************************************************************
c                            RMVS2_STEP_OUT.F
c*************************************************************************
c This subroutine takes a full dt step in helio coord for test particles
c in the outer region of an encounter.  It will also remember
c where the planets are for the planocentric integration if necessary.  
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
c             xtmpo,ytmpo,ztmpo  ==>  position of planet wrt time
c                                       (2d real arrays)
c             xht,yht,zht        ==>  initial tp position in helio coord 
c                                      (real arrays)
c             vxht,vyht,vzht     ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c             istat              ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c             ienco              ==> ienco(j) = 0 if tp j not involved in enc 
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
c
c
c Remarks: Adopted from martin's nbwh.f program
c Authors:  Hal Levison 
c Date:    8/25/94
c Last revision: 

      subroutine rmvs2_step_out(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,xtmpo,ytmpo,ztmpo,xht,yht,zht,
     &     vxht,vyht,vzht,istat,ienco,dt)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 xtmpo(NPLMAX,NTPENC),ytmpo(NPLMAX,NTPENC)
      real*8 ztmpo(NPLMAX,NTPENC)
      integer ienco(NTPMAX)
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer i1sttp,i,j,ic,istattmp(NTPMAX,NSTAT)
      real*8 dto
      real*8 xbegi(NPLMAX),ybegi(NPLMAX),zbegi(NPLMAX)
      real*8 xendi(NPLMAX),yendi(NPLMAX),zendi(NPLMAX)

c----
c...  Executable code 

      i1sttp = i1st

      dto = dt/float(NTENC)

c...  We only want to integrate the position of test particles in outer region
c...  make a temporary istat array ao only they are active
      do i=1,ntp
         if(ienco(i).eq.0) then
            istattmp(i,1) = 1       ! don't integrate it
          else
            istattmp(i,1) = 0       ! integrate it
          endif
          do j=2,NSTAT
             istattmp(i,j) = 0
          enddo
       enddo

c...  do integration of outer loop
       ic = 0
       do i=1,NTENC

c...      remember the current position of the planets
          if(i.eq.1) then
             do j=1,nbod
                xbegi(j) = xbeg(j)
                ybegi(j) = ybeg(j)
                zbegi(j) = zbeg(j)
             enddo
          else
             do j=1,nbod
                xbegi(j) = xtmpo(j,i-1)
                ybegi(j) = ytmpo(j,i-1)
                zbegi(j) = ztmpo(j,i-1)
             enddo
          endif

          do j=1,nbod
             xendi(j) = xtmpo(j,i)
             yendi(j) = ytmpo(j,i)
             zendi(j) = ztmpo(j,i)
          enddo

          call step_kdk_tp(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbegi,ybegi,zbegi,xendi,yendi,zendi,
     &              xht,yht,zht,vxht,vyht,vzht,istattmp,dto)	

       enddo

c...   Have to update istat just in case damby had problems
      do i=1,ntp
         if( (ienco(i).ne.0) .and. (istattmp(i,1).ne.0) ) then
            istat(i,1) = 1       ! it had problems
            do j=2,NSTAT
                 istat(i,j) = istattmp(i,j)
            enddo
         endif
       enddo
      
       return
       end     ! rmvs2_step_out
c----------------------------------------------------------------
