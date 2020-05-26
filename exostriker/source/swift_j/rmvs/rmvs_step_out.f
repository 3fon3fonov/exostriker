c*************************************************************************
c                            RMVS_STEP_OUT.F
c*************************************************************************
c This subroutine takes a full dt step in helio coord for particles
c in the outer region of an encounter.  It will also remember
c where the planets are for the planocentric integration if necessary.  
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
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
c                 icflg         ==> ecounters? = 1 Yes, in outer region only
c                                              = -1 in inner region
c                 ienco         ==> ienco(j) = 0 if tp j not involved in enc 
c                                   in outer region: = planet# if it is. 
c                                     (integer array)
c                                              =  0 No (integer scalar)  
c                 dt            ==>  time step (real sclar)
c             Output:
c                 xh,yh,zh      ==>  final planet position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final planet velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c                                      NOTE: only the tp in the outer region
c                                            will have their x and v's changed 
c                 xtmp,ytmp,ztmp  ==>  position of planet wrt time
c                                       (2d real arrays)
c
c
c Remarks: Adopted from martin's nbwh.f program
c Authors:  Hal Levison 
c Date:    2/19/93
c Last revision: 2/24/94

      subroutine rmvs_step_out(i1st,nbod,ntp,mass,j2rp2,j4rp4,xh,yh,zh,
     &              vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &              icflg,ienco,dt,xtmp,ytmp,ztmp)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      integer icflg,ienco(NTPMAX)

c...  Outputs only
      real*8 xtmp(NPLMAX,NTPENC),ytmp(NPLMAX,NTPENC)
      real*8 ztmp(NPLMAX,NTPENC)

c...  Internals
      integer i1sttp,nstep,i,j,k,ic,istattmp(NTPMAX,NSTAT)
      real*8 dto,dti
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)

c----
c...  Executable code 

      i1sttp = i1st

      dto = dt/float(NTENC)
      if(icflg.le.-1) then    ! need high t resolution of planocentric
        nstep = NTPHENC
        dti = dt/float(NTPENC)
      else
        nstep = 1
        dti = dt/float(NTENC)
      endif

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
          do j=1,nbod
             xbeg(j) = xh(j)
             ybeg(j) = yh(j)
             zbeg(j) = zh(j)
          enddo

          do j=1,nstep
            call step_kdk_pl(i1st,nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &            vxh,vyh,vzh,dti)	
            if(icflg.le.-1) then 
               ic = ic + 1
               do k=1,nbod
                  xtmp(k,ic) = xh(k)
                  ytmp(k,ic) = yh(k)
                  ztmp(k,ic) = zh(k)
               enddo
            endif
          enddo

c...      remember these positions
          do j=1,nbod
             xend(j) = xh(j)
             yend(j) = yh(j)
             zend(j) = zh(j)
          enddo

          call step_kdk_tp(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
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
       end     ! rmvs_step_out
c----------------------------------------------------------------
