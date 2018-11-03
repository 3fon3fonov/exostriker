c*************************************************************************
c                            RMVS3_STEP_OUT.F
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
c             xtmp,ytmp,ztmp     ==>  position of planet wrt time
c                                       (2d real arrays)
c             vxtmp,vytmp,vztmp   ==>  velocity of planet wrt time
c                                       (2d real arrays)
c             xht,yht,zht        ==>  initial tp position in helio coord 
c                                      (real arrays)
c             vxht,vyht,vzht     ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c             istat              ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c             ienc               ==> ienc(j) = 0 if tp j not involved in enc 
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
c             ienc               ==> Multiplied by -1 if particle entered inner
c                                     region (integer array)
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real array)
c
c Remarks: Adopted from rmvs2_step_out.f
c Authors:  Hal Levison 
c Date:    7/10/96
c Last revision: 

      subroutine rmvs3_step_out(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,xtmp,ytmp,ztmp,vxtmp,
     &     vytmp,vztmp,xht,yht,zht,vxht,vyht,vzht,istat,ienc,dt,
     &     isperi,peri)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 xtmp(NPLMAX,NTPENC),ytmp(NPLMAX,NTPENC)
      real*8 ztmp(NPLMAX,NTPENC)
      real*8 vxtmp(NPLMAX,NTPENC),vytmp(NPLMAX,NTPENC)
      real*8 vztmp(NPLMAX,NTPENC)
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      integer ienc(NTPMAX)

c...  Outputs:
      real*8 peri(NTPMAX)
      integer isperi(NTPMAX)

c...  Internals
      integer i1sttp,i,j,ic,istattmp(NTPMAX,NSTAT)
      real*8 dto
      real*8 xbegi(NPLMAX),ybegi(NPLMAX),zbegi(NPLMAX)
      real*8 xendi(NPLMAX),yendi(NPLMAX),zendi(NPLMAX)
      real*8 vxbegi(NPLMAX),vybegi(NPLMAX),vzbegi(NPLMAX)
      real*8 vxendi(NPLMAX),vyendi(NPLMAX),vzendi(NPLMAX)
      real*8 peril(NTPMAX)
      integer isperil(NTPMAX)

c----
c...  Executable code 

      i1sttp = i1st

      dto = dt/float(NTENC)

c...  We only want to integrate the position of test particles in outer region
c...  make a temporary istat array ao only they are active
      do i=1,ntp
         if(ienc(i).eq.0) then
            istattmp(i,1) = 1       ! don't integrate it
          else
            istattmp(i,1) = 0       ! integrate it
          endif
          do j=2,NSTAT
             istattmp(i,j) = istat(i,j)
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
                vxbegi(j) = vxbeg(j)
                vybegi(j) = vybeg(j)
                vzbegi(j) = vzbeg(j)
             enddo
          else
             do j=1,nbod
                xbegi(j) = xtmp(j,i-1)
                ybegi(j) = ytmp(j,i-1)
                zbegi(j) = ztmp(j,i-1)
                vxbegi(j) = vxtmp(j,i-1)
                vybegi(j) = vytmp(j,i-1)
                vzbegi(j) = vztmp(j,i-1)
             enddo
          endif

          do j=1,nbod
             xendi(j) = xtmp(j,i)
             yendi(j) = ytmp(j,i)
             zendi(j) = ztmp(j,i)
             vxendi(j) = vxtmp(j,i)
             vyendi(j) = vytmp(j,i)
             vzendi(j) = vztmp(j,i)
          enddo

          call rmvs3_step_out2(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &         xbegi,ybegi,zbegi,xendi,yendi,zendi,vxbegi,vybegi,
     &         vzbegi,vxendi,vyendi,vzendi,xht,yht,zht,
     &         vxht,vyht,vzht,istattmp,dto,ienc,isperil,peril)	

          do j=1,ntp
             if(isperil(j).eq.0) then
                peri(j) = peril(j)
                isperi(j) = isperil(j)
             endif
             if( (i.eq.1) .or. (isperi(j).ne.0) ) then 
                isperi(j) = isperil(j)
             endif
          enddo

       enddo


c...   Have to update istat just in case damby had problems
      do i=1,ntp
         if( (ienc(i).ne.0) .and. (istattmp(i,1).ne.0) ) then
            istat(i,1) = 1       ! it had problems
            do j=2,NSTAT
                 istat(i,j) = istattmp(i,j)
            enddo
         endif
         if(ienc(i).lt.0) then   ! Now update the encounter stuff
            do j=4,NSTAT
                 istat(i,j) = istattmp(i,j)
              enddo
           endif
       enddo

       return
       end     ! rmvs3_step_out
c----------------------------------------------------------------
