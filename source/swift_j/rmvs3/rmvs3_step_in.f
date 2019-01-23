c*************************************************************************
c                            RMVS3_STEP_IN.F
c*************************************************************************
c This subroutine takes a full dt step in helio coord for test particles
c in the inner region of an encounter.
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xpl,ypl,zpl    ==>  position of planet wrt time
c                                       (2d real arrays)
c                 xbeg,ybeg,zbeg ==>  position of planet @ beginning of dt
c                                       (real arrays)
c                 vxbeg,vybeg,vzbeg ==>  vel of planet @ beginning of dt
c                                       (2d real arrays)
c                 xend,yend,zend ==>  position of planet @ end of dt
c                                       (2d real arrays)
c                 vxend,vyend,vzend ==>  vel of planet @ end of dt
c                                       (2d real arrays)
c                 xht,yht,zht    ==>  initial tp position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c                 istat          ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 nenci          ==> nenci(i) is the number of tp enc planet i
c                                    in inner region (integer array)
c                 itpenci       ==> itpenci(*,i) is a list of tp enc planet i
c                                   in inner region (2d integer array)
c                 dt            ==>  time step (real sclar)
c             Output:
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c                                      NOTE: only the tp in the inner region
c                                            will have their x and v's changed 
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real array)
c
c
c Remarks: Same as rmvs_step_in.f except for the timesteps
c Authors:  Hal Levison 
c Date:    7/10/96
c Last revision: 10/3/96

      subroutine rmvs3_step_in(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &           xpl,ypl,zpl,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &           xend,yend,zend,vxend,vyend,vzend,xht,yht,zht,
     &           vxht,vyht,vzht,istat,nenci,itpenci,isperi,peri,dt)


      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 xpl(NPLMAX,NTPHENC),ypl(NPLMAX,NTPHENC)
      real*8 zpl(NPLMAX,NTPHENC)
      Integer nenci(NPLMAX),itpenci(NTPMAX,NPLMAX)
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Outputs:
      real*8 peri(NTPMAX)
      integer isperi(NTPMAX)

c...  Internals
      integer i1sttp,i,j,k,link,istattmp(NTPMAX,NSTAT)
      real*8 xtpt(NTPMAX),ytpt(NTPMAX),ztpt(NTPMAX),dti
      real*8 vxtpt(NTPMAX),vytpt(NTPMAX),vztpt(NTPMAX)
      real*8 xpltb(NPLMAX),ypltb(NPLMAX),zpltb(NPLMAX)
      real*8 xplte(NPLMAX),yplte(NPLMAX),zplte(NPLMAX)
      real*8 masst(NPLMAX),perit(NTPMAX)
      integer isperit(NTPMAX)
      logical*2 lperit(NTPMAX)

      real*8 irh(NPLMAX),aoblx(NPLMAX,0:NTPHENC)
      real*8 aobly(NPLMAX,0:NTPHENC),aoblz(NPLMAX,0:NTPHENC)

c----
c...  Executable code 

      dti = dt/float(NTPHENC)

c...  First get the accel due to J2 and J4 in barycentric frame on planets
      if(j2rp2.ne.0.0d0) then
         do i=2,nbod
            irh(i) = 1.0d0/sqrt( xbeg(i)**2 + ybeg(i)**2 + zbeg(i)**2 )
         enddo
         call obl_acc(nbod,mass,j2rp2,j4rp4,xbeg,ybeg,zbeg,
     &        irh,aoblx(1,0),aobly(1,0),aoblz(1,0))
         do j=1,NTPHENC
            do i=2,nbod
               irh(i) = 1.0d0/sqrt( xpl(i,j)**2 + ypl(i,j)**2 + 
     &              zpl(i,j)**2 )
            enddo
            call obl_acc(nbod,mass,j2rp2,j4rp4,xpl(1,j),ypl(1,j),
     &        zpl(1,j),irh,aoblx(1,j),aobly(1,j),aoblz(1,j))
         enddo
      endif

c...  Must do each planet one at a time
        do i=2,nbod
           if(nenci(i).ne.0) then

c...          set up test particles
              do j = 1,nenci(i)
                 link = itpenci(j,i)
                 xtpt(j) = xht(link) - xbeg(i)
                 ytpt(j) = yht(link) - ybeg(i)
                 ztpt(j) = zht(link) - zbeg(i)
                 vxtpt(j) = vxht(link) - vxbeg(i)
                 vytpt(j) = vyht(link) - vybeg(i)
                 vztpt(j) = vzht(link) - vzbeg(i)
                 istattmp(j,1) = 0
                 lperit(j) = .false.
              enddo

C...          set up planets at t=0
              call rmvs_step_in_mvpl(i,nbod,mass,xbeg,ybeg,zbeg,
     &                masst,xpltb,ypltb,zpltb,xplte,yplte,zplte)
              i1sttp = 0

              call util_peri(0,nenci(i),xtpt,ytpt,ztpt,vxtpt,
     &             vytpt,vztpt,masst(1),isperit,perit,lperit)

              do j=1,NTPHENC
                 call rmvs_step_in_mvpl(i,nbod,mass,xpl(1,j),ypl(1,j),
     &                                 zpl(1,j),masst,xpltb,ypltb,
     &                                 zpltb,xplte,yplte,zplte)
                 call rmvs_step_in_tp(i1sttp,nbod,nenci(i),masst,
     &              j2rp2,j4rp4,xpltb,ypltb,zpltb,xplte,yplte,zplte,
     &              xtpt,ytpt,ztpt,vxtpt,vytpt,vztpt,istattmp,dti,
     &              i,aoblx(i,j-1),aobly(i,j-1),aoblz(i,j-1),	
     &              aoblx(i,j),aobly(i,j),aoblz(i,j))

                 call util_peri(1,nenci(i),xtpt,ytpt,ztpt,vxtpt,
     &                vytpt,vztpt,masst(1),isperit,perit,lperit)

              enddo

c...          put things back
              do j = 1,nenci(i)
                 link = itpenci(j,i)
                 xht(link) = xtpt(j) + xend(i)
                 yht(link) = ytpt(j) + yend(i)
                 zht(link) = ztpt(j) + zend(i)
                 vxht(link) = vxtpt(j) + vxend(i)
                 vyht(link) = vytpt(j) + vyend(i)
                 vzht(link) = vztpt(j) + vzend(i)
                 if(lperit(j)) then
                    isperi(link) = 0
                 else
                    isperi(link) = isperit(j)
                 endif
                 peri(link) = perit(j)
                 if(istattmp(j,1).ne.0) then
                    istat(link,1) = 1
                    do k=2,NSTAT
                       istat(link,k) = istattmp(j,k)
                    enddo
                 endif
              enddo

           endif
        enddo

        return
        end  ! rmvs2_step_in
c------------------------------------------------------------------




