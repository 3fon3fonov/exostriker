c*************************************************************************
c                            RMVS_PERI.F
c*************************************************************************
c This subroutine determines whether peri has taken place
c
c             Input:
c                 iflg           ==>  = 0 if first step; = 1 not (int scalar)
c                 nenci          ==>  number of massive bodies (int scalar)
c                 xt,yt,zt       ==>  planocantric position of tp's
c                                       (real arrays)
c                 vxt,vyt,vzt    ==>  planocantric velcocities of tp's
c                                       (real arrays)
c                 massp          ==>  mass of the planet (real scalar)
c
c             Output:
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real array)
c
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/25/94
c Last revision: 

      subroutine rmvs_peri(iflg,nenci,xt,yt,zt,vxt,vyt,vzt,
     &     massp,isperi,peri)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nenci,iflg
      real*8 xt(nenci),yt(nenci),zt(nenci),massp
      real*8 vxt(nenci),vyt(nenci),vzt(nenci)

c...  Outputs:
      real*8 peri(NTPMAX)
      integer isperi(NTPMAX)

c...  Internals
      integer i,ialpha
      real*8 vdotr,a,e

c----
c...  Executable code 

      if(iflg.eq.0) then    ! are we just setting thing up?

         do i=1,nenci
            vdotr = xt(i)*vxt(i) + yt(i)*vyt(i) + zt(i)*vzt(i)
            if (vdotr .gt. 0.d0) then
               isperi(i) = 1
            else 
               isperi(i) =-1
            endif
         enddo

      else

         do i=1,nenci
            if(isperi(i).eq.-1)  then   ! still coming in
               vdotr = xt(i)*vxt(i) + yt(i)*vyt(i) + zt(i)*vzt(i)
               if (vdotr .gt. 0.d0) then
                  isperi(i) = 0
                  call orbel_xv2aeq(xt(i),yt(i),zt(i),vxt(i),vyt(i),
     &                 vzt(i),massp,ialpha,a,e,peri(i))
               endif
            endif
         enddo

      endif

      return
      end    ! rmvs_peri
c------------------------------------------------------------------


