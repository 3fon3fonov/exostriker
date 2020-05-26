c*************************************************************************
c                            UTIL_PERI.F
c*************************************************************************
c This subroutine determines whether peri has taken place
c
c             Input:
c                 iflg           ==>  = 0 if first step; = 1 not (int scalar)
c                 ntp          ==>  number of bodies (int scalar)
c                 xt,yt,zt       ==>  planocantric position of tp's
c                                       (real arrays)
c                 vxt,vyt,vzt    ==>  planocantric velcocities of tp's
c                                       (real arrays)
c                 massc          ==>  mass of the central body (real scalar)
c
c             Output:
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real array)
c                lperi           ==> set to .true. if isperi=0
c                                         (logical*2 array)
c
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/25/94
c Last revision: 7/14/94

      subroutine util_peri(iflg,ntp,xt,yt,zt,vxt,vyt,vzt,
     &     massc,isperi,peri,lperi)

      include '../swift.inc'

c...  Inputs Only: 
      integer ntp,iflg
      real*8 xt(ntp),yt(ntp),zt(ntp),massc
      real*8 vxt(ntp),vyt(ntp),vzt(ntp)

c...  Outputs:
      real*8 peri(NTPMAX)
      integer isperi(NTPMAX)
      logical*2 lperi(NTPMAX)

c...  Internals
      integer i,ialpha
      real*8 vdotr,a,e

c----
c...  Executable code 

      if(iflg.eq.0) then    ! are we just setting thing up?

         do i=1,ntp
            vdotr = xt(i)*vxt(i) + yt(i)*vyt(i) + zt(i)*vzt(i)
            if (vdotr .gt. 0.d0) then
               isperi(i) = 1
            else 
               isperi(i) =-1
            endif
         enddo

      else

         do i=1,ntp
            vdotr = xt(i)*vxt(i) + yt(i)*vyt(i) + zt(i)*vzt(i)
            if(isperi(i).eq.-1) then         ! was coming in

               if (vdotr .lt. 0.d0) then     ! still coming in
                  isperi(i) = -1
               else                          ! turned around
                  isperi(i) = 0
                  lperi(i) = .true.
                  call orbel_xv2aeq(xt(i),yt(i),zt(i),vxt(i),vyt(i),
     &                 vzt(i),massc,ialpha,a,e,peri(i))
               endif

            else

               if (vdotr .lt. 0.d0) then     ! coming in
                  isperi(i) = -1
               else
                  isperi(i) = 1              ! going out
               endif

            endif
         enddo

      endif

      return
      end    ! util_peri
c------------------------------------------------------------------


