c*************************************************************************
c                            DISCARD_PERI.F
c*************************************************************************
c This subroutine checks to see if a partical should be discarded because
c of its perihelion distance gets too small
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 ntp           ==>  number of test bodies (int scalar)
c                 xht,yht,zht    ==>   part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>   part vel in helio coord 
c                                      (real arrays)
c                 qmin            ==>  Smallest perihelion distance 
c                                      (real scalar)
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                 nbod            ==>  Number of planets (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c	                               istat(i,2) =  -4
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,2) perihelion distance.
c                                      rstat(i,1) time of discard.
c
c Remarks: 
c
c Authors:  Hal Levison 
c Date:    5/10/94
c Last revision: 1/20/97

         subroutine discard_peri(time,ntp,xht,yht,zht,vxht,vyht,
     &       vzht,qmin,istat,rstat,nbod,mass,xh,yh,zh,vxh,vyh,vzh)

      include '../swift.inc'

c...  Inputs: 
      integer ntp,nbod
      real*8 mass(nbod),time,qmin
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Input and Output
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)

c...  internal
      integer i,isperi(NTPMAX),i1st,j
      real*8 peri(NTPMAX),ih,r2
      logical*2 lperi(NTPMAX)
      real*8 r2hill(NTPMAX)

      data i1st/0/
      save i1st,isperi,r2hill

c-----
c...  Executable code 

      if(i1st.eq.0) then     ! if first time through, set things up
         call util_hills(nbod,mass,xh,yh,zh,vxh,vyh,vzh,r2hill)
         call util_peri(0,ntp,xht,yht,zht,vxht,vyht,vzht,
     &     mass(1),isperi,peri,lperi)
         i1st = 1
         return                 !  <==== RETURN
      endif

      call util_peri(1,ntp,xht,yht,zht,vxht,vyht,vzht,
     &     mass(1),isperi,peri,lperi)

      do i=1,ntp
         if( (istat(i,1).eq.0).and. (isperi(i).eq.0) ) then
            
            ih = 0
            do j=2,nbod
               r2 = (xht(i)-xh(j))**2 + (yht(i)-yh(j))**2 + 
     &              (zht(i)-zh(j))**2
               if(r2.le.r2hill(j)) then
                  ih = 1
               endif
            enddo

            if(ih.eq.0) then
               rstat(i,2) = peri(i)
               if(peri(i).le.qmin) then
                  write(*,*) 'Particle',i,' perihelion distance too',
     &                 ' small at t=',time
                  istat(i,1) =  1
                  istat(i,2) = -4
               endif
            endif

         endif
      enddo

      return
      end       ! discard_peri










