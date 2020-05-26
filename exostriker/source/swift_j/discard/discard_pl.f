c*************************************************************************
c                            DISCARD_PL.F
c*************************************************************************
c This subroutine checks to see if a partical should be discarded because
c of its position or becuase it becomes unbound
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 dt            ==>  time step  (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>   part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>   part vel in helio coord 
c                                      (real arrays)
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,3) planet of last inner enc
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,3) closest approach to a planet
c                                          as determined by encounter routines.
c                                          only set if peri happened last dt.
c
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c	                               istat(i,2) =  n    too close to planet n
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,1) time of discard.
c
c
c Remarks: 
c
c Authors:  Hal Levison 
c Date:    3/2/93
c Last revision: 2/22/94

      subroutine discard_pl(time,dt,nbod,ntp,mass,xh,yh,zh,vxh,
     &       vyh,vzh,xht,yht,zht,vxht,vyht,vzht,rplsq,istat,rstat)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,ntp
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp),rplsq(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 time,dt

c...  Input and Output
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)

c...  internal
      integer i,j,iflg
      real*8 xr,yr,zr,vxr,vyr,vzr,peri2,r2min

c-----
c...  Executable code 

      do j=1,ntp
         if(istat(j,1).eq.0) then

            peri2 = rstat(j,3)*rstat(j,3)
            if( (istat(j,3).ne.0) .and. (rstat(j,3).gt.0.0) .and.
     &           (peri2.lt.rplsq(istat(j,3))) ) then

               istat(j,1) = 1
               istat(j,2) = istat(j,3)
               write(*,*) 'Particle',j,' q with respect to Planet ',
     &              istat(j,3),' is too small at t=',time
               rstat(j,1) = time
               rstat(j,2) = -1.0d0*rstat(j,3)

            else

               do i=2,nbod
                  xr = xht(j) - xh(i)
                  yr = yht(j) - yh(i)
                  zr = zht(j) - zh(i)
                  vxr = vxht(j) - vxh(i)
                  vyr = vyht(j) - vyh(i)
                  vzr = vzht(j) - vzh(i)
                  call discard_pl_close(xr,yr,zr,vxr,vyr,vzr,
     &                 dt,rplsq(i),iflg,r2min)
                  if(iflg.eq.1) then
                     istat(j,1) = 1
                     istat(j,2) = i
                     write(*,*) 'Particle',j,' too close to Planet ',i,
     &                    ' at t=',time
                     rstat(j,1) = time
                     rstat(j,2) = sqrt(r2min)
                  endif
               enddo

            endif
         endif
      enddo

      return
      end                       ! discard_pl
c-------------------------------------------------------------------------


