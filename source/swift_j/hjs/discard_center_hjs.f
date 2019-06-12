c*************************************************************************
c                            DISCARD_CENTER_HJS.F
c*************************************************************************
c This subroutine checks to see if a partical should be discarded because
c of its position or becuase it becomes unbound
c
c             Input:
c                 time          ==>  current time
c                 ntp,itp       ==>  number of test particles (int scalar)
c                 msys          ==>  Total Planet Mass (real scalar)
c                 xjt,yjt,zjt    ==>   part position in Jacobi coord 
c                                      (real arrays)
c                 xbt,ybt,zbt    ==>   part position in bary coord 
c                                      (real arrays)
c                 vxbt,vybt,vzbt ==>   part velocity in bary coord 
c                                        (real arrays)
c                 rmin,rmax      ==>  maximum and min distance from Sun
c                                     if <0  then don't check
c                                        (real scalar)
c                 rmaxu          ==>  maximum distance from Sun in not bound
c                                     if <0  then don't check
c                                        (real scalar)
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c                                      istat(i,2) = -2  -> a<0 & r_bary>rmaxu
c	                               istat(i,2) = -3     r_helio>rmax
c	                               istat(i,2) =  1     r_helio<rmin
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,1) time of discard.
c
c
c Remarks: Adapted from discard_sun
c Authors:  H. Beust
c Date:    Feb. 15, 2002

      subroutine discard_center_hjs(time,ntp,itp,msys,xjt,yjt,zjt,
     &     xbt,ybt,zbt,vxbt,vybt,vzbt,rmin,rmax,rmaxu,istat,rstat)

      include '../swift.inc'

c...  Inputs: 
      integer ntp,itp(ntp)
      real*8 msys
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 xbt(ntp),ybt(ntp),zbt(ntp)
      real*8 vxbt(ntp),vybt(ntp),vzbt(ntp)
      real*8 rmin,rmax,rmaxu,time

c...  Input and Output
      integer istat(ntp,NSTAT)
      real*8 rstat(ntp,NSTATR)

c...  Internal
      integer i
      real*8 energy,vb2,rb2,rj2
      real*8 rmin2,rmax2,rmaxu2

c-----
c...  Executable code 

      rmin2 = rmin*rmin
      rmax2 = rmax*rmax
      rmaxu2 = rmaxu*rmaxu

      do i=1,ntp
         if(istat(i,1).eq.0) then
            rj2 = xjt(i)**2 + yjt(i)**2 + zjt(i)**2
            if( (rmax.ge.0.0) .and. (rj2.gt.rmax2) ) then
               write(*,*) 'Particle',itp(i),
     &                     ' too far from Center at t=',time
               istat(i,1) = 1
               istat(i,2) = -3
               rstat(i,1) = time
            endif
            if( (rmin.ge.0.0) .and. (rj2.lt.rmin2) ) then
               write(*,*) 'Particle',itp(i),
     &                     ' too close to Center at t=',time
               istat(i,1) = 1
               istat(i,2) = 1
               rstat(i,1) = time
            endif
         endif
         if( (istat(i,1).eq.0) .and. (rmaxu.ge.0.0) ) then
            rb2 = xbt(i)**2 + ybt(i)**2 + zbt(i)**2
            vb2 = vxbt(i)**2 + vybt(i)**2 + vzbt(i)**2
            energy = 0.5*vb2 - msys/sqrt(rb2)
            if( (energy.gt.0.0) .and. (rb2.gt.rmaxu2) ) then
               write(*,*) 'Particle',itp(i),
     & ' is unbound and too far from barycenter at t=',time
               istat(i,1) = 1
               istat(i,2) = -2
               rstat(i,1) = time
            endif
         endif
      enddo

      return
      end     ! discard_center_hjs
c------------------------------------------------------------------



