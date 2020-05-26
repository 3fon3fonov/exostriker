c*************************************************************************
c                        KICKVH_TP.F
c*************************************************************************
c To kick the velocity components vxh(*) by axh(*)*dt for test particles 
c
c             Input:
c                 ntp          ==>  number of bodies (int scalar)
c                 vxh,vyh,vzh  ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 axh,ayh,azh  ==>  acceleration in helio coord
c                                    (real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c                          dt   ==>  time step
c             Output:
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                    (real arrays)
c
c     ALGORITHM: Obvious  
c       
c     AUTHOR:  M. Duncan.
c     DATE WRITTEN:  Feb. 2, 1993.
c     REVISIONS: 2/18/93   HFL

      subroutine kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dt) 

      include '../../swift.inc'

c...  Inputs Only: 
	integer ntp,istat(ntp)
	real*8 axht(ntp),ayht(ntp),azht(ntp)
	real*8 dt

c...   Inputs and Output:
	real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
	integer n

c----
c...  Executable code 

	do n= 1, ntp
           if(istat(n).eq.0) then
              vxht(n) = vxht(n) + axht(n)*dt
              vyht(n) = vyht(n) + ayht(n)*dt
              vzht(n) = vzht(n) + azht(n)*dt
           endif
	enddo

        return
        end    ! kickvh_tp
c-----------------------------------------------------------------------------
