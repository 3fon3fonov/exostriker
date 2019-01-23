c*************************************************************************
c                        LYAP2_DRIFT_TP.F
c*************************************************************************
c This subroutine loops thorugh the TEST particles and calls the danby routine
c
c             INPUT:
c                 ntp             ==>  number of test particles (int scalar)
c                 msun            ==>  mass of the sun (real scalar)
c                 xjt,yjt,zjt     ==>  initial position in jacobi coord 
c                                      (real arrays)
c                 vxjt,vyjt,vzjt  ==>  initial position in jacobi coord 
c                                      (real arrays)
c               dxjt,dyjt,dzjt    ==>  initial separation in position
c                                      (real arrays)
c               dvxjt,dvyjt,dvzjt ==>  initial separation in velocity
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt              ==>  time step
c             OUTPUT:
c                 xjt,yjt,zjt     ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxjt,vyjt,vzjt  ==>  final position in jacobi coord 
c                                      (real arrays)
c               dxjt,dyjt,dzjt    ==>  final separation in position
c                                      (real arrays)
c               dvxjt,dvyjt,dvzjt ==>  final separation in velocity
c                                       (real arrays)
c
c Comments: based on drift_tp.f
c Authors:  Hal Levison 
c Date:    7/11/95
c Last revision:

      subroutine lyap2_drift_tp(ntp,msun,xjt,yjt,zjt,vxjt,vyjt,
     &     vzjt,dxjt,dyjt,dzjt,dvxjt,dvyjt,dvzjt,dt,istat)	

      include '../swift.inc'

c...  Inputs Only: 
      integer ntp
      real*8 msun,dt

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)
      real*8 dxjt(ntp),dyjt(ntp),dzjt(ntp)
      real*8 dvxjt(ntp),dvyjt(ntp),dvzjt(ntp)

c...  Internals:
	integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth

	do j = 1,ntp
           if(istat(j,1).eq.0) then
	      call lyap2_drift_one(msun,xjt(j),yjt(j),zjt(j),
     &             vxjt(j),vyjt(j),vzjt(j),dxjt(j),dyjt(j),dzjt(j),
     &             dvxjt(j),dvyjt(j),dvzjt(j),dt,iflg)
              if(iflg.ne.0) then
                 istat(j,1) = 1
                 istat(j,2) = -1
              endif
           endif
	enddo

	return
	end
c--------------------------------------------------------------------------

