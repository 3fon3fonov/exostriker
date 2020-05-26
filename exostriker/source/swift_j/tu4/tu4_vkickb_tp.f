c*************************************************************************
c                         TU4_VKICKB_TP.F
c*************************************************************************
c VKICKB gives a 'kick' to the bary. velocities for use in the symplectic
c code based on Ham. of the form T(p) + V(x). TEST PARTICLES
c
c             Input:
c                 ntp              ==>  number of massive bodies (int scalar)
c                 vxbt,vybt,vzbt   ==>  initial velocity in beri coord 
c                                    (real arrays)
c                 axbt,aybt,azbt   ==>  accel in beri coord (real arrays) 
c                 dt               ==>  time step
c                 istat            ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c
c             Output:
c                vxbt,vybt,vzbt   ==>  final velocity in beri coord 
c                                       (real arrays)
c
c Remarks:  Based on Martin's NB4M routines
c Authors:  Martin Duncan 
c Date:    3/8/93
c Last revision: 

      subroutine tu4_vkickb_tp(ntp,vxbt,vybt,vzbt,axbt,aybt,
     &           azbt,dt,istat)

      include '../swift.inc'

c...  Inputs Only: 
      integer ntp,istat(ntp)
      real*8 axbt(ntp),aybt(ntp),azbt(ntp),dt

c...  Inputs and Outputs:
      real*8 vxbt(ntp),vybt(ntp),vzbt(ntp)

c...  Internals
      integer i

c---
c...  Executable code

	  do i=1,ntp
             if(istat(i).eq.0) then
                vxbt(i) = vxbt(i) + axbt(i)*dt
                vybt(i) = vybt(i) + aybt(i)*dt
                vzbt(i) = vzbt(i) + azbt(i)*dt
             endif
	  enddo

	  return
	  end   ! tu4_vkickb_tp
c____________________________________________________________________________
