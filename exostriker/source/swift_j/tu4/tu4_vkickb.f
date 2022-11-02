c*************************************************************************
c                         TU4_VKICKB.F
c*************************************************************************
c VKICKB gives a 'kick' to the bary. velocities for use in the symplectic
c code based on Ham. of the form T(p) + V(x).
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 vxb,vyb,vzb   ==>  initial velocity in beri coord 
c                                    (real arrays)
c                 axb,ayb,azb   ==>  accel in beri coord (real arrays) 
c                 dt            ==>  time step
c             Output:
c                 vxb,vyb,vzb   ==>  final velocity in beri coord 
c                                       (real arrays)
c
c Remarks:  Based on Martin's NB4M routines 
c Authors:  Martin Duncan 
c Date:    3/8/93
c Last revision: 

	  subroutine tu4_vkickb(nbod,vxb,vyb,vzb,axb,ayb,azb,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 axb(nbod),ayb(nbod),azb(nbod),dt

c...  Inputs and Outputs:
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Internals
      integer i

c---
c...  Executable code

	  do i=1,nbod
	     vxb(i) = vxb(i) + axb(i)*dt
	     vyb(i) = vyb(i) + ayb(i)*dt
	     vzb(i) = vzb(i) + azb(i)*dt
	  enddo

	  return
	  end   ! tu4_vkickb
c____________________________________________________________________________
