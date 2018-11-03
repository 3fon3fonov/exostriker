c*******************************************************************
c                           TU4_LDRIFT.F
c*********************************************************************
c LDRIFT gives a 'linear drift' to the bary. positions for use in 
c the symplectic code based on Ham. of the form T(p) + V(x).
c 
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 xb,yb,zb      ==>  initial position in beri coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial velocity in beri coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xb,yb,zb      ==>  final position in beri coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  final velocity in beri coord 
c                                       (real arrays)
c
c Remarks:  Based on Martin's NB4M routines 
c Authors:  Martin Duncan 
c Date:    3/8/93
c Last revision: 

      subroutine tu4_ldrift(nbod,xb,yb,zb,vxb,vyb,vzb,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 dt

c...  Inputs and Outputs:
      real*8 xb(nbod),yb(nbod),zb(nbod)
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Internals
      integer i

c---
c...  Executable code

      do i=1,nbod
         xb(i) = xb(i) + vxb(i)*dt
         yb(i) = yb(i) + vyb(i)*dt
         zb(i) = zb(i) + vzb(i)*dt
      enddo

      return
      end                       ! tu4_ldrift
c____________________________________________________________________________
