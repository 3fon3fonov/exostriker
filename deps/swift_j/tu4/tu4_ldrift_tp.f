c*******************************************************************
c                           TU4_LDRIFT_TP.F
c*********************************************************************
c LDRIFT gives a 'linear drift' to the bary. positions for use in 
c the symplectic code based on Ham. of the form T(p) + V(x).
c  TEST PARTICLES
c
c             Input:
c                 ntp              ==>  number of bodies (int scalar)
c                 xbt,ybt,zbt      ==>  initial position in beri coord 
c                                    (real arrays)
c                 vxbt,vybt,vzbt   ==>  initial velocity in beri coord 
c                                    (real arrays)
c                 dt               ==>  time step
c                 istat            ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c                 xbt,ybt,zbt      ==>  final position in beri coord 
c                                    (real arrays)
c                 vxbt,vybt,vzbt   ==>  final velocity in beri coord 
c                                       (real arrays)
c
c Remarks:  Based on Martin's NB4M routines 
c Authors:  Martin Duncan 
c Date:    3/8/93
c Last revision: 

      subroutine tu4_ldrift_tp(ntp,xbt,ybt,zbt,vxbt,vybt,vzbt,
     &            dt,istat)

      include '../swift.inc'

c...  Inputs Only: 
      integer ntp,istat(ntp)
      real*8 dt

c...  Inputs and Outputs:
      real*8 xbt(ntp),ybt(ntp),zbt(ntp)
      real*8 vxbt(ntp),vybt(ntp),vzbt(ntp)

c...  Internals
      integer i

c---
c...  Executable code

      do i=1,ntp
         if(istat(i).eq.0) then
            xbt(i) = xbt(i) + vxbt(i)*dt
            ybt(i) = ybt(i) + vybt(i)*dt
            zbt(i) = zbt(i) + vzbt(i)*dt
         endif
      enddo

      return
      end                       ! tu4_ldrift_tp
c____________________________________________________________________________

