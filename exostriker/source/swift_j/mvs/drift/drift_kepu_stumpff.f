c*************************************************************************
c                        DRIFT_KEPU_STUMPFF.F
c*************************************************************************
c subroutine for the calculation of stumpff functions
c see Danby p.172  equations 6.9.15
c
c             Input:
c                 x             ==>  argument
c             Output:
c                 c0,c1,c2,c3   ==>  c's from p171-172
c                                       (real scalors)
c Author:  Hal Levison  
c Date:    2/3/93
c Last revision: 2/3/93

      subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)

      include '../../swift.inc'

c...  Inputs: 
      real*8 x

c...  Outputs:
      real*8 c0,c1,c2,c3

c...  Internals:
      integer n,i
      real*8 xm

c----
c...  Executable code 

      n = 0
      xm = 0.1
      do while(abs(x).ge.xm)
         n = n + 1
         x = x/4.0
      enddo

      c2 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/182.)
     &       /132.)/90.)/56.)/30.)/12.)/2.
      c3 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/210.)
     &       /156.)/110.)/72.)/42.)/20.)/6.
      c1 = 1. - x*c3
      c0 = 1. - x*c2

      if(n.ne.0) then
         do i=n,1,-1
            c3 = (c2 + c0*c3)/4.
            c2 = c1*c1/2.
            c1 = c0*c1
            c0 = 2.*c0*c0 - 1.
            x = x * 4.
          enddo
       endif

       return
       end     !   drift_kepu_stumpff
c------------------------------------------------------------------
