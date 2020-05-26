c*************************************************************************
c                        DRIFT_KEPU_FCHK.F
c*************************************************************************
c Returns the value of the function f of which we are trying to find the root
c in universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalar)
c                 r0            ==>  Distance between `Sun' and particle
c                                     (real scalar)
c                 mu            ==>  Reduced mass of system (real scalar)
c                 alpha         ==>  Twice the binding energy (real scalar)
c                 u             ==>  Vel. dot radial vector (real scalar)
c                 s             ==>  Approx. root of f 
c             Output:
c                 f             ==>  function value ( = 0 if O.K.) (integer)
c
c Author:  Martin Duncan  
c Date:    March 12/93
c Last revision: March 12/93

      subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)

c...  Inputs: 
      real*8 dt,r0,mu,alpha,u,s

c...  Outputs:
      real*8 f

c...  Internals:
      real*8  x,c0,c1,c2,c3

c----
c...  Executable code 

        x=s*s*alpha
        call drift_kepu_stumpff(x,c0,c1,c2,c3)
        c1=c1*s
        c2 = c2*s*s
        c3 = c3*s*s*s
        f = r0*c1 + u*c2 + mu*c3 - dt

        return
        end     !   drift_kepu_fchk
c-------------------------------------------------------------------
