c*************************************************************************
c                        DRIFT_KEPU_LAG.F
c*************************************************************************
c subroutine for solving kepler's equation in universal variables.
c using LAGUERRE'S METHOD
c
c             Input:
c                 s             ==>  inital value of universal variable
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 s             ==>  final value of universal variable
c                 fp            ==>  f' from p170  
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison  
c Date:    2/3/93
c Last revision: 4/21/93

      subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

      include '../../swift.inc'

c...  Inputs: 
      real*8 s,dt,r0,mu,alpha,u

c...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflg

c...  Internals:
      integer nc,ncmax
      real*8 ln
      real*8 x,fpp,ds,c0,f
      real*8 fdt

      integer NTMP
      parameter(NTMP=NLAG2+1)

c----
c...  Executable code 

c...    To get close approch needed to take lots of iterations if alpha<0
        if(alpha.lt.0.0) then
           ncmax = NLAG2
        else
           ncmax = NLAG2
        endif

        ln = 5.0
c...    start laguere's method
        do nc =0,ncmax
           x = s*s*alpha
           call drift_kepu_stumpff(x,c0,c1,c2,c3)
           c1 = c1*s 
           c2 = c2*s*s 
           c3 = c3*s*s*s
           f = r0*c1 + u*c2 + mu*c3 - dt
           fp = r0*c0 + u*c1 + mu*c2
           fpp = (-40.0*alpha + mu)*c1 + u*c0
           ds = - ln*f/(fp + dsign(1.d0,fp)*sqrt(abs((ln - 1.0)*
     &       (ln - 1.0)*fp*fp - (ln - 1.0)*ln*f*fpp)))
           s = s + ds

           fdt = f/dt

c..        quartic convergence
           if( fdt*fdt.lt.DANBYB*DANBYB) then 
             iflg = 0
             return
           endif
c...      Laguerre's method succeeded
        enddo

        iflg = 2

        return

        end    !    drift_kepu_leg
c-----------------------------------------------------------------------
