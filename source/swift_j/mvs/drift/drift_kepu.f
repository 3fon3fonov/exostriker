c*************************************************************************
c                        DRIFT_KEPU.F
c*************************************************************************
c subroutine for solving kepler's equation using universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 fp            ==>  f' from p170  
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflg          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison  
c Date:    2/3/93
c Last revision: 2/3/93

      subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

      include '../../swift.inc'

c...  Inputs: 
      real*8 dt,r0,mu,alpha,u

c...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflg

c...  Internals:
      real*8 s,st,fo,fn

c----
c...  Executable code 

        call drift_kepu_guess(dt,r0,mu,alpha,u,s)
         
        st = s
c..     store initial guess for possible use later in
c..     laguerre's method, in case newton's method fails.

        call drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
        if(iflg.ne.0) then
           call drift_kepu_fchk(dt,r0,mu,alpha,u,st,fo)
           call drift_kepu_fchk(dt,r0,mu,alpha,u,s,fn)
           if(abs(fo).lt.abs(fn)) then
               s = st 
           endif
           call drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
        endif

        return
        end    ! drift_kepu
c----------------------------------------------------------------------
