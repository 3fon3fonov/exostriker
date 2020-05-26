c*************************************************************************
c                        LYAP2_DRIFT_KEPU.F
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
c                 c4,c5         ==>  Mikkola's extra c'c
c                                       (real scalors)
c                    s          ==> Value of universal variable
c                                       (real scalor)
c                 iflg          ==>  =0 if converged; !=0 if not
c                                       (int scalor)
c
c Comments: Based on drift_kepu
c Author:  Hal Levison  
c Date:    7/11/95
c Last revision: 

      subroutine lyap2_drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,c4,
     &                      c5,s,iflg)

      include '../swift.inc'

c...  Inputs: 
      real*8 dt,r0,mu,alpha,u

c...  Outputs:
      real*8 fp,c1,c2,c3,c4,c5,s
      integer iflg

c...  Internals:
      real*8 st,fo,fn

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

c...    Now call Mikkola's routines
        call lyap2_kepu_mikk(s,alpha,c1,c2,c3,c4,c5)

        return
        end    ! lyap2_drift_kepu
c----------------------------------------------------------------------
