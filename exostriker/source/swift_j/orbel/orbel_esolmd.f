***********************************************************************
c                    ORBEL_ESOLMD.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           m ==> mean anomaly. (real scalar)
*             Returns:
*                orbel_esolmd ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Some sort of quartic convergence from Wisdom. 
*     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
*         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
*      	  ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI 
*     INCLUDES: needs SCGET.F
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 7, 1992.
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_esolmd(e,m)

      include '../swift.inc'

c...  Inputs Only: 
      real*8 e,m

c...  Internals:
	real*8 x,sm,cm,sx,cx
	real*8 es,ec,f,fp,fpp,fppp,dx

c----
c...  Executable code 

c...    Function to solve Kepler's eqn for E (here called
c...    x) for given e and M. returns value of x.

	  call orbel_scget(m,sm,cm)
	  x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))

	  call orbel_scget(x,sx,cx)
	  es = e*sx
	  ec = e*cx
	  f = x - es  - m
	  fp = 1.d0 - ec 
	  fpp = es 
	  fppp = ec 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)

	  orbel_esolmd = x + dx

	return   ! orbel_esolmd
	end
c--------------------------------------------------------------------
