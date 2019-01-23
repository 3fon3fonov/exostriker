***********************************************************************
c                    ORBEL_EHYBRID.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           m ==> mean anomaly. (real scalar)
*             Returns:
*              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: For e < 0.18 uses fast routine ESOLMD 
*	         For larger e but less than 0.8, uses EGET
*	         For e > 0.8 uses EHIE
*     REMARKS: Only EHIE brings M and E into range (0,TWOPI)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 25,1992.
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_ehybrid(e,m)

      include '../swift.inc'

c...  Inputs Only: 
	real*8 e,m

c...  Internals:
	real*8 orbel_esolmd,orbel_eget,orbel_ehie

c----
c...  Executable code 

	if(e .lt. 0.18d0) then
	  orbel_ehybrid = orbel_esolmd(e,m)
	else 
	  if( e .le. 0.8d0) then
	     orbel_ehybrid = orbel_eget(e,m)
	  else
	     orbel_ehybrid = orbel_ehie(e,m)
	  endif
	endif   

	return
	end     ! orbel_ehybrid

c--------------------------------------------------------------------




