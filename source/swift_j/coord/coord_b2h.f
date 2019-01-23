c***********************************************************************
c	                    COORD_B2H.F
c***********************************************************************
*     PURPOSE: Converts from Barycentric to Helio coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
*                                 SYMMETRY IN SUBROUTINE CALLS
*		     xb(*),yb(*),zb(*) ==> Barycentric particle coords
*                                          (real array)
*		     vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    xh(*),yh(*),zh(*) ==> Helio particle positions
*                                          (real array)
*                    vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
*                                            (real array)
*       
*     ALGORITHM: Obvious 
*     REMARKS:  Can of course use this to get coords. relative to any body.
*              by changing the one subtracted off.
*
*     Authors:  Martin Duncan
*     WRITTEN:  Jan 27/93
*     REVISIONS: 2/17/95  HFL

	subroutine coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &      xh,yh,zh,vxh,vyh,vzh)


      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c...  Outputs:
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c...  Internals:
	integer n

c----
c...  Executable code 

	do n=1,nbod
	  xh(n) = xb(n) - xb(1)
	  yh(n) = yb(n) - yb(1)
	  zh(n) = zb(n) - zb(1)
	  vxh(n) = vxb(n) - vxb(1)
	  vyh(n) = vyb(n) - vyb(1)
	  vzh(n) = vzb(n) - vzb(1)
	enddo

	return
	end     ! coord_b2h

c--------------------------------------------------------------------------

