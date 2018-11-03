***********************************************************************
c                    ORBEL_EHIE.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           m ==> mean anomaly. (real scalar)
*             Returns:
*              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Use Danby's quartic for 3 iterations. 
*                Eqn. is f(x) = x - e*sin(x+M). Note  that
*	         E = x + M. First guess is very good for e near 1.
*	         Need to first get M between 0. and PI and use
*		 symmetry to return right answer if M between PI and 2PI
*     REMARKS: Modifies M so that both E and M are in range (0,TWOPI)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 25,1992.
*     REVISIONS: 
***********************************************************************

      real*8 function orbel_ehie(e,m)

      include '../swift.inc'

c...  Inputs Only: 
	real*8 e,m

c...  Internals:
      integer iflag,nper,niter,NMAX
      real*8 dx,x,sa,ca,esa,eca,f,fp

      parameter (NMAX = 3)

c----
c...  Executable code 

c In this section, bring M into the range (0,TWOPI) and if
c the result is greater than PI, solve for (TWOPI - M).
	iflag = 0
	nper = m/TWOPI
	m = m - nper*TWOPI
	if (m .lt. 0.d0) m = m + TWOPI

	if (m.gt.PI) then
	   m = TWOPI - m
	   iflag = 1
	endif

c Make a first guess that works well for e near 1.
	x = (6.d0*m)**(1.d0/3.d0) - m
	niter =0

c Iteration loop
	do niter =1,NMAX
	    call orbel_scget(x + m,sa,ca)
	    esa = e*sa
	    eca = e*ca
	    f = x - esa
	    fp = 1.d0 -eca
	    dx = -f/fp
	    dx = -f/(fp + 0.5d0*dx*esa)
	    dx = -f/(fp + 0.5d0*dx*(esa+0.3333333333333333d0*eca*dx))
	    x = x + dx
	enddo

	orbel_ehie = m + x

	if (iflag.eq.1) then
	  orbel_ehie = TWOPI - orbel_ehie
	  m = TWOPI - m
	endif

	return
	end         !orbel_ehie
c------------------------------------------------------------------

