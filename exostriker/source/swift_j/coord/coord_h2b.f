c***********************************************************************
c	                    COORD_H2B.F
c***********************************************************************
*     PURPOSE: Converts from Heliocentric to Barycentric coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     xh(*),yh(*),zh(*) ==> heliocentric particle coords
*                                          (real array)
*		     vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
*                                             (real array)
*                 Returned are
*                    xb(*),yb(*),zb(*) ==> bary. particle positions
*                                          (real array)
*                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
*                                            (real array)
*                    msys              ==>  Total mass of of system
*                                            (real scalar)       
*     Authors:  Martin Duncan
*     ALGORITHM: Obvious 
*     WRITTEN:  Jan 27/93
*     REVISIONS: 2/22/94  HFL

	subroutine coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &      xb,yb,zb,vxb,vyb,vzb,msys)

      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c...  Outputs:
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c...  Internals:
	real*8 msys,xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp
	integer n

c----
c...  Executable code 

	msys = mass(1)
	xtmp =0.d0
	ytmp =0.d0
	ztmp =0.d0
	vxtmp =0.d0
	vytmp =0.d0
	vztmp =0.d0

	do n=2,nbod
	   msys = msys +mass(n)
	   xtmp = xtmp + mass(n)*xh(n)
	   ytmp = ytmp + mass(n)*yh(n)
	   ztmp = ztmp + mass(n)*zh(n)
	   vxtmp = vxtmp + mass(n)*vxh(n)
	   vytmp = vytmp + mass(n)*vyh(n)
	   vztmp = vztmp + mass(n)*vzh(n)
	enddo

	xb(1) = -xtmp/msys
	yb(1) = -ytmp/msys
	zb(1) = -ztmp/msys
	vxb(1) = -vxtmp/msys
	vyb(1) = -vytmp/msys
	vzb(1) = -vztmp/msys

	do n=2,nbod
	  xb(n) = xh(n) + xb(1)
	  yb(n) = yh(n) + yb(1)
	  zb(n) = zh(n) + zb(1)
	  vxb(n) = vxh(n) + vxb(1)
	  vyb(n) = vyh(n) + vyb(1)
	  vzb(n) = vzh(n) + vzb(1)
	enddo

	return
	end     ! coord_h2b
c--------------------------------------------------------------------------

