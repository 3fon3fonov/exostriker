c***********************************************************************
c	                    COORD_B2J.F
c***********************************************************************
*     PURPOSE: Converts from Barycentric to Jacobi coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     xb(*),yb(*),zb(*) ==> Barycentric particle coords
*                                          (real array)
*		     vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    xj(*),yj(*),zj(*) ==> jacobi particle positions
*                                          (real array)
*                    vxj(*),vyj(*),vzj(*) ==> jacobi particle velocities
*                                            (real array)
*       
*     ALGORITHM:  See e.g. Wisdom and Holman 
*     Authors:  Martin Duncan
*     WRITTEN:  Jan 29/93
*     REVISIONS: 2/17/95  HFL

	subroutine coord_b2j(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &      xj,yj,zj,vxj,vyj,vzj)


      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c...  Outputs:
	real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)

c...  Internals:
	real*8 eta(NPLMAX)
	real*8 sumx,sumy,sumz,sumvx,sumvy,sumvz	
	real*8 capx,capy,capz,capvx,capvy,capvz
	integer n

c----
c...  Executable code 


c First calc. the array eta(*) then convert to jacobi coords

	eta(1) = mass(1)
	do n = 2,nbod
	  eta(n) = eta(n-1) + mass(n)
        enddo

	xj(1) = 0.d0
	yj(1) = 0.d0
	zj(1) = 0.d0
	vxj(1) = 0.d0
	vyj(1) = 0.d0
	vzj(1) = 0.d0

	sumx = mass(1)*xb(1)
	sumy = mass(1)*yb(1)
	sumz = mass(1)*zb(1)
	sumvx = mass(1)*vxb(1)
	sumvy = mass(1)*vyb(1)
	sumvz = mass(1)*vzb(1)

	capx = xb(1)
	capy = yb(1)
	capz = zb(1)
	capvx = vxb(1)
	capvy = vyb(1)
	capvz = vzb(1)

	do n=2,nbod - 1
	  xj(n) = xb(n) - capx
	  yj(n) = yb(n) - capy
	  zj(n) = zb(n) - capz
	  vxj(n) = vxb(n) - capvx
	  vyj(n) = vyb(n) - capvy
	  vzj(n) = vzb(n) - capvz

	  sumx = sumx + mass(n)*xb(n)
	  sumy = sumy + mass(n)*yb(n)
	  sumz = sumz + mass(n)*zb(n)
	  sumvx = sumvx + mass(n)*vxb(n)
	  sumvy = sumvy + mass(n)*vyb(n)
	  sumvz = sumvz + mass(n)*vzb(n)

	  capx = sumx/eta(n)
	  capy = sumy/eta(n)
	  capz = sumz/eta(n)
	  capvx =sumvx/eta(n)
	  capvy =sumvy/eta(n)
	  capvz =sumvz/eta(n)

	enddo

	  xj(nbod) = xb(nbod) - capx
	  yj(nbod) = yb(nbod) - capy
	  zj(nbod) = zb(nbod) - capz
	  vxj(nbod) = vxb(nbod) - capvx
	  vyj(nbod) = vyb(nbod) - capvy
	  vzj(nbod) = vzb(nbod) - capvz

	return
	end   ! coord_b2j

c--------------------------------------------------------------------------

