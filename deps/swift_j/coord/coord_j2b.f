c***********************************************************************
c	                    COORD_J2B.F
c***********************************************************************
*     PURPOSE: Converts from Jacobi to Barycentric coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     xj(*),yj(*),zj(*) ==> Jacobi particle coords
*                                          (real array)
*		     vxj(*),vyj(*),vzj(*) ==> Jacobi particle velocities
*                                             (real array)
*                 Returned are
*                    xb(*),yb(*),zb(*) ==>  Barycentric particle positions
*                                          (real array)
*                    vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                            (real array)
*       
*     ALGORITHM: See my notes on Nov 21. 
*
*     Authors:  Martin Duncan
*     WRITTEN:  Jan 27/93
*     REVISIONS: 2/17/95  HFL

	subroutine coord_j2b(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xb,yb,zb,vxb,vyb,vzb)

      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)

c...  Outputs:
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c...  Internals:
	integer n
	real*8 xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp
	real*8 capx,capy,capz,capvx,capvy,capvz
	real*8 mtot,rat,rat2
	real*8 eta(NPLMAX)
c----
c...  Executable code 

c First compute the necessary auxiliary vbles, then compute the bary. positions

	eta(1) = mass(1)
	do n = 2,nbod
	  eta(n) = eta(n-1) + mass(n)
        enddo

	mtot = eta(nbod)
	xb(nbod) = eta(nbod-1)*xj(nbod)/mtot
	yb(nbod) = eta(nbod-1)*yj(nbod)/mtot
	zb(nbod) = eta(nbod-1)*zj(nbod)/mtot
	vxb(nbod) = eta(nbod-1)*vxj(nbod)/mtot
	vyb(nbod) = eta(nbod-1)*vyj(nbod)/mtot
	vzb(nbod) = eta(nbod-1)*vzj(nbod)/mtot
	
	capx = mass(nbod)*xj(nbod)/mtot
	capy = mass(nbod)*yj(nbod)/mtot
	capz = mass(nbod)*zj(nbod)/mtot
	capvx = mass(nbod)*vxj(nbod)/mtot
	capvy = mass(nbod)*vyj(nbod)/mtot
	capvz = mass(nbod)*vzj(nbod)/mtot

	do n = nbod-1,2,-1
	  rat = eta(n-1)/eta(n)
	  xb(n) = rat*xj(n) - capx
	  yb(n) = rat*yj(n) - capy
	  zb(n) = rat*zj(n) - capz
	  vxb(n) = rat*vxj(n) - capvx
	  vyb(n) = rat*vyj(n) - capvy
	  vzb(n) = rat*vzj(n) - capvz

	  rat2 = mass(n)/eta(n)
	  capx = capx + rat2*xj(n)
	  capy = capy + rat2*yj(n)
	  capz = capz + rat2*zj(n)
	  capvx = capvx + rat2*vxj(n)
	  capvy = capvy + rat2*vyj(n)
	  capvz = capvz + rat2*vzj(n)
	enddo
	
c Now compute the Sun's barycentric position
	xtmp =0.d0
	ytmp =0.d0
	ztmp =0.d0
	vxtmp =0.d0
	vytmp =0.d0
	vztmp =0.d0

	do n=2,nbod
	   xtmp = xtmp + mass(n)*xb(n)
	   ytmp = ytmp + mass(n)*yb(n)
	   ztmp = ztmp + mass(n)*zb(n)
	   vxtmp = vxtmp + mass(n)*vxb(n)
	   vytmp = vytmp + mass(n)*vyb(n)
	   vztmp = vztmp + mass(n)*vzb(n)
	enddo

	xb(1) = -xtmp/mass(1)
	yb(1) = -ytmp/mass(1)
	zb(1) = -ztmp/mass(1)
	vxb(1) = -vxtmp/mass(1)
	vyb(1) = -vytmp/mass(1)
	vzb(1) = -vztmp/mass(1)

	return
	end    ! coord_j2b

c--------------------------------------------------------------------------
