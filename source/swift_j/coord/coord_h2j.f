***********************************************************************
*	                    COORD_H2J.F
***********************************************************************
*     PURPOSE: Converts from Heliocentric to Jacobi coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of  objects (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==> planetary masses (real array)
*		     xh(*),yh(*),zh(*) ==> barycentric particle coords
*                                            (real array)
*		     vxh(*),vyh(*),vzh(*) ==>barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    xj(*),yj(*),zj(*) ==> jacobi. particle positions
*                                             (real array)
*                    vxj(*),vyj(*),vzj(*) ==> jacobi. particle velocities
*                                              (real array)
*       
*     ALGORITHM: See my notes Nov 21/92 
*     REMARKS:  Note that we set the Jacobi coord of the Sun = 0
*               This is not in accord with the definition, but
*               since we never use the Sun's Jacobi coord, it was the
*               fastest thing to do.
*       
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  Jan 27, 1993.
*     REVISIONS:  2/17/93 HFL

	subroutine coord_h2j(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &      xj,yj,zj,vxj,vyj,vzj)

      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

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

	xj(1) =0.d0
	yj(1) =0.d0
	zj(1) =0.d0
	vxj(1) =0.d0
	vyj(1) =0.d0
	vzj(1) =0.d0

	xj(2) = xh(2)
	yj(2) = yh(2)
	zj(2) = zh(2)
	vxj(2) = vxh(2)
	vyj(2) = vyh(2)
	vzj(2) = vzh(2)

	sumx = mass(2)*xh(2)
	sumy = mass(2)*yh(2)
	sumz = mass(2)*zh(2)
	sumvx = mass(2)*vxh(2)
	sumvy = mass(2)*vyh(2)
	sumvz = mass(2)*vzh(2)

	capx = sumx/eta(2)
	capy = sumy/eta(2)
	capz = sumz/eta(2)
	capvx = sumvx/eta(2)
	capvy = sumvy/eta(2)
	capvz = sumvz/eta(2)

	do n=3,nbod
	  xj(n) = xh(n) - capx
	  yj(n) = yh(n) - capy
	  zj(n) = zh(n) - capz
	  vxj(n) = vxh(n) - capvx
	  vyj(n) = vyh(n) - capvy
	  vzj(n) = vzh(n) - capvz

	  if(n.lt.nbod) then
             sumx = sumx + mass(n)*xh(n)
             sumy = sumy + mass(n)*yh(n)
             sumz = sumz + mass(n)*zh(n)
             sumvx = sumvx + mass(n)*vxh(n)
             sumvy = sumvy + mass(n)*vyh(n)
             sumvz = sumvz + mass(n)*vzh(n)

             capx = sumx/eta(n)
             capy = sumy/eta(n)
             capz = sumz/eta(n)
             capvx =sumvx/eta(n)
             capvy =sumvy/eta(n)
             capvz =sumvz/eta(n)
          endif

	enddo

  	return
	end    ! coord_h2j
c--------------------------------------------------------------------
