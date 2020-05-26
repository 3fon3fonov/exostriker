***********************************************************************
*	                    COORD_VH2VJ.F
***********************************************************************
*     PURPOSE: Converts from Heliocentric to Jacobi coords. VELOCITIES ONLY.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of  objects (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==> planetary masses (real array)
*		     vxh(*),vyh(*),vzh(*) ==>barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    vxj(*),vyj(*),vzj(*) ==> jacobi. particle velocities
*                                              (real array)
*       
*     ALGORITHM: See my notes Nov 21/92 
*     REMARKS:  
*       
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  Jan 29, 1993.
*     REVISIONS:  2/17/93 HFL

	subroutine coord_vh2vj(nbod,mass,vxh,vyh,vzh,vxj,vyj,vzj)


      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c...  Outputs:
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)

c...  Internals:
	real*8 eta(NPLMAX)
	real*8 sumvx,sumvy,sumvz	
	real*8 capvx,capvy,capvz
	integer n

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi velocities

	eta(1) = mass(1)
	do n = 2,nbod
	  eta(n) = eta(n-1) + mass(n)
	enddo

	vxj(1) =0.d0
	vyj(1) =0.d0
	vzj(1) =0.d0

	vxj(2) = vxh(2)
	vyj(2) = vyh(2)
	vzj(2) = vzh(2)

	sumvx = mass(2)*vxh(2)
	sumvy = mass(2)*vyh(2)
	sumvz = mass(2)*vzh(2)

	capvx = sumvx/eta(2)
	capvy = sumvy/eta(2)
	capvz = sumvz/eta(2)

	do n=3,nbod
	  vxj(n) = vxh(n) - capvx
	  vyj(n) = vyh(n) - capvy
	  vzj(n) = vzh(n) - capvz

	  if(n.lt.nbod) then
             sumvx = sumvx + mass(n)*vxh(n)
             sumvy = sumvy + mass(n)*vyh(n)
             sumvz = sumvz + mass(n)*vzh(n)

             capvx =sumvx/eta(n)
             capvy =sumvy/eta(n)
             capvz =sumvz/eta(n)
          endif

	enddo

123	return
	end    ! coord_vh2vj
c--------------------------------------------------------------------------
