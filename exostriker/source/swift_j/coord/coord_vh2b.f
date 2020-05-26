c***********************************************************************
c	                    COORD_VH2B.F
c***********************************************************************
*     PURPOSE: Converts from Heliocentric to Barycentric coords. 
*              Velocity only
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
*                                             (real array)
*                 Returned are
*                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
*                                            (real array)
*                    msys              ==>  Total mass of of system
*                                            (real scalar)       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  11/14/96
*     REVISIONS: 11/21/96

	subroutine coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)

      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(nbod)
	real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Outputs:
	real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Internals:
	real*8 msys,vxtmp,vytmp,vztmp
	integer n

c----
c...  Executable code 

	msys = mass(1)
	vxtmp =0.d0
	vytmp =0.d0
	vztmp =0.d0

	do n=2,nbod
	   msys = msys +mass(n)
	   vxtmp = vxtmp + mass(n)*vxh(n)
	   vytmp = vytmp + mass(n)*vyh(n)
	   vztmp = vztmp + mass(n)*vzh(n)
	enddo

	vxb(1) = -vxtmp/msys
	vyb(1) = -vytmp/msys
	vzb(1) = -vztmp/msys

	do n=2,nbod
	  vxb(n) = vxh(n) + vxb(1)
	  vyb(n) = vyh(n) + vyb(1)
	  vzb(n) = vzh(n) + vzb(1)
	enddo

	return
	end     ! coord_vh2b
c--------------------------------------------------------------------------

