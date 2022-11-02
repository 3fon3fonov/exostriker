c***********************************************************************
c	                    COORD_VB2H.F
c***********************************************************************
*     PURPOSE: Converts from Barycentric to Helio coords.
*               Velocity only
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
*                                 SYMMETRY IN SUBROUTINE CALLS
*		     vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
*                                            (real array)
*       
*     ALGORITHM: Obvious 
*     Authors:  Hal Levison
*     WRITTEN:  11/14/96
*     REVISIONS: 11/21/96

      subroutine coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod)
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Outputs:
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Internals:
      integer i

c----
c...  Executable code 

	vxb(1) = - mass(2)*vxb(2)
	vyb(1) = - mass(2)*vyb(2)
	vzb(1) = - mass(2)*vzb(2)

	do i=3,nbod
	   vxb(1) = vxb(1) - mass(i)*vxb(i)
	   vyb(1) = vyb(1) - mass(i)*vyb(i)
	   vzb(1) = vzb(1) - mass(i)*vzb(i)
	enddo

	vxb(1) = vxb(1)/mass(1)
	vyb(1) = vyb(1)/mass(1)
	vzb(1) = vzb(1)/mass(1)

	do i=2,nbod
	   vxh(i) = vxb(i)  - vxb(1)  
	   vyh(i) = vyb(i)  - vyb(1)  
	   vzh(i) = vzb(i)  - vzb(1)  
	enddo

	return
	end     ! coord_vb2h

c--------------------------------------------------------------------------

