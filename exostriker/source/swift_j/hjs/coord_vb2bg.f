c***********************************************************************
c	                    COORD_VB2VG.F
c***********************************************************************
*     PURPOSE: Converts from Ganeralized Jacobi to Barycentric coords.
*                 (velocities only, but works for accels. also)
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*                    umat(*,*) ==> Passage matrix (inverse)
*                    vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                            (real array)
*                 Returned are
*		     vxj(*),vyj(*),vzj(*) ==> Gen. Jacobi particle velocities
*                                             (real array)
*     Authors:  Hervé Beust
*     WRITTEN:  Jan 24, 2002

        subroutine coord_vb2vg(nbod,mat,mass,vxb,vyb,vzb,vxj,vyj,vzj)

        include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(nbod),mat(NPLMAX,NPLMAX)
	real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Outputs:
	real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Internals:
	integer i,j
c----
c...  Executable code 

c...  Matrix product      VJ=MAT*VB

	vxj(1) = 0.d0
	vyj(1) = 0.d0
	vzj(1) = 0.d0

c... The first Jacobi coordinates are necessarily zero 
        do i = 2,nbod
          vxj(i) = 0.0d0
          vyj(i) = 0.0d0
          vzj(i) = 0.0d0
          do j = 1,nbod
            if (mat(i,j).ne.0.0d0) then
              vxj(i) = vxj(i) + mat(i,j)*vxb(j)
              vyj(i) = vyj(i) + mat(i,j)*vyb(j)
              vzj(i) = vzj(i) + mat(i,j)*vzb(j)
            end if
          end do
        end do

	return
	end    ! coord_vb2vg

c--------------------------------------------------------------------------
