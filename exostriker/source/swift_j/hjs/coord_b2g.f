c***********************************************************************
c	                    COORD_B2G.F
c***********************************************************************
*     PURPOSE: Converts from Ganeralized Jacobi to Barycentric coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*                    umat(NPLMAX,NPLMAX) ==> Passage matrix (2D real array)
*	             mass(*) ==>  masses (real array)
*                    xb(*),yb(*),zb(*) ==>  Barycentric particle positions
*                                          (real array)
*                    vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                            (real array)
*                 Returned are
*		     xj(*),yj(*),zj(*) ==> Gen. Jacobi particle coords
*                                          (real array)
*		     vxj(*),vyj(*),vzj(*) ==> Gen. Jacobi particle velocities
*                                             (real array)
*     Authors:  Hervé Beust
*     WRITTEN:  Feb. 11, 2002

	subroutine coord_b2g(nbod,mat,mass,xb,yb,zb,vxb,vyb,vzb,
     &           xj,yj,zj,vxj,vyj,vzj)

        include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(nbod),mat(NPLMAX,NPLMAX)
	real*8 xb(nbod),yb(nbod),zb(nbod)
	real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Outputs:
	real*8 xj(nbod),yj(nbod),zj(nbod)
	real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Internals:
	integer i,j
c----
c...  Executable code 

c...  Matrix product      XJ=MAT*XB

	xj(1) = 0.d0
	yj(1) = 0.d0
	zj(1) = 0.d0
	vxj(1) = 0.d0
	vyj(1) = 0.d0
	vzj(1) = 0.d0

c... The first Generalized Jacobi coordinates are necessarily zero 
        do i = 2,nbod
          xj(i) = 0.0d0
          yj(i) = 0.0d0
          zj(i) = 0.0d0
          vxj(i) = 0.0d0
          vyj(i) = 0.0d0
          vzj(i) = 0.0d0
          do j = 1,nbod
            if (mat(i,j).ne.0.0d0) then
              xj(i) = xj(i) + mat(i,j)*xb(j)
              yj(i) = yj(i) + mat(i,j)*yb(j)
              zj(i) = zj(i) + mat(i,j)*zb(j)
              vxj(i) = vxj(i) + mat(i,j)*vxb(j)
              vyj(i) = vyj(i) + mat(i,j)*vyb(j)
              vzj(i) = vzj(i) + mat(i,j)*vzb(j)
            end if
          end do
        end do

	return
	end    ! coord_b2g

c--------------------------------------------------------------------------
