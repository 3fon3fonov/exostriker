c***********************************************************************
c	                    COORD_G2B.F
c***********************************************************************
*     PURPOSE: Converts from GEneralized Jacobi to Barycentric coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*                    umat(NPLMAX,NPLMAX) ==> Passage matrix (2D real array)
*	             mass(*) ==>  masses (real array)
*		     xj(*),yj(*),zj(*) ==> Gen. Jacobi particle coords
*                                          (real array)
*		     vxj(*),vyj(*),vzj(*) ==> Gen. Jacobi particle velocities
*                                             (real array)
*                 Returned are
*                    xb(*),yb(*),zb(*) ==>  Barycentric particle positions
*                                          (real array)
*                    vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                            (real array)
*     Authors:  Hervé Beust
*     WRITTEN:  Jan. 24, 2002

	subroutine coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xb,yb,zb,vxb,vyb,vzb)

      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(nbod),umat(NPLMAX,NPLMAX)
	real*8 xj(nbod),yj(nbod),zj(nbod)
	real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Outputs:
	real*8 xb(nbod),yb(nbod),zb(nbod)
	real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Internals:
	integer i,j
c----
c...  Executable code 

c...  Matrix product      XB=UMAT*XJ

        do i = 1,nbod
          xb(i) = 0.0d0
          yb(i) = 0.0d0
          zb(i) = 0.0d0
          vxb(i) = 0.0d0
          vyb(i) = 0.0d0
          vzb(i) = 0.0d0
          do j = 2,nbod   ! 2 because first jac. coord is zero
            if (umat(i,j).ne.0.0d0) then
              xb(i) = xb(i) + umat(i,j)*xj(j)
              yb(i) = yb(i) + umat(i,j)*yj(j)
              zb(i) = zb(i) + umat(i,j)*zj(j)
              vxb(i) = vxb(i) + umat(i,j)*vxj(j)
              vyb(i) = vyb(i) + umat(i,j)*vyj(j)
              vzb(i) = vzb(i) + umat(i,j)*vzj(j)
            end if
          end do
        end do

	return
	end    ! coord_g2b

c--------------------------------------------------------------------------
