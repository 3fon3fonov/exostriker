c***********************************************************************
c	                    COORD_G2B_TP.F
c***********************************************************************
*     PURPOSE: Converts from Generalized Jacobi to BARY coords for ONE tp.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*                    umatp  ==> The reverse conversion vector for tp
*		     xj(*),yj(*),zj(*) ==> Jacobi particle coords
*                                          (real array)
*		     vxj(*),vyj(*),vzj(*) ==> Jacobi particle velocities
*                                             (real array)
*                    xjt,yjt,zjt,vzjt,vyjt,vzjt ==> Gen. Jacobi  pos. and vel.
*                                                of the tp.
*                 Returned are
*                    xbt,ybt,zbt,vzbt,vybt,vzbt ==> Bary pos. and vel.
*                                                of the tp.
*       
*     Authors:  Herve Beust
*     WRITTEN:  Jan. 24, 2002
*     Adapted from coord_j2h.f

	subroutine coord_g2b_tp(nbod,umatp,xj,yj,zj,vxj,vyj,vzj,
     &      xjt,yjt,zjt,vxjt,vyjt,vzjt,xbt,ybt,zbt,vxbt,vybt,vzbt)


        include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 umatp(nbod)
	real*8 xj(nbod),yj(nbod),zj(nbod)
	real*8 vxj(nbod),vyj(nbod),vzj(nbod)
        real*8 xjt,yjt,zjt,vxjt,vyjt,vzjt

c...  Outputs:
	real*8 xbt,ybt,zbt,vxbt,vybt,vzbt

c...  Internals:
	integer j
	real*8 sumx,sumy,sumz,sumvx,sumvy,sumvz

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi coords

	sumx = 0.0d0
	sumy = 0.0d0
	sumz = 0.0d0
	sumvx = 0.0d0
	sumvy = 0.0d0
	sumvz = 0.0d0

        do j = 2,nbod   ! 2 because first jac. coord is zero
          if (umatp(j).ne.0.0d0) then
            sumx = sumx + umatp(j)*xj(j)
            sumy = sumy + umatp(j)*yj(j)
            sumz = sumz + umatp(j)*zj(j)
            sumvx = sumvx + umatp(j)*vxj(j)
            sumvy = sumvy + umatp(j)*vyj(j)
            sumvz = sumvz + umatp(j)*vzj(j)
          end if
        end do

	xbt = xjt + sumx
	ybt = yjt + sumy
	zbt = zjt + sumz
	vxbt = vxjt + sumvx
	vybt = vyjt + sumvy
	vzbt = vzjt + sumvz

123	return
	end     ! coord_g2b_tp

c--------------------------------------------------------------------------

