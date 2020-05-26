c***********************************************************************
c	                    COORD_VB2VG_TP.F
c***********************************************************************
*     PURPOSE: Converts BARY coords to Generalized Jacobi for ONE tp.
*                (vels. only)
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*                    umatp  ==> The reverse conversion vector for tp
*		     vxb(*),vyb(*),vzb(*) ==> BARY particle velocities
*                                             (real array)
*                    vzbt,vybt,vzbt ==> Bary pos. and vel.
*                                                of the tp.
*                 Returned are
*                    vzjt,vyjt,vzjt ==> Gen. Jacobi  pos. and vel.
*                                                of the tp.
*       
*     Authors:  Herve Beust
*     WRITTEN:  Feb. 07, 2002
*     Adapted from coord_b2g.f

	subroutine coord_vb2vg_tp(nbod,matp,vxb,vyb,vzb,
     &      vxbt,vybt,vzbt,vxjt,vyjt,vzjt)


        include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 matp(nbod)
	real*8 vxb(nbod),vyb(nbod),vzb(nbod)
        real*8 vxbt,vybt,vzbt

c...  Outputs:
	real*8 vxjt,vyjt,vzjt

c...  Internals:
	integer j
	real*8 sumvx,sumvy,sumvz

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi coords

	sumvx = 0.0d0
	sumvy = 0.0d0
	sumvz = 0.0d0

        do j = 1,nbod 
          if (matp(j).ne.0.0d0) then
            sumvx = sumvx + matp(j)*vxb(j)
            sumvy = sumvy + matp(j)*vyb(j)
            sumvz = sumvz + matp(j)*vzb(j)
          end if
        end do

	vxjt = vxbt + sumvx
	vyjt = vybt + sumvy
	vzjt = vzbt + sumvz

123	return
	end     ! coord_vb2vg_tp

c--------------------------------------------------------------------------

