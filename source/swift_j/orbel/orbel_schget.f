***********************************************************************
c	                  ORBEL_SCHGET.F
***********************************************************************
*     PURPOSE:  Given an angle, efficiently compute sinh and cosh.
*
*        Input:
*             angle ==> angle in radians (real scalar)
*        
*        Output:
*             shx    ==>  sinh(angle)  (real scalar)
*             chx    ==>  cosh(angle)  (real scalar)
*
*     ALGORITHM: Obvious from the code 
*     REMARKS: Based on the routine SCGET for sine's and cosine's.
*       We use the sqrt rather than cosh (it's faster)
*       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
*       OR OVERFLOWS WILL OCCUR!
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 6, 1992.
*     REVISIONS: 
***********************************************************************

	subroutine orbel_schget(angle,shx,chx)

      include '../swift.inc'

c...  Inputs Only: 
        real*8 angle

c...  Output:
	real*8 shx,chx

c----
c...  Executable code 

	shx = sinh(angle)
	chx= sqrt(1.d0 + shx*shx)

	return
	end   ! orbel_schget
c---------------------------------------------------------------------
