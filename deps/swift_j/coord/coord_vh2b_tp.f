c***********************************************************************
c	                    COORD_VH2B_TP.F
c***********************************************************************
*     PURPOSE: Converts test part from Heliocentric to Barycentric coords.
*              Velocity only
*     ARGUMENTS:  Input is 
*                              ntp ==> number of test part (<= NTPMAX)
*                                              (integer)
*		     vxht,vyht,vzht ==> heliocentric particle velocities
*                                             (real array)
*		     vxhs,vyhs,vzhs ==> bary vel of the Sun
*                                          (real scalar)
*                 Returned are
*                    vxbt,vybt,vzbt ==> bary. particle velocities
*                                            (real array)
*       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  11/14/96
*     REVISIONS:

	subroutine coord_vh2b_tp(ntp,vxht,vyht,vzht,vxs,vys,vzs,
     &      vxbt,vybt,vzbt)

      include '../swift.inc'

c...  Inputs: 
	integer ntp
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
        real*8 vxs,vys,vzs

c...  Outputs:
	real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)

c...  Internals:
	integer i

c----
c...  Executable code 
	do i=1,ntp
	  vxbt(i) = vxht(i) + vxs
	  vybt(i) = vyht(i) + vys
	  vzbt(i) = vzht(i) + vzs
	enddo

	return
	end     ! coord_vh2b_tp
c--------------------------------------------------------------------------

