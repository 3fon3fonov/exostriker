c***********************************************************************
c	                    COORD_H2B_TP.F
c***********************************************************************
*     PURPOSE: Converts test part from Heliocentric to Barycentric coords.
*     ARGUMENTS:  Input is 
*                              ntp ==> number of test part (<= NTPMAX)
*                                              (integer)
*		        xht,yht,zht ==> heliocentric particle coords
*                                          (real array)
*		     vxht,vyht,vzht ==> heliocentric particle velocities
*                                             (real array)
*		        xhs,yhs,zhs ==> bary coords of the Sun
*                                          (real scalar)
*		     vxhs,vyhs,vzhs ==> bary vel of the Sun
*                                          (real scalar)
*                 Returned are
*                       xbt,ybt,zbt ==> bary. particle positions
*                                            (real array)
*                    vxbt,vybt,vzbt ==> bary. particle velocities
*                                            (real array)
*       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  2/18/92
*     REVISIONS:

	subroutine coord_h2b_tp(ntp,xht,yht,zht,vxht,vyht,vzht,
     &      xs,ys,zs,vxs,vys,vzs,
     &      xbt,ybt,zbt,vxbt,vybt,vzbt)

      include '../swift.inc'

c...  Inputs: 
	integer ntp
	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
        real*8 xs,ys,zs,vxs,vys,vzs

c...  Outputs:
	real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)
	real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)

c...  Internals:
	integer i

c----
c...  Executable code 
	do i=1,ntp
	  xbt(i) = xht(i) + xs
	  ybt(i) = yht(i) + ys
	  zbt(i) = zht(i) + zs
	  vxbt(i) = vxht(i) + vxs
	  vybt(i) = vyht(i) + vys
	  vzbt(i) = vzht(i) + vzs
	enddo

	return
	end     ! coord_h2b_tp
c--------------------------------------------------------------------------

