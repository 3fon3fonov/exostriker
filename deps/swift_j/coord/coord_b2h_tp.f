c***********************************************************************
c	                    COORD_B2H_TP.F
c***********************************************************************
*     PURPOSE: Converts test part from Barycentric to Heliocentric coords.
*     ARGUMENTS:  Input is 
*                              ntp ==> number of test part (<= NTPMAX)
*                                              (integer)
*                       xbt,ybt,zbt ==> bary. particle positions
*                                            (real array)
*                    vxbt,vybt,vzbt ==> bary. particle velocities
*                                            (real array)
*		        xhs,yhs,zhs ==> bary coords of the Sun
*                                          (real scalar)
*		     vxhs,vyhs,vzhs ==> bary vel of the Sun
*                                          (real scalar)
*                 Returned are
*		        xht,yht,zht ==> heliocentric particle coords
*                                          (real array)
*		     vxht,vyht,vzht ==> heliocentric particle velocities
*                                             (real array)
*       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  2/18/92
*     REVISIONS:

	subroutine coord_b2h_tp(ntp,xbt,ybt,zbt,vxbt,vybt,vzbt,
     &      xs,ys,zs,vxs,vys,vzs,
     &      xht,yht,zht,vxht,vyht,vzht)


      include '../swift.inc'

c...  Inputs: 
	integer ntp
	real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)
	real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
        real*8 xs,ys,zs,vxs,vys,vzs

c...  Outputs:
	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

c...  Internals:
	integer i

c----
c...  Executable code 
	do i=1,ntp
	  xht(i) = xbt(i) - xs
	  yht(i) = ybt(i) - ys
	  zht(i) = zbt(i) - zs
	  vxht(i) = vxbt(i) - vxs
	  vyht(i) = vybt(i) - vys
	  vzht(i) = vzbt(i) - vzs
	enddo

	return
	end     ! coord_b2h_tp
c--------------------------------------------------------------------------

