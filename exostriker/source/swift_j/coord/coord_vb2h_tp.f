c***********************************************************************
c	                    COORD_VB2H_TP.F
c***********************************************************************
*     PURPOSE: Converts test part from Barycentric to Heliocentric coords.
*              Velocity only
*     ARGUMENTS:  Input is 
*                              ntp ==> number of test part (<= NTPMAX)
*                                              (integer)
*                             istat ==>  Status flag
*                    vxbt,vybt,vzbt ==> bary. particle velocities
*                                            (real array)
*		     vxhs,vyhs,vzhs ==> bary vel of the Sun
*                                          (real scalar)
*                 Returned are
*		     vxht,vyht,vzht ==> heliocentric particle velocities
*                                             (real array)
*       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  11/14/96
*     REVISIONS: 11/15/96

	subroutine coord_vb2h_tp(ntp,istat,vxbt,vybt,vzbt,vxs,
     &      vys,vzs,vxht,vyht,vzht)


      include '../swift.inc'

c...  Inputs: 
	integer ntp
	real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
        real*8 vxs,vys,vzs
        integer istat(ntp)

c...  Outputs:
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

c...  Internals:
	integer i

c----
c...  Executable code 
	do i=1,ntp
           if(istat(i).eq.0) then
              vxht(i) = vxbt(i) - vxs
              vyht(i) = vybt(i) - vys
              vzht(i) = vzbt(i) - vzs
           endif
        enddo

	return
	end     ! coord_vb2h_tp
c--------------------------------------------------------------------------

