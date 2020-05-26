c*************************************************************************
c                            DISCARD_PL_CLOSE.F
c*************************************************************************
c  Subroutine to check if a test particle and planet
c  are having or **will** have a encounter w/ r^2<r2crit within
c  in the next timestep.  Uses LINEAR exploitation.
c
c             Input:
c                 xr,yr,zr     ==>  relative position of tp wrt planet
c                                   (real scalar)
c                 vxr,vyr,vzr  ==>  relative velocity of tp wrt planet
c                                   (real scalor)
c                 dt           ==>  time step (real scalor)
c                 r2crit       ==> boundary of encounter region
c                                   (real scalor)
c             Output:
c                 iflg        ==> encounter?  =  0 no
c                                              =  1 yes
c                 r2min       ==> square of smallest predicted distance
c
c Remarks: Based on rmvs_chk_ind.f
c Authors:  Hal Levison 
c Date:    2/21/94
c Last revision: 7/14/94

      subroutine discard_pl_close(xr,yr,zr,vxr,vyr,vzr,dt,r2crit,iflg,
     &     r2min)


      include '../swift.inc'

c...  Inputs: 
      real*8 xr,yr,zr,vxr,vyr,vzr,dt,r2crit

c...  Outputs
	integer iflg
        real*8 r2min

c...  Internals
	real*8 r2,v2,vdotr,tmin

c-----
c...  Executable code 

c...    First check if we're already in the encounter region. If so return
c.             with flag set to one.
	r2 = xr**2 + yr**2 + zr**2
	if (r2 .le. r2crit) then
	   iflg = 1
	   return
	endif

c...    If we're heading outward, then we are done
	vdotr = xr*vxr + yr*vyr + zr*vzr
	if (vdotr . gt. 0.d0) then
           iflg = 0
           return
	endif	

c...    We're not yet inside and are converging so we need to calc. the
c.           minimum separation attained in time dt.
	v2 = vxr**2 + vyr**2 + vzr**2
	tmin = -vdotr/v2

	if(tmin.lt.dt) then
	   r2min = r2 -(vdotr**2)/v2
	else
	   r2min = r2 + 2.d0*vdotr*dt + v2*dt*dt
	endif

        r2min = dmin1(r2min,r2)     ! really make sure

        if (r2min.le.r2crit)then
	   iflg = 1
	else 
	   iflg = 0
	endif

	return
	end  ! discard_pl_close
c--------------------------------------------------------------

