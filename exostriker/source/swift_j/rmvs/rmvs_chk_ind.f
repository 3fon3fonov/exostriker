c*************************************************************************
c                            RMVS_CHK_IND.F
c*************************************************************************
c  Subroutine to check if a test particle and planet
c  are having or **will** have an encounter 
c  in the next timestep. 
c
c             Input:
c                 xr,yr,zr     ==>  relative position of tp wrt planet
c                                   (real scalar)
c                 vxr,vyr,vzr  ==>  relative velocity of tp wrt planet
c                                   (real scalor)
c                 dt           ==>  time step (real scalor)
c                 r2crit       ==> boundary of outer enc region
c                                   (real scalor)
c                 r2critp      ==> boundary of inner (planocentric) enc region
c                                   (real scalor)
c             Output:
c                 iflag        ==> encounter?  =  0 no
c                                              =  1 yes, in outer region
c                                              = -1 yes, in inner region
c
c
c Remarks: Based on Hal's wiscl_fk.f' but origonaly written by Martin Duncan
c Authors:  Hal Levison 
c Date:    2/19/93
c Last revision: 

      subroutine rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,
     &                         r2crit,r2critp,iflag)


      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs: 
      real*8 xr,yr,zr,vxr,vyr,vzr,dt,r2crit,r2critp

c...  Outputs
	integer iflag

c...  Internals
	real*8 r2,v2,vdotr,tmin,r2min

c-----
c...  Executable code 

c...    First check if we're already in the encounter region. If so return
c.             with flag set to one.
	r2 = xr**2 + yr**2 + zr**2
	if (r2 .le. r2critp) then
	   iflag = -1
	   return
	endif

c...    If we're heading outward, use r2 to calc iflag
	vdotr = xr*vxr + yr*vyr + zr*vzr
	if (vdotr . gt. 0.d0) then
           if (r2 .ge. r2crit) then
	      iflag = 0
           else
	      iflag = 1
           endif
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

        if (r2min.le.r2critp)then
	   iflag = -1
	else if(r2min.le.r2crit)then
	   iflag = 1
	else 
	   iflag = 0
	endif

	return
	end  ! rmvs_chk_ind
c--------------------------------------------------------------


