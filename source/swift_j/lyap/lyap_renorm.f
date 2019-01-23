c******************************************************************
c                            LYAP_RENORM.F
c*****************************************************************
c   Computes the distance and adds to the cum. sum for computing
c   Lyap exponent for a tp and its shadow. Moves the shadow back to
c   its original distance in phase space along the separation vector.
c             Input:
c              xt,yt,zt      ==>  current position of a TP in Helio coord 
c                                    (real scalars)
c              vxt,vyt,vzt   ==>  current position of TP in Helio coord 
c                                    (real scalars)
c              xs,ys,zs      ==>  current position of shawdow part in Helio coord 
c                                    (real scalars)
c              vxs,vys,vzs   ==>  current position of SP in Helio coord 
c                                    (real scalars)
c                  distorig  ==>  original distance between TP and SP
c                  logsum    ==>  current values of the log of gamma
c                                    (real scalar)
c
c             Output:
c              xs,ys,zs      ==>  renormalized position of SP in Helio coord 
c                                    (real scalars)
c              vxs,vys,vzs   ==>  renormalized position of SP in Helio coord 
c                                    (real scalars)
c               logsum       ==>  renormalized values of the log of gamma
c                                    (real scalar)
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    5/5/93 
c Last revision:  5/11/93  HFL

	subroutine lyap_renorm(xt,yt,zt,vxt,vyt,vzt,xs,ys,zs,
     &                    vxs,vys,vzs,distorig,logsum)

	include '../swift.inc'

c...    Input Only
        real*8 xt,yt,zt,vxt,vyt,vzt
	real*8 distorig

c....   Input and Output
	real*8 xs,ys,zs,vxs,vys,vzs
	real*8 logsum

c...    Internal
	real*8 dx,dy,dz,dvx,dvy,dvz,dist,fac

c-----
c...  Executable code


	dx = xt - xs
	dy = yt - ys
	dz = zt - zs
	dvx = vxt - vxs
	dvy = vyt - vys
	dvz = vzt - vzs

	dist = sqrt(dx**2 + dy**2 + dz**2 + dvx**2 + dvy**2 + dvz**2)
	fac = distorig/dist

	xs = xt + fac*dx
	ys = yt + fac*dy
	zs = zt + fac*dz
	vxs = vxt + fac*dvx
	vys = vyt + fac*dvy
	vzs = vzt + fac*dvz

	logsum = logsum + log(1.d0/fac)
	return
	end ! lyap_renorm
c---------------------------------------------------------------------------
