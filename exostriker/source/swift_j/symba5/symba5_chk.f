c*************************************************************************
c                            SYMBA5_CHK.F
c*************************************************************************
c This subroutine checks to see if there are encounters
c
c             Input:
c                 rhill         ==>  Radius of hill sphere (real array)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ip1,ip2       ==>  The two bodies to check (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xr,yr,zr      ==>  relative position of the two bodies
c                                    (real scalar)
c                 vxr,vyr,vzr   ==>  relative velocity of the two bodies
c                                    (real scalar)
c                 dt            ==>  time step  (real scalar)
c                 irec          ==>  current recursion level (int scalar)
c             Output:
c                 icflg         ==> ecounter?  = 1 Yes
c                                              =  0 No (integer scalar)  
c                 svdotr        ==> = .true. if i,j are receding
c                                   = .false is approaching
c                                     (logical*1 scalar)
c
c Remarks: Based on plh_chk.f.  Same as symba_chk.f
c Authors:  Hal Levison
c Date:   3/20/97
c Last revision: 5/28/99
c

      subroutine symba5_chk(rhill,nbod,ip1,ip2,mass,xr,yr,zr,
     &     vxr,vyr,vzr,dt,irec,icflg,svdotr)

      include '../swift.inc'
      include 'symba5.inc'

c...  Inputs: 
      integer nbod,irec,ip1,ip2
      real*8 mass(nbod),rhill(nbod)
      real*8 xr,yr,zr,vxr,vyr,vzr,dt

c...  Outputs
      integer icflg
      logical*1 svdotr

c...  Internals
      real*8 r2crit,r2critp,rcrit
      real*8 vdotr

c-----
c...  Executable code 

      rcrit = (rhill(ip1)+rhill(ip2)) * RHSCALE * (RSHELL**(irec))
      r2crit = rcrit*rcrit
      r2critp = -1.0d0          ! not used here
      
      call rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,
     &     r2crit,r2critp,icflg)
      
      vdotr = xr*vxr + yr*vyr + zr*vzr
      svdotr = (vdotr.lt.0.0d0)

      return
      end                       ! symba5_chk
c------------------------------------------------------
