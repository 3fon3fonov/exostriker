c*************************************************************************
c                            SYMBA5_MERGE.F
c*************************************************************************
c This subroutine checks to see if there are encounters
c
c             Input:
c                 t             ==>  current time (real scalar)
c                 dt            ==>  time step (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ip1,ip2       ==>  The two bodies to check (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 ireci         ==>  current recursion level (int scalar)
c                 irecl         ==>  maximum recursion level (int scalar)
c                 svdotrold     ==>  old radial velocity test
c                                   = .true. if i,j are receding
c                                   = .false is approaching
c                                     (logical*1 scalar)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c             mergelst          ==>  list of mergers (int array)
c             mergecnt          ==>  count of mergers (int array)
c             rhill             ==>  Hill sphere of planet (real Scalar)
c             eoff              ==>  Energy offset (real scalar)
c                ielc           ==>  number of encounters (integer*2 scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c
c             Output:  Changed only if a Megrer happens
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c             mergelst          ==>  list of mergers (int array)
c             mergecnt          ==>  count of mergers (int array)
c             rhill             ==>  Hill sphere of planet (real Scalar)
c             eoff              ==>  Energy offset (real scalar)
c                ielc           ==>  number of encounters (integer*2 scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c
c Remarks: 
c Authors:  Hal Levison
c Date:   1/2/97
c Last revision: 1/24/97
c

      subroutine symba5_merge(t,dt,nbod,ip1,ip2,mass,xh,yh,zh,vxb,
     &     vyb,vzb,ireci,irecl,svdotrold,iecnt,rpl,
     &     mergelst,mergecnt,rhill,eoff,ielc,ielst)


      include '../swift.inc'
      include 'symba5.inc'

c...  Inputs: 
      integer nbod,ireci,irecl,ip1,ip2
      real*8 t,dt
      logical*1 svdotrold

c...  Inputs and Outputs:
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod),eoff
      real*8 vxb(nbod),vyb(nbod),vzb(nbod),rpl(nbod),rhill(nbod)
      integer*2 iecnt(NTPMAX)
      integer mergelst(2,NTPMAX),mergecnt
      integer ip1l,ip2l
      integer*2 ielst(2,NENMAX),ielc

c...  Outputs

c...  Internals
      integer ialpha
      real*8 xr,yr,zr,vxr,vyr,vzr,vdotr,tcross2
      real*8 rlim,rlim2,rr2,massc,a,e,peri,dt2

c-----
c...  Executable code 

      xr = xh(ip2) - xh(ip1)
      yr = yh(ip2) - yh(ip1)
      zr = zh(ip2) - zh(ip1)
      rr2 = xr**2 + yr**2 + zr**2

      rlim = rpl(ip1)+rpl(ip2)

      if(rlim.eq.0.0d0) RETURN  ! <======  NOTE !!!!!

      rlim2 = rlim*rlim

      if(rlim2.ge.rr2) then
         ip1l = ip1
         ip2l = ip2 
         call discard_mass_merge5(t,nbod,ip1l,ip2l,mass,xh,yh,
     &        zh,vxb,vyb,vzb,rpl,eoff,ielc,ielst,NENMAX)
         mergecnt = mergecnt + 1
         mergelst(1,mergecnt) = ip1l
         mergelst(2,mergecnt) = ip2l
         rhill(ip2l) = 0.0d0
         call util_hills1(mass(1),mass(ip1l),xh(ip1l),yh(ip1l),
     &        zh(ip1l),vxb(ip1l),vyb(ip1l),vzb(ip1l),rhill(ip1l))
         return      !   <=== NOTE !!!!!!!!!
      endif

      vxr = vxb(ip2) - vxb(ip1)
      vyr = vyb(ip2) - vyb(ip1)
      vzr = vzb(ip2) - vzb(ip1)
      vdotr = xr*vxr + yr*vyr + zr*vzr

      if( svdotrold .and. (vdotr.gt.0.0d0)) then

         tcross2 = rr2/(vxr**2+vyr**2+vzr**2)
         dt2 = dt*dt

         if(tcross2.le.dt2) then
            massc = mass(ip1) + mass(ip2)
            call orbel_xv2aeq(xr,yr,zr,vxr,vyr,vzr,massc,
     &           ialpha,a,e,peri)
            if( peri.lt.rlim) then
               ip1l = ip1
               ip2l = ip2 
               call discard_mass_merge5(t,nbod,ip1l,ip2l,mass,xh,
     &              yh,zh,vxb,vyb,vzb,rpl,eoff,ielc,
     &              ielst,NENMAX)
               mergecnt = mergecnt + 1
               mergelst(1,mergecnt) = ip1l
               mergelst(2,mergecnt) = ip2l
               rhill(ip2l) = 0.0d0
               call util_hills1(mass(1),mass(ip1l),xh(ip1l),yh(ip1l),
     &              zh(ip1l),vxb(ip1l),vyb(ip1l),vzb(ip1l),rhill(ip1l))
            endif
         endif
      endif

      return
      end                       ! symba5_merge
c------------------------------------------------------

