c*************************************************************************
c                            DISCARD_HJS.F
c*************************************************************************
c This subroutine checks to see if a particle should be discarded because
c of its position or because it becomes unbound
c  (HJS, Generalized Jacobi coordinates only)
c             Input:
c                 time          ==>  current time (real scalar)
c                 dt            ==>  time step  (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 itp           ==>  # of tp (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xj,yj,zj      ==>  Planet positions in Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  Planet vels. in Jacobi coord 
c                                    (real arrays)
c                 xjt,yjt,zjt    ==> Particle position in Jacobi coord 
c                                      (real scalars)
c                 vxjt,vyjt,vzjt ==> Particle  velocity in Jacobi coord 
c                                        (real scalars)
c                 umatp          ==> Conversion vector Jacobi => Barycentric
c                                      (real array)
c                 etatpt         ==> Mass of centers fot tp
c                                        (real scalar)
c                 rmin,rmax      ==>  maximum and min distance from Sun
c                                     if <0  then don't check
c                                        (real scalar)
c                 rmaxu          ==>  maximum distance from Sun in not bound
c                                     if <0  then don't check
c                                        (real scalar)
c                  qmin          ==> Smallest perihelion distance
c                                      if <0  then don't check
c                                          (real scalar)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,1) closest approach to a planet
c
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c                                      istat(i,2) = -2  -> a<0 & r_cent>rmaxu
c	                               istat(i,2) = -3     r_cent>rmax
c	                               istat(i,2) = -4     q<qmin
c	                               istat(i,2) =  1     r_cent<rmin
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,2) closest approach to a planet
c                                          as determined by encounter routines.
c
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c	                               istat(i,2) =  n    too close to planet n
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,1) time of discard.
c
c Remarks: Adapted from discard_mjs.f
c Authors:  Hervé Beust
c Date:    Feb. 08, 2002

      subroutine discard_hjs(time,dt,nbod,itp,mass,xj,yj,zj,vxj,
     &     vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,umatp,etatpt,rmin,rmax,
     &     rmaxu,qmin,lclose,rplsq,istat,rstat)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,itp
      real*8 mass(nbod),xj(nbod),yj(nbod),zj(nbod),time,rplsq(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xjt,yjt,zjt,umatp(nbod)
      real*8 vxjt,vyjt,vzjt
      real*8 etatpt
      real*8 rmin,rmax,rmaxu,dt,qmin
      logical*2 lclose

c...  Input and Output
      integer istat(NSTAT)
      real*8 rstat(NSTATR)
c...  Internals
      integer ialpha,i
      real*8  a,e,q,gm
      real*8 xbt,ybt,zbt
      real*8 vxbt,vybt,vzbt

c-----
c...  Executable code 

c...  Compute barycentric coordinates
      call coord_g2b_tp(nbod,umatp,xj,yj,zj,
     &          vxj,vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,
     &          vzjt,xbt,ybt,zbt,vxbt,vybt,vzbt)

      if( (rmin.ge.0.0) .or. (rmax.ge.0.0) .or. (rmaxu.ge.0.0) ) then
         gm = 0.0d0
         do i = 1,nbod
           gm = gm+mass(i)
         end do
         call discard_center_hjs(time,1,itp,gm,xjt,yjt,zjt,
     &       xbt,ybt,zbt,vxbt,vybt,vzbt,rmin,rmax,rmaxu,istat,rstat)
      endif

      if (qmin.ge.0.0) then
         call discard_peri(time,1,itp,etatpt,xjt,yjt,zjt,
     &       vxjt,vyjt,vzjt,0.0d0,qmin,istat,rstat)
      endif

      return
      end    ! discard_hjs.f
c-----------------------------------------------------------------------


