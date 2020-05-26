c*************************************************************************
c                            STEP_KDK_PL_HJS.F
c*************************************************************************
c This subroutine takes a step in generalized Jacobi coord (HJS)  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 umat,mat      ==>  Conversion matrixes (Jacobi - Barycentric)
c                                     (2D real arrays)
c                 oloc          ==>  Link Bodies <--> Orbits (2D int array)
c                 eta,mu        ==>  Masses for centers & sats (real arrays)
c                 xj,yj,zj      ==>  initial position in Gen. Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial velocity in Gen. Jacobi coord 
c                                    (real arrays)
c                 dt            ==>  time step
c            Input & Output
c                 xbbeg,ybbeg,zbbeg   ==>  initial position in bary coord 
c                                    (real arrays)
c                 axbbeg,aybbeg,azbbeg   ==>  initial accel. in bary coord 
c                                    (real arrays)
c                 xj,yj,zj      ==>  final position in Gen. Jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in Gen. Jacobi coord 
c                                       (real arrays)
c                 ir3jbeg,ir3jend ==> Inverse radii^3 beg & end (real arrays)
c                 xbend,ybend,zbend ==> Final position in bary coord 
c                                    (real arrays)
c                 axbend,aybend,azbend ==> Final accel in bary coord 
c                                    (real arrays)
c                 vxjh,vyjh,vzjh==>  Middle velocity in Gen. Jacobi coord 
c                                       (real arrays)
c
c Remarks: Adapted from step_kdk_pl.f
c Authors:  Herve Beust
c Date:    Feb 11, 2002
c Last revision: Dec 10, 2002

      subroutine step_kdk_pl_hjs(i1st,nbod,oloc,umat,mat,mass,eta,mu,
     &     xj,yj,zj,vxj,vyj,vzj,xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,
     &     azbbeg,ir3jbeg,vxjh,vyjh,vzjh,xbend,ybend,zbend,
     &     axbend,aybend,azbend,ir3jend,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod),dt,eta(nbod),mu(nbod)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 vxjh(nbod),vyjh(nbod),vzjh(nbod)
      real*8 xbbeg(nbod),ybbeg(nbod),zbbeg(nbod)
      real*8 axbbeg(nbod),aybbeg(nbod),azbbeg(nbod)
      real*8 ir3jbeg(nbod),ir3jend(nbod)
      real*8 xbend(nbod),ybend(nbod),zbend(nbod)
      real*8 axbend(nbod),aybend(nbod),azbend(nbod)

c...  Internals:
      integer i
      real*8 dth 
      real*8 axj(NPLMAX),ayj(NPLMAX),azj(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

      save axj,ayj,azj     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if (i1st.eq.0) then
c...     Compute barycentric coords
        call  coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xbbeg,ybbeg,zbbeg,vxb,vyb,vzb)
c...     Get the accelerations in bary frame. if frist time step
c        call getacch_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
c     &     xbbeg,ybbeg,zbbeg,ir3jbeg,axbbeg,aybbeg,azbbeg)
c...     Convert bary accels to Jacobi accels (use vel. procedure)
c        call coord_vb2vg(nbod,mat,mass,axbbeg,aybbeg,azbbeg,
c     &                                               axj,ayj,azj)
c...     Get the Jacobi accels if first time step
         call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbbeg,ybbeg,zbbeg,ir3jbeg,axj,ayj,azj)
         call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axbbeg,aybbeg,azbbeg)

        i1st = 1    ! turn this off
      endif
c...       Apply a Jacobi kick for a half dt 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dth)

c...   Drift in Jacobi coords for the full step 
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dt)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                   xbend,ybend,zbend,vxb,vyb,vzb)
c...  Save Jacobi velocities at middle point (for tp's after)
      do i=1,nbod
        vxjh(i) = vxj(i)
        vyjh(i) = vyj(i)
        vzjh(i) = vzj(i)
      end do
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbend,ybend,zbend,ir3jend,axj,ayj,azj)
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axbend,aybend,azbend)

c...  Compute now accelerations in barycentric coordinates
c      call getacch_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
c     &           xbend,ybend,zbend,ir3jend,axbend,aybend,azbend)
c...  Convert bary accels to Jacobi accels (use vel. procedure)
c      call coord_vb2vg(nbod,mat,mass,axbend,aybend,azbend,
c     &                                               axj,ayj,azj)
c...  Apply a Jacobi kick for a half dt 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dth)

c... Compute barycentric coordinates
c      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
c     &             xb,yb,zb,vxb,vyb,vzb)

      return
      end   ! step_kdk_pl_hjs
c---------------------------------------------------------------------

