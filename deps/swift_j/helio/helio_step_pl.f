c*************************************************************************
c                            HELIO_STEP_PL.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xbeg,ybeg,zbeg ==>  position at beginning of drift
c                                       (real arrays)
c                 xend,yend,zend ==>  position at end of drift
c                                       (real arrays)
c              vxbeg,vybeg,vzbeg ==>  velocity at beginning of drift
c                                       (real arrays)
c              vxend,vyend,vzend ==>  velocity at end of drift
c                                       (real arrays)
c               ptxb,ptyb,ptzb  ==> beg momentum of sun: tp's need this   
c                                       (real scalars)
c               ptxe,ptye,ptze  ==> end momentum of sun: tp's need this   
c                                       (real scalars)
c              vxsb,vxsb,vxsb  ==> Initial vel of the Sun: tp's need this
c                                       (real scalars)
c              vxse,vxse,vxse  ==> final vel of the Sun: tp's need this
c                                       (real scalars)
c
c Remarks: Based on step_kdk_pl.f 
c Authors:  Hal Levison 
c Date:    10/14/96
c Last revision: 

      subroutine helio_step_pl(i1st,nbod,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,dt,xbeg,ybeg,zbeg,     
     &     xend,yend,zend,vxbeg,vybeg,vzbeg,     
     &     vxend,vyend,vzend,ptxb,ptyb,ptzb,ptxe,ptye,
     &     ptze,vxsb,vysb,vzsb,vxse,vyse,vzse)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Outputs Only: 
      real*8 ptxb,ptyb,ptzb
      real*8 ptxe,ptye,ptze
      real*8 vxsb,vysb,vzsb,vxse,vyse,vzse
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)

c...  Internals:
      integer i
      real*8 dth 
      real*8 axh(NPLMAX),ayh(NPLMAX),azh(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX),msys

      save axh,ayh,azh     ! Note this !!
      save vxb,vyb,vzb     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then
c...      Convert vel to bery to jacobi coords
          call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)
         i1st = 1    ! turn this off
      endif

      vxsb = vxb(1)
      vysb = vyb(1)
      vzsb = vzb(1)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxb,ptyb,ptzb)

c...  Get the accelerations in helio frame. if frist time step
      call helio_getacch(nbod,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  remember the current position of the planets
      do i=1,nbod
         xbeg(i) = xh(i)
         ybeg(i) = yh(i)
         zbeg(i) = zh(i)
         vxbeg(i) = vxb(i)
         vybeg(i) = vyb(i)
         vzbeg(i) = vzb(i)
      enddo

c..   Drift in Jacobi coords for the full step 
      call helio_drift(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)

c...  now remember these positions
      do i=1,nbod
         xend(i) = xh(i)
         yend(i) = yh(i)
         zend(i) = zh(i)
         vxend(i) = vxb(i)
         vyend(i) = vyb(i)
         vzend(i) = vzb(i)
      enddo

c...  Get the accelerations in helio frame. if frist time step
      call helio_getacch(nbod,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxe,ptye,ptze)

c...  convert back to helio velocities
      call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)

      vxse = vxb(1)
      vyse = vyb(1)
      vzse = vzb(1)

      return
      end   ! helio_step_pl
c---------------------------------------------------------------------

