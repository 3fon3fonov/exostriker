c*************************************************************************
c                            SPIN_PLANETS_NOENC.F
c*************************************************************************
c This subroutine evolves spin axis due to torque from other planets for dt
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 preconst,aref ==>  Precession constant and its reference
c                                    semimajor axis (real arrays)
c                 sx,sy,sz      ==>  initial direction of spin (real arrays)
c                 dt            ==>  time step
c                 ielc          ==>  number of encounters (int scalar)
c                 ielst         ==>  list of ecnounters (2D integer*2 array)
c             Output:
c                 sx,sy,sz      ==>  final direction of spin (real arrays)
c Remarks: Based on spin_planets.f, but does not include torques between
c          bodies having encounters.
c
c Authors:  Man Hoi Lee
c Date:    7/19/06
c Last revision: 

      subroutine spin_planets_noenc(nbod,nbodm,mass,xh,yh,zh,preconst,
     &                              aref,sx,sy,sz,dt,ielc,ielst)

      include '../swift.inc'
      include '../symba5/symba5.inc'

c...  Inputs Only: 
      integer nbod,nbodm,ielc
      real*8 mass(nbod),dt
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 preconst(NTPMAX),aref(NTPMAX)
      integer*2 ielst(2,NENMAX)

c...  Inputs and Outputs:
      real*8 sx(NTPMAX),sy(NTPMAX),sz(NTPMAX)

c...  Internals:
      integer i,j,ie
      real*8 dth
      real*8 xhi,yhi,zhi,sxi,syi,szi
      real*8 dx,dy,dz,rsq,r,rdots,faci,facj
      real*8 sxtmp,sytmp,sztmp,torqxi,torqyi,torqzi
      real*8 rh(3,NTPMAX),stmp(3,NTPMAX),torq(3,NTPMAX)

c----
c...  Executable code 

      do i=1,nbod
         rh(1,i) = xh(i)
         rh(2,i) = yh(i)
         rh(3,i) = zh(i)
         torq(1,i) = 0.d0
         torq(2,i) = 0.d0
         torq(3,i) = 0.d0
      enddo

c     Calculate the torques
      do i=2,nbodm

         xhi = rh(1,i)
         yhi = rh(2,i)
         zhi = rh(3,i)

         sxi = sx(i)
         syi = sy(i)
         szi = sz(i)

         torqxi = torq(1,i)
         torqyi = torq(2,i)
         torqzi = torq(3,i)

         do j=i+1,nbod

               dx = rh(1,j) - xhi
               dy = rh(2,j) - yhi
               dz = rh(3,j) - zhi
               rsq = dx*dx + dy*dy + dz*dz
               r = dsqrt(rsq)

               rdots = dx*sxi + dy*syi + dz*szi
               facj = 2.d0*preconst(i)*(mass(j)/mass(1))*(aref(i)/r)**3*
     &                rdots/rsq
               torqxi = torqxi + facj*dx
               torqyi = torqyi + facj*dy
               torqzi = torqzi + facj*dz

               rdots = -(dx*sx(j) + dy*sy(j) + dz*sz(j))
               faci = 2.d0*preconst(j)*(mass(i)/mass(1))*(aref(j)/r)**3*
     &                rdots/rsq
               torq(1,j) = torq(1,j) - faci*dx
               torq(2,j) = torq(2,j) - faci*dy
               torq(3,j) = torq(3,j) - faci*dz

         enddo

         torq(1,i) = torqxi
         torq(2,i) = torqyi
         torq(3,i) = torqzi

      enddo

c...  Now subtract off anyone in an encounter
      do ie=1,ielc

         i = ielst(1,ie)
         j = ielst(2,ie)

         dx = rh(1,j) - rh(1,i)
         dy = rh(2,j) - rh(2,i)
         dz = rh(3,j) - rh(3,i)
         rsq = dx*dx + dy*dy + dz*dz
         r = dsqrt(rsq)

         rdots = dx*sx(i) + dy*sy(i) + dz*sz(i)
         facj = 2.d0*preconst(i)*(mass(j)/mass(1))*(aref(i)/r)**3*
     &          rdots/rsq
         torq(1,i) = torq(1,i) - facj*dx
         torq(2,i) = torq(2,i) - facj*dy
         torq(3,i) = torq(3,i) - facj*dz

         rdots = -(dx*sx(j) + dy*sy(j) + dz*sz(j))
         faci = 2.d0*preconst(j)*(mass(i)/mass(1))*(aref(j)/r)**3*
     &          rdots/rsq
         torq(1,j) = torq(1,j) + faci*dx
         torq(2,j) = torq(2,j) + faci*dy
         torq(3,j) = torq(3,j) + faci*dz

      enddo

c     Take a half step
      dth = 0.5d0*dt

      do i=2,nbod
         sxi = sx(i)
         syi = sy(i)
         szi = sz(i)
         stmp(1,i) = sxi + (torq(2,i)*szi - torq(3,i)*syi)*dth
         stmp(2,i) = syi + (torq(3,i)*sxi - torq(1,i)*szi)*dth
         stmp(3,i) = szi + (torq(1,i)*syi - torq(2,i)*sxi)*dth
      enddo

c     Recalculate the torques
      do i=1,nbod
         torq(1,i) = 0.d0
         torq(2,i) = 0.d0
         torq(3,i) = 0.d0
      enddo

      do i=2,nbodm

         xhi = rh(1,i)
         yhi = rh(2,i)
         zhi = rh(3,i)

         sxi = stmp(1,i)
         syi = stmp(2,i)
         szi = stmp(3,i)

         torqxi = torq(1,i)
         torqyi = torq(2,i)
         torqzi = torq(3,i)

         do j=i+1,nbod

               dx = rh(1,j) - xhi
               dy = rh(2,j) - yhi
               dz = rh(3,j) - zhi
               rsq = dx*dx + dy*dy + dz*dz
               r = dsqrt(rsq)

               rdots = dx*sxi + dy*syi + dz*szi
               facj = 2.d0*preconst(i)*(mass(j)/mass(1))*(aref(i)/r)**3*
     &                rdots/rsq
               torqxi = torqxi + facj*dx
               torqyi = torqyi + facj*dy
               torqzi = torqzi + facj*dz

               rdots = -(dx*stmp(1,j) + dy*stmp(2,j) + dz*stmp(3,j))
               faci = 2.d0*preconst(j)*(mass(i)/mass(1))*(aref(j)/r)**3*
     &                rdots/rsq
               torq(1,j) = torq(1,j) - faci*dx
               torq(2,j) = torq(2,j) - faci*dy
               torq(3,j) = torq(3,j) - faci*dz

         enddo

         torq(1,i) = torqxi
         torq(2,i) = torqyi
         torq(3,i) = torqzi

      enddo

c...  Now subtract off anyone in an encounter
      do ie=1,ielc

         i = ielst(1,ie)
         j = ielst(2,ie)

         dx = rh(1,j) - rh(1,i)
         dy = rh(2,j) - rh(2,i)
         dz = rh(3,j) - rh(3,i)
         rsq = dx*dx + dy*dy + dz*dz
         r = dsqrt(rsq)

         rdots = dx*stmp(1,i) + dy*stmp(2,i) + dz*stmp(3,i)
         facj = 2.d0*preconst(i)*(mass(j)/mass(1))*(aref(i)/r)**3*
     &          rdots/rsq
         torq(1,i) = torq(1,i) - facj*dx
         torq(2,i) = torq(2,i) - facj*dy
         torq(3,i) = torq(3,i) - facj*dz

         rdots = -(dx*stmp(1,j) + dy*stmp(2,j) + dz*stmp(3,j))
         faci = 2.d0*preconst(j)*(mass(i)/mass(1))*(aref(j)/r)**3*
     &          rdots/rsq
         torq(1,j) = torq(1,j) + faci*dx
         torq(2,j) = torq(2,j) + faci*dy
         torq(3,j) = torq(3,j) + faci*dz

      enddo

c     Now take a full step
      do i=2,nbod
         sxtmp = stmp(1,i)
         sytmp = stmp(2,i)
         sztmp = stmp(3,i)
         sx(i) = sx(i) + (torq(2,i)*sztmp - torq(3,i)*sytmp)*dt
         sy(i) = sy(i) + (torq(3,i)*sxtmp - torq(1,i)*sztmp)*dt
         sz(i) = sz(i) + (torq(1,i)*sytmp - torq(2,i)*sxtmp)*dt
      enddo

      return
      end   ! spin_planets_noenc
c---------------------------------------------------------------------
