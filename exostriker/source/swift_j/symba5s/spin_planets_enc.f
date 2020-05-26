c*************************************************************************
c                             SPIN_PLANETS_ENC.F
c*************************************************************************
c Spin evolution during encounters
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 irec          ==>  recursion level  (integer scalar)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 ielev         ==>  The level that this particle should go
c                                             (int*2 array)
c                 rhill         ==>  Hill sphere of planet (real Scalar)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                dt             ==>  timestep  (real scalar)
c                ielc           ==>  number of encounters (int scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c                 preconst,aref ==>  Precession constant and its reference
c                                    semimajor axis (real arrays)
c                 sx,sy,sz      ==>  initial direction of spin (real arrays)
c            Output:
c                 sx,sy,sz      ==>  final direction of spin (real arrays)
c
c Remarks: Based on spin_planets.f and symba5_kick.f
c Authors:  Man Hoi Lee
c Date:   7/20/06
c Last revision: 

      subroutine spin_planets_enc(nbod,mass,irec,iecnt,ielev,rhill,
     &        xh,yh,zh,dt,ielc,ielst,preconst,aref,sx,sy,sz)

      include '../swift.inc'
      include '../symba5/symba5.inc'

c...  Inputs Only: 
      integer nbod,irec,ielc
      real*8 mass(nbod),dt,rhill(nbod)
      integer*2 iecnt(NTPMAX),ielev(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      integer*2 ielst(2,NENMAX)
      real*8 preconst(NTPMAX),aref(NTPMAX)

c...  Inputs and Outputs:
      real*8 sx(NTPMAX),sy(NTPMAX),sz(NTPMAX)

c...  Internals: 
      real*8 dth,sxi,syi,szi,sxtmp,sytmp,sztmp,sgn
      real*8 dx,dy,dz,fac,ris,r,rsq,rdots
      real*8 ri,rr,faci,facj,rim1
      integer i,j,k,irm1,irecl,ie,nsgn
      real*8 stmp(3,NTPMAX),torq(3,NTPMAX)

c----
c...  Executable code 

      do i=1,nbod
         stmp(1,i) = sx(i)
         stmp(2,i) = sy(i)
         stmp(3,i) = sz(i)
      enddo

c     Take a half step

c     Loop through once with sgn=+1 if irec = 1 or twice with sgn=+-1.
      if (irec.eq.1) then
         nsgn = 1
      else
         nsgn = 2
      endif

      sgn = -1.0d0

      do k=1,nsgn

         sgn = -sgn

         irm1 = irec - 1
         if (sgn.lt.0.0d0) then
            irecl = irec - 1
         else
            irecl = irec
         endif

         do i=1,nbod
            torq(1,i) = 0.d0
            torq(2,i) = 0.d0
            torq(3,i) = 0.d0
         enddo

c        Calculate the torques
         do ie=1,ielc

            i = ielst(1,ie)
            j = ielst(2,ie)

            if ((ielev(i).ge.irm1) .and. (ielev(j).ge.irm1)) then

               ri = (rhill(i)+rhill(j))**2 * RHSCALE*RHSCALE * 
     &              (RSHELL**(2*irecl))
               rim1 = ri*RSHELL*RSHELL
            
               dx = xh(j) - xh(i)
               dy = yh(j) - yh(i)
               dz = zh(j) - zh(i)
               rsq = dx*dx + dy*dy + dz*dz
               r = dsqrt(rsq)

               if (rsq.lt.rim1) then
                  fac = 0.0d0
               else if (rsq.lt.ri) then
                  ris = sqrt(ri)
                  rr = (ris-r)/(ris*(1.0d0-RSHELL))
                  fac = (1.0d0 - 3.0d0*rr*rr + 2.0d0*(rr**3))
               else
                  fac = 1.0d0
               endif

               rdots = dx*sx(i) + dy*sy(i) + dz*sz(i)
               facj = fac*2.d0*preconst(i)*(mass(j)/mass(1))*
     &                (aref(i)/r)**3*rdots/rsq
               torq(1,i) = torq(1,i) + facj*dx
               torq(2,i) = torq(2,i) + facj*dy
               torq(3,i) = torq(3,i) + facj*dz

               rdots = -(dx*sx(j) + dy*sy(j) + dz*sz(j))
               faci = fac*2.d0*preconst(j)*(mass(i)/mass(1))*
     &                (aref(j)/r)**3*rdots/rsq
               torq(1,j) = torq(1,j) - faci*dx
               torq(2,j) = torq(2,j) - faci*dy
               torq(3,j) = torq(3,j) - faci*dz
            
            endif

         enddo

c        Apply the torques for a half step
         dth = 0.5d0*dt

         do i=2,nbod
            if ((iecnt(i).ne.0) .and. (ielev(i).ge.irm1)) then
               sxi = sx(i)
               syi = sy(i)
               szi = sz(i)
               stmp(1,i) = stmp(1,i) + (torq(2,i)*szi - torq(3,i)*syi)*
     &                                 dth*sgn
               stmp(2,i) = stmp(2,i) + (torq(3,i)*sxi - torq(1,i)*szi)*
     &                                 dth*sgn
               stmp(3,i) = stmp(3,i) + (torq(1,i)*syi - torq(2,i)*sxi)*
     &                                 dth*sgn
            endif
         enddo

      enddo

c     Now take a full step

c     Loop through once with sgn=+1 if irec = 1 or twice with sgn=+-1.
      if (irec.eq.1) then
         nsgn = 1
      else
         nsgn = 2
      endif

      sgn = -1.0d0

      do k=1,nsgn

         sgn = -sgn

         irm1 = irec - 1
         if (sgn.lt.0.0d0) then
            irecl = irec - 1
         else
            irecl = irec
         endif

         do i=1,nbod
            torq(1,i) = 0.d0
            torq(2,i) = 0.d0
            torq(3,i) = 0.d0
         enddo

c        Recalculate the torques
         do ie=1,ielc

            i = ielst(1,ie)
            j = ielst(2,ie)

            if ((ielev(i).ge.irm1) .and. (ielev(j).ge.irm1)) then

               ri = (rhill(i)+rhill(j))**2 * RHSCALE*RHSCALE * 
     &              (RSHELL**(2*irecl))
               rim1 = ri*RSHELL*RSHELL
            
               dx = xh(j) - xh(i)
               dy = yh(j) - yh(i)
               dz = zh(j) - zh(i)
               rsq = dx*dx + dy*dy + dz*dz
               r = dsqrt(rsq)

               if (rsq.lt.rim1) then
                  fac = 0.0d0
               else if (rsq.lt.ri) then
                  ris = sqrt(ri)
                  rr = (ris-r)/(ris*(1.0d0-RSHELL))
                  fac = (1.0d0 - 3.0d0*rr*rr + 2.0d0*(rr**3))
               else
                  fac = 1.0d0
               endif

               rdots = dx*stmp(1,i) + dy*stmp(2,i) + dz*stmp(3,i)
               facj = fac*2.d0*preconst(i)*(mass(j)/mass(1))*
     &                (aref(i)/r)**3*rdots/rsq
               torq(1,i) = torq(1,i) + facj*dx
               torq(2,i) = torq(2,i) + facj*dy
               torq(3,i) = torq(3,i) + facj*dz

               rdots = -(dx*stmp(1,j) + dy*stmp(2,j) + dz*stmp(3,j))
               faci = fac*2.d0*preconst(j)*(mass(i)/mass(1))*
     &                (aref(j)/r)**3*rdots/rsq
               torq(1,j) = torq(1,j) - faci*dx
               torq(2,j) = torq(2,j) - faci*dy
               torq(3,j) = torq(3,j) - faci*dz
            
            endif

         enddo

c        Apply the torques for a full step
         do i=2,nbod
            if ((iecnt(i).ne.0) .and. (ielev(i).ge.irm1)) then
               sxtmp = stmp(1,i)
               sytmp = stmp(2,i)
               sztmp = stmp(3,i)
               sx(i) = sx(i) + (torq(2,i)*sztmp - torq(3,i)*sytmp)*
     &                         dt*sgn
               sy(i) = sy(i) + (torq(3,i)*sxtmp - torq(1,i)*sztmp)*
     &                         dt*sgn
               sz(i) = sz(i) + (torq(1,i)*sytmp - torq(2,i)*sxtmp)*
     &                         dt*sgn
            endif
         enddo

      enddo

      return
      end      ! spin_planets_enc.f
c--------------------------------------------------------------
