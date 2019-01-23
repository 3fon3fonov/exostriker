c*************************************************************************
c                            SPIN_SUN.F
c*************************************************************************
c This subroutine evolves spin axis due to torque from Sun for dt
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
c             Output:
c                 sx,sy,sz      ==>  final direction of spin (real arrays)
c Remarks: 
c
c Authors:  Man Hoi Lee
c Date:    7/18/06
c Last revision: 

      subroutine spin_sun(nbod,nbodm,mass,xh,yh,zh,preconst,aref,
     &                    sx,sy,sz,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,nbodm
      real*8 mass(nbod),dt
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 preconst(NTPMAX),aref(NTPMAX)

c...  Inputs and Outputs:
      real*8 sx(NTPMAX),sy(NTPMAX),sz(NTPMAX)

c...  Internals:
      integer n
      real*8 rsq,r,rdots,rxsx,rxsy,rxsz,phi,cosphi,sinphi,fac1,fac2

c----
c...  Executable code 

      do n=2,nbod
         rsq = xh(n)**2 + yh(n)**2 + zh(n)**2
         r = dsqrt(rsq)
         rdots = xh(n)*sx(n) + yh(n)*sy(n) + zh(n)*sz(n)
         rxsx = yh(n)*sz(n) - zh(n)*sy(n)
         rxsy = zh(n)*sx(n) - xh(n)*sz(n)
         rxsz = xh(n)*sy(n) - yh(n)*sx(n)
         phi = 2.d0*preconst(n)*(aref(n)/r)**3*(rdots/r)*dt
         cosphi = dcos(phi)
         sinphi = dsin(phi)
         fac1 = (1.d0 - cosphi)*rdots/rsq
         fac2 = sinphi/r
         sx(n) = cosphi*sx(n) + fac1*xh(n) + fac2*rxsx
         sy(n) = cosphi*sy(n) + fac1*yh(n) + fac2*rxsy
         sz(n) = cosphi*sz(n) + fac1*zh(n) + fac2*rxsz
      enddo

      return
      end   ! spin_sun
c---------------------------------------------------------------------
