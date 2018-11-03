c***************************************************************************
c			TU4_GETACCB_TP.F
c*************************************************************************
c GETACCB_TP returns the bary. acc. on EACH of ntp test particles
c due to nbod massive objects by direct summation
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of planets (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xb,yb,zb      ==>  position of planets in beri coord 
c                                    (real arrays)
c                 ntp           ==> number of test particles (int scalar)
c                 xbt,ybt,zbt   ==>  position of test part in beri coord 
c                                    (real arrays)
c                 istat         ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c
c             Output:
c               axbt,aybt,azbt   ==>  accel in beri coord (real arrays) 
c
c Remarks:  Based on Martin's NB4M routines
c Authors:  Martin Duncan 
c Date:    3/8/93
c Last revision: 2/24/94

      subroutine tu4_getaccb_tp(nbod,mass,j2rp2,j4rp4,xb,yb,zb,
     &	              ntp,xbt,ybt,zbt,istat,axbt,aybt,azbt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,istat(ntp)
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xb(nbod),yb(nbod),zb(nbod)
      real*8 xbt(ntp),ybt(ntp),zbt(ntp)

c...  Output
      real*8 axbt(ntp),aybt(ntp),azbt(ntp)

c...  Internals
      real*8 xx,yy,zz,rr2,fac,fac1
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX),irht(NTPMAX)
      real*8 aoblxt(NTPMAX),aoblyt(NTPMAX),aoblzt(NTPMAX)
      integer i,j

c----
c...  executable code

c...  first do the Sun

      do i = 1,ntp
         if(istat(i).eq.0.0) then
            xx = xbt(i) - xb(1)
            yy = ybt(i) - yb(1)
            zz = zbt(i) - zb(1)
            rr2 = xx**2 + yy**2 + zz**2 
            fac1 = 1.d0/sqrt(rr2)
            fac = mass(1)*fac1/rr2

c..         save for the J2 and J4 calculations
            xht(i) = xx
            yht(i) = yy
            zht(i) = zz
            irht(i) = fac1

            axbt(i) = - fac*xx
            aybt(i) = - fac*yy
            azbt(i) = - fac*zz
	 endif
      enddo

c...  do the rest of the planets

      do i = 1,ntp
         if(istat(i).eq.0.0) then
            do j = 2,nbod
               xx = xbt(i) - xb(j)
               yy = ybt(i) - yb(j)
               zz = zbt(i) - zb(j)
               rr2 = xx**2 + yy**2 + zz**2 
               fac = mass(j)/(rr2*sqrt(rr2))
               axbt(i) = axbt(i) - fac*xx
               aybt(i) = aybt(i) - fac*yy
               azbt(i) = azbt(i) - fac*zz
            enddo
	 endif
      enddo

      if(j2rp2.ne.0.0d0) then
         call obl_acc_tp(ntp,istat,mass(1),j2rp2,j4rp4,xht,yht,zht,
     &        irht,aoblxt,aoblyt,aoblzt)
         do i = 1,ntp
            axbt(i) = axbt(i) + aoblxt(i)
            aybt(i) = aybt(i) + aoblyt(i)
            azbt(i) = azbt(i) + aoblzt(i)
         enddo
      endif

      return	
      end                       !  tu4_getaccb_tp
c____________________________________________________________________________

