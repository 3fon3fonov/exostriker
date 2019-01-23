c***************************************************************************
c			OBL_ACC_TP.F
c*************************************************************************
c OBL_ACC returns the BARYCENTRIC x,y,z components of the acc. on NTP
c test particles due to the oblateness of mass(1) using  
c the values of J2RP2 and J4RP4 passed into the routine.
c (J2RP2 for example is the product of 
c J_2 times the square of the central body's radius)
c Here we return the net acc. produced
c only by the J2 and J4 terms (i.e. including
c neither the monopole nor higher order terms).
c	
c
c             Input:
c                 ntp      ==>  number of massive bodies (incl. central one)
c                 istat    ==>  status of the test paricles
c                                    (integer array)
c                                istat(i) = 0 ==> active:  = 1 not
c                                NOTE: it is really a 2d array but 
c                                      we only use the 1st row
c                 msun     ==>  mass of the central particle (real*8 scalar)
c                 j2rp2    ==>  scaled value of j2 moment (real*8 scalar)
c                 j4rp4    ==>  scaled value of j4 moment (real*8 scalar)
c   xht(*),yht(*),zht(*)   ==>  HELIO. positions of particles
c                irht(*)   ==> 1./ magnitude of radius vector (real*8 vector)
c                                (passed in to save calcs.)
c             Output:
c  aoblxt(*),aoblyt(*),aoblzt(*)  ==>  BARY. components of accel 
c                                        (real*8 vectors) 
c
c Remarks:  Based on Martin's OBL_ACC.F
c Authors:  Hal Levison
c Date:    3/4/94
c Last revision: 

      subroutine obl_acc_tp(ntp,istat,msun,j2rp2,j4rp4,xht,yht,zht,
     &        irht,aoblxt,aoblyt,aoblzt)

      include '../swift.inc'

c...  Inputs Only: 
      integer ntp,istat(ntp)
      real*8 msun,j2rp2,j4rp4
      real*8 xht(ntp),yht(ntp),zht(ntp),irht(NPLMAX)

c...  Output
      real*8 aoblxt(ntp),aoblyt(ntp),aoblzt(ntp)

c...  Internals
      integer n
      real*8 rinv2,t0,t1,t2,t3
      real*8 fac1,fac2

c----
c...  executable code

c First get the bary acc. of each tp due to central oblate "sun"

      do n =1,ntp

         rinv2= irht(n)**2
         t0 = -msun*rinv2*rinv2*irht(n)
         t1 = 1.5d0 *j2rp2
         t2 = zht(n)*zht(n)*rinv2
         t3 = 1.875d0 *j4rp4*rinv2
         
         fac1 = t0*(t1 - t3 - (5.d0*t1 - (14.d0 - 21.d0*t2)*t3)*t2)
         fac2 = 2.d0*t0*(t1 - (2.d0 - (14.d0*t2/3.d0))*t3)
         
         aoblxt(n) = fac1*xht(n)
         aoblyt(n) = fac1*yht(n)
         aoblzt(n) = (fac1 + fac2)*zht(n)
         
      enddo
      
      return	
      end                       !  obl_acc_tp.f
c____________________________________________________________________________
