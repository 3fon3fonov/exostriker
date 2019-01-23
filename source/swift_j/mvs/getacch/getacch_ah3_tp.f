c*************************************************************************
c                        GETACCH_A3_TP.F
c*************************************************************************
c This subroutine calculates the 3rd term acceleration on the test particles
c in the HELIOCENTRIC frame. This term is the direct cross terms
c             Input:
c                 nbod          ==>  number of massive bodies (int scalor)
c                 ntp           ==>  number of test particles (int scalor)
c                 mass          ==>  mass of massive bodies (real array)
c                 xh,yh,zh      ==>  massive position in heliocentric coord 
c                                   (real arrays)
c                 xht,yht,zht   ==>  tp position in heliocentric coord 
c                                   (real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c                 axh3,ayh3,azh3 ==>  3rd term of acceleration in helio coord 
c                                     (real arrays)
c
c Author:  Hal Levison  
c Date:    2/18/93
c Last revision: 

      subroutine getacch_ah3_tp(nbod,ntp,mass,xh,yh,zh,xht,yht,zht,
     &                     istat,axh3,ayh3,azh3)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX)
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)

c...  Outputs:
      real*8 axh3(NTPMAX),ayh3(NTPMAX),azh3(NTPMAX)

c...  Internals:
      integer i,j
      real*8 dx,dy,dz,rji2,fac,irij3

c------
c...  Executable code

      do i=1,ntp
         axh3(i) = 0.0
         ayh3(i) = 0.0
         azh3(i) = 0.0
      enddo

      do j=1,ntp
         if(istat(j).eq.0) then
            do i=2,nbod

               dx = xht(j) - xh(i)
               dy = yht(j) - yh(i)
               dz = zht(j) - zh(i)
               rji2 = dx*dx + dy*dy + dz*dz

               irij3 = 1.0d0/(rji2*sqrt(rji2))
               fac = mass(i)*irij3

               axh3(j) = axh3(j) - fac*dx
               ayh3(j) = ayh3(j) - fac*dy
               azh3(j) = azh3(j) - fac*dz

            enddo
         endif
      enddo

      return
      end     ! getacch_ah3_tp
c--------------------------------------------------------------------------

