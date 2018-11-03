c*************************************************************************
c                        GETACCH_A3.F
c*************************************************************************
c This subroutine calculates the 3rd term acceleration on the massive particles
c in the HELIOCENTRIC frame. This term is the direct cross terms
c             Input:
c                 nbod          ==>  number of massive bodies (int scalor)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  position in heliocentric coord 
c                                   (real arrays)
c             Output:
c                 axh3,ayh3,azh3 ==>  3rd term of acceleration in helio coord 
c                                     (real arrays)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 11/21/96

      subroutine getacch_ah3(nbod,mass,xh,yh,zh,axh3,ayh3,azh3)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axh3(nbod),ayh3(nbod),azh3(nbod)

c...  Internals:
      integer i,j
      real*8 dx,dy,dz,rji2,faci,facj,irij3

c------
c...  Executable code

      do i=1,nbod
         axh3(i) = 0.0
         ayh3(i) = 0.0
         azh3(i) = 0.0
      enddo

      do i=2,nbod-1
         do j=i+1,nbod

             dx = xh(j) - xh(i)
             dy = yh(j) - yh(i)
             dz = zh(j) - zh(i)
             rji2 = dx*dx + dy*dy + dz*dz

             irij3 = 1.0d0/(rji2*sqrt(rji2))
             faci = mass(i)*irij3
             facj = mass(j)*irij3

             axh3(j) = axh3(j) - faci*dx
             ayh3(j) = ayh3(j) - faci*dy
             azh3(j) = azh3(j) - faci*dz

             axh3(i) = axh3(i) + facj*dx
             ayh3(i) = ayh3(i) + facj*dy
             azh3(i) = azh3(i) + facj*dz

         enddo
      enddo

      return
      end     ! getacch_ah3

