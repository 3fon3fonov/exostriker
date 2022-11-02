c*************************************************************************
c                        GETACCH_AH0.F
c*************************************************************************
c This subroutine calculates the 0th term of acceleration 
c on the massive particles in the HELIOCENTRIC frame. 
c             Input:
c                 istart      ==>  1st planet in the loop
c                                  istart=3 to get accel for planets
c                                  istart=2 to get accel for tp
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 xh,yh,zh    ==>  position in heliocentric coord (real array)
c                 ir3h         ==> inv radii in heliocentric coord (real array)
c             Output:
c                 axh0,ayh0,azh0 ==>  0th term acceleration in helio coord 
c                                    (real scalor)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 2/18/93

      subroutine getacch_ah0(istart,nbod,mass,xh,yh,zh,ir3h,
     &                        axh0,ayh0,azh0)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod, istart
      real*8 mass(nbod),ir3h(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axh0,ayh0,azh0
                
c...  Internals:
      integer i
      real*8 fac

c----
c...  Executable code 

      axh0 = 0.0d0
      ayh0 = 0.0d0
      azh0 = 0.0d0
      do i=istart,nbod

         fac = mass(i)*ir3h(i)

         axh0 = axh0 - fac*xh(i)
         ayh0 = ayh0 - fac*yh(i)
         azh0 = azh0 - fac*zh(i)

      enddo

      return
      end   ! getacch_ah0

c---------------------------------------------------------------------
