c*************************************************************************
c                        GETACCH_AH1.F
c*************************************************************************
c This subroutine calculates the 1st term of acceleration 
c on the massive particles in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 xh,yh,zh    ==>  position in heliocentric coord (real array)
c                 xj,yj,zj    ==>  position in jacobi coord (real array)
c                 ir3h        ==> inv radii in heliocentric coord (real array)
c                 ir3j        ==> inv radii in jacobi coord (real array)
c             Output:
c                 axh1,ayh1,azh1 ==>  1st term acceleration in helio coord 
c                                    (real array)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 2/2/93

      subroutine getacch_ah1(nbod,mass,xh,yh,zh,xj,yj,zj,ir3h,ir3j,
     &                 axh1,ayh1,azh1)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),ir3h(nbod),ir3j(nbod)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axh1(nbod),ayh1(nbod),azh1(nbod)

                
c...  Internals:
      integer i
      real*8 ah1h,ah1j

c----
c...  Executable code 

      axh1(1) = 0.0
      ayh1(1) = 0.0
      azh1(1) = 0.0

      axh1(2) = 0.0     ! because xj=xh
      ayh1(2) = 0.0
      azh1(2) = 0.0

      do i=3,nbod

         ah1j = xj(i)*ir3j(i)
         ah1h = xh(i)*ir3h(i)
         axh1(i) = mass(1)*(ah1j - ah1h)

         ah1j = yj(i)*ir3j(i)
         ah1h = yh(i)*ir3h(i)
         ayh1(i) = mass(1)*(ah1j - ah1h)

         ah1j = zj(i)*ir3j(i)
         ah1h = zh(i)*ir3h(i)
         azh1(i) = mass(1)*(ah1j - ah1h)

      enddo

      return
      end   ! getacch_ah1

c---------------------------------------------------------------------


