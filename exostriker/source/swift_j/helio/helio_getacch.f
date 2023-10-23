c*************************************************************************
c                        HELIO_GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Author:  Hal Levison  
c Date:    11/14/96
c Last revision: 11/21/96

      subroutine helio_getacch(nbod,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axh(nbod),ayh(nbod),azh(nbod)
                
c...  Internals:
      integer i
      real*8 aoblx(NTPMAX),aobly(NTPMAX),aoblz(NTPMAX) 
      real*8 ir3h(NTPMAX),irh(NTPMAX) ! use NTPMAX so symba can use 

c----
c...  Executable code 

c...  get thr r^-3's
      call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)

      axh(1) = 0.0
      ayh(1) = 0.0
      azh(1) = 0.0

c...  now the third terms
      call getacch_ah3(nbod,mass,xh,yh,zh,axh,ayh,azh)

c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
         do i = 2,nbod
            axh(i) = axh(i) + aoblx(i) - aoblx(1)
            ayh(i) = ayh(i) + aobly(i) - aobly(1)
            azh(i) = azh(i) + aoblz(i) - aoblz(1)
         enddo
      endif

      return
      end      ! helio_getacch

c---------------------------------------------------------------------




