c*************************************************************************
c                        GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xj,yj,zj    ==>  position in jacobi coord (real arrays)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 2/18/93

      subroutine getacch(nbod,mass,j2rp2,j4rp4,xj,yj,zj,
     &     xh,yh,zh,axh,ayh,azh)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(NPLMAX),xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),j2rp2,j4rp4
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)

c...  Outputs:
      real*8 axh(NPLMAX),ayh(NPLMAX),azh(NPLMAX)
                
c...  Internals:
      integer i
      real*8 ir3h(NPLMAX),ir3j(NPLMAX)
      real*8 irh(NPLMAX),irj(NPLMAX)
      real*8 axh1(NPLMAX),ayh1(NPLMAX),azh1(NPLMAX)
      real*8 axh2(NPLMAX),ayh2(NPLMAX),azh2(NPLMAX)
      real*8 axh3(NPLMAX),ayh3(NPLMAX),azh3(NPLMAX)
      real*8 axh0,ayh0,azh0
      real*8 aoblx(NPLMAX),aobly(NPLMAX),aoblz(NPLMAX)

c----
c...  Executable code 

c...  get thr r^-3's
      call getacch_ir3(nbod,2,xj,yj,zj,ir3j,irj)
      call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)

c...  calc the ah0's:  recall that they are the same for all particles
      call getacch_ah0(3,nbod,mass,xh,yh,zh,ir3h,axh0,ayh0,azh0) 

c...  now the first terms
      call getacch_ah1(nbod,mass,xh,yh,zh,xj,yj,zj,ir3h,ir3j,
     &                 axh1,ayh1,azh1)

c...  now the second terms
      call getacch_ah2(nbod,mass,xj,yj,zj,ir3j,axh2,ayh2,azh2)

c...  now the third terms
      call getacch_ah3(nbod,mass,xh,yh,zh,axh3,ayh3,azh3)

c...  add them all together
      axh(1) = 0.0
      ayh(1) = 0.0
      azh(1) = 0.0
      do i=2,nbod
        axh(i) = axh0 + axh1(i) + axh2(i) + axh3(i)
        ayh(i) = ayh0 + ayh1(i) + ayh2(i) + ayh3(i)
        azh(i) = azh0 + azh1(i) + azh2(i) + azh3(i)
      enddo

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
      end      ! getacch

c---------------------------------------------------------------------




