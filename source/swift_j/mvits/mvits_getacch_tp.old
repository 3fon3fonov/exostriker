c*************************************************************************
c                        MVITS_GETACCH_TP.F
c*************************************************************************
c This subroutine calculates the acceleration on the test particles
c in the HELIOCENTRIC frame. Used with MVITS. 
c             Input:
c                  nbod        ==>  number of massive bodies (int scalor)
c                  ntp         ==>  number of tp bodies (int scalor)
c                  mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4  ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                          (real scalars)
c                  xh,yh,zh    ==>  massive part position in helio coord 
c                                     (real arrays)
c                  xht,yht,zht ==>  test part position in heliocentric coord 
c                                     (real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c               axht,ayht,azht ==>  tp acceleration in helio coord 
c                                   (2d real arrays)
c
c Author:  Hal Levison  
c Date:    12/16/93
c Last revision: 2/24/94

      subroutine mvits_getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xh,yh,zh,
     &     xht,yht,zht,istat,axht,ayht,azht)

      include '../swift.inc'
      include 'mvits.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX)
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX),j2rp2,j4rp4

c...  Outputs:
      real*8 axht(NTPMAX,NPLMAX),ayht(NTPMAX,NPLMAX),azht(NTPMAX,NPLMAX)
                
c...  Internals:
      integer i,j
      real*8 ir3h(NPLMAX),irh(NPLMAX)
      real*8 ir3ht(NPLMAX),irht(NPLMAX)
      real*8 fac,dx,dy,dz,rji2,irij3
      real*8 axh0(NPLMAX),ayh0(NPLMAX),azh0(NPLMAX)
      real*8 aoblx(NPLMAX),aobly(NPLMAX),aoblz(NPLMAX)
      real*8 aoblxt(NTPMAX),aoblyt(NTPMAX),aoblzt(NTPMAX)

c----
c...  Executable code 

c...  get thr r^-3's  for the planets
      call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
      call getacch_ir3(ntp,1,xht,yht,zht,ir3ht,irht)

c...  calc the ah0's:  recall that they are the same for all particles
      do i=2,nbod
         fac = mass(i)*ir3h(i)
         axh0(i) =  - fac*xh(i)
         ayh0(i) =  - fac*yh(i)
         azh0(i) =  - fac*zh(i)
      enddo

c...  the first terms are 0

c...  the second terms are 0

c...  now the third terms
      do j=1,ntp
         if(istat(j).eq.0) then
            do i=2,nbod

               dx = xht(j) - xh(i)
               dy = yht(j) - yh(i)
               dz = zht(j) - zh(i)
               rji2 = dx*dx + dy*dy + dz*dz

               irij3 = 1.0d0/(rji2*sqrt(rji2))
               fac = mass(i)*irij3

               axht(j,i) = axh0(i) - fac*dx
               ayht(j,i) = ayh0(i) - fac*dy
               azht(j,i) = azh0(i) - fac*dz

            enddo
         endif
      enddo

c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
         call obl_acc_tp(ntp,istat,mass(1),j2rp2,j4rp4,xht,yht,zht,
     &        irht,aoblxt,aoblyt,aoblzt)
         do i = 1,ntp
            axh0(i) = axh0(i) + aoblxt(i) - aoblx(1)
            ayh0(i) = ayh0(i) + aoblyt(i) - aobly(1)
            azh0(i) = azh0(i) + aoblzt(i) - aoblz(1)
         enddo
      endif

      return
      end      ! mvits_getacch_tp

c---------------------------------------------------------------------




