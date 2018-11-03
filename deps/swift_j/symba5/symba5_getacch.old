c*************************************************************************
c                        SYMBA5_GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 nbodm       ==>  Location of last massive body(int scalar)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c                 mtiny       ==>  Small mass  (real array)
c                ielc           ==>  number of encounters (integer*2 scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Remarks: Based on helio_getacch.f, but does not include the forces of
c          an body B on body A, if body B and A are having an encounter.
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 11/22/97

      subroutine symba5_getacch(nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)

      include '../swift.inc'
      include 'symba5.inc'

c...  Inputs: 
      integer nbod,nbodm
      real*8 mass(nbod),j2rp2,j4rp4,mtiny
      real*8 xh(nbod),yh(nbod),zh(nbod)
      integer*2 ielst(2,NENMAX),ielc

c...  Outputs:
      real*8 axh(nbod),ayh(nbod),azh(nbod)
                
c...  Internals:
      real*8 aoblx(NTPMAX),aobly(NTPMAX),aoblz(NTPMAX) 
      real*8 ir3h(NTPMAX),irh(NTPMAX) 
      integer i,j,ie
      real*8 dx,dy,dz,rji2,faci,facj,irij3

c---
c...  Executable code 

c...  Zero things
      do i=1,nbod
         axh(i) = 0.0
         ayh(i) = 0.0
         azh(i) = 0.0
      enddo

c...  now the third terms
      do i=2,nbodm
         do j=i+1,nbod

            dx = xh(j) - xh(i)
            dy = yh(j) - yh(i)
            dz = zh(j) - zh(i)
            rji2 = dx*dx + dy*dy + dz*dz
            
            irij3 = 1.0d0/(rji2*sqrt(rji2))
            faci = mass(i)*irij3
            facj = mass(j)*irij3
            
            axh(j) = axh(j) - faci*dx
            ayh(j) = ayh(j) - faci*dy
            azh(j) = azh(j) - faci*dz
            
            axh(i) = axh(i) + facj*dx
            ayh(i) = ayh(i) + facj*dy
            azh(i) = azh(i) + facj*dz

         enddo
      enddo
      

c...  Now subtract off anyone in an encounter
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)

         dx = xh(j) - xh(i)
         dy = yh(j) - yh(i)
         dz = zh(j) - zh(i)
         rji2 = dx*dx + dy*dy + dz*dz
         
         irij3 = 1.0d0/(rji2*sqrt(rji2))
         faci = mass(i)*irij3
         facj = mass(j)*irij3
            
         axh(j) = axh(j) + faci*dx
         ayh(j) = ayh(j) + faci*dy
         azh(j) = azh(j) + faci*dz
            
         axh(i) = axh(i) - facj*dx
         ayh(i) = ayh(i) - facj*dy
         azh(i) = azh(i) - facj*dz

      enddo


c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
         do i = 2,nbod
            if(mass(i).ne.0.0d0) then
               axh(i) = axh(i) + aoblx(i) - aoblx(1)
               ayh(i) = ayh(i) + aobly(i) - aobly(1)
               azh(i) = azh(i) + aoblz(i) - aoblz(1)
            endif
         enddo
      endif

      return
      end      ! symba5_getacch

c---------------------------------------------------------------------




