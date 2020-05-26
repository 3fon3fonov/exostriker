c*************************************************************************
c                        SYMBA5_GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalar)
c                 nbodm       ==>  Location of last massive body(int scalar)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c                 mtiny       ==>  Small mass  (real array)
c                ielc           ==>  number of encounters (int scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Remarks: Based on helio_getacch.f, but does not include the forces of
c          an body B on body A, if body B and A are having an encounter.
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 5/28/99

      subroutine symba5_getacch(nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)

      include '../swift.inc'
      include 'symba5.inc'

c...  Inputs: 
      integer nbod,nbodm,ielc
      real*8 mass(nbod),j2rp2,j4rp4,mtiny
      real*8 xh(nbod),yh(nbod),zh(nbod)
      integer*2 ielst(2,NENMAX)

c...  Outputs:
      real*8 axh(nbod),ayh(nbod),azh(nbod)
                
c...  Internals:
      real*8 aoblx(NTPMAX),aobly(NTPMAX),aoblz(NTPMAX) 
      real*8 ir3h(NTPMAX),irh(NTPMAX) 
      integer i,j,ie
      real*8 dx,dy,dz,rji2,faci,facj,irij3
      real*8 massi,xhi,yhi,zhi,axhi,ayhi,azhi
      real*8 rh(3,NTPMAX),ah(3,NTPMAX)

c---
c...  Executable code 

c...  Initialize things
      do i=1,nbod
         rh(1,i) = xh(i)
         rh(2,i) = yh(i)
         rh(3,i) = zh(i)
         ah(1,i) = 0.0
         ah(2,i) = 0.0
         ah(3,i) = 0.0
      enddo

c...  now the third terms
      do i=2,nbodm

         massi = mass(i)
         xhi = rh(1,i)
         yhi = rh(2,i)
         zhi = rh(3,i)
         axhi = ah(1,i)
         ayhi = ah(2,i)
         azhi = ah(3,i)

         do j=i+1,nbod

            dx = rh(1,j) - xhi
            dy = rh(2,j) - yhi
            dz = rh(3,j) - zhi
            rji2 = dx*dx + dy*dy + dz*dz
            
            irij3 = 1.0d0/(rji2*sqrt(rji2))
            faci = massi*irij3
            facj = mass(j)*irij3
            
            ah(1,j) = ah(1,j) - faci*dx
            ah(2,j) = ah(2,j) - faci*dy
            ah(3,j) = ah(3,j) - faci*dz
            
            axhi = axhi + facj*dx
            ayhi = ayhi + facj*dy
            azhi = azhi + facj*dz

         enddo

         ah(1,i) = axhi
         ah(2,i) = ayhi
         ah(3,i) = azhi

      enddo
      

c...  Now subtract off anyone in an encounter
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)

         dx = rh(1,j) - rh(1,i)
         dy = rh(2,j) - rh(2,i)
         dz = rh(3,j) - rh(3,i)
         rji2 = dx*dx + dy*dy + dz*dz
         
         irij3 = 1.0d0/(rji2*sqrt(rji2))
         faci = mass(i)*irij3
         facj = mass(j)*irij3
            
         ah(1,j) = ah(1,j) + faci*dx
         ah(2,j) = ah(2,j) + faci*dy
         ah(3,j) = ah(3,j) + faci*dz
            
         ah(1,i) = ah(1,i) - facj*dx
         ah(2,i) = ah(2,i) - facj*dy
         ah(3,i) = ah(3,i) - facj*dz

      enddo


c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
         do i = 2,nbod
            if(mass(i).ne.0.0d0) then
               ah(1,i) = ah(1,i) + aoblx(i) - aoblx(1)
               ah(2,i) = ah(2,i) + aobly(i) - aobly(1)
               ah(3,i) = ah(3,i) + aoblz(i) - aoblz(1)
            endif
         enddo
      endif


      do i=1,nbod
         axh(i) = ah(1,i)
         ayh(i) = ah(2,i)
         azh(i) = ah(3,i)
      enddo

      return
      end      ! symba5_getacch

c---------------------------------------------------------------------
