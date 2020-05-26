c*************************************************************************
c                        SYMBA3_GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 mass        ==>  mass of bodies (real array)
c                 lemat         ==>  Matrix of encounters (logical*1 2d array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Remarks: Based on helio_getacch.f, but does not include the forces of
c          an body B on body A, if body B and A are having an encounter.
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 

      subroutine symba3_getacch(nbod,nbodm,mass,lemat,j2rp2,
     &     j4rp4,xh,yh,zh,axh,ayh,azh)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,nbodm
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(nbod),yh(nbod),zh(nbod)
      logical*1 lemat(NTPMAX,NTPMAX)

c...  Outputs:
      real*8 axh(nbod),ayh(nbod),azh(nbod)
                
c...  Internals:
      real*8 aoblx(NTPMAX),aobly(NTPMAX),aoblz(NTPMAX) 
      real*8 ir3h(NTPMAX),irh(NTPMAX) 
      integer i,j
      real*8 dx,dy,dz,rji2,faci,facj,irij3

c----
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

            if( (.not.lemat(i,j)) .and. 
     &        ( (mass(i).gt.TINY) .or. (mass(j).gt.TINY) )) then

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

            endif

         enddo
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
      end      ! symba3_getacch

c---------------------------------------------------------------------




