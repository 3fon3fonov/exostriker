c*************************************************************************
c                            DISCARD_MASS_REORDER.F
c*************************************************************************
c Remove a massive body
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ip            ==>  planets to remove (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 rhill         ==>  size of a planet's hill's sphere.
c                                    (real array)
c                 isperip       ==>  Perihelion check array (integer array)
c                 lemat         ==> encounter matrix (logical 2D array)
c                 isperih       ==> heliocentric peri flags. (real array)
c             Output:
c                 ip            ==>  planets to remove (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 rhill         ==>  size of a planet's hill's sphere.
c                                    (real array)
c                 isperip       ==>  Perihelion check array (integer array)
c                 lemat         ==> encounter matrix (logical 2D array)
c                 isperih       ==> heliocentric peri flags. (real array)
c
c Remarks: 
c
c Authors:  Hal Levison 
c Date:    1/2/97
c Last revision: 1/8/97

      subroutine discard_mass_reorder(ip,nbod,mass,xh,yh,zh,
     &           vxh,vyh,vzh,rpl,rhill,isperip,lemat,isperih)

      include '../swift.inc'

c...  Inputs: 
      integer ip

c...  Input and Output
      integer nbod
      logical*1 lemat(NTPMAX,NTPMAX)
      integer*2 isperip(NTPMAX,NTPMAX)
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod),rpl(nbod)
      real*8 rhill(nbod)
      integer isperih(nbod)

c...  internal
      integer i,j

c-----
c...  Executable code 

      do i=ip,nbod-1
         xh(i) = xh(i+1)
         yh(i) = yh(i+1)
         zh(i) = zh(i+1)
         vxh(i) = vxh(i+1)
         vyh(i) = vyh(i+1)
         vzh(i) = vzh(i+1)
         mass(i) = mass(i+1)
         rpl(i) = rpl(i+1)
         rhill(i) = rhill(i+1)
         isperih(i) = isperih(i+1)
         do j=1,nbod
            lemat(i,j) = lemat(i+1,j)
            isperip(i,j) = isperip(i+1,j)
         enddo
      enddo
      do i=ip,nbod-1
         do j=1,nbod
            lemat(j,i) = lemat(j,i+1)
            isperip(j,i) = isperip(j,i+1)
         enddo
      enddo
      nbod = nbod - 1

      return
      end
