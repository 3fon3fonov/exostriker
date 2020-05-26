c*************************************************************************
c                            MVITS_ROTATE.F
c*************************************************************************
c Rotates a coordinate frame
c
c             Input:
c                 sgn                ==> direction of rotation, must be +/-1
c                                             (real scalar)
c                 ntp                ==>  number of test bodies (int scalar)
c                 cosrot,sinrot      ==>  cos and sin of rotation angle of tp
c                                          (real arrays)
c                 xht,yht        ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht      ==>  initial velocity in helio coord 
c                                        (real arrays)
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    12/16/93
c Last revision: 

      subroutine mvits_rotate(sgn,ntp,cosrot,sinrot,xht,yht,vxht,vyht)

      include '../swift.inc'
      include 'mvits.inc'

c...  Inputs Only: 
      integer ntp
      real*8 sgn
      real*8 cosrot(NTPMAX),sinrot(NTPMAX)

c...  Inputs and Outputs:
      real*8 xht(ntp),yht(ntp)
      real*8 vxht(ntp),vyht(ntp)

c...  Internals:
      integer i
      real*8 x,y

c----
c...  Executable code 

      do i=1,ntp

         x = xht(i)*cosrot(i) - sgn*yht(i)*sinrot(i)
         y = yht(i)*cosrot(i) + sgn*xht(i)*sinrot(i)
         xht(i) = x
         yht(i) = y

         x = vxht(i)*cosrot(i) - sgn*vyht(i)*sinrot(i)
         y = vyht(i)*cosrot(i) + sgn*vxht(i)*sinrot(i)
         vxht(i) = x
         vyht(i) = y

      enddo

      return
      end      ! mvits_rotate
c---------------------------------------------------------------------
