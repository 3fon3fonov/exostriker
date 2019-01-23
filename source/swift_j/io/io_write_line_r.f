c*************************************************************************
c                            IO_WRITE_LINE_R
c*************************************************************************
c write out one line to real*4 binary file.
c
c      Input:
c            iu       ==> unit number to write to
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      subroutine io_write_line_r(iu,id,a,e,inc,capom,omega,capm) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu,id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      integer*2 id2
      real*4 a4,e4,inc4,capom4,omega4,capm4


c----
c...  Executable code 

      id2 = id

      a4 = a
      e4 = e
      inc4 = inc
      capom4 = capom
      capm4 = capm
      omega4 = omega

      write(iu) id2,a4,e4,inc4,capom4,omega4,capm4

      return
      end      ! io_write_line_r
c--------------------------------------------------------------------------
