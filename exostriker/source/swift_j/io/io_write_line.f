c*************************************************************************
c                            IO_WRITE_LINE
c*************************************************************************
c write out one line to integer*2 binary file.
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
c Date:    3/1/93 (HAPPY BIRTHDAY TO ME!!!)
c Last revision: 2/21/94

      subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu,id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      integer*2 id2,ia,ie,iinc,inode,iarg,imean
      real*8 scale,ascal,escal
      real*8 asm,abg,esm,ebg

      data scale,ascal,escal/32767.0,6.0,3.0/

c----
c...  Executable code 

      id2 = id

      abg = 10.0**ascal
      if(a.gt.abg) then
         a = abg
      endif
      asm = 10.0**(-1.0*ascal)
      if(a.lt.asm) then
         a = asm
      endif

      ebg = 10.0**escal
      if(e.gt.ebg) then
         e = ebg
      endif
      esm = 10.0**(-1.0*escal)
      if(e.lt.esm) then
         e = esm
      endif

      ia = nint(dlog10(a)*scale/ascal)
      ie = nint(dlog10(e)*scale/escal)
      iinc = nint(inc*scale/TWOPI)
      inode = nint(capom*scale/TWOPI)
      iarg = nint(omega*scale/TWOPI)
      imean = nint(capm*scale/TWOPI)
      write(iu) id2,ia,ie,iinc,inode,iarg,imean

      return
      end      ! io_write_line
c--------------------------------------------------------------------------

