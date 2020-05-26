c*************************************************************************
c                            IO_READ_LINE
c*************************************************************************
c read one line from integer*2 binary file.
c
c      Input:
c            iu       ==> unit number to write to
c      Output:
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c       Returns:
c      io_read_line    ==>   =0 read ok
c                           !=0 read failed is set to iostat variable
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    3/1/93 (HAPPY BIRTHDAY TO ME!!!)
c Last revision: 2/21/94

      integer function io_read_line(iu,id,a,e,inc,capom,omega,capm) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output: 
      integer id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      integer*2 id2,ia,ie,iinc,inode,iarg,imean
      real*8 scale,ascal,escal
      integer ierr

      data scale,ascal,escal/32767.0,6.0,3.0/

c----
c...  Executable code 

      read(iu,iostat=ierr) id2,ia,ie,iinc,inode,iarg,imean
      io_read_line = ierr
      if(ierr.ne.0) then
         return
      endif

      id = id2

      a = 10.0**(float(ia)*ascal/scale)
      e = 10.0**(float(ie)*escal/scale)

      inc = float(iinc)*TWOPI/scale
      capom = float(inode)*TWOPI/scale
      omega = float(iarg)*TWOPI/scale
      capm = float(imean)*TWOPI/scale

      return
      end      ! io_read_line
c--------------------------------------------------------------------------

