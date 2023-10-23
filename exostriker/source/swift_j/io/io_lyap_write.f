c*************************************************************************
c                            IO_LYAP_WRITE
c*************************************************************************
c write to lyap output file
c
c             Input:
c                 iul           ==>  unit number (int scalar)
c                 tin           ==>  current time (real scaler)
c                 logpr         ==>  gamma (real array)
c                 ntp           ==>  number of objects (int scalar)
c
c Remarks: 
c Authors:  Hal Levison
c Date:    8/12/93  
c Last revision: 3/3/94

      subroutine io_lyap_write(iul,tin,logpr,ntp)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iul,ntp
      real*8 logpr(ntp),tin

c...  Internals
      integer i,ierr

c----
c...  Executable code 

      call io_open(iul,'lyap.out','append','formatted',ierr)

      write(iul,234) tin,(logpr(i),i=1,ntp)
 234  format(1p1e12.5,1x,1000(1p1e10.3,1x))

      close(iul)

      return
      end      ! io_lyap_write

c----------------------------------------------------------------------------
