c*************************************************************************
c                            IO_LYAP2_WRITE
c*************************************************************************
c write to lyap2 output file
c
c             Input:
c                 iul           ==>  unit number (int scalar)
c                 time          ==>  current time (real scaler)
c                 dist          ==>  phase space distance (real array)
c                 ntp           ==>  number of objects (int scalar)
c
c Remarks: Based on io_lyap_write.f
c Authors:  Hal Levison
c Date:    7/11/95
c Last revision: 

      subroutine io_lyap2_write(iul,time,dist,ntp)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iul,ntp
      real*8 dist(ntp),time

c...  Internals
      integer i,ierr,iwrite

      data iwrite/0/
      save iwrite

c----
c...  Executable code 


      if(iwrite.eq.0) then
c         call io_open(iul,'lyap2.out','unknown','formatted',ierr)
         call io_open(iul,'lyap2.out','new','formatted',ierr)
         iwrite = 1
      else
         call io_open(iul,'lyap2.out','append','formatted',ierr)
      endif

      if(ierr.ne.0) then
         write(*,*) ' ERROR opening lyap2.out '
         call util_exit(1)
      endif

      write(iul,234) time,(dist(i),i=1,ntp)
 234  format(1p1e12.5,1x,1000(1p1e10.3,1x))

      close(iul)

      return
      end      ! io_lyap2_write

c----------------------------------------------------------------------------
