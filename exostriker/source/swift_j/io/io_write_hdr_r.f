c*************************************************************************
c                            IO_WRITE_HDR_R
c*************************************************************************
c write out header part of the real*4 binary file
c
c             Input:
c                 iu              ==> unit number to write to
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 istat           ==>  status of the test paricles
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      subroutine io_write_hdr_r(iu,time,nbod,ntp,istat) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX,NSTAT),iu
      real*8 time

c...  Internals
      integer i
      real*4 ttmp
      integer*2 nleft,nbod2

c----
c...  Executable code 


c...  calculate number of remaining test particles
      nleft = 0
      do i=1,ntp
         if(istat(i,1).eq.0) then
            nleft = nleft + 1
         endif
      enddo

      nbod2 = nbod

      ttmp = time
 
      write(iu) ttmp,nbod2,nleft

      return
      end     ! io_write_hdr_r.f
c---------------------------------------------------------------------------

