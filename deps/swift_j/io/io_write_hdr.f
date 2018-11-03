c*************************************************************************
c                            IO_WRITE_HDR
c*************************************************************************
c write out header part of the integer*2 binary file
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
c Date:    3/1/93 (HAPPY BIRTHDAY TO ME!!!)
c Last revision: 2/21/94

      subroutine io_write_hdr(iu,time,nbod,ntp,istat) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,ntp,istat(ntp),iu
      real*8 time

c...  Internals
      integer i
      real*8 ttmp,t1,t2,scale,tscal
      integer*2 it1,it2,it3,nleft,nbod2

      data scale,tscal/32767.0,100.0/

c----
c...  Executable code 


c...  calculate number of remaining test particles
      nleft = 0
      do i=1,ntp
         if(istat(i).eq.0) then
            nleft = nleft + 1
         endif
      enddo

      nbod2 = nbod

      ttmp = time*tscal
      it1 = nint(ttmp/(scale*scale))
      t1 = float(it1)*scale*scale
      it2 = nint( (ttmp-t1)/scale )
      t2 = float(it2)*scale
      it3 = nint(ttmp-t1-t2)
 
      write(iu) it1,it2,it3,nbod2,nleft

      return
      end     ! io_write_hdr.f
c---------------------------------------------------------------------------

