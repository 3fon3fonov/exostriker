c*************************************************************************
c                            IO_READ_HDR_R
c*************************************************************************
c read in header part of the real*4 file
c
c             Input:
c                 iu            ==> unit number to write to
c             Output:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nleft         ==>  number of active tp (int scalar)
c
c             Returns:
c               io_read_hdr_r     ==>   =0 read ok
c                                    !=0 read failed is set to iostat variable
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      integer function io_read_hdr_r(iu,time,nbod,nleft) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output
      integer nbod,nleft
      real*8 time

c...  Internals
      real*4 ttmp
      integer*2 nleft2,nbod2
      integer ierr

c----
c...  Executable code 


      read(iu,iostat=ierr) ttmp,nbod2,nleft2
      io_read_hdr_r = ierr
      if(ierr.ne.0) then
         return
      endif

      nbod = nbod2
      nleft = nleft2
      time = ttmp

      return
      end     ! io_read_hdr_r.f
c---------------------------------------------------------------------------

