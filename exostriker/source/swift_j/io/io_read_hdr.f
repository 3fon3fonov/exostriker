c*************************************************************************
c                            IO_READ_HDR
c*************************************************************************
c read in header part of the integer*2 binary file
c
c             Input:
c                 iu            ==> unit number to write to
c             Output:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nleft         ==>  number of active tp (int scalar)
c
c             Returns:
c               io_read_hdr     ==>   =0 read ok
c                                    !=0 read failed is set to iostat variable
c Remarks: 
c Authors:  Hal Levison 
c Date:    3/1/93 (HAPPY BIRTHDAY TO ME!!!)
c Last revision: 2/21/94

      integer function io_read_hdr(iu,time,nbod,nleft) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output
      integer nbod,nleft
      real*8 time

c...  Internals
      real*8 t3,t1,t2,scale,tscal
      integer*2 it1,it2,it3,nleft2,nbod2
      integer ierr

      data scale,tscal/32767.0,100.0/

c----
c...  Executable code 


      read(iu,iostat=ierr) it1,it2,it3,nbod2,nleft2
      io_read_hdr = ierr
      if(ierr.ne.0) then
         return
      endif

      nbod = nbod2
      nleft = nleft2

      t1 = float(it1)*scale*scale
      t2 = float(it2)*scale
      t3 = float(it3)
      time = (t1 + t2 + t3)/tscal

      return
      end     ! io_read_hdr.f
c---------------------------------------------------------------------------

