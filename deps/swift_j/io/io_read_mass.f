c*************************************************************************
c                            IO_READ_MASS
c*************************************************************************
c read in the mass file.
c
c             Output:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 iu              ==> unit number to read to
c
c             Returns:
c               io_read_mass     ==>   =0 read ok
c                                    !=0 read failed is set to iostat variable
c
c Remarks: Based on io_read_frame
c Authors:  Hal Levison 
c Date:    1/9/97
c Last revision: 

      integer function io_read_mass(time,nbod,mass,iu)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Outputs
      integer nbod
      real*8 mass(nbod),time

c...  Internals
      real*4 mass4(NTPMAX)
      real*4 ttmp
      integer*2 nbod2
      integer i,ierr

c----
c...  Executable code 

      read(iu,iostat=ierr) ttmp,nbod2
      io_read_mass = ierr
      if(ierr.ne.0) then
         return
      endif

      read(iu,iostat=ierr) (mass4(i),i=1,nbod2)
      io_read_mass = ierr
      if(ierr.ne.0) then
         return
      endif

      do i=1,nbod
         mass(i) = mass4(i)
      enddo
      nbod = nbod2
      time = ttmp

      return
      end      ! io_read_mass
c----------------------------------------------------------------------
