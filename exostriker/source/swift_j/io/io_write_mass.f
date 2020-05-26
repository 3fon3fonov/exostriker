c*************************************************************************
c                            IO_WRITE_MASS
c*************************************************************************
c write out masses
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Hal Levison 
c Date:    1/9/97
c Last revision: 3/19/97 

      subroutine io_write_mass(time,nbod,mass,oname,iu,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,iu
      real*8 mass(nbod),time
      character*80 oname,fopenstat

c...  Internals
      real*4 mass4(NTPMAX)
      real*4 ttmp
      integer*2 nbod2
      integer ierr,i,ldir,lfile
      character*80 dirname,filename

      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st,dirname,filename,ldir,lfile

c----
c...  Executable code 

c...  if first time through open file
      if(i1st.eq.0) then
         call io_splitname(oname,dirname,ldir,filename,lfile)
         call io_open(iu,dirname(1:ldir)//'mass.'//filename(1:lfile),
     &        fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_mass: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1st = 1
      else
        call io_open(iu,dirname(1:ldir)//'mass.'//filename(1:lfile),
     &        'append','UNFORMATTED',ierr)
      endif

      do i=1,nbod
         mass4(i) = mass(i)
      enddo
      nbod2 = nbod
      ttmp = time

      write(iu) ttmp,nbod2
      write(iu) (mass4(i),i=1,nbod)

      close(iu)
      return
      end      ! io_write_mass
c----------------------------------------------------------------------
