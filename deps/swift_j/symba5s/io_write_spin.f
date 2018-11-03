c*************************************************************************
c                            IO_WRITE_SPIN
c*************************************************************************
c write out spin directions
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 sx,sy,sz      ==>  vector of spin direction (real arrays)
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
c Last revision: 7/18/06 MHL

      subroutine io_write_spin(time,nbod,sx,sy,sz,oname,iu,fopenstat)

      include '../swift.inc'
      include '../io/io.inc'

c...  Inputs: 
      integer nbod,iu
      real*8 sx(nbod),sy(nbod),sz(nbod),time
      character*80 oname,fopenstat

c...  Internals
      real*4 sx4(NTPMAX),sy4(NTPMAX),sz4(NTPMAX),ds4(NTPMAX)
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
         call io_open(iu,dirname(1:ldir)//'spin.'//filename(1:lfile),
     &        fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_spin: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1st = 1
      else
        call io_open(iu,dirname(1:ldir)//'spin.'//filename(1:lfile),
     &        'append','UNFORMATTED',ierr)
      endif

      do i=1,nbod
         sx4(i) = sx(i)
         sy4(i) = sy(i)
         sz4(i) = sz(i)
         ds4(i) = dsqrt(sx(i)**2+sy(i)**2+sz(i)**2) - 1.d0
      enddo
      nbod2 = nbod
      ttmp = time

      write(iu) ttmp,nbod2
      write(iu) (sx4(i),sy4(i),sz4(i),ds4(i),i=1,nbod)

      close(iu)
      return
      end      ! io_write_spin
c----------------------------------------------------------------------
