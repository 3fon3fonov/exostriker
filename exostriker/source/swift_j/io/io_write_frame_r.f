c*************************************************************************
c                            IO_WRITE_FRAME_R
c*************************************************************************
c write out a whole frame to an real*4 binary file.
c both massive and test particles
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  current position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  current velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  current part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  current velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      subroutine io_write_frame_r(time,nbod,ntp,mass,xh,yh,zh,vxh,
     &           vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,oname,
     &           iu,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,ntp,iu
      real*8 mass(nbod),time
      integer istat(NTPMAX,NSTAT)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      character*80 oname,fopenstat

c...  Internals
      integer i,id
      integer ialpha,ierr
      real*8 a,e,inc,capom,omega,capm
      real*8 gm
      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st

c----
c...  Executable code 

c...  if first time through open file
      if(i1st.eq.0) then
         call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_frame: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1st = 1
      else
        call io_open(iu,oname,'append','UNFORMATTED',ierr)
      endif

      call io_write_hdr_r(iu,time,nbod,ntp,istat)
      
c...  write out planets
      do i=2,nbod
         gm = mass(1)+mass(i)
         id = -1*i
 	 call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),gm,
     &          ialpha,a,e,inc,capom,omega,capm)
         call io_write_line_r(iu,id,a,e,inc,capom,omega,capm)
      enddo

c...  write out test particles
      gm = mass(1)
      do i=1,ntp
         if(istat(i,1).eq.0) then
            call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &          vzht(i),gm,ialpha,a,e,inc,capom,omega,capm)
            call io_write_line_r(iu,i,a,e,inc,capom,omega,capm)
         endif
      enddo

      close(iu)
      return
      end      ! io_write_frame_r
c----------------------------------------------------------------------
