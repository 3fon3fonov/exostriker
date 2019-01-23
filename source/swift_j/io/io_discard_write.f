c*************************************************************************
c                            IO_DISCARD_WRITE
c*************************************************************************
c checks to see if a particle was removed since te last time the subroutine
c was called.  If so write out the position and vol of tp and planets
c
c             Input:
c                 init          ==>  initiize flag if = 0 initialize and return
c                                                     = 1 run through 
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp           ==>  number of test bodies (int scalar)
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
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,1) time of discard.
c                                      rstat(i,2) closest approach to a planet
c                                          as determined by encounter routines.
c                 iu              ==> unit number to write to
c                 rname           ==> output file name (character string) 
c
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c                 nleft           ==>  number of active test bodies(int scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    5/7/93
c Last revision: 7/11/94

      subroutine io_discard_write(init,time,nbod,ntp,xh,yh,zh,
     &     vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &     rstat,iu,rname,fopenstat,nleft)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer init,nbod,ntp,iu
      real*8 time
      real*8 rstat(NTPMAX,NSTATR)
      integer istat(NTPMAX,NSTAT)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      character*(*) rname,fopenstat

c...  Outputs: 
      integer nleft

c...  Internals
      integer istold(NTPMAX),i,in,irem,j
      integer iwrite,ierr,izero(NSTAT)
      real*8 rzero(NSTATR)

      save istold,izero,iwrite,rzero

c----
c...  Executable code 


c...  if  initialize flag=0
      if(init.eq.0) then

         do i=1,ntp
            istold(i) = istat(i,1)
         enddo

         do i=1,NSTAT
           izero(i) = 0
         enddo
         do i=1,NSTATR
           rzero(i) = 0.0d0
         enddo
         iwrite = 0
         return      ! NOTE !!!!!!
      endif

c...  if  initialize flag=1

      irem = 0
      nleft = 0
      do i=1,ntp
         if(istat(i,1).eq.0) then
            nleft = nleft + 1
         endif
         if(istold(i).ne.istat(i,1)) then
            if(irem.eq.0) then
               if(iwrite.eq.0) then
                  call io_open(iu,rname,fopenstat,'FORMATTED',ierr)

c...              if there was an error and fopenstat='append' then
c...              try to open as new
                  if(ierr.ne.0) then  
                     if( (fopenstat(1:6).eq.'append') .or. 
     &                    (fopenstat(1:6).eq.'APPEND') ) then
                        call io_open(iu,rname,'new','FORMATTED',ierr)
                     endif
                  endif

                  if(ierr.ne.0) then
                     write(*,*) ' SWIFT ERROR: in io_discard_write: '
                     write(*,*) '    Could not open discard output file'
                     call util_exit(1)
                  endif
               else
                  call io_open(iu,rname,'append','FORMATTED',ierr)
               endif
               write(iu,1000) time
 1000          format(1x,1p1e23.16)
            endif
            iwrite = 1
            irem = irem + 1
            write(iu,2000) i,(istat(i,j),j=1,NSTAT)
 2000       format(3x,i5,100(1x,i5))
            write(iu,3000) (rstat(i,j),j=1,NSTATR)
            write(iu,3000) xht(i),yht(i),zht(i)
            write(iu,3000) vxht(i),vyht(i),vzht(i)
 3000       format(4(1p1e23.16,1x))
         endif
         istold(i) = istat(i,1)
      enddo

      if(irem.ne.0) then
         do i=2,nbod       ! don't do the Sun
            in = -1*i
            write(iu,2000) in,(izero(j),j=1,NSTAT)
            write(iu,3000) (rzero(j),j=1,NSTATR)
            write(iu,3000) xh(i),yh(i),zh(i)
            write(iu,3000) vxh(i),vyh(i),vzh(i)
         enddo
         in = 0
         write(iu,2000) in,(izero(j),j=1,NSTAT)
      endif

      close(iu)

      return   
      end     !io_discard_write
