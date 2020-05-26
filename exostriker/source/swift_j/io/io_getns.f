c**********************************************************************
c			IO_GETNS.F
c**********************************************************************
c Determines that number of istat variables there are in a tp.in file
c
c             Input:
c                 iu            ==> unit number (integer)
c
c             Output:
c                 ns            ==>  number of istat variable (int scalar)
c
c
c Remarks: 
c Authors:  Hal Levison
c Date:    10/1/96
c Last revision:  3/18/97

      subroutine io_getns(iu,ns)

      include '../swift.inc'
      include 'io.inc'

c...  Input
      integer iu

c...  Output
      integer ns
      
c...  Internal
      character*1024 line
      integer i,i1,ib
      real*8 xht,yht,zht    
      real*8 vxht,vyht,vzht

c-----
c...  Executable code      

c...  get the irrelavant stuff
      read(7,*) xht,yht,zht    
      read(7,*) vxht,vyht,vzht

      ns = 0
      do while(.true.)
         read(7,fmt='(a)') line

c...     if there are `.' then it is not a istat line. Therefore leave
         do i = 1,1024
            if(line(i:i).eq.'.') goto 99
         enddo


c...     Find the first non-blank character
         i1 = 0
         do i=1,1024
            if( (i1.eq.0) .and. (line(i:i).ne.' ') ) then
               i1 = i
            endif
         enddo

         ib = 1
         do i=i1+1,1024
            if( (ib.eq.1) .and. (line(i:i).eq.' ') ) then
               ns = ns + 1
            endif
            if (line(i:i).eq.' ') then
               ib = 0
            else
               ib = 1
            endif
         enddo

      enddo
 99   continue

      rewind(7)

      return
      end     ! io_getns
c----------------------------------------------------


