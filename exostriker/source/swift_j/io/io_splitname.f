c*************************************************************************
c                            IO_SPLITNAME
c*************************************************************************
c splits the directory from the filename in a string
c
c             Input:
c                 oname        ==> string with the full path (character*80)
c   
c             Output:
c                 dirname     ==> string with the path (character*80)
c                 ldir        ==> length of dirname (integer scalar)
c                 filename     ==> string with the file name (character*80)
c                 lfile        ==> length of filename (integer scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:   3/19/97 
c Last revision: 

      subroutine io_splitname(oname,dirname,ldir,filename,lfile)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      character*80 oname

c...  Outputs:
      integer ldir,lfile
      character*80 dirname,filename

c...  Internals
      integer i,il,is

c----
c...  Executable code 

c... Find the last character
      il = 0
      do i=1,80
         if(oname(i:i).eq.' ') then
            il = i - 1
            goto 99
         endif
      enddo
 99   continue

      if(il.eq.0) then
         il = 80
      endif

c... Find the last /
      is = 0
      do i=1,il
         if(oname(i:i).eq.'/') then
            is = i
         endif
      enddo

      if(is.eq.0) then          ! there is no path
         do i=1,il
            filename(i:i) = oname(i:i)
         enddo
         lfile = il
         write(dirname,1000)
 1000    format('./')
         ldir = 2
      else
         do i=1,is
            dirname(i:i) = oname(i:i)
         enddo
         ldir = is
         do i=is+1,il
            filename(i-is:i-is) = oname(i:i)
         enddo
         lfile = il - is
      endif

      return
      end            ! io_splitname
c---------------------------------------------------------

