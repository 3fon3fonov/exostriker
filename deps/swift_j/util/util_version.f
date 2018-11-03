c*************************************************************************
c                            UTIL_VERSION.F
c*************************************************************************
c Prints version of Swift
c
c             NO IO!!!!!!
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/21/94
c Last revision: 

      subroutine util_version

      include '../swift.inc'

c-----
c...  Executable code 

      
      write(*,1000) VER_NUM
 1000 format('************* SWIFT: Version ',f3.1,' *************')

      write(*,*) ' '
      write(*,*) 'Authors:'
      write(*,*) '   Martin Duncan: Queen''s University '
      write(*,*) '   Hal Levison: Southwest Research Institute '
      write(*,*) ' '
      write(*,*) ' Please address any comments or questions to:'
      write(*,*) '   Hal Levison '
      write(*,*) '   Geophysical, Astrophysical, & Planetary Sciences'
      write(*,*) '   Southwest Research Institute'
      write(*,*) '   1050 Walnut St.'
      write(*,*) '   Suite 429 '
      write(*,*) '   Boulder, Co  80302 '
      write(*,*) '   (303) 546-0290 '
      write(*,*) '   Fax: (303) 546-9687 '
      write(*,*) '   (D)  swri::levison '
      write(*,*) '   (I)  hal@gort.space.swri.edu '
      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '

      return

      end  ! util_exit

c---------------------------------------------------
