c*************************************************************************
c                            UTIL_EXIT.F
c*************************************************************************
c Exits program
c
c             Input:
c                 iflg          ==>  status of exit
c                                       = 0 if normal exit
c                                       = 1 if exit because error
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    8/6/93
c Last revision: MD : change calc. of rhil Apr. 25

      subroutine util_exit(iflg)

      include '../swift.inc'

c...  Inputs: 
      integer iflg


c-----
c...  Executable code 

      write(*,*) ' '

      if(iflg.eq.0) then
        write(*,1000) VER_NUM 
 1000   format('Normal termination of SWIFT (version ',f3.1,')')
      else
        write(*,2000) VER_NUM 
 2000   format('Terminating SWIFT (version',f3.1,') due to ERROR!!! ')
      endif

      write(*,*) '----------------------------------------------------'

      stop
      end  ! util_exit

c---------------------------------------------------


