c*************************************************************************
c                            RMVS_ELOG.F
c*************************************************************************
c This subroutine keeps a log of the close encounters using istat
c
c             Input:
c                 ntp           ==>  number of massive bodies (int scalar)
c                 icflg         ==> ecounters? = 1 Yes, in outer region only
c                                              = -1 in inner region
c                                              =  0 No (integer scalar)  
c                 ienci         ==> ienci(j) = 0 if tp j not involved in enc 
c                                   in inner region: = planet# if it is. 
c                                     (integer array)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,NSTATP+n-1) is number of
c                                          inner encounters that tp i had
c                                          with planet n 
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,NSTATP+n-1) is number of
c                                          inner encounters that tp i had
c                                          with planet n (updated).
c
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/23/94
c Last revision: 

      subroutine rmvs_elog(ntp,icflg,ienci,istat)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer ntp,icflg,ienci(NTPMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)

c...  Internals
      integer i,i1st,iencio(NTPMAX),np

      data i1st/0/

      save i1st,iencio

c----
c...  Executable code 

c...  set up the old storage location
      if(i1st.eq.0) then
         do i=1,ntp
            iencio(i) = 0
         enddo
         i1st = 1
      endif

c...  check for new encounters
      if(icflg.eq.-1) then
         do i=1,ntp
            if( (iencio(i).eq.0) .and. (ienci(i).ne.0) ) then 
               np = NSTATP + ienci(i) - 1
               istat(i,np) = istat(i,np) + 1
            endif
         enddo
      endif

c.... save the current values
      do i=1,ntp
         iencio(i) = ienci(i) 
      enddo

      return
      end                       ! rmvs_elog
c------------------------------------------------------

