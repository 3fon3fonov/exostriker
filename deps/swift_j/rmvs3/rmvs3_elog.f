c*************************************************************************
c                            RMVS3_ELOG.F
c*************************************************************************
c This subroutine keeps a log of the close encounters using istat
c
c             Input:
c                 ntp           ==>  number of massive bodies (int scalar)
c                 icflg         ==> ecounters? = 1 Yes, in inner region only
c                                              =  0 No (integer scalar)  
c                 ienc          ==> ienc(j) = 0 if tp j not involved in enc 
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
c Remarks: BAsed on rmvs_elog
c Authors:  Hal Levison 
c Date:    7/10/96
c Last revision: 

      subroutine rmvs3_elog(ntp,icflg,ienc,istat)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer ntp,icflg,ienc(NTPMAX)

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
      if(icflg.eq.1) then
         do i=1,ntp
            if( (iencio(i).eq.0) .and. (ienc(i).ne.0) ) then 
               np = NSTATP + ienc(i) - 1
               istat(i,np) = istat(i,np) + 1
            endif
         enddo
      endif

c.... save the current values
      do i=1,ntp
         iencio(i) = ienc(i) 
      enddo

      return
      end                       ! rmvs3_elog
c------------------------------------------------------

