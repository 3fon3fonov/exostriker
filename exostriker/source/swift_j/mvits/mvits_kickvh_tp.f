c*************************************************************************
c                        MVITS_KICKVH_TP.F
c*************************************************************************
c To kick the velocity components vxh(*) by axh(*)*dt for test particles 
c
c             Input:
c                 ntp          ==>  number of bodies (int scalar)
c                 nbod         ==>  number of massive bodies (int scalar)
c                 vxh,vyh,vzh  ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 axh,ayh,azh  ==>  acceleration in helio coord
c                                    (2d real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c                        dtpli  ==>  time step for each planet (real array)
c                        idpl   ==>  link list of planet order  (real array)
c             Output:
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                    (real arrays)
c
c Remarks: kick_tp
c Authors:  Hal Levison 
c Date:    12/16/93
c Last revision: 

      subroutine mvits_kickvh_tp(nbod,ntp,vxht,vyht,vzht,axht,
     &     ayht,azht,istat,dtpli,idpl) 

      include '../swift.inc'
      include 'mvits.inc'

c...  Inputs Only: 
      integer ntp,istat(ntp),nbod,idpl(NPLMAX)
      real*8 axht(NTPMAX,NPLMAX),ayht(NTPMAX,NPLMAX),azht(NTPMAX,NPLMAX)
      real*8 dtpli(NTPMAX)

c...   Inputs and Output:
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
      integer n,i

c----
c...  Executable code 

      do n= 1, ntp
         if(istat(n).eq.0) then
            do i=2,nbod
               vxht(n) = vxht(n) + axht(n,i)*dtpli(idpl(i))/2.0
               vyht(n) = vyht(n) + ayht(n,i)*dtpli(idpl(i))/2.0
               vzht(n) = vzht(n) + azht(n,i)*dtpli(idpl(i))/2.0
            enddo
         endif
      enddo

      return
      end                       ! mvits_kickvh_tp
c-----------------------------------------------------------------------------
