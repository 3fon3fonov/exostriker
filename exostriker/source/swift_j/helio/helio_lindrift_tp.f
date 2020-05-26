c*************************************************************************
c                            HELIO_LINDRIFT_TP.F
c*************************************************************************
c This subroutine takes a linear drift due to mometum of Sun
c
c             Input:
c                 ntp              ==>  number of TPs (int scalar)
c                 ptx,pty,ptz      ==> momentum of sun
c                                          (real scalars)
c                 dt               ==>  time step
c                 istat            ==>  status of the test paricles
c                                         (integer array)
c                                         istat(i) = 0 ==> active:  = 1 not
c                                        NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c                 xht,yht,zht      ==>  initial position in helio coord 
c                                          (real arrays)
c             Output:
c                 xht,yht,zht      ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Bases on Martin's code h2.f
c Authors:  Hal Levison 
c Date:    11/14/96
c Last revision: 11/15/96

      subroutine helio_lindrift_tp(ntp,ptx,pty,ptz,dt,
     &     istat,xht,yht,zht)

      include '../swift.inc'

c...  Inputs Only: 
      integer ntp,istat(NTPMAX)
      real*8 ptx,pty,ptz,dt

c...  Inputs and Outputs:
      real*8 xht(ntp),yht(ntp),zht(ntp)

c...  Internals:
      integer n

c----
c...  Executable code 

      do n=1,ntp
         if(istat(n).eq.0) then
            xht(n) = xht(n) + ptx*dt
            yht(n) = yht(n) + pty*dt
            zht(n) = zht(n) + ptz*dt
         endif
      enddo

      return
      end       ! helio_lindrift_tp
c---------------------------------------------------
