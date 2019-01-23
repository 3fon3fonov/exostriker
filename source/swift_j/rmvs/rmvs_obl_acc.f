c*************************************************************************
c                            RMVS_OBL_ACC.F
c*************************************************************************
c This subroutine adds the differential accel due to J2 and J4
c
c             Input:
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 ipl            ==>  the planet that is currently in the 
c                                      center (integer scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xpl,ypl,zpl    ==>  massive part position at 
c                                       (real arrays)
c              aoblx,aobly,aoblz ==> acceleration of the Sun on the central pl
c                                     due to J2 anf J4
c                                      (real scalars)
c                 xpt,ypt,zpt    ==>  initial part position in planet-coord 
c                                      (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c              axpt,aypt,azpt    ==>  accel on tp w/o J2 and J4
c                                      (real arrays)
c
c             Output:
c                 axpt,aypt,azpt ==>  accel on tp WITH J2 and J4 added.
c                                       (real arrays)
c
c Remarks: Taken from step_kdk_tp.f
c Authors:  Hal Levison 
c Date:    2/24/94
c Last revision: 

      subroutine rmvs_obl_acc(nbod,ntp,ipl,mass,j2rp2,j4rp4,xpl,ypl,
     &     zpl,aoblx,aobly,aoblz,xpt,ypt,zpt,istat,axpt,aypt,azpt)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,ipl
      integer istat(NTPMAX,NSTAT)
      real*8 mass(nbod),j2rp2,j4rp4  
      real*8 xpl(NPLMAX),ypl(NPLMAX),zpl(NPLMAX)
      real*8 xpt(ntp),ypt(ntp),zpt(ntp)
      real*8 aoblx,aobly,aoblz

c...  Inputs and Outputs:
      real*8 axpt(ntp),aypt(ntp),azpt(ntp)

c...  Internals:
      integer i
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 irht(NPLMAX)
      real*8 aoblxt(NTPMAX),aoblyt(NTPMAX),aoblzt(NTPMAX)


c----
c...  Executable code 

c...  Do we need to do this?
      if(j2rp2.eq.0.0d0) then
         return                  !!!!!!! NOTE 
      endif

c...  first get barycentric accel 
      do i=1,ntp
         xht(i) = xpt(i)-xpl(ipl)
         yht(i) = ypt(i)-ypl(ipl)
         zht(i) = zpt(i)-zpl(ipl)
         irht(i) = 1.0d0/sqrt( xht(i)**2 + yht(i)**2 + zht(i)**2 )
      enddo

      call obl_acc_tp(ntp,istat,mass(ipl),j2rp2,j4rp4,xht,yht,zht,
     &     irht,aoblxt,aoblyt,aoblzt)

      do i = 1,ntp
         axpt(i) = axpt(i) + aoblxt(i) - aoblx
         aypt(i) = aypt(i) + aoblyt(i) - aobly
         azpt(i) = azpt(i) + aoblzt(i) - aoblz
      enddo

      return
      end     ! rmvs_obl_acc
c------------------------------------------------------------------------
