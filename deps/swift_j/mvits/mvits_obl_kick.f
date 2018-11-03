c*************************************************************************
c                        MVITS_OBL_KICK.F
c*************************************************************************
c This subroutine kicks the tp's for J2 and J4 only
c
c             Input:
c                 i1st         ==>  = 0 if first step; = 1 not (int scalar)
c                  nbod        ==>  number of massive bodies (int scalor)
c                  ntp         ==>  number of tp bodies (int scalor)
c                  mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4  ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                          (real scalars) NOT USED!!!
c                  xh,yh,zh    ==>  massive part position in helio coord 
c                                     (real arrays)
c                  xht,yht,zht ==>  test part position in heliocentric coord 
c                                     (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                  dt            ==>  time step
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c Author:  Hal Levison  
c Date:    6/11/98
c Last revision: 

      subroutine mvits_obl_kick(i1st,nbod,ntp,mass,j2rp2,j4rp4,xh,
     &     yh,zh,xht,yht,zht,vxht,vyht,vzht,istat,dt)

      include '../swift.inc'
      include 'mvits.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX),i1st
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),dt
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX),j2rp2,j4rp4

c...  Inputs and Output 
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
      integer i,j
      real*8 ir3h(NPLMAX),irh(NPLMAX)
      real*8 ir3ht(NPLMAX),irht(NPLMAX)
      real*8 fac,dx,dy,dz,rji2,irij3
      real*8 axh0(NPLMAX),ayh0(NPLMAX),azh0(NPLMAX)
      real*8 aoblx(NPLMAX),aobly(NPLMAX),aoblz(NPLMAX)
      real*8 aoblxt(NTPMAX),aoblyt(NTPMAX),aoblzt(NTPMAX)

      save aoblx,aobly,aoblz,aoblxt,aoblyt,aoblzt

c----
c...  Executable code 

      if(j2rp2.ne.0.0d0) then
         if(i1st.eq.0) then 
c...        get the r^-3's  for the planets
            call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
            call getacch_ir3(ntp,1,xht,yht,zht,ir3ht,irht)

            call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &           aoblx,aobly,aoblz)
            call obl_acc_tp(ntp,istat,mass(1),j2rp2,j4rp4,xht,
     &           yht,zht,irht,aoblxt,aoblyt,aoblzt)
         endif
         do i = 1,ntp
            if(istat(i).eq.0) then
               vxht(i) = vxht(i) + (aoblxt(i)-aoblx(1))*dt
               vyht(i) = vyht(i) + (aoblyt(i)-aobly(1))*dt
               vzht(i) = vzht(i) + (aoblzt(i)-aoblz(1))*dt
            endif
         enddo
      endif

      return
      end


