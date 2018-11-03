c*************************************************************************
c                        LYAP2_ACC_TP.F
c*************************************************************************
c This subroutine calculates the delta acceleration on the test particles
c for difference equations
c             Input:
c                  nbod        ==>  number of massive bodies (int scalor)
c                  ntp         ==>  number of tp bodies (int scalor)
c                  mass        ==>  mass of bodies (real array)
c                  j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                  xh,yh,zh    ==>  massive part position in helio coord 
c                                     (real arrays)
c                  xht,yht,zht ==>  test part position in heliocentric coord 
c                                     (real arrays)
c             dxht,dyht,dzht    ==>  separation in position
c                                     (real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c               daxht,dayht,dazht ==>  tp delta acceleration in helio coord 
c                                   (real arrays)
c
c Comments: Based on getacch_tp
c Author:  Hal Levison  
c Date:    7/11/95
c Last revision: 

      subroutine lyap2_acc_tp(nbod,ntp,mass,j2rp2,j4rp4,xh,yh,zh,
     &     xht,yht,zht,dxht,dyht,dzht,istat,daxht,dayht,dazht)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX)
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX),j2rp2,j4rp4
      real*8 dxht(NTPMAX),dyht(NTPMAX),dzht(NTPMAX)

c...  Outputs:
      real*8 daxht(NTPMAX),dayht(NTPMAX),dazht(NTPMAX)
                
c...  Internals:
      integer i,j
      real*8 rx,ry,rz,rji2,irij3,irij5
      real*8 t1,t2,t3

c----
c...  Executable code 

c...  the 0th term is zero

c...  the first terms are 0

c...  the second terms are 0

c...  now the third terms

      do i=1,ntp
         daxht(i) = 0.0
         dayht(i) = 0.0
         dazht(i) = 0.0
      enddo

      do j=1,ntp
         if(istat(j).eq.0) then
            do i=2,nbod

               rx = xht(j) - xh(i)
               ry = yht(j) - yh(i)
               rz = zht(j) - zh(i)
               rji2 = rx*rx + ry*ry + rz*rz
               irij3 = 1.0d0/(rji2*sqrt(rji2))
               irij5 = irij3/(rji2)

c...           x 
               t1 = dxht(j)*(-irij3 + 3.0*rx*rx*irij5)
               t2 = 3.0*rx*ry*irij5*dyht(j)
               t3 = 3.0*rx*rz*irij5*dzht(j)
               daxht(j) = daxht(j) + mass(i)*(t1+t2+t3)

c...           y
               t1 = 3.0*rx*ry*irij5*dxht(j)
               t2 = (-irij3 + 3.0*ry*ry*irij5)*dyht(j)
               t3 = 3.0*ry*rz*irij5*dzht(j)
               dayht(j) = dayht(j) + mass(i)*(t1+t2+t3)

c...           z
               t1 = 3.0*rx*rz*irij5*dxht(j)
               t2 = 3.0*ry*rz*irij5*dyht(j)
               t3 = (-irij3 + 3.0*rz*rz*irij5)*dzht(j)
               dazht(j) = dazht(j) + mass(i)*(t1+t2+t3)

            enddo
         endif
      enddo



c...  Now do j2 and j4 stuff
c...    Not included in this version !!!!!!

      return
      end      ! lyap2_acc_tp

c---------------------------------------------------------------------




