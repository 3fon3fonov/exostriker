c*************************************************************************
c                            MVITS_INIT_ROT.F
c*************************************************************************
c Calculates the rotation sin and cos for the test particles.
c
c             Input:
c                 ntp            ==>  number of test particles (int scalar)
c                 dtpl           ==>  time step of planet int (real scalar)
c                 nstep          ==> number of planet steps per big step
c                                          (int scalar)
c                 gm             ==> G * mass of the sun  (real scalar)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c             Output:
c                 cosrot,sinrot  ==>  cosine and sine of rotation angles
c                                      cosrot(i,1) is rotation after dtpl
c                                      cosrot(i,2) ==> dtpl*(nstep-1)
c                                        (real arrays)
c
c Authors:  Hal Levison 
c Date:    12/16/93
c Last revision: 

      subroutine mvits_init_rot(ntp,dtpl,nstep,gm,xht,yht,zht,
     &       vxht,vyht,vzht,cosrot,sinrot)


      include '../swift.inc'
      include 'mvits.inc'


c...  Inputs Only: 
      integer ntp,nstep
      real*8 dtpl,gm
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Output Only: 
      real*8 cosrot(NTPMAX,2),sinrot(NTPMAX,2)

c...  Internals
      integer i,ialpha
      real*8 a,e,inc,capom,omega,capm,meanmo,phi

c----
c...  Executable code 

      do i=1,ntp
         call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &          vzht(i),gm,ialpha,a,e,inc,capom,omega,capm)
         meanmo = sqrt(gm/(a**3))
         phi = dtpl*meanmo
         call orbel_scget(phi,sinrot(i,1),cosrot(i,1))
         phi = dtpl*meanmo*float(nstep-1)
         call orbel_scget(phi,sinrot(i,2),cosrot(i,2))
      enddo

      return
      end             !  mvits_init_rot
c----------------------------------------------------------
