C**************************************************************************
C	    		        BS_DER
C**************************************************************************
c This is the subroutine that calculates the derivatives of the independant var
c
c             Input:
c              nbod  ==> number of planets  (int scalar)
c              ntp   ==> number of test particles  (int scalar)
c              mass  ==>  mass of bodies (real array)
c      j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c             istat  ==>  status of the test paricles
c                           (2d integer array)
c                            istat(i,1) = 0 ==> active:  = 1 not
c                            istat(i,2) = -1 ==> Danby did not work
c 	        ybs  ==> values dependent variables  (real array)
c
c             Output:
c 	         dy  ==> derivatives of the independant var (real array)
c
c Remarks:  This used TU4 routines !!  
c Authors:  Hal Levison
c Date:    5/17/93
c Last revision: 2/24/94

      subroutine bs_der(ntp,nbod,mass,j2rp2,j4rp4,istat,ybs,dy)

      include '../swift.inc'
      include 'bs.inc'

c...  Inputs Only: 
      integer nbod,ntp
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 ybs(6,(NTPMAX+NPLMAX))

c...  Input and Outputs
      integer istat(NTPMAX,NSTAT)

c...  Output
      real*8 dy(6,(NTPMAX+NPLMAX))

c...  Internals
      integer i,j
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      real*8 axb(NPLMAX),ayb(NPLMAX),azb(NPLMAX)
      real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
      real*8 axbt(NTPMAX),aybt(NTPMAX),azbt(NTPMAX)

c----
c...  Executable code 

c...  move things so that I can deal with it
      do i=1,nbod
         xb(i) = ybs(1,i)
         yb(i) = ybs(2,i)
         zb(i) = ybs(3,i)
         vxb(i) = ybs(4,i)
         vyb(i) = ybs(5,i)
         vzb(i) = ybs(6,i)
      enddo

      do i=1,ntp
         j = i + nbod
         xbt(i) = ybs(1,j)
         ybt(i) = ybs(2,j)
         zbt(i) = ybs(3,j)
         vxbt(i) = ybs(4,j)
         vybt(i) = ybs(5,j)
         vzbt(i) = ybs(6,j)
      enddo

      call tu4_getaccb(nbod,mass,j2rp2,j4rp4,xb,yb,zb,axb,ayb,azb)
      call tu4_getaccb_tp(nbod,mass,j2rp2,j4rp4,xb,yb,zb,
     &     ntp,xbt,ybt,zbt,istat,axbt,aybt,azbt)


c.... moves things intp dy array
      do i=1,nbod
         dy(1,i) = ybs(4,i)
         dy(2,i) = ybs(5,i)
         dy(3,i) = ybs(6,i)
         dy(4,i) = axb(i)
         dy(5,i) = ayb(i)
         dy(6,i) = azb(i)
      enddo

      do i=1,ntp
         j = i + nbod
         dy(1,j) = ybs(4,j)
         dy(2,j) = ybs(5,j)
         dy(3,j) = ybs(6,j)
         dy(4,j) = axbt(i)
         dy(5,j) = aybt(i)
         dy(6,j) = azbt(i)
      enddo

      return
      end     ! bs_der
c-------------------------------------------------------------------------------
