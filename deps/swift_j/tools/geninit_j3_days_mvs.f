c GENINIT_J3 reads Jacobi orbital elements of nbod planets and generates
c initial position and velocity in Jacobi coords.
c This version outputs rpl and rhill.

c Last modified by Man Hoi Lee, Aug 16, 2003.

      include 'swift.inc'

      real*8 SMASSYR,MSUN,AU
      parameter (SMASSYR=TWOPI*TWOPI)
      parameter (MSUN=1.989d33, AU=1.49597892d13)

      integer nbod
      real*8 mass(NTPMAX),j2rp2,j4rp4
      real*8 xj(NTPMAX),yj(NTPMAX),zj(NTPMAX)
      real*8 vxj(NTPMAX),vyj(NTPMAX),vzj(NTPMAX)
      real*8 rpl(NTPMAX),rhill(NTPMAX)
      integer iflgchk
      logical*2 lclose

      real*8 mstar0,mpl0,frho3

      real*8 gm,a,ecc,inc,capom,omega,capm

      integer ialpha,i,iuflg

      character*80 outfile

 100  write (*,*) ' Units Menu:'
      write (*,*) '       0 ==> Solar masses, and AU'
      write (*,*) '       1 ==> AU, and Years'
      read (*,*) iuflg
      if ((iuflg.ne.0).and.(iuflg.ne.1)) goto 100
      write (*,*) ' Mass of central star (in solar mass):'
      read (*,*) mstar0
      write (*,*) ' Number of planets:'
      read (*,*) nbod
      write (*,*) ' f/rho^(1/3) for radius (rho in g/cm^3):'
      read (*,*) frho3
      write (*,*) ' Enter name of data file for symba5_j:'
      read (*,201) outfile
 201  format (a)

      lclose = .true.
      iflgchk = 0
      nbod = nbod + 1
      if (iuflg.eq.0) then
          mass(1) = mstar0
      else
          mass(1) = mstar0*(SMASSYR/(365.25d0**2.d0))
      endif
      j2rp2 = 0.d0
      j4rp4 = 0.d0
      xj(1) = 0.d0
      yj(1) = 0.d0
      zj(1) = 0.d0
      vxj(1) = 0.d0
      vyj(1) = 0.d0
      vzj(1) = 0.d0

      open (unit=12,file='geninit_j.out')
      write (12,*) ' m     a     e     i     omega capom M     rpl   rh'

      gm = mass(1)

      do 300 i=2,nbod
          write (*,*) ' Mass of planet',i-1,' (in solar mass):'
          read (*,*) mpl0
          if (iuflg.eq.0) then
              mass(i) = mpl0
          else
              mass(i) = mpl0*(SMASSYR/(365.25d0**2.d0))
          endif
          rpl(i) = frho3*(1.5d0*mpl0*MSUN/TWOPI)**0.3333333333d0/AU

          ialpha = -1
          gm = gm + mass(i)

          write (*,*)
     &        ' a (AU), e, i (deg), omega (deg), capom (deg), M (deg):'
          read (*,*) a,ecc,inc,omega,capom,capm

          inc = inc/DEGRAD
          capom = capom/DEGRAD
          omega = omega/DEGRAD
          capm = capm/DEGRAD

          rhill(i) = a*(mass(i)/(3.d0*mass(1)))**0.3333333333d0

          write (12,251) mass(i),a,ecc,inc,omega,capom,capm,
     &                   rpl(i),rhill(i)
 251      format (1p9e12.4)

          call ORBEL_EL2XV (gm,ialpha,a,ecc,inc,capom,omega,capm,
     &                      xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i))
 300  continue

      close (unit=12)

      call IO_DUMP_PL (outfile,nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                       lclose,iflgchk,rpl,j2rp2,j4rp4)


      stop
      end
