c ERROR computes the fractional change in energy and angular momentum using
c output in energy.out from swift_symba5.

c Last modified by Man Hoi Lee, Dec 7, 2001.

      implicit none

      real*8 t0,energy0,lx0,ly0,lz0,l0,t,energy,lx,ly,lz,de,dlx,dly,dlz
      character*80 infile

      write (*,*) ' Name of input file:'
      read (*,101) infile
 101  format (a)

      open (unit=10,file=infile,status='old')
      open (unit=20,file='error.out',status='new')

      read (10,*) t0,energy0,lx0,ly0,lz0

      l0 = sqrt(lx0**2 + ly0**2 + lz0**2)
      if (lx0.eq.0.d0) lx0 = l0
      if (ly0.eq.0.d0) ly0 = l0
      if (lz0.eq.0.d0) lz0 = l0

 200  continue
          read (10,*,end=300) t,energy,lx,ly,lz
          de = (energy - energy0)/energy0
          dlx = (lx - lx0)/lx0
          dly = (ly - ly0)/ly0
          dlz = (lz - lz0)/lz0
          write (20,201) t,de,dlx,dly,dlz
 201      format (1p5e12.4)
      goto 200

 300  close (unit=10)
      close (unit=20)

      stop
      end
