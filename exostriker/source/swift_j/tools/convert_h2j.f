c CONVERT_H2J converts astrocentric orbital elements in ascii files
c from follow_symba(2) to Jacobi elements. It assumes that there are
c only two planets. 

c Last modified by Man Hoi Lee, Aug. 13, 2003.

	include 'swift.inc'

	real*8 SMASSYR
	parameter (SMASSYR=TWOPI*TWOPI)

	integer nbod

	real*8 mass(NTPMAX)
	real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
	real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)
	real*8 xj(NTPMAX),yj(NTPMAX),zj(NTPMAX)
	real*8 vxj(NTPMAX),vyj(NTPMAX),vzj(NTPMAX)
	real*8 t

        real*8 a(3),e(3),inc(3),capom(3),omega(3),capm(3),peri,apo,gm

	integer i,lfile,ifoln(3),ialpha

	character*80 filecode

	write (*,*) ' Enter code of input/output files: '
	read (*,100) filecode
 100    format (a)
	write (*,*) ' Mass of central star (in solar mass):'
	read (*,*) mass(1)

	do i=1,80
	  if (filecode(i:i).ne.' ') lfile = i
	enddo

	nbod = 3

	mass(1) = mass(1)*SMASSYR
	xh(1) = 0.d0
	yh(1) = 0.d0
	zh(1) = 0.d0
	vxh(1) = 0.d0
	vyh(1) = 0.d0
	vzh(1) = 0.d0

	open (unit=10, file='follow_symba.'//filecode(1:lfile)//'.1',
     &        status='old')
	open (unit=20, file='follow_symba.'//filecode(1:lfile)//'.2',
     &        status='old')
	open (unit=30, file='convert_h2j.'//filecode(1:lfile)//'.1',
     &        status='new')
	open (unit=40, file='convert_h2j.'//filecode(1:lfile)//'.2',
     &        status='new')

 200	continue

	read (10,*,err=900) t,ifoln(2),a(2),e(2),inc(2),capom(2),
     &	                    omega(2),capm(2),peri,apo,mass(2)
	read (20,*) t,ifoln(3),a(3),e(3),inc(3),capom(3),
     &	            omega(3),capm(3),peri,apo,mass(3)

	do i=2,nbod
	    ialpha = -1
	    mass(i) = mass(i)*mass(1)
	    gm = mass(1) + mass(i)
	    inc(i) = inc(i)/DEGRAD
	    capom(i) = capom(i)/DEGRAD
	    omega(i) = omega(i)/DEGRAD
	    capm(i) = capm(i)/DEGRAD
	    call ORBEL_EL2XV(gm,ialpha,a(i),e(i),inc(i),capom(i),
     &	         omega(i),capm(i),xh(i),yh(i),zh(i),vxh(i),vyh(i),
     &	         vzh(i))
	enddo

	call COORD_H2J(nbod,mass,xh,yh,zh,vxh,vyh,vzh,xj,yj,zj,
     &                 vxj,vyj,vzj)

	gm = mass(1)
	do i=2,nbod
	    gm = gm + mass(i)
	    call ORBEL_XV2EL(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),
     &	         gm,ialpha,a(i),e(i),inc(i),capom(i),omega(i),capm(i))
	    inc(i) = inc(i)*DEGRAD
	    capom(i) = capom(i)*DEGRAD
	    omega(i) = omega(i)*DEGRAD
	    capm(i) = capm(i)*DEGRAD
	enddo

	write (30,800) t,ifoln(2),a(2),e(2),inc(2),capom(2),
     &	               omega(2),capm(2)
	write (40,800) t,ifoln(3),a(3),e(3),inc(3),capom(3),
     &	               omega(3),capm(3)
 800	format (1x,e15.7,1x,i3,1x,e13.5,1x,f7.5,4(1x,f9.4))

	goto 200

 900    close (unit=10)
	close (unit=20)
	close (unit=30)
	close (unit=40)

	stop
	end
