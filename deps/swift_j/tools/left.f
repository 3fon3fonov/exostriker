c**********************************************************************
c		      LEFT.F
c**********************************************************************
c

     
	include '/users/hal/SWIFT/swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX)
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT)
	integer nbod,ntp,i
	integer iflgchk,iact

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 tfrac,rstat(NTPMAX,NSTATR)

	real*8 rmin,rmax,rmaxu,qmin,gm,rplsq(NPLMAX),j2rp2,j4rp4
        real*8 a,e,inc,capom,omega,capm,peri,apo
        integer ialpha,nwhy(-4:NPLMAX)
        logical*2 lclose

	character*80 outfile,fopenstat

c Get parameters 
	call io_init_param('dump_param.dat',t0,tstop,dt,dtout,
     &       dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,
     &       outfile,fopenstat)

c Get data for sun and planetary orbits.
	call io_init_pl('dump_pl.dat',lclose,iflgchk,nbod,mass,xh,yh,
     &                  zh,vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	call io_init_tp('dump_tp.dat',ntp,xht,yht,zht,
     &            vxht,vyht,vzht,istat,rstat)


c 
       open(unit=7,file='left.out')

       do i=-4,nbod
          nwhy(i) = 0
       enddo

       tfrac = t0/tstop
       if(t0.lt.10000.0) then
          write(*,*) ' Time = ',t0
          write(7,*) ' Time = ',t0
       else
          write(*,10) t0
          write(7,10) t0
 10       format(' Time = ',1p1e13.5)
       endif
       write(*,*) ' Fraction done = ',tfrac
       write(7,*) ' Fraction done = ',tfrac

       iact = 0
       do i=1,ntp
          if(istat(i,1).eq.0) then
             iact = iact + 1
          endif
       enddo
       write(*,*) iact,' out of ',ntp,' still active '
       write(7,*) iact,' out of ',ntp,' still active '

       write(7,*) ' '
       write(7,*) ' Planet : '
       write(7,*) '       #          a                e             i '
       do i=2,nbod
          gm = mass(1) + mass(i)
          call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),
     &          vzh(i),gm,ialpha,a,e,inc,capom,omega,capm)
          inc = inc*180.0d0/PI
          write(7,1000) i,a,e,inc
 1000     format(5x,i4,3(5x,f10.4))
       enddo

       write(7,*) ' '
       write(7,*) ' Active Particles : '
       write(7,1011) 
 1011  format(8x,'#',10x,'a',17x,'e',12x,'i',13x,'q              Q ')

       do i=1,ntp
          if(istat(i,1).eq.0) then
             call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &          vzht(i),mass(1),ialpha,a,e,inc,capom,omega,capm)
             inc = inc*180.0d0/PI
             apo = a*(1.0d0+e)
             peri = a*(1.0d0-e)
             write(7,1001) i,a,e,inc,peri,apo
 1001        format(5x,i4,5(5x,f10.4))
          else
             nwhy(istat(i,2)) = nwhy(istat(i,2)) + 1
             nwhy(0) = nwhy(0) + 1
          endif
       enddo

       write(7,*) ' '
       write(7,*) ' Discarded Particles : ',nwhy(0),' out of ',ntp

       write(7,*) ' '
c       write(7,*) '     #     why?   Last pl       a         q        Q'
       write(7,1010) 
 1010  format(6x,'#',5x,'why?  Last pl',6x,'a',11x,'q',11x,'Q',13x,'i')
       write(7,*) '               encountered '
       do i=1,ntp
          if(istat(i,1).ne.0) then
             call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &          vzht(i),mass(1),ialpha,a,e,inc,capom,omega,capm)
             apo = a*(1.0d0+e)
             peri = a*(1.0d0-e)
             write(7,1008) i,istat(i,2),istat(i,3),a,peri,apo,
     &            inc*DEGRAD
 1008        format(3(3x,i4),5(3x,f10.4))
          endif
       enddo

       write(7,*) ' '
       write(7,*) ' Fate of Discarded Particles : '

       write(7,1006) nwhy(-1)
 1006  format('   ',i3,' have istat(i,2) = -1',
     &      '  ==> Danby did not converge ')
       write(7,1002) nwhy(-2)
 1002  format('   ',i3,' have istat(i,2) = -2',
     &      '  ==> Ejected from the system ')
       write(7,1003) nwhy(-3)
 1003  format('   ',i3,' have istat(i,2) = -3',
     &      '  ==> Too far from the Sun ')
       write(7,1007) nwhy(-4)
 1007  format('   ',i3,' have istat(i,2) = -4',
     &      '  ==> Too small a perihelion ')
       write(7,1004) nwhy(1)
 1004  format('   ',i3,' have istat(i,2) =  1',
     &      '  ==> Too close to the Sun ')

       do i=2,nbod
          write(7,1005) nwhy(i),i,i
 1005     format('   ',i3,' have istat(i,2) = ',i2,
     &         '  ==> Too close to planet',i2)
       enddo


       write(7,*) '   --- '
       write(7,1009) nwhy(0),ntp
 1009  format('   ',i3,'/',i4)

       stop
       end                      ! left.f
c---------------------------------------------------------------------






