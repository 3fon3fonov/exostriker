c converts binary file to ascii file

	include 'swift.inc'

	real*8 mass(NTPMAX),dr
	real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
	real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

	integer nbodm,nbod,nbod0,ierr,ifol,istep
	integer iflgchk,iu,i,id,nleft,ium
        integer io_read_hdr,io_read_line,io_read_mass
        integer io_read_hdr_r,io_read_line_r

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tmax

	real*8 rmin,rmax,rmaxu,qmin,rpl(NTPMAX),rhill(NTPMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4
        real*8 peri,apo,tg

        integer plist(NTPMAX),ifoln

	character*80 outfile,inparfile,inplfile,fopenstat

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl_symba(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

        iu = 20
        ium = 30

        dr = 180.0/PI

        if(btest(iflgchk,0)) then
           write(*,*) ' Reading an integer*2 binary file '
        else if(btest(iflgchk,1)) then
           write(*,*) ' Reading an real*4 binary file '
        else
           write(*,*) ' ERROR: no binary file format specified '
           write(*,*) '        in param file '
           stop
        endif

        write(*,*) ' Input the particle number to follow '
        read(*,*) ifol
        ifol = abs(ifol)
        write(*,*) ' Following particle ',ifol

        open(unit=iu, file=outfile, status='old',form='unformatted')
        open(unit=ium, file='mass.'//outfile, status='old',
     *       form='unformatted')
        open(unit=7,file='follow_symba.out')

        nbod0 = nbod
        do i=1,nbod
           plist(i) = i
        enddo
        call follow_plist(0,ifol,ifoln,plist,nbod0,tg,tstop)

        ifoln = ifol

        write(*,*) '1  2 3 4  5    6     7    8    9   10  11 '
        write(*,*) 't,id,a,e,inc,capom,omega,capm,peri,apo, M '

        tmax = t0
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
                ierr = io_read_hdr(iu,t,nbod,nleft) 
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft) 
             endif
             if(ierr.ne.0) goto 2

             ierr = io_read_mass(t,nbodm,mass,ium)

             if(nbodm.ne.nbod) then
                write(*,*) ' Error 1:',nbod,nbodm
                stop
             endif

              do while(t.ge.tg)
                call follow_plist(1,ifol,ifoln,plist,nbod0,
     &               tg,tstop)
             enddo

             istep = 0
             do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2

                if(abs(id).eq.ifoln) then
                   istep = 1
                   inc = inc*dr
                   capom = capom*dr
                   omega = omega*dr
                   capm = capm*dr
                   peri = a*(1.0d0-e)
                   apo = a*(1.0d0+e)
                   write(7,1000) t,ifoln,a,e,inc,capom,omega,capm,
     &                  peri,apo,mass(abs(id))/mass(1)
 1000              format(1x,e15.7,1x,i3,1x,f10.4,1x,f7.5,4(1x,f9.4),
     &                  2(1x,f10.4),1e13.5)
                   tmax = t
                endif
             enddo

             if(istep.eq.0) goto 2     ! did not find particle this times step

        goto 1

 2      continue

        write(*,*) ' Tmax = ',tmax

        stop
        end
c-------------------------------------------
        subroutine follow_plist(iflg,ifol,ifoln,plist,nbod,tg,tstop)

	include 'swift.inc'
        real*8 tg,tstop
	integer nbod
        integer iflg,ifol,ifoln,plist(nbod),iwhy
        integer ig,im,idum,i,ierr
        save iwhy

        if(iflg.eq.0) then
           open(2,file='discard_mass.out',status='old',iostat=ierr)
           if(ierr.ne.0) then
              write(*,*) 'Could not open discard_mass.out'
              tg = 5.0*tstop
              return            ! <====== NOTE 
           endif
           read(2,*,iostat=ierr) tg,iwhy
           if(ierr.ne.0) then
              write(*,*) 'Could not read discard_mass.out'
              tg = 5.0*tstop
              return            ! <====== NOTE 
           endif
           ifoln = ifol
           return               ! <====== NOTE 
        endif

        if(iwhy.eq.2) then
           read(2,*) idum,im
           read(2,fmt='(1x)')
           read(2,fmt='(1x)')
           read(2,*) idum,ig
           call left_reorder(ig,im,nbod,plist)
           do i=1,5
              read(2,fmt='(1x)')
           enddo
        else
           read(2,*) idum,ig
           im = -1
           call left_reorder(ig,im,nbod,plist)
           read(2,fmt='(1x)')
           read(2,fmt='(1x)')
        endif

        read(2,*,iostat=ierr) tg,iwhy
        if(ierr.ne.0) then
           tg = 5.0 * tstop
        endif

        ifoln = plist(ifol)

        return
        end

c---------------------------------------------------------------------
       subroutine left_reorder(ig,im,nbod,plist)

       include 'swift.inc'

       integer ig,nbod,plist(nbod),im,i

       do i=1,nbod
          if(plist(i).eq.ig) then
             if(im.gt.0) then
                plist(i) = im
             else
                plist(i) = -1
             endif
          endif
       enddo

       do i=1,nbod
          if(plist(i).gt.ig) then
             plist(i) = plist(i) - 1
          endif
       enddo

       return
       end
