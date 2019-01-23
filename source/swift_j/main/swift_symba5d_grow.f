c**********************************************************************
c		      SWIFT_SYMBA5D_GROW.F
c**********************************************************************
c
c                 To run, need 2 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c
c  NOTE:  No test particles in this code and the massive bodies 
c         are dimensioned at NTPMAX
c
c  This version allows the first planet to grow from its initial mass to a
c  final mass in time tgrow.
c  This version also allows gas drag on all but the first planet.
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    11/21/96
c Last revision: 12/22/99 MHL

     
      include 'swift.inc'

      real*8 mass(NTPMAX),j2rp2,j4rp4
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

      real*8 xht(1),yht(1),zht(1)       ! Dummy for the io
      real*8 vxht(1),vyht(1),vzht(1)
      integer ntp,istat(1)

      integer nbod,i1st,i,nbodm,nbodo
      integer iflgchk,iub,iuj,iud,iue,ium
      
      real*8 t0,tstop,dt,dtout,dtdump
      real*8 t,tout,tdump,tfrac,eoff
      real*8 rpl(NTPMAX),rhill(NTPMAX)

      real*8 rmin,rmax,rmaxu,qmin,mtiny
      logical*2 lclose 
      integer isenc,ihills
      integer mergelst(2,NTPMAX),mergecnt
      integer*2 iecnt(NTPMAX)

      real*8 mass2i,mass2f,rpl2i,tgrow,ttmp

      real*8 kdrag0,nkdrag,eta0,neta

      character*80 outfile,inparfile,inplfile,fopenstat

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      real*8 ke,pot,energy,eltot(3)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c-----
c...  Executable code 

      ntp = 0

c...  print version number
      call util_version

c Get data for the run and the test particles
      write(*,*) 'Enter name of parameter data file : '
      read(*,999) inparfile
      call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &     iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
      write(*,*) ' '
      write(*,*) 'Enter name of planet data file : '
      read(*,999) inplfile
 999  format(a)
      call io_init_pl_symba(inplfile,lclose,iflgchk,nbod,mass,
     &     xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

      write(*,*) 'Enter the smallest mass to self gravitate :'
      read(*,*) mtiny
      write(*,*) ' mtiny = ',mtiny

c Prompt and read initial mass and radius, final mass, and growth time of
c first planet
      write(*,*) 'Enter initial mass and radius of first planet :'
      read(*,*) mass2i,rpl2i
      write(*,*) 'Enter final mass and growth time of first planet :'
      read(*,*) mass2f,tgrow
      write(*,*) ' mass2i,rpl2i,mass2f,tgrow = ',
     &     mass2i,rpl2i,mass2f,tgrow

c Prompt and read the "drag constant" K_drag,0 at R=1, the exponent nkdrag
c of its radial dependence, the gas velocity parameter eta_0 at r=1,
c and its exponent neta.
      write(*,*) 'Enter kdrag0 and nkdrag :'
      read(*,*) kdrag0,nkdrag
      write(*,*) 'Enter eta0 and neta :'
      read(*,*) eta0,neta
      write(*,*) ' kdrag0,nkdrag,eta0,neta = ',
     &     kdrag0,nkdrag,eta0,neta

c Initialize initial time and times for first output and first dump
      t = t0
      tout = t0 + dtout
      tdump = t0 + dtdump

      iub = 20
      iuj = 30
      iud = 40
      iue = 60
      ium = 21

c...    Do the initial io write
      if(btest(iflgchk,0))  then ! bit 0 is set
         call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &        xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
         call io_write_mass(t0,nbod,mass,outfile,ium,fopenstat)
      endif
      if(btest(iflgchk,1))  then ! bit 1 is set
         call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &        xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
         call io_write_mass(t0,nbod,mass,outfile,ium,fopenstat)
      endif
      if(btest(iflgchk,2))  then ! bit 2 is set
         eoff = 0.0d0
         call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &        vyh,vzh,iue,fopenstat,eoff)
      endif

c...  must initize discard io routine
      if(btest(iflgchk,4))  then ! bit 4 is set
         call io_discard_mass(0,t,0,mass(1),rpl(1),xh(1),yh(1),zh(1),
     &        vxh(1),vyh(1),vzh(1),iud,-1,fopenstat)
      endif

c...  Calculate the location of the last massive particle
      call symba5_nbodm(nbod,mass,mtiny,nbodm)

      ihills = 0
      i1st = 0
c***************here's the big loop *************************************
      write(*,*) ' ************** MAIN LOOP ****************** '

      do while ( (t .le. tstop) .and. (nbod.gt.1) )

         ttmp = t + 0.5d0*dt

c Update first planet's mass by half dt
         if(ttmp.le.tgrow) then
            mass(2) = mass2i + (mass2f - mass2i)*(ttmp/tgrow)**2*
     &                         (3.0d0 - 2.0d0*(ttmp/tgrow))
         else
            mass(2) = mass2f
         endif
         rpl(2) = rpl2i*(mass(2)/mass2i)**(1.0d0/3.0d0)
         call util_hills1(mass(1),mass(2),xh(2),yh(2),zh(2),
     &        vxh(2),vyh(2),vzh(2),rhill(2))

         call symba5d_step_pl(i1st,t,nbod,nbodm,mass,j2rp2,j4rp4,
     &        xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,
     &        mergelst,mergecnt,iecnt,eoff,rhill,mtiny,
     &        kdrag0,nkdrag,eta0,neta)

         t = t + dt

c Update first planet's mass by half dt
         if(t.le.tgrow) then
            mass(2) = mass2i + (mass2f - mass2i)*(t/tgrow)**2*
     &                         (3.0d0 - 2.0d0*(t/tgrow))
         else
            mass(2) = mass2f
         endif
         rpl(2) = rpl2i*(mass(2)/mass2i)**(1.0d0/3.0d0)
         call util_hills1(mass(1),mass(2),xh(2),yh(2),zh(2),
     &        vxh(2),vyh(2),vzh(2),rhill(2))

         if(btest(iflgchk,4))  then ! bit 4 is set
            nbodo = nbod
            call discard_massive5(t,dt,nbod,mass,xh,yh,zh,
     &           vxh,vyh,vzh,rmin,rmax,rmaxu,qmin,lclose,
     &           rpl,rhill,isenc,mergelst,mergecnt,
     &           iecnt,eoff,i1st)
            if(nbodo.ne.nbod) then
               call symba5_nbodm(nbod,mass,mtiny,nbodm)
            endif
         endif


c if it is time, output orb. elements, 
         if(t .ge. tout) then 

            if(btest(iflgchk,0))  then ! bit 0 is set
               call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &              vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &              iub,fopenstat)
               call io_write_mass(t,nbod,mass,outfile,ium,fopenstat)
            endif
            if(btest(iflgchk,1))  then ! bit 1 is set
               call  io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &              vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &              iub,fopenstat)
               call io_write_mass(t,nbod,mass,outfile,ium,fopenstat)
            endif

	    tout = tout + dtout
         endif

c If it is time, do a dump
         if(t.ge.tdump) then

            tfrac = (t-t0)/(tstop-t0)
            write(*,998) t,tfrac,nbod
 998        format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of bodies =',i4)
            call io_dump_pl_symba('dump_pl.dat',nbod,mass,xh,yh,zh,
     &           vxh,vyh,vzh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
            call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &           dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
            tdump = tdump + dtdump

            if(btest(iflgchk,2))  then ! bit 2 is set
               call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &              xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
            endif
            
	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

	call io_dump_pl_symba('dump_pl.dat',nbod,mass,xh,yh,zh,
     &            vxh,vyh,vzh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
	call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

        call util_exit(0)
        end    ! swift_symba5d_grow.f
c---------------------------------------------------------------------
