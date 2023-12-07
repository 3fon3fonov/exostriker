c**********************************************************************
c		      SWIFT_SYMBA5_J.F
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
c  This version inputs/outputs Jacobi coords and orbital elements.
c
c Authors: Man Hoi Lee
c Date:    1/20/05
c Last revision: 

     
      include 'swift.inc'

      real*8 mass(NTPMAX),j2rp2,j4rp4
      real*8 xj(NTPMAX),yj(NTPMAX),zj(NTPMAX)
      real*8 vxj(NTPMAX),vyj(NTPMAX),vzj(NTPMAX)

      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

      real*8 xjt(1),yjt(1),zjt(1)       ! Dummy for the io
      real*8 vxjt(1),vyjt(1),vzjt(1)
      integer ntp,istat(1,NSTAT)

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

      character*80 outfile,inparfile,inplfile,fopenstat


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
     &     xj,yj,zj,vxj,vyj,vzj,rpl,rhill,j2rp2,j4rp4)

      write(*,*) 'Enter the smallest mass to self gravitate :'
      read(*,*) mtiny
      write(*,*) ' mtiny = ',mtiny

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
      if(btest(iflgchk,1))  then ! bit 1 is set
         call io_write_frame_j(t0,nbod,ntp,mass,xj,yj,zj,vxj,vyj,vzj,
     &        xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,outfile,iub,fopenstat)
         call io_write_mass(t0,nbod,mass,outfile,ium,fopenstat)
      endif
      if(btest(iflgchk,2))  then ! bit 2 is set
         eoff = 0.0d0
         call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,xh,yh,zh,
     &        vxh,vyh,vzh)
         call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &        vyh,vzh,iue,fopenstat,eoff)
      endif

c...  must initize discard io routine
      if(btest(iflgchk,4))  then ! bit 4 is set
         call io_discard_mass(0,t,0,mass(1),rpl(1),xj(1),yj(1),zj(1),
     &        vxj(1),vyj(1),vzj(1),iud,-1,fopenstat)
      endif

c...  Calculate the location of the last massive particle
      call symba5_nbodm(nbod,mass,mtiny,nbodm)

      ihills = 0
      i1st = 0
c***************here's the big loop *************************************
      write(*,*) ' ************** MAIN LOOP ****************** '

      do while ( (t .le. tstop) .and. (nbod.gt.1) )

         i1st = 0

         call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,xh,yh,zh,
     &        vxh,vyh,vzh)

         call symba5_step_pl(i1st,t,nbod,nbodm,mass,j2rp2,j4rp4,
     &        xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,
     &        mergelst,mergecnt,iecnt,eoff,rhill,mtiny)

         call coord_h2j(nbod,mass,xh,yh,zh,vxh,vyh,vzh,xj,yj,zj,
     &        vxj,vyj,vzj)

         t = t + dt

         if(btest(iflgchk,4))  then ! bit 4 is set
            nbodo = nbod
            call discard_massive5(t,dt,nbod,mass,xj,yj,zj,
     &           vxj,vyj,vzj,rmin,rmax,rmaxu,qmin,lclose,
     &           rpl,rhill,isenc,mergelst,mergecnt,
     &           iecnt,eoff,i1st)
            if(nbodo.ne.nbod) then
               call symba5_nbodm(nbod,mass,mtiny,nbodm)
            endif
         endif


c if it is time, output orb. elements, 
         if(t .ge. tout) then 

            if(btest(iflgchk,1))  then ! bit 1 is set
               call  io_write_frame_j(t,nbod,ntp,mass,xj,yj,zj,vxj,
     &              vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,outfile,
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
            call io_dump_pl_symba('dump_pl.dat',nbod,mass,xj,yj,zj,
     &           vxj,vyj,vzj,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
            call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &           dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
            tdump = tdump + dtdump

            if(btest(iflgchk,2))  then ! bit 2 is set
               call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,xh,yh,zh,
     &              vxh,vyh,vzh)
               call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &              xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
            endif
            
	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

	call io_dump_pl_symba('dump_pl.dat',nbod,mass,xj,yj,zj,
     &            vxj,vyj,vzj,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
	call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

        call util_exit(0)
        end    ! swift_symba5_j_migrate6.f
c---------------------------------------------------------------------

c*************************************************************************
c                            IO_WRITE_FRAME_J
c*************************************************************************
c write out a whole frame to an real*4 binary file.
c both massive and test particles
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xj,yj,zj      ==>  current position in Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  current velocity in Jacobi coord 
c                                    (real arrays)
c                 xjt,yjt,zjt    ==>  current part position in Jacobi coord 
c                                      (real arrays)
c                 vxjt,vyjt,vzjt ==>  current velocity in Jacobi coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      subroutine io_write_frame_j(time,nbod,ntp,mass,xj,yj,zj,vxj,
     &           vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,oname,
     &           iu,fopenstat)

      include 'swift.inc'
      include '../io/io.inc'

c...  Inputs: 
      integer nbod,ntp,iu
      real*8 mass(nbod),time
      integer istat(1,NSTAT)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)
      character*80 oname,fopenstat

c...  Internals
      integer i,id
      integer ialpha,ierr
      real*8 a,e,inc,capom,omega,capm
      real*8 gm
      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st

c----
c...  Executable code 

c...  if first time through open file
      if(i1st.eq.0) then
         call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_frame: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1st = 1
      else
        call io_open(iu,oname,'append','UNFORMATTED',ierr)
      endif

      call io_write_hdr_r(iu,time,nbod,ntp,istat)
      
c...  write out planets
      gm = mass(1)
      do i=2,nbod
         gm = gm + mass(i)
         id = -1*i
 	 call orbel_xv2el(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),gm,
     &          ialpha,a,e,inc,capom,omega,capm)
         call io_write_line_r(iu,id,a,e,inc,capom,omega,capm)
      enddo

c...  write out test particles
      gm = mass(1)
      do i=1,ntp
         if(istat(i,1).eq.0) then
            call orbel_xv2el(xjt(i),yjt(i),zjt(i),vxjt(i),vyjt(i),
     &          vzjt(i),gm,ialpha,a,e,inc,capom,omega,capm)
            call io_write_line_r(iu,i,a,e,inc,capom,omega,capm)
         endif
      enddo

      close(iu)
      return
      end      ! io_write_frame_j
c----------------------------------------------------------------------
