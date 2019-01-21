c**********************************************************************
c		      SWIFT_MVS_J.F
c**********************************************************************
c
c                 NO CLOSE ENCOUNTERS
c                 To run, need 2 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c
c  This version inputs/outputs Jacobi coords and orbital elements and
c  groups terms so that hierarchical systems with comparable masses can
c  be integrated.
c  NOTE:  No test particles in this code
c
c Author:  Man Hoi Lee
c Date:    12/6/01
c Last revision: 

     
	include 'swift.inc'

	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)

	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	real*8 xjt(1),yjt(1),zjt(1)       ! Dummy for the io
	real*8 vxjt(1),vyjt(1),vzjt(1)

	integer istat(1,NSTAT),i1st
	integer nbod,ntp,nleft
	integer iflgchk,iub,iuj,iud,iue
        real*8 rstat(1,NSTATR)

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tout,tdump,tfrac,eoff

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose 

	character*80 outfile,inparfile,inplfile,fopenstat


c-----
c...    Executable code 

	ntp = 0

c...    print version number
        call util_version

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
	call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xj,yj,zj,
     &       vxj,vyj,vzj,rplsq,j2rp2,j4rp4)

c Initialize initial time and times for first output and first dump
	t = t0
	tout = t0 + dtout
	tdump = t0 + dtdump

        iub = 20
        iuj = 30
        iud = 40
        iue = 60

c...    Do the initial io write
        if(btest(iflgchk,1))  then ! bit 1 is set
           call io_write_frame_j(t0,nbod,ntp,mass,xj,yj,zj,vxj,vyj,vzj,
     &         xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,2))  then    ! bit 2 is set
           eoff = 0.0d0
           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,xh,yh,zh,
     &          vxh,vyh,vzh)
           call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &          vyh,vzh,iue,fopenstat,eoff)
        endif

c...    must initize discard io routine
        if(btest(iflgchk,4))  then ! bit 4 is set
           call io_discard_write(0,t,nbod,ntp,xj,yj,zj,vxj,vyj,
     &          vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,iud,
     &          'discard.out',fopenstat,nleft)
        endif

        nleft = ntp
	i1st = 0
c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '

	  do while ( (t .le. tstop) .and. 
     &       ((ntp.eq.0).or.(nleft.gt.0)) )

 	     call step_kdk_pl_j(i1st,nbod,mass,j2rp2,j4rp4,
     &            xj,yj,zj,vxj,vyj,vzj,dt)

	     t = t + dt

             if(btest(iflgchk,4))  then    ! bit 4 is set
                call discard(t,dt,nbod,ntp,mass,xj,yj,zj,vxj,vyj,vzj,
     &               xjt,yjt,zjt,vxjt,vyjt,vzjt,rmin,rmax,rmaxu,
     &               qmin,lclose,rplsq,istat,rstat)
                call io_discard_write(1,t,nbod,ntp,xj,yj,zj,vxj,vyj,
     &               vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,iud,
     &               'discard.out',fopenstat,nleft)
             else
                nleft = ntp
             endif

c if it is time, output orb. elements, 
	  if(t .ge. tout) then 

             if(btest(iflgchk,1))  then    ! bit 1 is set
                call  io_write_frame_j(t,nbod,ntp,mass,xj,yj,zj,vxj,
     &               vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,outfile,
     &               iub,fopenstat)
              endif

	    tout = tout + dtout
	  endif

c If it is time, do a dump
          if(t.ge.tdump) then

             tfrac = (t-t0)/(tstop-t0)
             write(*,998) t,tfrac,nleft
 998         format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of active tp =',i4)
	     call io_dump_pl('dump_pl.dat',nbod,mass,xj,yj,zj,
     &                 vxj,vyj,vzj,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	     call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &           dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
	     tdump = tdump + dtdump

             if(btest(iflgchk,2))  then    ! bit 2 is set
                call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,xh,yh,zh,
     &               vxh,vyh,vzh)
                call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &               xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
             endif

	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

	call io_dump_pl('dump_pl.dat',nbod,mass,xj,yj,zj,
     &            vxj,vyj,vzj,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

        call util_exit(0)
        end    ! swift_mvs.f
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

      include '../swift.inc'
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

c*************************************************************************
c                            STEP_KDK_PL_J.F
c*************************************************************************
c This subroutine takes a step in Jacobi coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xj,yj,zj      ==>  initial position in Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial velocity in Jacobi coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xj,yj,zj      ==>  final position in Jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in Jacobi coord 
c                                       (real arrays)
c
c Remarks: Adopted from step_kdk_pl
c Authors:  Man Hoi Lee
c Date:    12/6/01
c Last revision: 

      subroutine step_kdk_pl_j(i1st,nbod,mass,j2rp2,j4rp4,
     &     xj,yj,zj,vxj,vyj,vzj,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Internals:
      real*8 dth 
      real*8 axj(NPLMAX),ayj(NPLMAX),azj(NPLMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

      save axj,ayj,azj,xh,yh,zh     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then
c...      Convert to helio coords
          call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &          xh,yh,zh,vxh,vyh,vzh)
c...     Get the accelerations in Jacobi frame. if frist time step
         call getaccj(nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &         xj,yj,zj,axj,ayj,azj)
         i1st = 1    ! turn this off
      endif

c...  Apply a Jacobi kick for a half dt 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dth)

c..   Drift in Jacobi coords for the full step 
      call drift2(nbod,mass,xj,yj,zj,vxj,vyj,vzj,dt)

c...  After drift, compute helio. xh and vh for acceleration calculations
      call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &             xh,yh,zh,vxh,vyh,vzh)

c...  Get the accelerations in Jacobi frame.
      call getaccj(nbod,mass,j2rp2,j4rp4,xh,yh,zh,xj,yj,zj,axj,ayj,azj)

c...  Apply a Jacobi kick for a half dt 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dth)

      return
      end   ! step_kdk_pl_j
c---------------------------------------------------------------------

c*************************************************************************
c                        DRIFT2.F
c*************************************************************************
c This subroutine loops thorugh the particles and calls the danby routine
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xj,yj,zj      ==>  initial position in jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial position in jacobi coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xj,yj,zj      ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final position in jacobi coord 
c                                       (real arrays)
c
c Remarks: Based on drift
c Authors:  Man Hoi Lee
c Date:    12/6/01
c Last revision:

      subroutine drift2(nbod,mass,xj,yj,zj,vxj,vyj,vzj,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Internals:
	real*8 etaj
	integer j,iflg

c----
c...  Executable code 

c Take a drift forward dt

	etaj = mass(1)
	do j = 2,nbod
	   etaj = etaj + mass(j)
	   call drift_one(etaj,xj(j),yj(j),zj(j),
     &             vxj(j),vyj(j),vzj(j),dt,iflg)
           if(iflg.ne.0) then
              write(*,*) ' Planet ',j,' is lost !!!!!!!!!'
              write(*,*) etaj,dt
              write(*,*) xj(j),yj(j),zj(j)
              write(*,*) vxj(j),vyj(j),vzj(j)
              write(*,*) ' STOPPING '
              call util_exit(1)
           endif
	enddo

	return
	end
c--------------------------------------------------------------------------

c*************************************************************************
c                        GETACCJ.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the JACOBI frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c                 xj,yj,zj    ==>  position in jacobi coord (real arrays)
c             Output:
c                 axj,ayj,azj ==>  acceleration in Jacobi coord (real arrays)
c
c Remarks: Based on getacch
c Author:  Man Hoi Lee
c Date:    12/6/01
c Last revision:

      subroutine getaccj(nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &     xj,yj,zj,axj,ayj,azj)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),j2rp2,j4rp4
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)

c...  Outputs:
      real*8 axj(NPLMAX),ayj(NPLMAX),azj(NPLMAX)
                
c...  Internals:
      integer i
      real*8 ir3h(NPLMAX),ir3j(NPLMAX)
      real*8 irh(NPLMAX),irj(NPLMAX)
      real*8 axj1(NPLMAX),ayj1(NPLMAX),azj1(NPLMAX)
      real*8 axj2(NPLMAX),ayj2(NPLMAX),azj2(NPLMAX)
      real*8 axj3(NPLMAX),ayj3(NPLMAX),azj3(NPLMAX)
      real*8 aoblx(NPLMAX),aobly(NPLMAX),aoblz(NPLMAX)

c----
c...  Executable code 

c...  get thr r^-3's
      call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
      call getacch_ir3(nbod,2,xj,yj,zj,ir3j,irj)

c...  now the first terms
      call getaccj_aj1(nbod,mass,xh,yh,zh,xj,yj,zj,ir3h,ir3j,
     &                 axj1,ayj1,azj1)

c...  now the second terms
      call getaccj_aj2(nbod,mass,xh,yh,zh,ir3h,axj2,ayj2,azj2)

c...  now the third terms
      call getaccj_aj3(nbod,mass,xh,yh,zh,axj3,ayj3,azj3)

c...  add them all together
      axj(1) = 0.0
      ayj(1) = 0.0
      azj(1) = 0.0
      do i=2,nbod
        axj(i) = axj1(i) + axj2(i) + axj3(i)
        ayj(i) = ayj1(i) + ayj2(i) + ayj3(i)
        azj(i) = azj1(i) + azj2(i) + azj3(i)
      enddo

c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call obl_acc(nbod,mass,j2rp2,j4rp4,xj,yj,zj,irh,
     &        aoblx,aobly,aoblz)
         do i = 2,nbod
            axj(i) = axj(i) + aoblx(i) - aoblx(1)
            ayj(i) = ayj(i) + aobly(i) - aobly(1)
            azj(i) = azj(i) + aoblz(i) - aoblz(1)
         enddo
      endif

      return
      end      ! getaccj

c---------------------------------------------------------------------

c*************************************************************************
c                        GETACCJ_AJ1.F
c*************************************************************************
c This subroutine calculates the 1st term of acceleration 
c on the massive particles in the JACOBI frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalar)
c                 mass        ==>  mass of bodies (real array)
c                 xh,yh,zh    ==>  position in heliocentric coord (real array)
c                 xj,yj,zj    ==>  position in jacobi coord (real array)
c                 ir3h        ==> inv radii in heliocentric coord (real array)
c                 ir3j        ==> inv radii in jacobi coord (real array)
c             Output:
c                 axj1,ayj1,azj1 ==>  1st term acceleration in Jacobi coord 
c                                    (real array)
c
c Author:  Man Hoi Lee
c Date:    12/6/01
c Last revision:

      subroutine getaccj_aj1(nbod,mass,xh,yh,zh,xj,yj,zj,ir3h,ir3j,
     &                 axj1,ayj1,azj1)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),ir3h(nbod),ir3j(nbod)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axj1(nbod),ayj1(nbod),azj1(nbod)

                
c...  Internals:
      integer i
      real*8 etaim1,etai,aj1h,aj1j

c----
c...  Executable code 

      axj1(1) = 0.0
      ayj1(1) = 0.0
      azj1(1) = 0.0

      axj1(2) = 0.0     ! because xj=xh
      ayj1(2) = 0.0
      azj1(2) = 0.0

      etaim1 = mass(1) + mass(2)

      do i=3,nbod

         etai = etaim1 + mass(i)

         aj1j = xj(i)*ir3j(i)
         aj1h = mass(1)*xh(i)*ir3h(i)/etaim1
         axj1(i) = etai*(aj1j - aj1h)

         aj1j = yj(i)*ir3j(i)
         aj1h = mass(1)*yh(i)*ir3h(i)/etaim1
         ayj1(i) = etai*(aj1j - aj1h)

         aj1j = zj(i)*ir3j(i)
         aj1h = mass(1)*zh(i)*ir3h(i)/etaim1
         azj1(i) = etai*(aj1j - aj1h)

         etaim1 = etai

      enddo

      return
      end   ! getaccj_aj1

c---------------------------------------------------------------------

c*************************************************************************
c                        GETACCJ_AJ2.F
c*************************************************************************
c This subroutine calculates the 2nd term of acceleration 
c on the massive particles in the JACOBI frame.
c             Input:
c                 nbod        ==>  number of massive bodies (int scalar)
c                 mass        ==>  mass of bodies (real array)
c                 xh,yh,zh    ==>  position in helio coord (real array)
c                 ir3h        ==> inv radii in helio coord (real array)
c             Output:
c                 axj2,ayj2,azj2 ==>  2nd term acceleration in Jacobi coord
c                                    (real array)
c
c Author:  Man Hoi Lee
c Date:    12/6/01
c Last revision:

      subroutine getaccj_aj2(nbod,mass,xh,yh,zh,ir3h,axj2,ayj2,azj2)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),ir3h(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axj2(nbod),ayj2(nbod),azj2(nbod)
                
c...  Internals:
      integer i
      real*8 etaim1, fac

c----
c...  Executable code 

      axj2(nbod) = 0.0
      ayj2(nbod) = 0.0
      azj2(nbod) = 0.0

      do i=nbod,3,-1
         fac = mass(i)*ir3h(i)
         axj2(i-1) = axj2(i) - fac*xh(i)
         ayj2(i-1) = ayj2(i) - fac*yh(i)
         azj2(i-1) = azj2(i) - fac*zh(i)
      enddo

      axj2(1) = 0.0
      ayj2(1) = 0.0
      azj2(1) = 0.0

      etaim1 = 0.0d0
      do i=2,nbod-1
         etaim1 = etaim1 + mass(i-1)
         axj2(i) = axj2(i)*mass(1)/etaim1
         ayj2(i) = ayj2(i)*mass(1)/etaim1
         azj2(i) = azj2(i)*mass(1)/etaim1
      enddo

      return
      end     ! getaccj_aj2
c---------------------------------------------------------------------

c*************************************************************************
c                        GETACCJ_AJ3.F
c*************************************************************************
c This subroutine calculates the 3rd term acceleration on the massive particles
c in the JACOBI frame. This term is the direct cross terms
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  position in heliocentric coord 
c                                   (real arrays)
c             Output:
c                 axj3,ayj3,azj3 ==>  3rd term of acceleration in Jacobi coord
c                                     (real arrays)
c
c Author:  Man Hoi Lee
c Date:    12/6/01
c Last revision:

      subroutine getaccj_aj3(nbod,mass,xh,yh,zh,axj3,ayj3,azj3)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axj3(nbod),ayj3(nbod),azj3(nbod)

c...  Internals:
      integer i,j
      real*8 axh3(NPLMAX),ayh3(NPLMAX),azh3(NPLMAX)
      real*8 dx,dy,dz,rji2,faci,facj,irij3
      real*8 etaim1,sumax,sumay,sumaz

c------
c...  Executable code

      do i=1,nbod
         axh3(i) = 0.0
         ayh3(i) = 0.0
         azh3(i) = 0.0
      enddo

      do i=2,nbod-1
         do j=i+1,nbod

             dx = xh(j) - xh(i)
             dy = yh(j) - yh(i)
             dz = zh(j) - zh(i)
             rji2 = dx*dx + dy*dy + dz*dz

             irij3 = 1.0d0/(rji2*sqrt(rji2))
             faci = mass(i)*irij3
             facj = mass(j)*irij3

             axh3(j) = axh3(j) - faci*dx
             ayh3(j) = ayh3(j) - faci*dy
             azh3(j) = azh3(j) - faci*dz

             axh3(i) = axh3(i) + facj*dx
             ayh3(i) = ayh3(i) + facj*dy
             azh3(i) = azh3(i) + facj*dz

         enddo
      enddo

      axj3(1) = 0.d0
      ayj3(1) = 0.d0
      azj3(1) = 0.d0

      axj3(2) = axh3(2)
      ayj3(2) = ayh3(2)
      azj3(2) = azh3(2)

      etaim1 = mass(1)

      sumax = 0.d0
      sumay = 0.d0
      sumaz = 0.d0

      do i=3,nbod
         etaim1 = etaim1 + mass(i-1)

         sumax = sumax + mass(i-1)*axh3(i-1)
         sumay = sumay + mass(i-1)*ayh3(i-1)
         sumaz = sumaz + mass(i-1)*azh3(i-1)

         axj3(i) = axh3(i) - sumax/etaim1
         ayj3(i) = ayh3(i) - sumay/etaim1
         azj3(i) = azh3(i) - sumaz/etaim1

      enddo

      return
      end     ! getaccj_aj3

c----------------------------------------------------------------------
