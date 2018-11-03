c*************************************************************************
c                            LYAP2_STEP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Uses the difference equations as discussed my Mikola
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp           ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c	                                rstat(i,3), rstat(i,3), rstat(i,5) = 
c                                             the dx(1)-dx(3) from the 
c                                             difference equations.
c	                                rstat(i,6), rstat(i,7), rstat(i,8) = 
c                                             the dv(1)-dv(3) from the 
c                                             difference equations
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adopted from step_kbk
c Authors:  Hal Levison
c Date:    7/11/95
c Last revision: 

      subroutine lyap2_step(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer i1sttp,i,i1stin,iul
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 dx0,dy0,dz0,dvx0,dvy0,dvz0
      real*8 dxht(NTPMAX),dyht(NTPMAX),dzht(NTPMAX)
      real*8 dvxht(NTPMAX),dvyht(NTPMAX),dvzht(NTPMAX)

      real*8 dist(NTPMAX)
      real*8 tin,tlout,dtlout

      data i1stin/0/

      save i1stin,tin,tlout,dtlout,iul
      save dxht,dyht,dzht,dvxht,dvyht,dvzht

c----
c...  Executable code 

c...  set things up if this is the initial call
      if(i1stin.eq.0) then

        if( (j2rp2.ne.0.0d0) .or. (j4rp4.ne.0.0d0) ) then
           write(*,*) 'LYAP2 routines must have J2,J4=0!'
           call util_exit(1)   !    <==== NOTE!
        endif                 

	write(*,*) 'Input how aften the  distance is written:'
	read(*,*) dtlout
        write(*,*) 'Input the initial separation in phase space'
        write(*,*) '  If <0 then use the values currently in rstat'
        read(*,*) dx0,dy0,dz0,dvx0,dvy0,dvz0
	iul = 50
        tin = 0.0
        i1stin = 1
        tlout = dtlout
        if(dx0.gt.0.0d0) then
           do i=1,ntp
              dxht(i) = dx0
              dyht(i) = dy0
              dzht(i) = dz0
              dvxht(i) = dvx0
              dvyht(i) = dvy0
              dvzht(i) = dvz0
           enddo
        else
           do i=1,ntp
              dxht(i) = rstat(i,3)
              dyht(i) = rstat(i,4)
              dzht(i) = rstat(i,5)
              dvxht(i) = rstat(i,6)
              dvyht(i) = rstat(i,7)
              dvzht(i) = rstat(i,8)
           enddo
        endif

        do i=1,ntp
           if(istat(i,1).ne.0) then
              dist(i) = 0.0d0
           else
              dist(i) = sqrt( dxht(i)**2 + dyht(i)**2 + 
     &             dzht(i)**2 + dvxht(i)**2 + dvyht(i)**2 + 
     &             dvzht(i)**2 )
           endif
        enddo
        call io_lyap2_write(iul,time,dist,ntp)

        write(*,*) ' CONTINUE: '
      endif

c...  for the normal step

      i1sttp = i1st

c...  remember the current position of the planets
      do i=1,nbod
         xbeg(i) = xh(i)
         ybeg(i) = yh(i)
         zbeg(i) = zh(i)
      enddo

c...  first do the planets
      call step_kdk_pl(i1st,nbod,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,dt)

c...  now remember these positions
      do i=1,nbod
         xend(i) = xh(i)
         yend(i) = yh(i)
         zend(i) = zh(i)
      enddo

c...  next the test particles
      call lyap2_step_tp(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,
     &              dxht,dyht,dzht,dvxht,dvyht,dvzht,istat,dt)	

      tin = tin + dt

c...  now do the lyap anal stuff
      if(tin .ge. tlout) then 

         do i=1,ntp
            if(istat(i,1).ne.0) then
               dist(i) = 0.0d0
            else
               dist(i) = sqrt( dxht(i)**2 + dyht(i)**2 + 
     &              dzht(i)**2 + dvxht(i)**2 + dvyht(i)**2 + 
     &              dvzht(i)**2 )
            endif
         enddo

         call io_lyap2_write(iul,time,dist,ntp)

         tlout = tlout + dtlout
      endif

      do i=1,ntp
         rstat(i,3) = dxht(i)
         rstat(i,4) = dyht(i)
         rstat(i,5) = dzht(i)
         rstat(i,6) = dvxht(i)
         rstat(i,7) = dvyht(i)
         rstat(i,8) = dvzht(i)
      enddo

      return
      end   ! lyap2_step
c------------------------------------------------------------------------

