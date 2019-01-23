c*************************************************************************
c                            LYAP_STEP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c both massive and test particles and shadow particles for test guys.
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
c Authors:  Martin Duncan
c Date:    5/5/93
c Last revision: 2/24/94

      subroutine lyap_step(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
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
      integer i1sttp,i1stsh,i,i1stin,iul
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 xsh(NTPMAX),ysh(NTPMAX),zsh(NTPMAX)
      real*8 vxsh(NTPMAX),vysh(NTPMAX),vzsh(NTPMAX)
      real*8 dist0(NTPMAX),lrsum(NTPMAX),logpr(NTPMAX)
      real*8 tin,tnorm,dtnorm
      character*80 inshfile

      data i1stin/0/

      save i1stin,tin,tnorm,dtnorm
      save xsh,ysh,zsh
      save vxsh,vysh,vzsh
      save dist0,lrsum,iul

c----
c...  Executable code 

c...  set things up if this is the initial call
      if(i1stin.eq.0) then
	write(*,*) 'Enter name of shadow data file : '
	read(*,999) inshfile
999 	format(a)
	iul = 50
	call io_lyap_init(inshfile,ntp,xht,yht,zht,vxht,vyht,vzht,
     &            xsh,ysh,zsh,vxsh,vysh,vzsh,dist0,lrsum,dtnorm,iul)
        tin = 0.0
        i1stin = 1
        tnorm = dtnorm
        write(*,*) ' CONTINUE: '
      endif

c...  for the normal step

      i1sttp = i1st
      i1stsh = i1st

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
      call step_kdk_tp(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istat,dt)	

c... finally the shadow particles 
      call lyap_step_sh(i1stsh,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xsh,ysh,zsh,vxsh,vysh,vzsh,istat,dt)	

      tin = tin + dt

c...  now do the lyap anal stuff
      if(tin .ge. tnorm) then 

         do i=1,ntp
            if(istat(i,1).ne.0) then
               logpr(i) = 0.d0
            else
               call lyap_renorm(xht(i),yht(i),zht(i),
     &              vxht(i),vyht(i),vzht(i),xsh(i),ysh(i),zsh(i),
     &              vxsh(i),vysh(i),vzsh(i),dist0(i),lrsum(i))
               logpr(i) = dlog10(lrsum(i)/tin) 
            endif
         enddo

         call io_lyap_write(iul,tin,logpr,ntp)

         tnorm = tnorm + dtnorm
      endif

      return
      end   ! lyap_step
c------------------------------------------------------------------------

