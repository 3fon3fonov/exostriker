c*************************************************************************
c                            RMVS3_INTERP.F
c*************************************************************************
c This subroutine interpolates between two kepler orbits.
c For outer region only
c
c             Input:
c                 nbod                ==>  number of massive bodies 
c                                          (int scalar)
c                 xbeg,ybeg,zbeg      ==>  initial planet position in helio 
c                                            (real arrays)
c                 vxbeg,vybeg,vzbeg   ==>  initial planet vel in helio 
c                                            (real arrays)
c                 xend,yend,zend      ==>  final planet position in helio 
c                                            (real arrays)
c                 vxend,vyend,vzend   ==>  final planet position in helio 
c                                            (real arrays)
c                 dt                   ==>  time step (real sclar)
c                 msun                 ==>  mass of sun (real sclar)
c                 nt                   ==>  the number of intermediate steps
c                                           (integer scalar)
c             Output:
c                 xtmp,ytmp,ztmp      ==>  position of planet wrt time 
c                                          for outer region
c                                            (2d real arrays)
c                 vxtmp,vytmp,vztmp   ==>  velocoty of planet wrt time 
c                                          for outer region
c                                            (2d real arrays)
c
c
c Remarks: Based on rmvs2_interp_o 
c Authors:  Hal Levison 
c Date:    7/10/96
c Last revision: 

      subroutine rmvs3_interp(nbod,xbeg,ybeg,zbeg,vxbeg,vybeg,
     &     vzbeg,xend,yend,zend,vxend,vyend,vzend,dt,msun,nt,
     &     xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,nt
      real*8 dt,msun
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)

c...  Outputs:
      real*8 xtmp(NPLMAX,nt),ytmp(NPLMAX,nt)
      real*8 ztmp(NPLMAX,nt)
      real*8 vxtmp(NPLMAX,nt),vytmp(NPLMAX,nt)
      real*8 vztmp(NPLMAX,nt)

c...  Internals
      integer i,iflg,ib
      real*8 xc2(NPLMAX),yc2(NPLMAX),zc2(NPLMAX)
      real*8 vxc2(NPLMAX),vyc2(NPLMAX),vzc2(NPLMAX)
      real*8 xc1(NPLMAX),yc1(NPLMAX),zc1(NPLMAX)
      real*8 vxc1(NPLMAX),vyc1(NPLMAX),vzc1(NPLMAX)
      real*8 dti,dtb,frac,onemf

c----
c...  Executable code 

      dti = dt/float(nt)
      dtb = -1.0d0*dt

c...  move the end positions to beginning
      do i=2,nbod

         xc1(i) = xbeg(i)
         yc1(i) = ybeg(i)
         zc1(i) = zbeg(i)
         vxc1(i) = vxbeg(i)
         vyc1(i) = vybeg(i)
         vzc1(i) = vzbeg(i)

         xc2(i) = xend(i)
         yc2(i) = yend(i)
         zc2(i) = zend(i)
         vxc2(i) = vxend(i)
         vyc2(i) = vyend(i)
         vzc2(i) = vzend(i)
         call drift_one(msun,xc2(i),yc2(i),zc2(i),
     &             vxc2(i),vyc2(i),vzc2(i),dtb,iflg)
         if(iflg.ne.0) then
            write(*,*) ' Planet ',i,' is lost in rmvs2_interp !!!!!!!!!'
            write(*,*) msun,dtb
            write(*,*) xc2(i),yc2(i),zc2(i)
            write(*,*) vxc2(i),vyc2(i),vzc2(i)
            write(*,*) ' STOPPING '
            call util_exit(1)
         endif

      enddo

      do i=2,nbod

         do ib = 1,nt
            
            call drift_one(msun,xc1(i),yc1(i),zc1(i),
     &           vxc1(i),vyc1(i),vzc1(i),dti,iflg)
            if(iflg.ne.0) then
               write(*,*) ' Planet ',i,' is lost in rmvs2_interp_o !!!'
               write(*,*) msun,dtb
               write(*,*) xc1(i),yc1(i),zc1(i)
               write(*,*) vxc1(i),vyc1(i),vzc1(i)
               write(*,*) ' STOPPING '
               call util_exit(1)
            endif

            call drift_one(msun,xc2(i),yc2(i),zc2(i),
     &           vxc2(i),vyc2(i),vzc2(i),dti,iflg)
            if(iflg.ne.0) then
               write(*,*) ' Planet ',i,' is lost in rmvs2_interp_o !!!'
               write(*,*) msun,dtb
               write(*,*) xc2(i),yc2(i),zc2(i)
               write(*,*) vxc2(i),vyc2(i),vzc2(i)
               write(*,*) ' STOPPING '
               call util_exit(1)
            endif

            frac = float(ib)/float(nt)
            onemf = 1.0d0 - frac

            xtmp(i,ib) = onemf*xc1(i) + frac*xc2(i)
            ytmp(i,ib) = onemf*yc1(i) + frac*yc2(i)
            ztmp(i,ib) = onemf*zc1(i) + frac*zc2(i)
            vxtmp(i,ib) = onemf*vxc1(i) + frac*vxc2(i)
            vytmp(i,ib) = onemf*vyc1(i) + frac*vyc2(i)
            vztmp(i,ib) = onemf*vzc1(i) + frac*vzc2(i)

         enddo
      enddo

c...  put zeros in position 1
      do ib = 1,nt
         xtmp(1,ib) = 0.0d0
         ytmp(1,ib) = 0.0d0
         ztmp(1,ib) = 0.0d0
         vxtmp(1,ib) = 0.0d0
         vytmp(1,ib) = 0.0d0
         vztmp(1,ib) = 0.0d0
      enddo

      return
      end      ! rmvs3_interp.f
c-----------------------------------------------------------------------

