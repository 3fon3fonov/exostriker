c*************************************************************************
c                            RMVS2_INTERP.F
c*************************************************************************
c This subroutine interpolates between two kepler orbits.
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
c             Output:
c                 xtmp,ytmp,ztmp       ==>  position of planet wrt time 
c                                          for inner region
c                                            (2d real arrays)
c                 xtmpo,ytmpo,ztmpo    ==>  position of planet wrt time 
c                                          for outer region
c                                            (2d real arrays)
c
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    8/25/94
c Last revision: 

      subroutine rmvs2_interp(nbod,xbeg,ybeg,zbeg,vxbeg,vybeg,
     &     vzbeg,xend,yend,zend,vxend,vyend,vzend,dt,msun,
     &     xtmp,ytmp,ztmp,xtmpo,ytmpo,ztmpo)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod
      real*8 dt,msun
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)

c...  Outputs:
      real*8 xtmp(NPLMAX,NTPENC),ytmp(NPLMAX,NTPENC)
      real*8 ztmp(NPLMAX,NTPENC)
      real*8 xtmpo(NPLMAX,NTENC),ytmpo(NPLMAX,NTENC)
      real*8 ztmpo(NPLMAX,NTENC)

c...  Internals
      integer i,iflg,ib,is,ic 
      real*8 xc2(NPLMAX),yc2(NPLMAX),zc2(NPLMAX)
      real*8 vxc2(NPLMAX),vyc2(NPLMAX),vzc2(NPLMAX)
      real*8 xc1(NPLMAX),yc1(NPLMAX),zc1(NPLMAX)
      real*8 vxc1(NPLMAX),vyc1(NPLMAX),vzc1(NPLMAX)
      real*8 dti,dtb,frac,onemf

c----
c...  Executable code 

      dti = dt/float(NTPENC)
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

         ic = 0
         do ib = 1,NTENC
            do is=1,NTPHENC
               ic = ic + 1
            
               call drift_one(msun,xc1(i),yc1(i),zc1(i),
     &              vxc1(i),vyc1(i),vzc1(i),dti,iflg)
               if(iflg.ne.0) then
                  write(*,*) ' Planet ',i,' is lost in rmvs2_interp !!!'
                  write(*,*) msun,dtb
                  write(*,*) xc1(i),yc1(i),zc1(i)
                  write(*,*) vxc1(i),vyc1(i),vzc1(i)
                  write(*,*) ' STOPPING '
                  call util_exit(1)
               endif

               call drift_one(msun,xc2(i),yc2(i),zc2(i),
     &              vxc2(i),vyc2(i),vzc2(i),dti,iflg)
               if(iflg.ne.0) then
                  write(*,*) ' Planet ',i,' is lost in rmvs2_interp !!!'
                  write(*,*) msun,dtb
                  write(*,*) xc2(i),yc2(i),zc2(i)
                  write(*,*) vxc2(i),vyc2(i),vzc2(i)
                  write(*,*) ' STOPPING '
                  call util_exit(1)
               endif

               frac = float(ic)/float(NTPENC)
               onemf = 1.0d0 - frac

               xtmp(i,ic) = onemf*xc1(i) + frac*xc2(i)
               ytmp(i,ic) = onemf*yc1(i) + frac*yc2(i)
               ztmp(i,ic) = onemf*zc1(i) + frac*zc2(i)

            enddo

            frac = float(ib)/float(NTENC)
            onemf = 1.0d0 - frac

            xtmpo(i,ib) = onemf*xc1(i) + frac*xc2(i)
            ytmpo(i,ib) = onemf*yc1(i) + frac*yc2(i)
            ztmpo(i,ib) = onemf*zc1(i) + frac*zc2(i)

         enddo
      enddo

c...  put zeros in position 1
      do ib = 1,NTENC
         xtmpo(1,ib) = 0.0d0
         ytmpo(1,ib) = 0.0d0
         ztmpo(1,ib) = 0.0d0
      enddo
      do is = 1,NTPENC
         xtmp(1,is) = 0.0d0
         ytmp(1,is) = 0.0d0
         ztmp(1,is) = 0.0d0
      enddo

      return
      end      ! rmvs2_interp.f
c-----------------------------------------------------------------------

