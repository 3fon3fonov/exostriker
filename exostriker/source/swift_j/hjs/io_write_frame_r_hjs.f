c*************************************************************************
c                            IO_WRITE_FRAME_R_HJS
c*************************************************************************
c write out a whole frame to an real*4 binary file.
c both massive and test particles (Generalized Jacobi coordinates)
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 oloc          ==> Link between bodies & orbits (2D int array)
c                 mass          ==>  mass of bodies (real array)
c                 eta,mu        ==> Masses of centers & satellites 
c                                    (real arrays)
c                 umat          ==> Conversion matrix Jacobi => Barycentric
c                                    (2D real array)
c                 xj,yj,zj      ==>  current position in Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  current velocity in Jacobi coord 
c                                    (real arrays)
c                 matp          ==> Conversion vectors for tp's
c                                    (2D real array)
c                 etatp          ==> Masses of centers for tp's (real array)
c                 xjt,yjt,zjt    ==>  current part position in Jacobi coord 
c                                      (real arrays)
c                 vxjt,vyjt,vzjt ==>  current velocity in Jacobi coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Adapted from io_write_frame_r.f
c Authors:  Herve Beust 
c Date:    Feb 12, 2002

      subroutine io_write_frame_r_hjs(time,nbod,ntp,oloc,matp,umat,
     &           mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,etatp,xjt,yjt,zjt,
     &           vxjt,vyjt,vzjt,istat,oname,iu,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,ntp,iu
      real*8 mass(nbod),eta(nbod),mu(nbod),time
      integer istat(ntp),oloc(NPLMAX,NPLMAX)
      real*8 matp(NPLMAX,NTPMAX),umat(NPLMAX,NPLMAX),etatp(ntp)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)
      character*(*) oname,fopenstat

c...  Internals
      integer i,id,nc,j,k,jj,kk
      integer ialpha,ierr
      logical ok
      real*8 xb(nbod),yb(nbod),zb(nbod)
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)
      real*8 xt,yt,zt,vxt,vyt,vzt,c1,c2,c3,c,matr(3,3)
      real*8 a,e,inc,capom,omega,capm
      real*8 gm,msys,bot
      real*8 masstemp(NPLMAX),xjtemp(NPLMAX),yjtemp(NPLMAX)
      real*8 zjtemp(NPLMAX)
      real*8 vxjtemp(NPLMAX),vyjtemp(NPLMAX),vzjtemp(NPLMAX)
      real*8 xbtemp(NPLMAX),ybtemp(NPLMAX),zbtemp(NPLMAX)
      real*8 vxbtemp(NPLMAX),vybtemp(NPLMAX),vzbtemp(NPLMAX)
      real*8 umpart(NPLMAX,NPLMAX)
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
      do i=2,nbod
         gm = eta(i) + mu(i)
         id = -1*i
 	 call orbel_xv2el(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),gm,
     &          ialpha,a,e,inc,capom,omega,capm)
         call io_write_line_r(iu,id,a,e,inc,capom,omega,capm)
      end do

c...  write out test particles.
      do i=1,ntp
        if (istat(i).eq.0) then
c... Here we need to calculate the plane perp. to the local angular momentum 
c...    of the centers of the tp. The orbital elements of the tp's are given 
c...    relative to this base plane.
          nc = 0
          do j = 1,nbod
            if (matp(j,i).ne.0.0d0) then
              nc = nc+1
              masstemp(nc) = mass(j)
              xjtemp(nc) = xj(j)
              yjtemp(nc) = yj(j)
              zjtemp(nc) = zj(j)
              vxjtemp(nc) = vxj(j)
              vyjtemp(nc) = vyj(j)
              vzjtemp(nc) = vzj(j)
            end if
          end do
    
c... Of course we calculate the midplane if there is more than one center...
          if (nc.gt.1) then 

c... Now we extract a submatrix from the umat matrix 
            do j = 1,nbod
              do k = 1,nbod
                umpart(k,j) = 0.0d0
              end do
            end do  
            do j = 1,nbod
              umpart(j,1) = 1.0d0
            end do
            jj = 1
            do j = 2,nbod
              ok = .true.
              do k = 1,nbod
                if ((matp(k,i).eq.0.0d0).and.(oloc(j,k).ne.0))
     &                                                ok=.false.
              end do
              if (ok) then
                jj = jj+1
                kk = 0
                do k = 1,nbod
                  if (matp(k,i).ne.0.0d0) then
                    kk = kk+1
                    umpart(kk,jj) = umat(k,j)
                  end if
                end do
              end if
            end do
                   
            call coord_g2b(nc,umpart,masstemp,xjtemp,yjtemp,zjtemp,
     &          vxjtemp,vyjtemp,vzjtemp,xbtemp,ybtemp,zbtemp,vxbtemp,
     &                                                 vybtemp,vzbtemp)
            c1 = 0.0d0
            c2 = 0.0d0
            c3 = 0.0d0

            do j = 1,nc
              c1 = c1 + masstemp(j)*
     &                 (ybtemp(j)*vzbtemp(j)-zbtemp(j)*vybtemp(j))
              c2 = c2 + masstemp(j)*
     &                 (zbtemp(j)*vxbtemp(j)-xbtemp(j)*vzbtemp(j))
              c3 = c3 + masstemp(j)*
     &                 (xbtemp(j)*vybtemp(j)-ybtemp(j)*vxbtemp(j))
            end do
            c = sqrt(c1*c1 + c2*c2 + c3*c3 )
            c1 = c1/c
            c2 = c2/c
            c3 = c3/c
            bot = 1.0d0/(1.0d0 + c3)

            matr(1,1) = 1.0d0 - c1*c1*bot
            matr(1,2) = -1.0d0*c1*c2*bot
            matr(1,3) = -1.0d0*c1
            matr(2,1) = -1.0d0*c1*c2*bot
            matr(2,2) = 1.0d0 - c2*c2*bot
            matr(2,3) = -1.0d0*c2
            matr(3,1) = c1
            matr(3,2) = c2
            matr(3,3) = c3
c... Compute barycentric coordinates restricted to the centers

            xt = matr(1,1)*xjt(i)+matr(1,2)*yjt(i)+matr(1,3)*zjt(i) 
            yt = matr(2,1)*xjt(i)+matr(2,2)*yjt(i)+matr(2,3)*zjt(i) 
            zt = matr(3,1)*xjt(i)+matr(3,2)*yjt(i)+matr(3,3)*zjt(i) 
            vxt = matr(1,1)*vxjt(i)+matr(1,2)*vyjt(i)+matr(1,3)*vzjt(i) 
            vyt = matr(2,1)*vxjt(i)+matr(2,2)*vyjt(i)+matr(2,3)*vzjt(i) 
            vzt = matr(3,1)*vxjt(i)+matr(3,2)*vyjt(i)+matr(3,3)*vzjt(i) 
          else
            xt = xjt(i)
            yt = yjt(i)
            zt = zjt(i)
            vxt = vxjt(i)
            vyt = vyjt(i)
            vzt = vzjt(i)
          end if
c...
          call orbel_xv2el(xt,yt,zt,vxt,vyt,vzt,
     &                  etatp(i),ialpha,a,e,inc,capom,omega,capm)
          call io_write_line_r(iu,i,a,e,inc,capom,omega,capm)
        endif
      enddo

      close(iu)
      return
      end      ! io_write_frame_r_hjs
c----------------------------------------------------------------------
