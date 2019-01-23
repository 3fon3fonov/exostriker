C**************************************************************************
C	    		        BS_INT
C**************************************************************************
c This is the subroutine that does the the bs integration step.
c
c             Input:
c              nbod          ==> number of planets  (int scalar)
c              ntp           ==> number of test particles  (int scalar)
c              mass          ==>  mass of bodies (real array)
c              j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c              istat         ==>  status of the test paricles
c                                       (2d integer array)
c                                   istat(i,1) = 0 ==> active:  = 1 not
c                                   istat(i,2) = -1 ==> Danby did not work
c               x            ==> initial value independent variable 
c                                        (real scalar)
c 	        h0           ==> stepsize  (real scalar)
c 	        y            ==> initial value dependent variables  
c                                     (real array)
c               eps          ==> local truncation error tolerance
c
c             Output:
c 	          y  ==> final value dependent variables  (real array)
c                 x  ==> final value independent variable (real scalar)
c 	          h0 ==> recommended stepsize for the next call (real scalar)
c
c Remarks:  Based on Renu's code: mass,j2rp2,j4rp4 and istat are 
c           just passed on to bs_der
c Authors:  Hal Levison
c Date:    5/17/93
c Last revision: 2/24/94

      subroutine bs_int(nbod,ntp,mass,j2rp2,j4rp4,istat,x,h0,y,eps)

      include '../swift.inc'
      include 'bs.inc'

c...  Inputs Only: 
      integer nbod,ntp
      real*8 mass(nbod),h0,eps,j2rp2,j4rp4

c...  Input & Output
      real*8 x,y(N6DBS)
      integer istat(NTPMAX,NSTAT)

c...  Internals
      real*8 tp(NTEMP),dy(N6DBS),d(6),alt(10),lt(10)
      integer idiv,i,ii,m,l,m1,k,mmk,i1,i1max,ik,n
      real*8 xa,xb,varm,fl,h,hd,flt,eta2,dta,yb,c,b1
      real*8 den,dtn,b,var,varma
      logical lbreak

      data lt/1,2,3,4,6,8,12,16,24,32/
      data alt/1.d0,2.d0,3.d0,4.d0,6.d0,8.d0,12.d0,16.d0,24.d0,32.d0/

      save lt,alt

c----
c...  Executable code 

      n = 6*(nbod+ntp)

      xa=x
      call bs_der(ntp,nbod,mass,j2rp2,j4rp4,istat,y,dy)
      do i=1,n
         ii=12*i
         tp(ii-1)=dabs(y(i))
c
         if(tp(ii-1).lt.eps) then
            tp(ii-1)=eps
         endif
c
         tp(ii-4)=dy(i)
         tp(ii)=y(i)
      enddo

      do idiv=0,NTRYS

         xb=h0+xa
c
c        successive extrapolations
c
c         do m=1,10
         m = 1
         lbreak = .true.
         do while( (m.le.10) .and. lbreak )

            l=lt(m)
            fl=alt(m)
            varm=0.d0
            m1=min0(m-1,6)
c
c           calculation of d(k)=(h(m-k)/h(m))**2 for equation (6)
c
            if(m1.ne.0) then
               do k=1,m1
                  mmk=m-k
                  flt=alt(mmk)
                  d(k)=(fl/flt)**2
               enddo
            endif
            h=h0/fl
            hd=0.5d0*h
c
c           integration
c
            do i=1,n
               ii=12*i
               tp(ii-3)=tp(ii) 
               y(i)=tp(ii)+hd*tp(ii-4) !    equation (3b)
            enddo
            i1max=2*l-1
            x=xa
         
            do i1=1,i1max
               x=x+hd
               call bs_der(ntp,nbod,mass,j2rp2,j4rp4,istat,y,dy)
               do i=1,n
                  ii=12*i
                  tp(ii-1)=dmax1(tp(ii-1),dabs(y(i)))
                  eta2=tp(ii-3)+h*dy(i)
                  tp(ii-3)=y(i)
                  y(i)=eta2
               enddo 
            enddo
         
            call bs_der(ntp,nbod,mass,j2rp2,j4rp4,istat,y,dy)
            do i=1,n
               ii=12*i
               dta=tp(ii-11)
               yb=(tp(ii-3)+y(i)+hd*dy(i))/2.d0 !    equation (3d)
c     
c              extrapolated values
c
               c=yb             ! equation (6b)
               tp(ii-11)=yb     ! equation (6a)

               if(m1.ne.0) then
                  do k=1,m1
                     b1=d(k)*dta
                     den=b1-c
                     dtn=dta
                     if(den.ne.0.d0) then
                        b=(c-dta)/den
                        dtn=c*b       !   equation (6c)
                        c=b1*b        !    equation (6d)
                     endif
                     ik=ii-11+k
                     dta=tp(ik)
                     tp(ik)=dtn       !   equation (6c)
                     yb=yb+dtn        !     equation (6f)
                  enddo
                  var=dabs(tp(ii-2)-yb)/tp(ii-1)
                  varm=dmax1(varm,var)
               endif
               tp(ii-2)=yb
            enddo
            
            if(m.gt.3) then
               if(varm.le.eps) then       !   return results to calling program
                  x=xb
                  do i=1,n
                     ii=12*i
                     y(i)=tp(ii-2)
                  enddo
                  h0=h0*1.5d0*0.6d0**(m-1-m1)    !   recommend a new step size
                  return
               endif

               if(varm.ge.varma) then !  calculation did not converge
c                                        start again with half the step size
                  lbreak = .false.
               endif
            endif
            varma=varm
            m = m + 1
         enddo                  ! m

         h0=h0/2.d0
      enddo      ! idiv

      write(*,*) ' ERROR (b_int): lack of convergence !!!! '

c
      end      !  bs_int
c-----------------------------------------------------------------------------




