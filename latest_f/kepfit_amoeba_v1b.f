ccc   By Trifon Trifonov trifon@hku.hk
ccc   You can modify if as you want, but please
ccc   do not distribute without permision. 
ccc   This is not a final version!!! 
ccc   The final version will be available in the Python RVMod lib
ccc   Trifonov et al. (in prep).    

      implicit none
      real*8 PI, twopi
      parameter (PI=3.14159265358979d0)
      integer npl,ndset,idset,ndata,ma,mfit,i,j,NDSMAX,NPLMAX,MMAX
      integer writeflag_best_par
      integer writeflag_RV,writeflag_fit, amoebastarts
      parameter (NDSMAX=20, NPLMAX=20, MMAX=200)
      integer idsmax(NDSMAX),ia(MMAX),nt, ts(5000),ii, iter
      real*8 x(5000),y(5000),sig(5000),y_in(5000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 rms,mstar, mass(NPLMAX),ap(NPLMAX)
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*8 chisq,alamda,ochisq,dchisq, epsil, deltat
      real*8 sigscale,x0,xmax, incl(NPLMAX),cap0m(NPLMAX)
      real*8 t0,t1,t2,dt,offset,t_max,loglik,dy,sig2i
      real*8 st_mass,sini,m1,a1,m2,a2,epoch,ftol
      real*8 ymod(5000),dyda(MMAX), p(MMAX+1,MMAX),yamoeba(MMAX+1)
      real*8 loglikk, ologlikk, dloglikk 
      external rvkep, compute_abs_loglik
      character*80 infile
      real*4 t_stop,when_to_kill, model_max
      
      
      common /DSBLK/ npl,ndset,idsmax,idset

      twopi=2.d0*PI
      ftol=0.000001d0

c     first two just for consistency with dynamical input, not really used
      read (*,*) epsil,deltat, amoebastarts,
     &          when_to_kill, nt, model_max   
     
     
c      write(*,*) 'Stellar mass'
      read (*,*) st_mass, writeflag_best_par, writeflag_RV,
     & writeflag_fit 
 
      call io_read_data (ndata,x,ts,y,sig,epoch,
     &               x0,t_max,a,ia,ma,incl,cap0m)
       
      mfit = 0
      do j = 1,ma
          if (ia(j).ne.0) mfit = mfit + 1
      enddo
      
c      call prepare_for_amoeba(p,MMAX+1,MMAX,yamoeba,a,ia,ma,mfit,
c     & compute_abs_loglik,ndata,x,y,ymod,dyda,ts,sig)
 

      i = 0
 500  continue
 
         if (i.eq.amoebastarts) then
             i = 0
             goto 502
         endif

 
         i = i + 1
         ologlikk = loglikk

         call prepare_for_amoeba(p,MMAX+1,MMAX,yamoeba,a,ia,ma,mfit,
     & compute_abs_loglik,ndata,x,y,ymod,dyda,ts,sig, i)
         call amoeba(p,yamoeba,MMAX+1,MMAX,mfit,ftol,compute_abs_loglik,
     & iter,ndata,x,y,ymod,dyda,ma,ts,sig,a,ia,loglikk)
 
 
         CALL SECOND(t_stop)
         if (t_stop.ge.when_to_kill) then
            write(*,*) 'Max. time=',when_to_kill, 'sec ', 
     &                 'exceeded t_stop =', t_stop, 'sec ' 
            goto 502
         endif 
 
 
         loglikk = yamoeba(1)
         dloglikk = ologlikk - loglikk

         j=0
         do ii=1,ma
           if (ia(ii).ne.0) then
              j=j+1
              a(ii)=p(1,j)  
           endif
         enddo

      if (dabs(dloglikk).ge.0.000001d0) goto 500

     
502   idset = 1
      chisq=0.d0
      loglik=0.d0
      
      do i = 1,npl
         j = 5*(i-1)
         
         if (a(j+2).lt.0.d0) then  ! if P<0, set P>0 
            a(j+2) = abs(a(j+2))
         endif         
         
         if (a(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
            a(j+4) = a(j+4) + PI
            a(j+1) = abs(a(j+1))
            if (a(j+4).gt.2.d0*PI) a(j+4) = a(j+4)-2.d0*PI
         endif
         if (a(j+3).lt.0.d0) then  ! if e<0, set e>0 and w=w+PI, M0=M0-PI
            a(j+3) = abs(a(j+3))
            a(j+4) = a(j+4) +  PI
            if (a(j+4).gt.2.d0*PI) a(j+4) = a(j+4)-2.d0*PI
            a(j+5) = a(j+5) - PI
            if (a(j+5).lt.0.d0) a(j+5) = a(j+5)+2.d0*PI
         endif
         if(a(j+3).ge.1.d0) then ! if e>=1 set it to 0.99 to prevent errors
             a(j+3)=0.99d0
         endif
         if (a(j+4).lt.0.d0) a(j+4) = dmod(a(j+4)+2.d0*PI,  2.d0*PI)  
         if (a(j+5).lt.0.d0) a(j+5) = dmod(a(j+5)+2.d0*PI,  2.d0*PI) 
         if (a(j+4).gt.2.d0*PI) a(j+4) = dmod(a(j+4),  2.d0*PI )  
         if (a(j+5).gt.2.d0*PI) a(j+5) = dmod(a(j+5),  2.d0*PI )         
                              
      enddo        
      


      do i = 1,ndata

          idset = ts(i)
          call RVKEP (x(i),a,ymod(i),dyda,ma,idset)

          y_in(i) = y(i) - a(5*npl+idset)- a(5*npl+2*ndset+1)*x(i)
	  ymod(i) = ymod(i) - a(5*npl+idset) - 
     &    a(5*npl +2*ndset + 1)*x(i)

          dy = y_in(i) - ymod(i)

          if (writeflag_RV.gt.0) then 
              write(*,*) x0 + x(i),
     &        ymod(i), y_in(i),
     &        dy, sig(i), idset
   
          endif

          sig2i = 1.d0/(sig(i)**2 + a(5*npl+ndset+idset)**2)

 	  chisq  = chisq + dy*dy*sig2i

	  loglik =  loglik - 0.5*dy*dy*sig2i -
     &               0.5*dlog(twopi*(sig(i)**2
     &                + a(5*npl+ndset+idset)**2))  
          rms = rms + dy**2
      enddo


      rms = dsqrt(rms/dble(ndata))

      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
      write(*,*) loglik, chisq/dble(ndata-mfit),chisq,rms


51    format(f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3)
52    format(a,f14.3)
53    format(a,i4,a,i4,a,f7.3,a,f7.3,a,f12.3)    
    
      call MA_J (a,ma,npl,st_mass,sini,mass,ap)    
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (writeflag_best_par.gt.0) then

          do j = 1,npl+1
	      j_mass(j) = mass(j)/1.2667d17 
c             s_mass(j) = mass(j)/1.32712497d20 
c              swift_mass(j) = (mass(j)/1.32712497d20)*((4.d0*PI*PI)
c     &           /(365.25*365.25))
          enddo

          write (*,*) 'Best-fit K [m/s], P [days], e, w [deg], 
     & M0 [deg], i[deg], cap0m[deg] and their errors'

          do j = 1,npl
              i = 5*(j-1)
              write (*,*) a(i+1),a(i+2),a(i+3),a(i+4)*180.d0/PI,
     &        a(i+5)*180.d0/PI,incl(j),cap0m(j)
              write (*,*) dsqrt(covar(i+1,i+1)),dsqrt(covar(i+2,i+2)),
     &                 dsqrt(covar(i+3,i+3)), 
     &                 dsqrt(covar(i+4,i+4))*180.d0/PI,
     &                 dsqrt(covar(i+5,i+5))*180.d0/PI, 0.0, 0.0

          enddo
          write (*,*) 'Best-fit V0 [m/s] and their error bars:'
          do j = 1,ndset
              i = 5*npl + j
              write (*,*) a(i)
              write (*,*) dsqrt(covar(i,i))
          enddo

          write (*,*) 'Jitters for each data set:'
          do j = 1,ndset
              write (*,*) a(5*npl+ndset+j)
              write (*,*) '0'
          enddo          
          
          write (*,*) 'linear trend [m/s per day]:'
          write (*,*) a(5*npl + 2*ndset + 1)  
          write (*,*) dsqrt(covar(5*npl + ndset + 1,5*npl + ndset + 1))              
          
          write (*,*) ' ndata =',ndata
          write (*,*) ' mfit =',mfit
          write (*,*) ' RMS =',rms
          write (*,*) ' Chi^2 =',chisq/dble(ndata-mfit)
          write (*,*) ' epoch = ', x0

          write (*,*) 'Jupiter mass'
          write (*,*) (j_mass(i+1),i=1,npl)
 
          write(*,*) 'semi-major axes in Jacobi'
          write(*,*)  (ap(i)/1.49597892d11,i=1,npl)
 
      endif

 
      if(writeflag_fit.gt.0) then
          dt = (x(ndata)+model_max)/dble(nt - 1)
          do i = 1,nt
	      x(i) = ((i-1)*dt)-0.00

              do j = 1,ndset
                  a(5*npl + j) = 0.0
              enddo
              call RVKEP (x(i),a,ymod(i),dyda,ma,1)
              write(*,*) x0 + x(i), ymod(i)
          enddo
      endif


      stop
      end




      subroutine compute_abs_loglik(ndata,x,y,a2,ymod,dyda,ma,mfit,ts,
     & sig,loglik, num,a,ia)
      implicit none
      
      integer MMAX,NDSMAX,npl,ndset,idset,num, mfit    
      parameter (MMAX=200, NDSMAX=20)
      real*8 twopi, loglik
      parameter (twopi=6.28318530717958d0)
      integer ndata, i, j, ma, ts(5000), ia(MMAX), idsmax(NDSMAX)
      real*8 dy, sig(5000), dyda(MMAX), x(5000), y(5000)
      real*8 ymod(5000), a(MMAX), a2(mfit), a3(MMAX),sig2i, y_in(5000)
     & , y2(5000)
      
   
      common /DSBLK/ npl,ndset,idsmax,idset      
      
      loglik=0.d0
      j=0
      do i=1,ma 
         if (ia(i).ne.0) then
             j=j+1
             a3(i)=a2(j)
         else
             a3(i)=a(i)
         endif
      enddo
        do i = 1,ndata
              idset = ts(i)
              call RVKEP (x(i),a3,y2(i),dyda,ma,idset)
              y_in(i) = y(i) - a3(5*npl+idset)- 
     &                 a3(5*npl+2*ndset+1)*x(i)
	      y2(i) = y2(i) - a3(5*npl+idset) - 
     &        a3(5*npl +2*ndset + 1)*x(i)

          dy = y_in(i) - y2(i)

	  sig2i = 1.d0/(sig(i)**2 + a3(5*npl+ndset+idset)**2)

	  loglik =  loglik + 0.5*dy*dy*sig2i +
     &               dlog(dsqrt(twopi*(sig(i)**2 + 
     &               a3(5*npl+ndset+idset)**2))) 
     &               - dlog(dsqrt(twopi)) 
        enddo
     
      return
      end      

      subroutine io_read_data (ndata,t,ts,ys,sigs,epoch,t0,t_max,
     &          ar,iar,ma,incl,cap0m)  

      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl,ma
      
   
      real*8 t(5000),y(5000),sig(5000),ys(5000),sigs(5000),PI
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      parameter(PI=3.14159265358979d0)
      real*8 ar(MMAX),incl(NPLMAX),cap0m(NPLMAX)
      integer iar(MMAX),u_off(NDSMAX),u_jit(NDSMAX)
      integer idsmax(NDSMAX),ts(5000), u_incl, u_cap0m
      real*8 jitt(NDSMAX),sigscale,t0,t_max, epoch
      real*8 off(NDSMAX),loglik
      integer i,k,j
      character*80 infile
   
      common /DSBLK/ npl,ndset,idsmax,idset

c      write (*,*)  ' Number of Data Sets: '
      read (*,*) ndset
      if (ndset.gt.NDSMAX) stop ' KEPFIT: ndset > NDSMAX.'

      ndata = 1
      do i = 1,ndset
          read (*,50) infile
 50       format (a)
          read (*,*) off(i)
          read (*,*) u_off(i)
          read (*,*) jitt(i)
          read (*,*) u_jit(i)
          open (unit=10,file=infile)
 100      continue
             read (10,*,err=200,end=200) t0,y(ndata),sig(ndata)

c              sig(ndata) = dsqrt(sig(ndata)**2 + jitt(i)**2)
              if (ndata.eq.1) then              ! make sequences for datapoints
                 t(ndata) = t0
                 ts(ndata) = i
                 ys(ndata) = y(ndata)
                 sigs(ndata) = sig(ndata)              
              else
                 j = 1
 1002            if (t0.lt.t(j)) then
                    do k = ndata-1,j,-1
                       t(k+1) = t(k)
                       ts(k+1) = ts(k)
                       ys(k+1) = ys(k)
                       sigs(k+1) = sigs(k)
                    enddo
                    t(j) = t0
                    ts(j) = i
                    ys(j) = y(ndata)
                    sigs(j) = sig(ndata)
                    goto 1001
                 else
                    j = j + 1
                    if (j.eq.ndata) then
                       t(j) = t0
                       ts(j) = i
                       ys(j) = y(ndata)
                       sigs(j) = sig(ndata)
                       goto 1001
                    endif
                    goto 1002
                 endif
              endif

 1001         ndata = ndata + 1
          goto 100
 200      idsmax(i) = ndata - 1
          close (unit=10)
                    
          
          enddo
          read (*,*) npl
      if (npl.gt.NPLMAX) stop ' KEPFIT: npl > NPLMAX.'

      do i = 1,ndset
          ar(5*npl+i)=off(i)
          iar(5*npl+i)=u_off(i)
          ar(5*npl+ndset+i)=jitt(i)
          iar(5*npl+ndset+i)=u_jit(i)
      enddo
      
      ma = 5*npl + 2*ndset + 1
      do j = 1,npl
          i = 5*(j-1)
          read (*,*) ar(i+1),ar(i+2),ar(i+3),ar(i+4),ar(i+5),incl(j)
     &    ,cap0m(j)
          read (*,*) iar(i+1),iar(i+2),iar(i+3),iar(i+4),iar(i+5)
     &    ,u_incl, u_cap0m

c         inclinations and cap0m are always ignored in the fit, just for consistency with dynamical input and output

          ar(i+4) = ar(i+4)*PI/180.d0
          ar(i+5) = ar(i+5)*PI/180.d0
      enddo

          
      read (*,*) ar(5*npl+ 2*ndset+1)
      read (*,*) iar(5*npl+2*ndset+1)    

      ndata = ndata - 1


c      write(*,*) 'for epoch :'
      read (*,*) epoch
 
      t_max = t(ndata) 

      if (epoch.eq.0) then 
         t0 = t(1)
      else
         t0 = epoch
      endif
         



      do i = 1,ndata
         t(i) = (t(i) - t0)             ! time unit is day
      enddo


      return
      end      
 

      

      subroutine RVKEP (x,a,y,dyda,ma,ts)
      implicit none
      real*8 PI,TWOPI
      parameter (PI=3.14159265358979d0)
      parameter (TWOPI=2.0d0*PI)
      integer npl,ndset,idset,ma,i,j,NDSMAX,ts
      parameter (NDSMAX=20)
      integer idsmax(NDSMAX)
      real*8 x,y,a(ma),a2(ma),dyda(ma)
      real*8 cosw,sinw,capm,cape,cose,sine,cosf,sinf,fac1,fac2,fac3
      real*8 orbel_ehybrid, f, coswf

      common /DSBLK/ npl,ndset,idsmax,idset

      y = 0.d0
 
      do i = 1,ma
          a2(i)=a(i)
      enddo
 
      do i = 1,npl
         j = 5*(i-1)
         
         if (a2(j+2).lt.0.d0) then  ! if P<0, set P>0 
            a2(j+2) = abs(a2(j+2))
         endif         
         
         if (a2(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
            a2(j+4) = a2(j+4) + PI
            a2(j+1) = abs(a2(j+1))
            if (a2(j+4).gt.2.d0*PI) a2(j+4) = a2(j+4)-2.d0*PI
         endif
         if (a2(j+3).lt.0.d0) then  ! if e<0, set e>0 and w=w+PI, M0=M0-PI
            a2(j+3) = abs(a2(j+3))
            a2(j+4) = a2(j+4) +  PI
            if (a2(j+4).gt.2.d0*PI) a2(j+4) = a2(j+4)-2.d0*PI
            a2(j+5) = a2(j+5) - PI
            if (a2(j+5).lt.0.d0) a2(j+5) = a2(j+5)+2.d0*PI
         endif
         if(a2(j+3).ge.1.d0) then ! if e>=1 set it to 0.99 to prevent errors
             a2(j+3)=0.99d0
         endif
         if (a2(j+4).lt.0.d0) a2(j+4) = dmod(a2(j+4)+2.d0*PI,  2.d0*PI)  
         if (a2(j+5).lt.0.d0) a2(j+5) = dmod(a2(j+5)+2.d0*PI,  2.d0*PI) 
         if (a2(j+4).gt.2.d0*PI) a2(j+4) = dmod(a2(j+4),  2.d0*PI )  
         if (a2(j+5).gt.2.d0*PI) a2(j+5) = dmod(a2(j+5),  2.d0*PI )         
         
                        
      enddo  

      do j = 1,npl

          i = 5*(j-1)
          cosw = dcos(a2(4+i))
          sinw = dsin(a2(4+i))

          capm = TWOPI*x/a2(2+i) + a2(5+i)
          capm = dmod(capm,  2.d0*PI )

          cape = ORBEL_EHYBRID (a2(3+i),capm)
          cose = dcos(cape)
          sine = dsin(cape)
          
          cosf = (cose - a2(3+i))/(1.d0 - a2(3+i)*cose)
          sinf = (dsqrt(1.d0 - a2(3+i)**2)*sine)/(1.d0 - a2(3+i)*cose)
c          f = 2.0d0*datan2( dsqrt(1.d0 - a(3+i))*dcos(cape/2.d0),
c     &                 dsqrt(1.d0 + a(3+i))*dsin(cape/2.d0))

c          coswf = dcos(a(4+i)+f)
c          fac1 = coswf + a(3+i)*cosw

          fac1 = cosw*cosf - sinw*sinf + a2(3+i)*cosw

          fac2 = (cosw*sinf + sinw*cosf)/(1.d0 - a2(3+i)*cose)**2
          fac3 = -a2(1+i)*dsqrt(1.d0 - a2(3+i)**2)*fac2

          y = y + a2(1+i)*fac1
          dyda(1+i) = fac1
          dyda(2+i) = -TWOPI*fac3*x/a2(2+i)**2
          dyda(3+i) = -a2(1+i)*sine*(2.d0 - a2(3+i)**2 - a2(3+i)*cose)*
     &                 fac2/dsqrt(1.d0 - a2(3+i)**2)
          dyda(4+i) = -a2(1+i)*(sinw*cosf + cosw*sinf)
          dyda(5+i) = fac3

      enddo

c      do i = 1,idset
      y = y + a2(5*npl+ts)
      dyda(5*npl+ts) = 1.d0
c      enddo

c      do i = 1,idset
c          y = y + a(5*npl+i)
c          dyda(5*npl+i) = 1.d0
c      enddo


      y = y + a2(5*npl +2*ndset + 1)*x  
      dyda(5*npl + ndset + 1) = x
   
      do i = ts+1,ndset
          dyda(5*npl+i) = 0.d0
      enddo

      return
      end

      subroutine prepare_for_amoeba(p,mp,np,y,a,ia,ma,mfit,funk,ndata,
     & x,z,ymod,dyda,ts,sig, it)
      integer MMAX,NDSMAX,ma,ts(5000), ndata,mp,np,mfit,it
      parameter(MMAX=200,NDSMAX=20)
      REAL*8 ftol,p(mp,np),y(mp),a(MMAX), a2(mfit),fr,frjitt
      real*8 x(5000),z(5000),ymod(5000)
      real*8 dyda(MMAX), sig(5000), loglik
      parameter(fr=0.05, frjitt=0.05)
      INTEGER i,j,k, ia(MMAX), idsmax(NDSMAX)
      external funk
    
      common /DSBLK/ npl,ndset,idsmax,idset

      k=0
      do j=1,ma
          if(ia(j).ne.0) then
          k=k+1
          p(1,k)=a(j)
          do i=2,mfit+1
              if (k.eq.(i-1)) then
                  if (j.gt.(5*npl+ndset)) then
                  p(i,k)=(1+frjitt)*(p(1,k)+0.1)
                  else 
                  if (mod(j,5).eq.2) then
                     p(i,k)=(1+fr)*(p(1,k) + 0.1)
                  else if (mod(j,5).eq.3) then
                     p(i,k)=(1+frjitt)*(p(1,k)+0.1)
                  else
                     p(i,k)=(1+fr)*(p(1,k)+0.1)
                  endif
                  endif
              else
                  p(i,k)=p(1,k)
              endif
          enddo
          endif
      enddo
      do i=1,mfit+1
          do j=1,mfit
              a2(j)=p(i,j)
          enddo
          call funk(ndata,x,z,a2,ymod,dyda,ma,mfit,ts,sig,loglik,i,
     &              a,ia)
          y(i)=loglik
      enddo
      return
      end
          
          
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,ytry)
      implicit none
      INTEGER iter,mp,ndim,np,NMAX,ITMAX, MMAX,ma,ts(5000), ndata
      REAL*8 ftol,p(mp,np),y(mp),x(5000),z(5000),ymod(5000)
      PARAMETER (NMAX=20,ITMAX=50000,MMAX=200)
      real*8 dyda(MMAX), sig(5000), loglik, a(MMAX)
      EXTERNAL funk
      INTEGER i,ihi,ilo,inhi,j,m,n, ia(MMAX)
      REAL*8 rtol,summ,swap,ysave,ytry,psum(ndim),amotry
      iter=0
1     do 12 n=1,ndim
        summ=0.d0
        do 11 m=1,ndim+1
          summ=summ+p(m,n)
11      continue
        psum(n)=summ
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
      if (iter.ge.ITMAX) then
          write (*,*) 'ITMAX exceeded in amoeba'
          return
      endif
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia)
 
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              call funk(ndata,x,z,psum,ymod,dyda,ma,ndim,ts,sig,loglik,
     &                 i,a,ia)
              y(i)=loglik
             endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2

      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0=M,173+9.

      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia)
      implicit none
      INTEGER ihi,mp,ndim,np,NMAX, MMAX, ma, ts(5000),ndata
      PARAMETER (NMAX=20, MMAX=200)
      REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),x(5000),z(5000),
     & ymod(5000)
      real*8 dyda(MMAX), sig(5000),loglik
      EXTERNAL funk
      INTEGER j, ia(MMAX)
      REAL*8 fac1,fac2,ytry,ptry(ndim), a(MMAX)
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      call funk(ndata,x,z,ptry,ymod,dyda,ma,ndim,ts,sig,loglik,ihi,
     & a,ia)

      ytry=loglik
C      WRITE(*,*) loglik
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0=M,173+9.




	subroutine MA_J (a,ma,npl,m0,sini,mass,ap)
        
	implicit none
	real*8 m0,PI,TWOPI,THIRD,GMSUN,dm,MSUN
        integer npl,ma,i,j,NPLMAX
        parameter (NPLMAX=7)
        real*8 sini,mm(NPLMAX)
        real*8 a(ma),mass(NPLMAX),ap(NPLMAX),mpold(NPLMAX),mtotal
	parameter (THIRD=1.d0/3.d0)
        parameter (PI=3.14159265358979d0,TWOPI=2.d0*PI)
	parameter (GMSUN=1.32712497d20,MSUN=1.32712497d20)

c*******G is set to be unit, and s, m, kg as unit of time, length and mass
c*******expectively.        
        
        do j = 1,npl
           i = 5*(j-1) 

           mm(j) = 2.d0*PI/(a(i+2)*8.64d4)
           
        enddo 
 
        do i = 0,npl-1

           mass(1) = m0
	   mpold(i+1) = 0.d0
 101       continue
           if (i.eq.0) then
           mtotal = m0
	   mass(i+2) = a(5*i+1)*(TWOPI/mm(i+1)*(m0 + mpold(i+1))**2/
     &               (TWOPI*GMSUN))**THIRD*
     &               dsqrt(1.d0 - a(5*i+3)**2)
           else
              mtotal = m0
              do j = 0, i-1
                 mtotal = mtotal + mass(j+2)
              enddo
              mass(i+2) = a(5*i+1)*(TWOPI/mm(i+1)*(mtotal
     &                  +mpold(i+1))**2/(TWOPI*GMSUN))**THIRD*
     &                  dsqrt(1.d0 - a(5*i+3)**2)
           endif
           
	   dm = dabs(mass(i+2)-mpold(i+1))/mass(i+2)
	   mpold(i+1) = mass(i+2)
           if (dm.gt.0) goto 101

	   ap(i+1) = (GMSUN*(mtotal + mass(i+2))*(1.d0/mm(i+1))
     &               **2)**THIRD
           
        enddo

        do i = 1,npl+1
           mass(i) = mass(i)*MSUN
        enddo

        
	return
	end
