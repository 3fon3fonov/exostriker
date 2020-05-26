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
      integer writeflag_best_par,hkl,gr_flag
      integer writeflag_RV,writeflag_fit, amoebastarts
      parameter (NDSMAX=20, NPLMAX=20, MMAX=200)
      integer idsmax(NDSMAX),ia(MMAX),nt, ts(20000),ii, iter
      real*8 x(20000),y(20000),sig(20000),y_in(20000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 rms,mstar, mass(NPLMAX),ap(NPLMAX)
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*8 chisq,alamda,ochisq,dchisq, epsil, deltat
      real*8 sigscale,x0,xmax, incl(NPLMAX),cap0m(NPLMAX)
      real*8 t0,t1,t2,dt,offset,t_max,loglik,dy,sig2i
      real*8 st_mass,sini,m1,a1,m2,a2,epoch,ftol
      real*8 ymod(20000),dyda(MMAX), p(MMAX+1,MMAX),yamoeba(MMAX+1)
      real*8 loglikk, ologlikk, dloglikk,best_w,best_we
      external rvkep, compute_abs_loglik
      character*80 infile
      character*80 version_input, version
      
      real*4 t_stop,when_to_kill, model_max,model_min
      
      
      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag


      version = "0.07"
       
      CALL getarg(1, version_input)     
      if(version_input.eq.'-version') then
          write(*,*) version
          goto 222
      endif

      twopi=2.d0*PI
      ftol=0.000001d0

c     first two just for consistency with dynamical input, not really used
      read (*,*) epsil,deltat, amoebastarts,
     &          when_to_kill, nt, model_max, model_min ,gr_flag 
     
     
c      write(*,*) 'Stellar mass'
      read (*,*) st_mass, writeflag_best_par, writeflag_RV,
     & writeflag_fit 
 
      call io_read_data (ndata,x,ts,y,sig,epoch,
     &               x0,t_max,a,ia,ma,incl,cap0m,hkl)
 
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
     & compute_abs_loglik,ndata,x,y,ymod,dyda,ts,sig, i,hkl)
         call amoeba(p,yamoeba,MMAX+1,MMAX,mfit,ftol,compute_abs_loglik,
     & iter,ndata,x,y,ymod,dyda,ma,ts,sig,a,ia,loglikk,hkl)
 
 
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
      
      
      if (hkl.eq.0) then

          do i = 1,npl
             j = 6*(i-1)
             
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
             if (a(j+4).lt.0.d0) a(j+4) = dmod(a(j+4)+2.d0*PI, 2.d0*PI)  
             if (a(j+5).lt.0.d0) a(j+5) = dmod(a(j+5)+2.d0*PI, 2.d0*PI) 
c             if (a(j+6).lt.0.d0) a(j+6) = dmod(a(j+6)+2.d0*PI, 2.d0*PI) 
             
             if (a(j+4).gt.2.d0*PI) a(j+4) = dmod(a(j+4), 2.d0*PI )  
             if (a(j+5).gt.2.d0*PI) a(j+5) = dmod(a(j+5), 2.d0*PI )   
c             if (a(j+6).gt.2.d0*PI) a(j+6) = dmod(a(j+6), 2.d0*PI )   
           
 
                                                   
          enddo  

      else   
            
          do i = 1,npl
             j = 6*(i-1)
             if (a(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
                a(j+4) = -1.d0*a(j+4)       !     which is h = -h, k = -k
                a(j+3) = -1.d0*a(j+3)
                a(j+1) = abs(a(j+1))    
             endif
          
              
             if (a(j+5).lt.0.d0) a(j+5) = dmod(a(j+5)+2.d0*PI, 2.d0*PI) 
             if (a(j+5).gt.2.d0*PI) a(j+5) = dmod(a(j+5), 2.d0*PI )   
             if (a(j+6).lt.0.d0) a(j+6) = dmod(a(j+6)+2.d0*PI, 2.d0*PI) 
             if (a(j+6).gt.2.d0*PI) a(j+6) = dmod(a(j+6), 2.d0*PI )
             
c             write(*,*) a(j+4),a(j+4),ecc(i) ,omega(i) ,capmm(i) 
          enddo        
         
      endif     
      
      
      do i = 1,ndata

          idset = ts(i)
          call RVKEP (x(i),a,ymod(i),dyda,ma,idset,hkl)

          y_in(i) = y(i) - a(6*npl+idset) - a(6*npl+2*ndset+1)*x(i) -
     &      a(6*npl +2*ndset + 2)*x(i)**2
	      ymod(i) = ymod(i) - a(6*npl+idset) 
     &    - a(6*npl +2*ndset + 1)*x(i) - 
     &      a(6*npl +2*ndset + 2)*x(i)**2

          dy = y_in(i) - ymod(i)

          if (writeflag_RV.gt.0) then 
              write(*,*) x0 + x(i),
     &        ymod(i), y_in(i) + a(6*npl+2*ndset+1)*x(i) +
     &        a(6*npl +2*ndset + 2)*x(i)**2,
     &        dy, sig(i), idset
   
          endif

          sig2i = 1.d0/(sig(i)**2 + a(6*npl+ndset+idset)**2)

 	      chisq  = chisq + dy*dy*sig2i

	      loglik =  loglik - 0.5*dy*dy*sig2i -
     &               0.5*dlog(twopi*(sig(i)**2
     &                + a(6*npl+ndset+idset)**2))  
          rms = rms + dy**2
      enddo


      rms = dsqrt(rms/dble(ndata))

      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
      write(*,*) loglik, chisq/dble(ndata-mfit),chisq,rms


51    format(f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3)
52    format(a,f14.3)
53    format(a,i4,a,i4,a,f7.3,a,f7.3,a,f12.3)    
    
      call MA_J (a,ma,npl,st_mass,sini,mass,ap,hkl,gr_flag)    
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (writeflag_best_par.gt.0) then

          do j = 1,npl+1
	      j_mass(j) = mass(j)/1.2667d17 
c             s_mass(j) = mass(j)/1.32712497d20 
c              swift_mass(j) = (mass(j)/1.32712497d20)*((4.d0*PI*PI)
c     &           /(365.25*365.25))
          enddo

          write (*,*) 'Best-fit K [m/s], P [days], e, w [deg], 
     & M0 [deg], i[deg], cap0m[deg], w dot [deg/yr], and their errors'

          do j = 1,npl
              i = 6*(j-1)
              
              if (hkl.eq.0) then
                  best_w = a(i+4)*180.d0/PI
                  best_we = dsqrt(covar(i+4,i+4))*180.d0/PI
              else    
                  best_w = a(i+4) 
                  best_we = dsqrt(covar(i+4,i+4)) 
              endif    
                                
              write (*,*) a(i+1),a(i+2),a(i+3),best_w,
     &        a(i+5)*180.d0/PI,incl(j),cap0m(j),a(i+6)*180.d0/PI
              write (*,*) dsqrt(covar(i+1,i+1)),dsqrt(covar(i+2,i+2)),
     &                 dsqrt(covar(i+3,i+3)), 
     &                 best_we,
     &                 dsqrt(covar(i+5,i+5))*180.d0/PI,0.d0, 0.d0,
     &                 dsqrt(covar(i+6,i+6))*180.d0/PI

          enddo
          
          write (*,*) 'Best-fit V0 [m/s] and their error bars:'
          
          do j = 1,ndset
              i = 6*npl + j
              write (*,*) a(i)
              write (*,*) dsqrt(covar(i,i))
          enddo

          write (*,*) 'Jitters for each data set:'
          do j = 1,ndset
              write (*,*) a(6*npl+ndset+j)
              write (*,*) '0'
          enddo          
          
          write (*,*) 'linear trend [m/s per day]:'
          write (*,*) a(6*npl + 2*ndset + 1)  
          write (*,*) dsqrt(covar(6*npl + ndset + 1,6*npl + ndset + 1))    
          
          write (*,*) 'quad. trend [m/s per day]:'
          write (*,*) a(6*npl + 2*ndset + 2)  
          write (*,*) dsqrt(covar(6*npl + ndset + 2,6*npl + ndset + 2))           
                    
          
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
          dt = (x(ndata)+model_max+model_min)/dble(nt - 1)
          do i = 1,nt
	      x(i) = ((i-1)*dt)-model_min

              do j = 1,ndset
                  a(6*npl + j) = 0.0
              enddo
              call RVKEP (x(i),a,ymod(i),dyda,ma,1,hkl)
              write(*,*) x0 + x(i), ymod(i)
          enddo
      endif


c      stop
222   end




      subroutine compute_abs_loglik(ndata,x,y,a2,ymod,dyda,ma,mfit,ts,
     & sig,loglik, num,a,ia,hkl)
      implicit none
      
      integer MMAX,NDSMAX,npl,ndset,idset,num, mfit,gr_flag
      parameter (MMAX=200, NDSMAX=20)
      real*8 loglik, PI, TWOPI
      parameter (PI=3.14159265358979d0)
      parameter (TWOPI=2.0*PI)  
      integer ndata, i, j, ma, ts(20000), ia(MMAX), idsmax(NDSMAX),hkl
      real*8 dy, sig(20000), dyda(MMAX), x(20000), y(20000)
      real*8 ymod(20000),a(MMAX),a2(mfit),a3(MMAX),sig2i,y_in(20000)
     & , y2(20000)
      
   
      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag      
      
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
              call RVKEP (x(i),a3,y2(i),dyda,ma,idset,hkl)
              y_in(i) = y(i) - a3(6*npl+idset)- 
     &                 a3(6*npl+2*ndset+1)*x(i) 
     &               - a3(6*npl+2*ndset+2)*x(i)**2
     
	          y2(i) = y2(i) - a3(6*npl+idset) - 
     &         a3(6*npl+2*ndset+1)*x(i) 
     &       - a3(6*npl+2*ndset+2)*x(i)**2

          dy = y_in(i) - y2(i)

	      sig2i = 1.d0/(sig(i)**2 + a3(6*npl+ndset+idset)**2)

	      loglik =  loglik + 0.5*dy*dy*sig2i +
     &               dlog(dsqrt(TWOPI*(sig(i)**2 + 
     &               a3(6*npl+ndset+idset)**2))) 
     &               - dlog(dsqrt(TWOPI)) 
c        write(*,*) loglik
        enddo
       
      return
      end      

      subroutine io_read_data(ndata,t,ts,ys,sigs,epoch,t0,t_max,
     &          ar,iar,ma,incl,cap0m,hkl)  

      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl,ma
      
   
      real*8 t(20000),y(20000),sig(20000),ys(20000),sigs(20000),PI
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      parameter(PI=3.14159265358979d0)
      real*8 ar(MMAX),incl(NPLMAX),cap0m(NPLMAX)
      integer iar(MMAX),u_off(NDSMAX),u_jit(NDSMAX),hkl
      integer idsmax(NDSMAX),ts(20000), u_incl, u_cap0m
      real*8 jitt(NDSMAX),sigscale,t0,t_max, epoch
      real*8 off(NDSMAX),loglik
      integer i,k,j,gr_flag
      character*80 infile
   
      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag

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
          ar(6*npl+i)=off(i)
          iar(6*npl+i)=u_off(i)
          ar(6*npl+ndset+i)=jitt(i)
          iar(6*npl+ndset+i)=u_jit(i)
      enddo
      
      ma = 6*npl + 2*ndset + 2
      do j = 1,npl
          i = 6*(j-1)
          read (*,*) ar(i+1),ar(i+2),ar(i+3),ar(i+4),ar(i+5),incl(j),
     &    cap0m(j),ar(i+6)
          read (*,*) iar(i+1),iar(i+2),iar(i+3),iar(i+4),iar(i+5),
     &    u_incl, u_cap0m,iar(i+6)

c         inclinations and cap0m are always ignored in the fit, just for consistency with dynamical input and output
 
      enddo

          
      read (*,*) ar(6*npl+ 2*ndset+1)
      read (*,*) iar(6*npl+2*ndset+1)   
      
      read (*,*) ar(6*npl+ 2*ndset+2)
      read (*,*) iar(6*npl+2*ndset+2)   
      
      ndata = ndata - 1
 
c      write(*,*) 'for epoch :'
      read (*,*) epoch
 
      t_max = t(ndata) 

      if (epoch.eq.0) then 
         t0 = t(1)
      else
         t0 = epoch
      endif

      read (*,*) hkl      
         
      do j = 1,npl
          i = 6*(j-1)
          if (hkl.eq.0) then 
              ar(i+4) = ar(i+4)*PI/180.d0             
          endif
          ar(i+5) = ar(i+5)*PI/180.d0               
          ar(i+6) = ar(i+6)*PI/180.d0
      enddo         
 
      do i = 1,ndata
         t(i) = (t(i) - t0)             ! time unit is day
      enddo
 
      return
      end      
 

      
      subroutine RVKEP (x,a,y,dyda,ma,ts,hkl)
      implicit none
      real*8 PI,TWOPI
      parameter (PI=3.14159265358979d0)
      parameter (TWOPI=2.0d0*PI)
      integer npl,ndset,idset,ma,i,j,NDSMAX,ts,hkl,gr_flag
      parameter (NDSMAX=20)
      integer idsmax(NDSMAX)
      real*8 x,y,a(ma),a2(ma),dyda(ma),mass(10),ap(10)
      real*8 cosw,sinw,capm,cape,cose,sine,cosf,sinf,fac1,fac2,fac3
      real*8 orbel_ehybrid, f, coswf,omega(10),capmm(10),ecc(10)
      real*8 ecc2,wm,sinwm,coswm,sin2wm,cos2wm,sin3wm,cos3wm,omegad(10)

      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag
      
      y = 0.d0
      
      do i = 1,ma
          a2(i)=a(i)
      enddo      
      
      if (hkl.eq.0) then

          do i = 1,npl
             j = 6*(i-1)
             
             if (a2(j+2).lt.0.d0) then  ! if P<0, set P>0 
                a2(j+2) = dabs(a2(j+2))
             endif         
             
             if (a2(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
                a2(j+4) = a2(j+4) + PI
                a2(j+1) = dabs(a2(j+1))
                if (a2(j+4).gt.2.d0*PI) a2(j+4) = a2(j+4)-2.d0*PI
             endif
             if (a2(j+3).lt.0.d0) then  ! if e<0, set e>0 and w=w+PI, M0=M0-PI
                a2(j+3) = dabs(a2(j+3))
                a2(j+4) = a2(j+4) +  PI
                if (a2(j+4).gt.2.d0*PI) a2(j+4) = a2(j+4)-2.d0*PI
                a2(j+5) = a2(j+5) - PI
                if (a2(j+5).lt.0.d0) a2(j+5) = a2(j+5)+2.d0*PI
             endif  
             if (a2(j+4).lt.0.d0) a2(j+4)=dmod(a2(j+4)+2.d0*PI,2.d0*PI)  
             if (a2(j+5).lt.0.d0) a2(j+5)=dmod(a2(j+5)+2.d0*PI,2.d0*PI) 
c             if (a2(j+6).lt.0.d0) a2(j+6)=dmod(a2(j+6)+2.d0*PI,2.d0*PI) 
            
             
             if (a2(j+4).gt.2.d0*PI) a2(j+4)=dmod(a2(j+4), 2.d0*PI)  
             if (a2(j+5).gt.2.d0*PI) a2(j+5)=dmod(a2(j+5), 2.d0*PI)   
c             if (a2(j+6).gt.2.d0*PI) a2(j+6)=dmod(a2(j+6), 2.d0*PI)   
            
             ecc(i) = a2(j+3) 
             omega(i) = a2(j+4) 
             capmm(i) = a2(j+5)   
c             omegad(i) = a2(j+6)      
             
             if(gr_flag.ne.0) call MA_J (a,ma,npl,1.0d0,1.0d0,
     & 	           mass,ap,hkl,gr_flag)               
        
             omegad(i) = a(j+6)             
                                  
c             write(*,*) ecc(i) ,omega(i) ,capmm(i),omegad(i) 
                                                   
          enddo  

      else   
            
          do i = 1,npl
             j = 6*(i-1)
             if (a2(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
                a2(j+4) = -1.d0*a2(j+4)       !     which is h = -h, k = -k
                a2(j+3) = -1.d0*a2(j+3)
                a2(j+1) = dabs(a2(j+1))
             endif
          
             ecc(i) = dsqrt(a2(j+3)**2 + a2(j+4)**2)
             omega(i) = datan2(a2(j+3),a2(j+4)) 
          
             if(omega(i).lt.0.d0)omega(i)=dmod(omega(i)+2.d0*PI,2.d0*PI)
             if(omega(i).gt.0.d0)omega(i)=dmod(omega(i),        2.d0*PI)
             if (a2(j+5).lt.0.d0) a2(j+5)=dmod(a2(j+5)+2.d0*PI, 2.d0*PI)
             if (a2(j+5).gt.2.d0*PI) a2(j+5) = dmod(a2(j+5), 2.d0*PI)
              
             capmm(i) = a2(j+5) - omega(i)
                
             if(capmm(i).lt.0.d0)capmm(i)=dmod(capmm(i)+2.d0*PI,2.d0*PI)
             if(capmm(i).gt.0.d0)capmm(i)=dmod(capmm(i),        2.d0*PI)
             
c             write(*,*) a2(j+4),a2(j+4),ecc(i) ,omega(i) ,capmm(i) 
          enddo        
         
      endif

      if (hkl.eq.0) then
      do j = 1,npl

          i = 6*(j-1)
          cosw = dcos(omega(j)+omegad(j)*x/365.25d0)
          sinw = dsin(omega(j)+omegad(j)*x/365.25d0)

          capm = TWOPI*x/a2(2+i) + capmm(j)
          capm = dmod(capm,  2.d0*PI )

          cape = ORBEL_EHYBRID (ecc(j),capm)
          cose = dcos(cape)
          sine = dsin(cape)
          
          cosf = (cose - ecc(j))/(1.d0 - ecc(j)*cose)
          sinf = (dsqrt(1.d0 - ecc(j)**2)*sine)/(1.d0 - ecc(j)*cose)
c          f = 2.0d0*datan2( dsqrt(1.d0 - ecc(j))*dcos(cape/2.d0),
c     &                 dsqrt(1.d0 + ecc(j))*dsin(cape/2.d0))

c          coswf = dcos(omega(j)+f)
c          fac1 = coswf + ecc(j)*cosw

          fac1 = cosw*cosf - sinw*sinf + ecc(j)*cosw

          fac2 = (cosw*sinf + sinw*cosf)/(1.d0 - ecc(j)*cose)**2
          fac3 = -a2(1+i)*dsqrt(1.d0 - ecc(j)**2)*fac2

          y = y + a2(1+i)*fac1
          dyda(1+i) = fac1
          dyda(2+i) = -TWOPI*fac3*x/a2(2+i)**2
          dyda(3+i) = -a2(1+i)*sine*(2.d0 - ecc(j)**2 - ecc(j)*cose)*
     &                 fac2/dsqrt(1.d0 - ecc(j)**2)
          dyda(4+i) = -a2(1+i)*(sinw*cosf + cosw*sinf + ecc(j)*sinw)
          dyda(5+i) = fac3
          dyda(6+i) = -a2(1+i)*(sinw*cosf + cosw*sinf 
     &                 + ecc(j)*sinw)*x/365.25d0
      enddo
      
      else
      
      do j = 1,npl

          i = 6*(j-1)
c          ecc2 = dsqrt(a2(3+i)**2 + a2(4+i)**2)

          if (ecc(j).gt.1.d-2) then
             cosw = dcos(omega(j)+omegad(i)*x/365.25d0)
             sinw = dsin(omega(j)+omegad(i)*x/365.25d0)

             capm = TWOPI*x/a2(2+i) + capmm(j)
             capm = dmod(capm,  2.d0*PI )
c             write(*,*) capm
             cape = ORBEL_EHYBRID (ecc(j),capm)
             cose = dcos(cape)
             sine = dsin(cape)
             cosf = (cose - ecc(j))/(1.d0 - ecc(j)*cose)
             sinf = (dsqrt(1.d0 - ecc(j)**2)*sine)/(1.d0 - ecc(j)*cose)

             fac1 = cosw*cosf - sinw*sinf + ecc(j)*cosw
             fac2 = cosw*sinf + sinw*cosf
             fac3 = -a2(1+i)*dsqrt(1.d0 - ecc(j)**2)*fac2/
     &              (1.d0 - ecc(j)*cose)**2
    

             y = y + a2(1+i)*fac1
             dyda(1+i) = fac1
             dyda(2+i) = -TWOPI*fac3*x/a2(2+i)**2
             dyda(3+i) = -a2(1+i)*fac2*((2.d0-ecc(j)**2-ecc(j)*cose)*
     &                   sinw*sine/dsqrt(1.d0-ecc(j)**2) -
     &                   dsqrt(1.d0 - ecc(j)*2)*cosw/ecc(j))/
     &                   (1.d0 - ecc(j)*cose)**2 -
     &                   a2(1+i)*fac2*cosw/ecc(j)
             dyda(4+i) = -a2(1+i)*fac2*((2.d0-ecc(j)**2-ecc(j)*cose)*
     &                   cosw*sine/dsqrt(1.d0 - ecc(j)**2) +
     &                   dsqrt(1.d0 - ecc(j)*2)*sinw/ecc(j))/
     &                   (1.d0 - ecc(j)*cose)**2 +
     &                   a2(1+i)*fac2*sinw/ecc(j)
             dyda(5+i) = fac3
             dyda(6+i) = dyda(4+i)*x/365.25d0
             
          else
             wm = TWOPI*x/a2(2+i) + a2(5+i)
             wm = dmod(wm,  2.d0*PI )
             
             coswm = dcos(wm)
             sinwm = dsin(wm)
             cos2wm = dcos(2.d0*wm)
             sin2wm = dsin(2.d0*wm)
             cos3wm = dcos(3.d0*wm)
             sin3wm = dsin(3.d0*wm)

             fac1 = coswm + a2(3+i)*sin2wm - a2(4+i)*(1.d0 - cos2wm) -
     &              a2(3+i)**2*(0.875d0*coswm + 1.125d0*cos3wm) -
     &              a2(3+i)*a2(4+i)*(0.25d0*sinwm - 2.25d0*sin3wm) -
     &              a2(4+i)**2*1.125d0*(coswm - cos3wm)
             fac3 = -sinwm + 2.d0*a2(3*i)*cos2wm - 2.d0*a2(4+i)*sin2wm +
     &              a2(3+i)**2*(0.875d0*sinwm - 3.375d0*sin3wm) -
     &              a2(3+i)*a2(4+i)*(0.25d0*coswm - 6.75d0*cos3wm) +
     &              a2(4+i)**2*(1.125d0*coswm - 3.375*sin3wm)

             y = y + a2(1+i)*fac1
             dyda(1+i) = fac1
             dyda(2+i) = -a2(1+i)*TWOPI*fac3*x/a2(2+i)**2
             dyda(3+i) = a2(1+i)*(sin2wm -
     &                 a2(3+i)*(1.75d0*coswm + 2.25d0*cos3wm) -
     &                 a2(4+i)*(0.25d0*sinwm - 2.25d0*sin3wm))
             dyda(4+i) = a2(1+i)*(-1.d0 + a2(4+i)*cos2wm -
     &              a2(3+i)**(0.25d0*sinwm - 2.25d0*sin3wm) -
     &              a2(4+i)*2.25d0*(coswm - cos3wm))
             dyda(5+i) = a2(1+i)*fac3
             dyda(6+i) = dyda(4+i)*x/365.25d0
          endif

      enddo      
      
      endif
      

c      do i = 1,idset
      y = y + a2(6*npl+ts)
      dyda(6*npl+ts) = 1.d0
c      enddo

c      do i = 1,idset
c          y = y + a(5*npl+i)
c          dyda(5*npl+i) = 1.d0
c      enddo


      y = y + a2(6*npl +2*ndset + 1)*x 
      y = y + a2(6*npl +2*ndset + 2) + a2(6*npl+2*ndset + 2)*x**2
c      write(*,*) a2(6*npl +2*ndset + 1), a2(6*npl+ndset + 2)
c      dyda(6*npl + ndset + 1) = x
c      dyda(6*npl + ndset + 2) = x**2         
      
   
      do i = ts+1,ndset
          dyda(6*npl+i) = 0.d0
      enddo

      return
      end
 

      subroutine prepare_for_amoeba(p,mp,np,y,a,ia,ma,mfit,funk,ndata,
     & x,z,ymod,dyda,ts,sig, it,hkl)
      integer MMAX,NDSMAX,ma,ts(20000), ndata,mp,np,mfit,it
      parameter(MMAX=200,NDSMAX=20)
      REAL*8 ftol,p(mp,np),y(mp),a(MMAX), a2(mfit),fr,frjitt
      real*8 x(20000),z(20000),ymod(20000)
      real*8 dyda(MMAX), sig(20000), loglik
      parameter(fr=0.05, frjitt=0.05)
      INTEGER i,j,k, ia(MMAX), idsmax(NDSMAX),hkl,gr_flag
      external funk
    
      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag

      k=0
      do j=1,ma
          if(ia(j).ne.0) then
          k=k+1
          p(1,k)=a(j)
          do i=2,mfit+1
              if (k.eq.(i-1)) then
                  if (j.gt.(6*npl+ndset)) then
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
     &              a,ia,hkl)
          y(i)=loglik
c          write(*,*) a2(1),a2(2),a2(3),a2(4),a2(5),a2(6),a2(7) 
c          write(*,*) a2(8),a2(9),a2(10),a2(11),a2(12),a2(13),a2(14)
c          write(*,*) a2(15),a2(16),a2(17),a2(18),a2(19),a2(20),a2(21)          
      enddo
      return
      end
          
          
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,ytry,hkl)
      implicit none
      INTEGER iter,mp,ndim,np,NMAX,ITMAX, MMAX,ma,ts(20000), ndata
      REAL*8 ftol,p(mp,np),y(mp),x(20000),z(20000),ymod(20000)
      PARAMETER (NMAX=20,ITMAX=200000,MMAX=200)
      real*8 dyda(MMAX), sig(20000), loglik, a(MMAX)
      EXTERNAL funk
      INTEGER i,ihi,ilo,inhi,j,m,n, ia(MMAX),hkl
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
     & dyda,ma,ts,sig,a,ia,hkl)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,hkl)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,hkl)
 
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              call funk(ndata,x,z,psum,ymod,dyda,ma,ndim,ts,sig,loglik,
     &                 i,a,ia,hkl)
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
     & dyda,ma,ts,sig,a,ia,hkl)
      implicit none
      INTEGER ihi,mp,ndim,np,NMAX, MMAX, ma, ts(20000),ndata
      PARAMETER (NMAX=20, MMAX=200)
      REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),x(20000),z(20000),
     & ymod(20000)
      real*8 dyda(MMAX), sig(20000),loglik
      EXTERNAL funk
      INTEGER j, ia(MMAX),hkl
      REAL*8 fac1,fac2,ytry,ptry(ndim), a(MMAX)
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      call funk(ndata,x,z,ptry,ymod,dyda,ma,ndim,ts,sig,loglik,ihi,
     & a,ia,hkl)

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




	    subroutine MA_J (a,ma,npl,m0,sini,mass,ap,hkl,gr_flag)
        
	    implicit none
	    real*8 m0,PI,TWOPI,THIRD,GMSUN,dm,MSUN
        integer npl,ma,i,j,NPLMAX,hkl,gr_flag
        parameter (NPLMAX=7)
        real*8 sini,mm(NPLMAX),ecc,corr
        real*8 a(ma),mass(NPLMAX),ap(NPLMAX),mpold(NPLMAX),mtotal
	    parameter (THIRD=1.d0/3.d0)
        parameter (PI=3.14159265358979d0,TWOPI=2.d0*PI)
	    parameter (GMSUN=1.32712497d20,MSUN=1.32712497d20)

c*******G is set to be unit, and s, m, kg as unit of time, length and mass
c*******expectively.        
        
        do j = 1,npl
           i = 6*(j-1) 

           mm(j) = 2.d0*PI/(a(i+2)*8.64d4)
           
        enddo 
 
        do i = 0,npl-1
        
           if (hkl.eq.0) then             
               ecc = a(6*i+3)     
           else
               ecc = dsqrt(a(6*i+3)**2+a(6*i+4)**2)    !! only for h, k
           endif

           mass(1) = m0
	       mpold(i+1) = 0.d0
 101       continue
           if (i.eq.0) then
           mtotal = m0
	       mass(i+2) = a(6*i+1)*(TWOPI/mm(i+1)*(m0 + mpold(i+1))**2/
     &               (TWOPI*GMSUN))**THIRD*
     &               dsqrt(1.d0 - ecc**2)
           else
              mtotal = m0
              do j = 0, i-1
                 mtotal = mtotal + mass(j+2)
              enddo
              mass(i+2) = a(6*i+1)*(TWOPI/mm(i+1)*(mtotal
     &                  +mpold(i+1))**2/(TWOPI*GMSUN))**THIRD*
     &                  dsqrt(1.d0 - ecc**2)
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

        if(gr_flag.ne.0) then
           do i = 1,npl
             j = 6*(i-1)
             call gr_corr(ap(i),a(j+3),mass(1),mass(i+1),i,corr,1.0d0) 
             a(j+6) = corr*365.25
           enddo
        endif      

        
	    return
	    end


        subroutine gr_corr(a,e,gmi,mass,n,corr,dt)	
        implicit none
	
        real*8  a ,e, T, PI, c, GMSUN, AU, st_mass
        integer n,NPLMAX
      	real*8 mass, corr,THIRD,gmi,dt,prec_frac
        parameter (PI=3.14159265358979d0)
        parameter (c = 299792458.0d0)
	    parameter (GMSUN=1.32712497d20, AU=1.49597892d13)
        parameter (THIRD=1.d0/3.d0)
 
               
        T = 2.0d0*PI * sqrt((a**3.0d0)/(gmi ) ) 
        
         corr = (24.0d0 * (PI**3.0d0) * (a**2.0d0)) / ( (T**2.0d0) *
     &    (c**2.0d0)*(1.0d0-e**2.0d0) )
         
        corr = corr/T 
      
        return
        end



***********************************************************************
c                    ORBEL_EHYBRID.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           m ==> mean anomaly. (real scalar)
*             Returns:
*              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: For e < 0.18 uses fast routine ESOLMD 
*	         For larger e but less than 0.8, uses EGET
*	         For e > 0.8 uses EHIE
*     REMARKS: Only EHIE brings M and E into range (0,TWOPI)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 25,1992.
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_ehybrid(e,m)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
	real*8 e,m

c...  Internals:
	real*8 orbel_esolmd,orbel_eget,orbel_ehie

c----
c...  Executable code 

	if(e .lt. 0.18d0) then
	  orbel_ehybrid = orbel_esolmd(e,m)
	else 
	  if( e .le. 0.8d0) then
	     orbel_ehybrid = orbel_eget(e,m)
	  else
	     orbel_ehybrid = orbel_ehie(e,m)
	  endif
	endif   

	return
	end     ! orbel_ehybrid

c--------------------------------------------------------------------



***********************************************************************
c                    ORBEL_EHIE.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           m ==> mean anomaly. (real scalar)
*             Returns:
*              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Use Danby's quartic for 3 iterations. 
*                Eqn. is f(x) = x - e*sin(x+M). Note  that
*	         E = x + M. First guess is very good for e near 1.
*	         Need to first get M between 0. and PI and use
*		 symmetry to return right answer if M between PI and 2PI
*     REMARKS: Modifies M so that both E and M are in range (0,TWOPI)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 25,1992.
*     REVISIONS: 
***********************************************************************

      real*8 function orbel_ehie(e,m)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
	real*8 e,m

c...  Internals:
      integer iflag,nper,niter,NMAX
      real*8 dx,x,sa,ca,esa,eca,f,fp

      parameter (NMAX = 3)

c----
c...  Executable code 

c In this section, bring M into the range (0,TWOPI) and if
c the result is greater than PI, solve for (TWOPI - M).
	iflag = 0
	nper = m/TWOPI
	m = m - nper*TWOPI
	if (m .lt. 0.d0) m = m + TWOPI

	if (m.gt.PI) then
	   m = TWOPI - m
	   iflag = 1
	endif

c Make a first guess that works well for e near 1.
	x = (6.d0*m)**(1.d0/3.d0) - m
	niter =0

c Iteration loop
	do niter =1,NMAX
	    call orbel_scget(x + m,sa,ca)
	    esa = e*sa
	    eca = e*ca
	    f = x - esa
	    fp = 1.d0 -eca
	    dx = -f/fp
	    dx = -f/(fp + 0.5d0*dx*esa)
	    dx = -f/(fp + 0.5d0*dx*(esa+0.3333333333333333d0*eca*dx))
	    x = x + dx
	enddo

	orbel_ehie = m + x

	if (iflag.eq.1) then
	  orbel_ehie = TWOPI - orbel_ehie
	  m = TWOPI - m
	endif

	return
	end         !orbel_ehie
c------------------------------------------------------------------


***********************************************************************
c	                  ORBEL_SCGET.F
***********************************************************************
*     PURPOSE:  Given an angle, efficiently compute sin and cos.
*
*        Input:
*             angle ==> angle in radians (real scalar)
*        
*        Output:
*             sx    ==>  sin(angle)  (real scalar)
*             cx    ==>  cos(angle)  (real scalar)
*
*     ALGORITHM: Obvious from the code 
*     REMARKS: The HP 700 series won't return correct answers for sin
*       and cos if the angle is bigger than 3e7. We first reduce it
*       to the range [0,2pi) and use the sqrt rather than cos (it's faster)
*       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 6, 1992.
*     REVISIONS: 
***********************************************************************

	subroutine orbel_scget(angle,sx,cx)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
        real*8 angle

c...  Output:
	real*8 sx,cx

c... Internals:
	integer nper
	real*8 x
	real*8 PI3BY2
	parameter(PI3BY2 = 1.5d0*PI)

c----
c...  Executable code 

        nper = angle/TWOPI
	x = angle - nper*TWOPI
	if(x.lt.0.d0) then
           x = x + TWOPI
        endif
	sx = sin(x)
	cx= sqrt(1.d0 - sx*sx)
	if( (x .gt. PIBY2) .and. (x .lt.PI3BY2)) then
           cx = -cx
        endif

	return
	end   ! orbel_scget
c-------------------------------------------------------------------




***********************************************************************
c                    ORBEL_EGET.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           m ==> mean anomaly. (real scalar)
*             Returns:
*                  orbel_eget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Quartic convergence from Danby
*     REMARKS: For results very near roundoff, give it M between
*           0 and 2*pi. One can condition M before calling EGET
*           by calling my double precision function MOD2PI(M). 
*           This is not done within the routine to speed it up
*           and because it works fine even for large M.
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 7, 1992.
*     REVISIONS: May 21, 1992.  Now have it go through EXACTLY two iterations
*                with the premise that it will only be called if
*	         we have an ellipse with e between 0.15 and 0.8
***********************************************************************

	real*8 function orbel_eget(e,m)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
	real*8 e,m

c...  Internals:
	real*8 x,sm,cm,sx,cx
	real*8 es,ec,f,fp,fpp,fppp,dx

c----
c...  Executable code 

c Function to solve Kepler's eqn for E (here called
c x) for given e and M. returns value of x.
c MAY 21 : FOR e < 0.18 use ESOLMD for speed and sufficient accuracy
c MAY 21 : FOR e > 0.8 use EHIE - this one may not converge fast enough.

	  call orbel_scget(m,sm,cm)

c  begin with a guess accurate to order ecc**3	
	  x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))

c  Go through one iteration for improved estimate
	  call orbel_scget(x,sx,cx)
	  es = e*sx
	  ec = e*cx
	  f = x - es  - m
	  fp = 1.d0 - ec 
	  fpp = es 
	  fppp = ec 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
	  orbel_eget = x + dx

c Do another iteration.
c For m between 0 and 2*pi this seems to be enough to
c get near roundoff error for eccentricities between 0 and 0.8

	  x = orbel_eget
	  call orbel_scget(x,sx,cx)
	  es = e*sx
	  ec = e*cx
	  f = x - es  - m
	  fp = 1.d0 - ec 
	  fpp = es 
	  fppp = ec 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)

	  orbel_eget = x + dx

	return
	end  ! orbel_eget
c---------------------------------------------------------------------



***********************************************************************
c                    ORBEL_ESOLMD.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           m ==> mean anomaly. (real scalar)
*             Returns:
*                orbel_esolmd ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Some sort of quartic convergence from Wisdom. 
*     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
*         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
*      	  ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI 
*     INCLUDES: needs SCGET.F
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 7, 1992.
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_esolmd(e,m)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
      real*8 e,m

c...  Internals:
	real*8 x,sm,cm,sx,cx
	real*8 es,ec,f,fp,fpp,fppp,dx

c----
c...  Executable code 

c...    Function to solve Kepler's eqn for E (here called
c...    x) for given e and M. returns value of x.

	  call orbel_scget(m,sm,cm)
	  x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))

	  call orbel_scget(x,sx,cx)
	  es = e*sx
	  ec = e*cx
	  f = x - es  - m
	  fp = 1.d0 - ec 
	  fpp = es 
	  fppp = ec 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)

	  orbel_esolmd = x + dx

	return   ! orbel_esolmd
	end
c--------------------------------------------------------------------



