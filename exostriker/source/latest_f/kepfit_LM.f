ccc   Innitually writen by Man Hoi Lee 
ccc   Heavlly modified by Trifon Trifonov trifon@hku.hk
ccc   The final version will be available in the Python RVMod lib.    


      implicit none
      real*8 PI
      parameter (PI=3.14159265358979d0)
      integer npl,ndset,idset,ndata,ma,mfit,i,j,NDSMAX,NPLMAX,MMAX
      integer writeflag_best_par
      integer writeflag_RV,writeflag_fit, amoebastarts 
      parameter (NDSMAX=20, NPLMAX=20, MMAX=200)
      integer idsmax(NDSMAX),ia(MMAX),nt, ts(20000),hkl,gr_flag
      real*8 x(20000),y(20000),sig(20000),ymod2(20000),y_in(20000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 rms,mstar,mass(NPLMAX),ap(NPLMAX)
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*8 chisq,alamda,ochisq,dchisq, epsil, deltat
      real*8 jitter,sigscale,x0,xmax,loglik,dy
      real*8 t0,t1,t2,dt,offset,t_max, incl(NPLMAX), cap0m(NPLMAX)
      real*8 st_mass,sini,m1,a1,m2,a2,jitt(NDSMAX),epoch
      real*8 ymod(20000),dyda(MMAX),best_w,best_we
      real*4 t_stop,when_to_kill,model_max,model_min
      
      external rvkep
      character*80 infile
      character*80 version_input, version

      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag


      version = "0.07"
       
      CALL getarg(1, version_input)     
      if(version_input.eq.'-version') then
          write(*,*) version
          goto 222
      endif
      
      
      
      loglik=0.d0
      
c     these two just for consistency with dynamical input and amoebastarts for consistency with loglik, not really used
      read (*,*) epsil, deltat,a moebastarts, 
     &          when_to_kill, nt, model_max,model_min,gr_flag
      
c      write(*,*) 'Stellar mass, sini'
      read (*,*) st_mass,
     &          writeflag_best_par, 
     &	             writeflag_RV,writeflag_fit 
c      write (*,*) ' Stellar Jitter: '
c      read (*,*) jitter


      call io_read_data (ndata,x,ts,y,sig,jitt,epoch,
     &               x0,t_max,a,ia,ma,incl,cap0m,hkl)
      

c Initilize 
c      write(*,*) loglik, hkl
      mfit = 0
      do j = 1,ma
          if (ia(j).ne.0) mfit = mfit + 1
      enddo
 
cccccc Hack to make it compatible with the loglik code
      if (amoebastarts.eq.0) then 
          alamda = 1.d0
          mfit = 0
          do j = 1,ma
              ia(j) = 1
              mfit = mfit + 1
          enddo
          call MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,MMAX,
     & 	           chisq,rvkep,alamda,loglik,jitt,hkl)
          goto 333
      endif




      alamda = -1.d0
      call MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,MMAX,
     & 	           chisq,rvkep,alamda,loglik,jitt,hkl)

      i = 0
c      gr_flag = 0
 500  continue
          i = i + 1
          ochisq = chisq
            
          call MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,MMAX,
     &                 chisq,rvkep,alamda,loglik,jitt,hkl)

          dchisq = chisq - ochisq
          if (i.eq.10) then
              i = 0
CC              pause
          endif

          CALL SECOND(t_stop)
          if (t_stop.ge.when_to_kill) then
            write(*,*) 'Max. time=',when_to_kill, 'sec ', 
     &                 'exceeded t_stop =', t_stop, 'sec ' 
            goto 333
           endif
          
          
      if ((chisq.ge.ochisq).or.(dchisq.lt.-1.d-2)) goto 500


333   alamda = 0.d0
      call MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,MMAX,
     &                 chisq,rvkep,alamda,loglik,jitt,hkl)

      idset = 1
      rms = 0.d0
c     y_in = y
 
         do i = 1,ndata
              idset = ts(i)
c             if (i.gt.idsmax(idset)) idset = idset + 1
              call RVKEP (x(i),a,ymod(i),dyda,ma,ts(i),hkl)
 
	      xmax = x0 + x(i)
 
 
          y_in(i) = y(i) - a(6*npl+idset) 
     &    - a(6*npl+ndset+1)*x(i)
     &    - a(6*npl+ndset+2)*x(i)**2
                   
	      ymod(i) = ymod(i) - a(6*npl+idset) 
     &    - a(6*npl+ndset+1)*x(i)
     &    - a(6*npl+ndset+2)*x(i)**2

              dy = y_in(i) - ymod(i)
              rms = rms + (y_in(i) - ymod(i))**2
              if (writeflag_RV.gt.0) then 
                  write(*,*) x0 + x(i),
     &            ymod(i), y_in(i) 
     &            + a(6*npl+ndset+1)*x(i)
     &            + a(6*npl+ndset+2)*x(i)**2,    
     &            dy, sig(i), idset
   
              endif
 
         enddo

      rms=dsqrt(rms/dble(ndata))
      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
      write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms


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
     &                 dsqrt(covar(i+5,i+5))*180.d0/PI, 0.d0, 0.d0,
     &                 dsqrt(covar(i+6,i+6))*180.d0/PI
          enddo
          write (*,*) 'Best-fit V0 [m/s] and their error bars:'
          do j = 1,ndset
              i = 6*npl + j
              write (*,*) a(i)
c              offset(j) = a(i)
              write (*,*) dsqrt(covar(i,i))
          enddo

          write (*,*) 'Jitters (not fitted in chi^2 method):'
          do j = 1,ndset
              write (*,*) jitt(j)
              write (*,*) '0'
          enddo          
          
          write (*,*) 'linear trend [m/s per day]:'
          write (*,*) a(6*npl + ndset + 1)  
          write (*,*) dsqrt(covar(6*npl + ndset + 1,6*npl + ndset + 1))
          
          write (*,*) 'quad. trend [m/s per day]:'
          write (*,*) a(6*npl + ndset + 2)  
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


c      nt = 20000
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
      

c*************************************************************************      
c**********************       read RV data      **************************
c*************************************************************************

C     Xianyu Tan 2011
                                                                            
      subroutine io_read_data(ndata,t,ts,ys,sigs,jitt,epoch,t0,t_max,
     &          ar,iar,ma,incl,cap0m,hkl)  

      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl,ma
      real*8 t(20000),y(20000),sig(20000),ys(20000),sigs(20000),PI
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      parameter(PI=3.14159265358979d0)
      real*8 ar(MMAX),incl(NPLMAX),cap0m(NPLMAX)
      integer iar(MMAX),u_off(NDSMAX),u_jit, hkl
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
          read (*,*) u_jit
          open (unit=10,file=infile)
 100      continue
             read (10,*,err=200,end=200) t0,y(ndata),sig(ndata)

c              sig(ndata) = dsqrt(sig(ndata)**2 + jitt(i)**2)
c      	      write (*,*) "TEST"
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
c          WRITE (*,*) idsmax(i)
          close (unit=10)
                    
          
          enddo
          read (*,*) npl
      if (npl.gt.NPLMAX) stop ' KEPFIT: npl > NPLMAX.'

      do i = 1,ndset
          ar(6*npl+i)=off(i)
          iar(6*npl+i)=u_off(i)
      enddo
      
      ma = 6*npl + ndset + 2
c      write (*,*) 'Initial K, P, e, w, M0 and their flags: '
      do j = 1,npl
          i = 6*(j-1)
          read (*,*) ar(i+1),ar(i+2),ar(i+3),ar(i+4),ar(i+5),incl(j),
     &    cap0m(j),ar(i+6)
          read (*,*) iar(i+1),iar(i+2),iar(i+3),iar(i+4),iar(i+5),
     &    u_incl, u_cap0m,iar(i+6)

c         inclinations and cap0m are always ignored in the fit, just for consistency with dynamical input and output
 
      enddo
c u_jit read for consistency with input/output in loglik case, but here we don't actually save this information and not use it, jitters cannot be used for fit in chi^2 minimization
          
      read (*,*) ar(6*npl+ ndset+1)
      read (*,*) iar(6*npl+ndset+1)    

      read (*,*) ar(6*npl+ ndset+2)
      read (*,*) iar(6*npl+ndset+2) 

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
      integer npl,ndset,idset,ma,i,j,NDSMAX,ts,hkl
      parameter (NDSMAX=20)
      integer idsmax(NDSMAX),gr_flag
      real*8 x,y,a(ma),dyda(ma),st_mass,mass(10),ap(10)
      real*8 cosw,sinw,capm,cape,cose,sine,cosf,sinf,fac1,fac2,fac3
      real*8 orbel_ehybrid, f, coswf,omega(10),capmm(10),ecc(10)
      real*8 ecc2,wm,sinwm,coswm,sin2wm,cos2wm,sin3wm,cos3wm,omegad(10)

      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag
      
      y = 0.d0

      
      
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
             
             
             if (a(j+4).gt.2.d0*PI) a(j+4) = dmod(a(j+4), 2.d0*PI)  
             if (a(j+5).gt.2.d0*PI) a(j+5) = dmod(a(j+5), 2.d0*PI)   
c             if (a(j+6).gt.2.d0*PI) a(j+6) = dmod(a(j+6), 2.d0*PI)   
             
             ecc(i) = a(j+3) 
             omega(i) = a(j+4) 
             capmm(i) = a(j+5)                               
c             gr_flag = 1
             if(gr_flag.ne.0) call MA_J (a,ma,npl,1.0d0,1.0d0,
     & 	           mass,ap,hkl,gr_flag)               
        
             omegad(i) = a(j+6)
                          
c             write(*,*) ecc(i) ,omega(i) ,capmm(i),omegad(i)
                                                   
          enddo  

      else   
            
          do i = 1,npl
             j = 6*(i-1)
             if (a(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
                a(j+4) = -1.d0*a(j+4)       !     which is h = -h, k = -k
                a(j+3) = -1.d0*a(j+3)
                a(j+1) = abs(a(j+1))    
             endif
          
             ecc(i) = dsqrt(a(j+3)**2 + a(j+4)**2)
             omega(i) = atan2(a(j+3),a(j+4)) 
          
             if(omega(i).lt.0.d0)omega(i)=dmod(omega(i)+2.d0*PI,2.d0*PI)  
             if(omega(i).gt.0.d0)omega(i)=dmod(omega(i),        2.d0*PI)              
             if (a(j+5).lt.0.d0) a(j+5) = dmod(a(j+5)+2.d0*PI, 2.d0*PI) 
             if (a(j+5).gt.2.d0*PI) a(j+5) = dmod(a(j+5), 2.d0*PI )   
              
             capmm(i) = a(j+5) - omega(i)        
                
             if(capmm(i).lt.0.d0)capmm(i)=dmod(capmm(i)+2.d0*PI,2.d0*PI)  
             if(capmm(i).gt.0.d0)capmm(i)=dmod(capmm(i),        2.d0*PI)                 
             
c             write(*,*) a(j+4),a(j+4),ecc(i) ,omega(i) ,capmm(i) 
          enddo        
         
      endif

      if (hkl.eq.0) then
      do j = 1,npl

          i = 6*(j-1)
c          cosw = dcos(omega(j))
c          sinw = dsin(omega(j))
          cosw = dcos(omega(j)+omegad(j)*x/365.25d0)
          sinw = dsin(omega(j)+omegad(j)*x/365.25d0)
c          write(*,*) omegad(j)
          capm = TWOPI*x/a(2+i) + capmm(j)
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
          fac3 = -a(1+i)*dsqrt(1.d0 - ecc(j)**2)*fac2

          y = y + a(1+i)*fac1
          dyda(1+i) = fac1
          dyda(2+i) = -TWOPI*fac3*x/a(2+i)**2
          dyda(3+i) = -a(1+i)*sine*(2.d0 - ecc(j)**2 - ecc(j)*cose)*
     &                 fac2/dsqrt(1.d0 - ecc(j)**2)
          dyda(4+i) = -a(1+i)*(sinw*cosf + cosw*sinf + ecc(j)*sinw)
          dyda(5+i) = fac3
          dyda(6+i) = -a(1+i)*(sinw*cosf + cosw*sinf
     &                 + ecc(j)*sinw)*x/365.25d0
          
      enddo
      
      else
      
      do j = 1,npl

          i = 6*(j-1)
c          ecc2 = dsqrt(a(3+i)**2 + a(4+i)**2)

          if (ecc(j).gt.1.d-2) then
             cosw = dcos(omega(j))
             sinw = dsin(omega(j))

             capm = TWOPI*x/a(2+i) + capmm(j)
             capm = dmod(capm,  2.d0*PI )
c             write(*,*) capm
             cape = ORBEL_EHYBRID (ecc(j),capm)
             cose = dcos(cape)
             sine = dsin(cape)
             cosf = (cose - ecc(j))/(1.d0 - ecc(j)*cose)
             sinf = (dsqrt(1.d0 - ecc(j)**2)*sine)/(1.d0 - ecc(j)*cose)

             fac1 = cosw*cosf - sinw*sinf + ecc(j)*cosw
             fac2 = cosw*sinf + sinw*cosf
             fac3 = -a(1+i)*dsqrt(1.d0 - ecc(j)**2)*fac2/
     &              (1.d0 - ecc(j)*cose)**2
    

             y = y + a(1+i)*fac1
             dyda(1+i) = fac1
             dyda(2+i) = -TWOPI*fac3*x/a(2+i)**2
             dyda(3+i) = -a(1+i)*fac2*((2.d0 - ecc(j)**2 - ecc(j)*cose)*
     &                   sinw*sine/dsqrt(1.d0 - ecc(j)**2) -
     &                   dsqrt(1.d0 - ecc(j)*2)*cosw/ecc(j))/
     &                   (1.d0 - ecc(j)*cose)**2 -
     &                   a(1+i)*fac2*cosw/ecc(j)
             dyda(4+i) = -a(1+i)*fac2*((2.d0 - ecc(j)**2 - ecc(j)*cose)*
     &                   cosw*sine/dsqrt(1.d0 - ecc(j)**2) +
     &                   dsqrt(1.d0 - ecc(j)*2)*sinw/ecc(j))/
     &                   (1.d0 - ecc(j)*cose)**2 +
     &                   a(1+i)*fac2*sinw/ecc(j)
             dyda(5+i) = fac3
             dyda(6+i) = dyda(4+i)*x/365.25d0             
          else
             wm = TWOPI*x/a(2+i) + a(5+i)
             wm = dmod(wm,  2.d0*PI )
             
             coswm = dcos(wm)
             sinwm = dsin(wm)
             cos2wm = dcos(2.d0*wm)
             sin2wm = dsin(2.d0*wm)
             cos3wm = dcos(3.d0*wm)
             sin3wm = dsin(3.d0*wm)

             fac1 = coswm + a(3+i)*sin2wm - a(4+i)*(1.d0 - cos2wm) -
     &              a(3+i)**2*(0.875d0*coswm + 1.125d0*cos3wm) -
     &              a(3+i)*a(4+i)*(0.25d0*sinwm - 2.25d0*sin3wm) -
     &              a(4+i)**2*1.125d0*(coswm - cos3wm)
             fac3 = -sinwm + 2.d0*a(3*i)*cos2wm - 2.d0*a(4+i)*sin2wm +
     &              a(3+i)**2*(0.875d0*sinwm - 3.375d0*sin3wm) -
     &              a(3+i)*a(4+i)*(0.25d0*coswm - 6.75d0*cos3wm) +
     &              a(4+i)**2*(1.125d0*coswm - 3.375*sin3wm)

             y = y + a(1+i)*fac1
             dyda(1+i) = fac1
             dyda(2+i) = -a(1+i)*TWOPI*fac3*x/a(2+i)**2
             dyda(3+i) = a(1+i)*(sin2wm -
     &                 a(3+i)*(1.75d0*coswm + 2.25d0*cos3wm) -
     &                 a(4+i)*(0.25d0*sinwm - 2.25d0*sin3wm))
             dyda(4+i) = a(1+i)*(-1.d0 + a(4+i)*cos2wm -
     &              a(3+i)**(0.25d0*sinwm - 2.25d0*sin3wm) -
     &              a(4+i)*2.25d0*(coswm - cos3wm))
             dyda(5+i) = a(1+i)*fac3
             dyda(6+i) = dyda(4+i)*x/365.25d0             
             
          endif

      enddo      
      
      endif
      

c      do i = 1,idset
      y = y + a(6*npl+ts)
c      write(*,*) ts
c      do i = 1,ts-1
c          dyda(6*npl+i) = 0.d0
c      enddo
      
      dyda(6*npl+ts) = 1.d0

      y = y + a(6*npl +ndset + 1)*x + a(6*npl + ndset + 2)*x**2
        
      dyda(6*npl + ndset + 1) = x
      dyda(6*npl + ndset + 2) = x**2
   
      do i = ts+1,ndset
          dyda(6*npl+i) = 0.d0
      enddo

      return
      end
      
      
      

c MRQMIN attempts to reduce the chi^2 of a fit by the Levenberg-Marquardt
c method. It uses COVSRT, GAUSSJ, and MRQCOF.
c
c From Numerical Recipes.

	subroutine MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,nca,
     & 	                   chisq,funcs,alamda,loglik,jitt,hkl)
	implicit none
	integer ma,nca,ndata,ia(ma),MMAX,NDSMAX,ts(20000)
	real*8 alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     &	       sig(ndata),x(ndata),y(ndata),loglik
	external funcs
	parameter (MMAX=200,NDSMAX=20)
	integer j,k,l,mfit,hkl
	real*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX),jitt(NDSMAX)
	save ochisq,atry,beta,da,mfit

c Initialization.
	if (alamda.lt.0.d0) then
	  mfit = 0
	  do j = 1,ma
	    if (ia(j).ne.0) mfit = mfit + 1
	  enddo
	  alamda = 0.001d0
CC	  alamda = 1000.d0
	  call MRQCOF (x,y,sig,ndata,a,ia,ma,ts,
     &	       alpha,beta,nca,chisq,funcs,loglik,jitt,hkl)

	  ochisq = chisq
	  do j = 1,ma
	    atry(j) = a(j)
	  enddo
	endif
	       

c Alter linearized fitting matrix by augmenting diagonal elements.
	do j = 1,mfit
	  do k = 1,mfit
	    covar(j,k) = alpha(j,k)
	  enddo
	  covar(j,j) = alpha(j,j)*(1.d0 + alamda)
	  da(j) = beta(j)
	enddo

c Matrix solution.
	call GAUSSJ (covar,mfit,nca,da,1,1)

c Evaluate covariance matrix once converged.
	if (alamda.eq.0.d0) then
	  call COVSRT (covar,nca,ma,ia,mfit)
	  call COVSRT (alpha,nca,ma,ia,mfit)
	  return
	endif

	j = 0
	do l = 1,ma
	  if (ia(l).ne.0) then
	    j = j + 1
	    atry(l) = a(l) + da(j)
	  endif
	enddo

	call MRQCOF (x,y,sig,ndata,atry,ia,ma,
     &	       ts,covar,da,nca,chisq,funcs,loglik,jitt,hkl)

	if (chisq.lt.ochisq) then
c	  Accept new solution.
	  alamda = 0.1d0*alamda
	  ochisq = chisq
	  do j = 1,mfit
	    do k = 1,mfit
	      alpha(j,k) = covar(j,k)
	    enddo
	    beta(j) = da(j)
	  enddo
	  do l = 1,ma
	    a(l) = atry(l)
	  enddo
	else
c	  Increase alamda and return.
	  alamda = 10.d0*alamda
	  chisq = ochisq
	endif

	return
	end

c MRQCOF evaluates the linearized fitting matrix alpha, the vector beta,
c and chi^2.
c
c From Numerical Recipes.

	subroutine MRQCOF (x,y,sig,ndata,a,ia,ma,ts,alpha,beta,nalp,
     &	                   chisq,funcs,loglik,jitt,hkl)

	implicit none
	integer npl,ndset,idset,ma,nalp,ndata,ia(ma),NDSMAX,MMAX
	parameter (NDSMAX=20, MMAX=200)
        integer idsmax(NDSMAX),ts(20000),hkl,gr_flag
	real*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),
     &	       x(ndata),y(ndata),loglik,TWOPI,jitt(NDSMAX)
        parameter (TWOPI=2.d0*3.14159265358979d0)
	external funcs
	integer mfit,i,j,k,l,m
	real*8 dy,sig2i,wt,ymod,dyda(MMAX)

      common /DSBLK/ npl,ndset,idsmax,idset,gr_flag

	mfit = 0
	do j = 1,ma
	  if (ia(j).ne.0) mfit = mfit + 1
	enddo

c Initialize (symmetric) alpha and beta.
	do j = 1,mfit
	  do k = 1,j
	    alpha(j,k) = 0.d0
	  enddo
	  beta(j) = 0.d0
	enddo
	chisq = 0.d0
	loglik=0.d0

c Loop over all data.
        idset = 1
	do i = 1,ndata
cc          if (i.gt.idsmax(idset)) idset = idset + 1
          idset = ts(i)
 
	  call FUNCS (x(i),a,ymod,dyda,ma,idset,hkl)
c	  sig2i = 1.d0/(sig(i)*sig(i))
c          write(*,*) sig(i),ts(i)
          sig2i = 1.d0/(sig(i)**2 + jitt(idset)**2)
	  dy = y(i) - ymod
	  j = 0


	  do l = 1,ma
	    if (ia(l).ne.0) then
	      j = j + 1
	      wt = dyda(l)*sig2i
c              write(*,*) wt
	      k = 0
	      do m = 1,l
	        if (ia(m).ne.0) then
	          k = k + 1
	          alpha(j,k) = alpha(j,k) + wt*dyda(m)
	        endif
	      enddo
	      beta(j) = beta(j) + dy*wt
	    endif
	  enddo

c	  Compute chi^2 and loglik.
	  chisq = chisq + dy*dy*sig2i
	  loglik =  loglik - 0.5*dy*dy*sig2i -
     &               0.5*dlog(TWOPI*(sig(i)**2 + 
     &               jitt(idset)**2))
	  
	  enddo


c Fill in the symmetric side.
	do j = 2,mfit
	  do k = 1,j-1
	    alpha(k,j) = alpha(j,k)
	  enddo
	enddo

	return
	end

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

c GAUSSJ solves linear equation by Gauss-Jordan elimination.
c
c From Numerical Recipes.

	subroutine GAUSSJ (a,n,np,b,m,mp)

	implicit none
	integer m,mp,n,np,NMAX
	real*8 a(np,np),b(np,mp)
	parameter (NMAX=50)
	integer i,icol,irow,j,k,l,ll,kkk,
     &  indxc(NMAX),indxr(NMAX),ipiv(NMAX)
	real*8 big,dum,pivinv

	do j = 1,n
	  ipiv(j) = 0
	enddo

c Main loop over the columns to be reduced.
	do i = 1,n

c	  Search for a pivot element.
	  big = 0.d0
	  do j = 1,n
	    if (ipiv(j).ne.1) then
	      do k = 1,n
	        if (ipiv(k).eq.0) then
	          if (dabs(a(j,k)).ge.big) then
	            big = dabs(a(j,k))
	            irow = j
	            icol = k
	          endif
	        endif
	      enddo
	    endif
	  enddo
	  ipiv(icol) = ipiv(icol) + 1

c	  Interchange rows to put the pivot element on the diagonal.
	  if (irow.ne.icol) then
	    do l = 1,n
	      dum = a(irow,l)
	      a(irow,l) = a(icol,l)
	      a(icol,l) = dum
	    enddo
	  endif
	  indxr(i) = irow
	  indxc(i) = icol

c	  if (a(icol,icol).eq.0.d0) pause ' singular matrix in GAUSSJ'

c	  Divide pivot row by the pivot element.
	  pivinv = 1.d0/a(icol,icol)
	  a(icol,icol) = 1.d0
	  do l = 1,n
	    a(icol,l) = a(icol,l)*pivinv
	  enddo
	  do l = 1,m
	    b(icol,l) = b(icol,l)*pivinv
	  enddo

c	  Reduce the rows, except for the pivot one.
	  do ll = 1,n
	    if (ll.ne.icol) then
	      dum = a(ll,icol)
	      a(ll,icol) = 0.d0
	      do l = 1,n
	        a(ll,l) = a(ll,l) - a(icol,l)*dum
	      enddo
	      do l = 1,m
	        b(ll,l) = b(ll,l) - b(icol,l)*dum
	      enddo
	    endif
	  enddo

	enddo

c Unscramble the solution in view of the column interchanges.
	do l = n,1,-1
	  if (indxr(l).ne.indxc(l)) then
	    do k = 1,n
	      dum = a(k,indxr(l))
	      a(k,indxr(l)) = a(k,indxc(l))
	      a(k,indxc(l)) = dum
	    enddo
	    do kkk=1,m
	      dum = b(indxr(l),kkk)
	      b(indxr(l),kkk) = b(indxc(l),kkk)
	      b(indxc(l),kkk) = dum
	    enddo
	  endif
	enddo

	return
	end

c COVSRT expands in storage the covariance matrix covar, so as to take
c into account parameters that are being held fixed. (For the latter,
c return zero covariances.)
c
c From Numerical Recipes.

	subroutine COVSRT (covar,npc,ma,ia,mfit)

	implicit none
	integer ma,mfit,npc,ia(ma)
	real*8 covar(npc,npc)
	integer i,j,k
	real*8 swap

	do i = mfit+1,ma
	  do j = 1,i
	    covar(i,j) = 0.d0
	    covar(j,i) = 0.d0
	  enddo
	enddo

	k = mfit
	do j = ma,1,-1
	  if (ia(j).ne.0) then
	    do i = 1,ma
	      swap = covar(i,k)
	      covar(i,k) = covar(i,j)
	      covar(i,j) = swap
	    enddo
	    do i = 1,ma
	      swap = covar(k,i)
	      covar(k,i) = covar(j,i)
	      covar(j,i) = swap
	    enddo
	    k = k - 1
	  endif
	enddo

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

