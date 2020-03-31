ccc   Innitually writen by Man Hoi Lee, Xianyu Tan, 
ccc   Heavlly modified by Trifon Trifonov trifon@hku.hk
ccc   The final version will be available in the Python RVMod lib.    


      implicit none
      real*8 PI
      parameter (PI=3.14159265358979d0)
      integer npl,ndset,idset,ndata,ma,mfit,i,j,NDSMAX,NPLMAX,MMAX,k
      real*8 inc,sini,mstar
      integer writeflag_best_par, nt,hkl,gr_flag
      integer writeflag_RV,writeflag_fit, amoebastarts
      parameter (NDSMAX=20, NPLMAX=20, MMAX=200)
      integer idsmax(NDSMAX),ia(MMAX),ts(20000)
      real*8 t(20000),y(20000),sig(20000),ys(20000),sigs(20000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 chisq,alamda,ochisq,dchisq,loglik,epsil,deltat
      real*8 jitter(NDSMAX),sigscale,t0,t_max,epoch
      real*8 rms,ymod(20000),dyda(20000,MMAX)
      real*8 a0(MMAX),asave(MMAX),ap(NPLMAX)
      real*8 wdot(NPLMAX),u_wdot(NPLMAX)
 
      real*8 tstop,dt,dtout
      real*4 t_stop,when_to_kill,model_max,model_min
       
      external rvkep_ewcop_fin
      character*80 infile
      character*80 version_input, version

      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      version = "0.07"
       
      CALL getarg(1, version_input)     
      if(version_input.eq.'-version') then
          write(*,*) version
          goto 222
      endif

      loglik=0.d0
      rms=0.d0
     
c amoebastarts for consistency with loglik input
     
      read (*,*) epsil,deltat, amoebastarts,
     &          when_to_kill, nt, model_max,model_min,gr_flag
     
c      write(*,*) 'Stellar mass: '
      read (*,*) mstar,
     &          writeflag_best_par, 
     &	             writeflag_RV,writeflag_fit 
   

      call io_read_data (ndata,t,ts,ys,sigs,jitter,
     & 	           epoch,t0,t_max,a,ia,ma,mfit,hkl,wdot,u_wdot)


cccccc Hack to make it compatible with the loglik code
c      if (amoebastarts.eq.0) then 
c          alamda = 1.d0
c          mfit = 0
c          do j = 1,ma
c              ia(j) = 1
c              mfit = mfit + 1
c          enddo
c         call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
c     & chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat)
c          goto 333
c      endif


c*****set alamda to be negtive for initializing******
      alamda = -1.d0
      call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
     & chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat,hkl)
      
c     write(*,*) ' alamda,chi_nu^2: ',alamda,chisq/dble(ndata-mfit)
 
      i = 0
 500  continue
          i = i + 1
          ochisq = chisq
          call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
     &  chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat,hkl)
 
 
          dchisq = chisq - ochisq
          
c          write(*,*) chisq

          if ((i.eq.200).or.(alamda.ge.1d6)) then
              i = 0
              goto 502
          endif
          
          CALL SECOND(t_stop)
          if (t_stop.ge.when_to_kill) then
               write(*,*) 'Max. time=',when_to_kill, 'sec ', 
     &                 'exceeded t_stop =', t_stop, 'sec ' 
               goto 502
          endif        

      if ((chisq.ge.ochisq).or.(dchisq.lt.-1d-3)) goto 500

 502  alamda = -1.d0
      i = 0
 501  continue

          i = i + 1
          ochisq = chisq
          call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
     &  chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat,hkl)

 
c          if (alamda.gt.1d12) alamda = -1.d0
          dchisq = chisq - ochisq

          if ((i.eq.200).or.(alamda.ge.1d6)) then
              i = 0
              alamda = -1.d0
              goto 333
          endif
      if ((chisq.ge.ochisq).or.(dchisq.lt.-1d-5)) then
      goto 501
      endif


c*******final output******************
333   alamda = 0.d0

      call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
     & chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat,hkl)

     
      call io_write_bestfitpa_ewcop_fin (a,covar,t,ys,ndata,ts,
     & 	           ma,mfit,t0,t_max,sigs,chisq,rms,writeflag_RV,
     &  writeflag_best_par,writeflag_fit,jitter,epsil,deltat,
     &  nt, model_max,model_min,hkl,wdot,u_wdot)

c      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
c      write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms
 
222   end



      subroutine io_read_data (ndata,t,ts,ys,sigs,jitter,epoch,
     &               t0,t_max,ar,iar,ma,mfit,hkl,wdot,u_wdot) 


      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl
      real*8 t(20000),y(20000),sig(20000),ys(20000),sigs(20000)
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      integer idsmax(NDSMAX),ts(20000),hkl
      real*8 jitter(NDSMAX),t0,t_max,epoch,ar(MMAX),off(NDSMAX), PI
      parameter(PI=3.14159265358979d0)
      integer i,k,j, iar(MMAX), u_off(NDSMAX), u_jit, ma, mfit
      character*80 infile
      real*8 wdot(NPLMAX),u_wdot(NPLMAX)
   
      common /DSBLK/ npl,ndset,idsmax,idset

c      write(*,*)  ' Number of Data Sets: '
      read (*,*) ndset
      if (ndset.gt.NDSMAX) stop ' KEPFIT: ndset > NDSMAX.'

      ndata = 1
      do i = 1,ndset
          read (*,50) infile
 50       format (a)
          read (*,*) off(i)
          read (*,*) u_off(i)
          read (*,*) jitter(i)
          read (*,*) u_jit
          open (unit=10,file=infile)
 100      continue
              read (10,*,err=200,end=200) t0,y(ndata),sig(ndata)

c              sig(ndata) = dsqrt(sig(ndata)**2 + jitter(i)**2)
 
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

      ndata = ndata - 1
      
      read (*,*) npl
      
      do i = 1,ndset
          ar(7*npl+i)=off(i)
          iar(7*npl+i)=u_off(i)
      enddo
      
c           write(*,*) ' Number of Planets: '
      if (npl.gt.NPLMAX) stop ' KEPFIT: npl > NPLMAX.'

      ma = 7*npl + ndset + 2
      
c      write(*,*) 'Initial K, P, e, w, M0,Inc,Capom and their flags: '
      do j = 1,npl
          i = 7*(j-1)
          read (*,*) ar(i+1),ar(i+2),ar(i+3),ar(i+4),ar(i+5),ar(i+6),
     &               ar(i+7),wdot(i)
          read (*,*) iar(i+1),iar(i+2),iar(i+3),iar(i+4),iar(i+5),
     &               iar(i+6),iar(i+7),u_wdot(i)
 
 
c          a(i+2) = 2.d0*PI/(a(i+2)*8.64d4)         ! mean motion 
          ar(i+2) = ar(i+2)*8.64d4         ! second as unit
          u_wdot(i) = 0 
      enddo


      
c      write (*,*) 'linear trend:'      
      read (*,*) ar(7*npl + ndset + 1)
      read (*,*) iar(7*npl + ndset + 1)   

      read (*,*) ar(7*npl + ndset + 2)
      read (*,*) iar(7*npl + ndset + 2)        

      read (*,*) epoch
 
      t_max = t(ndata) 

      if (epoch.eq.0) then 
         t0 = t(1)
      else
         t0 = epoch
      endif
      
      read (*,*) hkl      
         
      do j = 1,npl
          i = 7*(j-1)
          if (hkl.eq.0) then 
              ar(i+4) = ar(i+4)*PI/180.d0
          endif
              
          ar(i+5) = ar(i+5)*PI/180.d0
          ar(i+6) = ar(i+6)*PI/180.d0
          ar(i+7) = ar(i+7)*PI/180.d0          
      enddo          
      
      
      
      do i = 1,ndata
         t(i) = (t(i) - t0)*8.64d4               ! time unit is second
      enddo
      
      mfit = 0
      do j = 1,ma
          if (iar(j).ne.0) mfit = mfit + 1
      enddo
      
      return
      end
      

C**************************************************************************
C**********   output best-fit parameters and errorbars    *****************
C**************************************************************************

      subroutine io_write_bestfitpa_ewcop_fin (a,covar,t,ys,ndata,ts,
     &           ma,mfit,t0,t_max,sigs,chisq,rms,writeflag_RV,
     &           writeflag_best_par,writeflag_fit,jitter,epsil,
     &           deltat,nt, model_max,model_min,hkl,wdot,u_wdot)
   
      implicit none 
      real*8 PI
      integer MMAX,NDSMAX,NPLMAX 
      parameter (PI=3.14159265358979d0,MMAX=200  ,NDSMAX=20,NPLMAX=20)
      real*8 a(MMAX),ia(MMAX),t(20000),ymod(20000),ys(20000)
      real*8 covar(MMAX,MMAX),dyda(20000,MMAX),AU,day,loglik,dy,twopi
      real*8 rms,mstar,sini,mass(NPLMAX),ap(NPLMAX)
      integer ts(20000),nbod,nt,writeflag_RV,
     &           writeflag_best_par,writeflag_fit
      real*8 t0,t1,t2,dt,offset,t_max,chisq,deltat,epsil,sig2i 
      real*8 x(20000),y(20000),sigs(20000),jitter(NDSMAX)
      integer i,j,npl,ndset,ndata,idset,mfit,ma,idsmax(NDSMAX),hkl
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),vxh(NPLMAX),vyh(NPLMAX)
     &       ,vzh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),vxj(NPLMAX),vyj(NPLMAX)
     &       ,vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX)
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*4 model_max,model_min
      parameter (AU=1.49597892d11, day = 86400.d0)
      real*8 wdot(NPLMAX),u_wdot(NPLMAX),best_w,best_we


      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini
            
      nbod = npl+1
      rms = 0.d0
      twopi = 2.0d0*PI  
      chisq = 0
ccccccccccccccccccc t[JD], obs., cal., O-C   ccccccccccccc   

      call RVKEP_ewcop_fin (t,a,ymod,dyda,ma,ndata,ia,epsil,deltat,hkl)
      call MA_J_cop_fin (a,ma,npl,mstar,sini,mass,ap,hkl)





      do i = 1,npl
         j = 7*(i-1)
         
   
          if (hkl.eq.0) then
              
c             a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)
             
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
             if (a(j+4).lt.0.d0) a(j+4)=dmod(a(j+4)+2.d0*PI,2.d0*PI)  
             if (a(j+5).lt.0.d0) a(j+5)=dmod(a(j+5)+2.d0*PI,2.d0*PI) 
             if (a(j+4).gt.2.d0*PI) a(j+4)=dmod(a(j+4),2.d0*PI )  
             if (a(j+5).gt.2.d0*PI) a(j+5)=dmod(a(j+5),2.d0*PI )         
             if (a(j+6).lt.0.d0) a(j+6)=dmod(a(j+6)+2.d0*PI,2.d0*PI)  
             if (a(j+7).lt.0.d0) a(j+7)=dmod(a(j+7)+2.d0*PI,2.d0*PI) 
             if (a(j+6).gt.2.d0*PI) a(j+6)=dmod(a(j+6),2.d0*PI)  
             if (a(j+7).gt.2.d0*PI) a(j+7)=dmod(a(j+7),2.d0*PI)
     
          else
     
             if (a(j+2).lt.0.d0) then  ! if P<0, set P>0 
                a(j+2) = abs(a(j+2))
             endif                   
             
             if (a(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
                a(j+4) = -1.d0*a(j+4)       !     which is h = -h, k = -k
                a(j+3) = -1.d0*a(j+3)
                a(j+1) = abs(a(j+1))    
             endif

             if (a(j+5).lt.0.d0) a(j+5)=dmod(a(j+5)+2.d0*PI,2.d0*PI) 
             if (a(j+6).lt.0.d0) a(j+6)=dmod(a(j+6)+2.d0*PI,2.d0*PI)  
             if (a(j+7).lt.0.d0) a(j+7)=dmod(a(j+7)+2.d0*PI,2.d0*PI) 
             if (a(j+5).gt.2.d0*PI) a(j+5)=dmod(a(j+5),2.d0*PI)
             if (a(j+6).gt.2.d0*PI) a(j+6)=dmod(a(j+6),2.d0*PI)  
             if (a(j+7).gt.2.d0*PI) a(j+7)=dmod(a(j+7),2.d0*PI)    
     
          endif
      enddo



      do i = 1,ndata
          idset = ts(i)
 
	      ys(i) = ys(i) - a(7*npl +idset)
     &                -a(7*npl +ndset + 1)*(t(i)/86400.d0)
     &                -a(7*npl +ndset + 2)*(t(i)/86400.d0)**2
 
          if (writeflag_RV.gt.0) then
          write(*,*) t0 + t(i)/8.64d4 ,ymod(i),ys(i)
     &                +a(7*npl +ndset + 1)*(t(i)/86400.d0)
     &                +a(7*npl +ndset + 2)*(t(i)/86400.d0)**2,     
     &                ys(i) - ymod(i),sigs(i),ts(i)

          endif
     
          sig2i = 1.d0/(sigs(i)**2 + jitter(idset)**2)

          dy =  ys(i) -ymod(i)  
 

 	      chisq  = chisq + dy*dy*sig2i
 
c          write(*,*) "TEST:",loglik,dy,ymod(i),jitter(idset)
	      loglik =  loglik - 0.5*dy*dy*sig2i -
     &               0.5*dlog(twopi*(sigs(i)**2
     &                + jitter(idset)**2)) 
     
 

          rms = rms + dy**2     

c          rms = rms + (ys(i) - ymod(i))**2
      enddo
      rms = dsqrt(rms/dble(ndata))




      if(writeflag_best_par.gt.0) then
          write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
          write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms

          write (*,*) 'Best-fit K [m/s], P [days], e, w [deg], 
     & M0 [deg], i[deg], cap0m[deg], w dot [deg/yr], and their errors'
          do j = 1,npl
              i = 7*(j-1)

              if (hkl.eq.0) then
                  best_w = a(i+4)*180.d0/PI
                  best_we = dsqrt(covar(i+4,i+4))*180.d0/PI
              else    
                  best_w = a(i+4) 
                  best_we = dsqrt(covar(i+4,i+4)) 
              endif


              write(*,*) a(i+1),a(i+2)/8.64d4,a(i+3),
     &                best_w,a(i+5)*180.d0/PI,
     &                dmod(a(i+6)*180.d0/PI,180.d0),
     &                dmod(a(i+7)*180.d0/PI,360.d0),wdot(i)
              write(*,*) dsqrt(covar(i+1,i+1)),
     &                dsqrt(covar(i+2,i+2))/8.64d4,
c     &                2.d0*PI/a(i+2)**2*dsqrt(covar(i+2,i+2))/8.64d4,
     &                dsqrt(covar(i+3,i+3)),best_we,
     &                dsqrt(covar(i+5,i+5))*180.d0/PI,
     &                dmod(dsqrt(covar(i+6,i+6))*180.d0/PI,180.d0),
     &                dmod(dsqrt(covar(i+7,i+7))*180.d0/PI,360.d0),
     &                u_wdot(i)
          enddo
 
   
          write (*,*) 'Best-fit V0 [m/s] and their error bars:'
          do j = 1,ndset
              i = 7*npl + j
 
              write (*,*) a(i)
              write (*,*) dsqrt(covar(i,i))
          enddo   
      
          write (*,*) 'Jitters (not fitted in chi^2 method):'
          do j = 1,ndset
              write (*,*) jitter(j)
              write (*,*) '0.0'
          enddo 
       
          write (*,*) 'linear trend [m/s per day]:'
          write (*,*) a(7*npl + ndset + 1)
          write (*,*) dsqrt(covar(7*npl + ndset + 1,7*npl + ndset + 1))

          write (*,*) 'quad. trend [m/s per day]:'
          write (*,*) a(7*npl + ndset + 2)
          write (*,*) dsqrt(covar(7*npl + ndset + 2,7*npl + ndset + 2))  
         
          write(*,*) ' ndata =',ndata
          write(*,*) ' mfit =',mfit
          write(*,*) ' RMS =',rms
          write(*,*) ' Chi^2 =',chisq/dble(ndata-mfit)
          write(*,*) ' epoch = ', t0
        

          do j = 1,npl+1

	        j_mass(j) = mass(j)/1.2667d17 
c	  s_mass(j) = mass(j)/1.32712497d20 
 
c          swift_mass(j) = (mass(j)/1.32712497d20)*
c     &          ((4.d0*PI*PI)/(365.25*365.25))
           enddo
      

           write(*,*) 'Jupiter mass'
c           write(*,*) (j_mass(i),i=1,npl+1)
           write(*,*) (j_mass(i+1),i=1,npl) 
                    
           write (*,*) 'semi-major axes in Jacobi'
           write (*,*)  (ap(i)/1.49597892d11,i=1,npl)
      else
           write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
           write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms
      endif
 

      if(writeflag_fit.gt.0) then 

          dt = (t_max+model_max - t0)/dble(nt - 1)
          do i = 1,nt
             x(i) = (i-1)*dt*8.64d4
          enddo
          call RVKEP_ewcop_fin (x,a,ymod,dyda,ma,nt,ia,epsil,deltat,hkl)
          do i = 1,nt
             write(*,*) t0 + x(i)/8.64d4,
     &       ymod(i) + a(7*npl +ndset + 1)*(x(i)/86400.d0)
     &               + a(7*npl +ndset +2)*(x(i)/86400.d0)**2
     
          enddo

      endif

      return

      end      


c MRQMIN attempts to reduce the chi^2 of a fit by the Levenberg-Marquardt
c method. It uses COVSRT, GAUSSJ, and MRQCOF.
c
c From Numerical Recipes.

	subroutine MRQMIN (x,ts,y,sig,ndata,a,ia,ma,covar,alpha,nca,
     & 	         chisq,funcs,alamda,loglik,jitt,epsil,deltat,hkl)

	implicit none
	integer ma,nca,ndata,ia(ma),MMAX,NDSMAX,ts(ndata)
	parameter (MMAX=200,NDSMAX=20)
        integer npl,ndset,idset,idsmax(NDSMAX),hkl
	real*8 alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     &	       sig(ndata),x(ndata),y(ndata),loglik,jitt(NDSMAX)
	external funcs

	integer j,k,l,mfit
	real*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX),epsil,deltat
	save ochisq,atry,beta,da,mfit

      common /DSBLK/ npl,ndset,idsmax,idset
	
c Initialization.
	if (alamda.lt.0.d0) then
	  mfit = 0
	  do j = 1,ma
	    if (ia(j).ne.0) mfit = mfit + 1
	  enddo
	  alamda = 0.001d0
cc	  alamda = 0.00001d0
CC	  alamda = 20000.d0
      call MRQCOF (x,ts,y,sig,ndata,a,ia,ma,alpha,beta,
     &       nca,chisq,funcs,loglik,jitt,epsil,deltat,hkl)
   
	  ochisq = chisq
 
	  do j = 1,ma
	    atry(j) = a(j)

	  enddo
	endif
	

c        pause
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
	
       
	call MRQCOF (x,ts,y,sig,ndata,atry,ia,ma,covar,da,
     &               nca,chisq,funcs,loglik,jitt,epsil,deltat,hkl)
	
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

	subroutine MRQCOF (x,ts,y,sig,ndata,a,ia,ma,alpha,beta,nalp,
     &	             chisq,funcs,loglik,jitt,epsil,deltat,hkl)

	implicit none
	integer npl,ndset,idset,ma,nalp,ndata,ia(ma),NDSMAX,MMAX
	parameter (NDSMAX=20, MMAX=200)
        integer idsmax(NDSMAX),ts(ndata),hkl
	real*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),
     &	       x(ndata),y(ndata),loglik,TWOPI,epsil,deltat
        parameter (TWOPI=2.d0*3.14159265358979d0)
	external funcs
	integer mfit,i,j,k,l,m
	real*8 dy,sig2i,wt,ymod(ndata),dyda(ndata,MMAX),jitt(NDSMAX)

        common /DSBLK/ npl,ndset,idsmax,idset

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

        call FUNCS (x,a,ymod,dyda,ma,ndata,ia,epsil,deltat,hkl)
        
	do i = 1,ndata
          idset = ts(i)
c          do j = 1,idset

c             ymod(i) = ymod(i) + a(7*npl+j)
c             dyda(i,7*npl+j) = 1.d0
c          enddo
       
          ymod(i) = ymod(i) + a(7*npl+idset)
          dyda(i,7*npl+idset) = 1.d0 
          
cccccccccccccccccccccc     lin. trend:     cccccccccccccccccccccccc
           
          ymod(i) = ymod(i) + a(7*npl +ndset + 1)*(x(i)/86400.d0)  
     &                      + a(7*npl +ndset + 2)*(x(i)/86400.d0)**2  
          dyda(i,7*npl + ndset + 1) = (x(i)/86400.d0)
          dyda(i,7*npl + ndset + 2) = (x(i)/86400.d0)**2


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c          do j = 1,idset-1
c             dyda(i,7*npl+j) = 0.d0
c          enddo
          do j = idset+1,ndset
             dyda(i,7*npl+j) = 0.d0
          enddo

c	  sig2i = 1.d0/(sig(i)*sig(i))
          sig2i = 1.d0/(sig(i)**2 + jitt(idset)**2)	  
	  
 
	  dy = y(i) - ymod(i)
	  j = 0
	  do l = 1,ma
	    if (ia(l).ne.0) then
	      j = j + 1
	      wt = dyda(i,l)*sig2i
	      k = 0
	      do m = 1,l
	        if (ia(m).ne.0) then
	          k = k + 1
	          alpha(j,k) = alpha(j,k) + wt*dyda(i,m)        
	        endif
	      enddo

	      beta(j) = beta(j) + dy*wt
	    endif
	  enddo
c	  Compute chi^2

          chisq = chisq + dy*dy*sig2i

	  loglik =  loglik - 0.5*dy*dy*sig2i -
     &    0.5*dlog(TWOPI*(sig(i)**2 + 
     &    jitt(idset)**2))   
	  
	  enddo
	  
	  j = 0

c Fill in the symmetric side.
	do j = 2,mfit
	  do k = 1,j-1
	    alpha(k,j) = alpha(j,k)
	  enddo
	enddo

	j = 0
	
	return
	end





c GAUSSJ solves linear equation by Gauss-Jordan elimination.
c
c From Numerical Recipes.

	subroutine GAUSSJ (a,n,np,b,m,mp)

	implicit none
	integer m,mp,n,np,NMAX
	real*8 a(np,np),b(np,mp)
	parameter (NMAX=51)
	integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
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
	    do l = 1,m
	      dum = b(irow,l)
	      b(irow,l) = b(icol,l)
	      b(icol,l) = dum

	    enddo
	  endif


	  
	  indxr(i) = irow
	  indxc(i) = icol
	  


	  if (a(icol,icol).ne.0.d0) then 

          
c	  Divide pivot row by the pivot element.
	  pivinv = 1.d0/a(icol,icol)

c	  endif


	  a(icol,icol) = 1.d0
	  do l = 1,n
	    a(icol,l) = a(icol,l)*pivinv
	  enddo
	  do l = 1,m
	    b(icol,l) = b(icol,l)*pivinv
	  enddo
      endif

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


 

      subroutine RVKEP_ewcop_fin (t,a,ymod,dyda,ma,ndata,ia,epsil,
     & deltat,hkl)
      implicit none
      real*8 PI,TWOPI,eps,epsil,deltat
      parameter (PI=3.14159265358979d0,eps=1.d-6)
      parameter (TWOPI=2.0d0*PI)
      integer npl,ma,i,j,NPLMAX,na,ndset,NDSMAX,idset,ndata,nbod,z,p
      parameter (NPLMAX=20,NDSMAX=20)
      real*8 t(ndata),ymod(ndata),a(ma),dyda(ndata,ma),x(ma),ia(ma)
      real*8 mstar,ap(NPLMAX),mass(NPLMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),vxh(NPLMAX),vyh(NPLMAX)
     &       ,vzh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),vxj(NPLMAX),vyj(NPLMAX)
     &       ,vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX)
      real*8 ah(ma),ahh(ma),ymodhb(ndata),ymodha(ndata)
      real*8 sini,sinih,sinihh,dydsini(ndata),factor
      integer ts(ndata),correct, idsmax(NDSMAX),hkl

      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      nbod = npl + 1
      na = 7*npl



      if (hkl.eq.0) then

          do i = 1,npl
             j = 7*(i-1)

c             a2(j+2) = 2.d0*PI/(a2(j+2)*8.64d4)

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
             if (a(j+4).lt.0.d0) a(j+4)=dmod(a(j+4)+2.d0*PI,2.d0*PI)
             if (a(j+5).lt.0.d0) a(j+5)=dmod(a(j+5)+2.d0*PI,2.d0*PI)
             if (a(j+4).gt.2.d0*PI) a(j+4)=dmod(a(j+4),2.d0*PI)
             if (a(j+5).gt.2.d0*PI) a(j+5)=dmod(a(j+5),2.d0*PI)
             if (a(j+6).lt.0.d0) a(j+6)=dmod(a(j+6)+2.d0*PI,2.d0*PI)
             if (a(j+7).lt.0.d0) a(j+7)=dmod(a(j+7)+2.d0*PI,2.d0*PI)
             if (a(j+6).gt.2.d0*PI) a(j+6)=dmod(a(j+6),2.d0*PI)
             if (a(j+7).gt.2.d0*PI) a(j+7)=dmod(a(j+7),2.d0*PI)
          enddo


      else

          do i = 1,npl
             j = 7*(i-1)

c             a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)

             if (a(j+2).lt.0.d0) then  ! if P<0, set P>0 
                a(j+2) = abs(a(j+2))
             endif
             
             if (a(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
                a(j+4) = -1.d0*a(j+4)       !     which is h = -h, k = -k
                a(j+3) = -1.d0*a(j+3)
                a(j+1) = abs(a(j+1))
             endif

             if (a(j+5).lt.0.d0) a(j+5)=dmod(a(j+5)+2.d0*PI,2.d0*PI)
             if (a(j+6).lt.0.d0) a(j+6)=dmod(a(j+6)+2.d0*PI,2.d0*PI)
             if (a(j+7).lt.0.d0) a(j+7)=dmod(a(j+7)+2.d0*PI,2.d0*PI)
             if (a(j+5).gt.2.d0*PI) a(j+5)=dmod(a(j+5),2.d0*PI)
             if (a(j+6).gt.2.d0*PI) a(j+6)=dmod(a(j+6),2.d0*PI)
             if (a(j+7).gt.2.d0*PI) a(j+7)=dmod(a(j+7),2.d0*PI)
          enddo

      endif
 
c-----get ymod first


      call MA_J_cop_fin (a,ma,npl,mstar,sini,mass,ap,hkl)
           
      
      call GENINIT_J3_ewcop (nbod,ap,a,
     &                 mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)
      
      call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)
     
        do i = 1,ndata        !  initialize ymod 
           ymod(i) = 0.0
        enddo 	

       
      call integrate_cop_fin(ymod,t,nbod,na,ndata,mass,a,
     &  xh,yh,zh,vxh,vyh,vzh,epsil,deltat,hkl)

     
      do i = 1,ma
         do j = 1,ndata
            dyda(j,i) = 0.0
         enddo
      enddo

       
      factor = 1.0d0

      do i = 0,npl-1
         ahh(7*i+1) = 0.05d0*factor 
         ahh(7*i+2) = 8640d0 
         ahh(7*i+3) = 0.0001d0*factor
         ahh(7*i+4) = 0.01d0*factor
         ahh(7*i+5) = 0.01d0*factor
         ahh(7*i+6) = 0.01d0*factor     
         ahh(7*i+7) = 0.01d0*factor
c         ahh(7*npl + ndset + 1) = 0.02d0*factor
      enddo
      
c       do i = 0,npl-1
c         ahh(7*i+1) = a(7*i+1)/1000.0d0
c         ahh(7*i+2) = a(7*i+2)/100.0d0
c         ahh(7*i+3) = a(7*i+3)/1000.0d0
c         ahh(7*i+4) = a(7*i+4)/10.0d0
c         ahh(7*i+5) = a(7*i+5)/10.0d0
c         ahh(7*i+6) = a(7*i+6)/10.0d0  
c         ahh(7*i+7) = a(7*i+7)/10.0d0
c      enddo     
      

      

c     calculate dyda
      do i = 1,na

         do j = 1,na 
            ah(j) = a(j)
         enddo       

c        get ah of back
c         if (mod(i,7).ne.3) then
         if (mod(i,7).ne.44) then
        

c           if(ia(i).lt.1.0d0) then 
           
c             do j = 1,ndata    
c               dyda(j,i) = ymod(j)
c             enddo 
                 
c             goto 1001
c           endif         

           
           ah(i) = a(i) - ahh(i)

           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap,hkl)


           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                 mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)

           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)
     
           call integrate_cop_fin(ymodhb,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat,hkl)
    
             
c        get ah of ahead
    
           ah(i) = a(i) + ahh(i)
       
           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap,hkl)
           
           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                 mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)

           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)

           call integrate_cop_fin(ymodha,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat,hkl)
         
         
c        calculate the ith dyda 
           do j = 1,ndata
               dyda(j,i) = (ymodha(j) - ymodhb(j))/(2.d0*ahh(i))          

           enddo
           
         else
         
         
           ah(i) = a(i) - ahh(i)
c           ah(i) = dmod(a(i) - ahh(i),  0.0d0 ) 

 

c           if (ah(j+6).le.0.d0) ah(j+6) = dmod(ah(j+6) + PI,  PI )  
c           if (ah(j+6).ge.PI)   ah(j+6) = dmod(ah(j+6) +0.00001,  PI) 

                       
         
           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap,hkl)

           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                 mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)
     
           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)

           call integrate_cop_fin(ymodhb,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat,hkl)
             
c        get ah of ahead
    
           ah(i) = a(i) + ahh(i)
c           ah(i) = dmod(a(i) + ahh(i),  PI ) 


           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap,hkl)
           
           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                 mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)

           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)

           call integrate_cop_fin(ymodha,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat,hkl)
         
 
c        calculate the ith dyda 
           do j = 1,ndata
               dyda(j,i) = (ymodha(j) - ymodhb(j))/(2.d0*ahh(i))

           enddo
c           pause
1001     endif
 
      enddo   
 
      
      return
      end
      
   
c MA_J calculates the actual masses and Jacobi semimajor axes of a two-planet
c system for assumed sin(i) using the parameters P, K and e from a
c two-Kepler fit.

	subroutine MA_J_cop_fin (a,ma,npl,m0,sini,mass,ap,hkl)
        
	implicit none
	real*8 m0,sini,PI,TWOPI,THIRD,GMSUN,dm,MSUN,ecc
        integer npl,ma,i,j,NPLMAX,hkl
        parameter (NPLMAX=20)
        real*8 a(ma),mass(NPLMAX),ap(NPLMAX),mpold(NPLMAX),mtotal
	parameter (THIRD=1.d0/3.d0)
        parameter (PI=3.14159265358979d0,TWOPI=2.d0*PI)
	parameter (GMSUN=1.32712497d20,MSUN=1.32712497d20)

c*******G is set to be unit, and s, m, kg as unit of time, length and mass
c*******expectively.        
                                          
 
        do i = 0,npl-1
        
           if (hkl.eq.0) then
               ecc = a(7*i+3)
           else
               ecc = dsqrt(a(7*i+3)**2+a(7*i+4)**2)    !! only for h, k
           endif

           mass(1) = m0
           mpold(i+1) = 0.d0
 101       continue
           if (i.eq.0) then
           mtotal = m0
           mass(i+2) = abs(a(7*i+1))*(a(7*i+2)*(m0 + mpold(i+1))**2/
     &               (TWOPI*GMSUN))**THIRD*dsqrt(1.d0-ecc**2)
     &               /dabs(dsin(a(7*i+6)))
           else
              mtotal = m0
              do j = 0, i-1
                 mtotal = mtotal + mass(j+2)
              enddo
              mass(i+2) = abs(a(7*i+1))*(a(7*i+2)*(mtotal
     &                  +mpold(i+1))**2/(TWOPI*GMSUN))**THIRD
     &                  *dsqrt(1.d0-ecc**2)/dabs(dsin(a(7*i+6)))
           endif
           
           dm = dabs(mass(i+2)-mpold(i+1))/mass(i+2)
           mpold(i+1) = mass(i+2)
           if (dm.gt.0) goto 101

           ap(i+1) = (GMSUN*(mtotal + mass(i+2))*(a(7*i+2)/TWOPI)
     &               **2)**THIRD
           
        enddo

        do i = 1,npl+1
           mass(i) = mass(i)*MSUN
 

        enddo
 
	return
	end

 
c GENINIT_J3 reads Jacobi orbital elements of nbod planets and generates
c initial position and velocity in Jacobi coords.
c This version outputs rpl and rhill.

c Last modified by Man Hoi Lee, Aug 16, 2003.

      subroutine GENINIT_J3_ewcop (nbod,ap,a,
     &                       mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)

c      include 'swift.inc'

      real*8 SMASSYR,MSUN,PI,eps,THIRD
      parameter (PI=3.14159265358979d0,eps=1.d-7)
      parameter (SMASSYR=4.d0*PI*PI)
      parameter (MSUN=1.32712497d20)
      parameter (THIRD=1.d0/3.d0)

      integer nbod,NPLMAX,i,j,hkl
      parameter (NPLMAX=20)
      real*8 mass(NPLMAX),ecc
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
      real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX)
      real*8 mstar0,m1,m2,frho3,ap(NPLMAX),a(NPLMAX)
      real*8 gm,inc,capom,a1,e1,omega,capm

      integer ialpha

c SET F/RHO^(1/3) FOR RADIUS (RHO IN G/CM^3) TO 1.D0 FOR NOW.
      frho3 = 1.d0

      xj(1) = 0.d0
      yj(1) = 0.d0
      zj(1) = 0.d0
      vxj(1) = 0.d0
      vyj(1) = 0.d0
      vzj(1) = 0.d0

      ialpha = -1
      inc = 0.d0
      capom = 0.d0
      gm = mass(1)

      do i = 2,nbod
         j = 7*(i-2)

         if (hkl.eq.0) then  
             ecc = a(j+3)
             omega = a(j+4)
             capm = a(j+5)
         else
             ecc = dsqrt(a(j+3)**2 + a(j+4)**2)
             omega = datan2(a(j+3),a(j+4))
             capm = a(j+5) - omega
         endif
         
        if (ecc.lt.0.d0) ecc = 0.000001d0 

         gm = gm + mass(i)
         rpl(i) = frho3*(1.5d0*mass(i)/2.d0*PI)**THIRD
         rhill(i) = ap(i-1)*(mass(i)/(3.d0*mass(1)))**THIRD
          
         call ORBEL_EL2XV (gm,ialpha,ap(i-1),ecc,a(j+6),a(j+7),
     &       omega,capm,xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i))

      enddo
      return
      end






      subroutine integrate_cop_fin(ymod,t,nbod,ma,ndata,mass,a,
     &                     xh,yh,zh,vxh,vyh,vzh,epsil,deltat,hkl)

      include 'swift_Jakub.inc'
      
        integer IO_NBITS
        parameter(IO_NBITS=6)
	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
	real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)
        real*8 xht(NPLMAX),yht(NPLMAX),zht(NPLMAX)
        real*8 vxht(NPLMAX),vyht(NPLMAX),vzht(NPLMAX)


	integer istat(NTPMAX,NSTAT),i1st
	integer nbod,ntp,nleft
	integer iflgchk,iub,iuj,iud,iue,hkl
        real*8 rstat(NTPMAX,NSTATR)

	real*8 t0,tstop,deltat,dtout,dtdump
	real*8 time,tout,tdump,tfrac,eoff

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        logical*1 lflg(0:IO_NBITS-1)

        integer ndata,ma,i,j,nd,flag,ierr
        real*8 ymod(ndata),t(ndata),a(ma)
        real*8 ap(NPLMAX),h,epsil,mtotal
        real*8 mstar,sini
         

	character*80 outfile,inparfile,inplfile,intpfile,fopenstat

 

        ntp = 0

c Prompt and read name of planet data file
      j2rp2 = 0.0d0
      j4rp4 = 0.0d0




c Initialize initial time and times for first output and first dump
	time = 0.d0
	tstop = t(ndata)
	

        iub = 20
        iuj = 30
        iud = 40
        iue = 60


        nleft = ntp
	i1st = 0

        do i = 1,ndata        !  initialize ymod 
           ymod(i) = 0.0
        enddo


c---------here is the big loop--------------------------

        nd = 1
c------output the first ymod if time of dataset begins from 0
        if (t(1).lt.1.d-10) then
           mtotal = 0.d0
           do i = 1,nbod
              mtotal = mass(i) + mtotal
           enddo
           do i = 2,nbod     
               j = 7*(i-2)         
c              ymod(nd) = ymod(nd) + mass(i)/mtotal*vyh(i)*a(ma)
c              ymod(nd) = ymod(nd) + mass(i)/mtotal*vzh(i)*dsin(a(j+6))

              ymod(nd) = ymod(nd) + mass(i)/mtotal*vzh(i)              
        
           enddo          
           nd = nd + 1
        endif
c-------loop-----
        do while( time.le.tstop )
           h = deltat
           flag = 0        ! flag for controling output
           do i = 1,ndata    

              if ((time.lt.t(i)).and.((t(i)-time).le.deltat)) then
                 flag = 1
                 h = t(i) - time
                 goto 555
              endif
           enddo

           
 555       call bs_step2(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &         xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &         vzht,istat,rstat,h,epsil)
           
           if (flag.eq.1) then
              mtotal = 0.d0
              do i = 1,nbod
                 mtotal = mtotal + mass(i)
              enddo
              do i = 2,nbod
                 j = 7*(i-2) 
c                 ymod(nd) = ymod(nd) + mass(i)/mtotal*vyh(i)*a(ma)
c                 ymod(nd) = ymod(nd) + mass(i)/mtotal*vzh(i)*
c     &            dsin(a(j+6))

             ymod(nd) = ymod(nd) + mass(i)/mtotal*vzh(i)   
     
              enddo          
              nd = nd + 1

           endif

           time = time + h
         
c           if(btest(iflgchk,4))  then    ! bit 4 is set
c              call discard(t,dt,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
c     &               xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
c     &               qmin,lclose,rplsq,istat,rstat)
c              call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
c     &               vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
c     &               'discard.out',fopenstat,nleft)
c           else
c c             nleft = ntp
c           endif

	enddo


        return
        
        end
        






      subroutine bs_step2(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt,eps)	

      include 'swift_Jakub.inc'
      include 'bs.inc'

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
      integer j,i,ntpi,link(NTPMAX)
      integer istattmp(NTPMAX,NSTAT),jj,i1stin
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX),eps
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
      real*8 ybs(6,(NTPMAX+NPLMAX)),tfake,dttmp,msys

      

c----
c...  Executable code 

c...  set things up if this is the initial call

c...  Convert to barycentric coords
      call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &     xb,yb,zb,vxb,vyb,vzb,msys)
c      call coord_h2b_tp(ntp,xht,yht,zht,vxht,vyht,vzht,
c     &     xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
c     &     xbt,ybt,zbt,vxbt,vybt,vzbt)

c...  copy to the big array
      do i=1,nbod
         ybs(1,i) = xb(i)
         ybs(2,i) = yb(i)
         ybs(3,i) = zb(i)
         ybs(4,i) = vxb(i)
         ybs(5,i) = vyb(i)
         ybs(6,i) = vzb(i)
      enddo

      ntpi = 0
c      do i=1,ntp
c         if(istat(i,1).eq.0) then
c            ntpi = ntpi + 1
c            j = ntpi + nbod
c            link(ntpi) = i
c            ybs(1,j) = xbt(i)
c            ybs(2,j) = ybt(i)
c            ybs(3,j) = zbt(i)
c            ybs(4,j) = vxbt(i)
c            ybs(5,j) = vybt(i)
c            ybs(6,j) = vzbt(i)
c            do jj = 1,NSTAT
c               istattmp(ntpi,jj) = istat(i,jj)
c            enddo
c         endif
c      enddo

      tfake = 0.0d0
      dttmp = dt

c      do while(tfake.lt.dt)
      do while( (abs(tfake-dt)/dt) .gt. 1.0e-7 )    ! just to be real safe
         call bs_int_pl(nbod,ntpi,mass,j2rp2,j4rp4,istattmp,
     &        tfake,dttmp,ybs,eps)

         dttmp = dt - tfake
      enddo

c...  put things back
      do i=1,nbod
         xb(i) = ybs(1,i)
         yb(i) = ybs(2,i)
         zb(i) = ybs(3,i)
         vxb(i) = ybs(4,i)
         vyb(i) = ybs(5,i)
         vzb(i) = ybs(6,i)
      enddo

c      do i=1,ntpi
c         j = i + nbod
c         xbt(link(i)) = ybs(1,j)
c         ybt(link(i)) = ybs(2,j)
c         zbt(link(i)) = ybs(3,j)
c         vxbt(link(i)) = ybs(4,j)
c         vybt(link(i)) = ybs(5,j)
c         vzbt(link(i)) = ybs(6,j)
c         do jj = 1,NSTAT
c            istat(link(i),jj) = istattmp(i,jj)
c         enddo
c      enddo

c...  Convert back to helio. coords at the end of the step
	call coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &         xh,yh,zh,vxh,vyh,vzh)
c	call coord_b2h_tp(ntp,xbt,ybt,zbt,vxbt,vybt,vzbt,
c     &         xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
c     &         xht,yht,zht,vxht,vyht,vzht)

      return

      end   ! bs_step
c------------------------------------------------------------------------






*****************************************************************************
*                          ORBEL_EL2XV.F
*****************************************************************************
*     PURPOSE: To compute cartesian positions and velocities given
*               central mass, ialpha ( = +1 for hyp., 0 for para. and
*               -1 for ellipse), and orbital elements.
C       input:
c            gm       ==> G times central mass (real scalar)
c	     ialpha   ==> conic section type ( see PURPOSE, integer scalar)
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
*       
c       Output:
c            x,y,z    ==>  position of object (real scalars)
c            vx,vy,vz ==>  velocity of object (real scalars)
c
*     ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
*     REMARKS: All angles are in RADIANS
*       
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 11, 1992.
*     REVISIONS: May 26 - now use better Kepler solver for ellipses
*                 and hyperbolae called EHYBRID.F and FHYBRID.F
***********************************************************************

	subroutine orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &                x,y,z,vx,vy,vz)

	include 'swift_loglik_Jakub.inc'

c...  Inputs Only: 
	integer ialpha
	real*8 gm,a,e,inc,capom,omega,capm

c...  Outputs:
	real*8 x,y,z,vx,vy,vz

c...  Internals:
	real*8 cape,capf,zpara,em1
	real*8 sp,cp,so,co,si,ci
	real*8 d11,d12,d13,d21,d22,d23
	real*8 scap,ccap,shcap,chcap
	real*8 sqe,sqgma,xfac1,xfac2,ri,vfac1,vfac2
        real*8 orbel_ehybrid,orbel_fhybrid,orbel_zget

c----
c...  Executable code 

        if(e.lt.0.0) then
           write(*,*) ' ERROR in orbel_el2xv: e<0, setting e=0!!1'
           e = 0.0
        endif

c...    check for inconsistencies between ialpha and e
        em1 = e - 1.d0
        if(
     &     ((ialpha.eq.0) .and. (abs(em1).gt.TINY))  .or.
     &     ((ialpha.lt.0) .and. (e.gt.1.0d0))  .or.
     &     ((ialpha.gt.0) .and. (e.lt.1.0d0)) )  then
        write(*,*) 'ERROR in orbel_el2xv: ialpha and e inconsistent'
             write(*,*) '                       ialpha = ',ialpha
             write(*,*) '                            e = ',e
        endif

C Generate rotation matrices (on p. 42 of Fitzpatrick)
C
	call orbel_scget(omega,sp,cp)
	call orbel_scget(capom,so,co)
	call orbel_scget(inc,si,ci)
	d11 = cp*co - sp*so*ci
	d12 = cp*so + sp*co*ci
	d13 = sp*si
	d21 = -sp*co - cp*so*ci
	d22 = -sp*so + cp*co*ci
	d23 = cp*si

C--
C Get the other quantities depending on orbit type ( i.e. IALPHA)
C
	if (ialpha .eq. -1) then
	  cape = orbel_ehybrid(e,capm)
	  call orbel_scget(cape,scap,ccap)
	  sqe = sqrt(1.d0 -e*e)
	  sqgma = sqrt(gm*a)
	  xfac1 = a*(ccap - e)
	  xfac2 = a*sqe*scap
	  ri = 1.d0/(a*(1.d0 - e*ccap))
	  vfac1 = -ri * sqgma * scap
	  vfac2 = ri * sqgma * sqe * ccap
	endif
c--
	if (ialpha .eq. +1) then
	  capf = orbel_fhybrid(e,capm)
	  call orbel_schget(capf,shcap,chcap)
	  sqe = sqrt(e*e - 1.d0 )
	  sqgma = sqrt(gm*a)
	  xfac1 = a*(e - chcap)
	  xfac2 = a*sqe*shcap
	  ri = 1.d0/(a*(e*chcap - 1.d0))
	  vfac1 = -ri * sqgma * shcap
	  vfac2 = ri * sqgma * sqe * chcap
	endif
C--
	if (ialpha .eq. 0) then
	  zpara = orbel_zget(capm)
	  sqgma = sqrt(2.d0*gm*a)
	  xfac1 = a*(1.d0 - zpara*zpara)
	  xfac2 = 2.d0*a*zpara
	  ri = 1.d0/(a*(1.d0 + zpara*zpara))
	  vfac1 = -ri * sqgma * zpara
	  vfac2 = ri * sqgma 
	endif
C--
	x =  d11*xfac1 + d21*xfac2
	y =  d12*xfac1 + d22*xfac2
	z =  d13*xfac1 + d23*xfac2
	vx = d11*vfac1 + d21*vfac2
	vy = d12*vfac1 + d22*vfac2
	vz = d13*vfac1 + d23*vfac2

	return
	end    ! orbel_el2xv

c-----------------------------------------------------------------------


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
c                    ORBEL_FHYBRID.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           n ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
*	         For larger N, uses FGET
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26,1992.
*     REVISIONS: 
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_fhybrid(e,n)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
	real*8 e,n

c...  Internals:
	real*8 abn
        real*8 orbel_flon,orbel_fget

c----
c...  Executable code 

	abn = n
	if(n.lt.0.d0) abn = -abn

	if(abn .lt. 0.636d0*e -0.6d0) then
	  orbel_fhybrid = orbel_flon(e,n)
	else 
	  orbel_fhybrid = orbel_fget(e,n)
	endif   

	return
	end  ! orbel_fhybrid
c-------------------------------------------------------------------

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




c***************************************************************************
c			TU4_GETACCB.F
c*************************************************************************
c GETACCB returns the bary. acc. on each of n mutually
c interacting objects by direct pairwise summation
c	
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of planets (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xb,yb,zb      ==>  position of planets in beri coord 
c                                    (real arrays)
c             Output:
c               axb,ayb,azb   ==>  accel in beri coord (real arrays) 
c
c Remarks:  Based on Martin's NB4M routines
c Authors:  Martin Duncan 
c Date:    3/8/93
c Last revision: 4/5/95

      subroutine tu4_getaccb(nbod,mass,j2rp2,j4rp4,xb,yb,zb,axb,ayb,azb)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xb(nbod),yb(nbod),zb(nbod)

c...  Output
      real*8 axb(nbod),ayb(nbod),azb(nbod)

c...  Internals
      real*8 xx,yy,zz,rr2,fac,mi,mj
      real*8 axx,ayy,azz,fac1
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),irh(NPLMAX)
      real*8 aoblx(NPLMAX),aobly(NPLMAX),aoblz(NPLMAX)
      integer i,j

c----
c...  executable code

      xh(1) = 0.0d0
      yh(1) = 0.0d0
      zh(1) = 0.0d0

c...  do sun part first
      axb(1) = 0.0d0
      ayb(1) = 0.0d0
      azb(1) = 0.0d0
      i = 1
      do j = i+1,nbod
         mi = mass(i)
         mj = mass(j)
         xx = xb(i) - xb(j)
         yy = yb(i) - yb(j)
         zz = zb(i) - zb(j)
         rr2 = xx**2 + yy**2 + zz**2

         fac1 = 1.d0/sqrt(rr2)
         fac = fac1/rr2
c         if (j.eq.2) write(*,*) fac, rr2, xx , yy , zz

c..      save for the J2 and J4 calculations
         xh(j) = -xx
         yh(j) = -yy
         zh(j) = -zz
         irh(j) = fac1

         axx = xx*fac
         ayy = yy*fac
         azz = zz*fac
         axb(i) = axb(i) - axx*mj
         ayb(i) = ayb(i) - ayy*mj
         azb(i) = azb(i) - azz*mj
         axb(j) =   axx*mi
         ayb(j) =   ayy*mi
         azb(j) =   azz*mi
      enddo

      do i = 2,nbod-1
         do j = i+1,nbod
            mi = mass(i)
            mj = mass(j)
            xx = xb(i) - xb(j)
            yy = yb(i) - yb(j)
            zz = zb(i) - zb(j)
            rr2 = xx**2 + yy**2 + zz**2
            fac = 1.d0/(rr2*sqrt(rr2))
            axx = xx*fac
            ayy = yy*fac
            azz = zz*fac
            axb(i) = axb(i) - axx*mj
            ayb(i) = ayb(i) - ayy*mj
            azb(i) = azb(i) - azz*mj
            axb(j) = axb(j) + axx*mi
            ayb(j) = ayb(j) + ayy*mi
            azb(j) = azb(j) + azz*mi
         enddo
      enddo

      if(j2rp2.ne.0.0d0) then
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
         do i = 1,nbod
            axb(i) = axb(i) + aoblx(i)
            ayb(i) = ayb(i) + aobly(i)
            azb(i) = azb(i) + aoblz(i)
         enddo
      endif

      return	
      end                       !  tu4_getaccb
c____________________________________________________________________________



c***********************************************************************
c	                    COORD_H2B.F
c***********************************************************************
*     PURPOSE: Converts from Heliocentric to Barycentric coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     xh(*),yh(*),zh(*) ==> heliocentric particle coords
*                                          (real array)
*		     vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
*                                             (real array)
*                 Returned are
*                    xb(*),yb(*),zb(*) ==> bary. particle positions
*                                          (real array)
*                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
*                                            (real array)
*                    msys              ==>  Total mass of of system
*                                            (real scalar)       
*     Authors:  Martin Duncan
*     ALGORITHM: Obvious 
*     WRITTEN:  Jan 27/93
*     REVISIONS: 2/22/94  HFL

	subroutine coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &      xb,yb,zb,vxb,vyb,vzb,msys)

    include 'swift_loglik_Jakub.inc'


c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c...  Outputs:
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c...  Internals:
	real*8 msys,xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp
	integer n

c----
c...  Executable code 

	msys = mass(1)
	xtmp =0.d0
	ytmp =0.d0
	ztmp =0.d0
	vxtmp =0.d0
	vytmp =0.d0
	vztmp =0.d0

	do n=2,nbod
	   msys = msys +mass(n)
	   xtmp = xtmp + mass(n)*xh(n)
	   ytmp = ytmp + mass(n)*yh(n)
	   ztmp = ztmp + mass(n)*zh(n)
	   vxtmp = vxtmp + mass(n)*vxh(n)
	   vytmp = vytmp + mass(n)*vyh(n)
	   vztmp = vztmp + mass(n)*vzh(n)
	enddo

	xb(1) = -xtmp/msys
	yb(1) = -ytmp/msys
	zb(1) = -ztmp/msys
	vxb(1) = -vxtmp/msys
	vyb(1) = -vytmp/msys
	vzb(1) = -vztmp/msys

	do n=2,nbod
	  xb(n) = xh(n) + xb(1)
	  yb(n) = yh(n) + yb(1)
	  zb(n) = zh(n) + zb(1)
	  vxb(n) = vxh(n) + vxb(1)
	  vyb(n) = vyh(n) + vyb(1)
	  vzb(n) = vzh(n) + vzb(1)
	enddo

	return
	end     ! coord_h2b
c--------------------------------------------------------------------------


c***********************************************************************
c	                    COORD_B2H.F
c***********************************************************************
*     PURPOSE: Converts from Barycentric to Helio coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
*                                 SYMMETRY IN SUBROUTINE CALLS
*		     xb(*),yb(*),zb(*) ==> Barycentric particle coords
*                                          (real array)
*		     vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    xh(*),yh(*),zh(*) ==> Helio particle positions
*                                          (real array)
*                    vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
*                                            (real array)
*       
*     ALGORITHM: Obvious 
*     REMARKS:  Can of course use this to get coords. relative to any body.
*              by changing the one subtracted off.
*
*     Authors:  Martin Duncan
*     WRITTEN:  Jan 27/93
*     REVISIONS: 2/17/95  HFL

	subroutine coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &      xh,yh,zh,vxh,vyh,vzh)


    include 'swift_loglik_Jakub.inc'


c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c...  Outputs:
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c...  Internals:
	integer n

c----
c...  Executable code 

	do n=1,nbod
	  xh(n) = xb(n) - xb(1)
	  yh(n) = yb(n) - yb(1)
	  zh(n) = zb(n) - zb(1)
	  vxh(n) = vxb(n) - vxb(1)
	  vyh(n) = vyb(n) - vyb(1)
	  vzh(n) = vzb(n) - vzb(1)
	enddo

	return
	end     ! coord_b2h

c--------------------------------------------------------------------------


c***************************************************************************
c			OBL_ACC.F
c*************************************************************************
c OBL_ACC returns the BARYCENTRIC x,y,z components of the acc. on NBOD
c particles due to the oblateness of mass(1) using  
c the values of J2RP2 and J4RP4 passed into the routine.
c (J2RP2 for example is the product of 
c J_2 times the square of the central body's radius)
c Here we return the net acc. produced
c only by the J2 and J4 terms (i.e. including
c neither the monopole nor higher order terms).
c	
c
c             Input:
c                 nbod     ==>  number of massive bodies (incl. central one)
c                 mass(*)  ==>  masses of particles (real*8 array)
c                 j2rp2    ==>  scaled value of j2 moment (real*8 scalar)
c                 j4rp4    ==>  scaled value of j4 moment (real*8 scalar)
c                                    (real*8 vectors)
c                 xh(*),yh(*),zh(*)   ==>  HELIO. positions of particles
c                 irh(*)   ==> 1./ magnitude of radius vector (real*8 vector)
c                                (passed in to save calcs.)
c             Output:
c               aoblx(*),aobly(*),aoblz(*)  ==>  BARY. components of accel 
c                                        (real*8 vectors) 
c
c Remarks:  aoblx(1) (for example) contains x-component of
c           bary. acc. of central body
c Authors:  Martin Duncan 
c Date:    3/4/94
c Last revision: 

      subroutine obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &     aoblx,aobly,aoblz)


    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
      integer nbod
      real*8 j2rp2,j4rp4
      real*8 mass(NPLMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),irh(NPLMAX)

c...  Output
      real*8 aoblx(NPLMAX),aobly(NPLMAX),aoblz(NPLMAX)

c...  Internals
      integer n
      real*8 rinv2,t0,t1,t2,t3
      real*8 fac1,fac2

c----
c...  executable code

c First get the bary acc. of each "planet" due to central oblate "sun"

	do n =2,nbod

c Note that here we assume we know inverse of radius rather than calc. it
c from (x,y,z) to save the sqrt.
	  rinv2= irh(n)**2
	  t0 = -mass(1)*rinv2*rinv2*irh(n)
	  t1 = 1.5d0 *j2rp2
	  t2 = zh(n)*zh(n)*rinv2
	  t3 = 1.875d0 *j4rp4*rinv2

	  fac1 = t0*(t1 - t3 - (5.d0*t1 - (14.d0 - 21.d0*t2)*t3)*t2)
	  fac2 = 2.d0*t0*(t1 - (2.d0 - (14.d0*t2/3.d0))*t3)
      
          aoblx(n) = fac1*xh(n)
          aobly(n) = fac1*yh(n)
          aoblz(n) = (fac1 + fac2)*zh(n)

	enddo

c Now compute the bary. acc. of Sun due to all the planets

	aoblx(1) = 0.d0
	aobly(1) = 0.d0
	aoblz(1) = 0.d0
	do n=2,nbod
	   aoblx(1) = aoblx(1) - mass(n)*aoblx(n)/mass(1)
	   aobly(1) = aobly(1) - mass(n)*aobly(n)/mass(1)
	   aoblz(1) = aoblz(1) - mass(n)*aoblz(n)/mass(1)
	enddo

        return	
        end                       !  obl_acc.f
c____________________________________________________________________________


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
c                    ORBEL_FGET.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_fget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
*           Cel. Mech. ".  Quartic convergence from Danby's book.
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_fget(e,capn)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
	real*8 e,capn

c...  Internals:
	integer i,IMAX
	real*8 tmp,x,shx,chx
	real*8 esh,ech,f,fp,fpp,fppp,dx
	PARAMETER (IMAX = 10)

c----
c...  Executable code 

c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. 

c  begin with a guess proposed by Danby	
	if( capn .lt. 0.d0) then
	   tmp = -2.d0*capn/e + 1.8d0
	   x = -log(tmp)
	else
	   tmp = +2.d0*capn/e + 1.8d0
	   x = log( tmp)
	endif

	orbel_fget = x

	do i = 1,IMAX
	  call orbel_schget(x,shx,chx)
	  esh = e*shx
	  ech = e*chx
	  f = esh - x - capn
c	  write(6,*) 'i,x,f : ',i,x,f
	  fp = ech - 1.d0  
	  fpp = esh 
	  fppp = ech 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
	  orbel_fget = x + dx
c   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) RETURN
	  x = orbel_fget
	enddo	

	write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	return
	end   ! orbel_fget
c------------------------------------------------------------------

***********************************************************************
c                    ORBEL_ZGET.F
***********************************************************************
*     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
*          given Q (Fitz. notation.)
*
*             Input:
*                           q ==>  parabola mean anomaly. (real scalar)
*             Returns:
*                  orbel_zget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
*     REMARKS: For a parabola we can solve analytically.
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: May 27 - corrected it for negative Q and use power
*	      series for small Q.
***********************************************************************

	real*8 function orbel_zget(q)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
	real*8 q

c...  Internals:
	integer iflag
	real*8 x,tmp

c----
c...  Executable code 

	iflag = 0
	if(q.lt.0.d0) then
	  iflag = 1
	  q = -q
	endif

	if (q.lt.1.d-3) then
	   orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
	else
	   x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
	   tmp = x**(1.d0/3.d0)
	   orbel_zget = tmp - 1.d0/tmp
	endif

	if(iflag .eq.1) then
           orbel_zget = -orbel_zget
	   q = -q
	endif
	
	return
	end    ! orbel_zget
c----------------------------------------------------------------------


***********************************************************************
c	                  ORBEL_SCHGET.F
***********************************************************************
*     PURPOSE:  Given an angle, efficiently compute sinh and cosh.
*
*        Input:
*             angle ==> angle in radians (real scalar)
*        
*        Output:
*             shx    ==>  sinh(angle)  (real scalar)
*             chx    ==>  cosh(angle)  (real scalar)
*
*     ALGORITHM: Obvious from the code 
*     REMARKS: Based on the routine SCGET for sine's and cosine's.
*       We use the sqrt rather than cosh (it's faster)
*       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
*       OR OVERFLOWS WILL OCCUR!
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 6, 1992.
*     REVISIONS: 
***********************************************************************

	subroutine orbel_schget(angle,shx,chx)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
        real*8 angle

c...  Output:
	real*8 shx,chx

c----
c...  Executable code 

	shx = sinh(angle)
	chx= sqrt(1.d0 + shx*shx)

	return
	end   ! orbel_schget
c---------------------------------------------------------------------



***********************************************************************
c                    ORBEL_FLON.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_flon ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Uses power series for N in terms of F and Newton,s method
*     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26, 1992.
*     REVISIONS: 
***********************************************************************

	real*8 function orbel_flon(e,capn)

    include 'swift_loglik_Jakub.inc'


c...  Inputs Only: 
	real*8 e,capn

c...  Internals:
	integer iflag,i,IMAX
	real*8 a,b,sq,biga,bigb
	real*8 x,x2
	real*8 f,fp,dx
	real*8 diff
	real*8 a0,a1,a3,a5,a7,a9,a11
	real*8 b1,b3,b5,b7,b9,b11
	PARAMETER (IMAX = 10)
	PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
	PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
	PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
	PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

c----
c...  Executable code 


c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. Only good for smallish CAPN 

	iflag = 0
	if( capn .lt. 0.d0) then
	   iflag = 1
	   capn = -capn
	endif

	a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
	a0 = -6227020800.d0*capn/e
	b1 = a1

c  Set iflag nonzero if capn < 0., in which case solve for -capn
c  and change the sign of the final answer for F.
c  Begin with a reasonable guess based on solving the cubic for small F	


	a = 6.d0*(e-1.d0)/e
	b = -6.d0*capn/e
	sq = sqrt(0.25*b*b +a*a*a/27.d0)
	biga = (-0.5*b + sq)**0.3333333333333333d0
	bigb = -(+0.5*b + sq)**0.3333333333333333d0
	x = biga + bigb
c	write(6,*) 'cubic = ',x**3 +a*x +b
	orbel_flon = x
c If capn is tiny (or zero) no need to go further than cubic even for
c e =1.
	if( capn .lt. TINY) go to 100

	do i = 1,IMAX
	  x2 = x*x
	  f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
	  fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
	  dx = -f/fp
c	  write(6,*) 'i,dx,x,f : '
c	  write(6,432) i,dx,x,f
432	  format(1x,i3,3(2x,1p1e22.15))
	  orbel_flon = x + dx
c   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) go to 100
	  x = orbel_flon
	enddo	

c Abnormal return here - we've gone thru the loop 
c IMAX times without convergence
	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif
	write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	  diff = e*sinh(orbel_flon) - orbel_flon - capn
	  write(6,*) 'N, F, ecc*sinh(F) - F - N : '
	  write(6,*) capn,orbel_flon,diff
	return

c  Normal return here, but check if capn was originally negative
100	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif

	return
	end     ! orbel_flon
c------------------------------------------------------------------


c***********************************************************************
c	                    COORD_J2H.F
c***********************************************************************
*     PURPOSE: Converts from Jacobi to Helio coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     xj(*),yj(*),zj(*) ==> Jacobi particle coords
*                                          (real array)
*		     vxj(*),vyj(*),vzj(*) ==> Jacobi particle velocities
*                                             (real array)
*                 Returned are
*                    xh(*),yh(*),zh(*) ==> Helio particle positions
*                                          (real array)
*                    vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
*                                            (real array)
*       
*     ALGORITHM: See my notes on Nov 21. 
*
*     Authors:  Martin Duncan
*     WRITTEN:  Jan 27/93
*     REVISIONS: 2/17/95  HFL

	subroutine coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xh,yh,zh,vxh,vyh,vzh)


    include 'swift_loglik_Jakub.inc'


c...  Inputs: 
	integer nbod
	real*8 mass(NPLMAX)
	real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)

c...  Outputs:
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c...  Internals:
	integer n
	real*8 sumx,sumy,sumz,sumvx,sumvy,sumvz
	real*8 eta(NPLMAX)

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi coords

	eta(1) = mass(1)
	do n = 2,nbod
	   eta(n) = eta(n-1) + mass(n)
        enddo

        xh(1) = 0.d0
        yh(1) =  0.d0
	zh(1) =  0.d0
	vxh(1) =  0.d0
	vyh(1) =  0.d0
	vzh(1) =  0.d0

        xh(2) = xj(2) 
        yh(2) = yj(2)
	zh(2) = zj(2)
	vxh(2) = vxj(2)
	vyh(2) = vyj(2)
	vzh(2) = vzj(2)

	sumx = mass(2)*xj(2)/eta(2)
	sumy = mass(2)*yj(2)/eta(2)
	sumz = mass(2)*zj(2)/eta(2)
	sumvx = mass(2)*vxj(2)/eta(2)
	sumvy = mass(2)*vyj(2)/eta(2)
	sumvz = mass(2)*vzj(2)/eta(2)

	do n=3,nbod 
	  xh(n) = xj(n) + sumx
	  yh(n) = yj(n) + sumy
	  zh(n) = zj(n) + sumz
	  vxh(n) = vxj(n) + sumvx
	  vyh(n) = vyj(n) + sumvy
	  vzh(n) = vzj(n) + sumvz

	  if(n.lt.nbod) then
             sumx = sumx + mass(n)*xj(n)/eta(n)
             sumy = sumy + mass(n)*yj(n)/eta(n)
             sumz = sumz + mass(n)*zj(n)/eta(n)
             sumvx = sumvx + mass(n)*vxj(n)/eta(n)
             sumvy = sumvy + mass(n)*vyj(n)/eta(n)
             sumvz = sumvz + mass(n)*vzj(n)/eta(n)
          endif

	enddo


123	return
	end     ! coord_j2h

c--------------------------------------------------------------------------


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

      subroutine bs_int_pl(nbod,ntp,mass,j2rp2,j4rp4,istat,x,h0,y,eps)

      include 'swift_loglik_Jakub.inc'
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
      call bs_der_pl(ntp,nbod,mass,j2rp2,j4rp4,istat,y,dy)
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
               call bs_der_pl(ntp,nbod,mass,j2rp2,j4rp4,istat,y,dy)
               do i=1,n
                  ii=12*i
                  tp(ii-1)=dmax1(tp(ii-1),dabs(y(i)))
                  eta2=tp(ii-3)+h*dy(i)
                  tp(ii-3)=y(i)
                  y(i)=eta2
               enddo 
            enddo
         
            call bs_der_pl(ntp,nbod,mass,j2rp2,j4rp4,istat,y,dy)
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


C**************************************************************************
C	    		        BS_DER
C**************************************************************************
c This is the subroutine that calculates the derivatives of the independant var
c
c             Input:
c              nbod  ==> number of planets  (int scalar)
c              ntp   ==> number of test particles  (int scalar)
c              mass  ==>  mass of bodies (real array)
c      j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c             istat  ==>  status of the test paricles
c                           (2d integer array)
c                            istat(i,1) = 0 ==> active:  = 1 not
c                            istat(i,2) = -1 ==> Danby did not work
c 	        ybs  ==> values dependent variables  (real array)
c
c             Output:
c 	         dy  ==> derivatives of the independant var (real array)
c
c Remarks:  This used TU4 routines !!  
c Authors:  Hal Levison
c Date:    5/17/93
c Last revision: 2/24/94

      subroutine bs_der_pl(ntp,nbod,mass,j2rp2,j4rp4,istat,ybs,dy)

      include 'swift_loglik_Jakub.inc'
      include 'bs.inc'

c...  Inputs Only: 
      integer nbod,ntp
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 ybs(6,(NTPMAX+NPLMAX))

c...  Input and Outputs
      integer istat(NTPMAX,NSTAT)

c...  Output
      real*8 dy(6,(NTPMAX+NPLMAX))

c...  Internals
      integer i,j
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      real*8 axb(NPLMAX),ayb(NPLMAX),azb(NPLMAX)
      real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
      real*8 axbt(NTPMAX),aybt(NTPMAX),azbt(NTPMAX)

c----
c...  Executable code 

c...  move things so that I can deal with it
      do i=1,nbod
         xb(i) = ybs(1,i)
         yb(i) = ybs(2,i)
         zb(i) = ybs(3,i)
         vxb(i) = ybs(4,i)
         vyb(i) = ybs(5,i)
         vzb(i) = ybs(6,i)
      enddo

c      do i=1,ntp
c         j = i + nbod
c         xbt(i) = ybs(1,j)
c         ybt(i) = ybs(2,j)
c         zbt(i) = ybs(3,j)
c         vxbt(i) = ybs(4,j)
c         vybt(i) = ybs(5,j)
c         vzbt(i) = ybs(6,j)
c      enddo

      call tu4_getaccb(nbod,mass,j2rp2,j4rp4,xb,yb,zb,axb,ayb,azb)
c      call tu4_getaccb_tp(nbod,mass,j2rp2,j4rp4,xb,yb,zb,
c     &     ntp,xbt,ybt,zbt,istat,axbt,aybt,azbt)


c.... moves things intp dy array
      do i=1,nbod
         dy(1,i) = ybs(4,i)
         dy(2,i) = ybs(5,i)
         dy(3,i) = ybs(6,i)
         dy(4,i) = axb(i)
         dy(5,i) = ayb(i)
         dy(6,i) = azb(i)
      enddo

c      do i=1,ntp
c         j = i + nbod
c         dy(1,j) = ybs(4,j)
c         dy(2,j) = ybs(5,j)
c         dy(3,j) = ybs(6,j)
c         dy(4,j) = axbt(i)
c         dy(5,j) = aybt(i)
c         dy(6,j) = azbt(i)
c      enddo

      return
      end     ! bs_der
c-------------------------------------------------------------------------------


