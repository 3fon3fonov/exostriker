ccc   Innitually writen by Man Hoi Lee, Xianyu Tan, 
ccc   Heavlly modified by Trifon Trifonov trifon@hku.hk
ccc   The final version will be available in the Python RVMod lib.    


      implicit none
      real*8 PI
      parameter (PI=3.14159265358979d0)
      integer npl,ndset,idset,ndata,ma,mfit,i,j,NDSMAX,NPLMAX,MMAX,k
      real*8 inc,sini,mstar
      integer writeflag_best_par, nt
      integer writeflag_RV,writeflag_fit, amoebastarts
      parameter (NDSMAX=20, NPLMAX=20, MMAX=200)
      integer idsmax(NDSMAX),ia(MMAX),ts(10000)
      real*8 t(10000),y(10000),sig(10000),ys(10000),sigs(10000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 chisq,alamda,ochisq,dchisq,loglik,epsil,deltat
      real*8 jitter(NDSMAX),sigscale,t0,t_max,epoch
      real*8 rms,ymod(10000),dyda(10000,MMAX)
      real*8 a0(MMAX),asave(MMAX),ap(NPLMAX)
 
      real*8 tstop,dt,dtout
      real*4 t_stop,when_to_kill,model_max
       
      external rvkep_ewcop_fin
      character*80 infile

      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      loglik=0.d0
      rms=0.d0
     
c amoebastarts for consistency with loglik input
     
      read (*,*) epsil,deltat, amoebastarts,
     &          when_to_kill, nt, model_max
     
c      write(*,*) 'Stellar mass: '
      read (*,*) mstar,
     &          writeflag_best_par, 
     &	             writeflag_RV,writeflag_fit 
   

      call io_read_data (ndata,t,ts,ys,sigs,jitter,
     & 	           epoch,t0,t_max,a,ia,ma,mfit)



c      write(*,*) ' i_max,i_min,isteps: '
c      read (*,*) i_max,i_min,istep
c      write(*,*) i_max,i_min,istep

c     call io_read_initpa_ewcop_fin (a,ia,t,ma,ndata,mfit)
       
c      write(*,*)            
c*****set alamda to be negtive for initializing******
      alamda = -1.d0
   
      call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
     & chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat)
      
c     write(*,*) ' alamda,chi_nu^2: ',alamda,chisq/dble(ndata-mfit)

 

c      pause

      i = 0
 500  continue
          i = i + 1
          ochisq = chisq
          call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
     &  chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat)
 
 
c          write(*,*) ' alamda,chi_nu^2: ',alamda,chisq/dble(ndata-mfit)

          dchisq = chisq - ochisq

          if ((i.eq.200).or.(alamda.ge.1d6)) then
              i = 0
              goto 502
c              pause
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
     &  chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat)

      
c          write(*,*) ' alamda,chi_nu^2: ',alamda,chisq/dble(ndata-mfit)
          
          
c          if (alamda.gt.1d12) alamda = -1.d0
          dchisq = chisq - ochisq

          if ((i.eq.200).or.(alamda.ge.1d6)) then
c          if (i.eq.200) then
              i = 0
              alamda = -1.d0
              goto 333
c              pause
          endif
      if ((chisq.ge.ochisq).or.(dchisq.lt.-1d-5)) then
      goto 501
      endif


c*******final output******************
 333  alamda = 0.d0

      call MRQMIN (t,ts,ys,sigs,ndata,a,ia,ma,covar,alpha,MMAX,
     & chisq,rvkep_ewcop_fin,alamda,loglik,jitter,epsil,deltat)

     
      call io_write_bestfitpa_ewcop_fin (a,covar,t,ys,ndata,ts,
     & 	           ma,mfit,t0,t_max,sigs,chisq,rms,writeflag_RV,
     &  writeflag_best_par,writeflag_fit,jitter,epsil,deltat,
     &  nt, model_max)

      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
      write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms
     
      stop
      end



      subroutine io_read_data (ndata,t,ts,ys,sigs,jitter,epoch,
     &               t0,t_max,ar,iar,ma,mfit)  


      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl
      real*8 t(10000),y(10000),sig(10000),ys(10000),sigs(10000)
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      integer idsmax(NDSMAX),ts(10000)
      real*8 jitter(NDSMAX),t0,t_max,epoch,ar(MMAX),off(NDSMAX), PI
      parameter(PI=3.14159265358979d0)
      integer i,k,j, iar(MMAX), u_off(NDSMAX), u_jit, ma, mfit
      character*80 infile
   
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

      ma = 7*npl + ndset + 1
      
c      write(*,*) 'Initial K, P, e, w, M0,Inc,Capom and their flags: '
      do j = 1,npl
          i = 7*(j-1)
          read (*,*) ar(i+1),ar(i+2),ar(i+3),ar(i+4),ar(i+5),ar(i+6),
     &               ar(i+7)
          read (*,*) iar(i+1),iar(i+2),iar(i+3),iar(i+4),iar(i+5),
     &               iar(i+6),iar(i+7)
 
 
c          a(i+2) = 2.d0*PI/(a(i+2)*8.64d4)         ! mean motion 
          ar(i+2) = ar(i+2)*8.64d4         ! second as unit
          ar(i+4) = ar(i+4)*PI/180.d0                ! radians
          ar(i+5) = ar(i+5)*PI/180.d0
          ar(i+6) = ar(i+6)*PI/180.d0
          ar(i+7) = ar(i+7)*PI/180.d0
 
      enddo
      
c      write (*,*) 'linear trend:'      
      read (*,*) ar(7*npl + ndset + 1)
      read (*,*) iar(7*npl + ndset + 1)   
      

      read (*,*) epoch
 
      t_max = t(ndata) 

      if (epoch.eq.0) then 
         t0 = t(1)
      else
         t0 = epoch
      endif
      
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
C
Cc    Xianyu Tan 2011

      subroutine io_write_bestfitpa_ewcop_fin (a,covar,t,ys,ndata,ts,
     &           ma,mfit,t0,t_max,sigs,chisq,rms,writeflag_RV,
     &           writeflag_best_par,writeflag_fit,jitter,epsil,
     &           deltat,nt, model_max)
   
      implicit none 
      real*8 PI
      integer MMAX,NDSMAX,NPLMAX 
      parameter (PI=3.14159265358979d0,MMAX=200  ,NDSMAX=20,NPLMAX=20)
      real*8 a(MMAX),ia(MMAX),t(5000),ymod(5000),ys(5000)
      real*8 covar(MMAX,MMAX),dyda(5000,MMAX),AU,day
      real*8 rms,mstar,sini(7),mass(NPLMAX),ap(NPLMAX)
      integer ts(5000),nbod,nt,writeflag_RV,
     &           writeflag_best_par,writeflag_fit
      real*8 t0,t1,t2,dt,offset,t_max,chisq,deltat,epsil
      real*8 x(5000),y(5000),sigs(5000),jitter(NDSMAX)
      integer i,j,npl,ndset,ndata,idset,mfit,ma,idsmax(NDSMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),vxh(NPLMAX),vyh(NPLMAX)
     &       ,vzh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),vxj(NPLMAX),vyj(NPLMAX)
     &       ,vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX)
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*4 model_max
      parameter (AU=1.49597892d11, day = 86400.d0)


      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini
            
      nbod = npl+1
      rms = 0.d0
 
 
ccccccccccccccccccc t[JD], obs., cal., O-C   ccccccccccccc   

 
      call RVKEP_ewcop_fin (t,a,ymod,dyda,ma,ndata,ia,epsil,deltat)
      call MA_J_cop_fin (a,ma,npl,mstar,sini,mass,ap)

      do i = 1,ndata
          idset = ts(i)
 
	  ys(i) = ys(i) - a(7*npl+idset)
 
          if (writeflag_RV.gt.0) then
          write(*,*) t0 + t(i)/8.64d4 ,ymod(i),ys(i), 
     &                ys(i) - ymod(i),sigs(i),ts(i)

          endif
     
          rms = rms + (ys(i) - ymod(i))**2
      enddo
      rms = dsqrt(rms/dble(ndata))




      if(writeflag_best_par.gt.0) then

                write (*,*) 'Best-fit K [m/s], P [days], e, w [deg], 
     &          M0 [deg], i[deg], cap0m[deg] and their errors'
              do j = 1,npl
              i = 7*(j-1)


              write(*,*) a(i+1),a(i+2)/8.64d4,a(i+3),
     &                a(i+4)*180.d0/PI,a(i+5)*180.d0/PI,
     &                dmod(a(i+6)*180.d0/PI,180.d0),
     &                dmod(a(i+7)*180.d0/PI,360.d0)
              write(*,*) dsqrt(covar(i+1,i+1)),
     &                dsqrt(covar(i+2,i+2))/8.64d4,
c     &                2.d0*PI/a(i+2)**2*dsqrt(covar(i+2,i+2))/8.64d4,
     &                dsqrt(covar(i+3,i+3)),
     &                dsqrt(covar(i+4,i+4))*180.d0/PI,
     &                dsqrt(covar(i+5,i+5))*180.d0/PI,
     &                dmod(dsqrt(covar(i+6,i+6))*180.d0/PI,180.d0),
     &                dmod(dsqrt(covar(i+7,i+7))*180.d0/PI,360.d0)
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
              write (*,*) '0'
          enddo 
       
          write (*,*) 'linear trend  [m/s per day]:'
          write (*,*) a(7*npl + ndset + 1)
          write (*,*) dsqrt(covar(7*npl + ndset + 1,7*npl + ndset + 1))          
         
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
      
      endif
 

      if(writeflag_fit.gt.0) then 

          dt = (t_max+model_max - t0)/dble(nt - 1)
          do i = 1,nt
             x(i) = (i-1)*dt*8.64d4
          enddo
          call RVKEP_ewcop_fin (x,a,ymod,dyda,ma,nt,ia,epsil,deltat)
          do i = 1,nt
             write(*,*) t0 + x(i)/8.64d4,ymod(i)
          enddo

      endif

      return

      end      


c MRQMIN attempts to reduce the chi^2 of a fit by the Levenberg-Marquardt
c method. It uses COVSRT, GAUSSJ, and MRQCOF.
c
c From Numerical Recipes.

	subroutine MRQMIN (x,ts,y,sig,ndata,a,ia,ma,covar,alpha,nca,
     & 	                   chisq,funcs,alamda,loglik,jitt,epsil,deltat)

	implicit none
	integer ma,nca,ndata,ia(ma),MMAX,NDSMAX,ts(ndata)
	parameter (MMAX=200,NDSMAX=20)
        integer npl,ndset,idset,idsmax(NDSMAX)
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
CC	  alamda = 5000.d0
      call MRQCOF (x,ts,y,sig,ndata,a,ia,ma,alpha,beta,
     &                 nca,chisq,funcs,loglik,jitt,epsil,deltat)
   
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
     &               nca,chisq,funcs,loglik,jitt,epsil,deltat)
	
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
     &	                   chisq,funcs,loglik,jitt,epsil,deltat)

	implicit none
	integer npl,ndset,idset,ma,nalp,ndata,ia(ma),NDSMAX,MMAX
	parameter (NDSMAX=20, MMAX=200)
        integer idsmax(NDSMAX),ts(ndata)
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

        call FUNCS (x,a,ymod,dyda,ma,ndata,ia,epsil,deltat)
        
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
          dyda(i,7*npl + ndset + 1) = (x(i)/86400.d0)


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
     & deltat)
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
      integer ts(ndata),correct, idsmax(NDSMAX)

      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      nbod = npl + 1
      na = 7*npl


      correct = 1
      if(correct.gt.0) then      
      do i = 1,npl
         j = 7*(i-1)
         
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
         if (a(j+4).lt.0.d0) a(j+4) = dmod(a(j+4)+2.d0*PI,  2.d0*PI )  
         if (a(j+5).lt.0.d0) a(j+5) = dmod(a(j+5)+2.d0*PI,  2.d0*PI ) 
         if (a(j+4).gt.2.d0*PI) a(j+4) = dmod(a(j+4),  2.d0*PI )  
         if (a(j+5).gt.2.d0*PI) a(j+5) = dmod(a(j+5),  2.d0*PI )  
                
c         if (a(j+6).lt.0.0d0) a(j+6) = dmod(a(j+6) + PI,  PI )  
         if (a(j+7).lt.0.0d0) a(j+7) = dmod(a(j+7) + 2.d0*PI,  2.d0*PI ) 
         
c         if (a(j+6).ge.PI) a(j+6) = dmod(a(j+6)+ 0.001d0,  PI )  
         if (a(j+7).gt.2.d0*PI) a(j+7) = dmod(a(j+7),2.d0*PI)                         
      enddo  
      endif
 
c-----get ymod first


      call MA_J_cop_fin (a,ma,npl,mstar,sini,mass,ap)
           
      
      call GENINIT_J3_ewcop (nbod,ap,a,
     &                          mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)
      
      call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)
     
        do i = 1,ndata        !  initialize ymod 
           ymod(i) = 0.0
        enddo 	

       
      call integrate_cop_fin(ymod,t,nbod,na,ndata,mass,a,
     &  xh,yh,zh,vxh,vyh,vzh,epsil,deltat)

     
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

           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap)


           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                          mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)

           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)
     
           call integrate_cop_fin(ymodhb,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat)
    
             
c        get ah of ahead
    
           ah(i) = a(i) + ahh(i)
       
           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap)
           
           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                         mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)

           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)

           call integrate_cop_fin(ymodha,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat)
         
         
c        calculate the ith dyda 
           do j = 1,ndata
               dyda(j,i) = (ymodha(j) - ymodhb(j))/(2.d0*ahh(i))          

           enddo
           
         else
         
         
           ah(i) = a(i) - ahh(i)
c           ah(i) = dmod(a(i) - ahh(i),  0.0d0 ) 

 

c           if (ah(j+6).le.0.d0) ah(j+6) = dmod(ah(j+6) + PI,  PI )  
c           if (ah(j+6).ge.PI)   ah(j+6) = dmod(ah(j+6) +0.00001,  PI) 

                       
         
           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap)

           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                          mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)
     
           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)

           call integrate_cop_fin(ymodhb,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat)
             
c        get ah of ahead
    
           ah(i) = a(i) + ahh(i)
c           ah(i) = dmod(a(i) + ahh(i),  PI ) 


           call MA_J_cop_fin (ah,ma,npl,mstar,sini,mass,ap)
           
           call GENINIT_J3_ewcop (nbod,ap,ah,
     &                         mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)

           call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)

           call integrate_cop_fin(ymodha,t,nbod,na,ndata,mass,ah,
     &                  xh,yh,zh,vxh,vyh,vzh,epsil,deltat)
         
 
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

	subroutine MA_J_cop_fin (a,ma,npl,m0,sini,mass,ap)
        
	implicit none
	real*8 m0,sini,PI,TWOPI,THIRD,GMSUN,dm,MSUN
        integer npl,ma,i,j,NPLMAX
        parameter (NPLMAX=20)
        real*8 a(ma),mass(NPLMAX),ap(NPLMAX),mpold(NPLMAX),mtotal
	parameter (THIRD=1.d0/3.d0)
        parameter (PI=3.14159265358979d0,TWOPI=2.d0*PI)
	parameter (GMSUN=1.32712497d20,MSUN=1.32712497d20)

c*******G is set to be unit, and s, m, kg as unit of time, length and mass
c*******expectively.        
                                          
 
        do i = 0,npl-1

           mass(1) = m0
	   mpold(i+1) = 0.d0
 101       continue
           if (i.eq.0) then
           mtotal = m0
	   mass(i+2) = abs(a(7*i+1))*(a(7*i+2)*(m0 + mpold(i+1))**2/
     &               (TWOPI*GMSUN))**THIRD*dsqrt(1.d0-a(7*i+3)**2)
     &               /dabs(dsin(a(7*i+6)))
           else
              mtotal = m0
              do j = 0, i-1
                 mtotal = mtotal + mass(j+2)
              enddo
              mass(i+2) = abs(a(7*i+1))*(a(7*i+2)*(mtotal
     &                  +mpold(i+1))**2/(TWOPI*GMSUN))**THIRD
     &                  *dsqrt(1.d0-a(7*i+3)**2)/dabs(dsin(a(7*i+6)))
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
     &                       mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)

c      include 'swift.inc'

      real*8 SMASSYR,MSUN,PI,eps
      parameter (PI=3.14159265358979d0,eps=1.d-7)
      parameter (SMASSYR=4.d0*PI*PI)
      parameter (MSUN=1.32712497d20)

      integer nbod,NPLMAX,i,j
      parameter (NPLMAX=20)
      real*8 mass(NPLMAX)
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
      

c      pause 
      
      do i = 2,nbod
         j = 7*(i-2)


        gm = gm + mass(i)
        rpl(i) = frho3*(1.5d0*mass(i)/2.d0*PI)**0.3333333333d0
c        rhill(i) = ap(i-1)*(mass(i)/(3.d0*mass(1)))**0.3333333333d0
        rhill(i) = ap(i-1)*(1-a(j+3))*(mass(i)/(3.d0*mass(1)))**0.3333333333d0        


c        if (a(j+6).lt.0.d0) a(j+6) = dmod(a(j+6)+ PI,  PI )  
c        if (a(j+7).lt.0.d0) a(j+7) = dmod(a(j+7)+2.d0*PI,  2.d0*PI ) 
c        if (a(j+6).gt.PI) a(j+6) = dmod(a(j+6),  PI )  
c        if (a(j+7).gt.2.d0*PI) a(j+7) = dmod(a(j+7),  2.d0*PI ) 


c        if (a(j+6).lt.0.d0) a(j+6) = dmod(a(j+6) + PI,  PI )  
c        if (a(j+7).lt.0.d0) a(j+7) = dmod(a(j+7) + 2.d0*PI,  2.d0*PI ) 
         
c        if (a(j+6).gt.PI) a(j+6) = PI - dmod(a(j+6),  PI )  
c        if (a(j+7).gt.2.d0*PI) a(j+7) = dmod(a(j+7),2.d0*PI)   
        
        
               
        call ORBEL_EL2XV (gm,ialpha,ap(i-1),a(j+3),a(j+6),a(j+7)     
     &          ,a(j+4),a(j+5),xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i))
     
     
c       pause 


 
      enddo
c      pause
      return
      end




      subroutine integrate_cop_fin(ymod,t,nbod,ma,ndata,mass,a,
     &                     xh,yh,zh,vxh,vyh,vzh,epsil,deltat)

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
	integer iflgchk,iub,iuj,iud,iue
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
      call coord_h2b_tp(ntp,xht,yht,zht,vxht,vyht,vzht,
     &     xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &     xbt,ybt,zbt,vxbt,vybt,vzbt)

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
      do i=1,ntp
         if(istat(i,1).eq.0) then
            ntpi = ntpi + 1
            j = ntpi + nbod
            link(ntpi) = i
            ybs(1,j) = xbt(i)
            ybs(2,j) = ybt(i)
            ybs(3,j) = zbt(i)
            ybs(4,j) = vxbt(i)
            ybs(5,j) = vybt(i)
            ybs(6,j) = vzbt(i)
            do jj = 1,NSTAT
               istattmp(ntpi,jj) = istat(i,jj)
            enddo
         endif
      enddo

      tfake = 0.0d0
      dttmp = dt

c      do while(tfake.lt.dt)
      do while( (abs(tfake-dt)/dt) .gt. 1.0e-7 )    ! just to be real safe
         call bs_int(nbod,ntpi,mass,j2rp2,j4rp4,istattmp,
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

      do i=1,ntpi
         j = i + nbod
         xbt(link(i)) = ybs(1,j)
         ybt(link(i)) = ybs(2,j)
         zbt(link(i)) = ybs(3,j)
         vxbt(link(i)) = ybs(4,j)
         vybt(link(i)) = ybs(5,j)
         vzbt(link(i)) = ybs(6,j)
         do jj = 1,NSTAT
            istat(link(i),jj) = istattmp(i,jj)
         enddo
      enddo

c...  Convert back to helio. coords at the end of the step
	call coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &         xh,yh,zh,vxh,vyh,vzh)
	call coord_b2h_tp(ntp,xbt,ybt,zbt,vxbt,vybt,vzbt,
     &         xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &         xht,yht,zht,vxht,vyht,vzht)

      return

      end   ! bs_step
c------------------------------------------------------------------------



