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
      integer idsmax(NDSMAX),ia(MMAX),nt, ts(5000)
      real*8 x(5000),y(5000),sig(5000),ymod2(5000),y_in(5000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 rms,mstar,mass(NPLMAX),ap(NPLMAX)
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*8 chisq,alamda,ochisq,dchisq, epsil, deltat
      real*8 jitter,sigscale,x0,xmax,loglik
      real*8 t0,t1,t2,dt,offset,t_max, incl(NPLMAX), cap0m(NPLMAX)
      real*8 st_mass,sini,m1,a1,m2,a2,jitt(NDSMAX),epoch
      real*8 ymod(5000),dyda(MMAX)
      real*4 t_stop,when_to_kill,model_max
      
      external rvkep
      character*80 infile

      common /DSBLK/ npl,ndset,idsmax,idset

      loglik=0.d0
      
c     these two just for consistency with dynamical input and amoebastarts for consistency with loglik, not really used
      read (*,*) epsil, deltat,a moebastarts, 
     &          when_to_kill, nt, model_max
      
c      write(*,*) 'Stellar mass, sini'
      read (*,*) st_mass,
     &          writeflag_best_par, 
     &	             writeflag_RV,writeflag_fit 

c      write (*,*) ' Stellar Jitter: '
c      read (*,*) jitter


      call io_read_data (ndata,x,ts,y,sig,jitt,epoch,
     &               x0,t_max,a,ia,ma,incl,cap0m)
      
      mfit = 0
      do j = 1,ma
          if (ia(j).ne.0) mfit = mfit + 1
      enddo
 


      alamda = -1.d0
      call MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,MMAX,
     & 	           chisq,rvkep,alamda,loglik,jitt)

      i = 0
 500  continue
          i = i + 1
          ochisq = chisq
          call MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,MMAX,
     &                 chisq,rvkep,alamda,loglik,jitt)

          dchisq = chisq - ochisq
          if (i.eq.10) then
              i = 0
CC              pause
          endif
          
          CALL SECOND(t_stop)
          if (t_stop.ge.when_to_kill) then
            write(*,*) 'Max. time=',when_to_kill, 'sec ', 
     &                 'exceeded t_stop =', t_stop, 'sec ' 
            goto 502
           endif
          
          
      if ((chisq.ge.ochisq).or.(dchisq.lt.-1.d-2)) goto 500

502   alamda = 0.d0
      call MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,MMAX,
     &                 chisq,rvkep,alamda,loglik,jitt)

      idset = 1
      rms = 0.d0
c     y_in = y
  
 
         do i = 1,ndata
              idset = ts(i)
c             if (i.gt.idsmax(idset)) idset = idset + 1
              call RVKEP (x(i),a,ymod(i),dyda,ma,ts(i))
 
	      xmax = x0 + x(i)

              y_in(i) = y(i) - a(5*npl+idset) 
 
 
              rms = rms + (y(i) - ymod(i))**2
      if (writeflag_RV.gt.0) then 
              write(*,*) x0 + x(i),
     &                   ymod(i) - a(5*npl+idset) - 
     &                   a(5*npl+ndset+1)*x(i),
     &                   y_in(i) - a(5*npl+ndset+1)*x(i),
     &                   y(i) - ymod(i), sig(i), idset
   
        
      endif
         enddo

      rms=dsqrt(rms/dble(ndata))
      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
      write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms


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
c              offset(j) = a(i)
              write (*,*) dsqrt(covar(i,i))
          enddo

          write (*,*) 'Jitters (not fitted in chi^2 method):'
          do j = 1,ndset
              write (*,*) jitt(j)
              write (*,*) '0'
          enddo          
          
          write (*,*) 'linear trend [m/s per day]:'
          write (*,*) a(5*npl + ndset + 1)  
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


c      nt = 5000
      if(writeflag_fit.gt.0) then
 
              dt = (xmax+model_max - x0)/dble(nt - 1)
      
              do i = 1,nt
	            x(i) = ((i-1)*dt)-0.00

                    do j = 1,ndset
                         a(5*npl + j) = 0.0
                    enddo	  
	  
                    call RVKEP (x(i),a,ymod(i),dyda,ma,1)
 
                    write(*,*) x0 + x(i), ymod(i)
             enddo
 
      endif
      
      
      
      

c      stop
      end
      
      


c*************************************************************************      
c**********************       read RV data      **************************
c*************************************************************************

C     Xianyu Tan 2011
                                                                            
      subroutine io_read_data (ndata,t,ts,ys,sigs,jitt,epoch,t0,t_max,
     &          ar,iar,ma,incl,cap0m)  

      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl,ma
      real*8 t(5000),y(5000),sig(5000),ys(5000),sigs(5000),PI
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      parameter(PI=3.14159265358979d0)
      real*8 ar(MMAX),incl(NPLMAX),cap0m(NPLMAX)
      integer iar(MMAX),u_off(NDSMAX),u_jit
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
          ar(5*npl+i)=off(i)
          iar(5*npl+i)=u_off(i)
      enddo
      
      ma = 5*npl + ndset + 1
c      write (*,*) 'Initial K, P, e, w, M0 and their flags: '
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
c u_jit read for consistency with input/output in loglik case, but here we don't actually save this information and not use it, jitters cannot be used for fit in chi^2 minimization
          
      read (*,*) ar(5*npl+ ndset+1)
      read (*,*) iar(5*npl+ndset+1)    

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
      real*8 x,y,a(ma),dyda(ma)
      real*8 cosw,sinw,capm,cape,cose,sine,cosf,sinf,fac1,fac2,fac3
      real*8 orbel_ehybrid, f, coswf

      common /DSBLK/ npl,ndset,idsmax,idset
      
      y = 0.d0
 
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
         if (a(j+4).lt.0.d0) a(j+4) = dmod(a(j+4)+2.d0*PI,  2.d0*PI )  
         if (a(j+5).lt.0.d0) a(j+5) = dmod(a(j+5)+2.d0*PI,  2.d0*PI ) 
         if (a(j+4).gt.2.d0*PI) a(j+4) = dmod(a(j+4),  2.d0*PI )  
         if (a(j+5).gt.2.d0*PI) a(j+5) = dmod(a(j+5),  2.d0*PI )         
         
                        
      enddo  

      do j = 1,npl

          i = 5*(j-1)
          cosw = dcos(a(4+i))
          sinw = dsin(a(4+i))

          capm = TWOPI*x/a(2+i) + a(5+i)
          capm = dmod(capm,  2.d0*PI )

          cape = ORBEL_EHYBRID (a(3+i),capm)
          cose = dcos(cape)
          sine = dsin(cape)
          
          cosf = (cose - a(3+i))/(1.d0 - a(3+i)*cose)
          sinf = (dsqrt(1.d0 - a(3+i)**2)*sine)/(1.d0 - a(3+i)*cose)

          fac1 = cosw*cosf - sinw*sinf + a(3+i)*cosw

          fac2 = (cosw*sinf + sinw*cosf)/(1.d0 - a(3+i)*cose)**2
          fac3 = -a(1+i)*dsqrt(1.d0 - a(3+i)**2)*fac2

          y = y + a(1+i)*fac1
          dyda(1+i) = fac1
          dyda(2+i) = -TWOPI*fac3*x/a(2+i)**2
          dyda(3+i) = -a(1+i)*sine*(2.d0 - a(3+i)**2 - a(3+i)*cose)*fac2/dsqrt(1.d0 - a(3+i)**2)
          dyda(4+i) = -a(1+i)*(sinw*cosf + cosw*sinf)
          dyda(5+i) = fac3

      enddo

c      do i = 1,idset
      y = y + a(5*npl+ts)

      do i = 1,ts-1
          dyda(5*npl+i) = 0.d0
      enddo
      
      dyda(5*npl+ts) = 1.d0

      y = y + a(5*npl +ndset + 1)*x  
      dyda(5*npl + ndset + 1) = x
   
      do i = ts+1,ndset
          dyda(5*npl+i) = 0.d0
      enddo

      return
      end

c MRQMIN attempts to reduce the chi^2 of a fit by the Levenberg-Marquardt
c method. It uses COVSRT, GAUSSJ, and MRQCOF.
c
c From Numerical Recipes.

	subroutine MRQMIN (x,y,sig,ndata,a,ia,ma,ts,covar,alpha,nca,
     & 	                   chisq,funcs,alamda,loglik,jitt)
	implicit none
	integer ma,nca,ndata,ia(ma),MMAX,NDSMAX,ts(5000)
	real*8 alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     &	       sig(ndata),x(ndata),y(ndata),loglik
	external funcs
	parameter (MMAX=200,NDSMAX=20)
	integer j,k,l,mfit
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
     &	       alpha,beta,nca,chisq,funcs,loglik,jitt)

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
     &	       ts,covar,da,nca,chisq,funcs,loglik,jitt)

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
     &	                   chisq,funcs,loglik,jitt)

	implicit none
	integer npl,ndset,idset,ma,nalp,ndata,ia(ma),NDSMAX,MMAX
	parameter (NDSMAX=20, MMAX=200)
        integer idsmax(NDSMAX),ts(5000)
	real*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),
     &	       x(ndata),y(ndata),loglik,TWOPI,jitt(NDSMAX)
        parameter (TWOPI=2.d0*3.14159265358979d0)
	external funcs
	integer mfit,i,j,k,l,m
	real*8 dy,sig2i,wt,ymod,dyda(MMAX)

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
        idset = 1
	do i = 1,ndata
cc          if (i.gt.idsmax(idset)) idset = idset + 1
          idset = ts(i)
 
	  call FUNCS (x(i),a,ymod,dyda,ma,idset)
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


