
c*************************************************************************
C****************   Radial Velocity dynamical model      *****************
c*************************************************************************
 

      implicit none
      real*8 PI
      parameter (PI=3.14159265358979d0)
      integer npl,ndset,idset,ndata,ma,mfit,i,j,NDSMAX,NPLMAX,MMAX,k
      real*8 mstar,sini(7),sig2i,dy, loglik,sig2l,ftol
      parameter (NDSMAX=20, NPLMAX=20, MMAX=200)
      integer idsmax(NDSMAX),ia(MMAX),ts(20000) ,nt, iter, ii
      integer writeflag_best_par,hkl,gr_flag
      integer writeflag_RV,writeflag_fit, amoebastarts
      integer dynamical_planets(NPLMAX)
      real*8 t(20000),x(20000),y(20000),sig(20000),ys(20000),sigs(20000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 chisq,alamda,ochisq,dchisq,red_chisq, p(MMAX+1,MMAX),
     &       yamoeba(MMAX+1), loglikk, ologlikk, dloglikk
      real*8 sigscale,t0,t_max,twopi,dt, epoch, epsil,deltat
      real*8 rms,ymod(20000),dyda(20000,MMAX),jitt(NDSMAX)
      real*4 t_stop, when_to_kill,model_max,model_min
      real*8 wdot(NPLMAX),u_wdot(NPLMAX)
      
      
      external rvdyn, compute_abs_loglik
      character*80 infile
      character*80 version_input, version
      
      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      version = "0.06"
       
      CALL getarg(1, version_input)     
      if(version_input.eq.'-version') then
          write(*,*) version
          goto 222
      endif

      ftol=0.0015d0       
      read (*,*) epsil,deltat, amoebastarts,
     &          when_to_kill, nt, model_max,model_min,gr_flag      
      
      read (*,*) mstar,
     &          writeflag_best_par, 
     &	             writeflag_RV,writeflag_fit 
    
      call io_read_data (ndata,t,ts,ys,sigs,jitt,
     & 	           epoch,t0,t_max,a,ia,ma,mfit,dynamical_planets,
     & 	           hkl,wdot,u_wdot)

           i = 0
c      call timer(t_start)              
 500  continue
        if (i.eq.amoebastarts) then  
             i = 0
             goto 502
         endif
         i = i + 1
         ologlikk = loglikk

      
         call prepare_for_amoeba(p,MMAX+1,MMAX,yamoeba,a,ia,ma,mfit,
     & compute_abs_loglik,ndata,t,ys,ymod,dyda,ts,sigs,epsil,deltat,i,
     & dynamical_planets)
     
         call amoeba(p,yamoeba,MMAX+1,MMAX,mfit,ftol,compute_abs_loglik,
     & iter,ndata,t,ys,ymod,dyda,ma,ts,sigs,a,ia,epsil,deltat, 
     & dynamical_planets)
     
         CALL SECOND(t_stop)
         if (t_stop.ge.when_to_kill) then
            write(*,*) 'Max. time=',when_to_kill, 'sec ', 
     &                 'exceeded t_stop =', t_stop, 'sec ' 
            goto 502
         endif     
     
         loglikk = yamoeba(1)
 
         dloglikk = ologlikk - loglikk      

c         write (*,*) i, dloglikk, loglikk
         j=0
         do ii=1,ma
           if (ia(ii).ne.0) then
              j=j+1
              a(ii)=p(1,j)  
           endif
         enddo
c      call timer(t_stop)
c      CALL SECOND(t_stop)
c      if (t_stop.ge.10) then
c            write(*,*) t_stop 
c            goto 502
c      endif
         
      if (dabs(dloglikk).ge.0.000001d0) goto 500

502   j=0

      loglik = 0.0d0
      chisq = 0.0d0
      rms = 0.0d0
      
      call io_write_bestfitpa_ewcop_fin (a,covar,t,ys,ndata,ts,
     & 	           ma,mfit,t0,t_max,sigs,chisq,rms,loglik,writeflag_RV,
     &             writeflag_best_par,writeflag_fit,jitt,epsil,deltat,
     &             nt, model_max, model_min, dynamical_planets,
     &             hkl,wdot,u_wdot)


c      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
c      write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms
        

c      stop
222   end


      subroutine compute_abs_loglik(ndata,x,y,a2,ymod,dyda,ma,mfit,ts,
     & sig,loglik, num,a,ia,epsil,deltat,dynamical_planets)
      implicit none
      
      integer MMAX,npl,ndset,NDSMAX,NPLMAX,idset,num, mfit    
      parameter (MMAX=200, NDSMAX=20,NPLMAX=20)
      real*8 twopi, loglik
      parameter (twopi=6.28318530717958d0)
      integer ndata, i, j, ma, ts(20000), ia(MMAX), idsmax(NDSMAX)
      real*8 dy, sig(20000), dyda(20000,MMAX), x(20000), y(20000)
      real*8 ymod(20000), a(MMAX), a2(mfit), a3(MMAX),sig2i, y_in(20000)
     & , y2(20000), epsil,deltat
      integer dynamical_planets(NPLMAX)
   
      common /DSBLK/ npl,ndset,idsmax,idset      
   
      j=0
      do i=1,ma 
         if (ia(i).ne.0) then
             j=j+1
             a3(i)=a2(j)
         else
             a3(i)=a(i)
         endif
      enddo
      
      loglik=0.d0
        call RVDYN (x,a3,y2,dyda,ma,ndata,epsil,deltat, 
     &  dynamical_planets,ts)
c        write(*,*) a3(2),a2(2),a(2)
c        write(*,*) a3(1),a3(2),a3(3),a3(4),a3(5)
c        write(*,*) a3(8),a3(9),a3(10),a3(11),a3(12)
c        write(*,*) y_in(1),y_in(2)

      do i = 1,ndata
              idset = ts(i)
              y_in(i) = y(i) - a3(7*npl+idset) 
     &                 -a3(7*npl+2*ndset+1)*(x(i)/86400.d0) 
     &                 -a3(7*npl+2*ndset+2)*(x(i)/86400.d0)**2     
c	      y2(i) = y2(i) - a3(7*npl+idset) - 
c     &        a3(7*npl +2*ndset + 1)*x(i)

          dy = y_in(i) - y2(i)


          sig2i = 1.d0/(sig(i)**2 + a3(7*npl+ndset+idset)**2)
 
       
 

	  loglik =  loglik + 0.5*dy*dy*sig2i +
     &               dlog(dsqrt(twopi*(sig(i)**2 + 
     &               a3(7*npl+ndset+idset)**2))) 
     &               - dlog(dsqrt(twopi)) 


        enddo

c      write(*,*) "lnL", loglik
             
      return
      end            
      
      subroutine prepare_for_amoeba(p,mp,np,y,a,ia,ma,mfit,funk,ndata,
     & x,z,ymod,dyda,ts,sig,epsil,deltat,it,dynamical_planets)
      integer MMAX, NDSMAX, NPLMAX,ma,ts(20000), ndata,mp,np,mfit,it
      parameter(MMAX=200, NDSMAX=20, NPLMAX=20)
      REAL*8 ftol,p(mp,np),y(mp),a(MMAX), a2(mfit),fr,frjitt
      real*8 x(20000),z(20000),ymod(20000)
      real*8 dyda(MMAX), sig(20000), loglik,epsil,deltat
      parameter(fr=0.01, frjitt=0.05)
      integer i,j,k,ia(MMAX),idsmax(NDSMAX),dynamical_planets(NPLMAX)
      external funk
    
      common /DSBLK/ npl,ndset,idsmax,idset
      
      k=0
      do j=1,ma
          if(ia(j).ne.0) then
          k=k+1
          p(1,k)=a(j)
          do i=2,mfit+1
              if (k.eq.(i-1)) then
                  if (j.gt.(7*npl+ndset)) then
                  p(i,k)=(1+frjitt)*(p(1,k)+0.1)
                  else 
                  if (mod(j,7).eq.2) then
                     p(i,k)=(1+fr)*(p(1,k) + 0.1)
                  else if (mod(j,7).eq.3) then
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
c              write(*,*) a2(j)  
          enddo

          call funk(ndata,x,z,a2,ymod,dyda,ma,mfit,ts,sig,loglik,i,
     &              a,ia,epsil,deltat, dynamical_planets)
          y(i)=loglik
	  
      enddo
c      write(*,*) loglik 
      return
      end
      
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,epsil,deltat,dynamical_planets)
      implicit none
      INTEGER iter,mp,ndim,np,NMAX,ITMAX, MMAX,ma,ts(20000),ndata,NPLMAX
      REAL*8 ftol,p(mp,np),y(mp),x(20000),z(20000),ymod(20000)
      PARAMETER (NMAX=20,ITMAX=50000,MMAX=200,NPLMAX=20)
      real*8 dyda(20000,MMAX), sig(20000), loglik, a(MMAX)
      EXTERNAL funk
      INTEGER i,ihi,ilo,inhi,j,m,n, ia(MMAX), dynamical_planets(NPLMAX)
      REAL*8 rtol,summ,swap,ysave,ytry,psum(ndim),amotry,epsil,deltat
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
c      write (*,*) rtol,ftol
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
          return
      endif
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,epsil,deltat, dynamical_planets)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,epsil,deltat,dynamical_planets)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,epsil,deltat,dynamical_planets)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              call funk(ndata,x,z,psum,ymod,dyda,ma,ndim,ts,sig,loglik,
     &                 i,a,ia,epsil,deltat,dynamical_planets)
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
     & dyda,ma,ts,sig,a,ia,epsil,deltat,dynamical_planets)
      implicit none
      INTEGER ihi,mp,ndim,np,NMAX, MMAX, NPLMAX, ma, ts(20000),ndata
      PARAMETER (NMAX=20, MMAX=200, NPLMAX=20)
      REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),x(20000),z(20000),
     & ymod(20000), epsil, deltat
      real*8 dyda(20000,MMAX), sig(20000),loglik
      EXTERNAL funk
      INTEGER j, ia(MMAX), dynamical_planets(NPLMAX)
      REAL*8 fac1,fac2,ytry,ptry(ndim), a(MMAX)
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      call funk(ndata,x,z,ptry,ymod,dyda,ma,ndim,ts,sig,loglik,ihi,
     & a,ia,epsil,deltat,dynamical_planets)   
      ytry=loglik
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
c      write(*,*) amotry
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0=M,173+9.      
      
c*************************************************************************      
c**********************       read RV data      **************************
c*************************************************************************

 
                                                                            
      subroutine io_read_data (ndata,t,ts,ys,sigs,jitter,epoch,
     &               t0,t_max,ar,iar,ma,mfit, dynamical_planets,
     & 	           hkl,wdot,u_wdot)  


      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl,hkl
      real*8 t(20000),y(20000),sig(20000),ys(20000),sigs(20000)
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      integer idsmax(NDSMAX),ts(20000), dynamical_planets(NPLMAX)
      real*8 jitter(NDSMAX),t0,t_max,epoch,ar(MMAX),off(NDSMAX), PI
      parameter(PI=3.14159265358979d0)
      integer i,k,j, iar(MMAX), u_off(NDSMAX), u_jit(NDSMAX), ma, mfit
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
          read (*,*) u_jit(i)
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
      
          do i = 1,npl
            read (*,*) dynamical_planets(i)
          enddo  
      
      do i = 1,ndset
          ar(7*npl+i)=off(i)
          iar(7*npl+i)=u_off(i)
          ar(7*npl+ndset+i)=jitter(i)
          iar(7*npl+ndset+i)=u_jit(i)
          enddo
      
c           write(*,*) ' Number of Planets: '
      if (npl.gt.NPLMAX) stop ' KEPFIT: npl > NPLMAX.'

      ma = 7*npl + 2*ndset + 2
      
c      write(*,*) 'Initial K, P, e, w, M0,Inc,Capom and their flags: '
      do j = 1,npl
          i = 7*(j-1)
          read (*,*) ar(i+1),ar(i+2),ar(i+3),ar(i+4),ar(i+5),ar(i+6),
     &               ar(i+7),wdot(i)
          read (*,*) iar(i+1),iar(i+2),iar(i+3),iar(i+4),iar(i+5),
     &               iar(i+6),iar(i+7),u_wdot(i)
 
 
c          ar(i+2) = 2.d0*PI/(ar(i+2)*8.64d4)         ! mean motion 
c          ar(i+2) = ar(i+2)*8.64d4         ! second as unit
c          ar(i+4) = ar(i+4)*PI/180.d0                ! radians
c          ar(i+5) = ar(i+5)*PI/180.d0
c          ar(i+6) = ar(i+6)*PI/180.d0
c          ar(i+7) = ar(i+7)*PI/180.d0
          u_wdot(i) = 0  
      enddo
  
c      write (*,*) 'linear trend:'      
      read (*,*) ar(7*npl + 2*ndset + 1)
      read (*,*) iar(7*npl + 2*ndset + 1)   
      
      read (*,*) ar(7*npl + 2*ndset + 2)
      read (*,*) iar(7*npl + 2*ndset + 2)  
      
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
C##########################################################################
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C##########################################################################
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C**************************************************************************
C**********   output best-fit parameters and errorbars    *****************
C**************************************************************************
C 
 

       subroutine io_write_bestfitpa_ewcop_fin (a,covar,t,ys,ndata,ts,
     &           ma,mfit,t0,t_max,sigs,chisq,rms,loglik,writeflag_RV,
     &           writeflag_best_par,writeflag_fit,jitter,epsil,deltat,
     &  nt, model_max,model_min, dynamical_planets,hkl,wdot,u_wdot)   
    
      implicit none 
      real*8 PI
      integer MMAX,NDSMAX,NPLMAX 
      parameter (PI=3.14159265358979d0,MMAX=200  ,NDSMAX=20,NPLMAX=20)
      real*8 a(MMAX),ia(MMAX),t(20000),ymod(20000),ys(20000)
      real*8 covar(MMAX,MMAX),dyda(20000,MMAX),AU,day
      real*8 rms,mstar,sini(7),mass(NPLMAX),ap(NPLMAX)
      integer ts(20000),nbod,nt,writeflag_RV,hkl,
     &  writeflag_best_par,writeflag_fit, dynamical_planets(NPLMAX)
      real*8 t0,t1,t2,dt,offset,t_max,chisq,loglik,dy,sig2i,twopi
      real*8 x(20000),y(20000),sigs(20000),jitter(NDSMAX)
      integer i,j,npl,ndset,ndata,idset,mfit,ma,idsmax(NDSMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),vxh(NPLMAX),vyh(NPLMAX)
     &       ,vzh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),vxj(NPLMAX),vyj(NPLMAX)
     &       ,vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX),epsil,deltat
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*4 model_max,model_min,best_w,best_we
      parameter (AU=1.49597892d11, day = 86400.d0)
      real*8 wdot(NPLMAX),u_wdot(NPLMAX)      
 
      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini
            
      nbod = npl+1
      rms = 0.d0
      twopi = 2.0d0*PI 
 
ccccccccccccccccccc t[JD], obs., cal., O-C   ccccccccccccc   

      call RVDYN(t,a,ymod,dyda,ma,ndata,epsil,deltat,
     & dynamical_planets,ts)   
     
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
 
	  ys(i) = ys(i) - a(7*npl+idset) - 
     &               a(7*npl  + 2*ndset + 1)*(t(i)/8.64d4)- 
     &               a(7*npl  + 2*ndset + 2)*(t(i)/8.64d4)**2

c          write(*,*) a(7*npl+idset), a(7*npl  + 2*ndset + 1)

          if (writeflag_RV.gt.0) then
          write(*,*) t0 + t(i)/8.64d4 ,ymod(i),ys(i) +
     &                a(7*npl  + 2*ndset + 1)*(t(i)/8.64d4)+ 
     &                a(7*npl  + 2*ndset + 2)*(t(i)/8.64d4)**2, 
     &                ys(i) - ymod(i),sigs(i),ts(i)

          endif
     
          sig2i = 1.d0/(sigs(i)**2 + a(7*npl+ndset+idset)**2)

          dy =  ys(i) -ymod(i)  

 	  chisq  = chisq + dy*dy*sig2i
 
c          write(*,*) "TEST:",loglik,dy,ymod(i)

	  loglik =  loglik - 0.5*dy*dy*sig2i -
     &               0.5*dlog(twopi*(sigs(i)**2 + 
     &               a(7*npl+ndset+idset)**2)) 
c     &               + dlog(dsqrt(twopi))  

          rms = rms + dy**2
      enddo

      rms = dsqrt(rms/dble(ndata))

      write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
      write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms

      if(writeflag_best_par.gt.0) then



        write (*,*) 'Best-fit K [m/s], P [days], e, w [deg], 
     &                M0 [deg], i[deg], cap0m[deg] and their errors'
          do j = 1,npl
              i = 7*(j-1)
              
              if (hkl.eq.0) then
                  best_w = a(i+4)*180.d0/PI
                  best_we = dsqrt(covar(i+4,i+4))*180.d0/PI
              else    
                  best_w = a(i+4) 
                  best_we = dsqrt(covar(i+4,i+4)) 
              endif    
!         a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)
              write(*,*) a(i+1),
     &                a(i+2),a(i+3),best_w, a(i+5)*180.d0/PI,
     &                dmod(a(i+6)*180.d0/PI,180.d0),
     &                dmod(a(i+7)*180.d0/PI,360.d0),wdot(i)
              write(*,*) dsqrt(covar(i+1,i+1)),
     &                dsqrt(covar(i+2,i+2)),
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
      
          write (*,*) 'Jitters for each dataset:'
          do j = 1,ndset
              write (*,*) a(7*npl+ndset+j)
              write (*,*) '0'
          enddo 
       
          write (*,*) 'linear trend  [m/s per day]:'
          write (*,*) a(7*npl + 2*ndset + 1)
          write (*,*) dsqrt(covar(7*npl + 2*ndset + 1,
     &                7*npl + 2*ndset + 1))
     
          write (*,*) 'quad. trend  [m/s per day]:'
          write (*,*) a(7*npl + 2*ndset + 2)
          write (*,*) dsqrt(covar(7*npl + 2*ndset + 2,
     &                7*npl + 2*ndset + 2))     

          write(*,*) ' ndata =',ndata
          write(*,*) ' mfit =',mfit
          write(*,*) ' RMS =',rms
          write(*,*) ' Chi^2 =',chisq/dble(ndata-mfit)
          write(*,*) ' epoch = ', t0

          do i = 1,npl
             j = 7*(i-1)
         
             a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)   
             
          enddo        

          call MA_J (a,ma,npl,mstar,sini,mass,ap)
          
          call GENINIT_J3 (nbod,ap,a,
     &                         mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill) 
     
          do i = 1,npl
             j = 7*(i-1)
         
             a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)   
             
          enddo        




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
          dt = ((t_max- t0) + model_max )/dble(nt - 1)      

          do i = 1,nt

             x(i) = (i-1)*dt*8.64d4
          enddo
          call RVDYN (x,a,ymod,dyda,ma,nt,epsil,deltat,
     &    dynamical_planets,ts)
          do i = 1,nt
             write(*,*) t0 + x(i)/8.64d4, 
     &       ymod(i) + a(7*npl  + 2*ndset + 1)*(x(i)/8.64d4)+ 
     &                 a(7*npl  + 2*ndset + 2)*(x(i)/8.64d4)**2          
          enddo

       endif

 

       return
       end     

C##########################################################################
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       subroutine split_parameters(a,a_kep,a_dyn,dynamical_planets,k,d)
     
	implicit none
	integer MMAX, NPLMAX, NDSMAX
	parameter (MMAX=200, NPLMAX=20, NDSMAX=20)	
	real*8 a(MMAX), a_kep(MMAX), a_dyn(MMAX), PI
	parameter (PI=3.14159265358979d0)
	integer ia(MMAX), ia_kep(MMAX), ia_dyn(MMAX)
	integer dynamical_planets(NPLMAX), i, j,k,d
        integer npl,ndset,idset	
        integer idsmax(NDSMAX)       

       common /DSBLK/ npl,ndset,idsmax,idset
	
	k=0
	d=0
c	write (*,*) 'HELLO'
	do i=1,npl
c           write (*,*) dynamical_planets(i)
	   if (dynamical_planets(i).eq.1.d0) then
                do j=1,7
                    a_dyn(7*d+j)=a(7*(i-1)+j)

                enddo
                d=d+1
            else
                do j=1,7
                    a_kep(7*k+j)=a(7*(i-1)+j)
c                    write(*,*) a(7*(i-1)+j)                    
                enddo
                k=k+1           
            endif
        enddo
           
       do i = 1,d
         j = 7*(i-1)
         a_dyn(j+2) = 2.d0*PI/(a_dyn(j+2)*8.64d4)          
       enddo
           
      end

      subroutine RVKEP (x,a,y,dyda,ma,ts,k)
      implicit none
      real*8 PI,TWOPI
      parameter (PI=3.14159265358979d0)
      parameter (TWOPI=2.0d0*PI)
      integer npl,ndset,idset,ma,i,j,NDSMAX,ts,k
      parameter (NDSMAX=20)
      integer idsmax(NDSMAX)
      real*8 x,y,a(ma),dyda(ma)
      real*8 cosw,sinw,capm,cape,cose,sine,cosf,sinf,fac1,fac2,fac3
      real*8 orbel_ehybrid, f, coswf

       common /DSBLK/ npl,ndset,idsmax,idset
 
 
       do j = 1,k

          i = 7*(j-1)
          cosw = dcos(a(4+i))
          sinw = dsin(a(4+i))

          capm = TWOPI*x/(a(2+i)*86400.d0) + a(5+i)
          capm = dmod(capm,  2.d0*PI )

          cape = ORBEL_EHYBRID (a(3+i),capm)
          cose = dcos(cape)
          sine = dsin(cape)
          
          cosf = (cose - a(3+i))/(1.d0 - a(3+i)*cose)
          sinf = (dsqrt(1.d0 - a(3+i)**2)*sine)/(1.d0 - a(3+i)*cose)

          fac1 = cosw*cosf - sinw*sinf + a(3+i)*cosw

          fac2 = (cosw*sinf + sinw*cosf)/(1.d0 - a(3+i)*cose)**2
          fac3 = -a(1+i)*dsqrt(1.d0 - a(3+i)**2)*fac2

          y = y+a(1+i)*fac1
          dyda(1+i) = fac1
          dyda(2+i) = -TWOPI*fac3*x/(a(2+i)*86400.d0)**2
          dyda(3+i) = -a(1+i)*sine*(2.d0 - a(3+i)**2 - 
     &       a(3+i)*cose)*fac2/dsqrt(1.d0 - a(3+i)**2)
          dyda(4+i) = -a(1+i)*(sinw*cosf + cosw*sinf)
          dyda(5+i) = fac3

      enddo

      do i = 1,ts-1
          dyda(7*npl+i) = 0.d0
      enddo
      
      dyda(7*npl+ts) = 1.d0

      dyda(7*npl + ndset + 1) = x
      dyda(7*npl + ndset + 2) = x*x
   
      do i = ts+1,ndset
          dyda(7*npl+i) = 0.d0
      enddo
      return
      end 

      subroutine RVDYN (t,a,ymod,dyda,ma,ndata,epsil,deltat,
     & dynamical_planets,ts)  
      
      implicit none
      real*8 PI,TWOPI,GMSUN,AU
      parameter (GMSUN=1.32712497d26,AU=1.49597892d11)
      parameter (PI=3.14159265358979d0)
      parameter (TWOPI=2.0d0*PI)
      integer npl,nbod,ndata,ma,i,j,NPLMAX,na,ndset,NDSMAX,idset
      parameter (NPLMAX=20, NDSMAX=20)
      real*8 a_kep(ma), a_dyn(ma)
      real*8 t(ndata),ymod(ndata),a(ma),a2(ma),dyda(ndata,ma),x(ma)
      real*8 ymod_kep(ndata)
      real*8 mstar,sini(NPLMAX),ap(NPLMAX),mass(NPLMAX), dyda_kep(ma)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),vxh(NPLMAX),vyh(NPLMAX)
     &       ,vzh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),vxj(NPLMAX),vyj(NPLMAX)
     &       ,vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX),epsil,deltat

      integer correct, idsmax(NDSMAX),
     &         dynamical_planets(NPLMAX), k, d, ts(20000)
     
      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      do i = 1,ma
          a2(i)=a(i)
      enddo      

      correct = 1
      if(correct.gt.0) then      
      do i = 1,npl
         j = 7*(i-1)

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
         if (a2(j+4).lt.0.d0) a2(j+4) = dmod(a2(j+4)+2.d0*PI,  2.d0*PI )  
         if (a2(j+5).lt.0.d0) a(j+5) = dmod(a2(j+5)+2.d0*PI,  2.d0*PI ) 
         if (a2(j+4).gt.2.d0*PI) a2(j+4) = dmod(a2(j+4),  2.d0*PI )  
         if (a2(j+5).gt.2.d0*PI) a2(j+5) = dmod(a2(j+5),  2.d0*PI )  
                
         if (a2(j+6).lt.0.0d0) a2(j+6) = dmod(a2(j+6) + PI,  PI )  
         if (a2(j+7).lt.0.0d0) a2(j+7)=dmod(a2(j+7)+2.d0*PI,2.d0*PI) 
         
         if (a2(j+6).ge.PI) a2(j+6) = dmod(a2(j+6)+ 0.001d0,  PI )  
         if (a2(j+7).gt.2.d0*PI) a2(j+7) = dmod(a2(j+7),2.d0*PI)                         
      enddo  
      endif      
      
      k=0
      d=0
      call split_parameters(a2,a_kep,a_dyn,dynamical_planets,k,d)
      nbod = d + 1
      na = 7*d 
          do i = 1,ndata        
             ymod(i)=0.d0
             ymod_kep(i)=0.d0
          enddo
 
c-----get ymod first
c        if (a_kep(1).gt.1) then
       do i = 1,ndata        !  initialize ymod 
          call RVKEP(t(i),a_kep,ymod_kep(i),dyda_kep,ma-na,ts(i),k)
 
       enddo 	
c        endif

        if (d.eq.0.d0) then
         go to 888
        endif
c******calculating masses and radius of planets in Jacobian system********
       call MA_J (a_dyn,na,d,mstar,sini,mass,ap)
           
           
c******getting Xs and Vs in Jacobian system**********
       call GENINIT_J3 (nbod,ap,a_dyn,
     &                         mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)
     
c       write(*,*) a2(1),a2(2),a2(3),a2(4),a2(5)
c       write(*,*) a2(8),a2(9),a2(10),a2(11),a2(12)
c       write(*,*) a2(15),a2(16),a2(17)

c       write(*,*) mass(1),mass(2),mass(3)
c       write(*,*) ap(1),ap(2)
c       write(*,*)  xj(1),yj(1),zj(1),vxj(1),vyj(1),vzj(1)
c       write(*,*)  xj(2),yj(2),zj(2),vxj(2),vyj(2),vzj(2)

c******transform to heliocentric coordinates, in order for using BS integrator**

       call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)
          do i = 1,npl
          j = 7*(i-1)
          x(j+1) = mass(i+1)
          x(j+2) = xh(i+1)
          x(j+3) = yh(i+1)
          x(j+4) = zh(i+1)
          x(j+5) = vxh(i+1)
          x(j+6) = vyh(i+1)
          x(j+7) = vzh(i+1) 
          enddo
 
C******calculating the RVs of model and the derivatives of orbel elements****** 
       call integrate(x,mass,ymod,dyda,t,a_dyn,ap,nbod,na,ndata,epsil,
     & deltat)
c       write(*,*) ymod
c       pause     
    
c       write(*,*) 'mass(2) , mass(3) , mass(4):'
c       write(*,*) mass(2)/1.2667d17,mass(3)/1.2667d17,mass(4)/1.2667d17
                               ! compare to Jupiter mass
c      write(*,*) mass(1),mass(2)
c      write(*,*) x(1),x(2)
c      write(*,*) na,ndata
c      write(*,*) ap(1),ap(2)
c      write(*,*) t(1),t(2)
c      write(*,*) ymod(1),ymod(2)

888      do i=1,ndata
         ymod(i)=ymod(i)+ymod_kep(i)
        enddo

      return
      end

          

c MA_J calculates the actual masses and Jacobi semimajor axes of a two-planet
c system for assumed sin(i) using the parameters P, K and e from a
c two-Kepler fit.

	subroutine MA_J (a,ma,npl,m0,sini,mass,ap)
        
	implicit none
	real*8 m0,PI,TWOPI,THIRD,GMSUN,dm,MSUN
        integer npl,ma,i,j,NPLMAX
        parameter (NPLMAX=20)
        real*8 sini(ma)
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
	   mass(i+2) = a(7*i+1)*(TWOPI/a(7*i+2)*(m0 + mpold(i+1))**2/
     &               (TWOPI*GMSUN))**THIRD*
     &               dsqrt(1.d0 - a(7*i+3)**2)
           else
              mtotal = m0
              do j = 0, i-1
                 mtotal = mtotal + mass(j+2)
              enddo
              mass(i+2) = a(7*i+1)*(TWOPI/a(7*i+2)*(mtotal
     &                  +mpold(i+1))**2/(TWOPI*GMSUN))**THIRD*
     &                  dsqrt(1.d0 - a(7*i+3)**2)
           endif
           
	   dm = dabs(mass(i+2)-mpold(i+1))/mass(i+2)
	   mpold(i+1) = mass(i+2)
           if (dm.gt.0) goto 101

	   ap(i+1) = (GMSUN*(mtotal + mass(i+2))*(1.d0/a(7*i+2))
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

      subroutine GENINIT_J3 (nbod,ap,a,
     &                       mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill)

c      include 'swift_loglik_Jakub.inc'

      real*8 SMASSYR,MSUN,PI
      parameter (PI=3.14159265358979d0)
      parameter (SMASSYR=4.d0*PI*PI)
      parameter (MSUN=1.32712497d20)

      integer nbod,NPLMAX,i,j
      parameter (NPLMAX=20)
      real*8 mass(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
      real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX)

      real*8 mstar0,m1,m2,frho3,ap(NPLMAX),a(NPLMAX)

      real*8 gm,inc,capom,a1,e1,omega1,capm1,a2,e2,omega2,capm2

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
         
         gm = gm + mass(i)
         rpl(i) = frho3*(1.5d0*mass(i)/2.d0*PI)**0.3333333333d0
         rhill(i) = ap(i-1)*(mass(i)/(3.d0*mass(1)))**0.3333333333d0

         call ORBEL_EL2XV (gm,ialpha,ap(i-1),a(j+3),a(j+6),a(j+7),
     &          a(j+4),a(j+5),xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i))


      enddo

      return
      end
 





    
      subroutine parametric(x,mass,dfdx,ma,npl)
      implicit none
      integer i,j,npl,ma
      real*8 M,dfdx(ma,ma),x(ma),mass(51)
      real*8 ri,rj,rij
      
      M = mass(1)

      do i = 1, 7*npl
         do j = 1, 7*npl
            dfdx(i,j) = 0.0
         enddo
      enddo
      
      do i = 0, npl-1
         
         ri = dsqrt(x(7*i+2)**2 + x(7*i+3)**2 + x(7*i+4)**2)

         dfdx(7*i+5,7*i+2) = 1.d0
         dfdx(7*i+6,7*i+3) = 1.d0
         dfdx(7*i+7,7*i+4) = 1.d0
         dfdx(7*i+1,7*i+5) = -1.d0*x(7*i+2)/ri**3 
         
         dfdx(7*i+2,7*i+5) = -1.d0*(M+x(7*i+1))*(1.0/ri**3
     &             - 3.d0*x(7*i+2)**2/ri**5)
         dfdx(7*i+3,7*i+5) = (M+x(7*i+1))*3.d0*x(7*i+2)*x(7*i+3)/ri**5

         dfdx(7*i+4,7*i+5) = (M+x(7*i+1))*3.d0*x(7*i+2)*x(7*i+4)/ri**5
          
         dfdx(7*i+1,7*i+6) = -1.d0*x(7*i+3)/ri**3
         
         dfdx(7*i+2,7*i+6) = (M+x(7*i+1))*3.d0*x(7*i+3)*x(7*i+2)/ri**5
                   
         dfdx(7*i+3,7*i+6) = -1.d0*(M+x(7*i+1))*(1.0/ri**3
     &             - 3.d0*x(7*i+3)**2/ri**5)

         dfdx(7*i+4,7*i+6) = (M+x(7*i+1))*3.d0*x(7*i+4)*x(7*i+3)/ri**5

         dfdx(7*i+1,7*i+7) = -1.d0*x(7*i+4)/ri**3
         
         dfdx(7*i+2,7*i+7) = (M+x(7*i+1))*3.d0*x(7*i+4)*x(7*i+2)/ri**5
                   
         dfdx(7*i+3,7*i+7) = (M+x(7*i+1))*3.d0*x(7*i+4)*x(7*i+3)/ri**5

         dfdx(7*i+4,7*i+7) = -1.d0*(M+x(7*i+1))*(1.0/ri**3
     &             - 3.d0*x(7*i+4)**2/ri**5)
         
            do j = 0, npl-1
               if (j.ne.i) then
                  rij = dsqrt((x(7*i+2)-x(7*j+2))**2+(x(7*i+3)-x(7*j+3))
     &             **2 + (x(7*i+4)-x(7*j+4))**2)

                  rj = dsqrt(x(7*j+2)**2 + x(7*j+3)**2 + x(7*j+4)**2)

                  dfdx(7*i+2,7*i+5) = dfdx(7*i+2,7*i+5) + x(7*j+1)*(
     &             -1.d0/rij**3 + 3.d0*(x(7*j+2)-x(7*i+2))**2/rij**5)

                  dfdx(7*i+3,7*i+5) = dfdx(7*i+3,7*i+5) + 3.d0*x(7*j
     &             +1)*(x(7*j+2)-x(7*i+2))*(x(7*j+3)-x(7*i+3))/rij**5

                  dfdx(7*i+4,7*i+5) = dfdx(7*i+4,7*i+5) + 3.d0*x(7*j
     &             +1)*(x(7*j+2)-x(7*i+2))*(x(7*j+4)-x(7*i+4))/rij**5
                  
                  dfdx(7*i+2,7*i+6) = dfdx(7*i+2,7*i+6) + 3.d0*x(7*j
     &             +1)*(x(7*j+2)-x(7*i+2))*(x(7*j+3)-x(7*i+3))/rij**5

                  dfdx(7*i+3,7*i+6) = dfdx(7*i+3,7*i+6) + x(7*j+1)*(
     &             -1.d0/rij**3 + 3.d0*(x(7*j+3)-x(7*i+3))**2/rij**5)

                  dfdx(7*i+4,7*i+6) = dfdx(7*i+4,7*i+6) + 3.d0*x(7*j
     &             +1)*(x(7*j+3)-x(7*i+3))*(x(7*j+4)-x(7*i+4))/rij**5

                  dfdx(7*i+2,7*i+7) = dfdx(7*i+2,7*i+7) + 3.d0*x(7*j
     &             +1)*(x(7*j+4)-x(7*i+4))*(x(7*j+2)-x(7*i+2))/rij**5

                  dfdx(7*i+3,7*i+7) = dfdx(7*i+3,7*i+7) + 3.d0*x(7*j
     &             +1)*(x(7*j+4)-x(7*i+4))*(x(7*j+3)-x(7*i+3))/rij**5

                  dfdx(7*i+4,7*i+7) = dfdx(7*i+4,7*i+7) + x(7*j+1)*(
     &             -1.d0/rij**3 + 3.d0*(x(7*j+4)-x(7*i+4))**2/rij**5)

                  dfdx(7*j+1,7*i+5) = ((x(7*j+2)-x(7*i+2))/rij**3 -
     &             x(7*j+2)/rj**3)

                  dfdx(7*j+2,7*i+5) = x(7*j+1)*(1.d0/rij**3 - 3.d0*
     &             (x(7*j+2)-x(7*i+2))**2/rij**5 - 1.d0/rj**3 + 3.d0*
     &             x(7*j+2)**2/rj**5)

                  dfdx(7*j+3,7*i+5) = x(7*j+1)*(-3.d0*(x(7*j+2)-x(7*i
     &             +2))*(x(7*j+3)-x(7*i+3))/rij**5 + 3.d0*x(7*j+2)*
     &             x(7*j+3)/rj**5)

                  dfdx(7*j+4,7*i+5) = x(7*j+1)*(-3.d0*(x(7*j+2)-x(7*i
     &             +2))*(x(7*j+4)-x(7*i+4))/rij**5 + 3.d0*x(7*j+2)*
     &             x(7*j+4)/rj**5)

                  dfdx(7*j+1,7*i+6) = ((x(7*j+3)-x(7*i+3))/rij**3 -
     &             x(7*j+3)/rj**3)

                  dfdx(7*j+2,7*i+6) = x(7*j+1)*(-3.d0*(x(7*j+2)-x(7*i
     &             +2))*(x(7*j+3)-x(7*i+3))/rij**5 + 3.d0*x(7*j+2)*
     &             x(7*j+3)/rj**5)

                  dfdx(7*j+3,7*i+6) = x(7*j+1)*(1.d0/rij**3 - 3.d0*
     &             (x(7*j+3)-x(7*i+3))**2/rij**5 - 1.d0/rj**3 + 3.d0*
     &             x(7*j+3)**2/rj**5)

                  dfdx(7*j+4,7*i+6) = x(7*j+1)*(-3.d0*(x(7*j+4)-x(7*i
     &             +4))*(x(7*j+3)-x(7*i+3))/rij**5 + 3.d0*x(7*j+4)*
     &             x(7*j+3)/rj**5)

                  dfdx(7*j+1,7*i+7) = ((x(7*j+4)-x(7*i+4))/rij**3 -
     &             x(7*j+4)/rj**3)

                  dfdx(7*j+2,7*i+7) = x(7*j+1)*(-3.d0*(x(7*j+4)-x(7*i
     &             +4))*(x(7*j+2)-x(7*i+2))/rij**5 + 3.d0*x(7*j+4)*
     &             x(7*j+2)/rj**5)

                  dfdx(7*j+3,7*i+7) = x(7*j+1)*(-3.d0*(x(7*j+4)-x(7*i
     &             +4))*(x(7*j+3)-x(7*i+3))/rij**5 + 3.d0*x(7*j+4)*
     &             x(7*j+3)/rj**5)

                  dfdx(7*j+4,7*i+7) = x(7*j+1)*(1.d0/rij**3 - 3.d0*
     &             (x(7*j+4)-x(7*i+4))**2/rij**5 - 1.d0/rj**3 + 3.d0*
     &             x(7*j+4)**2/rj**5)
               endif
            enddo
      enddo

      return
      end



                                                                               


      subroutine j2h(x,tjh,ma,npl,mass)

      integer i,j,npl,ma,k,p      
      real*8 tjh(ma,ma),lamda(npl),x(ma),mass(51)
                
     
      do i = 1,ma
         do j = 1,ma
            tjh(i,j) = 0.0
            if ((j.eq.ma).and.(i.eq.ma)) tjh(i,j) = 1.d0
         enddo
      enddo 
      do i = 1, npl
         lamda(i) = mass(1)
      enddo

      do i = 0, npl-1
         do j = 0, i
         lamda(i+1) = lamda(i+1) + x(7*j + 1)
         enddo
      enddo

      do j = 0,npl-1
        do i = 0,npl-1
          if (i.eq.j) then
            do p = 1,7 
              do k = 1,7          
                if (p.eq.k) tjh(7*i+k,7*j+p) = 1.0
                if (i.eq.0) then
                   tjh(7*i+1,7*j+1) = lamda(i+1)/mass(1)
                else
                   tjh(7*i+1,7*j+1) = lamda(i+1)/lamda(i)
                endif
              enddo
            enddo
          elseif (i.lt.j) then
            do p = 2,7
              do k = 2,7
                if(p.eq.k) tjh(7*i+k,7*j+p) = x(7*i+1)/lamda(i+1)
              enddo
            enddo
          endif
        enddo 
      enddo
      
      return
      end





      subroutine h2j (x,thj,ma,npl,mass)
      
      integer i,j,npl,ma,k,p
      real*8 thj(ma,ma),lamda(npl),mass(51),x(ma)

      do i = 1,ma
         do j = 1,ma
            thj(i,j) = 0.0
            if ((j.eq.ma).and.(i.eq.ma)) thj(i,j) = 1.d0
         enddo
      enddo
      do i = 1,npl
         lamda(i) = mass(1)
      enddo

      do i = 0, npl-1
         do j = 0, i
         lamda(i+1) = lamda(i+1) + x(7*j + 1)
         enddo
      enddo

      do j = 0,npl-1
        do i = 0,npl-1
          if (i.eq.j) then
            do p = 1,7
              do k = 1,7
                if (p.eq.k) thj(7*i+k,7*j+p) = 1.0
                if (i.eq.0) then
                   thj(7*i+1,7*j+1) = mass(1)/lamda(i+1)
                else
                   thj(7*i+1,7*j+1) = lamda(i)/lamda(i+1)
                endif
              enddo
            enddo
          elseif (i.lt.j) then
            do p = 2,7
              do k = 2,7
                if (p.eq.k) thj(7*i+k,7*j+p) = -1.d0*x(7*i+1)/lamda(j)
              enddo
            enddo         
          endif 
        enddo
      enddo
   
      return
      end





      subroutine transfer (ma,npl,mass,a,ap,dxda)
      
      implicit none
      real*8 PI,M
      integer npl,ma,i,j,NPLMAX,k,l
      parameter (NPLMAX=20)
      parameter (PI=3.14159265358979d0)
      real*8 dxda(ma,ma),ap(NPLMAX),mp(NPLMAX),a(ma),mass(NPLMAX)
      real*8 cosw,sinw,capm,cape,sine,cose,cosf,sinf,capom,inc
      real*8 coscapom,sincapom,sininc,cosinc
      real*8 dfde,dfdecc,dedm,dedecc
      real*8 r,coswf,sinwf,drdm,drdecc,drdw
      real*8 dcoswfdm,dcoswfdw,dcoswfdecc,dsinwfdm,dsinwfdecc,dsinwfdw
      real*8 dadk,dadn,dv1,dv2,x1,x2,x3,v1,v2,v3 
      real*8 orbel_ehybrid
      real*8 mstar,sini(7),lamda(npl)

      common mstar,sini


      do i = 1,ma
         do j= 1, ma
            dxda(i,j) = 0.0
         enddo
      enddo 
      do i = 1,npl
         lamda(i) = 0.d0
         do j = 1,i+1
            lamda(i) = lamda(i) + mass(j)  
         enddo
      enddo

************* begin loop for each planet ***********************
      do j = 1,npl
         i = 7*(j - 1)

         if (j.eq.1) then                       !!!   Jacobi masses 
            mp(j) = mass(1)/lamda(j)*mass(j+1)
            M = mass(1) - mp(j)
         else
            mp(j) = lamda(j-1)/lamda(j)*mass(j+1)
            M = lamda(j-1) - mp(j)
         endif         

      cosw = dcos(a(4+i))
      sinw = dsin(a(4+i))
      
      inc = a(i+6)
      capom = a(i+7)
      coscapom = dcos(capom)
      sincapom = dsin(capom)
      cosinc = dcos(inc)
      sininc = dsin(inc)
      capm = a(i+5)
      cape = ORBEL_EHYBRID (a(i+3),capm)
      cose = dcos(cape)
      sine = dsin(cape)
      cosf = (cose - a(i+3))/(1.d0 - a(3+i)*cose)
      sinf = dsqrt(1.d0 - a(3+i)**2)*sine/(1.d0 - a(3+i)*cose)

      r = ap(j)*(1.d0 - a(i+3)*cose)
      coswf = cosw*cosf - sinw*sinf
      sinwf = sinw*cosf + cosw*sinf
      
c*********************************************************
      dfde = dsqrt((1.d0 + a(3+i))/(1.d0-a(i+3)))*(1.d0 + cosf)
     &       /(1.d0 + cose)
      dedm = 1.d0/(1.d0 - a(i+3)*cose)
      dedecc = sine/(1.d0 - a(i+3)*cose) 
      dfdecc = dfde*(sine/(1.d0 - a(i+3)**2))
      
      drdm = ap(j)*a(3+i)*sine*dedm
      drdecc = ap(j)*(-1.d0*cose + a(i+3)*sine*dedecc)

      drdw = 0
      dcoswfdm = -1.d0*sinwf*dfde*dedm

      dcoswfdecc = -1.d0*sinwf*(dfdecc + dfde*dedecc) 
      dcoswfdw = -1.d0*sinwf
      dsinwfdm = coswf*dfde*dedm

      dsinwfdecc = coswf*(dfdecc + dfde*dedecc)
      dsinwfdw = coswf

c*****************************************
      dadk = 3.d0*ap(j)*a(i+1)**2*(M+mp(j))**2*(1.d0-a(i+3)**2)**1.5
     &       /(3.d0*a(i+2)*mp(j)**2*(3.d0*M+mp(j)))
      dadn = -2.d0*ap(j)/(3.d0*a(i+2)) - ap(j)*a(i+1)**3*(M+mp(j))**2*
     &       (1.d0-a(i+3)**2)**1.5/(3.d0*a(i+2)**2*mp(j)**2*
     &       (3.d0*M+mp(j)))
      dv1 =  a(i+2)*ap(j)/dsqrt(1.d0-a(i+3)**2)
      dv2 =  a(i+2)*ap(j)*a(i+3)/(1.d0-a(i+3)**2)**1.5

      x1 = coscapom*coswf - sincapom*sinwf*cosinc
      x2 = sincapom*coswf + coscapom*sinwf*cosinc
      x3 = sinwf*sininc
      v1 = coscapom*sinwf + sincapom*cosinc*coswf + a(i+3)*(coscapom*
     &     sinw + sincapom*cosinc*cosw)
      v2 = -1.d0*sincapom*sinwf + coscapom*cosinc*coswf + a(i+3)*(-1.d0*
     &     sincapom*sinw + coscapom*cosinc*cosw)
      v3 = sininc*coswf + a(i+3)*sininc*cosw
*************************************************************
ccccccccccccccccccc     begin matrix      cccccccccccccccc

      dxda(i+1,i+1) = 3.d0*a(i+1)**2*(M + mp(j))**3*
     &        (1.d0-a(i+3)**2)**1.5
     &        /(a(i+2)*(3.d0*M+mp(j))*mp(j)**2)
      dxda(i+2,i+1) = -1.d0*a(i+1)**3*(M + mp(j))**3*
     &        (1.d0-a(i+3)**2)**1.5
     &        /(a(i+2)**2*(3.d0*M + mp(j))*mp(j)**2)
      dxda(i+3,i+1) = 0.0
      dxda(i+4,i+1) = 0.0
      dxda(i+5,i+1) = 0.0
      dxda(i+6,i+1) = 0.0
      dxda(i+7,i+1) = 0.0

      dxda(i+1,i+2) = dadk*(1.d0 - a(i+3)*cose)*x1
      dxda(i+2,i+2) = dadn*(1.d0 - a(i+3)*cose)*x1
      dxda(i+3,i+2) = drdecc*x1 + r*(coscapom*dcoswfdecc - sincapom*
     &                cosinc*dsinwfdecc)
      dxda(i+4,i+2) = drdw*x1 + r*(coscapom*dcoswfdw - sincapom*
     &                cosinc*dsinwfdw)
      dxda(i+5,i+2) = drdm*x1 + r*(coscapom*dcoswfdm - sincapom*
     &                cosinc*dsinwfdm)
      dxda(i+6,i+2) = r*sincapom*sinwf*sininc
      dxda(i+7,i+2) = r*(-1.d0*sincapom*coswf - coscapom*sinwf*cosinc)

      dxda(i+1,i+3) = dadk*(1.d0 - a(i+3)*cose)*x2
      dxda(i+2,i+3) = dadn*(1.d0 - a(i+3)*cose)*x2
      dxda(i+3,i+3) = drdecc*x2 + r*(sincapom*dcoswfdecc + coscapom*
     &                cosinc*dsinwfdecc)
      dxda(i+4,i+3) = drdw*x2 + r*(sincapom*dcoswfdw + coscapom*
     &                cosinc*dsinwfdw)
      dxda(i+5,i+3) = drdm*x2 + r*(sincapom*dcoswfdm + coscapom*
     &                cosinc*dsinwfdm)
      dxda(i+6,i+3) = -1.d0*r*coscapom*sinwf*sininc
      dxda(i+7,i+3) = r*(coscapom*coswf - sincapom*sinwf*cosinc)

      dxda(i+1,i+4) = dadk*(1.d0 - a(i+3)*cose)*x3
      dxda(i+2,i+4) = dadn*(1.d0 - a(i+3)*cose)*x3
      dxda(i+3,i+4) = drdecc*x3 + r*dsinwfdecc*sininc
      dxda(i+4,i+4) = drdw*x3 + r*dsinwfdw*sininc
      dxda(i+5,i+4) = drdm*x3 + r*dsinwfdm*sininc
      dxda(i+6,i+4) = r*sinwf*cosinc
      dxda(i+7,i+4) = 0.0

      dxda(i+1,i+5) = -1.d0*dadk*a(i+2)/dsqrt(1.d0 - a(i+3)**2)*v1
      dxda(i+2,i+5) = (a(i+2)*dadn + ap(j))*(-1.d0/
     &                dsqrt(1.d0 - a(i+3)**2)*v1)
      dxda(i+3,i+5) = -1.d0*dv2*v1 - dv1*(coscapom*dsinwfdecc + 
     &                sincapom*cosinc*dcoswfdecc + coscapom*sinw +
     &                sincapom*cosinc*cosw)
      dxda(i+4,i+5) = -1.d0*dv1*(coscapom*coswf - sincapom*cosinc*sinwf
     &                + a(i+3)*(coscapom*cosw - sincapom*cosinc*sinw))
c      dxda(i+5,i+4) = -1.d0*dv1*(dsinwfdm + sinw/sine)
      dxda(i+5,i+5) = -1.d0*dv1*(coscapom*dsinwfdm + sincapom*cosinc*
     &                dcoswfdm)
      dxda(i+6,i+5) = -1.d0*dv1*(-1.d0*sincapom*sininc*coswf - a(i+3)*
     &                sincapom*sininc*cosw)
      dxda(i+7,i+5) = -1.d0*dv1*(-1.d0*sincapom*sinwf + coscapom*cosinc*
     &                coswf + a(i+3)*(-1.d0*sincapom*sinw + coscapom*
     &                cosinc*cosw)) 

      dxda(i+1,i+6) = dadk*a(i+2)/dsqrt(1.d0 - a(i+3)**2)*v2
      dxda(i+2,i+6) = (a(i+2)*dadn + ap(j))*1.d0/dsqrt(1.d0
     &                - a(i+3)**2)*v2
      dxda(i+3,i+6) = dv2*v2 + dv1*(-1.d0*sincapom*dsinwfdecc + coscapom
     &                *cosinc*dcoswfdecc - sincapom*sinw + coscapom*
     &                cosinc*cosw)
      dxda(i+4,i+6) = dv1*(-1.d0*sincapom*coswf - coscapom*cosinc*sinwf
     &                + a(i+3)*(-1.d0*sincapom*cosw - coscapom*cosinc*
     &                sinw))
c      dxda(i+5,i+5) = dv1*(dcoswfdm + cosw/sine)
      dxda(i+5,i+6) = dv1*(-1.d0*sincapom*dsinwfdm + coscapom*cosinc*
     &                dcoswfdm)
      dxda(i+6,i+6) = dv1*(-1.d0*coscapom*sininc*coswf - a(i+3)*
     &                coscapom*sininc*cosw)
      dxda(i+7,i+6) = dv1*(-1.d0*coscapom*sinwf - sincapom*cosinc*coswf
     &                + a(i+3)*(-1.d0*coscapom*sinw - sincapom*cosinc
     &                *cosw))

      dxda(i+1,i+7) = dadk*a(i+2)/dsqrt(1.d0 - a(i+3)**2)*v3
      dxda(i+2,i+7) =  (a(i+2)*dadn + ap(j))*1.d0/dsqrt(1.d0
     &                - a(i+3)**2)*v3
      dxda(i+3,i+7) = dv2*v3 + dv1*(sininc*dcoswfdecc + sininc*cosw)
      dxda(i+4,i+7) = dv1*(-1.d0*sininc*sinwf - a(i+3)*sininc*sinw)
      dxda(i+5,i+7) = dv1*sininc*dcoswfdm
      dxda(i+6,i+7) = dv1*(cosinc*coswf + a(i+3)*cosinc*cosw)
      dxda(i+7,i+7) = 0.0


      enddo
 
      return
      end
      
      
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

      subroutine bs_der_parametric(nbod,ma,mass,j2rp2,j4rp4,ybs,dy,x)

      include 'swift_loglik_Jakub.inc'
      include 'bs.inc'

c...  Inputs Only: 
      integer nbod,ma,npl
      real*8 mass(nbod),ybs(N6DBS),x(ma)  !x(ma) is x and v
      real*8 j2rp2,j4rp4

c...  Output
      real*8 dy(N6DBS)
     
c...  Internals
      integer i,j,p,nbodm,l
      real*8 dfdx(ma,ma),sz(ma,ma),z(ma,ma)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
      real*8 axb(NPLMAX),ayb(NPLMAX),azb(NPLMAX)      
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c----
c...  Executable code 
c...  move things so that I can deal with it
c------------------------first for motion equations-------------------------
      do i=1,nbod
         j = 6*(i-1)
         xb(i) = ybs(j+1)
         yb(i) = ybs(j+2)
         zb(i) = ybs(j+3)
         vxb(i) = ybs(j+4)
         vyb(i) = ybs(j+5)
         vzb(i) = ybs(j+6)
      enddo
 
      call tu4_getaccb(nbod,mass,j2rp2,j4rp4,xb,yb,zb,axb,ayb,azb)      

 
c------moves things back by array
      do i = 1,nbod
         j = 6*(i-1)
         dy(j+1) = vxb(i)
         dy(j+2) = vyb(i)
         dy(j+3) = vzb(i)
         dy(j+4) = axb(i)
         dy(j+5) = ayb(i)
         dy(j+6) = azb(i)
      enddo

 
c-----------------------second, for parametric equations------------------      
      npl = nbod - 1

      call coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &         xh,yh,zh,vxh,vyh,vzh)

      do i = 2,nbod
         j = 7*(i-2)
         x(j+2) = xh(i)
         x(j+3) = yh(i)
         x(j+4) = zh(i)
         x(j+5) = vxh(i)
         x(j+6) = vyh(i)
         x(j+7) = vzh(i)
      enddo

      call parametric(x,mass,dfdx,ma,npl)

      do i = 1,ma
         do j = 1,ma
               sz(i,j) = 0.0
         enddo
      enddo
       
      l = 1
      do j = 1,ma
         do i = 1,ma
           z(i,j) = ybs(6*nbod+l)
           l = l + 1
         enddo
      enddo 

      do j = 1,ma
         do i = 1,ma
            do p = 1,ma
               sz(i,j) = sz(i,j) + z(i,p)*dfdx(p,j) 
            enddo
         enddo
      enddo

      l = 1
      do i = 1,ma
         do j = 1,ma
            dy(6*nbod+l) = sz(j,i)
            l = l + 1
         enddo
      enddo

      return
      end





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

      subroutine bs_int_parametric(nbod,ma,mass,j2rp2,j4rp4,t,h0,
     &                             y,x,eps)

      include 'swift_loglik_Jakub.inc'
      include 'bs.inc'

c...  Inputs Only: 
      integer nbod,npl,ma
      real*8 mass(nbod),h0,eps,j2rp2,j4rp4,x(ma)

c...  Input & Output
      real*8 t,y(N6DBS)
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

      npl = nbod - 1
      n = 6*nbod + ma*ma

      xa=t
      call bs_der_parametric(nbod,ma,mass,j2rp2,j4rp4,y,dy,x)

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
            t=xa
         
            do i1=1,i1max
               t=t+hd       
               call bs_der_parametric(nbod,ma,mass,j2rp2,j4rp4,y,dy,x)        
                  do i=1,n
                  ii=12*i
                  tp(ii-1)=dmax1(tp(ii-1),dabs(y(i)))
                  eta2=tp(ii-3)+h*dy(i)
                  tp(ii-3)=y(i)
                  y(i)=eta2
               enddo 
            enddo
         
            call bs_der_parametric(nbod,ma,mass,j2rp2,j4rp4,y,dy,x)

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
                  t=xb
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
      return
c
      end      !  bs_int
c-----------------------------------------------------------------------------





c*************************************************************************
c                            BS_STEP.F
c*************************************************************************
c This subroutine has same i/o as STEP_KDK but here we use a BS
c Steps both massive and test particles in the same subroutine.
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp           ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c
c Remarks:  
c Authors:  Hal Levison
c Date:    5/17/93
c Last revision: 2/24/94

      subroutine bs_step_parametric(i1st,time,nbod,ma,mass,j2rp2,j4rp4,
     &           xh,yh,zh,vxh,vyh,vzh,x,z,dt,eps)	

      include 'swift_loglik_Jakub.inc'
      include 'bs.inc'

c...  Inputs Only: 
      integer nbod,i1st,ma
      real*8 mass(nbod),dt,time,x(ma),j2rp2,j4rp4

c...  Inputs and Outputs:

      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 z(ma,ma)

c...  Internals
      integer j,i,ntpi,link(NTPMAX),l
      integer istattmp(NTPMAX,NSTAT),jj,i1stin
      real*8 eps
      real*8 tfake,dttmp,msys,y(N6DBS)
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      

c----
c...  Executable code 

c...  set things up if this is the initial call
      call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &     xb,yb,zb,vxb,vyb,vzb,msys)    

      do i = 1,nbod
        j = 6*(i-1)
        y(j+1) = xb(i)
        y(j+2) = yb(i)
        y(j+3) = zb(i)
        y(j+4) = vxb(i)
        y(j+5) = vyb(i)
        y(j+6) = vzb(i)
      enddo

      l = 1
      do i = 1,ma
       do j = 1,ma
         y(6*nbod + l) = z(j,i)
         l = l + 1
       enddo
      enddo

      tfake = 0.0d0
      dttmp = dt


c      do while(tfake.lt.dt)
      do while( (abs(tfake-dt)/dt) .gt. 1.0e-7 )    ! just to be real safe
      call bs_int_parametric(nbod,ma,mass,j2rp2,j4rp4,
     &                          tfake,dttmp,y,x,eps)
         dttmp = dt - tfake
      enddo

      do i = 1,nbod
       j = 6*(i-1) 
       xb(i) = y(j+1)
       yb(i) = y(j+2)
       zb(i) = y(j+3)
       vxb(i) = y(j+4)
       vyb(i) = y(j+5) 
       vzb(i) = y(j+6)
      enddo
  
      l = 1
      do i = 1,ma
       do j = 1,ma
        z(j,i) = y(6*nbod + l)
        l = l + 1
       enddo
      enddo
      call coord_b2h(nbod,mass,xb,yb,zb,vxb,vyb,vzb,
     &         xh,yh,zh,vxh,vyh,vzh)
      return

      end   ! bs_step
c------------------------------------------------------------------------





c**********************************************************************
c		      SWIFT_BS.F
c**********************************************************************
c
c                 INCLUDES CLOSE ENCOUNTERS
c                 To run, need 3 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c                   test particle file like   tp.in
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    8/5/93
c Last revision: 12/27/96

        subroutine integrate(x,mass,ymod,dyda,t,a,ap,nbod,ma,ndata,
     & epsil,dt)
	include 'swift_loglik_Jakub.inc'

        integer IO_NBITS
        parameter(IO_NBITS=6)
	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
        real*8 xht(NPLMAX),yht(NPLMAX),zht(NPLMAX)
        real*8 vxht(NPLMAX),vyht(NPLMAX),vzht(NPLMAX)


	integer istat(NTPMAX,NSTAT),i1st
	integer nbod,ntp,nleft
	integer iflgchk,iub,iuj,iud,iue
        real*8 rstat(NTPMAX,NSTATR)

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 time,tout,tdump,tfrac,eoff

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        logical*1 lflg(0:IO_NBITS-1)

        integer ndata,ma,i,j,nd,flag,ierr
        real*8 ymod(ndata),dyda(ndata,ma),t(ndata),a(ma),x(ma)
        real*8 z(ma,ma),ap(NPLMAX),h,eps,epsil,deltat
         

	character*80 outfile,inparfile,inplfile,intpfile,fopenstat


c-----
c...    Executable code 


 
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

        xh(1) = 0.0
        yh(1) = 0.0
        zh(1) = 0.0
        vxh(1) = 0.0
        vyh(1) = 0.0
        vzh(1) = 0.0
        do i = 2,nbod
           j = 7*(i-2)
           xh(i) = x(j+2)
           yh(i) = x(j+3)
           zh(i) = x(j+4)
           vxh(i) = x(j+5)
           vyh(i) = x(j+6)
           vzh(i) = x(j+7)
        enddo
        do i = 1,ndata        !  initialize ymod and dyda and z
           do j = 1,ma
              dyda(i,j) = 0.0
           enddo
        enddo

        do i = 1,ma
           do j = 1,ma
              if (i.eq.j) then
                 z(i,j) = 1.d0
              else 
                 z(i,j) = 0.0
              endif
           enddo
        enddo       
c***************here's the big loop *************************************

        nd = 1
c-------output the first round--------
        if (t(1).lt.1.d-10) then
           call output(ndata,nd,nbod,ma,x,a,z,ymod,dyda,mass,ap)
           nd = nd + 1
        endif
	  do while ( time .le. tstop )
             h = dt
             flag = 0                   ! flag for controling output
             do i = 1,ndata    
                if ((time.lt.t(i)).and.((t(i)-time).le.dt)) then
                   flag = 1
                   h = t(i) - time
                   goto 555
                endif
             enddo

 555   	     call bs_step_parametric(i1st,time,nbod,ma,mass,j2rp2,j4rp4,
     &            xh,yh,zh,vxh,vyh,vzh,x,z,h,epsil)
             do i = 2,nbod
                j = 7*(i-2)
                x(j+2) = xh(i)
                x(j+3) = yh(i)
                x(j+4) = zh(i)
                x(j+5) = vxh(i)
                x(j+6) = vyh(i)
                x(j+7) = vzh(i)
             enddo

             if (flag.eq.1) then
             call output(ndata,nd,nbod,ma,x,a,z,ymod,dyda,mass,ap)
c             write(*,*) "TEST:", time, t(1),dt
             nd = nd + 1
             endif

             time = time + h              

c             if(btest(iflgchk,4))  then    ! bit 4 is set
c                call discard(t,dt,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
c     &               xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
c     &               qmin,lclose,rplsq,istat,rstat)
c                call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
c     &               vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
c     &               'discard.out',fopenstat,nleft)
c             else
c                nleft = ntp
c             endif

	enddo

        return
        end    ! swift_bs.f
c---------------------------------------------------------------------




         subroutine output(ndata,nd,nbod,ma,x,a,z,ymod,dyda,mass,ap)
         implicit none
         integer ndata,nbod,nd,npl,i,j,p,ma,NPLMAX
         parameter (NPLMAX=20)
         real*8 mstar,sini(7)
         real*8 x(ma),z(ma,ma),mass(nbod),ymod(ndata),dyda(ndata,ma)
         real*8 tjh(ma,ma),thj(ma,ma),zout(ma,ma),ztest(ma,ma)
         real*8 mp(NPLMAX),ap(NPLMAX),dydx(ma),zxa(ma,ma),dxda(ma,ma)
         real*8 a(ma),mtotal,test(ma)
 	 real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
       	 real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
         
         common mstar,sini      

            npl = nbod - 1

            do i = 2,nbod       ! begin ymod
               j = 7*(i-2) 
               xh(i) = x(j+2)
               yh(i) = x(j+3)
               zh(i) = x(j+4)
               vxh(i) = x(j+5)
               vyh(i) = x(j+6)
               vzh(i) = x(j+7)              
            enddo
         
            mtotal = 0.d0         
            do i = 1,nbod
               mtotal = mtotal + mass(i)
            enddo 
            
            do i = 2,nbod
               ymod(nd) = ymod(nd) + mass(i)/mtotal*vzh(i)
            enddo  

          
            call j2h(x,tjh,ma,npl,mass)   ! begin dyda

            do i = 1,7*npl
               do j = 1,7*npl
                  zout(i,j) = 0.0
               enddo
            enddo
                 
            do j = 1,7*npl
               do i = 1,7*npl
                  do p = 1,7*npl
                     zout(i,j) = zout(i,j) +
     &                              z(p,j)*tjh(i,p)
                  enddo
               enddo
            enddo      

            do i = 1,npl
               mp(i) = mass(i+1)
            enddo 

            call transfer (ma,npl,mass,a,ap,dxda)       !    here need a dxda

           
            do i = 1,7*npl
               do j = 1,7*npl
                  zxa(i,j) = 0.0       !first two parts of the transfer equation
               enddo
            enddo

            do j = 1,7*npl
               do i = 1,7*npl
                  do p = 1,7*npl
                     zxa(i,j) = zxa(i,j) + dxda(i,p)*zout(p,j)
                  enddo
               enddo
            enddo

           
            do i = 1,7*npl
               dydx(i) = 0.0
            enddo
            mtotal = 0.d0
            do i = 1,nbod
               mtotal = mtotal + mass(i)
            enddo
            do i = nbod,2,-1
               dydx(7*(i-2)+7) = mass(i)/mtotal

               dydx(7*(i-2)+1) = (mtotal - mass(i))/mtotal**2
     &                          *vzh(i)
               
            enddo 

            do j = 1,7*npl
               do i = 1,7*npl
                  dyda(nd,j) = dyda(nd,j) + zxa(j,i)*dydx(i)
               enddo
            enddo
        
                        
         return
 
         end

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


