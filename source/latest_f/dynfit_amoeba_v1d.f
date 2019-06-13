
c*************************************************************************
C****************   Radial Velocity dynamical model      *****************
c*************************************************************************
 

      implicit none
      real*8 PI
      parameter (PI=3.14159265358979d0)
      integer npl,ndset,idset,ndata,ma,mfit,i,j,NDSMAX,NPLMAX,MMAX,k
      parameter (NDSMAX=20, NPLMAX=20, MMAX=200)
      real*8 mstar,sini(NPLMAX),sig2i,dy, loglik,sig2l,ftol      
      integer idsmax(NDSMAX),ia(MMAX),ts(10000) ,nt, iter, ii
      integer writeflag_best_par,hkl
      integer writeflag_RV,writeflag_fit, amoebastarts
      real*8 t(10000),x(10000),y(10000),sig(10000),ys(10000),sigs(10000)
      real*8 a(MMAX),covar(MMAX,MMAX),alpha(MMAX,MMAX)
      real*8 chisq,alamda,ochisq,dchisq,red_chisq, p(MMAX+1,MMAX),
     &       yamoeba(MMAX+1), loglikk, ologlikk, dloglikk
      real*8 sigscale,t0,t_max,twopi,dt, epoch, epsil,deltat
      real*8 rms,ymod(10000),dyda(10000,MMAX),jitt(NDSMAX)
      real*4 t_stop, when_to_kill,model_max,model_min
      
      
      external rvkep, compute_abs_loglik
      character*80 infile

      
      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      ftol=0.001d0
      
      read (*,*) epsil,deltat, amoebastarts,
     &          when_to_kill, nt, model_max,model_min      
      
      read (*,*) mstar,
     &          writeflag_best_par, 
     &	             writeflag_RV,writeflag_fit 
   

      call io_read_data (ndata,t,ts,ys,sigs,jitt,
     & 	           epoch,t0,t_max,a,ia,ma,mfit,hkl)

  
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
     &        compute_abs_loglik,ndata,t,ys,ymod,dyda,
     &        ts,sigs,epsil,deltat,i,hkl)
     
     
         call amoeba(p,yamoeba,MMAX+1,MMAX,mfit,ftol,compute_abs_loglik,
     & iter,ndata,t,ys,ymod,dyda,ma,ts,sigs,a,ia,epsil,deltat,hkl)
 
c         write (*,*) i, dloglikk, loglikk     
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
     &  nt, model_max, model_min,hkl)


      end


      subroutine compute_abs_loglik(ndata,x,y,a2,ymod,dyda,ma,mfit,ts,
     & sig,loglik, num,a,ia,epsil,deltat,hkl)
      implicit none
      
      integer MMAX,npl,ndset,NDSMAX,idset,num, mfit,hkl    
      parameter (MMAX=200, NDSMAX=20)
      real*8 twopi, loglik
      parameter (twopi=6.28318530717958d0)
      integer ndata, i, j, ma, ts(10000), ia(MMAX), idsmax(NDSMAX)
      real*8 dy, sig(10000), dyda(10000,MMAX), x(10000), y(10000)
      real*8 ymod(10000), a(MMAX),a2(mfit),a3(MMAX),sig2i,y_in(10000)
     & , y2(10000), epsil,deltat
      
   
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
      
c      write(*,*) a3(2)
      loglik=0.d0
 
        call RVKEP (x,a3,y2,dyda,ma,ndata,epsil,deltat,hkl)
        
        do i = 1,ndata
              idset = ts(i)
              y_in(i) = y(i) - a3(7*npl+idset)- 
     &                 a3(7*npl+2*ndset+1)*(x(i)/86400.d0)

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
     & x,z,ymod,dyda,ts,sig,epsil,deltat,it,hkl)
      integer MMAX, NDSMAX,ma,ts(10000), ndata,mp,np,mfit,it,hkl
      parameter(MMAX=200, NDSMAX=20)
      REAL*8 ftol,p(mp,np),y(mp),a(MMAX), a2(mfit),fr,frjitt
      real*8 x(10000),z(10000),ymod(10000)
      real*8 dyda(MMAX), sig(10000), loglik,epsil,deltat
      parameter(fr=0.01, frjitt=0.05)
      INTEGER i,j,k, ia(MMAX), idsmax(NDSMAX)
      external funk
    
      common /DSBLK/ npl,ndset,idsmax,idset

      k=0
      do j=1,ma
          if(ia(j).ne.0) then
          k=k+1
          p(1,k)=a(j)
          do i=2,mfit+1
          
              if (hkl.eq.0) then
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
              
              else
              
              if (k.eq.(i-1)) then
                  if (j.gt.(7*npl+ndset)) then
                  p(i,k)=(1+frjitt)*(p(1,k)+0.1)
                  else 
                  if (mod(j,7).eq.2) then
                     p(i,k)=(1+fr)*(p(1,k) + 0.1)
                  else if (mod(j,7).eq.3) then
                     p(i,k)=(1+frjitt)*(p(1,k)+0.0001)
                  else if (mod(j,7).eq.4) then
                     p(i,k)=(1+frjitt)*(p(1,k)+0.0001)                     
                  else
                     p(i,k)=(1+fr)*(p(1,k)+0.1)
                  endif
                  endif
              else
                  p(i,k)=p(1,k)
              endif              
 
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
     &              a,ia,epsil,deltat,hkl)
          y(i)=loglik
c          write(*,*) loglik 	  
      enddo
c      write(*,*) loglik 
      return
      end
      
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,epsil,deltat,hkl)
      implicit none
      INTEGER iter,mp,ndim,np,NMAX,ITMAX, MMAX,ma,ts(10000), ndata
      REAL*8 ftol,p(mp,np),y(mp),x(10000),z(10000),ymod(10000)
      PARAMETER (NMAX=20,ITMAX=50000,MMAX=200)
      real*8 dyda(10000,MMAX), sig(10000), loglik, a(MMAX)
      EXTERNAL funk
      INTEGER i,ihi,ilo,inhi,j,m,n, ia(MMAX),hkl
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
     & dyda,ma,ts,sig,a,ia,epsil,deltat,hkl)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,epsil,deltat,hkl)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0,ndata,x,z,ymod,
     & dyda,ma,ts,sig,a,ia,epsil,deltat,hkl)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              call funk(ndata,x,z,psum,ymod,dyda,ma,ndim,ts,sig,loglik,
     &                 i,a,ia,epsil,deltat,hkl)
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
     & dyda,ma,ts,sig,a,ia,epsil,deltat,hkl)
      implicit none
      INTEGER ihi,mp,ndim,np,NMAX, MMAX, ma, ts(10000),ndata
      PARAMETER (NMAX=20, MMAX=200)
      REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),x(10000),z(10000),
     & ymod(10000), epsil, deltat
      real*8 dyda(10000,MMAX), sig(10000),loglik
      EXTERNAL funk
      INTEGER j, ia(MMAX),hkl
      REAL*8 fac1,fac2,ytry,ptry(ndim), a(MMAX)
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      call funk(ndata,x,z,ptry,ymod,dyda,ma,ndim,ts,sig,loglik,ihi,
     & a,ia,epsil,deltat,hkl)

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
     &               t0,t_max,ar,iar,ma,mfit,hkl)  


      implicit none
      integer ndset,idset,ndata,NDSMAX,NPLMAX,MMAX,npl
      real*8 t(10000),y(10000),sig(10000),ys(10000),sigs(10000)
      parameter (NDSMAX=20,NPLMAX=20,MMAX=200)
      integer idsmax(NDSMAX),ts(10000),hkl
      real*8 jitter(NDSMAX),t0,t_max,epoch,ar(MMAX),off(NDSMAX), PI
      parameter(PI=3.14159265358979d0)
      integer i,k,j, iar(MMAX), u_off(NDSMAX), u_jit(NDSMAX), ma, mfit
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
      
      do i = 1,ndset
          ar(7*npl+i)=off(i)
          iar(7*npl+i)=u_off(i)
          ar(7*npl+ndset+i)=jitter(i)
          iar(7*npl+ndset+i)=u_jit(i)
          enddo
      
c           write(*,*) ' Number of Planets: '
      if (npl.gt.NPLMAX) stop ' KEPFIT: npl > NPLMAX.'

      ma = 7*npl + 2*ndset + 1
      
c      write(*,*) 'Initial K, P, e, w, M0,Inc,Capom and their flags: '
      do j = 1,npl
          i = 7*(j-1)
          read (*,*) ar(i+1),ar(i+2),ar(i+3),ar(i+4),ar(i+5),ar(i+6),
     &               ar(i+7)
          read (*,*) iar(i+1),iar(i+2),iar(i+3),iar(i+4),iar(i+5),
     &               iar(i+6),iar(i+7)
 
 
c          ar(i+2) = 2.d0*PI/(ar(i+2)*8.64d4)         ! mean motion 
          ar(i+2) = ar(i+2)*8.64d4         ! second as unit
 
      enddo
  
c      write (*,*) 'linear trend:'      
      read (*,*) ar(7*npl + 2*ndset + 1)
      read (*,*) iar(7*npl + 2*ndset + 1)   
      

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
C
 

      subroutine io_write_bestfitpa_ewcop_fin (a,covar,t,ys,ndata,ts,
     &           ma,mfit,t0,t_max,sigs,chisq,rms,loglik,writeflag_RV,
     &           writeflag_best_par,writeflag_fit,jitter,epsil,
     &           deltat,nt, model_max,model_min,hkl)
   
      implicit none 
      real*8 PI
      integer MMAX,NDSMAX,NPLMAX 
      parameter (PI=3.14159265358979d0,MMAX=200  ,NDSMAX=20,NPLMAX=20)
      real*8 a(MMAX),ia(MMAX),t(10000),ymod(10000),ys(10000)
      real*8 covar(MMAX,MMAX),dyda(10000,MMAX),AU,day
      real*8 rms,mstar,sini(NPLMAX),mass(NPLMAX),ap(NPLMAX)
      integer ts(10000),nbod,nt,writeflag_RV,
     &           writeflag_best_par,writeflag_fit,hkl
      real*8 t0,t1,t2,dt,offset,t_max,chisq,loglik,dy,sig2i,twopi
      real*8 x(10000),y(10000),sigs(10000),jitter(NDSMAX)
      integer i,j,npl,ndset,ndata,idset,mfit,ma,idsmax(NDSMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),vxh(NPLMAX),vyh(NPLMAX)
     &       ,vzh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),vxj(NPLMAX),vyj(NPLMAX)
     &       ,vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX),deltat,epsil
      real*8 swift_mass(NPLMAX),s_mass(NPLMAX),j_mass(NPLMAX)
      real*4 model_max,model_min,best_w,best_we
      parameter (AU=1.49597892d11, day = 86400.d0)


      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini
            
      nbod = npl+1
      rms = 0.d0
      twopi = 2.0d0*PI 
 
ccccccccccccccccccc t[JD], obs., cal., O-C   ccccccccccccc   


      call RVKEP(t,a,ymod,dyda,ma,ndata,epsil,deltat,hkl)

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
     &               a(7*npl  + 2*ndset + 1)*(t(i)/8.64d4)

c          write(*,*) a(7*npl+idset), a(7*npl  + 2*ndset + 1)

          if (writeflag_RV.gt.0) then
          write(*,*) t0 + t(i)/8.64d4 ,ymod(i),ys(i) +
     &                a(7*npl  + 2*ndset + 1)*(t(i)/8.64d4), 
     &                ys(i) - ymod(i),sigs(i),ts(i)

          endif
     
          sig2i = 1.d0/(sigs(i)**2 + a(7*npl+ndset+idset)**2)

          dy =  ys(i) -ymod(i)  
 

 	      chisq  = chisq + dy*dy*sig2i
 
c          write(*,*) "TEST:",loglik,dy,ymod(i),a(7*npl+ndset+idset)
	      loglik =  loglik - 0.5*dy*dy*sig2i -
     &               0.5*dlog(twopi*(sigs(i)**2
     &                + a(7*npl+ndset+idset)**2)) 
     
 

          rms = rms + dy**2
      enddo

      rms = dsqrt(rms/dble(ndata))

      if(writeflag_best_par.gt.0) then

        write(*,*) 'loglik, reduced chi^2, chi^2, rms:'
        write(*,*) loglik, chisq/dble(ndata-mfit),chisq, rms


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

!             a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)
!              a(i+2) = a(i+2)/8.64d4

              write(*,*) a(i+1),
     &                a(i+2)/8.64d4,a(i+3),best_w, a(i+5)*180.d0/PI,
     &                dmod(a(i+6)*180.d0/PI,180.d0),
     &                dmod(a(i+7)*180.d0/PI,360.d0)
              write(*,*) dsqrt(covar(i+1,i+1)),
     &                dsqrt(covar(i+2,i+2)),
c     &                2.d0*PI/a(i+2)**2*dsqrt(covar(i+2,i+2))/8.64d4,
     &                dsqrt(covar(i+3,i+3)),best_we,
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
      
          write (*,*) 'Jitters for each dataset:'
          do j = 1,ndset
              write (*,*) a(7*npl+ndset+j)
              write (*,*) '0'
          enddo 
       
          write (*,*) 'linear trend  [m/s per day]:'
          write (*,*) a(7*npl + 2*ndset + 1)
          write (*,*) dsqrt(covar(7*npl + 2*ndset + 1,
     &                7*npl + 2*ndset + 1))

          write(*,*) ' ndata =',ndata
          write(*,*) ' mfit =',mfit
          write(*,*) ' RMS =',rms
          write(*,*) ' Chi^2 =',chisq/dble(ndata-mfit)
          write(*,*) ' epoch = ', t0


 
          call MA_J_cop_fin (a,ma,npl,mstar,sini,mass,ap,hkl)
 
          call GENINIT_J3_ewcop (nbod,ap,a,
     &                          mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)

c          do i = 1,npl
c             j = 7*(i-1)
         
c             a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)   
             
c          enddo        


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
          call RVKEP (x,a,ymod,dyda,ma,nt,epsil,deltat,hkl)
          do i = 1,nt
             write(*,*) t0 + x(i)/8.64d4, 
     &       ymod(i) + a(7*npl  + 2*ndset + 1)*(x(i)/8.64d4)           
          enddo

      endif

 

      return

      end      

 

 
      subroutine RVKEP (t,a,ymod,dyda,ma,ndata,epsil,dt,hkl) 
c      subroutine RVKEP(t,a,ymod,dyda,ma,ndata,ia,hkl)
      implicit none
      real*8 PI,TWOPI,eps, dt 
      parameter (PI=3.14159265358979d0,eps=1.d-6)
      parameter (TWOPI=2.0d0*PI)
      integer npl,ma,i,j,NPLMAX,na,ndset,NDSMAX 
      parameter (NPLMAX=20,NDSMAX=20)
      real*8 t(ndata),ymod(ndata),a(ma),dyda(ndata,ma),x(ma),ia(ma)
      real*8 mstar,ap(NPLMAX),mass(NPLMAX),epsil,deltat,a2(ma)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),vxh(NPLMAX),vyh(NPLMAX)
     &       ,vzh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX),vxj(NPLMAX),vyj(NPLMAX)
     &       ,vzj(NPLMAX)
      real*8 rpl(NPLMAX),rhill(NPLMAX)
      real*8 ah(ma),ahh(ma),ymodhb(ndata),ymodha(ndata)
      real*8 sini(NPLMAX),sinih,sinihh,dydsini(ndata),factor
      integer ts(ndata),correct,hkl,idset,ndata,nbod,z,p
      
      integer idsmax(NDSMAX)

      common /DSBLK/ npl,ndset,idsmax,idset
      common mstar,sini

      nbod = npl + 1
      na = 7*npl
      
      
      do i = 1,ma
          a2(i)=a(i)
      enddo     
      
      if (hkl.eq.0) then

          do i = 1,npl
             j = 7*(i-1)
             
c             a2(j+2) = 2.d0*PI/(a2(j+2)*8.64d4)
             
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
             if (a2(j+4).lt.0.d0) a2(j+4)=dmod(a2(j+4)+2.d0*PI,2.d0*PI)  
             if (a2(j+5).lt.0.d0) a2(j+5)=dmod(a2(j+5)+2.d0*PI,2.d0*PI) 
             if (a2(j+4).gt.2.d0*PI) a2(j+4)=dmod(a2(j+4),2.d0*PI )  
             if (a2(j+5).gt.2.d0*PI) a2(j+5)=dmod(a2(j+5),2.d0*PI )         
             if (a2(j+6).lt.0.d0) a2(j+6)=dmod(a2(j+6)+2.d0*PI,2.d0*PI)  
             if (a2(j+7).lt.0.d0) a2(j+7)=dmod(a2(j+7)+2.d0*PI,2.d0*PI) 
             if (a2(j+6).gt.2.d0*PI) a2(j+6)=dmod(a2(j+6),2.d0*PI)  
             if (a2(j+7).gt.2.d0*PI) a2(j+7)=dmod(a2(j+7),2.d0*PI)                         
          enddo  


      else   
            
          do i = 1,npl
             j = 7*(i-1)
             
c             a2(j+2) = 2.d0*PI/(a2(j+2)*8.64d4)
             
             
             if (a2(j+2).lt.0.d0) then  ! if P<0, set P>0 
                a2(j+2) = abs(a2(j+2))
             endif                   
             
             if (a2(j+1).lt.0.d0) then  ! if K<0, set K>0 and w = w+PI 
                a2(j+4) = -1.d0*a2(j+4)       !     which is h = -h, k = -k
                a2(j+3) = -1.d0*a2(j+3)
                a2(j+1) = abs(a2(j+1))    
             endif
          
                         
             if (a2(j+5).lt.0.d0) a2(j+5)=dmod(a2(j+5)+2.d0*PI,2.d0*PI) 
             if (a2(j+6).lt.0.d0) a2(j+6)=dmod(a2(j+6)+2.d0*PI,2.d0*PI)  
             if (a2(j+7).lt.0.d0) a2(j+7)=dmod(a2(j+7)+2.d0*PI,2.d0*PI) 
             if (a2(j+5).gt.2.d0*PI) a2(j+5)=dmod(a2(j+5),2.d0*PI)               
             if (a2(j+6).gt.2.d0*PI) a2(j+6)=dmod(a2(j+6),2.d0*PI)  
             if (a2(j+7).gt.2.d0*PI) a2(j+7)=dmod(a2(j+7),2.d0*PI)    
          enddo        
         
      endif
 
c-----get ymod first



      call MA_J_cop_fin (a2,ma,npl,mstar,sini,mass,ap,hkl)
           
      
      call GENINIT_J3_ewcop (nbod,ap,a2,
     &                          mass,xj,yj,zj,vxj,vyj,vzj,rpl,rhill,hkl)
      
      call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &                 xh,yh,zh,vxh,vyh,vzh)
      
      call integrate_cop_fin(ymod,t,nbod,na,ndata,mass,a2,
     &              xh,yh,zh,vxh,vyh,vzh,epsil,dt,hkl)
     
     
 
 

c      pause
      return
      end
      
   
c MA_J calculates the actual masses and Jacobi semimajor axes of a two-planet
c system for assumed sin(i) using the parameters P, K and e from a
c two-Kepler fit.

	    subroutine MA_J_cop_fin (a,ma,npl,m0,sini,mass,ap,hkl)
        
	    implicit none
	    real*8 m0,PI,TWOPI,THIRD,GMSUN,dm,MSUN,ecc
        integer npl,ma,i,j,NPLMAX,hkl
        parameter (NPLMAX=20)
        real*8 sini(NPLMAX)        
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

      real*8 SMASSYR,MSUN,PI,eps
      parameter (PI=3.14159265358979d0,eps=1.d-7)
      parameter (SMASSYR=4.d0*PI*PI)
      parameter (MSUN=1.32712497d20)

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
      
c      write(*,*) mass(1),mass(2),mass(3)
c      pause 
      
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
                  
c         write(*,*) ecc, omega

          gm = gm + mass(i)
          rpl(i) = frho3*(1.5d0*mass(i)/2.d0*PI)**0.3333333333d0
          rhill(i) = ap(i-1)*(mass(i)/(3.d0*mass(1)))**0.3333333333d0
    
          call ORBEL_EL2XV (gm,ialpha,ap(i-1),ecc,a(j+6),a(j+7),
     &       omega,capm,xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i))
 
 
c       write(*,*)  mass(1),mass(2),mass(3), ap(1),ap(2)
c       write(*,*)  xj(2),yj(2),zj(2),vxj(2),vyj(2),vzj(2)
c       write(*,*)  xj(3),yj(3),zj(3),vxj(3),vyj(3),vzj(3)
       


 
      enddo
c      pause
      return
      end



c      subroutine integrate_cop_fin(x,mass,ymod,dyda,t,a,ap,nbod,ma,ndata,
c     & epsil,dt,hkl)
      subroutine integrate_cop_fin(ymod,t,nbod,ma,ndata,mass,a,
     &                     xh,yh,zh,vxh,vyh,vzh,eps,dt,hkl)

	    include 'swift_loglik_Jakub.inc'
	          
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

	    real*8 t0,tstop,dt,dtout,dtdump
	    real*8 time,tout,tdump,tfrac,eoff

	    real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        logical*1 lflg(0:IO_NBITS-1)

        integer ndata,ma,i,j,nd,flag,ierr
        real*8 ymod(ndata),t(ndata),a(ma)
        real*8 ap(NPLMAX),h,eps,mtotal
        real*8 mstar,sini(NPLMAX)
         

	    character*80 outfile,inparfile,inplfile,intpfile,fopenstat



        ntp = 0

c Get data for the run and the test particles
c      call io_open(7,'dynfit_param.in','old','formatted',ierr)

 
      
      
      

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
c        write(*,*) dt, eps
c        write(*,*) a(1),a(2),a(3),a(4),a(5),a(6),a(7)
c        write(*,*) a(8),a(9),a(10),a(11),a(12),a(13),a(14) 
c        write(*,*) a(15)
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
           h = dt
           flag = 0        ! flag for controling output
           do i = 1,ndata    
              if ((time.lt.t(i)).and.((t(i)-time).le.dt)) then
                 flag = 1
                 h = t(i) - time
                 goto 555
              endif
           enddo

555       call bs_step2(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &         xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &         vzht,istat,rstat,h,eps)
           
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
        






      subroutine bs_step2(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt,eps)	

      include 'swift_loglik_Jakub.inc'
      include 'bs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4

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






