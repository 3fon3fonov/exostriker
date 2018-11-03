c**********************************************************************
c		      SWIFT_MVS.F
c**********************************************************************
c
c                 NO CLOSE ENCOUNTERS
c                 To run, need 3 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c                   test particle file like   tp.in
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    5/7/93
c Last revision: 12/27/96

     
	include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT),i1st
	integer nbod,ntp,nleft,tidal_eff,z,i
	integer iflgchk,iub,iuj,iud,iue
        real*8 rstat(NTPMAX,NSTATR)

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tout,tdump,tfrac,eoff

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose 

	character*80 outfile,inparfile,inplfile
	character*80 intpfile,fopenstat, ssefile
        real*8 gmi,ai,ei,inci,capomi,omegai,capmi
        real*8 dadt(NTPMAX),dedt(NTPMAX),indays

c      real*8 or_dadt(NTPMAX),or_dedt(NTPMAX)
        real*8 or_dadt,or_dedt

        real*8 Mt(2000),t_sse(2000),L(2000)
        real*8 R(2000),Menv(2000),spin(2000),Renv(2000)

        integer steps, ll, p
        real*8  dM, dL, dR, dMenv,dRenv,dspin,dM1
        real*8  L_i, R_i, Menv_i, Renv_i, spin_i, M_i
      
        integer ialphai
 

c-----
c...    Executable code 

c...    print version number
        call util_version

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

        write(*,*) ' '
        write(*,*) 'Enter name of SSE .dat data file : '
        read(*,999) ssefile

        write(*,*) 'Tidal effects? 0 - no, 1 - yes :'
        read(*,*) tidal_eff
        write(*,*) ' tidal_eff = ',tidal_eff


c Initialize initial time and times for first output and first dump
	t = t0
	tout = t0 + dtout
	tdump = t0 + dtdump

        iub = 20
        iuj = 30
        iud = 40
        iue = 60

c...    Do the initial io write
        if(btest(iflgchk,0))  then ! bit 0 is set
           call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,1))  then ! bit 1 is set
           call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,2))  then    ! bit 2 is set
           eoff = 0.0d0
           call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &          vyh,vzh,iue,fopenstat,eoff)
        endif
        if(btest(iflgchk,3))  then    ! bit 3 is set
           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
        endif

c...    must initize discard io routine
        if(btest(iflgchk,4))  then ! bit 4 is set
           call io_discard_write(0,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &          vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &          'discard.out',fopenstat,nleft)
        endif

        nleft = ntp
	i1st = 0

        call ess_read(ssefile, t_sse, Mt,L,R,Menv,Renv,spin)

        z = 1
        dM1 = 0.d0 
      
      
        mass(1) = Mt(z)*39.478417604357354d0
      
        M_i  = Mt(z)          
        L_i  = L(z)               
        R_i  = R(z)  
        Menv_i = Menv(z)                         
        Renv_i = Renv(z)
        spin_i = spin(z) 


 

        p = 0
        ll = 1000
        indays = 365.25*365.25

c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '

	  do while ( (t .le. tstop) .and. 
     &       ((ntp.eq.0).or.(nleft.gt.0)) )


	 
c        ************* the mass loss (from SSE!!!) ****************************
             if(t.ge.(t_sse(z)-t_sse(1))*1000000.0d0) then

              steps=int(abs((( t_sse(z+1)-t_sse(z))*1000000.0d0)/dt))- 1
              
              
              dM1 = ( Mt(z+1) - Mt(z) )*39.478417604357354d0 / steps
              dM  = ( Mt(z+1) - Mt(z) ) / steps              
              dL  = (  L(z+1) - L(z) ) / steps              
              dR  = (  R(z+1) - R(z) ) / steps  
              dMenv = ( Menv(z+1) - Menv(z) ) / steps                          
              dRenv = ( Renv(z+1) - Renv(z) ) / steps 
              dspin = ( spin(z+1) - spin(z) ) / steps 
              
                                                     
 
              mass(1) = Mt(z)*(39.478417604357354d0/indays)
	      M_i  = Mt(z)          
              L_i  = L(z)               
              R_i  = R(z)  
              Menv_i = Menv(z)                         
              Renv_i = Renv(z)
              spin_i = spin(z)
              
              z = z+1              
c              write(*,*) or_dadt,or_dedt, ai,ei,mass(2)
              write(*,*) t_sse(z),mass(1), M_i,L_i, R_i,Menv_i,Renv_i
c              write(*,*)t_sse(z), dM1,dM,dL,dR,dMenv,dRenv,  steps  
 
c              write(*,*) steps
c              write(*,*)z,t_sse(z),mass(1),mass(1)/39.478417604357354d0,
c     &        Mt(z),dM, mass(1)/39.478417604357354d0 - Mt(z)
c              write(*,*) or_dadt,or_dedt

c         endif

             else
              mass(1) = mass(1) + dM1
         
	      M_i  = M_i + dM         
              L_i  = L_i + dL            
              R_i  = R_i + dR 
              Menv_i = Menv_i + dMenv                        
              Renv_i = Renv_i + dRenv
              spin_i = spin_i + dspin
c               write(*,*) t_sse(z),mass(1), M_i,L_i, R_i,Menv_i,Renv_i
c               write(*,*)t_sse(z), dM1,dM,dL,dR,dMenv,dRenv,  steps 

	     endif              



c         pause
 	     
	     if((tidal_eff.eq.1).and.(p.eq.ll)) then 
c	 if(tidal_eff.eq.1) then 	     
               do i=2,nbod
                gmi = gmi + mass(i)
            
                call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),
     &           gmi,ialphai,ai,ei,inci,capomi,omegai,capmi)           
            
c            write(*,*) ai,ei, gmi
            
                call orb_retard2 (mass(i),ai,ei,M_i, L_i,
     &      R_i, Renv_i, Menv_i, spin_i, 
     &      dt, or_dadt,or_dedt)
     
                ai = ai + or_dadt*ll
                ei = ei + or_dedt*ll
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c this is just workaround a symba problem when one of the 
c bodies is lost - it must be corrected at some point!

                if(ei.lt.0.d0) then
                  ei = 0.0001d0
                endif
                if((ai.lt.0.d0).or.(ei.lt.0.d0)) then
                          
                  ei = 0.0001d0
                  ai = 0.0001d0
                  inci = 0.0d0
                  capomi = 0.0d0
                  omegai = 0.0d0
                  capmi = 0.0d0
              
                 write(*,*) "a",i," < 0 !"
                endif         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            write(*,*) or_dadt,or_dedt
                call orbel_el2xv(gmi,ialphai,ai,ei,inci,capomi,omegai,
     &           capmi,xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i))
 
c            pause
               enddo   
               p = 0        
             endif

             p = p + 1
                     






 	     call step_kdk(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &         xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &         vzht,istat,rstat,dt)

	     t = t + dt

             if(btest(iflgchk,4))  then    ! bit 4 is set
                call discard(t,dt,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &               xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
     &               qmin,lclose,rplsq,istat,rstat)
                call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &               vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &               'discard.out',fopenstat,nleft)
             else
                nleft = ntp
             endif

c if it is time, output orb. elements, 
	  if(t .ge. tout) then 

             if(btest(iflgchk,0))  then    ! bit 0 is set
                call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &               iub,fopenstat)
              endif
             if(btest(iflgchk,1))  then    ! bit 1 is set
                call  io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &               iub,fopenstat)
              endif

	    tout = tout + dtout
	  endif

c If it is time, do a dump
          if(t.ge.tdump) then

             tfrac = (t-t0)/(tstop-t0)
             write(*,998) t,tfrac,nleft
 998         format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of active tp =',i4)
	     call io_dump_pl('dump_pl.dat',nbod,mass,xh,yh,zh,
     &                 vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	     call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     &                      vxht,vyht,vzht,istat,rstat)
	     call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &           dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
	     tdump = tdump + dtdump

             if(btest(iflgchk,2))  then    ! bit 2 is set
                call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &               xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
             endif
             if(btest(iflgchk,3))  then    ! bit 3 is set
                call anal_jacobi_write(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,
     &               iuj,fopenstat)
             endif

	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

	call io_dump_pl('dump_pl.dat',nbod,mass,xh,yh,zh,
     &            vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     &              vxht,vyht,vzht,istat,rstat)
	call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

        call util_exit(0)
        end    ! swift_mvs.f
c---------------------------------------------------------------------


c reads the SSE data file T.T.

      subroutine ess_read (ssefile, t,Mt,L, R,Menv,Renv,spin)
 
      character*(*) ssefile
      real*8 t(2000),typ(2000),Mo(2000),Mt(2000), L(2000),R(2000),
     &    Teff(2000),Mc(2000),Menv(2000),Renv(2000),
     &    epoch(2000),spin(2000)
      integer n
      
      n = 1
      
      
      open (unit=15, file=ssefile, status='old',
     &    access='sequential', form='formatted', action='read' )
 100  continue
 
 
         read (15, *,err=200,end=200) t(n),typ(n),Mo(n),Mt(n),
     &    L(n),R(n),
     &    Teff(n),Mc(n),Menv(n),Renv(n),epoch(n),spin(n)
     
 110     format (12(F12.8) )
 

 120     format (12(F12.3) )
 
         n = n + 1
 
         goto 100
         close (unit=15)
         
         
 200  return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
      subroutine  orb_retard (pl_mass_ii,semiM_ii, 
     &    ecc_ii, Mt_ii, L_ii,R_ii,Renv_ii,
     &    Menv_ii,spin_ii,dt, or_dadt,or_dedt)

c      include 'swift.inc'

      real*8 pl_mass,semiM, ecc, T,  prec_frac
      real*8 or_dadt,or_dedt, t_conv, R_env
      real*8 Mt_iii,L_iii,R_iii, Menv_iii,spin_iii
     &    pe_atail,pe_etail
      real*8 Mt_ii,L_ii,R_ii, Menv_ii,spin_ii,Renv_ii    
      real*8 dt,pl_mass_ii,Mrat,Rrat,semiM_ii, ecc_ii   
      real*8 pe(3), MSUN, RSUN, LSUN, GMSUN2

      parameter (PI=3.14159265358979d0,TWOPI=2.d0*PI,THIRD=1.d0/3.d0)
      parameter (GM=6.674d-11,AU=149597870700.d0) 
      parameter (RSUN=6.95700d8, SMASSYR=TWOPI*TWOPI)
      parameter (MSUN=1.98855d30, LSUN = 3.828d26 )

      integer j
    
      L_iii     = L_ii * LSUN

      R_iii     = R_ii * RSUN
      Renv_iii  = Renv_ii * RSUN
      
      pl_mass   = (pl_mass_ii/SMASSYR) 
      
      Menv_iii  = Menv_ii * MSUN
      Mt_iii    = Mt_ii * MSUN

      semiM     = semiM_ii * AU
      ecc       = ecc_ii
      

      t_conv  = ( (Menv_ii * (R_ii - R_env_ii)**2.0d0)/
     &   ( 3.0d0*L_ii) )**THIRD
     
      Mrat = (Menv_ii/Mt_ii) * (pl_mass/Mt_ii) *
     &            (1.d0 + (pl_mass / Mt_ii))
     
      Rrat = (R_iii / semiM)**8.0d0
      
c      Rrat = (R_ii / (semiM_ii/0.004650467d0))**8.0d0      

      do j=1,3        
        pe(j) = (4.d0*(PI*PI)*semiM**3.d0) / 
     &  (((j**2.0d0)*(GM*(Mt_iii+(pl_mass*MSUN))))*(t_conv**2.d0))

 
        pe(j) = 9.d0/2.d0 * dmin1(1.d0,pe(j)) 
      enddo

      pe_atail = 2.d0*pe(2) + (ecc*ecc)*( (7.d0/8.d0)*pe(1)
     &         - 10.d0*pe(2) + (441.d0/8.d0)*pe(3) )   

      pe_etail = ( (5.d0/4.d0)*pe(1) - 2.d0*pe(2) +
     &    (147.d0/4.d0)*pe(3) ) 
 
      or_dadt =  (((semiM_ii ) / (9.d0 * t_conv)) * Mrat * Rrat *
     &           pe_atail)*(-1.d0)
     
      or_dedt = ( (ecc_ii  / (36.d0 * t_conv)) * Mrat * Rrat * 
     &           pe_etail)*(-1.d0)    
     
     
      T = 2.0d0*PI * sqrt((semiM**3.0d0)/(GM*(Mt_iii+pl_mass*MSUN) ) ) 

      prec_frac =  (T /(86400.0d0 *365.25))  / dt  
 
      or_dadt =   or_dadt / prec_frac
      or_dedt =   or_dedt / prec_frac  
 
      return
      end


 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
      subroutine  orb_retard2 (pl_mass_ii,semiM_ii, 
     &    ecc_ii, Mt_ii, L_ii,R_ii,Renv_ii,
     &    Menv_ii,spin_ii,dt, or_dadt,or_dedt)

c      include 'swift.inc'

      real*8 pl_mass,semiM, ecc, T,  prec_frac
      real*8 or_dadt,or_dedt, t_conv, R_env
      real*8 Mt_iii,L_iii,R_iii, Menv_iii,spin_iii
     &    pe_atail,pe_etail
      real*8 Mt_ii,L_ii,R_ii, Menv_ii,spin_ii,Renv_ii    
      real*8 dt,pl_mass_ii,Mrat,Rrat,semiM_ii, ecc_ii   
      real*8 pe(3), MSUN, RSUN, LSUN, GMSUN2, meanM

      parameter (PI=3.14159265358979d0,TWOPI=2.d0*PI,THIRD=1.d0/3.d0)
      parameter (GM=6.674d-11,AU=149597870700.d0) 
      parameter (RSUN=6.95700d8, SMASSYR=TWOPI*TWOPI)
      parameter (MSUN=1.98855d30, LSUN = 3.828d26 )

      integer j
    
      L_iii     = L_ii * LSUN

      R_iii     = R_ii * RSUN
      Renv_iii  = Renv_ii * RSUN
      
      pl_mass   = (pl_mass_ii/SMASSYR) 
      
      Menv_iii  = Menv_ii * MSUN
      Mt_iii    = Mt_ii * MSUN

      semiM     = semiM_ii * AU
      ecc       = ecc_ii
      

      t_conv  = ( (Menv_ii * (R_ii - R_env_ii)**2.0d0)/
     &   ( 3.0d0*L_ii) )**THIRD
     
      Mrat = (Menv_ii/Mt_ii) * (pl_mass/Mt_ii) *
     &            (1.d0 + (pl_mass / Mt_ii))
     
      Rrat = (R_iii / semiM)**8.0d0
      
c      Rrat = (R_ii / (semiM_ii/0.004650467d0))**8.0d0  



      meanM    = sqrt((SMASSYR*(Mt_ii+pl_mass))/semiM_ii**3.0d0)

      do j=1,3        
        pe(j) = (TWOPI / (j*meanM*t_conv))**2.d0

 
        pe(j) = 9.d0/2.d0 * dmin1(1.d0,pe(j)) 
      enddo

      pe_atail = 2.d0*pe(2) + (ecc*ecc)*( (7.d0/8.d0)*pe(1)
     &         - 10.d0*pe(2) + (441.d0/8.d0)*pe(3) )   

      pe_etail = ( (5.d0/4.d0)*pe(1) - 2.d0*pe(2) +
     &    (147.d0/4.d0)*pe(3) ) 
 
      or_dadt =  (((semiM_ii ) / (9.d0 * t_conv)) * Mrat * Rrat *
     &           pe_atail)*(-1.d0)
     
      or_dedt = ( (ecc_ii  / (36.d0 * t_conv)) * Mrat * Rrat * 
     &           pe_etail)*(-1.d0)    
     
c      write(*,*) or_dadt, or_dedt
c      T = 2.0d0*PI * sqrt((semiM**3.0d0)/(GM*(Mt_iii+pl_mass*MSUN) ) ) 

c      write(*,*) T/86400.0

      T = 2.0d0*PI * sqrt((semiM_ii**3.0d0)/(SMASSYR*(Mt_ii+pl_mass)) ) 
      
      
c      write(*,*) or_dadt, or_dedt , T
c      pause
c      prec_frac =  (T /(86400.0d0 *365.25)) / dt  
 
c      or_dadt =   (or_dadt / (prec_frac ))
 
c      or_dedt =   (or_dedt / (prec_frac ))   
      
      or_dadt =   or_dadt/ (T * dt) 
      or_dedt =   or_dedt/ (T * dt)     
      
 
      return
      end




