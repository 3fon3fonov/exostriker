c*************************************************************************
c                            RMVS_CHK.F
c*************************************************************************
c This subroutine checks to see if there are encounters
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c                 dt            ==>  time step  (real scalor)
c             Output:
c                 icflg         ==> ecounters? = 1 Yes, in outer region only
c                                              = -1 in inner region
c                                              =  0 No (integer scalar)  
c                 nenco         ==> nenco(i) is the number of tp enc planet i
c                                   in outer region (integer array)
c                 nenci         ==> nenci(i) is the number of tp enc planet i
c                                   in inner region (integer array)
c                 itpenco       ==> itpenco(*,i) is a list of tp enc planet i
c                                   in outer region (2d integer array)
c                 itpenci       ==> itpenci(*,i) is a list of tp enc planet i
c                                   in inner region (2d integer array)
c                 ienco         ==> ienco(j) = 0 if tp j not involved in enc 
c                                   in outer region: = planet# if it is. 
c                                     (integer array)
c                 ienci         ==> same but for inner region.
c
c
c Remarks: Based on Hal's wiscl_fk.f
c Authors:  Hal Levison 
c Date:    2/19/93
c
c Last revision: 7/14/94,  Brett Gladman,  main loop changed to "do while" to 
c                                          improved efficiency.
c 

      subroutine rmvs_chk(nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,
     &       zht,vxht,vyht,vzht,istat,dt,icflg,nenco,itpenco,nenci,
     &       itpenci,ienci,ienco)


      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs: 
      integer nbod,ntp,istat(ntp)
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod),dt
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Outputs
       integer icflg
       Integer nenco(NPLMAX),itpenco(NTPMAX,NPLMAX),nenci(NPLMAX)
       integer itpenci(NTPMAX,NPLMAX),ienci(NTPMAX),ienco(NTPMAX)

c...  Internals
       integer i,j,iflag,icflgi
       real*8 r2hill(NPLMAX),rts,rps,r2crit,r2critp
       real*8 xr,yr,zr,vxr,vyr,vzr

       integer i1st    ! =0 first time through; =1 after
       data i1st/0/
       save i1st,r2hill

c-----
c...  Executable code 

c...  if first time through, calc hill's shere for the planets
        if(i1st.eq.0) then
           call util_hills(nbod,mass,xh,yh,zh,vxh,vyh,vzh,r2hill)
           i1st = 1
        endif

c...    clear everything out
        icflg = 0
        icflgi = 0
        do i=1,nbod
           nenco(i) = 0
           nenci(i) = 0
           do j=1,ntp
              itpenco(j,i) = 0
              itpenci(j,i) = 0
           enddo
        enddo

        do i=1,ntp
           ienco(i) = 0
           ienci(i) = 0
        enddo

        rts = RHSCALE*RHSCALE
        rps = RHPSCALE*RHPSCALE
        
        do j=1,ntp
           if(istat(j).eq.0) then

	      i = 2
	      iflag = 0				! precaution
 
c... Check for close encounters until we find one or run out of planets. BG

              do while ( (iflag .eq. 0) .and. (i .le. nbod) )

                 r2crit = r2hill(i)*rts
                 r2critp = r2hill(i)*rps

                 xr = xht(j) - xh(i)
                 yr = yht(j) - yh(i)
                 zr = zht(j) - zh(i)
                 vxr = vxht(j) - vxh(i)
                 vyr = vyht(j) - vyh(i)
                 vzr = vzht(j) - vzh(i)
                 call rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,
     &                         r2crit,r2critp,iflag)
                 if(iflag.gt.0) then
                    icflg = 1
                    ienco(j) = i
                    nenco(i) = nenco(i) + 1
                    itpenco(nenco(i),i) = j
                 endif
                 if(iflag.lt.0) then
                    icflgi = 1
                    ienci(j) = i
                    nenci(i) = nenci(i) + 1
                    itpenci(nenci(i),i) = j
                 endif
		 
 		 i = i + 1    		! next planet
              enddo
           endif
        enddo

        if(icflgi.eq.1) then
           icflg = -1
        endif


        return
        end  ! rmvs_chk
c------------------------------------------------------

