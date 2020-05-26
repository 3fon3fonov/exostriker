c*************************************************************************
c                            RMVS3_CHK.F
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
c                 rts           ==>  The fraction of hill sphere for encounter
c                                               (real scalor)
c             Output:
c                 icflg         ==> ecounters? = 1 Yes
c                                              =  0 No (integer scalar)  
c                 nenc          ==> nenc(i) is the number of tp enc planet i
c                                   (integer array)
c                 itpenc        ==> itpenc(*,i) is a list of tp enc planet i
c                                   (2d integer array)
c                 ienc          ==> ienc(j) = 0 if tp j not involved in enc 
c                                           = planet# if it is. 
c                                           (integer array)
c
c Remarks: Based on RMVS_CHK.F
c Authors:  Hal Levison 
c Date:    7/10/96
c
c Last revision: 

      subroutine rmvs3_chk(nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,
     &       zht,vxht,vyht,vzht,istat,dt,rts,icflg,nenc,itpenc,ienc)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs: 
      integer nbod,ntp,istat(ntp)
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod),dt
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 rts

c...  Outputs
       integer icflg
       Integer nenc(NPLMAX),itpenc(NTPMAX,NPLMAX)
       integer ienc(NTPMAX)

c...  Internals
       integer i,j,iflag
       real*8 r2hill(NPLMAX),r2crit,r2critp
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
        do i=1,nbod
           nenc(i) = 0
           do j=1,ntp
              itpenc(j,i) = 0
           enddo
        enddo

        do i=1,ntp
           ienc(i) = 0
        enddo

        do j=1,ntp
           if(istat(j).eq.0) then

	      i = 2
	      iflag = 0				! precaution
 
c... Check for close encounters until we find one or run out of planets. BG

              do while ( (iflag .eq. 0) .and. (i .le. nbod) )

                 r2crit = r2hill(i)*rts
                 r2critp = 0.0           ! dummy so we can use rmvs routine

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
                    ienc(j) = i
                    nenc(i) = nenc(i) + 1
                    itpenc(nenc(i),i) = j
                 endif
		 
 		 i = i + 1    		! next planet
              enddo
           endif
        enddo

        return
        end  ! rmvs3_chk
c------------------------------------------------------

