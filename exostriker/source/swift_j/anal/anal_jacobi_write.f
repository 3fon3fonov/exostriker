c*************************************************************************
c                          ANAL_JACOBI_WRITE.F
c*************************************************************************
c Writes the mean and maximum absolute value of the change in jacobi
c to the screen as a function of time.  Writes the value of the first 10 
c test particles to a file called jacobi.out (unit=iu)
c
c      Input:
c            t              ==>  current time
c            nbod           ==>  number of massive bodies (int scalar)
c            ntp            ==>  number of tp (int scalar)
c            mass           ==>  mass of bodies (real array)
c            xh,yh,zh       ==>  current position in helio coord 
c                               (real arrays)
c            vxh,vyh,vzh    ==>  current velocity in helio coord 
c                               (real arrays)
c            xht,yht,zht    ==>  current tp position in helio coord 
c                               (real arrays)
c            vxht,vyht,vzht ==>  current tp velocity in helio coord 
c                               (real arrays)
c            istat          ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c            ipl            ==>  Planet to take jacobi with respect to 
c            iu             ==>  unit to write to
c            fopenstat      ==>  The status flag for the open 
c                                statements of the output files.  
c                                          (character*80)
c
c Remarks: If the particle is not active a value of -100 is written out 
c Authors:  Hal Levison 
c Date:    3/4/93
c Last revision: 10/3/96

      subroutine anal_jacobi_write(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &    vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,ipl,iu,fopenstat)


      include '../swift.inc'

c...  Inputs: 
      integer nbod,ntp,ipl,iu
      integer istat(ntp)
      real*8 mass(nbod),t
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      character*80 fopenstat

c...  Internals
      integer i1st,i,icnt,nw
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
      real*8 jac,jac0(NTPMAX),djmean,djmax,dj(NTPMAX)
      real*8 gmsum,energy,aplh,omega,fac
      real*8 omegax,omegay,omegaz,msys

      data i1st/0/
      save i1st,jac0

c----
c...  Executable code 

      nw = min0(ntp,10)

c...   Compute ang. mom. vector for Sun-planet relative orbit
      gmsum = mass(1) + mass(ipl)
      energy = 0.5d0*(vxh(ipl)**2 + vyh(ipl)**2 + vzh(ipl)**2) 
      energy = energy - gmsum/sqrt(xh(ipl)**2 + yh(ipl)**2 + zh(ipl)**2)
      aplh = -0.5d0*gmsum/energy
      omega = sqrt(gmsum/(aplh**3))
      omegax = yh(ipl)*vzh(ipl) - zh(ipl)*vyh(ipl)
      omegay = zh(ipl)*vxh(ipl) - xh(ipl)*vzh(ipl)
      omegaz = xh(ipl)*vyh(ipl) - yh(ipl)*vxh(ipl)
      fac = omega/sqrt(omegax**2 + omegay**2 + omegaz**2)
      omegax = fac*omegax
      omegay = fac*omegay
      omegaz = fac*omegaz

c...  put things in bary
      call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &           xb,yb,zb,vxb,vyb,vzb,msys)   

      call coord_h2b_tp(ntp,xht,yht,zht,vxht,vyht,vzht,
     &      xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &      xbt,ybt,zbt,vxbt,vybt,vzbt)


      if(i1st.eq.0) then

         do i=1,ntp
            call anal_jacobi(mass(1),mass(ipl),omegax,omegay,omegaz,
     &           xbt(i),ybt(i),zbt(i),vxbt(i),vybt(i),vzbt(i),xb(1),
     &           yb(1),zb(1),xb(ipl),yb(ipl),zb(ipl),jac0(i))
            dj(i) = 0.0
         enddo

         call io_jacobi_write(i1st,t,jac0,dj,nw,iu,fopenstat)

         i1st = 1

      else

         icnt = 0
         djmean = 0.0
         djmax = 0.0
         do i=1,ntp
            if(istat(i).eq.0) then
               call anal_jacobi(mass(1),mass(ipl),omegax,omegay,omegaz,
     &              xbt(i),ybt(i),zbt(i),vxbt(i),vybt(i),vzbt(i),xb(1),
     &              yb(1),zb(1),xb(ipl),yb(ipl),zb(ipl),jac)
               icnt = icnt + 1
               dj(i) = jac/jac0(i) - 1.d0
               djmean = djmean + abs(dj(i))
               djmax = dmax1(djmax,abs(dj(i)))
            else
               dj(i) = -100.0
            endif
         enddo
         djmean = djmean/float(icnt)
         write(*,1) djmean,djmax
 1       format(5x,'mean |dj/j|, max |dj/j|,',2(2x,1p1e12.5))

         call io_jacobi_write(i1st,t,jac0,dj,nw,iu,fopenstat)
         
      endif

      return
      end                       ! anal_jacobi_write
c-------------------------------------------------------------------------
