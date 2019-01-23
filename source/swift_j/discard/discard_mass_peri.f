c*************************************************************************
c                            DISCARD_MASS_PERI.F
c*************************************************************************
c This subroutine checks to see if a partical should be discarded because
c of its perihelion distance gets too small
c
c             Input:
c                 time           ==>  current time (real scalar)
c                 nbod           ==>  number of test bodies (int scalar)
c                 iecnt          ==>  Number of encounters (int*2 array)
c                 mass           ==>  mass of bodies (real array)
c                 xh,yh,zh       ==>   part position in helio coord 
c                                      (real arrays)
c                 vxh,vyh,vzh    ==>   part vel in helio coord 
c                                      (real arrays)
c                 qmin           ==>  Smallest perihelion distance 
c                                      (real scalar)
c                 iwhy           ==>  status of the object
c                                      (integer array)
c             Output:
c                 iwhy           ==>  status of the object
c                                      (integer array)
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c Remarks: Based on discard_peri
c Authors:  Hal Levison 
c Date:    12/30/96
c Last revision: 

      subroutine discard_mass_peri(time,nbod,iecnt,mass,xh,yh,zh,
     &     vxh,vyh,vzh,qmin,iwhy,isperi)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),time,qmin
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      integer*2 iecnt(NTPMAX)

c...  Input and Output
      integer iwhy(nbod)
      integer isperi(nbod)

c...  internal
      integer i,i1st
      real*8 peri(NTPMAX)
      logical*2 lperi(NTPMAX)

      data i1st/0/

      save i1st

c-----
c...  Executable code 

      if(i1st.eq.0) then     ! if first time through, set things up
         call util_mass_peri(0,nbod,xh,yh,zh,vxh,vyh,vzh,
     &     mass,isperi,peri,lperi)
         i1st = 1
         return                 !  <==== RETURN
      endif

      call util_mass_peri(1,nbod,xh,yh,zh,vxh,vyh,vzh,
     &     mass,isperi,peri,lperi)

      do i=2,nbod
         if((isperi(i).eq.0) .and. (iecnt(i).eq.0)) then
            if(peri(i).le.qmin) then
               write(*,*) 'Particle',i,' perihelion distance too',
     &               ' small at t=',time
               iwhy(i) = -4
            endif
         endif
      enddo

      return
      end       ! discard_mass_peri
c------------------------------------------------------










