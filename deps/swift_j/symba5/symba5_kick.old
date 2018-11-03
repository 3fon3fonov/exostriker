c*************************************************************************
c                             SYMBA5_KICK.F
c*************************************************************************
c Do a symba5 kick
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 irec          ==>  recursion level  (integer scalar)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 ielev         ==>  The level that this particle should go
c                                             (int*2 array)
c                 rhill         ==>  Hill sphere of planet (real Scalar)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial velocity in bari coord 
c                                    (real arrays)
c                dt             ==>  timestep  (real scalar)
c                sgn            ==>  add or subtract the force (real scalar)
c                ielc           ==>  number of encounters (integer*2 scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c            Output:
c                 vxb,vyb,vzb   ==>  final velocity in bari coord 
c                                    (real arrays)
c
c Remarks: Uses Man Hoi's force
c Authors:  Hal Levison 
c Date:   3/20/97
c Last revision: 

      subroutine symba5_kick(nbod,mass,irec,iecnt,ielev,
     &        rhill,xh,yh,zh,vxb,vyb,vzb,dt,sgn,ielc,ielst)

      include '../swift.inc'
      include 'symba5.inc'

c...  Inputs Only: 
      integer nbod,irec
      real*8 mass(nbod),dt,rhill(nbod),sgn
      integer*2 iecnt(NTPMAX),ielev(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      integer*2 ielst(2,NENMAX),ielc

c...  Inputs and Outputs:
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Internals: 
      real*8 ax(NTPMAX),ay(NTPMAX),az(NTPMAX)
      real*8 dx,dy,dz,fac,ris,r
      real*8 ri,rr,r2,faci,facj,ir3,rim1
      integer i,j,irm1,irecl,ie

c----
c...  Executable code 

      irm1 = irec - 1
      if(sgn.lt.0.0d0) then
         irecl = irec - 1
      else
         irecl = irec
      endif

c...  Zero everthing
      do i=1,nbod
         ax(i) = 0.0d0
         ay(i) = 0.0d0
         az(i) = 0.0d0
      enddo

c...  calculate the accelerations

      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)

         if((ielev(i).ge.irm1) .and. (ielev(j).ge.irm1) ) then

            ri = (rhill(i)+rhill(j))**2 * RHSCALE*RHSCALE * 
     &           (RSHELL**(2*irecl))
            rim1 = ri*RSHELL*RSHELL
            
            dx = xh(j) - xh(i)
            dy = yh(j) - yh(i)
            dz = zh(j) - zh(i)
            r2 = dx*dx + dy*dy + dz*dz
            ir3 = 1.0d0/(r2*sqrt(r2))

            if (r2.lt.rim1) then
               fac = 0.0d0
            else if (r2.lt.ri) then
               ris = sqrt(ri)
               r = sqrt(r2)
               rr = (ris-r)/(ris*(1.0-RSHELL))
               fac = (r2**(-1.5d0)) * 
     &              ( 1.0d0 - 3.0d0*rr*rr + 2.0d0*(rr**3))
            else
               fac = ir3
            endif
            
            faci = mass(i)*fac
            facj = mass(j)*fac
            
            ax(j) = ax(j) - faci*dx
            ay(j) = ay(j) - faci*dy
            az(j) = az(j) - faci*dz
            
            ax(i) = ax(i) + facj*dx
            ay(i) = ay(i) + facj*dy
            az(i) = az(i) + facj*dz
            
         endif

      enddo
      
c...  apply the kick

      do i=2,nbod
         if( (iecnt(i).ne.0) .and. (ielev(i).ge.irm1) ) then
            vxb(i) = vxb(i) + ax(i)*dt*sgn
            vyb(i) = vyb(i) + ay(i)*dt*sgn
            vzb(i) = vzb(i) + az(i)*dt*sgn
         endif
      enddo

      return
      end      ! symba5_kick.f
c--------------------------------------------------------------
