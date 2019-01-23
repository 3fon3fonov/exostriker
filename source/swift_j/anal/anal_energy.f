c*************************************************************************
c                          ANAL_ENERGY.F
c*************************************************************************
c Calculates the energy of the total system (massive bodies) wrt time.
c returns the total energy of n objects by direct pairwise summation
c G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
c
c      Input:
c            t             ==>  current time
c            nbod          ==>  number of massive bodies (int scalar)
c            mass          ==>  mass of bodies (real array)
c            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
c            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
c            xh,yh,zh      ==>  current position in heliocentric coord 
c                               (real arrays)
c            vxh,vyh,vzh   ==>  current velocity in heliocentric coord 
c                               (real arrays)
c
c      Output:
c            ke            ==>  kinetic energy
c            pot           ==>  potential energy
c            energy        ==>  Total energy
c            eltot         ==>  components of total angular momentum
c                               (real array)
c
c Remarks: 
c Authors:  Martin Duncan
c Date:  ?
c Last revision:  1/24/97 HFL 

      subroutine anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &           vxh,vyh,vzh,ke,pot,energy,eltot)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Output
      real*8 energy,eltot(3),ke,pot

c...  Internals
      real*8 elx,ely,elz
      real*8 xx,yy,zz,rr2,oblpot,msys,irh(NTPMAX),ir3h(NTPMAX)
      real*8 xb(NTPMAX),yb(NTPMAX),zb(NTPMAX)       ! Used NTPMAX for symba
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)
      integer i,j

c----
c...  Executable code 

      call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &           xb,yb,zb,vxb,vyb,vzb,msys)   

      eltot(1)=(yb(nbod)*vzb(nbod)-zb(nbod)*vyb(nbod))*mass(nbod)
      eltot(2)=(zb(nbod)*vxb(nbod)-xb(nbod)*vzb(nbod))*mass(nbod)
      eltot(3)=(xb(nbod)*vyb(nbod)-yb(nbod)*vxb(nbod))*mass(nbod)

      ke = 0.5*mass(nbod)*(vxb(nbod)**2 + vyb(nbod)**2 + vzb(nbod)**2)
      pot= 0.d0

      do i = 1,nbod-1

         elx=(yb(i)*vzb(i)-zb(i)*vyb(i))*mass(i)
         ely=(zb(i)*vxb(i)-xb(i)*vzb(i))*mass(i)
         elz=(xb(i)*vyb(i)-yb(i)*vxb(i))*mass(i)
         eltot(1) = eltot(1) + elx
         eltot(2) = eltot(2) + ely
         eltot(3) = eltot(3) + elz
         
         ke = ke + 0.5*mass(i)*(vxb(i)**2 + vyb(i)**2 + vzb(i)**2)
         do j = i+1,nbod

	    xx = xb(i) - xb(j)
	    yy = yb(i) - yb(j)
	    zz = zb(i) - zb(j)
	    rr2 = xx**2 + yy**2 + zz**2 
            if((mass(i).ne.0.0d0).and.(mass(j).ne.0.0d0)) then
               pot = pot - mass(i)*mass(j)/(sqrt(rr2))
            endif
         enddo
      enddo

      if(j2rp2.ne.0.0d0) then
         call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
         call obl_pot(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &                oblpot)
         pot = pot + oblpot
      endif

      energy = ke + pot

      return	
      end      ! anal_energy
c-----------------------------------------------------------------------

