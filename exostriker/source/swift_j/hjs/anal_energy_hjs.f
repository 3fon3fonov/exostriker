c*************************************************************************
c                          ANAL_ENERGY_HJS.F
c*************************************************************************
c Returns the energy of the total system (massive bodies) in
c   the hierarchical Jacobi case
c
c      Input:
c            nbod          ==>  number of massive bodies (int scalar)
c            mass          ==>  mass of bodies (real array)
c            umat          ==>  Conversion matrix Jacobi => Barycentric
c                                (2D real array)
c            xj,yj,zj      ==>  current position in Jacobi coord 
c                               (real arrays)
c            vxj,vyj,vzj   ==>  current velocity in Jacobi coord 
c                               (real arrays)
c            eoff          ==> An energy offset that is added to the energy
c                                      (real*8 scalar)
c      Output:
c            energy        ==> the energy
c
c Remarks: Adapted from anal_energy_write_hjs
c Authors:  Herve Beust
c Date:    Sep 4, 2006

      subroutine anal_energy_hjs(nbod,umat,mass,xj,yj,zj,
     &                                       vxj,vyj,vzj,eoff,energy)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,iu
      real*8 mass(nbod),t,umat(NPLMAX,NPLMAX),eoff
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Internals
      integer i1st
      real*8 energy,ke,pot
      real*8 xx,yy,zz,rr2
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      integer i,j

c----
c...  Executable code 

c Compute initial ke,pot,energy

      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xb,yb,zb,vxb,vyb,vzb) 

      ke = 0.5d0*mass(nbod)*(vxb(nbod)**2+vyb(nbod)**2+vzb(nbod)**2)
      pot = 0.d0

      do i = 1,nbod-1

         ke = ke + 0.5d0*mass(i)*(vxb(i)**2+vyb(i)**2+vzb(i)**2)
         do j = i+1,nbod

	    xx = xb(i) - xb(j)
	    yy = yb(i) - yb(j)
	    zz = zb(i) - zb(j)
	    rr2 = xx**2 + yy**2 + zz**2 
	    pot = pot - mass(i)*mass(j)/(sqrt(rr2))

         enddo
      enddo

      energy = ke + pot + eoff

      return
      end     !anal_energy_hjs
c--------------------------------------------------------------------------

