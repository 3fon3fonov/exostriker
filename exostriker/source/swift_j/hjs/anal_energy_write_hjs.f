c*************************************************************************
c                          ANAL_ENERGY_WRITE_HJS.F
c*************************************************************************
c Writes the energy of the total system (massive bodies) wrt time in
c   the hierarchical Jacobi case
c
c      Input:
c            t             ==>  current time
c            nbod          ==>  number of massive bodies (int scalar)
c            mass          ==>  mass of bodies (real array)
c            umat          ==>  Conversion matrix Jacobi => Barycentric
c                                (2D real array)
c            xj,yj,zj      ==>  current position in Jacobi coord 
c                               (real arrays)
c            vxj,vyj,vzj   ==>  current velocity in Jacobi coord 
c                               (real arrays)
c            iu            ==>  unit to write to (int scalar)
c            fopenstat     ==>  The status flag for the open 
c                                statements of the output files.  
c                                          (character*80)
c            eoff          ==> An energy offset that is added to the energy
c                                      (real*8 scalar)
c
c Remarks: Adapted from anal_energy_write & anal_energy
c Authors:  Herve Beust
c Date:    Sep 4, 2006

      subroutine anal_energy_write_hjs(t,nbod,umat,mass,xj,yj,zj,
     &     vxj,vyj,vzj,iu,fopenstat,diro,eoff)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,iu
      real*8 mass(nbod),t,umat(NPLMAX,NPLMAX),eoff
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      character*80 fopenstat,diro

c...  Internals
      integer i1st
      real*8 energy,eltot(3),ke,pot
      real*8 elx,ely,elz
      real*8 xx,yy,zz,rr2
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      integer i,j
      data i1st/0/
      save i1st

c----
c...  Executable code 

c Compute and print initial ke,pot,energy and ang. mom.

      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xb,yb,zb,vxb,vyb,vzb) 

      eltot(1)=(yb(nbod)*vzb(nbod)-zb(nbod)*vyb(nbod))*mass(nbod)
      eltot(2)=(zb(nbod)*vxb(nbod)-xb(nbod)*vzb(nbod))*mass(nbod)
      eltot(3)=(xb(nbod)*vyb(nbod)-yb(nbod)*vxb(nbod))*mass(nbod)

      ke = 0.5*mass(nbod)*(vxb(nbod)**2 + vyb(nbod)**2 + vzb(nbod)**2)
      pot = 0.d0

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
	    pot = pot - mass(i)*mass(j)/(sqrt(rr2))

         enddo
      enddo

      energy = ke + pot + eoff

      call io_energy_write(i1st,t,energy,eltot,iu,diro,fopenstat)

      if(i1st.eq.0) then
         i1st=1
      endif

      return
      end     !anal_energy_write_hjs
c--------------------------------------------------------------------------

