c*************************************************************************
c                        GETACCH_HJS.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in a hierarchical Genealized Jacobi case 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 oloc        ==> Link between bodies and orbits
c                                  (2D integer array)
c                 eta         ==> Masses of centers for orbits (real array)
c                 mu          ==> Masses of satellites for orbits (real arr.)
c                 mass        ==>  mass of bodies (real array)
c                 xj,yj,zj    ==>  position in jacobi coord (real arrays)
c                 xb,yb,zb    ==>  position in bary centric coord (real arrays)
c             Output:
c                 ir3j         ==> Inverse Jacobi radii^3 (real arrays)
c                 axb,ayb,azb ==>  acceleration in bary. coord (real arrays)
c
c Author:  H. Beust  
c Date:    Jan. 24, 2002
c Remarks : Adapted from getacch.f

      subroutine getacch_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xb,yb,zb,ir3j,axb,ayb,azb)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      integer oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod),eta(nbod),mu(nbod)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 xb(nbod),yb(nbod),zb(nbod)

c...  Outputs:
      real*8 axb(nbod),ayb(nbod),azb(nbod),ir3j(nbod)
                
c...  Internals:
      integer i,j
      real*8 irj(NPLMAX),faccen,facsat,faci,facj,dx,dy,dz
      real*8 rji2,irij3

c----
c...  Executable code 

c...  get thr r^-3's
      call getacch_ir3(nbod,2,xj,yj,zj,ir3j,irj)
      do i=1,nbod
        axb(i) = 0.0d0
        ayb(i) = 0.0d0
        azb(i) = 0.0d0
      end do

c...  now the jacobi terms
      do j=2,nbod
        facsat = eta(j)*ir3j(j)
        faccen = -mu(j)*ir3j(j)
        do i=1,nbod
          if (oloc(j,i).eq.1) then   !  body #i is a satellite in orbit #j
            axb(i) = axb(i) + facsat*xj(j)
            ayb(i) = ayb(i) + facsat*yj(j)
            azb(i) = azb(i) + facsat*zj(j)
          else if (oloc(j,i).eq.-1) then ! body #i is a center in orbit #j
            axb(i) = axb(i) + faccen*xj(j)
            ayb(i) = ayb(i) + faccen*yj(j)
            azb(i) = azb(i) + faccen*zj(j)
          end if
        end do
      end do

c      do i=1,nbod
c        axb(i) = 0.0d0
c        ayb(i) = 0.0d0
c        azb(i) = 0.0d0
c      end do

c...  now the third terms

      do i=1,nbod-1
        do j=i+1,nbod
          dx = xb(j) - xb(i)
          dy = yb(j) - yb(i)
          dz = zb(j) - zb(i)
          rji2 = dx*dx + dy*dy + dz*dz
          irij3 = 1.0d0/(rji2*sqrt(rji2))
          faci = mass(i)*irij3
          facj = mass(j)*irij3

          axb(j) = axb(j) - faci*dx
          ayb(j) = ayb(j) - faci*dy
          azb(j) = azb(j) - faci*dz

          axb(i) = axb(i) + facj*dx
          ayb(i) = ayb(i) + facj*dy
          azb(i) = azb(i) + facj*dz
        end do
      end do

      return
      end      ! getacch_hjs

c---------------------------------------------------------------------




