c*************************************************************************
c                        GETACCJ_HJS.F
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
c                 axj,ayj,azj ==>  acceleration in Jac. coord (real arrays)
c
c Author:  H. Beust  
c Date:    Sep. 12, 2006
c Remarks : Adapted from getacch_hjs.f

      subroutine getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xb,yb,zb,ir3j,axj,ayj,azj)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      integer oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod),eta(nbod),mu(nbod)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 xb(nbod),yb(nbod),zb(nbod)

c...  Outputs:
      real*8 axj(nbod),ayj(nbod),azj(nbod),ir3j(nbod)
                
c...  Internals:
      integer i,j,k
      real*8 r2(NPLMAX),r3(NPLMAX),imas(NPLMAX)
      real*8 dx(NPLMAX,NPLMAX),dy(NPLMAX,NPLMAX),dz(NPLMAX,NPLMAX)
      real*8 accijx(NPLMAX,NPLMAX),accijy(NPLMAX,NPLMAX)
      real*8 accijz(NPLMAX,NPLMAX),rij2,irij3
      real*8 epsx,epsy,epsz,eps2,xx,yy,zz,ps,mm

c----
c...  Executable code 

c...  get thr r^-3's

      do i = 2,nbod
         r2(i) = xj(i)*xj(i) + yj(i)*yj(i) + zj(i)*zj(i)
         r3(i) = r2(i)*sqrt(r2(i))
         imas(i) = 1.d0/eta(i)+1.d0/mu(i)
         ir3j(i) = 1.d0/r3(i)
         axj(i) = 0.0d0
         ayj(i) = 0.0d0
         azj(i) = 0.0d0
      end do
      
c...  First compute the cross terms

      do i=1,nbod-1
        do j=i+1,nbod
          dx(i,j) = xb(j) - xb(i)
          dy(i,j) = yb(j) - yb(i)
          dz(i,j) = zb(j) - zb(i)
          rij2 = dx(i,j)*dx(i,j) + dy(i,j)*dy(i,j) + dz(i,j)*dz(i,j)
          irij3 = 1.0d0/(rij2*sqrt(rij2))
          accijx(i,j) = mass(i)*mass(j)*dx(i,j)*irij3
          accijy(i,j) = mass(i)*mass(j)*dy(i,j)*irij3
          accijz(i,j) = mass(i)*mass(j)*dz(i,j)*irij3

          dx(j,i) = -dx(i,j)
          dy(j,i) = -dy(i,j)
          dz(j,i) = -dz(i,j)
          accijx(j,i) = -accijx(i,j)
          accijy(j,i) = -accijy(i,j)
          accijz(j,i) = -accijz(i,j)
        end do
      end do

c... Compute accel for each orbit k
      do k=2,nbod

c...  First compute the contribution of the outer bodies
        do j=1,nbod
          if (oloc(k,j).eq.0) then  ! Body #j is an outer body
            do i=1,nbod
              if (oloc(k,i).eq.1) then !  body #i is a satellite in orbit #k
                axj(k) = axj(k) + accijx(i,j)/mu(k)
                ayj(k) = ayj(k) + accijy(i,j)/mu(k)
                azj(k) = azj(k) + accijz(i,j)/mu(k)
              else if (oloc(k,i).eq.-1) then ! body #i is a center in orbit #k
                axj(k) = axj(k) - accijx(i,j)/eta(k)
                ayj(k) = ayj(k) - accijy(i,j)/eta(k)
                azj(k) = azj(k) - accijz(i,j)/eta(k)
              end if
            end do
          end if
        end do

c... Now compute the internal contributions
        do i=1,nbod
          if (oloc(k,i).eq.1) then  
            do j=1,nbod
              if (oloc(k,j).eq.-1) then ! i = satellite, j = center
                epsx = dx(i,j) - xj(k)
                epsy = dy(i,j) - yj(k)
                epsz = dz(i,j) - zj(k)
                eps2 = epsx*epsx + epsy*epsy + epsz*epsz
                ps = epsx*xj(k) + epsy*yj(k) + epsz*zj(k)
                xx = eps2 + 2.d0*ps
                yy = (r2(k)+xx)*sqrt(r2(k)+xx)
                zz = -xx*ir3j(k)*(3.d0*r2(k)*r2(k)+3.d0*r2(k)*xx+xx*xx)
     &                        /(r3(k)+yy)/yy
                mm = imas(k)*mass(i)*mass(j)
                axj(k) = axj(k) + mm*dx(i,j)*zz
                ayj(k) = ayj(k) + mm*dy(i,j)*zz
                azj(k) = azj(k) + mm*dz(i,j)*zz
              end if
            end do
          end if
        end do
      end do

      return
      end      ! getaccj_hjs

c---------------------------------------------------------------------




