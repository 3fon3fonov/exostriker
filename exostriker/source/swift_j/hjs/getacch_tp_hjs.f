c*************************************************************************
c                        GETACCH_TP_HJS.F
c*************************************************************************
c This subroutine calculates the acceleration on ONE TP
c in the HJS case 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 eta         ==> Masses of centers for orbits (real array)
c                 mu          ==> Masses of satellites for orbits (real arr.)
c                 xj,yj,zj    ==> Pl. positions in jacobi coord (real arrays)
c                 xb,yb,zb    ==> Pl. positions in bary coord (real arrays)
c                 ir3j        ==> Inverse Jacobi radii^3 (real arrays)
c                 oloct       ==> Link between tp and orbits (int array)
c                 etatp       ==> Masses of center for tp (real array)
c                 xjt,yjt,zjt ==>  tp position in jacobi coord (real arrays)
c                 xbt,ybt,zbt ==>  tp position in bary coord (real arrays)
c             Output:
c                 axbt,aybt,azbt ==> tp accel. in bary. coord (real arrays)
c
c Author:  H. Beust  
c Date:    Feb. 05, 2002
c Remarks : Adapted from getacch_hjs.f


      subroutine getacch_tp_hjs(nbod,mass,eta,mu,xj,yj,zj,
     &     xb,yb,zb,ir3j,oloct,etatp,xbt,ybt,zbt,xjt,yjt,zjt,
     &     axbt,aybt,azbt)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      integer oloct(nbod)
      real*8 mass(nbod),eta(nbod),mu(nbod)
      real*8 xj(nbod),yj(nbod),zj(nbod),ir3j(nbod)
      real*8 xb(nbod),yb(nbod),zb(nbod)
      real*8 xjt,yjt,zjt,xbt,ybt,zbt,etatp

c...  Outputs:
      real*8 axbt,aybt,azbt
                
c...  Internals:
      integer i,j
      real*8 fac,dx,dy,dz
      real*8 rji2,ir3jt,irij3,irjt

c----
c...  Executable code 

c...  get thr r^-3
      call getacch_ir3(1,1,xjt,yjt,zjt,ir3jt,irjt)

c... The first term relative to the orbit of the tp
      fac = etatp*ir3jt
      axbt = fac*xjt
      aybt = fac*yjt
      azbt = fac*zjt

c...  now the jacobi terms
      do j=2,nbod    !  Check all the orbits
        if (oloct(j).eq.-1) then !  tp is a center in orbit #i
          fac = -mu(j)*ir3j(j)
          axbt = axbt + fac*xj(j)
          aybt = aybt + fac*yj(j)
          azbt = azbt + fac*zj(j)
        else if (oloct(j).eq.1) then ! tp is a satellite in orbit #i
          fac = eta(j)*ir3j(j)
          axbt = axbt + fac*xj(j)
          aybt = aybt + fac*yj(j)
          azbt = azbt + fac*zj(j)
        end if
      end do

c...  now the third terms. We need to consider the bodies

      do i=1,nbod
        dx = xbt - xb(i)
        dy = ybt - yb(i)
        dz = zbt - zb(i)
        rji2 = dx*dx + dy*dy + dz*dz
        irij3 = 1.0d0/(rji2*sqrt(rji2))
        fac = mass(i)*irij3

        axbt = axbt - fac*dx
        aybt = aybt - fac*dy
        azbt = azbt - fac*dz
      end do

      return
      end      ! getacch_tp_hjs

c---------------------------------------------------------------------

