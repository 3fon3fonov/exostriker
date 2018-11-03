c*************************************************************************
c                        HELIO_DRAG.F
c*************************************************************************
c This subroutine loops through the particles (EXCEPT j=2) and applies
c kick due to gas drag
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord
c                                    (real arrays)
c                 dt            ==>  time step
c                 rpl           ==>  physical size of a planet
c                                    (real array)
c                 kdrag0,nkdrag ==>  Drag constant and exponent
c                                    (real scalars)
c                 eta0,neta     ==>  Gas velocity parameter and exponent
c                                    (real scalars)
c             Output:
c                 vxh,vyh,vzh   ==>  final velocity in helio coord
c                                    (real arrays)
c
c Remarks:
c Authors:  Man Hoi Lee
c Date:    12/22/99
c Last revision:

      subroutine helio_drag(nbod,mass,xh,yh,zh,vxh,vyh,vzh,dt,rpl,
     &     kdrag0,nkdrag,eta0,neta)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 rpl(nbod),kdrag0,nkdrag,eta0,neta

c...  Inputs and Outputs:
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Internals:
      integer j
      real*8 varpi,varpi2,rj,rj2,kdrag,vgas,vgasx,vgasy,u,fac

c----
c...  Executable code 

      do j = 3,nbod
         if(mass(j).ne.0.0d0) then
            varpi2 = xh(j)**2 + yh(j)**2
            varpi = sqrt(varpi2)
            rj2 = varpi2 + zh(j)**2
            rj = sqrt(rj2)
            kdrag = kdrag0*(varpi**nkdrag)*(rpl(j)**2)/mass(j)
            vgas = sqrt((mass(1)/rj)*
     &                  (varpi2/rj2 -2.d0*eta0*(rj**neta)))
            vgasx = (-yh(j)/varpi)*vgas
            vgasy =  (xh(j)/varpi)*vgas
            u = sqrt((vxh(j) - vgasx)**2 + (vyh(j) - vgasy)**2 +
     &               vzh(j)**2)
            fac = u*kdrag*dt
            vxh(j) = (vxh(j) + fac*vgasx)/(1 + fac)
            vyh(j) = (vyh(j) + fac*vgasy)/(1 + fac)
            vzh(j) = vzh(j)/(1 + fac)
         endif
      enddo

      return
      end
c--------------------------------------------------------------------------
