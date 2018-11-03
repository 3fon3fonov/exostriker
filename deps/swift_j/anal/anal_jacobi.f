c*************************************************************************
c                          ANAL_JACOBI.F
c*************************************************************************
c Calculates the Jacobi integral of a particle given the position and velocity
c and ang. velocity vector of the Sun and Planet
c
c      Input:
c            msun          ==>   Mass of the Sun (real scalar)
c            mpl           ==>   Mass of the Planet (real scalar)
c            omegax,..y,z  ==>   vector components of planet-sun orbit's
c                                 angular velocity vector(real scalars)
c            xb,yb,zb      ==>   Barycentric position of particle 
c                                   (real scalars)
c            vxb,vyb,vzb   ==>   Barycentric vel of particle (real scalars)
c            xsb,ysb,zsb   ==>   Barycentric position of Sun (real scalars)
c            xplb,yplb,zplb ==>  Barycentric position of planet (real scalars)
c
c       Output:
c            jacobi         ==>  Value of the jacobi constant
c
c Remarks: 
c Authors:  Martin Duncan
c Date:   April 13/93 
c Last revision:  3/4/93 HFL 

      subroutine anal_jacobi(msun,mpl,omegax,omegay,omegaz,xb,yb,zb,
     &      vxb,vyb,vzb,xsb,ysb,zsb,xplb,yplb,zplb,jacobi)

      include '../swift.inc'

c...  Inputs: 
      real*8 msun,mpl,omegax,omegay,omegaz,xb,yb,zb,vxb,vyb,vzb
      real*8 xsb,ysb,zsb,xplb,yplb,zplb

c...  Outputs:
      real*8 jacobi

c...  Internals
      real*8 rr,jx,jy,jz

c----
c...  Executable code 

      jacobi = 0.5d0*(vxb**2 + vyb**2 + vzb**2)

      rr = sqrt((xb-xsb)**2 + (yb-ysb)**2 + (zb-zsb)**2)
      jacobi = jacobi - msun/rr

      rr = sqrt((xb-xplb)**2 + (yb-yplb)**2 + (zb-zplb)**2)
      jacobi = jacobi - mpl/rr

      jx = yb*vzb - zb*vyb
      jy = zb*vxb - xb*vzb
      jz = xb*vyb - yb*vxb

      jacobi = jacobi - omegax*jx - omegay*jy - omegaz*jz

      return
      end       ! anal_jacobi
c-------------------------------------------------------------------------





