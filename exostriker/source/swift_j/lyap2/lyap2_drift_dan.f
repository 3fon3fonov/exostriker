c*************************************************************************
c                        LYAP2_DRIFT_DAN.F
c*************************************************************************
c This subroutine does the Danby and decides which vbles to use
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 x0,y0,z0         ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 vx0,vy0,vz0      ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 dx0,dy0,dz0      ==>  initial difference in position 
c                                    (real scalar)
c              dvx0,dvy0,dvz0      ==>  initial difference in velocity
c                                    (real scalar)
c                 dt0            ==>  time step
c             Output:
c                 x0,y0,z0         ==>  final position in jacobi coord 
c                                       (real scalars)
c                 vx0,vy0,vz0      ==>  final position in jacobi coord 
c                                       (real scalars)
c                 dx0,dy0,dz0      ==>  final difference in position 
c                                    (real scalar)
c              dvx0,dvy0,dvz0      ==>  final difference in velocity
c                                    (real scalar)
c                 iflg             ==>  integer flag (zero if satisfactory)
c					      (non-zero if nonconvergence)
c
c Comments: Based on drift_dan.f
c Authors:  Hal Levison 
c Date:    7/11/95
c Last revision: 

      subroutine lyap2_drift_dan(mu,x0,y0,z0,vx0,vy0,vz0,dx0,
     &     dy0,dz0,dvx0,dvy0,dvz0,dt0,iflg)

      include '../swift.inc'

c...  Inputs Only: 
      real*8 mu,dt0

c...  Inputs and Outputs:
      real*8 x0,y0,z0
      real*8 vx0,vy0,vz0
      real*8 dx0,dy0,dz0
      real*8 dvx0,dvy0,dvz0

c...  Output
      integer iflg

c...  Internals:
      real*8 x,y,z,vx,vy,vz,dt
      real*8 sqrtmu,sn,alfa,eta,zeta,x0dx,w0dx,w0dv,r0pr,alfapr
      real*8 x0dv,etapr,zetapr,g3a,g2a,g1a,xpr,g0,g1pr,g2pr,g3pr
      real*8 rpr,fpr,gpr,dfpr,dgpr
      real*8 wx,wy,wz
      real*8 dwx0,dwy0,dwz0
      real*8 dwx,dwy,dwz
      real*8 dx,dy,dz
      real*8 f,g,fdot,gdot
      real*8 u,alpha,fp,r0,v0s,s
      real*8 cn1,cn2,cn3,cn4,cn5
      real*8 c1,c2,c3,c4,c5

c----
c...  Executable code 

c...  Set dt = dt0 to be sure timestep is not altered while solving
c...  for new coords.
      dt = dt0
      iflg = 0
      r0 = sqrt(x0*x0 + y0*y0 + z0*z0)
      v0s = vx0*vx0 + vy0*vy0 + vz0*vz0
      u = x0*vx0 + y0*vy0 + z0*vz0
      alpha = 2.0*mu/r0 - v0s
        
      call lyap2_drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,c4,
     &     c5,s,iflg)

      if(iflg .eq.0) then
         f = 1.0 - (mu/r0)*c2
         g = dt - mu*c3
         fdot = -(mu/(fp*r0))*c1
         gdot = 1. - (mu/fp)*c2

         x = x0*f + vx0*g
         y = y0*f + vy0*g
         z = z0*f + vz0*g
         vx = x0*fdot + vx0*gdot
         vy = y0*fdot + vy0*gdot
         vz = z0*fdot + vz0*gdot

c--------------------------------------
c...     now do the differentials from Mikkola's code

         sqrtmu = sqrt(mu)
         wx = vx0/sqrtmu
         wy = vy0/sqrtmu
         wz = vz0/sqrtmu
         dwx0 = dvx0/sqrtmu
         dwy0 = dvy0/sqrtmu
         dwz0 = dvz0/sqrtmu

         cn1 = c1*sqrtmu
         cn2 = c2*mu
         cn3 = c3*(sqrtmu*mu)
         cn4 = c4*(mu**2)
         cn5 = c5*(sqrtmu*mu**2)
         sn = s*sqrtmu

         alfa = 2.d0/r0 - (wx**2 + wy**2 + wz**2)
         eta = x0*wx + y0*wy + z0*wz
         zeta = 1.0d0 - alfa*r0

         x0dx = x0*dx0 + y0*dy0 + z0*dz0
         x0dv = x0*dwx0 + y0*dwy0 + z0*dwz0
         w0dx = wx*dx0 + wy*dy0 + wz*dz0
         w0dv = wx*dwx0 + wy*dwy0 + wz*dwz0
         r0pr=x0dx/r0
         alfapr = -2.0d0*(r0pr/(r0**2) + w0dv)
         etapr = w0dx + x0dv
         zetapr = -alfa*r0pr - r0*alfapr
         g3a = 0.5d0 * (3.0d0*cn5 - sn*cn4)
         g2a = 0.5d0 * (2.0d0*cn4 - sn*cn3)
         g1a = 0.5d0 * (      cn3 - sn*cn2)
         xpr = -( sn*r0pr + cn2*etapr + cn3*zetapr +
     &      (eta*g2a+zeta*g3a)*alfapr )/fp       
         g0 = 1.d0 - alfa*cn2
         g1pr = g0*xpr + g1a*alfapr
         g2pr = cn1*xpr + g2a*alfapr
         g3pr = cn2*xpr + g3a*alfapr
         rpr = r0pr + cn1*etapr + cn2*zetapr + eta*g1pr + zeta*g2pr
         fpr = cn2*r0pr/r0**2 - g2pr/r0
         gpr = -g3pr 
         dfpr = -g1pr/(r0*fp) + cn1*rpr/(fp*fp*r0) + 
     &        cn1*r0pr/(r0*r0*fp)
         dgpr = -g2pr/fp + cn2*rpr/(fp*fp)
 
         dx = f*dx0 + g*dwx0*sqrtmu + x0*fpr + wx*gpr
         dy = f*dy0 + g*dwy0*sqrtmu + y0*fpr + wy*gpr
         dz = f*dz0 + g*dwz0*sqrtmu + z0*fpr + wz*gpr

         dwx = fdot*dx0/sqrtmu + gdot*dwx0 + x0*dfpr + wx*dgpr
         dwy = fdot*dy0/sqrtmu + gdot*dwy0 + y0*dfpr + wy*dgpr
         dwz = fdot*dz0/sqrtmu + gdot*dwz0 + z0*dfpr + wz*dgpr

         dx0 = dx
         dy0 = dy
         dz0 = dz
         dvx0 = dwx*sqrtmu
         dvy0 = dwy*sqrtmu
         dvz0 = dwz*sqrtmu

c-----------------------------------------

         x0 = x
         y0 = y
         z0 = z
         vx0 = vx
         vy0 = vy
         vz0 = vz
      endif
        
      return
      end                       ! lyap2_drift_dan


