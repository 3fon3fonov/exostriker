c********************************************************************#
c                  DRIFT_KEPMD
c********************************************************************#
c  Subroutine for solving kepler's equation in difference form for an
c  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
c  for the criteria.
c  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
c  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
c  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.
c
c	Input:
c	    dm		==> increment in mean anomaly M (real*8 scalar)
c	    es,ec       ==> ecc. times sin and cos of E_0 (real*8 scalars)
c
c       Output:
c            x          ==> solution to Kepler's difference eqn (real*8 scalar)
c            s,c        ==> sin and cosine of x (real*8 scalars)
c

        subroutine drift_kepmd(dm,es,ec,x,s,c)

	implicit none

c...    Inputs
	real*8 dm,es,ec
	
c...	Outputs
	real*8 x,s,c

c...    Internals
	real*8 A0, A1, A2, A3, A4
        parameter(A0 = 39916800.d0, A1 = 6652800.d0, A2 = 332640.d0)
	parameter(A3 = 7920.d0, A4 = 110.d0)
	real*8 dx
	real*8 fac1,fac2,q,y
	real*8 f,fp,fpp,fppp


c...    calc initial guess for root
	fac1 = 1.d0/(1.d0 - ec)
	q = fac1*dm
	fac2 = es*es*fac1 - ec/3.d0
	x = q*(1.d0 -0.5d0*fac1*q*(es -q*fac2))

c...  excellent approx. to sin and cos of x for small x.
	y = x*x
	s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
        c = sqrt(1.d0 - s*s)

c...    Compute better value for the root using quartic Newton method
        f = x - ec*s + es*(1.-c) - dm
        fp = 1. - ec*c + es*s
        fpp = ec*s + es*c
        fppp = ec*c - es*s
        dx = -f/fp
        dx = -f/(fp + 0.5*dx*fpp)
        dx = -f/(fp + 0.5*dx*fpp + 0.16666666666666666*dx*dx*fppp)
        x = x + dx
     
c...  excellent approx. to sin and cos of x for small x.
	y = x*x
	s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
        c = sqrt(1.d0 - s*s)

	return
	end
c-----------------------------------------------------------------------------
