c***************************************************************************
c			OBL_POT.F
c*************************************************************************
c OBL_POT returns the total potential in the barycentric frame for NBOD
c particles due to the oblateness of mass(1) using  
c the values of J2RP2 and J4RP4 passed into the routine.
c (J2RP2 for example is the product of 
c J_2 times the square of the central body's radius)
c Here we return the potential produced
c only by the J2 and J4 terms (i.e. including
c neither the monopole nor higher order terms).
c	
c
c             Input:
c                 nbod     ==>  number of massive bodies (incl. central one)
c                 mass(*)  ==>  masses of particles (real*8 array)
c                 j2rp2    ==>  scaled value of j2 moment (real*8 scalar)
c                 j4rp4    ==>  scaled value of j4 moment (real*8 scalar)
c                 xh(*),yh(*),zh(*)   ==>  HELIO. positions of particles
c                                    (real*8 vectors)
c                 irh(*)   ==> 1./ magnitude of radius vector (real*8 vector)
c                                (passed in to save calcs.)
c
c             Output:
c                 oblpot  ==>  BARY. potential
c                                        (real*8 scalar) 
c
c Remarks:  
c Authors:  Martin Duncan 
c Date:    3/4/94
c Last revision: 

      subroutine obl_pot(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &                oblpot)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(NPLMAX)
      real*8 j2rp2,j4rp4
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),irh(NPLMAX)

c...  Output
      real*8 oblpot

c...  Internals
      integer n
      real*8 rinv2,t0,t1,t2,t3
      real*8 p2,p4

c----
c...  executable code

c Sum all the the bary terms for each "planet" due to central oblate "sun"

	oblpot = 0.d0
	do n =2,nbod

c Note that here we assume we know inverse of radius rather than calc. it
c from (x,y,z) to save the sqrt.
	  rinv2= irh(n)**2
	  t0 = mass(1)*mass(n)*rinv2*irh(n)
	  t1 = j2rp2
	  t2 = zh(n)*zh(n)*rinv2
	  t3 = j4rp4*rinv2

	  p2 = 0.5d0*(3.d0*t2 - 1.d0)
	  p4 = 0.125d0*( (35.d0*t2 -30.d0)*t2 +3.d0)
      
	  oblpot = oblpot + t0*(t1*p2 + t3*p4)

	enddo

        return	
        end                       !  obl_pot.f
c____________________________________________________________________________

