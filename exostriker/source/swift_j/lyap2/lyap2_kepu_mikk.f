c*************************************************************************
c                        LYAP2_KEPU_MIKK
c*************************************************************************
c Given an s, returns mikkola's c1,c2,c3,c4,c5
c
c             Input:
c                 alpha         ==>  Twice the binding energy (real scalar)
c                 s             ==>  Approx. root of f 
c             Output:
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 c4,c5         ==>  Mikkola's extra c'c
c                                       (real scalors)
c
c Based on Mikkola's code
c Author:  Hal Levison
c Date:    7/11/95
c Last revision:

      subroutine lyap2_kepu_mikk(s,alpha,c1,c2,c3,c4,c5)

      include '../swift.inc'

c...  Inputs: 
      real*8 alpha,s

c...  Outputs:
      real*8 c1,c2,c3,c4,c5

c...  Internals:
      real*8  x,h
      real*8 cc6,cc132,cc56,cc30,cc24,cc156,cc90,cc110,cc16
      real*8 cc8,cc72,cc42,cc120,cc2,cc3
      integer k,i

      parameter(cc6=1.d0/6.d0,cc132=1.d0/132.d0,cc56=1.d0/56.d0,
     &  cc30=1.d0/30d0,cc24=1.d0/24.d0,cc156=1.d0/156.d0,
     &  cc90=1.d0/90.d0,cc110=1.d0/110.d0,cc16=1.d0/16.d0,
     &  cc8=1.d0/8.d0,cc72=1.d0/72.d0,cc42=1.d0/42.d0,
     &  cc120=1.d0/120.d0)

c----
c...  Executable code 

      x=s*s*alpha

c...  From Mikkola's code

      h=x
      do k=0,10
         if(abs(h) .ge. 0.1) goto 1
         h=0.25d0*h
      enddo
 1    continue

      c4 = (1.0d0-h*(1.0d0-h*(1.0d0-h*cc90/(1.0d0+h*cc132))*
     &     cc56)*cc30)*cc24
      c5 = (1.0d0-h*(1.0d0-h*(1.0d0-h*cc110/(1.0d0+h*cc156))*
     &     cc72)*cc42)*cc120

      do i=1,k
         cc3 =  cc6 - h*c5
         cc2 = 0.5D0 - h*c4
         c5 = (c5 + c4+cc2*cc3)*cc16
         c4 = cc3*(2.d0-h*cc3)*cc8
         h=4.d0*h
      enddo

      c3 = cc6 - x*c5
      c2 = 0.5d0 - x*c4
      c1 = 1.d0 - x*c3

c...  End of Mikkola's code

      c1=c1*s
      c2 = c2*s*s
      c3 = c3*s*s*s
      c4 = c4*s*s*s*s
      c5 = c5*s*s*s*s*s

      return
      end                       !   lyap2_kepu_mikk
c-------------------------------------------------------------------
