c*************************************************************************
c                        LYAP2_DRIFT_ONE.F
c*************************************************************************
c This subroutine does the danby-type drift for one particle, using 
c appropriate vbles and redoing a drift if the accuracy is too poor 
c (as flagged by the integer iflg).  Also integrates the difference equations.
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mu            ==>  mass of central body (real scalar) 
c                 x,y,z         ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 vx,vy,vz      ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 dx,dy,dz      ==>  initial difference in position 
c                                    (real scalar)
c              dvx,dvy,dvz      ==>  initial difference in velocity
c                                    (real scalar)
c                 dt            ==>  time step
c             Output:
c                 x,y,z         ==>  final position in jacobi coord 
c                                       (real scalars)
c                 vx,vy,vz      ==>  final position in jacobi coord 
c                                       (real scalars)
c                 dx,dy,dz      ==>  final difference in position 
c                                    (real scalar)
c              dvx,dvy,dvz      ==>  final difference in velocity
c                                    (real scalar)
c                 iflg          ==>  integer (zero for successful step)
c
c Comment: Based on drift_one.f
c Authors:  Hal Levison 
c Date:    7/11/95
c Last revision: 
c

      subroutine lyap2_drift_one(mu,x,y,z,vx,vy,vz,dx,dy,dz,
     &     dvx,dvy,dvz,dt,iflg)

      include '../swift.inc'

c...  Inputs Only: 
      real*8 mu,dt

c...  Inputs and Outputs:
      real*8 x,y,z
      real*8 vx,vy,vz
      real*8 dx,dy,dz
      real*8 dvx,dvy,dvz

c...  Output
	integer iflg
	
c...  Internals:
	integer i
	real*8 dttmp

c----
c...  Executable code 

           call lyap2_drift_dan(mu,x,y,z,vx,vy,vz,dx,dy,dz,
     &     dvx,dvy,dvz,dt,iflg)

	   if(iflg .ne. 0) then
	    
	     do i = 1,10
	       dttmp = dt/10.d0
               call lyap2_drift_dan(mu,x,y,z,vx,vy,vz,dx,dy,dz,
     &              dvx,dvy,dvz,dttmp,iflg)
	       if(iflg .ne. 0) return
	     enddo

	   endif

        return
        end    ! lyap2_drift_one
c-------------------------------------------------------------------
