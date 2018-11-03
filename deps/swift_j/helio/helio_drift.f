c*************************************************************************
c                        HELIO_DRIFT.F
c*************************************************************************
c This subroutine loops thorugh the particles and calls the danby routine
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial position in bary coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxb,vyb,vzb   ==>  final position in bary coord 
c                                       (real arrays)
c
c Remarks:  Based on drift.f
c Authors:  Hal Levison 
c Date:    11/14/96
c Last revision: 1/8/97  for symba

      subroutine helio_drift(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Internals:
      integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth

      do j = 2,nbod
         if(mass(j).ne.0.0d0) then
            call drift_one(mass(1),xh(j),yh(j),zh(j),
     &           vxb(j),vyb(j),vzb(j),dt,iflg)
            if(iflg.ne.0) then
               write(*,*) ' Planet ',j,' is lost !!!!!!!!!'
               write(*,*) mass(1),dt
               write(*,*) xh(j),yh(j),zh(j),' H '
               write(*,*) vxb(j),vyb(j),vzb(j),' B '
               write(*,*) ' STOPPING '
               call util_exit(1)
            endif
         endif
      enddo

      return
      end
c--------------------------------------------------------------------------
