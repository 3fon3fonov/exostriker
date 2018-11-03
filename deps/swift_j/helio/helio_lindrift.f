c*************************************************************************
c                            HELIO_LINDRIFT.F
c*************************************************************************
c This subroutine takes a linear drift due to mometum of Sun
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 vxb,vyb,vzb   ==>  velocity in bary coord 
c                                    (real arrays)
c                 dt            ==>  time step
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                       (real arrays)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 ptx,pty,ptz  ==> momentum of sun: tp's need this   
c                                       (real scalars)
c
c Remarks: Bases on Martin's code h2.f
c Authors:  Hal Levison 
c Date:    11/14/96
c Last revision: 1/8/97

      subroutine helio_lindrift(nbod,mass,vxb,vyb,vzb,dt,
     &     xh,yh,zh,ptx,pty,ptz)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Outputs Only: 
      real*8 ptx,pty,ptz

c...  Internals:
      integer n

c----
c...  Executable code 

      ptx = mass(2)*vxb(2)
      pty = mass(2)*vyb(2)
      ptz = mass(2)*vzb(2)
      do n=3,nbod
         ptx = ptx + mass(n)*vxb(n)
         pty = pty + mass(n)*vyb(n)
         ptz = ptz + mass(n)*vzb(n)
      enddo

      ptx = ptx/mass(1)
      pty = pty/mass(1)
      ptz = ptz/mass(1)

      do n=2,nbod
         if(mass(n).ne.0.0d0) then
            xh(n) = xh(n) + ptx*dt
            yh(n) = yh(n) + pty*dt
            zh(n) = zh(n) + ptz*dt
         endif
      enddo

      return
      end       ! helio_lindrift
c---------------------------------------------------
