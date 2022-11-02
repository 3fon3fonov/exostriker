c*************************************************************************
c                            STEP_KDK_TP_HJS.F
c*************************************************************************
c This subroutine takes a step in Generalized Jacobi coords (HJS)  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 matp,umatp     ==>  Conversion vectors for tp's
c                                      (2D real arrays)
c                 oloct          ==>  Link between tp's and orbits
c                                      (2D integer array)
c                 eta,mu,etatp   ==>  Masses of centers & sats (real arrays)
c                 xjbeg,yjbeg,zjbeg ==> Bodies Jac. position at beginning
c                                       (real arrays)
c                 xbbeg,ybbeg,zbbeg ==> Bodies bary position at beginning
c                                       (real arrays)
c                 vxjbeg,vyjbeg,vzjbeg ==> Bodies Jac. veloc. at beginning
c                                       (real arrays)
c                 ir3jbeg,ir3jend ==> 1/rj^3, beginnig and end (real arrays)
c                 axbbeg,aybbeg,azbbeg ==> Bodies bary accs. at beginning
c                                       (real arrays)
c                 vxjh,vyjh,vzjh    ==> Bod. velocities at middle point
c                                       (real arrays)
c                 xjend,yjend,zjend ==> Bodies Jac. position at end
c                                       (real arrays)
c                 xbend,ybend,zbend ==> Bodies bary position at end
c                                       (real arrays)
c                 axbend,aybend,azbend ==> Bodies bary accs. at end
c                                       (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adapted from step_kdk_tp.f
c Authors:  Herve Beust
c Date:    Feb 12, 2002 
c Last revision: Feb 06, 2003

      subroutine step_kdk_tp_hjs(i1st,nbod,ntp,matp,umatp,oloct,mass,
     &              eta,mu,etatp,xjbeg,yjbeg,zjbeg,vxjbeg,vyjbeg,vzjbeg,
     &              ir3jbeg,xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,
     &              vxjh,vyjh,vzjh,xjend,yjend,zjend,ir3jend,
     &              xbend,ybend,zbend,axbend,aybend,azbend,
     &              xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st,oloct(NPLMAX,NTPMAX)
      real*8 umatp(NPLMAX,NTPMAX),matp(NPLMAX,NTPMAX)
      real*8 mass(nbod),eta(nbod),mu(nbod),dt
      real*8 xjbeg(nbod),yjbeg(nbod),zjbeg(nbod)
      real*8 vxjbeg(nbod),vyjbeg(nbod),vzjbeg(nbod)
      real*8 xbbeg(nbod),ybbeg(nbod),zbbeg(nbod)
      real*8 axbbeg(nbod),aybbeg(nbod),azbbeg(nbod)
      real*8 ir3jbeg(nbod),ir3jend(nbod),etatp(ntp)
      real*8 xjend(nbod),yjend(nbod),zjend(nbod)
      real*8 xbend(nbod),ybend(nbod),zbend(nbod)
      real*8 axbend(nbod),aybend(nbod),azbend(nbod)
      real*8 vxjh(nbod),vyjh(nbod),vzjh(nbod)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)

c...  Internals:
      integer j,k,istati(NSTAT)
      real*8 dth,axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp 
      real*8 axjt(NTPMAX),ayjt(NTPMAX),azjt(NTPMAX)
      real*8 xbttp,zbttp,ybttp,vxbttp,vybttp,vzbttp 

      save axjt,ayjt,azjt     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  loop over all tp's
      do j = 1,ntp
        if (istat(j,1).eq.0) then
          if (i1st.eq.0) then 
c...      Convert to barycentric coordinates
            call coord_g2b_tp(nbod,umatp(1,j),xjbeg,yjbeg,zjbeg,
     &          vxjbeg,vyjbeg,vzjbeg,xjt(j),yjt(j),zjt(j),vxjt(j),
     &          vyjt(j),vzjt(j),xbttp,ybttp,zbttp,vxbttp,vybttp,vzbttp)

c...      Get the acceleration in bary frame.
            call getacch_tp_hjs(nbod,mass,eta,mu,xjbeg,yjbeg,
     &           zjbeg,xbbeg,ybbeg,zbbeg,ir3jbeg,oloct(1,j),
     &           etatp(j),xbttp,ybttp,zbttp,xjt(j),yjt(j),zjt(j),
     &           axbttp,aybttp,azbttp)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
            call coord_vb2vg_tp(nbod,matp(1,j),axbbeg,aybbeg,azbbeg,
     &                 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp)
          else
            axjttp = axjt(j)
            ayjttp = ayjt(j)
            azjttp = azjt(j)
          end if
c...  Apply a Jacobi kick for a half dt 
          call kickvh_tp(1,vxjt(j),vyjt(j),vzjt(j),
     &                            axjttp,ayjttp,azjttp,istat(j,1),dth) 

c...  Take a drift forward full step
          do k=1,NSTAT
            istati(k) = istat(j,k)
          end do
          call drift_tp(1,etatp(j),xjt(j),yjt(j),zjt(j),vxjt(j),
     &                    vyjt(j),vzjt(j),dt,istati)	

c...  After drift, compute bary pos. and vels. 
          call coord_g2b_tp(nbod,umatp(1,j),xjend,yjend,zjend,
     &          vxjh,vyjh,vzjh,xjt(j),yjt(j),zjt(j),vxjt(j),vyjt(j),
     &          vzjt(j),xbttp,ybttp,zbttp,vxbttp,vybttp,vzbttp)

c...  Get the acceleration in bary frame.
          call getacch_tp_hjs(nbod,mass,eta,mu,xjend,yjend,
     &           zjend,xbend,ybend,zbend,ir3jend,oloct(1,j),
     &           etatp(j),xbttp,ybttp,zbttp,xjt(j),yjt(j),zjt(j),
     &           axbttp,aybttp,azbttp)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
          call coord_vb2vg_tp(nbod,matp(1,j),axbend,aybend,azbend,
     &                 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp)
c...  Apply a Jacobi kick for a half dt 
          call kickvh_tp(1,vxjt(j),vyjt(j),vzjt(j),
     &                           axjttp,ayjttp,azjttp,istat(j,1),dth) 
          axjt(j) = axjttp
          ayjt(j) = ayjttp
          azjt(j) = azjttp
        end if
      end do
      i1st = 1

      return
      end   ! step_kdk_tp_hjs
c---------------------------------------------------------------------

