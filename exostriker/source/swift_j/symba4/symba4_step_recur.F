c*************************************************************************
c                            SYMBA4_STEP_RECUR.F
c*************************************************************************
c                      THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c
c             Input:
c                 t             ==>  time (real Scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 ireci         ==>  Input recursion level  (integer scalar)
c                 ilevl         ==>  largest recursion level used 
c                                    (integer array)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 ielev         ==>  The level that this particle should go
c                                             (int*2 array)
c                 lemat         ==>  Matrix of encounters (logical*1 2d array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 rhill         ==>  Hill sphere of planet (real Scalar)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial velocity in bari coord 
c                                    (real arrays)
c                dt0            ==>  Global timestep  (real scalar)
c                lclose         ==> .true. --> marge particles if they
c                                    get too close. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                rpl            ==>  physical size of a planet.
c                                    (real array)
c                eoff           ==>  Energy offset (real scalar)
c                ielc           ==>  number of encounters (integer*2 scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxb,vyb,vzb   ==>  final velocity in bari coord 
c                                       (real arrays)
c             mergelst          ==>  list of mergers (int array)
c             mergecnt          ==>  count of mergers (int array)
c                 rpl           ==>  Recalculated physical size of a planet.
c                                    if merger happened (real array)
c                 mass          ==>  Recalculated mass of bodies 
c                                    if merger happened (real array)
c                eoff           ==>  Energy offset (real scalar)
c                svdotr         ==> vdotr relative flag
c                                   = .true. if i,j are receding
c                                   = .false is approaching
c                                     (2D logical*1 array)
c                                   Used internally, but only need 1 copy.
c
c Remarks: If marger occurs, does not change nbod and puts the mass
c          of one of the particles to zero.
c Authors:  Hal Levison 
c Date:    3/20/97
c Last revision: 1/30/98

#ifdef  _RECUR_SUB
      recursive subroutine
#else
      subroutine
#endif
     &   symba4_step_recur(t,nbod,nbodm,mass,ireci,ilevl,lemat,
     &     iecnt,ielev,rhill,xh,yh,zh,vxb,vyb,vzb,lclose,rpl,mergelst,
     &     mergecnt,dt0,eoff,svdotr,ielc,ielst)

      include '../swift.inc'
      include 'symba4.inc'

c...  Inputs Only: 
      integer nbod,ireci,nbodm
      real*8 mass(nbod),dt0,rhill(nbod),t
      integer*2 iecnt(NTPMAX),ielev(nbod)
      logical*1 lemat(NTPMAX,NTPMAX)
      logical*2 lclose 
      integer*2 ielst(2,NTPMAXSQ),ielc

c...  Inputs and Outputs:
      integer ilevl(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod),eoff
      real*8 vxb(nbod),vyb(nbod),vzb(nbod),rpl(nbod)
      integer mergelst(2,NTPMAX),mergecnt
      logical*1 svdotr(NTPMAX,NTPMAX)

c...  Internals:
      integer i,j,ie
      integer icflg,it,irecp,ieflg
      real*8 dtl,dth,sgn

c----
c...  Executable code 

      dtl = dt0/(float(NTENC)**(ireci))
      dth = dtl/2.0d0

      if( (dtl/dt0) .le. TINY ) then
         write(*,*) ' Warning in SYMBA_STEP_RECUR: '
         write(*,*) '         Local timestep too small '
         write(*,*) '         Roundoff will be important!!!! '
      endif

      if(ireci.eq.0) then

         irecp = ireci + 1

c...     Do we need to go deeper?
         icflg = 0
         do ie=1,ielc
            i = ielst(1,ie)
            j = ielst(2,ie)
            if((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
               call symba4_chk(rhill,nbod,i,j,mass,xh,yh,zh,
     &              vxb,vyb,vzb,dtl,irecp,ieflg,svdotr(i,j))
               if(ieflg.ne.0) then
                  icflg = 1
                  ielev(i) = irecp
                  ielev(j) = irecp
                  ilevl(i) = max(irecp,ilevl(i))
                  ilevl(j) = max(irecp,ilevl(j))
               endif
            endif
         enddo

         sgn = 1.0d0
         call symba4_kick(nbod,mass,irecp,iecnt,ielev,
     &        lemat,rhill,xh,yh,zh,vxb,vyb,vzb,dth,sgn,ielc,
     &        ielst)

         call symba4_helio_drift(nbod,ielev,ireci,mass,xh,
     &        yh,zh,vxb,vyb,vzb,dtl)
         if(icflg.ne.0) then
            call symba4_step_recur(t,nbod,nbodm,mass,irecp,ilevl,
     &           lemat,iecnt,ielev,rhill,xh,yh,zh,vxb,vyb,vzb,lclose,
     &           rpl,mergelst,mergecnt,dt0,eoff,svdotr,ielc,ielst)
         endif

         sgn = 1.0d0
         call symba4_kick(nbod,mass,irecp,iecnt,ielev,lemat,
     &        rhill,xh,yh,zh,vxb,vyb,vzb,dth,sgn,ielc,ielst)

         if( lclose ) then       ! look for mergers
            do ie=1,ielc
               i = ielst(1,ie)
               j = ielst(2,ie)
               if((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
                  call symba4_merge(t,dtl,nbod,i,j,mass,xh,yh,zh,
     &                 vxb,vyb,vzb,ireci,ilevl,svdotr(i,j),lemat,
     &                 iecnt,rpl,mergelst,mergecnt,rhill,eoff,ielc,
     &                 ielst)
               endif
            enddo
         endif

         do i=2,nbod
            if(ielev(i).eq.irecp) then
               ielev(i) = ireci
            endif
         enddo


      else

         irecp = ireci + 1
         do it=1,NTENC

c...        Do we need to go deeper?
            icflg = 0

            do ie=1,ielc
               i = ielst(1,ie)
               j = ielst(2,ie)
               if((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
                  call symba4_chk(rhill,nbod,i,j,mass,xh,yh,
     &                 zh,vxb,vyb,vzb,dtl,irecp,ieflg,
     &                 svdotr(i,j))
                  if(ieflg.ne.0) then
                     icflg = 1
                     ielev(i) = irecp
                     ielev(j) = irecp
                     ilevl(i) = max(irecp,ilevl(i))
                     ilevl(j) = max(irecp,ilevl(j))
                  endif
               endif
            enddo

            sgn = 1.0d0
            call symba4_kick(nbod,mass,irecp,iecnt,ielev,lemat,
     &           rhill,xh,yh,zh,vxb,vyb,vzb,dth,sgn,ielc,ielst)
            sgn = -1.0d0
            call symba4_kick(nbod,mass,irecp,iecnt,ielev,lemat,
     &           rhill,xh,yh,zh,vxb,vyb,vzb,dth,sgn,ielc,ielst)


            call symba4_helio_drift(nbod,ielev,ireci,mass,xh,
     &           yh,zh,vxb,vyb,vzb,dtl)
            if(icflg.ne.0) then
               call symba4_step_recur(t,nbod,nbodm,mass,irecp,ilevl,
     &             lemat,iecnt,ielev,rhill,xh,yh,zh,vxb,vyb,vzb,lclose,
     &             rpl,mergelst,mergecnt,dt0,eoff,svdotr,ielc,ielst)
            endif

            sgn = 1.0d0
            call symba4_kick(nbod,mass,irecp,iecnt,ielev,lemat,
     &           rhill,xh,yh,zh,vxb,vyb,vzb,dth,sgn,ielc,ielst)
            sgn = -1.0d0
            call symba4_kick(nbod,mass,irecp,iecnt,ielev,lemat,
     &           rhill,xh,yh,zh,vxb,vyb,vzb,dth,sgn,ielc,ielst)

            if( lclose ) then ! look for mergers
               do ie=1,ielc
                  i = ielst(1,ie)
                  j = ielst(2,ie)
                  if((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
                     call symba4_merge(t,dtl,nbod,i,j,mass,xh,yh,zh,
     &                    vxb,vyb,vzb,ireci,ilevl,svdotr(i,j),lemat,
     &                    iecnt,rpl,mergelst,mergecnt,rhill,eoff,ielc,
     &                    ielst)
                  endif
               enddo
            endif
            
            do i=2,nbod
               if(ielev(i).eq.irecp) then
                  ielev(i) = ireci
               endif
            enddo

         enddo
         
      endif

      return
      end     ! symba4_step_recur.f
c--------------------------------------------------------------


