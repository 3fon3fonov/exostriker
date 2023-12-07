c************************************************************************
c                              IO_INIT_PL_HJS.F
c************************************************************************
c IO_INIT_PL_HJS reads in the data for the massive bodies for HJS 
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 iflgchk        
c
c             Output:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 oloc          ==>  Link matrix between bodies & orbits
c                         oloc(j,i)=1  : body #i is a satellite in orbit #j
c                         oloc(j,i)=-1 : body #i is a center in orbit #j
c                                    (2D integer array)
c                 mass          ==>  mass of bodies (real array)
c                 eta,mu        ==>  masses of centers & satellites 
c                                    (real arrays)
c                 mat           ==>  Conversion matrix Bary => Gen. Jacobi
c                 umat          ==>  Conversion matrix Gen. Jacobi => Bary
c                                     (2D real arrays)
c                 xj,yj,zj      ==>  initial position in Gen. Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial position in Gen. Jacobi coord 
c                                    (real arrays)
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c
c Remarks: Adapted from io_init_pl.f
c Authors:  Herve Beust
c Date:    Feb. 11, 2002
c Last revision: Dec. 9, 2002

	subroutine io_init_pl_hjs(infile,lclose,iflgchk,nbod,oloc,
     &     mass,eta,mu,mat,umat,xj,yj,zj,vxj,vyj,vzj,rplsq)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile
	integer iflgchk
        logical*2 lclose

c...    Output
	real*8 mass(NPLMAX),rplsq(NPLMAX)
	real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)
        real*8 eta(NPLMAX),mu(NPLMAX)
        real*8 mat(NPLMAX,NPLMAX),umat(NPLMAX,NPLMAX)
	integer nbod,oloc(NPLMAX,NPLMAX)

c...   Internal
	integer j,i,ierr
        real*8 rpl,mtot,vsat,vcen
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c-----
c...  Executable code      

	write(*,*) 'Planet data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)

c Read number of planets
	read(7,*) nbod

        if(nbod.gt.NPLMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_pl: '
           write(*,*) '   The number of massive bodies,',nbod,','
           write(*,*) '   is too large, it must be less than',NPLMAX
           call util_exit(1)
        endif

	write(*,23) nbod
23	format(/,'Number of bodies is ',i3,/,
     &   'For each, list mass ',/,
     &   'Followed by x,y,z,vx,vy,vz : '/)

c... Read information relative to bodies: masses, bary pos & vels.
	do j=1,nbod
           if(lclose) then
              read(7,*) mass(j),rpl
              write(*,*) mass(j),rpl
              rplsq(j) = rpl*rpl
           else
              read(7,*) mass(j)
              write(*,*) mass(j)
           endif
	  read(7,*) xb(j),yb(j),zb(j)
	  read(7,*) vxb(j),vyb(j),vzb(j)
	  write(*,*) xb(j),yb(j),zb(j)
	  write(*,*) vxb(j),vyb(j),vzb(j)
	enddo

c... Read information relative to orbits: link matrix oloc
        do i = 1,NPLMAX
          do j = 1,NPLMAX
            oloc(j,i) = 0
          end do
        end do
        do j = 2,nbod
          read(7,*)(oloc(j,i),i=1,nbod)          
        end do

c... Compute eta's and mu's:
        do j = 1,NPLMAX
          eta(j) = 0
          mu(j) = 0
        end do
        do j = 2,nbod
          do i = 1,nbod
            if (oloc(j,i).eq.1) mu(j) = mu(j)+mass(i)
            if (oloc(j,i).eq.-1) eta(j) = eta(j)+mass(i)
          end do
        end do

c... Build transform matrix Barycentric --> Generalized Jacobi
        do j = 1,NPLMAX
          do i = 1,NPLMAX
            mat(i,j) = 0.0d0
          end do
        end do
        mtot = 0.0d0
        do i = 1,nbod
          mtot = mtot+mass(i)
        end do
        do i = 1,nbod
          mat(1,i) = mass(i)/mtot
        end do
        do j = 2,nbod
          do i = 1,nbod
            if (oloc(j,i).eq.1) mat(j,i) = mass(i)/mu(j)
            if (oloc(j,i).eq.-1) mat(j,i) = -mass(i)/eta(j)
          end do
        end do    

c...    Build inverse transform matrix Generalized Jacobi --> Barycentric
        do j = 1,NPLMAX
          do i = 1,NPLMAX
            umat(i,j) = 0.0d0
          end do
        end do
        do i = 1,nbod
          umat(i,i) = 1.0d0
        end do
        do j = 2,nbod
          vsat = eta(j)/(mu(j)+eta(j))
          vcen = -mu(j)/(mu(j)+eta(j))
          do i = 1,nbod
            if (oloc(j,i).eq.1) umat(i,j) = vsat
            if (oloc(j,i).eq.-1) umat(i,j) = vcen
          end do
        end do         

c... Compute Generalized jacobi coordinates
        call coord_b2g(nbod,mat,mass,xb,yb,zb,vxb,vyb,vzb,
     &                               xj,yj,zj,vxj,vyj,vzj) 
	close(unit = 7)
	return
	end     ! io_init_pl_hjs.f

c--------------------------------------------------------------------------

