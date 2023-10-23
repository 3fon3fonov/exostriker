c**********************************************************************
c			IO_INIT_TP_HJS.F
c**********************************************************************
c Read in test particle data for HJS
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c                 oloc          ==>  Link matrix between bodies & orbits
c                         oloc(j,i)=1  : body #i is a satellite in orbit #j
c                         oloc(j,i)=-1 : body #i is a center in orbit #j
c                 nbod          ==>  Number of massive bodies (int scalar)
c                 mass          ==> Masses of massive bodies (real array)
c                 umat          ==> Conversion matrix Jacobi=>Bary
c                                   (2D real array)
c
c             Output:
c                 ntp           ==>  number of massive bodies (int scalar)
c                 oloct         ==>  Position of tp relative to orbits
c                          oloct(j,i)=1 : tp #i is a satellite in orbit #j
c                          oloct(j,i)=-1: tp #i is a center in orbit #j
c                 etatp         ==>  mass of centers for tp's (real array)
c                 matp          ==> Conversion vectors bary=>Jacobi for tp's
c                 umatp         ==> Conversion vectors Jacobi=>Bary for tp's
c                                       (2D real arrays)
c              xjt,yjt,zjt      ==>  initial position in Jacobi coord 
c                                    (real arrays)
c              vxjt,vyjt,vzjt   ==>  initial position in Jacobi coord 
c                                    (real arrays)
c               istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0  active
c                                      istat(i,1) = 1 not
c               rstat           ==>  status of the test paricles
c                                      (2d  real array)
c
c
c
c Remarks: Adapted from io_init_tp.f
c Authors:  Herve Beust
c Date:    Feb. 11, 2002

	subroutine io_init_tp_hjs(infile,nbod,ntp,oloc,oloct,mass,
     &       umat,etatp,matp,umatp,xjt,yjt,zjt,vxjt,vyjt,vzjt,
     &                                                 istat,rstat)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile
        integer iflgchk
        integer nbod,oloc(NPLMAX,NPLMAX)
        real*8 mass(NPLMAX)
        real*8 umat(NPLMAX,NPLMAX)

c...    Output
	real*8 xjt(NTPMAX),yjt(NTPMAX),zjt(NTPMAX)
	real*8 vxjt(NTPMAX),vyjt(NTPMAX),vzjt(NTPMAX)
        real*8 rstat(NTPMAX,NSTATR),etatp(NTPMAX)
        real*8 matp(NPLMAX,NTPMAX),umatp(NPLMAX,NTPMAX)
	integer istat(NTPMAX,NSTAT),oloct(NPLMAX,NTPMAX)

	integer ntp

c...   Internal
	integer i,j,k,ierr,orbct(NPLMAX)
        logical sat,cen

c-----
c...  Executable code      

	write(*,*) 'Test particle file called ',infile
        call io_open(7,infile,'old','formatted',ierr)

	read(7,*) ntp

        if(ntp.gt.NTPMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_tp: '
           write(*,*) '   The number of test bodies,',ntp,','
           write(*,*) '   is too large, it must be less than',NTPMAX
           call util_exit(1)
        endif

	write(*,*) ' '
	write(*,*) 'ntp : ',ntp

c Read in the x's and v's and istat(*,*)
	  write(*,*) ' '
	  do  i=1,ntp
            do j=1,NPLMAX
              oloct(j,i) = 0
            end do
	    read(7,*) xjt(i),yjt(i),zjt(i)
	    read(7,*) vxjt(i),vyjt(i),vzjt(i)
            read(7,*) (orbct(j),j=1,nbod)
	    read(7,*) (istat(i,j),j=1,NSTAT)
	    read(7,*) (rstat(i,j),j=1,NSTATR)
            etatp(i) = 0.0d0   
            do j=1,NPLMAX
              matp(j,i) = 0.0d0
              umatp(j,i) = 0.0d0
            end do
            do j = 1,nbod
              if (orbct(j).eq.-1) etatp(i) = etatp(i)+mass(j)
            end do
            do j = 1,nbod
              if (orbct(j).eq.-1) matp(j,i) = -mass(j)/etatp(i)
            end do
            do k = 1,nbod
              if (orbct(k).eq.-1) then
                do j=1,nbod     
                  umatp(j,i) = umatp(j,i) - matp(k,i)*umat(k,j)
                end do
              end if
            end do

c... Computation of the oloct array on the basis or the centers
c... listed in the orbct array

            do j=1,NPLMAX
              oloct(j,i) = 0
            end do
            do j = 2,nbod
              sat = .true.
              cen = .true.
              do k = 1,nbod
                if (orbct(k).eq.-1) then                  
                  sat = sat.and.(oloc(j,k).eq.1)
                  cen = cen.and.(oloc(j,k).eq.-1)
                end if
              end do
              if (sat) oloct(j,i) = 1
              if (cen) oloct(j,i) = -1
            end do

	  enddo

	close(unit = 7)
        write(*,*) ' '

	return
	end    ! io_init_tp_hjs.f
c-----------------------------------------------------------------

