C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C 
        PROGRAM GEN_PL_HJS

        include '../../swift.inc'

        real*8 SMASSYR,DR,EPOCH,YEAR,AU
        PARAMETER( SMASSYR = TWOPI*TWOPI,
     &          DR = 1.7453292519943294d-2,      ! 180/pi 
     &          EPOCH = 2449101.0d0,
     &          YEAR = 365.2422d0, 
     &          AU = 1.495978707d8 )               ! km

        REAL*8 MASS(NPLMAX),PERIOD,JDPERI,Q,RPLSQ(NPLMAX),
     &         ETA(NPLMAX),MU(NPLMAX),MAT(NPLMAX,NPLMAX),
     &         UMAT(NPLMAX,NPLMAX),MTOT,VSAT,VCEN,
     &         XB(NPLMAX),YB(NPLMAX),ZB(NPLMAX),
     &         VXB(NPLMAX),VYB(NPLMAX),VZB(NPLMAX),
     &         XCM,YCM,ZCM,VXCM,VYCM,VZCM,
     &         XJ(NPLMAX),YJ(NPLMAX),ZJ(NPLMAX),
     &         VXJ(NPLMAX),VYJ(NPLMAX),VZJ(NPLMAX),
     &         XJO(NPLMAX),YJO(NPLMAX),ZJO(NPLMAX),
     &         VXJO(NPLMAX),VYJO(NPLMAX),VZJO(NPLMAX),
     &         GM,A,E,INC,CAPOM,OMEGA,CAPM,MSYS,
     &         FACHILL,R2HILL(NPLMAX)
        LOGICAL*2 LCLOSE
        logical OK,SSYS
        INTEGER IFLGCHK,OLOC(NPLMAX,NPLMAX)
        CHARACTER*1 REP
        CHARACTER*32 NAME
        CHARACTER*80 INPLFILE
        CHARACTER*65 FPLANETS
        INTEGER NBOD,IALPHA,I,J,IP1,IP2,IUFLG,ICFLG,IRFLG
        INTEGER NNN

        OK=.FALSE.
        DO WHILE(.NOT.OK)
          WRITE(*,*) ' Units menu (G=1) : '
          WRITE(*,*) '       0 ==> Solar masses and AU '
          WRITE(*,*) '       1 ==> AU and years '
          READ(*,*) IUFLG
          OK = ((IUFLG.EQ.0).OR.(IUFLG.EQ.1))
        END DO

        OK=.FALSE.
        DO WHILE(.NOT.OK)
          WRITE(*,*) ' Coordinate Menu: '
          WRITE(*,*) '       0 ==> Initial plane '
          WRITE(*,*) '       1 ==> Invariable plane '
          READ(*,*) ICFLG
          OK = ((ICFLG.EQ.0).OR.(ICFLG.EQ.1))
        END DO

        WRITE(*,*)' Give number of massive bodies in the system '
        READ(*,*)NBOD
        DO I=1,NBOD
          WRITE(*,*)' Give mass of body #',i,' in solar masses '
          READ(*,*)MASS(I)
        END DO

        IF(IUFLG.EQ.1) then
          DO I=1,NBOD
            MASS(I) = MASS(I)*SMASSYR
          END DO
        END IF

        XJ(1) = 0.0d0
        YJ(1) = 0.0d0
        ZJ(1) = 0.0d0
        VXJ(1) = 0.0d0
        VYJ(1) = 0.0d0
        VZJ(1) = 0.0d0

        DO I=2,NBOD
          WRITE(*,*)' Information relative to orbit #',I-1,':'
          WRITE(*,*)' Give situation of bodies relative to that orbit.'
          WRITE(*,*)' 0 = foreign,   1 = satellite,  -1 = center'
          READ(*,*)(OLOC(I,J),J=1,NBOD)

          WRITE(*,*)
     &'Give now the orbital parameters (a,e,i,Omega,omega,M)'
          READ(*,*)A,E,INC,CAPOM,OMEGA,CAPM

          CAPM = DR*CAPM
          OMEGA = DR*OMEGA
          CAPOM = DR*CAPOM
          INC = DR*INC
         
          GM = 0.d0
          DO J = 1,NBOD
            IF (OLOC(I,J).NE.0) GM = GM+MASS(J)
          END DO
          IALPHA = -1
          CALL ORBEL_EL2XV(GM,IALPHA,A,E,INC,CAPOM,OMEGA,CAPM,
     &          XJ(I),YJ(I),ZJ(I),VXJ(I),VYJ(I),VZJ(I))
        END DO


        call hierarcheck(nbod,oloc,ok)
        write(*,*)' '
        if (.not.ok) then
          write(*,*),'stopping...'
          stop
        end if
        write(*,*)'Hierarchical structure ok !'
        write(*,*)' '

          LCLOSE = .FALSE.
          do i=1,nbod
            RPLSQ(i) = 0.d0
          end do
c... Compute eta's and mu's: center and satellite masses for orbits
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
          umat(i,1) = 1.0d0
        end do
        do j = 2,nbod
          vsat = eta(j)/(mu(j)+eta(j))
          vcen = -mu(j)/(mu(j)+eta(j))
          do i = 1,nbod
            if (oloc(j,i).eq.1) umat(i,j) = vsat
            if (oloc(j,i).eq.-1) umat(i,j) = vcen
          end do
        end do         

        CALL COORD_G2B(NBOD,UMAT,MASS,XJ,YJ,ZJ,VXJ,VYJ,VZJ,XB,YB,ZB,VXB,
     &     VYB,VZB)

        IF(ICFLG.EQ.1) CALL INVAR(NBOD,MASS,XB,YB,ZB,VXB,VYB,VZB)
        print*,xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1)

        CALL COORD_B2G(NBOD,MAT,MASS,XB,YB,ZB,VXB,VYB,
     &     VZB,XJO,YJO,ZJO,VXJO,VYJO,VZJO)

	WRITE(*,*) ' '
	WRITE(*,*) 'Enter name of planet data file : '
	READ(*,'(a)') INPLFILE
        IFLGCHK = 0
         print*,vzjo(1),vzjo(2),vzjo(3)
	CALL IO_DUMP_PL_HJS(INPLFILE,NBOD,OLOC,MASS,UMAT,
     &           XJO,YJO,ZJO,VXJO,VYJO,VZJO,LCLOSE,IFLGCHK,RPLSQ)


        END
C
C--------------------------------------------------------------------
C               For invariable plane
C--------------------------------------------------------------------
C
        SUBROUTINE INVAR(NBOD,MASS,X,Y,Z,VX,VY,VZ)

        include '../../swift.inc'

        INTEGER I,NBOD

        REAL*8 MASS(NBOD),X(NBOD),Y(NBOD),Z(NBOD),
     &                    VX(NBOD),VY(NBOD),VZ(NBOD),
     &                    XT,YT,ZT,VXT,VYT,VZT,
     &                    C1,C2,C3,C,BOT,A(3,3)

        C1 = 0.0d0
        C2 = 0.0d0
        C3 = 0.0d0
        DO I = 1,NBOD
          C1 = C1 + MASS(I)*(Y(I)*VZ(I)-Z(I)*VY(I))
          C2 = C2 + MASS(I)*(Z(I)*VX(I)-X(I)*VZ(I))
          C3 = C3 + MASS(I)*(X(I)*VY(I)-Y(I)*VX(I))
        END DO

        WRITE(*,*) ' C1,C2,C3 ',C1,C2,C3

        C = SQRT( C1*C1 + C2*C2 + C3*C3 )
        C1 = C1/C
        C2 = C2/C
        C3 = C3/C
        BOT = 1.0D0/(1.0D0 + C3)

        A(1,1) = 1.0D0 - C1*C1*BOT
        A(1,2) = -1.0D0*C1*C2*BOT
        A(1,3) = -1.0D0*C1

        A(2,1) = -1.0D0*C1*C2*BOT
        A(2,2) = 1.0D0 - C2*C2*BOT
        A(2,3) = -1.0D0*C2

        A(3,1) = C1
        A(3,2) = C2
        A(3,3) = C3

        DO I=1,NBOD
          XT = A(1,1)*X(I) + A(1,2)*Y(I) + A(1,3)*Z(I) 
          YT = A(2,1)*X(I) + A(2,2)*Y(I) + A(2,3)*Z(I) 
          ZT = A(3,1)*X(I) + A(3,2)*Y(I) + A(3,3)*Z(I) 
          X(I) = XT
          Y(I) = YT
          Z(I) = ZT
          VXT = A(1,1)*VX(I) + A(1,2)*VY(I) + A(1,3)*VZ(I) 
          VYT = A(2,1)*VX(I) + A(2,2)*VY(I) + A(2,3)*VZ(I) 
          VZT = A(3,1)*VX(I) + A(3,2)*VY(I) + A(3,3)*VZ(I) 
          VX(I) = VXT
          VY(I) = VYT
          VZ(I) = VZT
        END DO

C..    Check

        C1 = 0.0d0
        C2 = 0.0d0
        C3 = 0.0d0
        DO I = 1,NBOD
          C1 = C1 + MASS(I)*(Y(I)*VZ(I)-Z(I)*VY(I))
          C2 = C2 + MASS(I)*(Z(I)*VX(I)-X(I)*VZ(I))
          C3 = C3 + MASS(I)*(X(I)*VY(I)-Y(I)*VX(I))
        END DO

        WRITE(*,*) ' C1,C2,C3 ',C1,C2,C3

        END

C
C-----------------------------------------------------------------------
C        Subroutine dedicated to verify the hierarchical structure
C        Returns logical ok=true if the hierarchical structure is valid
C-----------------------------------------------------------------------
C
        subroutine hierarcheck(NBOD,oloc,ok)

        include '../../swift.inc'

        integer i,j,k,l,nnn,nbod,oloc(NPLMAX,NPLMAX)
        logical test,ok

        write(*,*)' '
        OK = .true.

c... We test over all pairs of orbits. Either they are foreign to
c... each other or one is inner to the other one
        DO I = 2,nbod-1
          do j = i+1,nbod
            test = .true.
c... We first test whether orbits i and j are foreign
            do k = 1,nbod
              test = test.and.(oloc(i,k)*oloc(j,k).eq.0)
            end do
            if (test) then
              write(*,*)'Orbits',i-1,' and',j-1,' are foreign'
            else
c... From here on they are not foreign. Now check whether orbit i
c... is inner to orbit j, i.e. all centers and satellites of orbit i
c... must fall in the centers or the satellites of orbit j.
              l = 1
              do while(oloc(i,l).eq.0)
                l = l+1
              end do
              nnn = oloc(j,l)
              if (nnn.ne.0) then
                test = .true.
                do k = l+1,nbod
                  if (oloc(i,k).ne.0) test =
     &                               test.and.(oloc(j,k).eq.nnn)
                end do
                if (test) write(*,*)'Orbit',i-1,
     &                                     ' is inner to orbit',j-1
              end if
              if ((nnn.eq.0).or.(.not.test)) then
c... Now orbit i is not inner to orbit j. We check whether orbit j
c... is inner to orbit i
                l = 1
                do while(oloc(j,l).eq.0)
                  l = l+1
                end do
                nnn = oloc(i,l)
                if (nnn.ne.0) then
                  test = .true.
                  do k = l+1,nbod
                    if (oloc(j,k).ne.0) test = 
     &                               test.and.(oloc(i,k).eq.nnn)
                  end do
                  if (test) then
                    write(*,*)'Orbit',j-1,' is inner to orbit',i-1
                  else
c... If all fails then we are sure that orbits i and j do not
c... fulfill the hierarchy rules
                    write(*,*)'Structure problem with orbits',
     &                                             i-1,' and',j-1
                    ok = .false.
                  end if
                else
                  write(*,*)'Structure problem with orbits',
     &                                             i-1,' and',j-1
                  ok=.false. 
                end if
              end if
            end if
          end do
        end do
        
        end
