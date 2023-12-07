C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C 
        PROGRAM GEN_AE_HJS

        include '../../swift.inc'

        REAL*8  SMASSYR
        PARAMETER (SMASSYR = TWOPI*TWOPI )
	INTEGER*4	NDATA          ! Taille du tableau de donnees.
        INTEGER IFLGCHK,NBOD,NTP,NTPT
        LOGICAL*2 LCLOSE
        REAL*8 XJT(NTPMAX),YJT(NTPMAX),ZJT(NTPMAX),
     &         VXJT(NTPMAX),VYJT(NTPMAX),VZJT(NTPMAX),
     &         RPLSQ(NPLMAX),MATR(3,3),UMAT(NPLMAX,NPLMAX),
     &         XJ(NPLMAX),YJ(NPLMAX),ZJ(NPLMAX),
     &         VXJ(NPLMAX),VYJ(NPLMAX),VZJ(NPLMAX),
     &         ETA(NPLMAX),MU(NPLMAX),MAT(NPLMAX,NPLMAX),
     &         VXH(NPLMAX),VYH(NPLMAX),VZH(NPLMAX),
     &         XB(NPLMAX),YB(NPLMAX),ZB(NPLMAX),
     &         VXB(NPLMAX),VYB(NPLMAX),VZB(NPLMAX),
     &         A,E,INC,CAPOM,OMEGA,CAPM,AMIN,AMAX,DA,
     &         RSTAT(NTPMAX,NSTATR),MSUN,MASS(NPLMAX),
     &         CIMAX,IMAX,DE,EMAX,EMIN,XX,AMINX,AMAXX,
     &         XT,YT,ZT,VXT,VYT,VZT
        REAL*8 ETATPT,MATPT(NPLMAX),UMATPT(NPLMAX),XJTEMP(NPLMAX),
     &         YJTEMP(NPLMAX),ZJTEMP(NPLMAX),VXJTEMP(NPLMAX),
     &         VYJTEMP(NPLMAX),VZJTEMP(NPLMAX),MASSTEMP(NPLMAX),
     &         XBTEMP(NPLMAX),YBTEMP(NPLMAX),ZBTEMP(NPLMAX),
     &         VXBTEMP(NPLMAX),VYBTEMP(NPLMAX),VZBTEMP(NPLMAX),
     &         UMPART(NPLMAX,NPLMAX),MATP(NPLMAX,NTPMAX),
     &         UMATP(NPLMAX,NTPMAX)
        real*8 xbt,ybt,zbt,vxbt,vybt,vzbt
        REAL*4 SEED

        INTEGER ISTAT(NTPMAX,NSTAT),OLOC(NPLMAX,NPLMAX)
        INTEGER ORBCT(NPLMAX),NC,KK,JJ,OLOCTT(NPLMAX),
     &          OLOCT(NPLMAX,NTPMAX)
        INTEGER IALPHA,I,J,II,K,ISEED,IUFLG,LLM,IND
        LOGICAL OK,OKC,SAT,CEN
        CHARACTER*80 INTPFILE,INPLFILE
        CHARACTER*1 PR

        OK=.FALSE.
        DO WHILE(.NOT.OK)
          WRITE(*,*) ' Units menu (G=1) : '
          WRITE(*,*) '       0 ==> Solar masses and AU '
          WRITE(*,*) '       1 ==> AU and years '
          READ(*,*) IUFLG
          OK = ((IUFLG.EQ.0).OR.(IUFLG.EQ.1))
        END DO

c Prompt and read name of planet data file
	WRITE(*,*) ' '
	WRITE(*,*) 'Enter name of planet data file : '
	READ(*,'(a)') INPLFILE
	CALL IO_INIT_PL_HJS(INPLFILE,LCLOSE,IFLGCHK,NBOD,OLOC,
     &       MASS,ETA,MU,MAT,UMAT,XJ,YJ,ZJ,VXJ,VYJ,VZJ,RPLSQ)

       
        DO I=1,NBOD
          ORBCT(I)=0.0d0
        END DO
        DO I=1,2**NBOD-1
          ORBCT(1)=ORBCT(1)-1
          DO J=1,NBOD-1
            IF (ORBCT(J).EQ.-2) THEN
              ORBCT(J)=0
              ORBCT(J+1)=ORBCT(J+1)-1
            END IF
          END DO
          call hierarchktp(nbod,oloc,orbct,ok)
          IF (OK) WRITE(*,*) (ORBCT(J),J=1,NBOD)
        END DO            
          
        stop


        OKC = .TRUE.
        WRITE(*,*) ' Input ISEED: large and odd '
        READ(*,*) ISEED
        SEED=FLOAT(ISEED)
        CALL SRAND(ISEED)
        NTPT = 0
        II = 0
        OPEN(12,FILE='gen_ae.out',status='UNKNOWN')
        DO WHILE(OKC)
          WRITE(*,*) '------------ Set of test particles --------------',
          WRITE(*,*) ' Input size of set '
          READ(*,*) NTP
          OKC = (NTP.NE.0)
          IF (OKC) THEN
            NTPT = NTPT+NTP
            WRITE(*,*)' Give orbct : centers for tp '
            READ(*,*) (ORBCT(I),I=1,NBOD)
            write(*,*),'Centers : ',(orbct(j),j=1,nbod)

c... Check the validity
            call hierarchktp(nbod,oloc,orbct,ok)
            write(*,*)' '
            if (.not.ok) then
              write(*,*),'stopping...'
              stop
            end if
            write(*,*)'Hierarchical structure ok !'
            write(*,*)' '

            etatpt = 0.0d0   
            
            do j = 1,nbod
              matpt(j) = 0.0d0
              umatpt(j) = 0.0d0
            end do
            do j = 1,nbod
              if (orbct(j).eq.-1) etatpt = etatpt+mass(j)
            end do
            do j = 1,nbod
              if (orbct(j).eq.-1) matpt(j) = -mass(j)/etatpt
            end do
            do k = 1,nbod
              if (orbct(k).eq.-1) then
                do j=1,nbod     
                  umatpt(j) = umatpt(j) - matpt(k)*umat(k,j)
                end do
              end if
            end do


c... Computation of the oloct array
c            do j=1,nbod
c              oloctt(j) = 0
c            end do
c            do j = 2,nbod
c              sat = .true.
c              cen = .true.
c              do i = 1,nbod
c                if (orbct(i).eq.-1) then                  
c                  sat = sat.and.(oloc(j,i).eq.1)
c                  cen = cen.and.(oloc(j,i).eq.-1)
c                end if
c              end do
c              if (sat) oloctt(j) = 1
c              if (cen) oloctt(j) = -1
c            end do
c            print*,'oloct',(oloctt(j),j=1,nbod)



            WRITE(*,*) ' Input ecc of particles '
            READ(*,*) EMIN,EMAX
            WRITE(*,*) ' Input max inc of particles (in deg) '
            READ(*,*) IMAX
            WRITE(*,*) ' Input AMIN and AMAX '
            READ(*,*) AMIN,AMAX
            WRITE(*,*) ' Input power index '
            READ(*,*) XX
            AMINX=AMIN**(XX+1.)
            AMAXX=AMAX**(XX+1.)
            E=EMIN
            CIMAX = COS(IMAX/DEGRAD)

c... Here we need to calculate the plane perp. to the local angular momentum 
c...    of the centers of the tp's. The orbital elements of the tp's are given 
c...    relative to this base plane.
            nc = 0
            do j = 1,nbod
              if (orbct(j).eq.-1) then
                nc = nc+1
                masstemp(nc) = mass(j)
                xjtemp(nc) = xj(j)
                yjtemp(nc) = yj(j)
                zjtemp(nc) = zj(j)
                vxjtemp(nc) = vxj(j)
                vyjtemp(nc) = vyj(j)
                vzjtemp(nc) = vzj(j)
              end if
            end do
    
c... Of course we calculate the midplane if there is more than one center...
            if (nc.gt.1) then 

c... Now we extract a submatrix from the umat matrix 
              do j = 1,nbod
                do i = 1,nbod
                  umpart(i,j) = 0.0d0
                end do
              end do
              do j = 1,nc
                umpart(j,1) = 1.0d0
              end do
              jj = 1
              do j = 2,nbod
                ok = .true.
                do k = 1,nbod
                  if ((matpt(k).eq.0.0d0).and.(oloc(j,k).ne.0))
     &                                                ok=.false.
                end do
                if (ok) then
                  jj = jj+1
                  kk = 0
                  do k = 1,nbod
                    if (matpt(k).ne.0.0d0) then
                      kk = kk+1
                      umpart(kk,jj) = umat(k,j)
                    end if
                  end do
                end if
              end do
                   
              call coord_g2b(nc,umpart,masstemp,xjtemp,yjtemp,zjtemp,
     &                  vxjtemp,vyjtemp,vzjtemp,xbtemp,ybtemp,zbtemp,
     &                  vxbtemp,vybtemp,vzbtemp)
              CALL INVAR(NC,MASSTEMP,XBTEMP,YBTEMP,ZBTEMP,VXBTEMP,
     &                 VYBTEMP,VZBTEMP,MATR)
            END IF

c... Now do the particles

            WRITE(*,*) ' Test particles: '
            WRITE(*,*) '   A      E      I    CAPOM  OMEGA    M '
            DO I=1,NTP
              II = II+1
              IALPHA = -1
              INC = ACOS( 1.0D0 - (1.0D0-CIMAX)*RAND())
              CAPM = TWOPI*RAND()
              OMEGA = TWOPI*RAND()
              CAPOM = TWOPI*RAND()
              A = (AMINX + (AMAXX - AMINX)*RAND())**(1./(XX+1.))
              E = EMIN + (EMAX - EMIN)*RAND()
              WRITE(*,*) II,(ORBCT(J),J=1,NBOD),SNGL(A),SNGL(E),
     &                   SNGL(INC),SNGL(CAPOM),SNGL(OMEGA),SNGL(CAPM)
              WRITE(12,*) II,(ORBCT(J),J=1,NBOD),SNGL(A),SNGL(E),
     &                   SNGL(INC),SNGL(CAPOM),SNGL(OMEGA),SNGL(CAPM)
              CALL ORBEL_EL2XV(ETATPT,IALPHA,A,E,INC,CAPOM,OMEGA,
     &       CAPM,XJT(II),YJT(II),ZJT(II),VXJT(II),VYJT(II),VZJT(II))
              IF (NC.GT.1) THEN
                XT = MATR(1,1)*XJT(II)+MATR(2,1)*YJT(II)
     &                                      +MATR(3,1)*ZJT(II) 
                YT = MATR(1,2)*XJT(II)+MATR(2,2)*YJT(II)
     &                                      +MATR(3,2)*ZJT(II) 
                ZT = MATR(1,3)*XJT(II)+MATR(2,3)*YJT(II)
     &                                      +MATR(3,3)*ZJT(II) 
                XJT(II) = XT
                YJT(II) = YT
                ZJT(II) = ZT
                VXT = MATR(1,1)*VXJT(II)+MATR(2,1)*VYJT(II)
     &                                      +MATR(3,1)*VZJT(II) 
                VYT = MATR(1,2)*VXJT(II)+MATR(2,2)*VYJT(II)
     &                                      +MATR(3,2)*VZJT(II) 
                VZT = MATR(1,3)*VXJT(II)+MATR(2,3)*VYJT(II)
     &                                      +MATR(3,3)*VZJT(II) 
                VXJT(II) = VXT
                VYJT(II) = VYT
                VZJT(II) = VZT
              END IF

              do j = 1,nbod
c                oloct(j,II) = OLOCTT(j)
                matp(j,II) = MATPT(j)
                umatp(j,II) = UMATPT(j)
              end do

              DO J=1,NSTAT
                ISTAT(II,J) = 0
              END DO
              DO J=1,NSTATR
                RSTAT(II,J) = 0.0D0
              END DO
            END DO
          END IF
        END DO

        close(12)
        WRITE(*,*) 'Enter name of test particle data file : '
        READ(*,'(a)') INTPFILE
        print*,intpfile

        call io_dump_tp_hjs(INTPFILE,NBOD,NTPT,matp,xjt,yjt,zjt,
     &                       vxjt,vyjt,vzjt,istat,rstat)

        END PROGRAM GEN_AE_HJS

C
C--------------------------------------------------------------------
C               For invariable plane 
C--------------------------------------------------------------------
C
        SUBROUTINE INVAR(NBOD,MASS,X,Y,Z,VX,VY,VZ,A)

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

        END

C
C-----------------------------------------------------------------------
C        Subroutine dedicated to verify the hierarchical structure
C        of the tp's we add. Returns logical ok=true if the hierarchical
C        structure is valid
C-----------------------------------------------------------------------
C
        subroutine hierarchktp(nbod,oloc,orbct,ok)

        include '../../swift.inc'

        integer i,k,l,nnn,nbod,oloc(NPLMAX,NPLMAX),orbct(NPLMAX)
        logical test,ok

c        write(*,*)' '
        OK = .true.
c... First check that the orbct array is made only of 0's and -1's
        do i=1,nbod
          ok = ok.and.((orbct(i).eq.0).or.(orbct(i).eq.-1))
        end do
        if (.not.ok) then
c          write(*,*)'The orbit given is not a valid tp orbit'
          return
        end if
c... Adding a tp means adding an orbit. We check the validity of the
c... the new orbit in the hierarchical structure of the system, i.e.
c... with all orbits i between massive bodies
        DO I = 2,nbod
          test = .true.
c... We first test whether tp orbit and orbit i are foreign
          do k = 1,nbod
            test = test.and.(oloc(i,k)*orbct(k).eq.0)
          end do
          if (test) then
c            write(*,*)'tp orbit and orbit',i-1,' are foreign'
          else
c... From here on they are not foreign. Now check whether tp orbit
c... is inner to orbit i, i.e. all centers of tp orbit
c... must fall in the centers or the satellites of orbit i.
            l = 1
            do while(orbct(l).eq.0)
              l = l+1
            end do
            nnn = oloc(i,l)
            if (nnn.ne.0) then
              test = .true.
              do k = l+1,nbod
                if (orbct(k).eq.-1) test =
     &                               test.and.(oloc(i,k).eq.nnn)
              end do
              if (test) then
c                write(*,*)'tp orbit is inner to orbit',i-1
              end if
            end if
            if ((nnn.eq.0).or.(.not.test)) then
c... Now tp orbit is not inner to orbit i. We check whether orbit i
c... is inner to tp orbit
              l = 1
              do while(oloc(i,l).eq.0)
                l = l+1
              end do
              nnn = orbct(l)
              if (nnn.ne.0) then
                test = .true.
                do k = l+1,nbod
                  if (oloc(i,k).ne.0) test = 
     &                               test.and.(orbct(k).eq.nnn)
                end do
                if (test) then
c                  write(*,*)'Orbit',i-1,' is inner to tp orbit'
                else
c... If all fails then we are sure that orbits i and j do not
c... fulfill the hierarchy rules
c                  write(*,*)'Structure problem relative to orbit',i-1
                  ok = .false.
                end if
              else
c                write(*,*)'Structure problem relative to orbit',i-1
                ok=.false. 
              end if
            end if
          end if
        end do
        
        end
