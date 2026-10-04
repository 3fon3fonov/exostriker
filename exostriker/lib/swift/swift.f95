! In-memory SWIFT interface. Numerical initialization and Jacobi/MVS helpers
! are derived from the original swift.f95 prototype; file-based drivers and
! followers are replaced by direct state sampling. Python serializes sessions.
module swift_memory_state
    use, intrinsic :: ieee_arithmetic
    use swift_memory_support
    include 'swift.inc'
    integer :: nbod=0, nbodm=0, method=0, first_mvs=0, gr_counter=0, gr_stride=1
    real(8) :: t=0.d0, stop_time=0.d0, dt=0.d0, interval=0.d0, tout=0.d0, eoff=0.d0
    real(8) :: mass(NTPMAX),xj(NTPMAX),yj(NTPMAX),zj(NTPMAX)
    real(8) :: vxj(NTPMAX),vyj(NTPMAX),vzj(NTPMAX)
    real(8) :: xh(NTPMAX),yh(NTPMAX),zh(NTPMAX),vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)
    real(8) :: rpl(NTPMAX),rhill(NTPMAX)
    logical :: initial=.false., finished=.true.
contains
    subroutine sample(energy,orbit,valid)
        real(8), intent(out) :: energy(:),orbit(:,:)
        integer, intent(out) :: valid(:)
        real(8) ke,pot,total,angular(3),gm,a,e,inc,node,omega,anomaly
        real(8) elements(NPLMAX,11)
        integer j,k,ialpha
        energy=0.d0
        orbit=0.d0
        valid=0
        if(nbod == 1) then
            energy=[t,eoff,0.d0,0.d0,0.d0]
            if(.not.ieee_is_finite(eoff)) swift_status=4
            return
        endif
        call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,xh,yh,zh,vxh,vyh,vzh)
        if(swift_status /= 0) return
        call anal_energy(nbod,mass,0.d0,0.d0,xh,yh,zh,vxh,vyh,vzh,ke,pot,total,angular)
        if(swift_status /= 0) return
        energy = [t,total+eoff,angular]
        gm=mass(1)
        do j=2,nbod
            gm=gm+mass(j)
            call orbel_xv2el(xj(j),yj(j),zj(j),vxj(j),vyj(j),vzj(j),gm, &
                            ialpha,a,e,inc,node,omega,anomaly)
            if(swift_status /= 0) return
            elements(j,:) = [t,real(j,8),a,e,inc*DEGRAD,node*DEGRAD, &
                            omega*DEGRAD,anomaly*DEGRAD,a*(1.d0-e),a*(1.d0+e),mass(j)/mass(1)]
            if(method /= 1) elements(j,2)=-real(j,8)
        enddo
        do k=1,ntracked
            j=tracked_slot(k)
            if(j >= 2 .and. j <= nbod) then
                valid(k)=1
                orbit(k,:)=elements(j,:)
            endif
        enddo
        if(.not.all(ieee_is_finite(energy)) .or. .not.all(ieee_is_finite(orbit))) swift_status=4
    end subroutine

    subroutine step()
        integer first,old_nbod,isenc,mergelst(2,NTPMAX),mergecnt,j,ialpha
        integer(2) iecnt(NTPMAX)
        real(8) gm,a,e,inc,node,omega,anomaly,corr
        logical(2) close_encounters
        close_encounters=.true.
        if(method == 1) then
            first=0
            call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,xh,yh,zh,vxh,vyh,vzh)
            if(swift_status /= 0) return
            call symba5_step_pl(first,t,nbod,nbodm,mass,0.d0,0.d0, &
                xh,yh,zh,vxh,vyh,vzh,dt,close_encounters,rpl,isenc, &
                mergelst,mergecnt,iecnt,eoff,rhill,1.d-40)
            if(swift_status /= 0) return
            call coord_h2j(nbod,mass,xh,yh,zh,vxh,vyh,vzh,xj,yj,zj,vxj,vyj,vzj)
            if(swift_status /= 0) return
            t=t+dt
            old_nbod=nbod
            ! Preserve the original SyMBA/Jacobi driver's discard convention.
            call discard_massive5(t,dt,nbod,mass,xj,yj,zj,vxj,vyj,vzj, &
                0.0001d0,50.d0,50.d0,-1.d0,close_encounters,rpl,rhill, &
                isenc,mergelst,mergecnt,iecnt,eoff,first)
            if(swift_status /= 0) return
            if(nbod /= old_nbod) call symba5_nbodm(nbod,mass,1.d-40,nbodm)
        else
            if(method == 3 .and. gr_counter == gr_stride) then
                gm=mass(1)
                do j=2,nbod
                    gm=gm+mass(j)
                    call orbel_xv2el(xj(j),yj(j),zj(j),vxj(j),vyj(j),vzj(j),gm, &
                                    ialpha,a,e,inc,node,omega,anomaly)
                    if(swift_status /= 0) return
                    if(ialpha /= -1) then
                        ! The orbit-averaged GR precession requires a bound
                        ! ellipse. Return completed histories as an early stop.
                        swift_status=5
                        return
                    endif
                    call gr_corr(a,e,gm,corr,dt)
                    omega=omega+corr*gr_stride
                    call orbel_el2xv(gm,ialpha,a,e,inc,node,omega,anomaly, &
                                    xj(j),yj(j),zj(j),vxj(j),vyj(j),vzj(j))
                    if(swift_status /= 0) return
                enddo
                ! GR changes positions as well as velocities. The next MVS
                ! half-kick must recompute its saved acceleration cache.
                first_mvs=0
                gr_counter=0
            endif
            if(method == 3) gr_counter=gr_counter+1
            call step_kdk_pl_j(first_mvs,nbod,mass,0.d0,0.d0,xj,yj,zj,vxj,vyj,vzj,dt)
            if(swift_status /= 0) return
            t=t+dt
            ! MVS's original discard() handles test particles only (ntp=0).
        endif
        if(.not.all(ieee_is_finite(xj(:nbod))) .or. .not.all(ieee_is_finite(yj(:nbod))) .or. &
           .not.all(ieee_is_finite(zj(:nbod))) .or. .not.all(ieee_is_finite(vxj(:nbod))) .or. &
           .not.all(ieee_is_finite(vyj(:nbod))) .or. .not.all(ieee_is_finite(vzj(:nbod)))) swift_status=4
    end subroutine
end module

subroutine swift_begin(stellar_mass,planets,nplan,integrator,tmax,timestep,output_step,gr_step,timeout,status)
    use swift_memory_state
    implicit none
    integer, intent(in) :: nplan,integrator,gr_step
    real(8), intent(in) :: stellar_mass,planets(nplan,7),tmax,timestep,output_step,timeout
    integer, intent(out) :: status
    real(8) coordinates(3*nplan+3,3),params(nplan,9)
    integer j
    status=1
    finished=.true.
    if(nplan < 1 .or. nplan > 50) return
    if(integrator < 1 .or. integrator > 3 .or. gr_step < 1) return
    if(.not.all(ieee_is_finite([stellar_mass,tmax,timestep,output_step,timeout]))) return
    if(min(stellar_mass,tmax,timestep,output_step,timeout) <= 0.d0) return
    if(.not.all(ieee_is_finite(planets))) return
    if(any(planets(:,1:2) <= 0.d0) .or. any(planets(:,3)<0.d0) .or. any(planets(:,3)>=1.d0)) return
    call memory_reset(nplan,timeout)
    mass=0.d0
    xj=0.d0; yj=0.d0; zj=0.d0
    vxj=0.d0; vyj=0.d0; vzj=0.d0
    rpl=0.d0; rhill=0.d0
    call geninit_j3_in_days(1,stellar_mass,nplan,1.d0,planets,params,coordinates)
    status=swift_status
    if(status /= 0) return
    nbod=nplan+1
    do j=1,nbod
        mass(j)=coordinates(3*j-2,1)
        rhill(j)=coordinates(3*j-2,2)
        rpl(j)=coordinates(3*j-2,3)
        xj(j)=coordinates(3*j-1,1); yj(j)=coordinates(3*j-1,2); zj(j)=coordinates(3*j-1,3)
        vxj(j)=coordinates(3*j,1); vyj(j)=coordinates(3*j,2); vzj(j)=coordinates(3*j,3)
    enddo
    call symba5_nbodm(nbod,mass,1.d-40,nbodm)
    method=integrator
    t=0.d0; stop_time=tmax; dt=timestep; interval=output_step; tout=interval; eoff=0.d0
    first_mvs=0; gr_counter=0; gr_stride=gr_step
    initial=.true.; finished=.false.
    call swift_poll()
    status=swift_status
end subroutine

subroutine swift_advance(nplan,max_steps,capacity,energy,orbits,valid,nfilled,done,status)
    use swift_memory_state
    implicit none
    integer, intent(in) :: nplan,max_steps,capacity
    real(8), intent(out) :: energy(capacity,5),orbits(capacity,nplan,11)
    integer, intent(out) :: valid(capacity,nplan),nfilled,done,status
    integer steps
    energy=0.d0; orbits=0.d0; valid=0
    nfilled=0; steps=0; done=0
    if(nplan /= ntracked .or. capacity < 1 .or. max_steps < 1) then
        status=1
        return
    endif
    call swift_poll()
    if(initial .and. swift_status == 0) then
        nfilled=1
        call sample(energy(1,:),orbits(1,:,:),valid(1,:))
        initial=.false.
    endif
    do while(.not.finished .and. nfilled < capacity .and. steps < max_steps .and. swift_status == 0)
        if(t > stop_time .or. (method == 1 .and. nbod <= 1)) then
            finished=.true.
            exit
        endif
        if(mod(steps,64) == 0) call swift_poll()
        if(swift_status /= 0) exit
        if(t+dt <= t) then
            swift_status=4
            exit
        endif
        call step()
        if(swift_status /= 0) exit
        steps=steps+1
        if(t >= tout) then
            nfilled=nfilled+1
            call sample(energy(nfilled,:),orbits(nfilled,:,:),valid(nfilled,:))
            tout=tout+interval
        endif
    enddo
    call swift_poll()
    if(finished) done=1
    status=swift_status
end subroutine

! Corrected initializer from the original prototype: nbod means planets and
! is never mutated. Coordinates include the star and stay double precision.
subroutine geninit_j3_in_days(iuflg,mstar0,nbod,frho3,planets_params,planets_results,coordinates)
    use swift_memory_support, only: swift_status
    include 'swift.inc'
    integer, intent(in) :: iuflg,nbod
    real(8), intent(in) :: mstar0,frho3,planets_params(nbod,7)
    real(8), intent(out) :: planets_results(nbod,9),coordinates(3*nbod+3,3)
    real(8) factor,mpl,gm,a,e,inc,omega,node,anomaly,radius,hill
    integer j
    coordinates=0.d0
    planets_results=0.d0
    factor=1.d0
    if(iuflg == 1) factor=TWOPI*TWOPI/(365.25d0*365.25d0)
    coordinates(1,1)=mstar0*factor
    gm=coordinates(1,1)
    do j=1,nbod
        mpl=planets_params(j,1)
        a=planets_params(j,2); e=planets_params(j,3)
        inc=planets_params(j,4)/DEGRAD; omega=planets_params(j,5)/DEGRAD
        node=planets_params(j,6)/DEGRAD; anomaly=planets_params(j,7)/DEGRAD
        radius=frho3*(1.5d0*mpl*1.989d33/TWOPI)**0.3333333333d0/1.49597892d13
        hill=a*(mpl/(3.d0*mstar0))**0.3333333333d0
        gm=gm+mpl*factor
        coordinates(3*j+1,:)=[mpl*factor,hill,radius]
        call orbel_el2xv(gm,-1,a,e,inc,node,omega,anomaly, &
            coordinates(3*j+2,1),coordinates(3*j+2,2),coordinates(3*j+2,3), &
            coordinates(3*j+3,1),coordinates(3*j+3,2),coordinates(3*j+3,3))
        if(swift_status /= 0) return
        planets_results(j,:)=[mpl*factor,a,e,inc,omega,node,anomaly,radius,hill]
    enddo
end subroutine

subroutine step_kdk_pl_j(i1st, nbod, mass, j2rp2, j4rp4, &
        xj, yj, zj, vxj, vyj, vzj, dt)

    use swift_memory_support, only: swift_status
    include 'swift.inc'

    !...  Inputs Only:
    integer nbod, i1st
    real(8) mass(nbod), dt, j2rp2, j4rp4

    !...  Inputs and Outputs:
    real(8) xj(nbod), yj(nbod), zj(nbod)
    real(8) vxj(nbod), vyj(nbod), vzj(nbod)

    !...  Internals:
    real(8) dth
    real(8) axj(NPLMAX), ayj(NPLMAX), azj(NPLMAX)
    real(8) xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)

    save axj, ayj, azj, xh, yh, zh     ! Note this !!

    !----
    !...  Executable code

    dth = 0.5d0 * dt

    if(i1st==0) then
        !...      Convert to helio coords
        call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
                xh, yh, zh, vxh, vyh, vzh)
    if(swift_status /= 0) return
        !...     Get the accelerations in Jacobi frame. if frist time step
        call getaccj(nbod, mass, j2rp2, j4rp4, xh, yh, zh, &
                xj, yj, zj, axj, ayj, azj)
    if(swift_status /= 0) return
        i1st = 1    ! turn this off
    endif

    !...  Apply a Jacobi kick for a half dt
    call kickvh(nbod, vxj, vyj, vzj, axj, ayj, azj, dth)
    if(swift_status /= 0) return

    !..   Drift in Jacobi coords for the full step
    call drift2(nbod, mass, xj, yj, zj, vxj, vyj, vzj, dt)
    if(swift_status /= 0) return

    !...  After drift, compute helio. xh and vh for acceleration calculations
    call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
            xh, yh, zh, vxh, vyh, vzh)
    if(swift_status /= 0) return

    !...  Get the accelerations in Jacobi frame.
    call getaccj(nbod, mass, j2rp2, j4rp4, xh, yh, zh, xj, yj, zj, axj, ayj, azj)
    if(swift_status /= 0) return

    !...  Apply a Jacobi kick for a half dt
    call kickvh(nbod, vxj, vyj, vzj, axj, ayj, azj, dth)
    if(swift_status /= 0) return

    return
end
! step_kdk_pl_j
!---------------------------------------------------------------------

!*************************************************************************
!                        DRIFT2.F
!*************************************************************************
! This subroutine loops thorugh the particles and calls the danby routine
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xj,yj,zj      ==>  initial position in jacobi coord
!                                    (real arrays)
!                 vxj,vyj,vzj   ==>  initial position in jacobi coord
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 xj,yj,zj      ==>  final position in jacobi coord
!                                       (real arrays)
!                 vxj,vyj,vzj   ==>  final position in jacobi coord
!                                       (real arrays)
!
! Remarks: Based on drift
! Authors:  Man Hoi Lee
! Date:    12/6/01
! Last revision:

subroutine drift2(nbod, mass, xj, yj, zj, vxj, vyj, vzj, dt)

    use swift_memory_support, only: swift_status
    include 'swift.inc'

    !...  Inputs Only:
    integer nbod
    real(8) mass(nbod), dt

    !...  Inputs and Outputs:
    real(8) xj(nbod), yj(nbod), zj(nbod)
    real(8) vxj(nbod), vyj(nbod), vzj(nbod)

    !...  Internals:
    real(8) etaj
    integer j, iflg

    !----
    !...  Executable code

    ! Take a drift forward dt

    etaj = mass(1)
    do j = 2, nbod
        etaj = etaj + mass(j)
        call drift_one(etaj, xj(j), yj(j), zj(j), &
                vxj(j), vyj(j), vzj(j), dt, iflg)
    if(swift_status /= 0) return
        if(iflg/=0) then
            swift_status=3
            return
        endif
    enddo

    return
end
!--------------------------------------------------------------------------

!*************************************************************************
!                        GETACCJ.F
!*************************************************************************
! This subroutine calculates the acceleration on the massive particles
! in the JACOBI frame.
!             Input:
!                 nbod        ==>  number of massive bodies (int scalor)
!                 mass        ==>  mass of bodies (real array)
!                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
!                 xj,yj,zj    ==>  position in jacobi coord (real arrays)
!             Output:
!                 axj,ayj,azj ==>  acceleration in Jacobi coord (real arrays)
!
! Remarks: Based on getacch
! Author:  Man Hoi Lee
! Date:    12/6/01
! Last revision:

subroutine getaccj(nbod, mass, j2rp2, j4rp4, xh, yh, zh, &
        xj, yj, zj, axj, ayj, azj)

    use swift_memory_support, only: swift_status
    include 'swift.inc'

    !...  Inputs:
    integer nbod
    real(8) mass(NPLMAX), xh(NPLMAX), yh(NPLMAX), zh(NPLMAX), j2rp2, j4rp4
    real(8) xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)

    !...  Outputs:
    real(8) axj(NPLMAX), ayj(NPLMAX), azj(NPLMAX)

    !...  Internals:
    integer i
    real(8) ir3h(NPLMAX), ir3j(NPLMAX)
    real(8) irh(NPLMAX), irj(NPLMAX)
    real(8) axj1(NPLMAX), ayj1(NPLMAX), azj1(NPLMAX)
    real(8) axj2(NPLMAX), ayj2(NPLMAX), azj2(NPLMAX)
    real(8) axj3(NPLMAX), ayj3(NPLMAX), azj3(NPLMAX)
    real(8) aoblx(NPLMAX), aobly(NPLMAX), aoblz(NPLMAX)

    !----
    !...  Executable code

    !...  get thr r^-3's
    call getacch_ir3(nbod, 2, xh, yh, zh, ir3h, irh)
    if(swift_status /= 0) return
    call getacch_ir3(nbod, 2, xj, yj, zj, ir3j, irj)
    if(swift_status /= 0) return

    !...  now the first terms
    call getaccj_aj1(nbod, mass, xh, yh, zh, xj, yj, zj, ir3h, ir3j, &
            axj1, ayj1, azj1)
    if(swift_status /= 0) return

    !...  now the second terms
    call getaccj_aj2(nbod, mass, xh, yh, zh, ir3h, axj2, ayj2, azj2)
    if(swift_status /= 0) return

    !...  now the third terms
    call getaccj_aj3(nbod, mass, xh, yh, zh, axj3, ayj3, azj3)
    if(swift_status /= 0) return

    !...  add them all together
    axj(1) = 0.0
    ayj(1) = 0.0
    azj(1) = 0.0
    do i = 2, nbod
        axj(i) = axj1(i) + axj2(i) + axj3(i)
        ayj(i) = ayj1(i) + ayj2(i) + ayj3(i)
        azj(i) = azj1(i) + azj2(i) + azj3(i)
    enddo

    !...  Now do j2 and j4 stuff
    if(j2rp2/=0.0d0) then
        call obl_acc(nbod, mass, j2rp2, j4rp4, xj, yj, zj, irh, &
                aoblx, aobly, aoblz)
    if(swift_status /= 0) return
        do i = 2, nbod
            axj(i) = axj(i) + aoblx(i) - aoblx(1)
            ayj(i) = ayj(i) + aobly(i) - aobly(1)
            azj(i) = azj(i) + aoblz(i) - aoblz(1)
        enddo
    endif

    return
end
! getaccj

!---------------------------------------------------------------------

!*************************************************************************
!                        GETACCJ_AJ1.F
!*************************************************************************
! This subroutine calculates the 1st term of acceleration
! on the massive particles in the JACOBI frame.
!             Input:
!                 nbod        ==>  number of massive bodies (int scalar)
!                 mass        ==>  mass of bodies (real array)
!                 xh,yh,zh    ==>  position in heliocentric coord (real array)
!                 xj,yj,zj    ==>  position in jacobi coord (real array)
!                 ir3h        ==> inv radii in heliocentric coord (real array)
!                 ir3j        ==> inv radii in jacobi coord (real array)
!             Output:
!                 axj1,ayj1,azj1 ==>  1st term acceleration in Jacobi coord
!                                    (real array)
!
! Author:  Man Hoi Lee
! Date:    12/6/01
! Last revision:

subroutine getaccj_aj1(nbod, mass, xh, yh, zh, xj, yj, zj, ir3h, ir3j, &
        axj1, ayj1, azj1)

    use swift_memory_support, only: swift_status
    include 'swift.inc'

    !...  Inputs:
    integer nbod
    real(8) mass(nbod), ir3h(nbod), ir3j(nbod)
    real(8) xj(nbod), yj(nbod), zj(nbod)
    real(8) xh(nbod), yh(nbod), zh(nbod)

    !...  Outputs:
    real(8) axj1(nbod), ayj1(nbod), azj1(nbod)


    !...  Internals:
    integer i
    real(8) etaim1, etai, aj1h, aj1j

    !----
    !...  Executable code

    axj1(1) = 0.0
    ayj1(1) = 0.0
    azj1(1) = 0.0

    axj1(2) = 0.0     ! because xj=xh
    ayj1(2) = 0.0
    azj1(2) = 0.0

    etaim1 = mass(1) + mass(2)

    do i = 3, nbod

        etai = etaim1 + mass(i)

        aj1j = xj(i) * ir3j(i)
        aj1h = mass(1) * xh(i) * ir3h(i) / etaim1
        axj1(i) = etai * (aj1j - aj1h)

        aj1j = yj(i) * ir3j(i)
        aj1h = mass(1) * yh(i) * ir3h(i) / etaim1
        ayj1(i) = etai * (aj1j - aj1h)

        aj1j = zj(i) * ir3j(i)
        aj1h = mass(1) * zh(i) * ir3h(i) / etaim1
        azj1(i) = etai * (aj1j - aj1h)

        etaim1 = etai

    enddo

    return
end
! getaccj_aj1

!---------------------------------------------------------------------

!*************************************************************************
!                        GETACCJ_AJ2.F
!*************************************************************************
! This subroutine calculates the 2nd term of acceleration
! on the massive particles in the JACOBI frame.
!             Input:
!                 nbod        ==>  number of massive bodies (int scalar)
!                 mass        ==>  mass of bodies (real array)
!                 xh,yh,zh    ==>  position in helio coord (real array)
!                 ir3h        ==> inv radii in helio coord (real array)
!             Output:
!                 axj2,ayj2,azj2 ==>  2nd term acceleration in Jacobi coord
!                                    (real array)
!
! Author:  Man Hoi Lee
! Date:    12/6/01
! Last revision:

subroutine getaccj_aj2(nbod, mass, xh, yh, zh, ir3h, axj2, ayj2, azj2)

    use swift_memory_support, only: swift_status
    include 'swift.inc'

    !...  Inputs:
    integer nbod
    real(8) mass(nbod), ir3h(nbod)
    real(8) xh(nbod), yh(nbod), zh(nbod)

    !...  Outputs:
    real(8) axj2(nbod), ayj2(nbod), azj2(nbod)

    !...  Internals:
    integer i
    real(8) etaim1, fac

    !----
    !...  Executable code

    axj2(nbod) = 0.0
    ayj2(nbod) = 0.0
    azj2(nbod) = 0.0

    do i = nbod, 3, -1
        fac = mass(i) * ir3h(i)
        axj2(i - 1) = axj2(i) - fac * xh(i)
        ayj2(i - 1) = ayj2(i) - fac * yh(i)
        azj2(i - 1) = azj2(i) - fac * zh(i)
    enddo

    axj2(1) = 0.0
    ayj2(1) = 0.0
    azj2(1) = 0.0

    etaim1 = 0.0d0
    do i = 2, nbod - 1
        etaim1 = etaim1 + mass(i - 1)
        axj2(i) = axj2(i) * mass(1) / etaim1
        ayj2(i) = ayj2(i) * mass(1) / etaim1
        azj2(i) = azj2(i) * mass(1) / etaim1
    enddo

    return
end
! getaccj_aj2
!---------------------------------------------------------------------

!*************************************************************************
!                        GETACCJ_AJ3.F
!*************************************************************************
! This subroutine calculates the 3rd term acceleration on the massive particles
! in the JACOBI frame. This term is the direct cross terms
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh,yh,zh      ==>  position in heliocentric coord
!                                   (real arrays)
!             Output:
!                 axj3,ayj3,azj3 ==>  3rd term of acceleration in Jacobi coord
!                                     (real arrays)
!
! Author:  Man Hoi Lee
! Date:    12/6/01
! Last revision:

subroutine getaccj_aj3(nbod, mass, xh, yh, zh, axj3, ayj3, azj3)

    use swift_memory_support, only: swift_status
    include 'swift.inc'

    !...  Inputs:
    integer nbod
    real(8) mass(nbod), xh(nbod), yh(nbod), zh(nbod)

    !...  Outputs:
    real(8) axj3(nbod), ayj3(nbod), azj3(nbod)

    !...  Internals:
    integer i, j
    real(8) axh3(NPLMAX), ayh3(NPLMAX), azh3(NPLMAX)
    real(8) dx, dy, dz, rji2, faci, facj, irij3
    real(8) etaim1, sumax, sumay, sumaz

    !------
    !...  Executable code

    do i = 1, nbod
        axh3(i) = 0.0
        ayh3(i) = 0.0
        azh3(i) = 0.0
    enddo

    do i = 2, nbod - 1
        do j = i + 1, nbod

            dx = xh(j) - xh(i)
            dy = yh(j) - yh(i)
            dz = zh(j) - zh(i)
            rji2 = dx * dx + dy * dy + dz * dz

            irij3 = 1.0d0 / (rji2 * sqrt(rji2))
            faci = mass(i) * irij3
            facj = mass(j) * irij3

            axh3(j) = axh3(j) - faci * dx
            ayh3(j) = ayh3(j) - faci * dy
            azh3(j) = azh3(j) - faci * dz

            axh3(i) = axh3(i) + facj * dx
            ayh3(i) = ayh3(i) + facj * dy
            azh3(i) = azh3(i) + facj * dz

        enddo
    enddo

    axj3(1) = 0.d0
    ayj3(1) = 0.d0
    azj3(1) = 0.d0

    axj3(2) = axh3(2)
    ayj3(2) = ayh3(2)
    azj3(2) = azh3(2)

    etaim1 = mass(1)

    sumax = 0.d0
    sumay = 0.d0
    sumaz = 0.d0

    do i = 3, nbod
        etaim1 = etaim1 + mass(i - 1)

        sumax = sumax + mass(i - 1) * axh3(i - 1)
        sumay = sumay + mass(i - 1) * ayh3(i - 1)
        sumaz = sumaz + mass(i - 1) * azh3(i - 1)

        axj3(i) = axh3(i) - sumax / etaim1
        ayj3(i) = ayh3(i) - sumay / etaim1
        azj3(i) = azh3(i) - sumaz / etaim1

    enddo

    return
end
! getaccj_aj3
!----------------------------------------------------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gr_corr(a, e, gmi, corr, dt)

    implicit none

    real(8)  a, e, T
    real(8) corr, THIRD, gmi, dt, prec_frac, a_corr
    real(8) PI, c, GMSUN, AU, st_mass
    parameter (PI = 3.14159265358979d0)
    parameter (c = 299792458.0d0)
    parameter (GMSUN = 1.32712497d20, AU = 1.49597892d13)
    parameter (THIRD = 1.d0 / 3.d0)

    st_mass = gmi / ((4.d0 * PI * PI) / (365.25 * 365.25))

    a_corr = a * 149597870700.0d0

    T = 2.0d0 * PI * sqrt((a_corr**3.0d0) / (GMSUN * (st_mass)))

    prec_frac = (T / 86400.0d0) / dt

    corr = (24.0d0 * (PI**3.0d0) * (a_corr**2.0d0)) / ((T**2.0d0) * &
            (c**2.0d0) * (1.0d0 - e**2.0d0))

    corr = corr / prec_frac

    return
end








