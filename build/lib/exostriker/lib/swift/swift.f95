! GENINIT_J3 reads Jacobi orbital elements of nbod planets and generates
! initial position and velocity in Jacobi coords.
! This version outputs rpl and rhill.

! Last modified by Man Hoi Lee, Aug 16, 2003.
subroutine geninit_j3_in_days(iuflg, mstar0, nbod, frho3, planets_params, planets_results, coordinates)
    include 'swift.inc'

    real(8) SMASSYR, MSUN, AU
    parameter (SMASSYR = TWOPI * TWOPI)
    parameter (MSUN = 1.989d33, AU = 1.49597892d13)

    integer nbod
    real(8) mass(NTPMAX), j2rp2, j4rp4
    real(8) xj(NTPMAX), yj(NTPMAX), zj(NTPMAX)
    real(8) vxj(NTPMAX), vyj(NTPMAX), vzj(NTPMAX)
    real(8) rpl(NTPMAX), rhill(NTPMAX)
    integer iflgchk
    logical(8) lclose

    real(8) mstar0, mpl0, frho3

    real(8) gm, a, ecc, inc, capom, omega, capm

    integer ialpha, i, iuflg, j

!    character(80) outfile

    real(8) planets_params(nbod, 7), planets_results(nbod, 9), coordinates(nbod*3+3, 3)
!f2py intent(in) iuflg, mstar0, nbod, frho3, outfile
!f2py intent(in) planets_params
!f2py intent(out) planets_results, coordinates
!f2py depend(nbod) planets_params, planets_results, coordinates

!100  write (*, *) ' Units Menu:'
!    write (*, *) '       0 ==> Solar masses, and AU'
!    write (*, *) '       1 ==> AU, and Years'
!    read (*, *) iuflg
!    if ((iuflg.ne.0).and.(iuflg.ne.1)) goto 100
!    write (*, *) ' Mass of central star (in solar mass):'
!    read (*, *) mstar0
!    write (*, *) ' Number of planets:'
!    read (*, *) nbod
!    write (*, *) ' f/rho^(1/3) for radius (rho in g/cm^3):'
!    read (*, *) frho3
!    write (*, *) ' Enter name of data file for symba5_j:'
!    read (*, 201) outfile
!201  format (a)

    if ((iuflg/=0).and.(iuflg/=1)) then
        write(*,*) "Units Menu 0 ==> Solar masses, and AU'; and 1 ==> AU, and Years'; value received: ", iuflg
        return
    endif

    lclose = .true.
    iflgchk = 0
    nbod = nbod + 1
    if (iuflg==0) then
        mass(1) = mstar0
    else
        mass(1) = mstar0 * SMASSYR / (365.25 * 365.25)
    endif
    j2rp2 = 0.d0
    j4rp4 = 0.d0
    xj(1) = 0.d0
    yj(1) = 0.d0
    zj(1) = 0.d0
    vxj(1) = 0.d0
    vyj(1) = 0.d0
    vzj(1) = 0.d0

!    open (unit = 12, file = 'geninit_j.out')
!    write (12, *) ' m     a     e     i     omega capom M     rpl   rh'

    gm = mass(1)

    do i = 2, nbod
!        write (*, *) ' Mass of planet', i - 1, ' (in solar mass):'
!        read (*, *) mpl0
        mpl0 = planets_params(i-1, 1)
        if (iuflg==0) then
            mass(i) = mpl0
        else
            mass(i) = mpl0 * SMASSYR / (365.25 * 365.25)
        endif
        rpl(i) = frho3 * (1.5d0 * mpl0 * MSUN / TWOPI)**0.3333333333d0 / AU

        ialpha = -1
        gm = gm + mass(i)

!        write (*, *)&
!                ' a (AU), e, i (deg), omega (deg), capom (deg), M (deg):'
!        read (*, *) a, ecc, inc, omega, capom, capm
        a = planets_params(i-1, 2)
        ecc = planets_params(i-1, 3)
        inc = planets_params(i-1, 4)
        omega = planets_params(i-1, 5)
        capom = planets_params(i-1, 6)
        capm = planets_params(i-1, 7)

        inc = inc / DEGRAD
        capom = capom / DEGRAD
        omega = omega / DEGRAD
        capm = capm / DEGRAD

        rhill(i) = a * (mass(i) / (3.d0 * mass(1)))**0.3333333333d0

!        write (12, 251) mass(i), a, ecc, inc, omega, capom, capm, &
!                rpl(i), rhill(i)
!251     format (1p9e12.4)
        planets_results(i-1, :) = (/mass(i), a, ecc, inc, omega, capom, capm, &
                rpl(i), rhill(i)/)

        call ORBEL_EL2XV (gm, ialpha, a, ecc, inc, capom, omega, capm, &
                xj(i), yj(i), zj(i), vxj(i), vyj(i), vzj(i))
    enddo

!    close (unit = 12)

!    call IO_DUMP_PL_SYMBA (outfile, nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
!            lclose, iflgchk, rpl, rhill, j2rp2, j4rp4)

    if(btest(iflgchk,5))  then ! bit 5 is set
           coordinates(1, :) = (/ mass(1),j2rp2,j4rp4 /)
        else
            coordinates(1, :) = (/ mass(1), 0.d0, 0.d0 /)
        endif
        coordinates(2, :) = (/ xj(1),yj(1),zj(1) /)
        coordinates(3, :) = (/ vxj(1),vyj(1),vzj(1) /)

    do j=2,nbod
           if(lclose) then
              coordinates((j-1)*3+1, :) = (/ mass(j),rhill(j),rpl(j) /)
           else
              coordinates((j-1)*3+1, :) = (/ mass(j),rhill(j), 0.d0 /)
           endif
           coordinates((j-1)*3+2, :) = (/ xj(j),yj(j),zj(j) /)
           coordinates((j-1)*3+3, :) = (/ vxj(j),vyj(j),vzj(j) /)
    enddo

    return
end

!**********************************************************************
!		      SWIFT_SYMBA5_J.F
!**********************************************************************
!
!                 To run, need 2 input files. The code prompts for
!                 the file names, but examples are :
!
!                   parameter file like       param.in
!		    planet file like          pl.in
!
!  NOTE:  No test particles in this code and the massive bodies
!         are dimensioned at NTPMAX
!
!  This version inputs/outputs Jacobi coords and orbital elements.
!
! Authors: Man Hoi Lee
! Date:    1/20/05
! Last revision:
subroutine swift_symba5_j(t0, tstop, dt, dtout, dtdump, lflg, &
        rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat, &
        nbod, coordinates, mtiny, energy_out_size, energy_out)
    include 'swift.inc'
    include 'io.inc'

    real(8) mass(NTPMAX), j2rp2, j4rp4
    real(8)  xj(NTPMAX), yj(NTPMAX), zj(NTPMAX)
    real(8)  vxj(NTPMAX), vyj(NTPMAX), vzj(NTPMAX)

    real(8)  xh(NTPMAX), yh(NTPMAX), zh(NTPMAX)
    real(8)  vxh(NTPMAX), vyh(NTPMAX), vzh(NTPMAX)

    real(8)  xjt(1), yjt(1), zjt(1)       ! Dummy for the io
    real(8)  vxjt(1), vyjt(1), vzjt(1)
    integer ntp, istat(1, NSTAT)

    integer nbod, i1st, i, nbodm, nbodo
    integer iflgchk, iub, iuj, iud, iue, ium

    real(8)  t0, tstop, dt, dtout, dtdump
    real(8)  t, tout, tdump, tfrac, eoff
    real(8)  rpl(NTPMAX), rhill(NTPMAX)

    real(8)  rmin, rmax, rmaxu, qmin, mtiny
    logical(2)  lclose
    integer isenc, ihills
    integer mergelst(2, NTPMAX), mergecnt
    integer(2)  iecnt(NTPMAX)

    character(80)  outfile, fopenstat

    logical(1) lflg(0:IO_NBITS-1)
    integer energy_out_size, iter_fill, j
    real(8) coordinates(nbod*3+3, 3), ke, pot, energy, eltot(3), energy_out(energy_out_size, 5)
!f2py intent(in) t0, tstop, dt, dtout, dtdump, lflg, nbod, mtiny
!f2py intent(in) rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat
!f2py intent(in) energy_out_size
!f2py intent(in) coordinates
!f2py intent(out) energy_out
!f2py depend(nbod) coordinates
!f2py depend(energy_out_size) energy_out

    !-----
    !...  Executable code

    ntp = 0

    !...  print version number
!    call util_version

    ! Get data for the run and the test particles
!    write(*, *) 'Enter name of parameter data file : '
!    read(*, 999) inparfile
!    call io_init_param(inparfile, t0, tstop, dt, dtout, dtdump, &
!            iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat)

    iflgchk=0
    do i=0,IO_NBITS-1
       if(lflg(i)) then
          iflgchk = ibset(iflgchk,i)
       endif
    enddo

    if(btest(iflgchk,0) .and. btest(iflgchk,1))  then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid logical flags '
       write(*,*) '    You cannot request that both a real and ',&
                  '       an integer binary file be written '
       call util_exit(1)
    endif

    if(btest(iflgchk,4))  then ! bit 4 is set
    else
       rmin = -1.0
       rmax = -1.0
       rmaxu = -1.0
       qmin = -1.0
       lclose = .false.
    endif

    if((fopenstat(1:3)/='new') .and. &
       (fopenstat(1:3)/='NEW') .and. &
       (fopenstat(1:7)/='unknown') .and. &
       (fopenstat(1:7)/='UNKNOWN') .and. &
       (fopenstat(1:6)/='append') .and. &
       (fopenstat(1:6)/='APPEND') ) then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid status flag:',fopenstat,':'
       call util_exit(1)
    endif

    ! Prompt and read name of planet data file
!    write(*, *) ' '
!    write(*, *) 'Enter name of planet data file : '
!    read(*, 999) inplfile
!    999  format(a)
!    call io_init_pl_symba(inplfile, lclose, iflgchk, nbod, mass, &
!            xj, yj, zj, vxj, vyj, vzj, rpl, rhill, j2rp2, j4rp4)

    if(nbod>NTPMAX) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   The number of massive bodies,',nbod,','
         write(*,*) '   is too large, it must be less than',NTPMAX
         call util_exit(1)
    endif

    do j=1,nbod
       mass(j) =  coordinates((j-1)*3+1, 1)
       if(j==1) then
           j2rp2 = coordinates(1, 2)
           j4rp4 = coordinates(1, 3)
           rpl(1) = 0.0d0
           rhill(1) = 0.0d0
       else
           rhill(j) =  coordinates((j-1)*3+1, 2)
           rpl(j) =  coordinates((j-1)*3+1, 3)
       endif

       xj(j) = coordinates((j-1)*3+2, 1)
       yj(j) = coordinates((j-1)*3+2, 2)
       zj(j) = coordinates((j-1)*3+2, 3)
       vxj(j) = coordinates((j-1)*3+3, 1)
       vyj(j) = coordinates((j-1)*3+3, 2)
       vzj(j) = coordinates((j-1)*3+3, 3)
    enddo

    if((xj(1)/=0.0d0) .or. &
       (yj(1)/=0.0d0) .or. &
       (zj(1)/=0.0d0) .or. &
       (vxj(1)/=0.0d0) .or. &
       (vyj(1)/=0.0d0) .or. &
       (vzj(1)/=0.0d0) ) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   Input MUST be in heliocentric coordinates '
         write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
         call util_exit(1)
    endif

!    write(*, *) 'Enter the smallest mass to self gravitate :'
!    read(*, *) mtiny
!    write(*, *) ' mtiny = ', mtiny

    ! Initialize initial time and times for first output and first dump
    t = t0
    tout = t0 + dtout
    tdump = t0 + dtdump

    iub = 20
    iuj = 30
    iud = 40
    iue = 60
    ium = 21

    !...    Do the initial io write
    if(btest(iflgchk, 1))  then ! bit 1 is set
        call io_write_frame_j(t0, nbod, ntp, mass, xj, yj, zj, vxj, vyj, vzj, &
                xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, outfile, iub, fopenstat)
        call io_write_mass(t0, nbod, mass, outfile, ium, fopenstat)
    endif
    if(btest(iflgchk, 2))  then ! bit 2 is set
        eoff = 0.0d0
        call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, xh, yh, zh, &
                vxh, vyh, vzh)
!        call anal_energy_write(t0, nbod, mass, j2rp2, j4rp4, xh, yh, zh, vxh, &
!                vyh, vzh, iue, fopenstat, eoff)
        call anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,&
                                vxh,vyh,vzh,ke,pot,energy,eltot)
        energy = energy + eoff
        energy_out(1, :) = (/ t0, energy, eltot /)
        if(i1st==0) then
            i1st=1
        endif
    endif

    !...  must initize discard io routine
    if(btest(iflgchk, 4))  then ! bit 4 is set
        call io_discard_mass(0, t, 0, mass(1), rpl(1), xj(1), yj(1), zj(1), &
                vxj(1), vyj(1), vzj(1), iud, -1, fopenstat)
    endif

    !...  Calculate the location of the last massive particle
    call symba5_nbodm(nbod, mass, mtiny, nbodm)

    ihills = 0
    i1st = 0
    iter_fill = 0
    !***************here's the big loop *************************************
    write(*, *) ' ************** MAIN LOOP ****************** '

    do while ((t <= tstop) .and. (nbod>1))
        i1st = 0

        call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, xh, yh, zh, &
                vxh, vyh, vzh)

        call symba5_step_pl(i1st, t, nbod, nbodm, mass, j2rp2, j4rp4, &
                xh, yh, zh, vxh, vyh, vzh, dt, lclose, rpl, isenc, &
                mergelst, mergecnt, iecnt, eoff, rhill, mtiny)

        call coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, &
                vxj, vyj, vzj)

        t = t + dt

        if(btest(iflgchk, 4))  then ! bit 4 is set
            nbodo = nbod
            call discard_massive5(t, dt, nbod, mass, xj, yj, zj, &
                    vxj, vyj, vzj, rmin, rmax, rmaxu, qmin, lclose, &
                    rpl, rhill, isenc, mergelst, mergecnt, &
                    iecnt, eoff, i1st)
            if(nbodo/=nbod) then
                call symba5_nbodm(nbod, mass, mtiny, nbodm)
            endif
        endif

        ! if it is time, output orb. elements,
        if(t >= tout) then

            if(btest(iflgchk, 1))  then ! bit 1 is set
                call  io_write_frame_j(t, nbod, ntp, mass, xj, yj, zj, vxj, &
                        vyj, vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, outfile, &
                        iub, fopenstat)
                call io_write_mass(t, nbod, mass, outfile, ium, fopenstat)
            endif

            tout = tout + dtout
        endif

        ! If it is time, do a dump
        if(t>=tdump) then

            tfrac = (t - t0) / (tstop - t0)
            write(*, 998) t, tfrac, nbod
998         format(' Time = ', 1p1e12.5, ': fraction done = ', 0pf5.3, &
                    ': Number of bodies =', i4)
!            call io_dump_pl_symba('dump_pl.dat', nbod, mass, xj, yj, zj, &
!                    vxj, vyj, vzj, lclose, iflgchk, rpl, rhill, j2rp2, j4rp4)
!            call io_dump_param('dump_param.dat', t, tstop, dt, dtout, &
!                    dtdump, iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile)
            tdump = tdump + dtdump

            if(btest(iflgchk, 2))  then ! bit 2 is set
                call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, xh, yh, zh, &
                        vxh, vyh, vzh)
!                call anal_energy_write(t, nbod, mass, j2rp2, j4rp4, &
!                        xh, yh, zh, vxh, vyh, vzh, iue, fopenstat, eoff)
                call anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,&
                                vxh,vyh,vzh,ke,pot,energy,eltot)
                energy = energy + eoff
                iter_fill = iter_fill + 1
                energy_out(iter_fill, :) = (/ t0, energy, eltot /)
                if(i1st==0) then
                    i1st=1
                endif
            endif
        endif
    enddo
    !********** end of the big loop from time 't0' to time 'tstop'

    ! Do a final dump for possible resumption later

!    call io_dump_pl_symba('dump_pl.dat', nbod, mass, xj, yj, zj, &
!            vxj, vyj, vzj, lclose, iflgchk, rpl, rhill, j2rp2, j4rp4)
!    call io_dump_param('dump_param.dat', t, tstop, dt, dtout, &
!            dtdump, iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile)

!    call util_exit(0)
     return
end

!**********************************************************************
!		      SWIFT_MVS_J.F
!**********************************************************************
!
!                 NO CLOSE ENCOUNTERS
!                 To run, need 2 input files. The code prompts for
!                 the file names, but examples are :
!
!                   parameter file like       param.in
!		    planet file like          pl.in
!
!  This version inputs/outputs Jacobi coords and orbital elements and
!  groups terms so that hierarchical systems with comparable masses can
!  be integrated.
!  NOTE:  No test particles in this code
!
! Author:  Man Hoi Lee
! Date:    12/6/01
! Last revision:

subroutine swift_mvs_j(t0, tstop, dt, dtout, dtdump, lflg, &
        rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat, &
        nbod, coordinates, energy_out_size, energy_out)
    include 'swift.inc'
    include 'io.inc'

    real(8) mass(NPLMAX), j2rp2, j4rp4
    real(8) xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
    real(8) vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)

    real(8) xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)

    real(8) xjt(1), yjt(1), zjt(1)       ! Dummy for the io
    real(8) vxjt(1), vyjt(1), vzjt(1)

    integer istat(1, NSTAT), i1st
    integer nbod, ntp, nleft
    integer iflgchk, iub, iuj, iud, iue
    real(8) rstat(1, NSTATR)

    real(8) t0, tstop, dt, dtout, dtdump
    real(8) t, tout, tdump, tfrac, eoff

    real(8) rmin, rmax, rmaxu, qmin, rplsq(NPLMAX)
    logical(2) lclose

    character(80) outfile, fopenstat

    logical(1) lflg(0:IO_NBITS-1)
    integer energy_out_size, iter_fill, j, i
    real(8) coordinates(nbod*3+3, 3), ke, pot, energy, eltot(3), energy_out(energy_out_size, 5)
!f2py intent(in) t0, tstop, dt, dtout, dtdump, lflg, nbod
!f2py intent(in) rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat
!f2py intent(in) energy_out_size
!f2py intent(in) coordinates
!f2py intent(out) energy_out
!f2py depend(nbod) coordinates
!f2py depend(energy_out_size) energy_out

    !-----
    !...    Executable code

    ntp = 0

!    !...    print version number
!    call util_version
!
!    ! Get data for the run and the test particles
!    write(*, *) 'Enter name of parameter data file : '
!    read(*, 999) inparfile
!    call io_init_param(inparfile, t0, tstop, dt, dtout, dtdump, &
!            iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat)

    iflgchk=0
    do i=0,IO_NBITS-1
       if(lflg(i)) then
          iflgchk = ibset(iflgchk,i)
       endif
    enddo

    if(btest(iflgchk,0) .and. btest(iflgchk,1))  then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid logical flags '
       write(*,*) '    You cannot request that both a real and ',&
                  '       an integer binary file be written '
       call util_exit(1)
    endif

    if(btest(iflgchk,4))  then ! bit 4 is set
    else
       rmin = -1.0
       rmax = -1.0
       rmaxu = -1.0
       qmin = -1.0
       lclose = .false.
    endif

    if((fopenstat(1:3)/='new') .and. &
       (fopenstat(1:3)/='NEW') .and. &
       (fopenstat(1:7)/='unknown') .and. &
       (fopenstat(1:7)/='UNKNOWN') .and. &
       (fopenstat(1:6)/='append') .and. &
       (fopenstat(1:6)/='APPEND') ) then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid status flag:',fopenstat,':'
       call util_exit(1)
    endif

!    ! Prompt and read name of planet data file
!    write(*, *) ' '
!    write(*, *) 'Enter name of planet data file : '
!    read(*, 999) inplfile
!    999    format(a)
!    call io_init_pl(inplfile, lclose, iflgchk, nbod, mass, xj, yj, zj, &
!            vxj, vyj, vzj, rplsq, j2rp2, j4rp4)

    if(nbod>NTPMAX) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   The number of massive bodies,',nbod,','
         write(*,*) '   is too large, it must be less than',NTPMAX
         call util_exit(1)
    endif

    do j=1,nbod
       mass(j) =  coordinates((j-1)*3+1, 1)
       if(j==1) then
           j2rp2 = coordinates(1, 2)
           j4rp4 = coordinates(1, 3)
           rplsq(1) = 0.0d0
       else
           rplsq(j) =  coordinates((j-1)*3+1, 3)**2
       endif

       xj(j) = coordinates((j-1)*3+2, 1)
       yj(j) = coordinates((j-1)*3+2, 2)
       zj(j) = coordinates((j-1)*3+2, 3)
       vxj(j) = coordinates((j-1)*3+3, 1)
       vyj(j) = coordinates((j-1)*3+3, 2)
       vzj(j) = coordinates((j-1)*3+3, 3)
    enddo

    if((xj(1)/=0.0d0) .or. &
       (yj(1)/=0.0d0) .or. &
       (zj(1)/=0.0d0) .or. &
       (vxj(1)/=0.0d0) .or. &
       (vyj(1)/=0.0d0) .or. &
       (vzj(1)/=0.0d0) ) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   Input MUST be in heliocentric coordinates '
         write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
         call util_exit(1)
    endif

    ! Initialize initial time and times for first output and first dump
    t = t0
    tout = t0 + dtout
    tdump = t0 + dtdump

    iub = 20
    iuj = 30
    iud = 40
    iue = 60

    !...    Do the initial io write
    if(btest(iflgchk, 1))  then ! bit 1 is set
        call io_write_frame_j(t0, nbod, ntp, mass, xj, yj, zj, vxj, vyj, vzj, &
                xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, outfile, iub, fopenstat)
    endif
    if(btest(iflgchk, 2))  then    ! bit 2 is set
        eoff = 0.0d0
        call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, xh, yh, zh, &
                vxh, vyh, vzh)
!        call anal_energy_write(t0, nbod, mass, j2rp2, j4rp4, xh, yh, zh, vxh, &
!                vyh, vzh, iue, fopenstat, eoff)
        call anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,&
                                vxh,vyh,vzh,ke,pot,energy,eltot)
        energy = energy + eoff
        energy_out(1, :) = (/ t0, energy, eltot /)
        if(i1st==0) then
            i1st=1
        endif
    endif

    !...    must initize discard io routine
    if(btest(iflgchk, 4))  then ! bit 4 is set
        call io_discard_write(0, t, nbod, ntp, xj, yj, zj, vxj, vyj, &
                vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, iud, &
                'discard.out', fopenstat, nleft)
    endif

    nleft = ntp
    i1st = 0
    iter_fill = 0
    !***************here's the big loop *************************************
    write(*, *) ' ************** MAIN LOOP ****************** '

    do while ((t <= tstop) .and.&
            ((ntp==0).or.(nleft>0)))

        call step_kdk_pl_j(i1st, nbod, mass, j2rp2, j4rp4, &
                xj, yj, zj, vxj, vyj, vzj, dt)

        t = t + dt

        if(btest(iflgchk, 4))  then    ! bit 4 is set
            call discard(t, dt, nbod, ntp, mass, xj, yj, zj, vxj, vyj, vzj, &
                    xjt, yjt, zjt, vxjt, vyjt, vzjt, rmin, rmax, rmaxu, &
                    qmin, lclose, rplsq, istat, rstat)
            call io_discard_write(1, t, nbod, ntp, xj, yj, zj, vxj, vyj, &
                    vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, iud, &
                    'discard.out', fopenstat, nleft)
        else
            nleft = ntp
        endif

        ! if it is time, output orb. elements,
        if(t >= tout) then

            if(btest(iflgchk, 1))  then    ! bit 1 is set
                call  io_write_frame_j(t, nbod, ntp, mass, xj, yj, zj, vxj, &
                        vyj, vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, outfile, &
                        iub, fopenstat)
            endif

            tout = tout + dtout
        endif

        ! If it is time, do a dump
        if(t>=tdump) then

            tfrac = (t - t0) / (tstop - t0)
            write(*, 9980) t, tfrac, nleft
9980        format(' Time = ', 1p1e12.5, ': fraction done = ', 0pf5.3, &
                    ': Number of active tp =', i4)
!            call io_dump_pl('dump_pl.dat', nbod, mass, xj, yj, zj, &
!                    vxj, vyj, vzj, lclose, iflgchk, rplsq, j2rp2, j4rp4)
!            call io_dump_param('dump_param.dat', t, tstop, dt, dtout, &
!                    dtdump, iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile)
            tdump = tdump + dtdump

            if(btest(iflgchk, 2))  then    ! bit 2 is set
                call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, xh, yh, zh, &
                        vxh, vyh, vzh)
!                call anal_energy_write(t, nbod, mass, j2rp2, j4rp4, &
!                        xh, yh, zh, vxh, vyh, vzh, iue, fopenstat, eoff)
                call anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,&
                                vxh,vyh,vzh,ke,pot,energy,eltot)
                energy = energy + eoff
                iter_fill = iter_fill + 1
                energy_out(iter_fill, :) = (/ t0, energy, eltot /)
                if(i1st==0) then
                    i1st=1
                endif
            endif
        endif
    enddo
    !********** end of the big loop from time 't0' to time 'tstop'

    ! Do a final dump for possible resumption later

!    call io_dump_pl('dump_pl.dat', nbod, mass, xj, yj, zj, &
!            vxj, vyj, vzj, lclose, iflgchk, rplsq, j2rp2, j4rp4)
!    call io_dump_param('dump_param.dat', t, tstop, dt, dtout, &
!            dtdump, iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile)

!    call util_exit(0)
    return
end
! swift_mvs.f
!---------------------------------------------------------------------

!**********************************************************************
!		      SWIFT_MVS_J_GR.F
!**********************************************************************
!
!                 NO CLOSE ENCOUNTERS
!                 To run, need 2 input files. The code prompts for
!                 the file names, but examples are :
!
!                   parameter file like       param.in
!		    planet file like          pl.in
!
!  This version inputs/outputs Jacobi coords and orbital elements and
!  groups terms so that hierarchical systems with comparable masses can
!  be integrated. GR precesion included
!  NOTE:  No test particles in this code
!
! Author:  Man Hoi Lee/Trifon Trifonov
! Date:    12/6/01
! Last revision: 01/06/2015

subroutine swift_mvs_j_gr(t0, tstop, dt, dtout, dtdump, lflg, &
        rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat, &
        nbod, coordinates, ll, energy_out_size, energy_out)
    include 'swift.inc'
    include 'io.inc'

    real(8) mass(NPLMAX), j2rp2, j4rp4, corr
    real(8) xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
    real(8) vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)

    real(8) xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)

    real(8) xjt(1), yjt(1), zjt(1)       ! Dummy for the io
    real(8) vxjt(1), vyjt(1), vzjt(1)

    integer istat(1, NSTAT), i1st, l, ll
    integer nbod, ntp, nleft, i, ialphai
    integer iflgchk, iub, iuj, iud, iue
    real(8) rstat(1, NSTATR)

    real(8) t0, tstop, dt, dtout, dtdump
    real(8) t, tout, tdump, tfrac, eoff

    real(8) rmin, rmax, rmaxu, qmin, rplsq(NPLMAX)
    logical(2) lclose

    character(80) outfile,fopenstat
    real(8) gmi, ai, ei, inci, capomi, omegai, capmi

    logical(1) lflg(0:IO_NBITS-1)
    integer energy_out_size, iter_fill, j
    real(8) coordinates(nbod*3+3, 3), ke, pot, energy, eltot(3), energy_out(energy_out_size, 5)
!f2py intent(in) t0, tstop, dt, dtout, dtdump, lflg, nbod, ll
!f2py intent(in) rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat
!f2py intent(in) energy_out_size
!f2py intent(in) coordinates
!f2py intent(out) energy_out
!f2py depend(nbod) coordinates
!f2py depend(energy_out_size) energy_out

    !-----
    !...    Executable code

    ntp = 0

!    !...    print version number
!    call util_version
!
!    ! Get data for the run and the test particles
!    write(*, *) 'Enter name of parameter data file : '
!    read(*, 999) inparfile
!    call io_init_param(inparfile, t0, tstop, dt, dtout, dtdump, &
!            iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat)

    iflgchk=0
    do i=0,IO_NBITS-1
       if(lflg(i)) then
          iflgchk = ibset(iflgchk,i)
       endif
    enddo

    if(btest(iflgchk,0) .and. btest(iflgchk,1))  then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid logical flags '
       write(*,*) '    You cannot request that both a real and ',&
                  '       an integer binary file be written '
       call util_exit(1)
    endif

    if(btest(iflgchk,4))  then ! bit 4 is set
    else
       rmin = -1.0
       rmax = -1.0
       rmaxu = -1.0
       qmin = -1.0
       lclose = .false.
    endif

    if((fopenstat(1:3)/='new') .and. &
       (fopenstat(1:3)/='NEW') .and. &
       (fopenstat(1:7)/='unknown') .and. &
       (fopenstat(1:7)/='UNKNOWN') .and. &
       (fopenstat(1:6)/='append') .and. &
       (fopenstat(1:6)/='APPEND') ) then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid status flag:',fopenstat,':'
       call util_exit(1)
    endif

!    ! Prompt and read name of planet data file
!    write(*, *) ' '
!    write(*, *) 'Enter name of planet data file : '
!    read(*, 999) inplfile
!    999    format(a)
!    call io_init_pl(inplfile, lclose, iflgchk, nbod, mass, xj, yj, zj, &
!            vxj, vyj, vzj, rplsq, j2rp2, j4rp4)

    if(nbod>NTPMAX) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   The number of massive bodies,',nbod,','
         write(*,*) '   is too large, it must be less than',NTPMAX
         call util_exit(1)
    endif

    do j=1,nbod
       mass(j) =  coordinates((j-1)*3+1, 1)
       if(j==1) then
           j2rp2 = coordinates(1, 2)
           j4rp4 = coordinates(1, 3)
           rplsq(1) = 0.0d0
       else
           rplsq(j) =  coordinates((j-1)*3+1, 3)**2
       endif

       xj(j) = coordinates((j-1)*3+2, 1)
       yj(j) = coordinates((j-1)*3+2, 2)
       zj(j) = coordinates((j-1)*3+2, 3)
       vxj(j) = coordinates((j-1)*3+3, 1)
       vyj(j) = coordinates((j-1)*3+3, 2)
       vzj(j) = coordinates((j-1)*3+3, 3)
    enddo

    if((xj(1)/=0.0d0) .or. &
       (yj(1)/=0.0d0) .or. &
       (zj(1)/=0.0d0) .or. &
       (vxj(1)/=0.0d0) .or. &
       (vyj(1)/=0.0d0) .or. &
       (vzj(1)/=0.0d0) ) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   Input MUST be in heliocentric coordinates '
         write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
         call util_exit(1)
    endif

!    write(*, *) 'On which step to apply GR precession (int) ?'
!    read(*, *) ll

    ! Initialize initial time and times for first output and first dump
    t = t0
    tout = t0 + dtout
    tdump = t0 + dtdump

    iub = 20
    iuj = 30
    iud = 40
    iue = 60

    !...    Do the initial io write
    if(btest(iflgchk, 1))  then ! bit 1 is set
        call io_write_frame_j(t0, nbod, ntp, mass, xj, yj, zj, vxj, vyj, vzj, &
                xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, outfile, iub, fopenstat)
    endif
    if(btest(iflgchk, 2))  then    ! bit 2 is set
        eoff = 0.0d0
        call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, xh, yh, zh, &
                vxh, vyh, vzh)
!        call anal_energy_write(t0, nbod, mass, j2rp2, j4rp4, xh, yh, zh, vxh, &
!                vyh, vzh, iue, fopenstat, eoff)
        call anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,&
                                vxh,vyh,vzh,ke,pot,energy,eltot)
        energy = energy + eoff
        energy_out(1, :) = (/ t0, energy, eltot /)
        if(i1st==0) then
            i1st=1
        endif
    endif

    !...    must initize discard io routine
    if(btest(iflgchk, 4))  then ! bit 4 is set
        call io_discard_write(0, t, nbod, ntp, xj, yj, zj, vxj, vyj, &
                vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, iud, &
                'discard.out', fopenstat, nleft)
    endif

    nleft = ntp
    i1st = 0
    l = 0
    iter_fill = 0
    !	ll = 1000
    !***************here's the big loop *************************************
    write(*, *) ' ************** MAIN LOOP ****************** '

    do while ((t <= tstop) .and.&
            ((ntp==0).or.(nleft>0)))

        if (l==ll) then
            gmi = mass(1)
            do i = 2, nbod
                gmi = gmi + mass(i)
                !            if ((dadt(i).ne.0.d0).or.(dedt(i).ne.0.d0)) then

                call orbel_xv2el(xj(i), yj(i), zj(i), vxj(i), vyj(i), &
                        vzj(i), gmi, ialphai, ai, ei, inci, capomi, omegai, capmi)

                call gr_corr(ai, ei, gmi, corr, dt)
                omegai = omegai + corr * ll
                !               write(*,*) ai

                call orbel_el2xv(gmi, ialphai, ai, ei, inci, capomi, omegai, &
                        capmi, xj(i), yj(i), zj(i), vxj(i), vyj(i), vzj(i))

            enddo
            l = 0
        endif

        l = l + 1
        !             write(*,*) l

        call step_kdk_pl_j(i1st, nbod, mass, j2rp2, j4rp4, &
                xj, yj, zj, vxj, vyj, vzj, dt)

        t = t + dt

        if(btest(iflgchk, 4))  then    ! bit 4 is set
            call discard(t, dt, nbod, ntp, mass, xj, yj, zj, vxj, vyj, vzj, &
                    xjt, yjt, zjt, vxjt, vyjt, vzjt, rmin, rmax, rmaxu, &
                    qmin, lclose, rplsq, istat, rstat)
            call io_discard_write(1, t, nbod, ntp, xj, yj, zj, vxj, vyj, &
                    vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, iud, &
                    'discard.out', fopenstat, nleft)
        else
            nleft = ntp
        endif

        ! if it is time, output orb. elements,
        if(t >= tout) then

            if(btest(iflgchk, 1))  then    ! bit 1 is set
                call  io_write_frame_j(t, nbod, ntp, mass, xj, yj, zj, vxj, &
                        vyj, vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, outfile, &
                        iub, fopenstat)
            endif

            tout = tout + dtout
        endif

        ! If it is time, do a dump
        if(t>=tdump) then

            tfrac = (t - t0) / (tstop - t0)
            write(*, 99800) t, tfrac, nleft
99800       format(' Time = ', 1p1e12.5, ': fraction done = ', 0pf5.3, &
                    ': Number of active tp =', i4)
!            call io_dump_pl('dump_pl.dat', nbod, mass, xj, yj, zj, &
!                    vxj, vyj, vzj, lclose, iflgchk, rplsq, j2rp2, j4rp4)
!            call io_dump_param('dump_param.dat', t, tstop, dt, dtout, &
!                    dtdump, iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile)
            tdump = tdump + dtdump

            if(btest(iflgchk, 2))  then    ! bit 2 is set
                call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, xh, yh, zh, &
                        vxh, vyh, vzh)
!                call anal_energy_write(t, nbod, mass, j2rp2, j4rp4, &
!                        xh, yh, zh, vxh, vyh, vzh, iue, fopenstat, eoff)
                call anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,&
                                vxh,vyh,vzh,ke,pot,energy,eltot)
                energy = energy + eoff
                iter_fill = iter_fill + 1
                energy_out(iter_fill, :) = (/ t0, energy, eltot /)
                if(i1st==0) then
                    i1st=1
                endif
            endif

        endif

    enddo
    !********** end of the big loop from time 't0' to time 'tstop'

    ! Do a final dump for possible resumption later

!    call io_dump_pl('dump_pl.dat', nbod, mass, xj, yj, zj, &
!            vxj, vyj, vzj, lclose, iflgchk, rplsq, j2rp2, j4rp4)
!    call io_dump_param('dump_param.dat', t, tstop, dt, dtout, &
!            dtdump, iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile)

!    call util_exit(0)
    return
end
! swift_mvs.f

subroutine follow_symba2(t0, tstop, dt, dtout, dtdump, lflg, &
                         rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat, &
                         nbod, coordinates, ifol, symba_out_size, symba_out)
    ! converts binary file to ascii file
    include 'swift.inc'
    include 'io.inc'

    real(8) mass(NTPMAX), dr
    real(8) xh(NTPMAX), yh(NTPMAX), zh(NTPMAX)
    real(8) vxh(NTPMAX), vyh(NTPMAX), vzh(NTPMAX)

    integer nbodm, nbod, nbod0, ierr, ifol, istep
    integer iflgchk, iu, i, id, nleft, ium
    integer io_read_hdr, io_read_line, io_read_mass
    integer io_read_hdr_r, io_read_line_r

    real(8) t0, tstop, dt, dtout, dtdump
    real(8) t, tmax

    real(8) rmin, rmax, rmaxu, qmin, rpl(NTPMAX), rhill(NTPMAX)
    logical(2) lclose
    real(8) a, e, inc, capom, omega, capm, j2rp2, j4rp4
    real(8) peri, apo, tg

    integer plist(NTPMAX), ifoln, j

    character(80) outfile, fopenstat

    logical(1) lflg(0:IO_NBITS-1)
    integer symba_out_size, iter_fill
    real(8) coordinates(nbod*3+3, 3), symba_out(symba_out_size, 11)
!f2py intent(in) t0, tstop, dt, dtout, dtdump, lflg, nbod, mtiny
!f2py intent(in) rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat
!f2py intent(in) symba_out_size
!f2py intent(in) coordinates
!f2py intent(out) symba_out
!f2py depend(nbod) coordinates
!f2py depend(symba_out_size) symba_out

    dtout = 0
    dtdump = 0
    dt = 0

    ! Get data for the run and the test particles
!    write(*, *) 'Enter name of parameter data file : '
!    read(*, 999) inparfile
!    call io_init_param(inparfile, t0, tstop, dt, dtout, dtdump, &
!            iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat)

    iflgchk=0
    do i=0,IO_NBITS-1
       if(lflg(i)) then
          iflgchk = ibset(iflgchk,i)
       endif
    enddo

    if(btest(iflgchk,0) .and. btest(iflgchk,1))  then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid logical flags '
       write(*,*) '    You cannot request that both a real and ',&
                  '       an integer binary file be written '
       call util_exit(1)
    endif

    if(btest(iflgchk,4))  then ! bit 4 is set
    else
       rmin = -1.0
       rmax = -1.0
       rmaxu = -1.0
       qmin = -1.0
       lclose = .false.
    endif

    if((fopenstat(1:3)/='new') .and. &
       (fopenstat(1:3)/='NEW') .and. &
       (fopenstat(1:7)/='unknown') .and. &
       (fopenstat(1:7)/='UNKNOWN') .and. &
       (fopenstat(1:6)/='append') .and. &
       (fopenstat(1:6)/='APPEND') ) then
       write(*,*) ' SWIFT ERROR: in io_init_param:'
       write(*,*) '    Invalid status flag:',fopenstat,':'
       call util_exit(1)
    endif

    ! Prompt and read name of planet data file
!    write(*, *) ' '
!    write(*, *) 'Enter name of planet data file : '
!    read(*, 999) inplfile
!    999    format(a)
!    call io_init_pl_symba(inplfile, lclose, iflgchk, nbod, mass, xh, yh, zh, &
!            vxh, vyh, vzh, rpl, rhill, j2rp2, j4rp4)

    if(nbod>NTPMAX) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   The number of massive bodies,',nbod,','
         write(*,*) '   is too large, it must be less than',NTPMAX
         call util_exit(1)
    endif

    do j=1,nbod
       mass(j) =  coordinates((j-1)*3+1, 1)
       if(j==1) then
           j2rp2 = coordinates(1, 2)
           j4rp4 = coordinates(1, 3)
           rpl(1) = 0.0d0
           rhill(1) = 0.0d0
       else
           rhill(j) =  coordinates((j-1)*3+1, 2)
           rpl(j) =  coordinates((j-1)*3+1, 3)
       endif

       xh(j) = coordinates((j-1)*3+2, 1)
       yh(j) = coordinates((j-1)*3+2, 2)
       zh(j) = coordinates((j-1)*3+2, 3)
       vxh(j) = coordinates((j-1)*3+3, 1)
       vyh(j) = coordinates((j-1)*3+3, 2)
       vzh(j) = coordinates((j-1)*3+3, 3)
    enddo

    if((xh(1)/=0.0d0) .or. &
       (yh(1)/=0.0d0) .or. &
       (zh(1)/=0.0d0) .or. &
       (vxh(1)/=0.0d0) .or. &
       (vyh(1)/=0.0d0) .or. &
       (vzh(1)/=0.0d0) ) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   Input MUST be in heliocentric coordinates '
         write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
         call util_exit(1)
    endif

    iu = 20
    ium = 30

    dr = 180.0 / PI

    if(btest(iflgchk, 0)) then
        write(*, *) ' Reading an integer*2 binary file '
    else if(btest(iflgchk, 1)) then
        write(*, *) ' Reading an real*4 binary file '
    else
        write(*, *) ' ERROR: no binary file format specified '
        write(*, *) '        in param file '
        stop
    endif

!    write(*, *) ' Input the particle number to follow '
!    read(*, *) ifol
    ifol = abs(ifol)
    write(*, *) ' Following particle ', ifol

    open(unit = iu, file = outfile, status = 'old', form = 'unformatted')
    open(unit = ium, file = 'mass.' // outfile, status = 'old', &
            form = 'unformatted')
    open(unit = 7, file = 'follow_symba.out')

    nbod0 = nbod
    do i = 1, nbod
        plist(i) = i
    enddo
    call follow_plist(0, ifol, ifoln, plist, nbod0, tg, tstop)

    ifoln = ifol

!    write(*, *) '1  2 3 4  5    6     7    8    9   10  11 '
!    write(*, *) 't,id,a,e,inc,capom,omega,capm,peri,apo, M '

    tmax = t0
    iter_fill = 0
1   continue
    if(btest(iflgchk, 0))  then ! bit 0 is set
        ierr = io_read_hdr(iu, t, nbod, nleft)
    else
        ierr = io_read_hdr_r(iu, t, nbod, nleft)
    endif
    if(ierr/=0) goto 2

    ierr = io_read_mass(t, nbodm, mass, ium)

    if(nbodm/=nbod) then
        write(*, *) ' Error 1:', nbod, nbodm
        stop
    endif

    do while(t>=tg)
        call follow_plist(1, ifol, ifoln, plist, nbod0, &
                tg, tstop)
    enddo

    istep = 0
    do i = 2, nbod
        if(btest(iflgchk, 0))  then ! bit 0 is set
            ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm)
        else
            ierr = io_read_line_r(iu, id, a, e, inc, capom, omega, capm)
        endif
        if(ierr/=0) goto 2

        if(abs(id)==ifoln) then
            istep = 1
            inc = inc * dr
            capom = capom * dr
            omega = omega * dr
            capm = capm * dr
            peri = a * (1.0d0 - e)
            apo = a * (1.0d0 + e)
!            write(7, 1000) t, ifoln, a, e, inc, capom, omega, capm, &
!                    peri, apo, mass(abs(id)) / mass(1)
!1000        format(1x, e15.7, 1x, i3, 1x, e13.5, 1x, f7.5, 4(1x, f9.4), &
!                    2(1x, e13.5), 1e13.5)
            iter_fill = iter_fill +1
            symba_out(iter_fill, :) = (/ real(t, 8), real(ifoln, 8), a, &
                                         e, inc, capom, omega, capm, &
                                         peri, apo, mass(abs(id)) / mass(1) /)
            tmax = t
        endif
    enddo

    if(istep==0) goto 2     ! did not find particle this times step
    goto 1
2   continue

    write(*, *) ' Tmax = ', tmax

    return
end
end

! swift_symba5_j_migrate6.f
!---------------------------------------------------------------------

!*************************************************************************
!                            IO_WRITE_FRAME_J
!*************************************************************************
! write out a whole frame to an real*4 binary file.
! both massive and test particles
!
!             Input:
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ntp            ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xj,yj,zj      ==>  current position in Jacobi coord
!                                    (real arrays)
!                 vxj,vyj,vzj   ==>  current velocity in Jacobi coord
!                                    (real arrays)
!                 xjt,yjt,zjt    ==>  current part position in Jacobi coord
!                                      (real arrays)
!                 vxjt,vyjt,vzjt ==>  current velocity in Jacobi coord
!                                        (real arrays)
!                 istat           ==>  status of the test paricles
!                                      (2d integer array)
!                                      istat(i,1) = 0 ==> active:  = 1 not
!                                      istat(i,2) = -1 ==> Danby did not work
!                 oname           ==> output file name (character string)
!                 iu              ==> unit number to write to
!                 fopenstat       ==>  The status flag for the open
!                                      statements of the output files.
!                                          (character(80) )
!
!
! Remarks: Based on io_write_frame
! Authors:  Hal Levison
! Date:    2/22/94
! Last revision:

subroutine io_write_frame_j(time, nbod, ntp, mass, xj, yj, zj, vxj, &
        vyj, vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, oname, &
        iu, fopenstat)

    include 'swift.inc'
    include 'io.inc'

    !...  Inputs:
    integer nbod, ntp, iu
    real(8)  mass(nbod), time
    integer istat(1, NSTAT)
    real(8)  xj(nbod), yj(nbod), zj(nbod)
    real(8)  vxj(nbod), vyj(nbod), vzj(nbod)
    real(8)  xjt(ntp), yjt(ntp), zjt(ntp)
    real(8)  vxjt(ntp), vyjt(ntp), vzjt(ntp)
    character(80)  oname, fopenstat

    !...  Internals
    integer i, id
    integer ialpha, ierr
    real(8)  a, e, inc, capom, omega, capm
    real(8)  gm
    integer i1st    ! =0 first time through; =1 after
    data i1st/0/
    save i1st

    !----
    !...  Executable code

    !...  if first time through open file
    if(i1st==0) then
        call io_open(iu, oname, fopenstat, 'UNFORMATTED', ierr)
        if(ierr/=0) then
            write(*, *) ' SWIFT ERROR: in io_write_frame: '
            write(*, *) '     Could not open binary output file:'
            call util_exit(1)
        endif
        i1st = 1
    else
        call io_open(iu, oname, 'append', 'UNFORMATTED', ierr)
    endif

    call io_write_hdr_r(iu, time, nbod, ntp, istat)

    !...  write out planets
    gm = mass(1)
    do i = 2, nbod
        gm = gm + mass(i)
        id = -1 * i
        call orbel_xv2el(xj(i), yj(i), zj(i), vxj(i), vyj(i), vzj(i), gm, &
                ialpha, a, e, inc, capom, omega, capm)
        call io_write_line_r(iu, id, a, e, inc, capom, omega, capm)
    enddo

    !...  write out test particles
    gm = mass(1)
    do i = 1, ntp
        if(istat(i, 1)==0) then
            call orbel_xv2el(xjt(i), yjt(i), zjt(i), vxjt(i), vyjt(i), &
                    vzjt(i), gm, ialpha, a, e, inc, capom, omega, capm)
            call io_write_line_r(iu, i, a, e, inc, capom, omega, capm)
        endif
    enddo

    close(iu)
    return
end
! io_write_frame_j
!----------------------------------------------------------------------

subroutine follow_plist(iflg, ifol, ifoln, plist, nbod, tg, tstop)

    include 'swift.inc'
    real(8) tg, tstop
    integer nbod
    integer iflg, ifol, ifoln, plist(nbod), iwhy
    integer ig, im, idum, i, ierr
    save iwhy

    if(iflg==0) then
        open(2, file = 'discard_mass.out', status = 'old', iostat = ierr)
        if(ierr/=0) then
            write(*, *) 'Could not open discard_mass.out'
            tg = 5.0 * tstop
            return            ! <====== NOTE
        endif
        read(2, *, iostat = ierr) tg, iwhy
        if(ierr/=0) then
            write(*, *) 'Could not read discard_mass.out'
            tg = 5.0 * tstop
            return            ! <====== NOTE
        endif
        ifoln = ifol
        return               ! <====== NOTE
    endif

    if(iwhy==2) then
        read(2, *) idum, im
        read(2, fmt = '(1x)')
        read(2, fmt = '(1x)')
        read(2, *) idum, ig
        call left_reorder(ig, im, nbod, plist)
        do i = 1, 5
            read(2, fmt = '(1x)')
        enddo
    else
        read(2, *) idum, ig
        im = -1
        call left_reorder(ig, im, nbod, plist)
        read(2, fmt = '(1x)')
        read(2, fmt = '(1x)')
    endif

    read(2, *, iostat = ierr) tg, iwhy
    if(ierr/=0) then
        tg = 5.0 * tstop
    endif

    ifoln = plist(ifol)

    return
end

!---------------------------------------------------------------------
subroutine left_reorder(ig, im, nbod, plist)

    include 'swift.inc'

    integer ig, nbod, plist(nbod), im, i

    do i = 1, nbod
        if(plist(i)==ig) then
            if(im>0) then
                plist(i) = im
            else
                plist(i) = -1
            endif
        endif
    enddo

    do i = 1, nbod
        if(plist(i)>ig) then
            plist(i) = plist(i) - 1
        endif
    enddo

    return
end

!*************************************************************************
!                            STEP_KDK_PL_J.F
!*************************************************************************
! This subroutine takes a step in Jacobi coord.
! Does a KICK than a DRIFT than a KICK.
! ONLY DOES MASSIVE PARTICLES
!
!             Input:
!                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xj,yj,zj      ==>  initial position in Jacobi coord
!                                    (real arrays)
!                 vxj,vyj,vzj   ==>  initial velocity in Jacobi coord
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 xj,yj,zj      ==>  final position in Jacobi coord
!                                       (real arrays)
!                 vxj,vyj,vzj   ==>  final velocity in Jacobi coord
!                                       (real arrays)
!
! Remarks: Adopted from step_kdk_pl
! Authors:  Man Hoi Lee
! Date:    12/6/01
! Last revision:

subroutine step_kdk_pl_j(i1st, nbod, mass, j2rp2, j4rp4, &
        xj, yj, zj, vxj, vyj, vzj, dt)

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
        !...     Get the accelerations in Jacobi frame. if frist time step
        call getaccj(nbod, mass, j2rp2, j4rp4, xh, yh, zh, &
                xj, yj, zj, axj, ayj, azj)
        i1st = 1    ! turn this off
    endif

    !...  Apply a Jacobi kick for a half dt
    call kickvh(nbod, vxj, vyj, vzj, axj, ayj, azj, dth)

    !..   Drift in Jacobi coords for the full step
    call drift2(nbod, mass, xj, yj, zj, vxj, vyj, vzj, dt)

    !...  After drift, compute helio. xh and vh for acceleration calculations
    call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
            xh, yh, zh, vxh, vyh, vzh)

    !...  Get the accelerations in Jacobi frame.
    call getaccj(nbod, mass, j2rp2, j4rp4, xh, yh, zh, xj, yj, zj, axj, ayj, azj)

    !...  Apply a Jacobi kick for a half dt
    call kickvh(nbod, vxj, vyj, vzj, axj, ayj, azj, dth)

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
        if(iflg/=0) then
            write(*, *) ' Planet ', j, ' is lost !!!!!!!!!'
            write(*, *) etaj, dt
            write(*, *) xj(j), yj(j), zj(j)
            write(*, *) vxj(j), vyj(j), vzj(j)
            write(*, *) ' STOPPING '
            call util_exit(1)
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
    call getacch_ir3(nbod, 2, xj, yj, zj, ir3j, irj)

    !...  now the first terms
    call getaccj_aj1(nbod, mass, xh, yh, zh, xj, yj, zj, ir3h, ir3j, &
            axj1, ayj1, azj1)

    !...  now the second terms
    call getaccj_aj2(nbod, mass, xh, yh, zh, ir3h, axj2, ayj2, azj2)

    !...  now the third terms
    call getaccj_aj3(nbod, mass, xh, yh, zh, axj3, ayj3, azj3)

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
