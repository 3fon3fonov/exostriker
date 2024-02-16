!cc   By Trifon Trifonov trifon@hku.hk
!cc   You can modify if as you want, but please
!cc   do not distribute without permision. 
!cc   This is not a final version!!! 
!cc   The final version will be available in the Python RVMod lib
!cc   Trifonov et al. (in prep).

subroutine kepfit_amoeba(epsil, deltat, amoebastarts, &
        when_to_kill, nt, &
        model_max, model_min, gr_flag_in, &
        st_mass, writeflag_best_par, writeflag_RV, &
        writeflag_fit, &
        ndset_in, ndata, data_array, files_param, &
        npl_in, array_npl, &
        final_params, &
        res_array, fit_return, &
        bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
        fit_array, dynamical_planets, coplar_inc)
    implicit none


    integer, intent(in) :: epsil, deltat, amoebastarts,npl_in
    integer, intent(in) :: writeflag_best_par, writeflag_RV
    integer, intent(in) :: writeflag_fit, gr_flag_in, ndset_in
    integer, intent(in) :: coplar_inc, ndata, nt
       
    real(4), intent(in) :: when_to_kill, model_max, model_min    

    real(8), intent(in) :: st_mass

    integer :: hkl, gr_flag
 
       
    integer :: npl, ndset, idset,  ma, mfit, i, j, NDSMAX, NPLMAX, MMAX      
    parameter (NDSMAX = 20, NPLMAX = 10, MMAX = 200)
    
    integer :: ii, iter, idsmax(NDSMAX), ia(MMAX) 
    
    integer, allocatable, dimension(:) :: ts
    real(8), allocatable, dimension(:) :: x, y, sig, y_in, ymod
    
!    real(8) :: y(20000), sig(20000), y_in(20000)
    real(8) :: a(MMAX), covar(MMAX, MMAX)
    real(8) :: rms, mass(NPLMAX), ap(NPLMAX)
    real(8) :: j_mass(NPLMAX), chisq
    real(8) :: x0, incl(NPLMAX), cap0m(NPLMAX)
    real(8) :: dt, t_max, loglik, dy, sig2i
    real(8) :: epoch, ftol, jitt(NDSMAX)
    real(8) :: dyda(MMAX), p(MMAX + 1, MMAX), yamoeba(MMAX + 1)
    real(8) :: ymod_pl(npl_in,ndata),ymod_pl2(npl_in,nt) 
    real(8) :: loglikk, ologlikk, dloglikk, best_w, best_we

!    character(80) version_input, version
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX)

    real(4) :: t_stop, t_init
    real(8) :: PI, twopi    
      
    integer, intent(in) :: dynamical_planets(npl_in)
    
    real(8), intent(in) :: data_array(ndata, 4)
    real(8), intent(in) :: files_param(ndset_in, 4)
    real(8), intent(in) :: array_npl(npl_in, 17, 2)
    real(8), intent(in) :: final_params(6)
    real(8), intent(out) :: res_array(ndata, 6+npl_in)
    real(8), intent(out) :: fit_return(4), fit_array(nt, 2+npl_in)
    real(8), intent(out) :: bestpar_1(npl_in, 17, 2), bestpar_2(ndset_in, 2)
    real(8), intent(out) :: bestpar_3(ndset_in, 2), bestpar_4(9 + 2 * npl_in)
    character(80) :: version
    character(20) :: mode
    

    external rvkep_kepamo, compute_abs_loglik_kep

!f2py intent(in) ndset_in, npl_in, ndata
!f2py intent(in) dynamical_planets
!f2py intent(out) res_array, fit_return, fit_array
!f2py intent(out) bestpar_1, bestpar_2, bestpar_3, bestpar_4
!f2py depend(ndset_in) files_array,files_param
!f2py depend(npl_in) array_npl, dynamical_planets
!f2py depend(ndata) data_array
!f2py depend(dt) fit_array

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    
    allocate(x(20000),y(20000),sig(20000),y_in(20000),ts(20000),ymod(20000))
!    allocate( t_stop,t_init)
    
    t_stop = 0
    t_init = 0
    gr_flag = gr_flag_in
    ndset = ndset_in
    npl = npl_in
    mode = "amoeba"
 
    rms = 0
    covar(:, :) = 0
!    dynamical_planets(:) = 0
!    coplar_inc = 0


    version = "1.05"

!    CALL getarg(1, version_input)
!    if(version_input=='-version') then
!        write(*, *) version
!        return
!    endif

    PI = 3.14159265358979d0
    twopi = 2.d0 * PI
    ftol = 0.000001d0
 

    call io_read_data(ndata, x, ts, y, sig, jitt, epoch, &
            x0, t_max, a, ia, ma, mfit, incl, cap0m, wdot, u_wdot, hkl, mode, &
            data_array, files_param, &
            array_npl, &
            final_params, coplar_inc)

    i = 0
    dloglikk = 10.d0
    loglikk = 0.d0
!    t_init = time()
    
    call cpu_time(t_init)    
    
    do while (dabs(dloglikk)>=0.000001d0)
        if (i.eq.amoebastarts) then
            i = 0
            exit
        endif

        i = i + 1
        ologlikk = loglikk

        call prepare_for_amoeba_kep(p, MMAX + 1, MMAX, yamoeba, a, ia, ma, &
                mfit, compute_abs_loglik_kep, ndata, x, y, dyda, ts, sig, hkl)
        call amoeba_kep(p, yamoeba, MMAX + 1, MMAX, mfit, ftol, &
                compute_abs_loglik_kep, &
                iter, ndata, x, y, dyda, ma, ts, sig, a, ia, loglikk, hkl)

        call cpu_time(t_stop)    
!        t_stop = time() - t_init

!        write(*,*) t_stop-t_init
        if ((t_stop-t_init)>=when_to_kill) then
            write(*, *) 'Max. time=', when_to_kill, 'sec ', &
                    'exceeded t_stop =', t_stop-t_init, 'sec '
            exit
        endif
        
!        CALL SECOND(t_stop)
!        if (t_stop.ge.when_to_kill) then
!           write(*,*) 'Max. time=',when_to_kill, 'sec ', &
!             'exceeded t_stop =', t_stop, 'sec ' 
!           exit
!        endif         
        

        loglikk = yamoeba(1)
        dloglikk = ologlikk - loglikk

        j = 0
        do ii = 1, ma
            if (ia(ii).ne.0) then
                j = j + 1
                a(ii) = p(1, j)
            endif
        enddo
    enddo

    idset = 1
    chisq = 0.d0
    loglik = 0.d0

    if (hkl.eq.0) then
        do i = 1, npl
            j = 6 * (i - 1)

            if (a(j + 2)<0.d0) then  !if P<0, set P>0
                a(j + 2) = abs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then  !if K<0, set K>0 and w = w+PI
                a(j + 4) = a(j + 4) + PI
                a(j + 1) = abs(a(j + 1))
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
            endif
            if (a(j + 3)<0.d0) then  !if e<0, set e>0 and w=w+PI, M0=M0-PI
                a(j + 3) = abs(a(j + 3))
                a(j + 4) = a(j + 4) + PI
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
                a(j + 5) = a(j + 5) - PI
                if (a(j + 5)<0.d0) a(j + 5) = a(j + 5) + 2.d0 * PI
            endif
            if (a(j + 4)<0.d0) a(j + 4) = dmod(a(j + 4) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            !if (a(j+6)<0.d0) a(j+6) = dmod(a(j+6)+2.d0*PI, 2.d0*PI)

            if (a(j + 4)>2.d0 * PI) a(j + 4) = dmod(a(j + 4), 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            !if (a(j+6)>2.d0*PI) a(j+6) = dmod(a(j+6), 2.d0*PI )
        enddo
    else
        do i = 1, npl
            j = 6 * (i - 1)
            if (a(j + 1)<0.d0) then  !if K<0, set K>0 and w = w+PI
                a(j + 4) = -1.d0 * a(j + 4) !which is h = -h, k = -k
                a(j + 3) = -1.d0 * a(j + 3)
                a(j + 1) = abs(a(j + 1))
            endif

            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            if (a(j + 6)<0.d0) a(j + 6) = dmod(a(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 6)>2.d0 * PI) a(j + 6) = dmod(a(j + 6), 2.d0 * PI)
        enddo
    endif

    do i = 1, ndata
        idset = ts(i)
        call RVKEP_kepamo (x(i), a, ymod(i), ymod_pl(:,i), dyda, ma, idset, hkl)

        y_in(i) = y(i) - a(6 * npl + idset) - a(6 * npl + 2 * ndset + 1) * x(i) - &
                a(6 * npl + 2 * ndset + 2) * x(i)**2
        ymod(i) = ymod(i) - a(6 * npl + idset)&
                - a(6 * npl + 2 * ndset + 1) * x(i) - &
                a(6 * npl + 2 * ndset + 2) * x(i)**2

        dy = y_in(i) - ymod(i)

        if (writeflag_RV>0) then

            res_array(i, :) = (/ x0 + x(i), &
                    y_in(i) + a(6 * npl + 2 * ndset + 1) * x(i) + &
                            a(6 * npl + 2 * ndset + 2) * x(i)**2, &
                     sig(i),  dble(idset) , dy, ymod(i), &
                            (ymod_pl(j,i), j=1,npl)/)                              
                            
        endif

        sig2i = 1.d0 / (sig(i)**2 + a(6 * npl + ndset + idset)**2)
        chisq = chisq + dy * dy * sig2i
        loglik = loglik - 0.5 * dy * dy * sig2i - &
                0.5 * dlog(twopi * (sig(i)**2&
                        + a(6 * npl + ndset + idset)**2))
        rms = rms + dy**2
    enddo

    rms = dsqrt(rms / dble(ndata))
    fit_return = (/ loglik, chisq / dble(ndata - mfit), chisq, rms /)

    call MA_J_kepamo (a, ma, npl, st_mass, mass, ap, hkl, gr_flag)

    if (writeflag_best_par>0) then
        do j = 1, npl + 1
            j_mass(j) = mass(j) / 1.26686534d17
        enddo

        do j = 1, npl
            i = 6 * (j - 1)

            if (hkl.eq.0) then
                best_w = a(i + 4) * 180.d0 / PI
                best_we = dsqrt(covar(i + 4, i + 4)) * 180.d0 / PI
            else
                best_w = a(i + 4)
                best_we = dsqrt(covar(i + 4, i + 4))
            endif

            bestpar_1(j, :, 1) = (/ a(i + 1), a(i + 2), a(i + 3), best_w, &
                    a(i + 5) * 180.d0 / PI, incl(j), cap0m(j), a(i + 6) * 180.d0 / PI, &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
            bestpar_1(j, :, 2) = (/ dsqrt(covar(i + 1, i + 1)), &
                    dsqrt(covar(i + 2, i + 2)), dsqrt(covar(i + 3, i + 3)), &
                    best_we, &
                    dsqrt(covar(i + 5, i + 5)) * 180.d0 / PI, 0.d0, 0.d0, &
                    dsqrt(covar(i + 6, i + 6)) * 180.d0 / PI, &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
        enddo

        do j = 1, ndset
            i = 6 * npl + j
            bestpar_2(j, :) = (/ a(i), dsqrt(covar(i, i)) /)

            bestpar_3(j, :) = (/ a(6 * npl + ndset + j), 0.d0 /)
        enddo

        bestpar_4 = (/ a(6 * npl + 2 * ndset + 1), &
                dsqrt(covar(6 * npl + ndset + 1, 6 * npl + ndset + 1)), &
                a(6 * npl + 2 * ndset + 2), &
                dsqrt(covar(6 * npl + ndset + 2, 6 * npl + ndset + 2)), &
                dble(ndata), dble(mfit), rms, &
                chisq / dble(ndata - mfit), x0, (j_mass(i + 1), i = 1, npl), &
                (ap(i) / 1.49597892d11, i = 1, npl) /)
    endif

    if(writeflag_fit>0) then
        dt = (x(ndata) + model_max + model_min) / dble(nt - 1)
        do i = 1, nt
            x(i) = ((i - 1) * dt) - model_min
            do j = 1, ndset
                a(6 * npl + j) = 0.0
            enddo
            call RVKEP_kepamo (x(i), a, ymod(i), ymod_pl2(:,i), dyda, ma, 1, hkl)
!            write(*,*)(ymod_pl(j,i),j=1,npl) 
            fit_array(i, :) = (/ x0 + x(i), ymod(i), (ymod_pl2(j,i),j=1,npl) /)
        enddo
    endif
end






subroutine kepfit_lm(epsil, deltat, amoebastarts, &
        when_to_kill, nt, &
        model_max, model_min, gr_flag_in, &
        st_mass, writeflag_best_par, writeflag_RV, &
        writeflag_fit, &
        ndset_in, &
        ndata, data_array, files_param, &
        npl_in, array_npl, &
        final_params, &
        res_array, fit_return, &
        bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
        fit_array, dynamical_planets, coplar_inc)
    implicit none
 


    integer, intent(in) :: epsil, deltat, amoebastarts,npl_in
    integer, intent(in) :: writeflag_best_par, writeflag_RV
    integer, intent(in) :: writeflag_fit, gr_flag_in, ndset_in
    integer, intent(in) :: coplar_inc, ndata, nt
       
    real(4), intent(in) :: when_to_kill, model_max, model_min    

    real(8), intent(in) :: st_mass

    integer :: hkl, gr_flag
 
       
    integer :: npl, ndset, idset,  ma, mfit, i, j, NDSMAX, NPLMAX, MMAX      
    parameter (NDSMAX = 20, NPLMAX = 10, MMAX = 200)
    
    integer :: ii, iter, idsmax(NDSMAX), ia(MMAX) 
    
    integer, allocatable, dimension(:) :: ts
    real(8), allocatable, dimension(:) :: x, y, sig, y_in, ymod

    real(8), allocatable, dimension(:,:) :: ymod_pl,ymod_pl2
     

    real(8) :: a(MMAX), covar(MMAX, MMAX), alpha(MMAX, MMAX)
    real(8) :: rms, mass(NPLMAX), ap(NPLMAX)
    real(8) :: j_mass(NPLMAX), chisq
    real(8) :: x0, incl(NPLMAX), cap0m(NPLMAX)
    real(8) :: dt, t_max, loglik, dy, sig2i
    real(8) :: epoch, ftol, jitt(NDSMAX)
    real(8) :: dyda(MMAX), p(MMAX + 1, MMAX), yamoeba(MMAX + 1)
    real(8) ::  best_w, best_we

!    character(80) version_input, version
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX)

    real(4) :: t_stop, t_init
    real(8) :: PI, twopi, alamda, ochisq, dchisq      
      
    integer, intent(in) :: dynamical_planets(npl_in)
    
    real(8), intent(in) :: data_array(ndata, 4)
    real(8), intent(in) :: files_param(ndset_in, 4)
    real(8), intent(in) :: array_npl(npl_in, 17, 2)
    real(8), intent(in) :: final_params(6)
    real(8), intent(out) :: res_array(ndata, 6+npl_in)
    real(8), intent(out) :: fit_return(4), fit_array(nt, 2+npl_in)
    real(8), intent(out) :: bestpar_1(npl_in, 17, 2), bestpar_2(ndset_in, 2)
    real(8), intent(out) :: bestpar_3(ndset_in, 2), bestpar_4(9 + 2 * npl_in)
    character(80) :: version
    character(20) :: mode  

    external rvkep_keplm

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    
    allocate(x(20000),y(20000),sig(20000),y_in(20000),ts(20000),ymod(20000))
    allocate(ymod_pl(npl_in,20000), ymod_pl2(npl_in,20000)) 
 
 




!f2py intent(in) ndset_in, npl_in, ndata
!f2py intent(in) dynamical_planets
!f2py intent(out) res_array, fit_return, fit_array
!f2py intent(out) bestpar_1, bestpar_2, bestpar_3, bestpar_4
!f2py depend(ndset_in) files_array,files_param
!f2py depend(npl_in) array_npl, dynamical_planets
!f2py depend(ndata) data_array
!f2py depend(dt) fit_array




    PI = 3.14159265358979d0
    twopi = 2.d0 * PI
        
    gr_flag = gr_flag_in
    ndset = ndset_in
    npl = npl_in
    mode = "lm"
!    epsil = 0
!    deltat = 0
!    dynamical_planets(:) = 0
!    coplar_inc = 0

    version = "1.05"

!    CALL getarg(1, version_input)
!    if(version_input=='-version') then
!        write(*, *) version
!        return
!    endif

    loglik = 0.d0

    call io_read_data(ndata, x, ts, y, sig, jitt, epoch, &
            x0, t_max, a, ia, ma, mfit, incl, cap0m, wdot, u_wdot, hkl, mode, &
            data_array, files_param, &
            array_npl, &
            final_params, coplar_inc)
    ! Hack to make it compatible with the loglik code
    if (amoebastarts.eq.0) then
        alamda = -1.d0
        mfit = 0
        do j = 1, ma
            ia(j) = 0
        enddo

        call MRQMIN_dynamo (x, y, sig, ndata, a, ia, ma, ts, covar, alpha, MMAX, &
                chisq, rvkep_keplm, alamda, loglik, jitt, hkl)
    else
        alamda = -1.d0
        call MRQMIN_dynamo (x, y, sig, ndata, a, ia, ma, ts, covar, alpha, MMAX, &
                chisq, rvkep_keplm, alamda, loglik, jitt, hkl)

        i = 0
        dchisq = 10.d0
        ochisq = chisq
!        t_init = time()
        
        call cpu_time(t_init)            
        do while ((chisq>=ochisq).or.(dchisq<-1.d-2))
            i = i + 1
!            write(*,*) alamda
            ochisq = chisq
            call MRQMIN_dynamo (x, y, sig, ndata, a, ia, ma, ts, covar, alpha, MMAX, &
                    chisq, rvkep_keplm, alamda, loglik, jitt, hkl)

            dchisq = chisq - ochisq
!            if (i.eq.10) then
!                i = 0
!            endif
            if ((i.eq.200).or.(alamda>=1d6)) then
                i = 0
                exit
            endif


!            t_stop = time() - t_init
            call cpu_time(t_stop)     
            if ((t_stop-t_init)>=when_to_kill) then
                write(*, *) 'Max. time=', when_to_kill, 'sec ', &
                        'exceeded t_stop =', t_stop, 'sec '
                exit
            endif
        enddo
    endif

    alamda = 0.d0
    call MRQMIN_dynamo (x, y, sig, ndata, a, ia, ma, ts, covar, alpha, MMAX, &
            chisq, rvkep_keplm, alamda, loglik, jitt, hkl)

    idset = 1
    rms = 0.d0

    do i = 1, ndata
        idset = ts(i)
        call RVKEP_keplm (x(i), a, ymod(i), ymod_pl(:,i), dyda, ma, ts(i), hkl)
!        xmax = x0 + x(i)

        y_in(i) = y(i) - a(6 * npl + idset)&
                - a(6 * npl + ndset + 1) * x(i)&
                - a(6 * npl + ndset + 2) * x(i)**2

        ymod(i) = ymod(i) - a(6 * npl + idset)&
                - a(6 * npl + ndset + 1) * x(i)&
                - a(6 * npl + ndset + 2) * x(i)**2

        dy = y_in(i) - ymod(i)
        rms = rms + (dy)**2
        if (writeflag_RV>0) then

                                                       
            res_array(i, :) = (/ x0 + x(i),  &
                    y_in(i) + a(6 * npl + ndset + 1) * x(i) + &
                            a(6 * npl + ndset + 2) * x(i)**2, &
                     sig(i),  dble(idset), dy ,ymod(i),  &
                            (ymod_pl(j,i), j=1,npl)/)                                
        endif
    enddo

    rms = dsqrt(rms / dble(ndata))
    fit_return = (/ loglik, chisq / dble(ndata - mfit), chisq, rms /)

    call MA_J_keplm (a, ma, npl, st_mass, mass, ap, hkl, gr_flag)

    if (writeflag_best_par>0) then
        do j = 1, npl + 1
            j_mass(j) = mass(j) / 1.26686534d17
        enddo

        do j = 1, npl
            i = 6 * (j - 1)

            if (hkl.eq.0) then
                best_w = a(i + 4) * 180.d0 / PI
                best_we = dsqrt(covar(i + 4, i + 4)) * 180.d0 / PI
            else
                best_w = a(i + 4)
                best_we = dsqrt(covar(i + 4, i + 4))
            endif

            bestpar_1(j, :, 1) = (/ a(i + 1), a(i + 2), a(i + 3), best_w, &
                    a(i + 5) * 180.d0 / PI, incl(j), cap0m(j), a(i + 6) * 180.d0 / PI, &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
            bestpar_1(j, :, 2) = (/ dsqrt(covar(i + 1, i + 1)), &
                    dsqrt(covar(i + 2, i + 2)), dsqrt(covar(i + 3, i + 3)), &
                    best_we, &
                    dsqrt(covar(i + 5, i + 5)) * 180.d0 / PI, 0.d0, 0.d0, &
                    dsqrt(covar(i + 6, i + 6)) * 180.d0 / PI, &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
        enddo

        do j = 1, ndset
            i = 6 * npl + j
            bestpar_2(j, :) = (/ a(i), dsqrt(covar(i, i)) /)

            bestpar_3(j, :) = (/ jitt(j), 0.d0 /)
        enddo

        bestpar_4 = (/ a(6 * npl + ndset + 1), &
                dsqrt(covar(6 * npl + ndset + 1, 6 * npl + ndset + 1)), &
                a(6 * npl + ndset + 2), &
                dsqrt(covar(6 * npl + ndset + 2, 6 * npl + ndset + 2)), &
                dble(ndata), dble(mfit), rms, &
                chisq / dble(ndata - mfit), x0, (j_mass(i + 1), i = 1, npl), &
                (ap(i) / 1.49597892d11, i = 1, npl) /)
    endif

    if(writeflag_fit>0) then
        dt = (x(ndata) + model_max + model_min) / dble(nt - 1)
        do i = 1, nt
            x(i) = ((i - 1) * dt) - model_min

            do j = 1, ndset
                a(6 * npl + j) = 0.0
            enddo

            call RVKEP_keplm (x(i), a, ymod(i), ymod_pl(:,i), dyda, ma, 1, hkl)
!            write(*,*)(ymod_pl(j,i),j=1,npl) 
            fit_array(i, :) = (/ x0 + x(i), ymod(i), (ymod_pl(j,i),j=1,npl) /)
        enddo
    endif
end



subroutine dynfit_amoeba__2(epsil, deltat, amoebastarts, &
        when_to_kill, nt, &
        model_max, model_min, gr_flag_in, &
        st_mass, writeflag_best_par, writeflag_RV, &
        writeflag_fit, &
        ndset_in, &
        ndata, data_array, files_param, &
        npl_in, array_npl, &
        final_params, &
        res_array, fit_return, &
        bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
        fit_array, dynamical_planets, coplar_inc)
    implicit none
    real(8) :: PI
    parameter (PI = 3.14159265358979d0)
    integer :: npl, ndset, idset, ndata, ma, mfit, i, j, NDSMAX, NPLMAX, MMAX
    parameter (NDSMAX = 20, NPLMAX = 10, MMAX = 200)
    integer :: idsmax(NDSMAX), ia(MMAX), ts(20000), nt, iter, ii
    integer :: writeflag_best_par, hkl
    integer :: writeflag_RV, writeflag_fit, amoebastarts
    real(8) :: mstar, sini(NPLMAX), loglik, ftol
    real(8) :: t(20000), ys(20000), sigs(20000)
    real(8) :: a(MMAX), covar(MMAX, MMAX)
    real(8) :: chisq, p(MMAX + 1, MMAX)
    real(8) :: yamoeba(MMAX + 1), loglikk, ologlikk, dloglikk
    real(8) :: t0, t_max, epoch, epsil, deltat, st_mass
    real(8) :: rms, jitt(NDSMAX)
    real(8) :: t_stop, t_init
    real(4) :: when_to_kill, model_max, model_min
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX)
    real(8) :: incl(NPLMAX), cap0m(NPLMAX)

    external rvkep_dynamo, compute_abs_loglik_dyn
    character(80) :: version
!    character(80) version_input, version

    integer :: ndset_in, npl_in, gr_flag_in, gr_flag, coplar_inc
    integer :: dynamical_planets(npl_in)
    real(8) :: data_array(ndata, 4)
    real(8) :: files_param(ndset_in, 4)
    real(8) :: array_npl(npl_in, 17, 2)
    real(8) :: final_params(6)
    real(8) :: res_array(ndata, 6+npl_in)
    real(8) :: fit_return(4), fit_array(nt, 2+npl_in)
    real(8) :: bestpar_1(npl_in, 17, 2), bestpar_2(ndset_in, 2)
    real(8) :: bestpar_3(ndset_in, 2), bestpar_4(9 + 2 * npl_in)
    character(20) :: mode
    
!f2py intent(in) ndset_in, npl_in, ndata
!f2py intent(in) dynamical_planets
!f2py intent(out) res_array, fit_return, fit_array
!f2py intent(out) bestpar_1, bestpar_2, bestpar_3, bestpar_4
!f2py depend(ndset_in) files_array,files_param
!f2py depend(npl_in) array_npl, dynamical_planets
!f2py depend(ndata) data_array
!f2py depend(dt) fit_array


    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    common mstar, sini

    gr_flag = gr_flag_in
    ndset = ndset_in
    npl = npl_in
    mode = "amoeba_dyn"
    mstar = st_mass
    covar(:, :) = 0

    version = "1.05"

!    CALL getarg(1, version_input)
!    if(version_input=='-version') then
!        write(*, *) version
!        return
!    endif

    loglikk = 0
    ftol = 0.001d0

    call io_read_data (ndata, t, ts, ys, sigs, jitt, &
            epoch, t0, t_max, a, ia, ma, mfit, incl, cap0m, wdot, u_wdot, hkl, mode, &
            data_array, files_param, &
            array_npl, &
            final_params, coplar_inc)
            
                    

    i = 0
    dloglikk = 10.d0
!    t_init = time()
    call cpu_time(t_init)         
    
    do while (dabs(dloglikk)>=0.000001d0)
        if (i.eq.amoebastarts) then
            i = 0
            exit
        endif
        i = i + 1
        ologlikk = loglikk

        call prepare_for_amoeba_dyn(p, MMAX + 1, MMAX, yamoeba, a, ia, ma, &
                mfit, compute_abs_loglik_dyn, ndata, t, ys, &
                ts, sigs, epsil, deltat, hkl, dynamical_planets, coplar_inc)

        call amoeba_dyn(p, yamoeba, MMAX + 1, MMAX, mfit, ftol, &
                compute_abs_loglik_dyn, &
                iter, ndata, t, ys, ma, ts, sigs, a, ia, epsil, deltat, hkl, &
                npl, dynamical_planets, coplar_inc)

!        t_stop = time() - t_init
        call cpu_time(t_stop)     
        if ((t_stop-t_init)>=when_to_kill) then
            write(*, *) 'Max. time=', when_to_kill, 'sec ', &
                    'exceeded t_stop =', t_stop, 'sec '
            exit
        endif

        loglikk = yamoeba(1)
        dloglikk = ologlikk - loglikk

        j = 0
        do ii = 1, ma
            if (ia(ii).ne.0) then
                j = j + 1
                a(ii) = p(1, j)
            endif
        enddo
 
     
    enddo

!    do ii = 1, ma
!        if (ia(ii).ne.0) then
!            j = j + 1
!            write(*,*) a(ii), p(1, j)
!        endif
!    enddo
!     do ii = 1, ma
! 
!             write(*,*) a(ii) 
! 
!    enddo

    j = 0
    loglik = 0.0d0
    chisq = 0.0d0
    rms = 0.0d0
     


    call io_write_bf_ewcop_fin_dynamo(a, covar, t, ys, ndata, ts, &
            ma, mfit, t0, t_max, sigs, chisq, rms, loglik, writeflag_RV, &
            writeflag_best_par, writeflag_fit, epsil, deltat, &
            nt, model_max, model_min, hkl, wdot, u_wdot, &
            dynamical_planets, &
            res_array, fit_return, &
            bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
            fit_array, coplar_inc)
return            
end





subroutine dynfit_amoeba(epsil, deltat, amoebastarts, &
        when_to_kill, nt, &
        model_max, model_min, gr_flag_in, &
        st_mass, writeflag_best_par, writeflag_RV, &
        writeflag_fit, &
        ndset_in, &
        ndata, data_array, files_param, &
        npl_in, array_npl, &
        final_params, &
        res_array, fit_return, &
        bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
        fit_array, dynamical_planets, coplar_inc)
    implicit none
    

    integer, intent(in) :: amoebastarts,npl_in
    integer, intent(in) :: writeflag_best_par, writeflag_RV
    integer, intent(in) :: writeflag_fit, gr_flag_in, ndset_in
    integer, intent(in) :: coplar_inc, ndata, nt
       
    real(4), intent(in) :: when_to_kill, model_max, model_min    

    real(8), intent(in) :: st_mass, epsil, deltat

    integer :: hkl, gr_flag
 
       
    integer :: npl, ndset, idset,  ma, mfit, i, j, NDSMAX, NPLMAX, MMAX      
    parameter (NDSMAX = 20, NPLMAX = 10, MMAX = 200)
    
    integer :: ii, iter, idsmax(NDSMAX), ia(MMAX),nbod,ts(20000)
    
!    integer, allocatable, dimension(:) :: ts
!    real(8), allocatable, dimension(:) :: x, t, ys, sigs, y_in, ymod
!    real(8), allocatable, dimension(:,:) :: ymod_pl,ymod_pl2

    real(8) :: t(20000),y(20000),ys(20000), sigs(20000)        
!    real(8) ::  x(20000), y_in(20000), ymod(20000)     
    real(8) :: a(MMAX), covar(MMAX, MMAX)
    real(8) :: rms, mass(NPLMAX), ap(NPLMAX)
    real(8) :: j_mass(NPLMAX), chisq
    real(8) :: x0, incl(NPLMAX), cap0m(NPLMAX)
    real(8) :: dt, t_max, loglik, dy, sig2i
    real(8) :: epoch, ftol, jitt(NDSMAX)
    real(8) :: dyda(MMAX), p(MMAX + 1, MMAX), yamoeba(MMAX + 1)

    real(8) :: loglikk, ologlikk, dloglikk, best_w, best_we, t0
!    real(8) :: xj(NPLMAX), yj(NPLMAX),  zj(NPLMAX) 
!    real(8) :: vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
!    real(8) :: rpl(NPLMAX), rhill(NPLMAX) 
    
!    character(80) version_input, version
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX)

    real(4) :: t_stop, t_init
    real(8) :: PI    
      
    integer, intent(in) :: dynamical_planets(npl_in)
    
    real(8), intent(in) :: data_array(ndata, 4)
    real(8), intent(in) :: files_param(ndset_in, 4)
    real(8), intent(in) :: array_npl(npl_in, 17, 2)
    real(8), intent(in) :: final_params(6)
    real(8), intent(out) :: res_array(ndata, 6+npl_in)
    real(8), intent(out) :: fit_return(4), fit_array(nt, 2+npl_in)
    real(8), intent(out) :: bestpar_1(npl_in, 17, 2), bestpar_2(ndset_in, 2)
    real(8), intent(out) :: bestpar_3(ndset_in, 2), bestpar_4(9 + 2 * npl_in)
    character(80) :: version
    character(20) :: mode
        
 
    real(8) :: mstar, sini(NPLMAX)  
!    character(80) version_input, version

 
    external RVKEP_dynamo, compute_abs_loglik_dyn 
    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    common mstar, sini
     
!f2py intent(in) ndset_in, npl_in, ndata
!f2py intent(in) dynamical_planets
!f2py intent(out) res_array, fit_return, fit_array
!f2py intent(out) bestpar_1, bestpar_2, bestpar_3, bestpar_4
!f2py depend(ndset_in) files_array,files_param
!f2py depend(npl_in) array_npl, dynamical_planets
!f2py depend(ndata) data_array


!    allocate(x(20000),t(20000),ys(20000),sigs(20000))
!    allocate(y_in(20000), ymod(20000))
!    allocate(ymod_pl(npl_in,20000), ymod_pl2(npl_in,20000)) 


    gr_flag = gr_flag_in
    ndset = ndset_in
    npl = npl_in
    mode = "amoeba_dyn"
    mstar = st_mass
    covar(:, :) = 0
    PI = 3.14159265358979d0
    version = "1.05"

!    CALL getarg(1, version_input)
!    if(version_input=='-version') then
!        write(*, *) version
!        return
!    endif

    loglikk = 0
    ftol = 0.001d0

    call io_read_data (ndata, t, ts, ys, sigs, jitt, &
            epoch, t0, t_max, a, ia, ma, mfit, incl, cap0m, wdot, u_wdot, hkl, mode, &
            data_array, files_param, &
            array_npl, &
            final_params, coplar_inc)

    i = 0
    dloglikk = 10.d0
!    t_init = time()
    call cpu_time(t_init)         
  

    
    do while (dabs(dloglikk)>=0.000001d0)
        if (i.eq.amoebastarts) then
            i = 0
            exit
        endif
        i = i + 1
        ologlikk = loglikk

        call prepare_for_amoeba_dyn(p, MMAX + 1, MMAX, yamoeba, a, ia, ma, &
                mfit, compute_abs_loglik_dyn, ndata, t, ys, &
                ts, sigs, epsil, deltat, hkl, dynamical_planets, coplar_inc)

        call amoeba_dyn(p, yamoeba, MMAX + 1, MMAX, mfit, ftol, &
                compute_abs_loglik_dyn, &
                iter, ndata, t, ys, ma, ts, sigs, a, ia, epsil, deltat, hkl, &
                npl, dynamical_planets, coplar_inc)

!        t_stop = time() - t_init
        call cpu_time(t_stop)     
        if ((t_stop-t_init)>=when_to_kill) then
            write(*, *) 'Max. time=', when_to_kill, 'sec ', &
                    'exceeded t_stop =', t_stop, 'sec '
            exit
        endif

        loglikk = yamoeba(1)
        dloglikk = ologlikk - loglikk

        j = 0
        do ii = 1, ma
            if (ia(ii).ne.0) then
                j = j + 1
                a(ii) = p(1, j)
            endif
        enddo
 
     
    enddo

 
    j = 0
    loglik = 0.0d0
    chisq = 0.0d0
    rms = 0.0d0
    

    call io_write_bf_ewcop_fin_dynamo(a, covar, t, ys, ndata, ts, &
            ma, mfit, t0, t_max, sigs, chisq, rms, loglik, writeflag_RV, &
            writeflag_best_par, writeflag_fit, epsil, deltat, &
            nt, model_max, model_min, hkl, wdot, u_wdot, &
            dynamical_planets, &
            res_array, fit_return, &
            bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
            fit_array, coplar_inc)

 
return
end

 




subroutine dynfit_lm(epsil, deltat, amoebastarts, &
        when_to_kill, nt, &
        model_max, model_min, gr_flag_in, &
        st_mass, writeflag_best_par, writeflag_RV, &
        writeflag_fit, &
        ndset_in, &
        ndata, data_array, files_param, &
        npl_in, array_npl, &
        final_params, &
        res_array, fit_return, &
        bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
        fit_array, dynamical_planets, coplar_inc)
    implicit none
    


    integer, intent(in) :: amoebastarts,npl_in
    integer, intent(in) :: writeflag_best_par, writeflag_RV
    integer, intent(in) :: writeflag_fit, gr_flag_in, ndset_in
    integer, intent(in) :: coplar_inc, ndata, nt
       
    real(4), intent(in) :: when_to_kill, model_max, model_min    

    real(8), intent(in) :: st_mass, epsil, deltat

    integer :: hkl, gr_flag
 
       
    integer :: npl, ndset, idset,  ma, mfit, i, j, NDSMAX, NPLMAX, MMAX      
    parameter (NDSMAX = 20, NPLMAX = 10, MMAX = 200)
    
    integer :: ii, iter, idsmax(NDSMAX), ia(MMAX),nbod,ts(20000)
    
!    integer, allocatable, dimension(:) :: ts
!    real(8), allocatable, dimension(:) :: x, t, ys, sigs, y_in, ymod
!    real(8), allocatable, dimension(:,:) :: ymod_pl,ymod_pl2

    real(8) :: t(20000),y(20000),ys(20000), sigs(20000)        
!    real(8) ::  x(20000), y_in(20000), ymod(20000)     
    real(8) :: a(MMAX), covar(MMAX, MMAX)
    real(8) :: rms, mass(NPLMAX), ap(NPLMAX)
    real(8) :: j_mass(NPLMAX), chisq,alamda, ochisq, dchisq
    real(8) :: x0, incl(NPLMAX), cap0m(NPLMAX)
    real(8) :: dt, t_max, loglik, dy, sig2i
    real(8) :: epoch, ftol, jitt(NDSMAX), alpha(MMAX, MMAX)
    real(8) :: dyda(MMAX), p(MMAX + 1, MMAX), yamoeba(MMAX + 1)

    real(8) ::  best_w, best_we, t0
!    real(8) :: xj(NPLMAX), yj(NPLMAX),  zj(NPLMAX) 
!    real(8) :: vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
!    real(8) :: rpl(NPLMAX), rhill(NPLMAX) 
    
!    character(80) version_input, version
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX)

    real(4) :: t_stop, t_init
    real(8) :: PI    
      
    integer, intent(in) :: dynamical_planets(npl_in)
    
    real(8), intent(in) :: data_array(ndata, 4)
    real(8), intent(in) :: files_param(ndset_in, 4)
    real(8), intent(in) :: array_npl(npl_in, 17, 2)
    real(8), intent(in) :: final_params(6)
    real(8), intent(out) :: res_array(ndata, 6+npl_in)
    real(8), intent(out) :: fit_return(4), fit_array(nt, 2+npl_in)
    real(8), intent(out) :: bestpar_1(npl_in, 17, 2), bestpar_2(ndset_in, 2)
    real(8), intent(out) :: bestpar_3(ndset_in, 2), bestpar_4(9 + 2 * npl_in)
    character(80) :: version
    character(20) :: mode
        
 
    real(8) :: mstar, sini(NPLMAX)  
!    character(80) version_input, version

    external rvkep_ewcop_fin 
    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    common mstar, sini    
 
 
 
!f2py intent(in) ndset_in, npl_in, ndata
!f2py intent(in) dynamical_planets
!f2py intent(out) res_array, fit_return, fit_array
!f2py intent(out) bestpar_1, bestpar_2, bestpar_3, bestpar_4
!f2py depend(ndset_in) files_array,files_param
!f2py depend(npl_in) array_npl, dynamical_planets
!f2py depend(ndata) data_array
!f2py depend(dt) fit_array
 
    PI = 3.14159265358979d0
    gr_flag = gr_flag_in
    ndset = ndset_in
    npl = npl_in
    mode = "lm_dyn"
    mstar = st_mass
    covar(:, :) = 0
!    dynamical_planets(:) = 0

    version = "1.05"

!    CALL getarg(1, version_input)
!    if(version_input=='-version') then
!        write(*, *) version
!        return
!    endif

    loglik = 0.d0
    rms = 0.d0

    call io_read_data (ndata, t, ts, ys, sigs, jitt, &
            epoch, t0, t_max, a, ia, ma, mfit, incl, cap0m, wdot, u_wdot, hkl, mode, &
            data_array, files_param, &
            array_npl, &
            final_params, coplar_inc)

    !*****set alamda to be negtive for initializing******
    alamda = -1.d0
    call MRQMIN_dynlm (t, ts, ys, sigs, ndata, a, ia, ma, covar, alpha, MMAX, &
            chisq, rvkep_ewcop_fin, alamda, loglik, jitt, epsil, deltat, hkl, coplar_inc)

    i = 0
    dchisq = 10.d0
    ochisq = chisq
!    t_init = time()
    call cpu_time(t_init)         
    do while ((chisq>=ochisq).or.(dchisq<-1d-5))
        if (i.eq.amoebastarts) then
            exit
        endif
        i = i + 1
        ochisq = chisq
!        write(*,*) alamda        
        call MRQMIN_dynlm (t, ts, ys, sigs, ndata, a, ia, ma, covar, alpha, &
                MMAX, chisq, rvkep_ewcop_fin, alamda, loglik, jitt, epsil, deltat, hkl, coplar_inc)

        dchisq = chisq - ochisq

        if ((i.eq.200).or.(alamda>=1d6)) then
            i = 0
            exit
        endif

!        t_stop = time() - t_init
        call cpu_time(t_stop)     
        if ((t_stop-t_init)>=when_to_kill) then
            write(*, *) 'Max. time=', when_to_kill, 'sec ', &
                    'exceeded t_stop =', t_stop, 'sec '
            exit
        endif
    enddo

    !*******final output******************
    alamda = 0.d0

    call MRQMIN_dynlm (t, ts, ys, sigs, ndata, a, ia, ma, covar, alpha, MMAX, &
            chisq, rvkep_ewcop_fin, alamda, loglik, jitt, epsil, deltat, hkl, coplar_inc)

    call io_write_bf_ewcop_fin_dynlm (a, covar, t, ys, ndata, ts, &
             ma, mfit, t0, t_max, sigs, chisq, rms, writeflag_RV, &
             writeflag_best_par, writeflag_fit, jitt, epsil, deltat, &
             nt, model_max, model_min, hkl, wdot, u_wdot, &
             res_array, fit_return, &
             bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
             fit_array, coplar_inc)
return
end
end

subroutine compute_abs_loglik_kep(ndata, x, y, a2, dyda, ma, mfit, ts, &
        sig, loglik, a, ia, hkl)
    implicit none

    integer :: MMAX, NDSMAX, npl, ndset, idset, mfit, gr_flag
    parameter (MMAX = 200, NDSMAX = 20)
    real(8) :: loglik, PI, TWOPI
    parameter (PI = 3.14159265358979d0)
    parameter (TWOPI = 2.0 * PI)
    integer :: ndata, i, j, ma, ts(ndata), ia(MMAX), idsmax(NDSMAX), hkl
    real(8) :: dy, sig(ndata), dyda(MMAX), x(ndata), y(ndata)
    real(8) :: a(MMAX), a2(mfit), a3(MMAX), sig2i, y_in(ndata),&
            y2(ndata), y2_pl(10,ndata)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    loglik = 0.d0
    j = 0
    do i = 1, ma
        if (ia(i).ne.0) then
            j = j + 1
            a3(i) = a2(j)
        else
            a3(i) = a(i)
        endif
    enddo

    do i = 1, ndata
        idset = ts(i)
        call RVKEP_kepamo (x(i), a3, y2(i), y2_pl(:,i), dyda, ma, idset, hkl)
        y_in(i) = y(i) - a3(6 * npl + idset) - &
                a3(6 * npl + 2 * ndset + 1) * x(i)&
                - a3(6 * npl + 2 * ndset + 2) * x(i)**2

        y2(i) = y2(i) - a3(6 * npl + idset) - &
                a3(6 * npl + 2 * ndset + 1) * x(i)&
                - a3(6 * npl + 2 * ndset + 2) * x(i)**2

        dy = y_in(i) - y2(i)

        sig2i = 1.d0 / (sig(i)**2 + a3(6 * npl + ndset + idset)**2)

        loglik = loglik + 0.5 * dy * dy * sig2i + &
                dlog(dsqrt(TWOPI * (sig(i)**2 + &
                        a3(6 * npl + ndset + idset)**2)))
    enddo
    loglik = loglik - ndata * dlog(dsqrt(TWOPI))
    return
end

subroutine compute_abs_loglik_dyn(ndata, x, y, a2, ma, mfit, ts, &
        sig, loglik, a, ia, epsil, deltat, hkl, dynamical_planets, &
        coplar_inc)

    implicit none
    integer :: MMAX, npl, ndset, NDSMAX, idset, mfit, hkl, gr_flag,ndata
    parameter (MMAX = 200, NDSMAX = 20)
    real(8) :: twopi, loglik
    parameter (twopi = 6.28318530717958d0)
    integer  :: i, j, ma, ts(20000), ia(MMAX), idsmax(NDSMAX)
    
    real(8) :: dy, sig(ndata), x(ndata), y(ndata), y_pl(npl,ndata)
    real(8) :: a(MMAX), a2(mfit), a3(MMAX), sig2i, y_in(ndata)
    real(8) :: y2(ndata), epsil, deltat
    integer dynamical_planets(npl), coplar_inc

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    j = 0
    do i = 1, ma
        if (ia(i).ne.0) then
            j = j + 1
            a3(i) = a2(j)
        else
            a3(i) = a(i)
        endif
    enddo
     

    loglik = 0.d0

    call RVKEP_dynamo (x, a3, y2, y_pl, ma, ndata, epsil, deltat, hkl, &
            dynamical_planets, ts, coplar_inc)


    do i = 1, ndata
        idset = ts(i)
        y_in(i) = y(i) - a3(7 * npl + idset)&
                - a3(7 * npl + 2 * ndset + 1) * (x(i) / 86400.d0)&
                - a3(7 * npl + 2 * ndset + 2) * (x(i) / 86400.d0)**2

        dy = y_in(i) - y2(i)

        sig2i = 1.d0 / (sig(i)**2 + a3(7 * npl + ndset + idset)**2)

        loglik = loglik + 0.5 * dy * dy * sig2i + &
                dlog(dsqrt(twopi * (sig(i)**2 + &
                        a3(7 * npl + ndset + idset)**2)))
    enddo
    loglik = loglik - ndata * dlog(dsqrt(TWOPI))
    return
end

subroutine io_read_data(ndata, t, ts, ys, sigs, jitt, epoch, t0, t_max, &
        ar, iar, ma, mfit, incl, cap0m, wdot, u_wdot, hkl, mode, &
        data_array, files_param, &
        array_npl, &
        final_params, coplar_inc)
    implicit none
    integer :: ndset, idset, ndata, NDSMAX, NPLMAX, MMAX, npl, ma
    real(8) :: t(ndata), ys(ndata), sigs(ndata), PI
    parameter (NDSMAX = 20, NPLMAX = 10, MMAX = 200)
    parameter(PI = 3.14159265358979d0)
    real(8) :: ar(MMAX), incl(NPLMAX), cap0m(NPLMAX)
    integer :: iar(MMAX), u_off(NDSMAX), u_jit(NDSMAX), hkl
    integer :: idsmax(NDSMAX), ts(ndata), u_incl, u_cap0m, mfit
    real(8) :: jitt(NDSMAX), t0, t_max, epoch
    real(8) :: off(NDSMAX)
    integer :: i, j, gr_flag
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX)

    real(8) :: data_array(ndata, 4)
    real(8) :: files_param(ndset, 4)
    real(8) :: array_npl(npl, 17, 2)
    real(8) :: final_params(6)
    character(20) :: mode
    integer :: ndsetmode, arr_pos, coplar_inc

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    if (ndset>NDSMAX) stop ' KEPFIT: ndset > NDSMAX.'

    off(1:ndset) = files_param(:, 1)
    u_off(1:ndset) = int(files_param(:, 2))
    jitt(1:ndset) = files_param(:, 3)
    u_jit(1:ndset) = int(files_param(:, 4))
    t(1:ndata) = data_array(:, 1)
    ys(1:ndata) = data_array(:, 2)
    sigs(1:ndata) = data_array(:, 3)
    ts(1:ndata) = int(data_array(:, 4))
    
!    write(*,*) ts

    if (npl>NPLMAX) stop ' KEPFIT: npl > NPLMAX.'

    ndsetmode = 0
    arr_pos = 0
    if (mode.eq."amoeba") then
        ndsetmode = 2 * ndset
        arr_pos = 6
    elseif (mode.eq."lm") then
        ndsetmode = ndset
        arr_pos = 6
    elseif (mode.eq."amoeba_dyn") then
        ndsetmode = 2 * ndset
        arr_pos = 7
    elseif (mode.eq."lm_dyn") then
        ndsetmode = ndset
        arr_pos = 7
    endif

    do i = 1, ndset
        ar(arr_pos * npl + i) = off(i)
        iar(arr_pos * npl + i) = u_off(i)
        if (mode.eq."amoeba" .or. mode.eq."amoeba_dyn") then
            ar(arr_pos * npl + ndset + i) = jitt(i)
            iar(arr_pos * npl + ndset + i) = u_jit(i)
        endif
    enddo

    ma = arr_pos * npl + ndsetmode + 2
    hkl = int(final_params(6))

    do j = 1, npl
        i = arr_pos * (j - 1)
        !array_npl order:
        !K P e w Ma incl ome wdot t0 RR aR a m tw h k lamb
        ar(i + 1) = array_npl(j, 1, 1)
        iar(i + 1) = int(array_npl(j, 1, 2))
        ar(i + 2) = array_npl(j, 2, 1)
        iar(i + 2) = int(array_npl(j, 2, 2))
        if (hkl.eq.0) then
            ar(i + 3) = array_npl(j, 3, 1)
            iar(i + 3) = int(array_npl(j, 3, 2))
            ar(i + 4) = array_npl(j, 4, 1)
            iar(i + 4) = int(array_npl(j, 4, 2))
            ar(i + 5) = array_npl(j, 5, 1)
            iar(i + 5) = int(array_npl(j, 5, 2))
        else
            ar(i + 3) = array_npl(j, 15, 1)
            iar(i + 3) = int(array_npl(j, 15, 2))
            ar(i + 4) = array_npl(j, 16, 1)
            iar(i + 4) = int(array_npl(j, 16, 2))
            ar(i + 5) = array_npl(j, 17, 1)
            iar(i + 5) = int(array_npl(j, 17, 2))
        endif

        if (mode.eq."amoeba" .or. mode.eq."lm") then
            incl(j) = array_npl(j, 6, 1)
            u_incl = int(array_npl(j, 6, 2))
            cap0m(j) = array_npl(j, 7, 1)
            u_cap0m = int(array_npl(j, 7, 2))
            ar(i + 6) = array_npl(j, 8, 1)
            iar(i + 6) = int(array_npl(j, 8, 2))

            wdot(i) = 0
            u_wdot(i) = 0
        elseif (mode.eq."amoeba_dyn".or.mode.eq."lm_dyn") then
            if (coplar_inc.eq.0) then
                ar(i + 6) = array_npl(j, 6, 1)
                iar(i + 6) = int(array_npl(j, 6, 2))
            else
                if (i.eq.0)  then
                    ar(i + 6) = array_npl(1, 6, 1)
                    iar(i + 6) = int(array_npl(1, 6, 2))
                else
                    ar(i + 6) = array_npl(1, 6, 1)
                    iar(i + 6) = 0    
                endif            
            endif

            ar(i + 7) = array_npl(j, 7, 1)
            iar(i + 7) = int(array_npl(j, 7, 2))
            wdot(i) = array_npl(j, 8, 1)
            u_wdot(i) = int(array_npl(j, 7, 2))

            incl(j) = 0
            u_incl = 0
            cap0m(j) = 0
            u_cap0m = 0

            if (mode.eq."amoeba_dyn".or. mode.eq."lm_dyn") then
                ar(i + 2) = ar(i + 2) * 8.64d4 !second as unit
            endif
        endif
        !inclinations and cap0m are always ignored in the fit, just for consistency with dynamical input and output
    enddo
    !u_jit read for consistency with input/output in loglik case, but here we don't actually save this information and not use it, jitters cannot be used for fit in chi^2 minimization

    ar(arr_pos * npl + ndsetmode + 1) = final_params(1)
    iar(arr_pos * npl + ndsetmode + 1) = int(final_params(2))
    ar(arr_pos * npl + ndsetmode + 2) = final_params(3)
    iar(arr_pos * npl + ndsetmode + 2) = int(final_params(4))

    epoch = final_params(5)
    t_max = t(ndata)

    if (epoch.eq.0) then
        t0 = t(1)
    else
        t0 = epoch
    endif

    do j = 1, npl
        i = arr_pos * (j - 1)
        if (hkl.eq.0) then
            ar(i + 4) = ar(i + 4) * PI / 180.d0
        endif
        ar(i + 5) = ar(i + 5) * PI / 180.d0
        ar(i + 6) = ar(i + 6) * PI / 180.d0
        if (mode.eq."amoeba_dyn".or.mode.eq."lm_dyn") then
            ar(i + 7) = ar(i + 7) * PI / 180.d0
        endif
    enddo

    do i = 1, ndata
        t(i) = (t(i) - t0)             ! time unit is day
        if (mode.eq."amoeba_dyn".or.mode.eq."lm_dyn") then
            t(i) = t(i) * 8.64d4          ! time unit is second
        endif
    enddo

    mfit = 0
    do j = 1, ma
        if (iar(j).ne.0) mfit = mfit + 1
    enddo
    return
end

subroutine RVKEP_kepamo (x, a, y, y_pl, dyda, ma, ts, hkl)
    implicit none
    real(8) :: PI, TWOPI
    parameter (PI = 3.14159265358979d0)
    parameter (TWOPI = 2.0d0 * PI)
    integer :: npl, ndset, idset, ma, i, j, NDSMAX, ts, hkl, gr_flag
    parameter (NDSMAX = 20)
    integer :: idsmax(NDSMAX)
    real(8) :: x, y, a(ma), a2(ma), dyda(ma), mass(10), ap(10)
    real(8) :: cosw, sinw, capm, cape, cose, sine, cosf, sinf, fac1, fac2, fac3
    real(8) :: orbel_ehybrid, omega(10), capmm(10), ecc(10)
    real(8) :: wm, sinwm, coswm, sin2wm, cos2wm
    real(8) :: sin3wm, cos3wm, omegad(10)
    real(8) :: y_pl(npl)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    y = 0.d0

    do i = 1, ma
        a2(i) = a(i)
    enddo

    if (hkl.eq.0) then

        do i = 1, npl
            j = 6 * (i - 1)

            if (a2(j + 2)<0.d0) then  ! if P<0, set P>0
                a2(j + 2) = dabs(a2(j + 2))
            endif

            if (a2(j + 1)<0.d0) then  ! if K<0, set K>0 and w = w+PI
                a2(j + 4) = a2(j + 4) + PI
                a2(j + 1) = dabs(a2(j + 1))
                if (a2(j + 4)>2.d0 * PI) a2(j + 4) = a2(j + 4) - 2.d0 * PI
            endif
            if (a2(j + 3)<0.d0) then  ! if e<0, set e>0 and w=w+PI, M0=M0-PI
                a2(j + 3) = dabs(a2(j + 3))
                a2(j + 4) = a2(j + 4) + PI
                if (a2(j + 4)>2.d0 * PI) a2(j + 4) = a2(j + 4) - 2.d0 * PI
                a2(j + 5) = a2(j + 5) - PI
                if (a2(j + 5)<0.d0) a2(j + 5) = a2(j + 5) + 2.d0 * PI
            endif
            if (a2(j + 4)<0.d0) a2(j + 4) = dmod(a2(j + 4) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 5)<0.d0) a2(j + 5) = dmod(a2(j + 5) + 2.d0 * PI, 2.d0 * PI)
            !if (a2(j+6)<0.d0) a2(j+6)=dmod(a2(j+6)+2.d0*PI,2.d0*PI)

            if (a2(j + 4)>2.d0 * PI) a2(j + 4) = dmod(a2(j + 4), 2.d0 * PI)
            if (a2(j + 5)>2.d0 * PI) a2(j + 5) = dmod(a2(j + 5), 2.d0 * PI)
            !if (a2(j+6)>2.d0*PI) a2(j+6)=dmod(a2(j+6), 2.d0*PI)

            ecc(i) = a2(j + 3)
            omega(i) = a2(j + 4)
            capmm(i) = a2(j + 5)

            if(gr_flag.ne.0) call MA_J_kepamo (a, ma, npl, 1.0d0, &
                    mass, ap, hkl, gr_flag)

            omegad(i) = a(j + 6)
        enddo
    else
        do i = 1, npl
            j = 6 * (i - 1)
            if (a2(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a2(j + 4) = -1.d0 * a2(j + 4) !which is h = -h, k = -k
                a2(j + 3) = -1.d0 * a2(j + 3)
                a2(j + 1) = dabs(a2(j + 1))
            endif

            ecc(i) = dsqrt(a2(j + 3)**2 + a2(j + 4)**2)
            omega(i) = datan2(a2(j + 3), a2(j + 4))

            if(omega(i)<0.d0)omega(i) = dmod(omega(i) + 2.d0 * PI, 2.d0 * PI)
            if(omega(i)>0.d0)omega(i) = dmod(omega(i), 2.d0 * PI)
            if (a2(j + 5)<0.d0) a2(j + 5) = dmod(a2(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 5)>2.d0 * PI) a2(j + 5) = dmod(a2(j + 5), 2.d0 * PI)

            capmm(i) = a2(j + 5) - omega(i)

            if(capmm(i)<0.d0)capmm(i) = dmod(capmm(i) + 2.d0 * PI, 2.d0 * PI)
            if(capmm(i)>0.d0)capmm(i) = dmod(capmm(i), 2.d0 * PI)

            if(gr_flag.ne.0) call MA_J_kepamo (a, ma, npl, 1.0d0, &
                    mass, ap, hkl, gr_flag)            
            omegad(i) = a(j + 6)
        enddo
    endif

    if (hkl.eq.0) then
        do j = 1, npl
            i = 6 * (j - 1)
            cosw = dcos(omega(j) + omegad(j) * x / 365.25d0)
            sinw = dsin(omega(j) + omegad(j) * x / 365.25d0)

            capm = TWOPI * x / a2(2 + i) + capmm(j)
            capm = dmod(capm, 2.d0 * PI)

            cape = ORBEL_EHYBRID (ecc(j), capm)
            cose = dcos(cape)
            sine = dsin(cape)

            cosf = (cose - ecc(j)) / (1.d0 - ecc(j) * cose)
            sinf = (dsqrt(1.d0 - ecc(j)**2) * sine) / (1.d0 - ecc(j) * cose)

            fac1 = cosw * cosf - sinw * sinf + ecc(j) * cosw

            fac2 = (cosw * sinf + sinw * cosf) / (1.d0 - ecc(j) * cose)**2
            fac3 = -a2(1 + i) * dsqrt(1.d0 - ecc(j)**2) * fac2

            y_pl(j) = a2(1 + i) * fac1
            y = y + a2(1 + i) * fac1
            dyda(1 + i) = fac1
            dyda(2 + i) = -TWOPI * fac3 * x / a2(2 + i)**2
            dyda(3 + i) = -a2(1 + i) * sine * (2.d0 - ecc(j)**2 - ecc(j) * cose) * &
                    fac2 / dsqrt(1.d0 - ecc(j)**2)
            dyda(4 + i) = -a2(1 + i) * (sinw * cosf + cosw * sinf + ecc(j) * sinw)
            dyda(5 + i) = fac3
            dyda(6 + i) = -a2(1 + i) * (sinw * cosf + cosw * sinf&
                    + ecc(j) * sinw) * x / 365.25d0
        enddo
    else
        do j = 1, npl
            i = 6 * (j - 1)
            !ecc2 = dsqrt(a2(3+i)**2 + a2(4+i)**2)
            if (ecc(j)>1.d-2) then
                cosw = dcos(omega(j) + omegad(i) * x / 365.25d0)
                sinw = dsin(omega(j) + omegad(i) * x / 365.25d0)

                capm = TWOPI * x / a2(2 + i) + capmm(j)
                capm = dmod(capm, 2.d0 * PI)

                cape = ORBEL_EHYBRID (ecc(j), capm)
                cose = dcos(cape)
                sine = dsin(cape)
                cosf = (cose - ecc(j)) / (1.d0 - ecc(j) * cose)
                sinf = (dsqrt(1.d0 - ecc(j)**2) * sine) / (1.d0 - ecc(j) * cose)

                fac1 = cosw * cosf - sinw * sinf + ecc(j) * cosw
                fac2 = cosw * sinf + sinw * cosf
                fac3 = -a2(1 + i) * dsqrt(1.d0 - ecc(j)**2) * fac2 / &
                        (1.d0 - ecc(j) * cose)**2

                y_pl(j) = a2(1 + i) * fac1
                y = y + a2(1 + i) * fac1
                dyda(1 + i) = fac1
                dyda(2 + i) = -TWOPI * fac3 * x / a2(2 + i)**2
                dyda(3 + i) = -a2(1 + i) * fac2 * ((2.d0 - ecc(j)**2 - ecc(j) * cose) * &
                        sinw * sine / dsqrt(1.d0 - ecc(j)**2) - &
                        dsqrt(1.d0 - ecc(j) * 2) * cosw / ecc(j)) / &
                        (1.d0 - ecc(j) * cose)**2 - &
                        a2(1 + i) * fac2 * cosw / ecc(j)
                dyda(4 + i) = -a2(1 + i) * fac2 * ((2.d0 - ecc(j)**2 - ecc(j) * cose) * &
                        cosw * sine / dsqrt(1.d0 - ecc(j)**2) + &
                        dsqrt(1.d0 - ecc(j) * 2) * sinw / ecc(j)) / &
                        (1.d0 - ecc(j) * cose)**2 + &
                        a2(1 + i) * fac2 * sinw / ecc(j)
                dyda(5 + i) = fac3
                dyda(6 + i) = dyda(4 + i) * x / 365.25d0
            else
                wm = TWOPI * x / a2(2 + i) + a2(5 + i)
                wm = dmod(wm, 2.d0 * PI)

                coswm = dcos(wm)
                sinwm = dsin(wm)
                cos2wm = dcos(2.d0 * wm)
                sin2wm = dsin(2.d0 * wm)
                cos3wm = dcos(3.d0 * wm)
                sin3wm = dsin(3.d0 * wm)

                fac1 = coswm + a2(3 + i) * sin2wm - a2(4 + i) * (1.d0 - cos2wm) - &
                        a2(3 + i)**2 * (0.875d0 * coswm + 1.125d0 * cos3wm) - &
                        a2(3 + i) * a2(4 + i) * (0.25d0 * sinwm - 2.25d0 * sin3wm) - &
                        a2(4 + i)**2 * 1.125d0 * (coswm - cos3wm)
                fac3 = -sinwm + 2.d0 * a2(3 * i) * cos2wm - 2.d0 * a2(4 + i) * sin2wm + &
                        a2(3 + i)**2 * (0.875d0 * sinwm - 3.375d0 * sin3wm) - &
                        a2(3 + i) * a2(4 + i) * (0.25d0 * coswm - 6.75d0 * cos3wm) + &
                        a2(4 + i)**2 * (1.125d0 * coswm - 3.375 * sin3wm)

                y_pl(j) = a2(1 + i) * fac1
                y = y + a2(1 + i) * fac1
                dyda(1 + i) = fac1
                dyda(2 + i) = -a2(1 + i) * TWOPI * fac3 * x / a2(2 + i)**2
                dyda(3 + i) = a2(1 + i) * (sin2wm - &
                        a2(3 + i) * (1.75d0 * coswm + 2.25d0 * cos3wm) - &
                        a2(4 + i) * (0.25d0 * sinwm - 2.25d0 * sin3wm))
                dyda(4 + i) = a2(1 + i) * (-1.d0 + a2(4 + i) * cos2wm - &
                        a2(3 + i)**(0.25d0 * sinwm - 2.25d0 * sin3wm) - &
                        a2(4 + i) * 2.25d0 * (coswm - cos3wm))
                dyda(5 + i) = a2(1 + i) * fac3
                dyda(6 + i) = dyda(4 + i) * x / 365.25d0
            endif
        enddo
    endif

    y = y + a2(6 * npl + ts)
    dyda(6 * npl + ts) = 1.d0

    y = y + a2(6 * npl + 2 * ndset + 1) * x
    y = y + a2(6 * npl + 2 * ndset + 2) + a2(6 * npl + 2 * ndset + 2) * x**2

    do i = ts + 1, ndset
        dyda(6 * npl + i) = 0.d0
    enddo

    return
end


subroutine prepare_for_amoeba_kep(p, mp, np, y, a, ia, ma, mfit, funk, &
        ndata, x, z, dyda, ts, sig, hkl)
    implicit none
    integer :: ndata, MMAX, NDSMAX, ma, ts(ndata), mp, np
    integer :: mfit,idset,ndset,npl
    parameter(MMAX = 200, NDSMAX = 20) 
    REAL(8) :: p(mp, np), y(mp), a(MMAX), a2(mfit), fr, frjitt
    real(8) :: x(ndata), z(ndata)
    real(8) :: dyda(MMAX), sig(ndata), loglik
    parameter(fr = 0.05, frjitt = 0.05)
    integer :: i, j, k, ia(MMAX), idsmax(NDSMAX), hkl, gr_flag
    external funk

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    k = 0
    do j = 1, ma
        if(ia(j).ne.0) then
            k = k + 1
            p(1, k) = a(j)
            do i = 2, mfit + 1
                if (k.eq.(i - 1)) then
                    if (j>(6 * npl + ndset)) then
                        p(i, k) = (1 + frjitt) * (p(1, k) + 0.1)
                    else
                        if (mod(j, 5).eq.2) then
                            p(i, k) = (1 + fr) * (p(1, k) + 0.1)
                        else if (mod(j, 5).eq.3) then
                            p(i, k) = (1 + frjitt) * (p(1, k) + 0.1)
                        else
                            p(i, k) = (1 + fr) * (p(1, k) + 0.1)
                        endif
                    endif
                else
                    p(i, k) = p(1, k)
                endif
            enddo
        endif
    enddo
    do i = 1, mfit + 1
        do j = 1, mfit
            a2(j) = p(i, j)
        enddo
        call funk(ndata, x, z, a2, dyda, ma, mfit, ts, sig, loglik, &
                a, ia, hkl)
        y(i) = loglik
    enddo
    return
end

SUBROUTINE amoeba_kep(p, y, mp, np, ndim, ftol, funk, iter, ndata, x, z, &
        dyda, ma, ts, sig, a, ia, ytry, hkl)
    implicit none
    integer :: ndata    
    integer :: iter, mp, ndim, np, NMAX, ITMAX, MMAX, ma, ts(ndata)
    REAL(8) :: ftol, p(mp, np), y(mp), x(ndata), z(ndata)
    PARAMETER (NMAX = 20, ITMAX = 200000, MMAX = 200)
    real(8) :: dyda(MMAX), sig(ndata), loglik, a(MMAX)
    EXTERNAL funk
    integer :: i, ihi, ilo, inhi, j, m, n, ia(MMAX), hkl
    REAL(8) :: rtol, summ, swap, ysave, ytry, psum(ndim), amotry_kep
    iter = 0
1   do n = 1, ndim
        summ = 0.d0
        do m = 1, ndim + 1
            summ = summ + p(m, n)
        enddo
        psum(n) = summ
    enddo
2   ilo = 1
    if (y(1)>y(2)) then
        ihi = 1
        inhi = 2
    else
        ihi = 2
        inhi = 1
    endif
    do i = 1, ndim + 1
        if(y(i)<=y(ilo)) ilo = i
        if(y(i)>y(ihi)) then
            inhi = ihi
            ihi = i
        else if(y(i)>y(inhi)) then
            if(i.ne.ihi) inhi = i
        endif
    enddo
    rtol = 2.d0 * abs(y(ihi) - y(ilo)) / (abs(y(ihi)) + abs(y(ilo)))
    if (rtol<ftol) then
        swap = y(1)
        y(1) = y(ilo)
        y(ilo) = swap
        do n = 1, ndim
            swap = p(1, n)
            p(1, n) = p(ilo, n)
            p(ilo, n) = swap
        enddo
        return
    endif
    if (iter>=ITMAX) then
        write (*, *) 'ITMAX exceeded in amoeba'
        return
    endif
    iter = iter + 2
    ytry = amotry_kep(p, y, psum, mp, np, ndim, funk, ihi, -1.0d0, ndata, x, z, &
            dyda, ma, ts, sig, a, ia, hkl)
    if (ytry<=y(ilo)) then
        ytry = amotry_kep(p, y, psum, mp, np, ndim, funk, ihi, 2.0d0, ndata, x, z, &
                dyda, ma, ts, sig, a, ia, hkl)
    else if (ytry>=y(inhi)) then
        ysave = y(ihi)
        ytry = amotry_kep(p, y, psum, mp, np, ndim, funk, ihi, 0.5d0, ndata, x, z, &
                dyda, ma, ts, sig, a, ia, hkl)

        if (ytry>=ysave) then
            do i = 1, ndim + 1
                if(i.ne.ilo)then
                    do j = 1, ndim
                        psum(j) = 0.5d0 * (p(i, j) + p(ilo, j))
                        p(i, j) = psum(j)
                    enddo
                    call funk(ndata, x, z, psum, dyda, ma, ndim, ts, sig, loglik, &
                            a, ia, hkl)
                    y(i) = loglik
                endif
            enddo
            iter = iter + ndim
            goto 1
        endif
    else
        iter = iter - 1
    endif
    goto 2
END
!(C) Copr. 1986-92 Numerical Recipes Software 0=M,173+9.

subroutine prepare_for_amoeba_dyn(p, mp, np, y, a, ia, ma, mfit, funk, &
        ndata, x, z, ts, sig, epsil, deltat, hkl, dynamical_planets, coplar_inc)
    implicit none        
    integer :: MMAX, NDSMAX, ma, ndset,  npl, ndata, idset
    integer :: mp, np, mfit, hkl, ts(ndata)
    parameter(MMAX = 200, NDSMAX = 20)
    REAL(8) :: p(mp, np), y(mp), a(MMAX), a2(mfit), fr, frjitt
    real(8) :: x(ndata), z(ndata)
    real(8) :: sig(ndata), loglik, epsil, deltat
    parameter(fr = 0.01, frjitt = 0.05)
    integer :: i, j, k, ia(MMAX), idsmax(NDSMAX), gr_flag, coplar_inc
    integer :: dynamical_planets(npl)
    external funk

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    k = 0
    do j = 1, ma
        if(ia(j).ne.0) then
            k = k + 1
            p(1, k) = a(j)
            do i = 2, mfit + 1
                if (hkl.eq.0) then
                    if (k.eq.(i - 1)) then
                        if (j>(7 * npl + ndset)) then
                            p(i, k) = (1 + frjitt) * (p(1, k) + 0.1)
                        else
                            if (mod(j, 7).eq.2) then
                                p(i, k) = (1 + fr) * (p(1, k) + 0.1)
                            else if (mod(j, 7).eq.3) then
                                p(i, k) = (1 + frjitt) * (p(1, k) + 0.1)
                            else if (mod(j, 7).eq.6) then
                                p(i, k) = (1 + 0.2) * (p(1, k) + 0.1)
                            else if (mod(j, 7).eq.7) then
                                p(i, k) = (1 + 0.4) * (p(1, k) + 0.1)
                            else
                                p(i, k) = (1 + fr) * (p(1, k) + 0.1)
                            endif
                        endif
                    else
                        p(i, k) = p(1, k)
                    endif
                else
                    if (k.eq.(i - 1)) then
                        if (j>(7 * npl + ndset)) then
                            p(i, k) = (1 + frjitt) * (p(1, k) + 0.1)
                        else
                            if (mod(j, 7).eq.2) then
                                p(i, k) = (1 + fr) * (p(1, k) + 0.1)
                            else if (mod(j, 7).eq.3) then
                                p(i, k) = (1 + frjitt) * (p(1, k) + 0.0001)
                            else if (mod(j, 7).eq.4) then
                                p(i, k) = (1 + frjitt) * (p(1, k) + 0.0001)
                            else if (mod(j, 7).eq.6) then
                                p(i, k) = (1 + 0.2) * (p(1, k) + 0.1)
                            else if (mod(j, 7).eq.7) then
                                p(i, k) = (1 + 0.4) * (p(1, k) + 0.1)
                            else
                                p(i, k) = (1 + fr) * (p(1, k) + 0.1)
                            endif
                        endif
                    else
                        p(i, k) = p(1, k)
                    endif
                endif
            enddo
        endif
    enddo

    do i = 1, mfit + 1
        do j = 1, mfit
            a2(j) = p(i, j)
        enddo
        call funk(ndata, x, z, a2, ma, mfit, ts, sig, loglik, &
                a, ia, epsil, deltat, hkl, dynamical_planets, coplar_inc)
        y(i) = loglik
    enddo

    return
end

SUBROUTINE amoeba_dyn(p, y, mp, np, ndim, ftol, funk, iter, ndata, x, z, &
        ma, ts, sig, a, ia, epsil, deltat, hkl, npl, dynamical_planets, coplar_inc)
    implicit none
    integer :: iter, mp, ndim, np, NMAX, ITMAX, MMAX, ma, ts(20000), ndata
    REAL(8) :: ftol, p(mp, np), y(mp), x(20000), z(20000)
    PARAMETER (NMAX = 20, ITMAX = 50000, MMAX = 200)
    real(8) :: sig(20000), loglik, a(MMAX), deltat
    EXTERNAL funk
    integer :: i, ihi, ilo, inhi, j, m, n, ia(MMAX), hkl, npl
    integer :: dynamical_planets(npl), coplar_inc
    REAL(8) :: rtol, summ, swap, ysave, ytry, psum(ndim), amotry_dyn, epsil
    iter = 0
11  do n = 1, ndim
        summ = 0.d0
        do m = 1, ndim + 1
            summ = summ + p(m, n)
        enddo
        psum(n) = summ
    enddo
22  ilo = 1
    if (y(1)>y(2)) then
        ihi = 1
        inhi = 2
    else
        ihi = 2
        inhi = 1
    endif
    do i = 1, ndim + 1
        if(y(i)<=y(ilo)) ilo = i
        if(y(i)>y(ihi)) then
            inhi = ihi
            ihi = i
        else if(y(i)>y(inhi)) then
            if(i.ne.ihi) inhi = i
        endif
    enddo
    rtol = 2.d0 * abs(y(ihi) - y(ilo)) / (abs(y(ihi)) + abs(y(ilo)))
    if (rtol<ftol) then
        swap = y(1)
        y(1) = y(ilo)
        y(ilo) = swap
        do n = 1, ndim
            swap = p(1, n)
            p(1, n) = p(ilo, n)
            p(ilo, n) = swap
        enddo
        return
    endif
    if (iter>=ITMAX) then
        write (*, *) 'ITMAX exceeded in amoeba'
        return
    endif
    
    
    
    
    iter = iter + 2
    ytry = amotry_dyn(p, y, psum, mp, np, ndim, funk, ihi, -1.0d0, ndata, x, z, &
            ma, ts, sig, a, ia, epsil, deltat, hkl, npl, dynamical_planets, coplar_inc)
    if (ytry<=y(ilo)) then
        ytry = amotry_dyn(p, y, psum, mp, np, ndim, funk, ihi, 2.0d0, ndata, x, z, &
                ma, ts, sig, a, ia, epsil, deltat, hkl, npl, dynamical_planets, coplar_inc)
    else if (ytry>=y(inhi)) then
        ysave = y(ihi)
        ytry = amotry_dyn(p, y, psum, mp, np, ndim, funk, ihi, 0.5d0, ndata, x, z, &
                ma, ts, sig, a, ia, epsil, deltat, hkl, npl, dynamical_planets, coplar_inc)
        if (ytry>=ysave) then
            do i = 1, ndim + 1
                if(i.ne.ilo)then
                    do j = 1, ndim
                        psum(j) = 0.5d0 * (p(i, j) + p(ilo, j))
                        p(i, j) = psum(j)
                    enddo
                    call funk(ndata, x, z, psum, ma, ndim, ts, sig, loglik, &
                            a, ia, epsil, deltat, hkl, dynamical_planets, coplar_inc)
                    y(i) = loglik
                endif
            enddo
            iter = iter + ndim
            goto 11
        endif
    else
        iter = iter - 1
    endif
    goto 22
END
!  (C) Copr. 1986-92 Numerical Recipes Software 0=M,173+9.

FUNCTION amotry_kep(p, y, psum, mp, np, ndim, funk, ihi, fac, ndata, x, z, &
        dyda, ma, ts, sig, a, ia, hkl)
    implicit none
    integer :: ihi, mp, ndim, np, NMAX, MMAX, ma, ts(20000), ndata
    PARAMETER (NMAX = 20, MMAX = 200)
    REAL(8) :: amotry_kep, fac, p(mp, np), psum(np), y(mp), x(20000), z(20000)
    real(8) :: dyda(MMAX), sig(20000), loglik
    EXTERNAL funk
    integer :: j, ia(MMAX), hkl
    REAL(8) :: fac1, fac2, ytry, ptry(ndim), a(MMAX)
    
    fac1 = (1.0d0 - fac) / ndim
    fac2 = fac1 - fac
    
    do j = 1, ndim
        ptry(j) = psum(j) * fac1 - p(ihi, j) * fac2
    enddo
    call funk(ndata, x, z, ptry, dyda, ma, ndim, ts, sig, loglik, &
            a, ia, hkl)

    ytry = loglik
    if (ytry<y(ihi)) then
        y(ihi) = ytry
        do j = 1, ndim
            psum(j) = psum(j) - p(ihi, j) + ptry(j)
            p(ihi, j) = ptry(j)
        enddo
    endif
    amotry_kep = ytry
    return
END
!  (C) Copr. 1986-92 Numerical Recipes Software 0=M,173+9.

FUNCTION amotry_dyn(p, y, psum, mp, np, ndim, funk, ihi, fac, ndata, x, z, &
        ma, ts, sig, a, ia, epsil, deltat, hkl, npl, dynamical_planets, coplar_inc)
    implicit none
    integer ::  ihi, mp, ndim, np, NMAX, MMAX, ma, ts(20000), ndata
    PARAMETER (NMAX = 20, MMAX = 200)
    REAL(8) :: amotry_dyn, fac, p(mp, np), psum(np), y(mp), x(20000), z(20000), &
            epsil, deltat
    real(8) :: sig(20000), loglik
    EXTERNAL funk
    integer ::  j, ia(MMAX), hkl, npl
    REAL(8) :: fac1, fac2, ytry, ptry(ndim), a(MMAX)
    integer :: dynamical_planets(npl), coplar_inc
    fac1 = (1.0d0 - fac) / ndim
    fac2 = fac1 - fac
    do j = 1, ndim
        ptry(j) = psum(j) * fac1 - p(ihi, j) * fac2
    enddo
    call funk(ndata, x, z, ptry, ma, ndim, ts, sig, loglik, &
            a, ia, epsil, deltat, hkl, dynamical_planets, coplar_inc)
           

    ytry = loglik
    if (ytry<y(ihi)) then
        y(ihi) = ytry
        do j = 1, ndim
            psum(j) = psum(j) - p(ihi, j) + ptry(j)
            p(ihi, j) = ptry(j)
        enddo
    endif
    amotry_dyn = ytry
    return
END

subroutine transit_tperi(per, ecc, om, ma, epoch, t_peri, t_transit)
    implicit none
    real(8) :: per, om, ma, E, ecc, epoch, t_peri, t_transit
    real(8) :: PI, TWOPI, per_days
    parameter (PI = 3.14159265358979d0)
    parameter (TWOPI = 2.0 * PI)

    per_days = per / 8.64d4

    write(*,*) per_days, om, ma, ecc, epoch
    E = 2.0 * datan(dsqrt(((1.0 - ecc) / (1.0 + ecc))) * &
            dtan((PI / 4.0) - (om / 2.0)))

    write(*,*) E
    t_peri = epoch - ((ma / TWOPI) * per_days)
    t_transit = t_peri + (E + ecc * dsin(E)) * (per_days / TWOPI)

    return
end

subroutine MA_J_kepamo (a, ma, npl, m0, mass, ap, hkl, gr_flag)
    implicit none
    real(8) :: m0, PI, TWOPI, THIRD, GMSUN, dm, MSUN
    integer :: npl, ma, i, j, NPLMAX, hkl, gr_flag
    parameter (NPLMAX = 10)
    real(8) :: mm(NPLMAX), ecc, corr
    real(8) :: a(ma), mass(NPLMAX), ap(NPLMAX), mpold(NPLMAX), mtotal
    parameter (THIRD = 1.d0 / 3.d0)
    parameter (PI = 3.14159265358979d0, TWOPI = 2.d0 * PI)
    parameter (GMSUN = 1.32712440018d20, MSUN = 1.32712440018d20)

    !G is set to be unit, and s, m, kg as unit of time, length and mass expectively.

    do j = 1, npl
        i = 6 * (j - 1)
        mm(j) = 2.d0 * PI / (a(i + 2) * 8.64d4)
    enddo

    do i = 0, npl - 1
        if (hkl.eq.0) then
            ecc = a(6 * i + 3)
        else
            ecc = dsqrt(a(6 * i + 3)**2 + a(6 * i + 4)**2)    !! only for h, k
        endif

        mass(1) = m0
        mpold(i + 1) = 0.d0
        dm = 10.d0
        do while (dm>0)
            if (i.eq.0) then
                mtotal = m0
                mass(i + 2) = a(6 * i + 1) * (TWOPI / mm(i + 1) * (m0 + mpold(i + 1))**2 / &
                        (TWOPI * GMSUN))**THIRD * &
                        dsqrt(1.d0 - ecc**2)
            else
                mtotal = m0
                do j = 0, i - 1
                    mtotal = mtotal + mass(j + 2)
                enddo
                mass(i + 2) = a(6 * i + 1) * (TWOPI / mm(i + 1) * (mtotal&
                        + mpold(i + 1))**2 / (TWOPI * GMSUN))**THIRD * &
                        dsqrt(1.d0 - ecc**2)
            endif

            dm = dabs(mass(i + 2) - mpold(i + 1)) / mass(i + 2)
            mpold(i + 1) = mass(i + 2)
        enddo

        ap(i + 1) = (GMSUN * (mtotal + mass(i + 2)) * (1.d0 / mm(i + 1))&
                **2)**THIRD
    enddo

    do i = 1, npl + 1
        mass(i) = mass(i) * MSUN
    enddo

    if(gr_flag.ne.0) then
        do i = 1, npl
            j = 6 * (i - 1)
            call gr_corr(ap(i), a(j + 3), mass(1), corr)
            a(j + 6) = corr * 365.25
        enddo
    endif
    return
end

subroutine gr_corr(a, e, gmi, corr)
    implicit none
    real(8) :: a, e, T, PI, c, GMSUN, AU
    real(8) :: corr, THIRD, gmi
    parameter (PI = 3.14159265358979d0)
    parameter (c = 299792458.0d0)
    parameter (GMSUN = 1.32712440018d20, AU = 1.49597892d13)
    parameter (THIRD = 1.d0 / 3.d0)

    T = 2.0d0 * PI * sqrt((a**3.0d0) / (gmi))
    corr = (24.0d0 * (PI**3.0d0) * (a**2.0d0)) / ((T**2.0d0) * &
            (c**2.0d0) * (1.0d0 - e**2.0d0))
    corr = corr / T

    return
end

subroutine io_write_bf_ewcop_fin_dynamo (a, covar, t, ys, &
        ndata, ts, ma, mfit, t0, t_max, sigs, chisq, rms, loglik, &
        writeflag_RV, writeflag_best_par, writeflag_fit, epsil, &
        deltat, nt, model_max, model_min, hkl, wdot, u_wdot, &
        dynamical_planets, &
        res_array, fit_return, &
        bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
        fit_array, coplar_inc)
    implicit none
    real(8) :: PI
    integer :: MMAX, NDSMAX, NPLMAX,ndata
    parameter (PI = 3.14159265358979d0, MMAX = 200, NDSMAX = 20, NPLMAX = 10)
    real(8) :: a(MMAX), t(20000), ymod(20000), ys(20000)
    real(8) :: covar(MMAX, MMAX), AU, day
    real(8) :: rms, mstar, sini(NPLMAX), mass(NPLMAX), ap(NPLMAX)
    integer :: ts(20000), nbod, nt, writeflag_RV, gr_flag
    integer :: writeflag_best_par, writeflag_fit, hkl
    real(8) :: t0, dt, t_max, chisq, loglik, dy, sig2i, twopi
    real(8) :: x(20000), sigs(20000)
    integer :: i, j, npl, ndset,  idset, mfit, ma, idsmax(NDSMAX)
    real(8) :: xj(NPLMAX), yj(NPLMAX),  zj(NPLMAX) 
    real(8) :: vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
    real(8) :: rpl(NPLMAX), rhill(NPLMAX), deltat, epsil
    real(8) :: j_mass(NPLMAX)
    real(4) :: model_max, model_min
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX), best_w, best_we
    integer :: dynamical_planets(npl), coplar_inc
    real(8) :: ymod_pl(npl,20000)     
    
    parameter (AU = 1.49597892d11, day = 86400.d0)

    real(8) :: res_array(ndata, 6+npl)
    real(8) :: fit_return(4), fit_array(nt, 2+npl)
    real(8) :: bestpar_1(npl, 17, 2), bestpar_2(ndset, 2)
    real(8) :: bestpar_3(ndset, 2), bestpar_4(9 + 2 * npl)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    common mstar, sini

    nbod = npl + 1
    rms = 0.d0
    twopi = 2.0d0 * PI


    if (coplar_inc.ne.0) then
        do i = 2, npl
             j = 7 * (i - 1)
             a(j + 6) = a(6)
        enddo
    endif

    call RVKEP_dynamo(t, a, ymod, ymod_pl, ma, ndata, epsil, deltat, hkl, &
            dynamical_planets, ts, coplar_inc)

    do i = 1, npl
        j = 7 * (i - 1)

        if (hkl.eq.0) then
            !a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)
            if (a(j + 2)<0.d0) then  !if P<0, set P>0
                a(j + 2) = dabs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then  !if K<0, set K>0 and w = w+PI
                a(j + 4) = a(j + 4) + PI
                a(j + 1) = dabs(a(j + 1))
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
            endif
            if (a(j + 3)<0.d0) then  !if e<0, set e>0 and w=w+PI, M0=M0-PI
                a(j + 3) = dabs(a(j + 3))
                a(j + 4) = a(j + 4) + PI
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
                a(j + 5) = a(j + 5) - PI
                if (a(j + 5)<0.d0) a(j + 5) = a(j + 5) + 2.d0 * PI
            endif
            if(a(j + 3)>=1.d0) then !if e>=1 set it to 0.99 to prevent errors
                a(j + 3) = 0.99d0
            endif
            if (a(j + 4)<0.d0) a(j + 4) = dmod(a(j + 4) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 4)>2.d0 * PI) a(j + 4) = dmod(a(j + 4), 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            if (a(j + 6)<0.d0) a(j + 6) = dmod(a(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 7)<0.d0) a(j + 7) = dmod(a(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 6)>2.d0 * PI) a(j + 6) = dmod(a(j + 6), 2.d0 * PI)
            if (a(j + 7)>2.d0 * PI) a(j + 7) = dmod(a(j + 7), 2.d0 * PI)
        else
            if (a(j + 2)<0.d0) then !if P<0, set P>0
                a(j + 2) = dabs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a(j + 4) = -1.d0 * a(j + 4) !which is h = -h, k = -k
                a(j + 3) = -1.d0 * a(j + 3)
                a(j + 1) = dabs(a(j + 1))
            endif

            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 6)<0.d0) a(j + 6) = dmod(a(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 7)<0.d0) a(j + 7) = dmod(a(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            if (a(j + 6)>2.d0 * PI) a(j + 6) = dmod(a(j + 6), 2.d0 * PI)
            if (a(j + 7)>2.d0 * PI) a(j + 7) = dmod(a(j + 7), 2.d0 * PI)
        endif
    enddo

    do i = 1, ndata
        idset = ts(i)

        ys(i) = ys(i) - a(7 * npl + idset)&
                - a(7 * npl + 2 * ndset + 1) * (t(i) / 8.64d4)&
                - a(7 * npl + 2 * ndset + 2) * (t(i) / 8.64d4)**2

        if (writeflag_RV>0) then
 
            res_array(i, :) = (/ t0 + t(i) / 8.64d4,  ys(i)&
                    + a(7 * npl + 2 * ndset + 1) * (t(i) / 8.64d4)&
                    + a(7 * npl + 2 * ndset + 2) * (t(i) / 8.64d4)**2, &
                    sigs(i), dble(ts(i)), ys(i) - ymod(i), ymod(i), &
                            (ymod_pl(j,i), j=1,npl)/)                              
                            
        endif

        sig2i = 1.d0 / (sigs(i)**2 + a(7 * npl + ndset + idset)**2)
        dy = ys(i) - ymod(i)
        chisq = chisq + dy * dy * sig2i
        loglik = loglik - 0.5 * dy * dy * sig2i - &
                0.5 * dlog(twopi * (sigs(i)**2&
                        + a(7 * npl + ndset + idset)**2))
        rms = rms + dy**2
    enddo

    rms = dsqrt(rms / dble(ndata))
    fit_return = (/ loglik, chisq / dble(ndata - mfit), chisq, rms /)

    if(writeflag_best_par>0) then
        do j = 1, npl
            i = 7 * (j - 1)

            if (hkl.eq.0) then
                best_w = a(i + 4) * 180.d0 / PI
                best_we = dsqrt(covar(i + 4, i + 4)) * 180.d0 / PI
            else
                best_w = a(i + 4)
                best_we = dsqrt(covar(i + 4, i + 4))
            endif

!            write(*,*)"incl1", dmod(a(i + 6) * 180.d0 / PI, 180.d0)
            bestpar_1(j, :, 1) = (/ a(i + 1), a(i + 2) / 8.64d4, a(i + 3), &
                    best_w, a(i + 5) * 180.d0 / PI, &
                    dmod(a(i + 6) * 180.d0 / PI, 180.d0), &
                    dmod(a(i + 7) * 180.d0 / PI, 360.d0), wdot(i), &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
            bestpar_1(j, :, 2) = (/ dsqrt(covar(i + 1, i + 1)), &
                    dsqrt(covar(i + 2, i + 2)), &
                    dsqrt(covar(i + 3, i + 3)), best_we, &
                    dsqrt(covar(i + 5, i + 5)) * 180.d0 / PI, &
                    dmod(dsqrt(covar(i + 6, i + 6)) * 180.d0 / PI, 180.d0), &
                    dmod(dsqrt(covar(i + 7, i + 7)) * 180.d0 / PI, 360.d0), &
                    u_wdot(i), &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
        enddo

        do j = 1, ndset
            i = 7 * npl + j
            bestpar_2(j, :) = (/ a(i), dsqrt(covar(i, i)) /)

            bestpar_3(j, :) = (/ a(7 * npl + ndset + j), 0.d0 /)
        enddo

        call MA_J_cop_fin (a, ma, npl, mstar, mass, ap, hkl)

        call GENINIT_J3_ewcop (nbod, ap, a, &
                mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)

        do j = 1, npl + 1
            j_mass(j) = mass(j) / 1.26686534d17
        enddo

        bestpar_4 = (/ a(7 * npl + 2 * ndset + 1), &
                dsqrt(covar(7 * npl + 2 * ndset + 1, 7 * npl + 2 * ndset + 1)), &
                a(7 * npl + 2 * ndset + 2), &
                dsqrt(covar(7 * npl + 2 * ndset + 2, 7 * npl + 2 * ndset + 2)), &
                dble(ndata), dble(mfit), rms, &
                chisq / dble(ndata - mfit), t0, (j_mass(i + 1), i = 1, npl), &
                (ap(i) / 1.49597892d11, i = 1, npl) /)
    endif

    if(writeflag_fit>0) then
        dt = ((t_max - t0) + model_max + model_min) / dble(nt - 1)

        do i = 1, nt
            x(i) = (i - 1) * dt * 8.64d4
        enddo
        do i = 1, npl
            j = 7 * (i - 1)
!            write(*,*)"incl", j, a(j+6)
        enddo
        call RVKEP_dynamo (x, a, ymod, ymod_pl, ma, nt, epsil, deltat, hkl, &
                dynamical_planets, ts, coplar_inc)
                
        do i = 1, nt
!            write(*,*)(ymod_pl(j,i),j=1,npl)                         
            fit_array(i, :) = (/ t0 + x(i) / 8.64d4, &
                    ymod(i) + a(7 * npl + 2 * ndset + 1) * (x(i) / 8.64d4)&
                            + a(7 * npl + 2 * ndset + 2) * (x(i) / 8.64d4)**2 &
                            ,(ymod_pl(j,i), j=1,npl)/)
        enddo
    endif

    return
end

subroutine io_write_bf_ewcop_fin_dynlm (a, covar, t, ys, ndata, ts, &
             ma, mfit, t0, t_max, sigs, chisq, rms, writeflag_RV, &
             writeflag_best_par, writeflag_fit, jitter, epsil, deltat, &
             nt, model_max, model_min, hkl, wdot, u_wdot, &
             res_array, fit_return, &
             bestpar_1, bestpar_2, bestpar_3, bestpar_4, &
             fit_array, coplar_inc)
    implicit none
    real(8) :: PI
    integer :: MMAX, NDSMAX, NPLMAX
    parameter (PI = 3.14159265358979d0, MMAX = 200, NDSMAX = 20, NPLMAX = 10)
    real(8) :: a(MMAX), t(20000), ymod(20000), ys(20000)
    real(8) :: covar(MMAX, MMAX), dyda(20000, MMAX), AU, day, loglik, dy, twopi
    real(8) :: rms, mstar, sini, mass(NPLMAX), ap(NPLMAX)
    integer :: ts(20000), nbod, nt, writeflag_RV, gr_flag
    integer :: writeflag_best_par, writeflag_fit, coplar_inc
    real(8) :: t0, dt, t_max, chisq, deltat, epsil, sig2i
    real(8) :: x(20000), sigs(20000), jitter(NDSMAX)
    integer :: i, j, npl, ndset, ndata, idset, mfit, ma, idsmax(NDSMAX), hkl
    real(8) :: j_mass(NPLMAX)
    real(4) :: model_max, model_min
    parameter (AU = 1.49597892d11, day = 86400.d0)
    real(8) :: wdot(NPLMAX), u_wdot(NPLMAX), best_w, best_we
    real(8) :: ymod_pl(npl, 20000)     
    real(8) :: res_array(ndata, 6 + npl)
    real(8) :: fit_return(4), fit_array(nt, 2 + npl)
    real(8) :: bestpar_1(npl, 17, 2), bestpar_2(ndset, 2)
    real(8) :: bestpar_3(ndset, 2), bestpar_4(9 + 2 * npl)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    common mstar, sini
    save dyda

    nbod = npl + 1
    rms = 0.d0
    twopi = 2.0d0 * PI 
    chisq = 0
    loglik = 0.d0

    call RVKEP_ewcop_fin (t, a, ymod, ymod_pl, dyda, ma, ndata, epsil, deltat, hkl, coplar_inc)
    call MA_J_cop_fin (a, ma, npl, mstar, mass, ap, hkl)

    do i = 1, npl
        j = 7 * (i - 1)

        if (hkl.eq.0) then
            !a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)
            if (a(j + 2)<0.d0) then  ! if P<0, set P>0
                a(j + 2) = abs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a(j + 4) = a(j + 4) + PI
                a(j + 1) = abs(a(j + 1))
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
            endif
            if (a(j + 3)<0.d0) then !if e<0, set e>0 and w=w+PI, M0=M0-PI
                a(j + 3) = abs(a(j + 3))
                a(j + 4) = a(j + 4) + PI
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
                a(j + 5) = a(j + 5) - PI
                if (a(j + 5)<0.d0) a(j + 5) = a(j + 5) + 2.d0 * PI
            endif
            if(a(j + 3)>=1.d0) then !if e>=1 set it to 0.99 to prevent errors
                a(j + 3) = 0.99d0
            endif
            if (a(j + 4)<0.d0) a(j + 4) = dmod(a(j + 4) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 4)>2.d0 * PI) a(j + 4) = dmod(a(j + 4), 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            if (a(j + 6)<0.d0) a(j + 6) = dmod(a(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 7)<0.d0) a(j + 7) = dmod(a(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 6)>2.d0 * PI) a(j + 6) = dmod(a(j + 6), 2.d0 * PI)
            if (a(j + 7)>2.d0 * PI) a(j + 7) = dmod(a(j + 7), 2.d0 * PI)
        else
            if (a(j + 2)<0.d0) then  ! if P<0, set P>0
                a(j + 2) = abs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a(j + 4) = -1.d0 * a(j + 4) !which is h = -h, k = -k
                a(j + 3) = -1.d0 * a(j + 3)
                a(j + 1) = abs(a(j + 1))
            endif

            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 6)<0.d0) a(j + 6) = dmod(a(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 7)<0.d0) a(j + 7) = dmod(a(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            if (a(j + 6)>2.d0 * PI) a(j + 6) = dmod(a(j + 6), 2.d0 * PI)
            if (a(j + 7)>2.d0 * PI) a(j + 7) = dmod(a(j + 7), 2.d0 * PI)
        endif
    enddo

    do i = 1, ndata
        idset = ts(i)

        ys(i) = ys(i) - a(7 * npl + idset)&
                - a(7 * npl + ndset + 1) * (t(i) / 86400.d0)&
                - a(7 * npl + ndset + 2) * (t(i) / 86400.d0)**2

        if (writeflag_RV>0) then

            res_array(i, :) = (/ t0 + t(i) / 8.64d4,  ys(i)&
                    + a(7 * npl + ndset + 1) * (t(i) / 8.64d4)&
                    + a(7 * npl + ndset + 2) * (t(i) / 8.64d4)**2, &
                     sigs(i), dble(ts(i)), ys(i) - ymod(i), ymod(i), &
                            (ymod_pl(j,i), j=1,npl)/)                              
                            
        endif

        sig2i = 1.d0 / (sigs(i)**2 + jitter(idset)**2)
        dy = ys(i) - ymod(i)
        chisq = chisq + dy * dy * sig2i
        loglik = loglik - 0.5 * dy * dy * sig2i - &
                0.5 * dlog(twopi * (sigs(i)**2&
                        + jitter(idset)**2))
        rms = rms + dy**2
    enddo

    rms = dsqrt(rms / dble(ndata))
    fit_return = (/ loglik, chisq / dble(ndata - mfit), chisq, rms /)

    if(writeflag_best_par>0) then
        do j = 1, npl
            i = 7 * (j - 1)

            if (hkl.eq.0) then
                best_w = a(i + 4) * 180.d0 / PI
                best_we = dsqrt(covar(i + 4, i + 4)) * 180.d0 / PI
            else
                best_w = a(i + 4)
                best_we = dsqrt(covar(i + 4, i + 4))
            endif

            bestpar_1(j, :, 1) = (/ a(i + 1), a(i + 2) / 8.64d4, a(i + 3), &
                    best_w, a(i + 5) * 180.d0 / PI, &
                    dmod(a(i + 6) * 180.d0 / PI, 180.d0), &
                    dmod(a(i + 7) * 180.d0 / PI, 360.d0), wdot(i), &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
            bestpar_1(j, :, 2) = (/ dsqrt(covar(i + 1, i + 1)), &
                    dsqrt(covar(i + 2, i + 2)) / 8.64d4, &
                    dsqrt(covar(i + 3, i + 3)), best_we, &
                    dsqrt(covar(i + 5, i + 5)) * 180.d0 / PI, &
                    dmod(dsqrt(covar(i + 6, i + 6)) * 180.d0 / PI, 180.d0), &
                    dmod(dsqrt(covar(i + 7, i + 7)) * 180.d0 / PI, 360.d0), &
                    u_wdot(i), &
                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
        enddo

        do j = 1, ndset
            i = 7 * npl + j
            bestpar_2(j, :) = (/ a(i), dsqrt(covar(i, i)) /)

            bestpar_3(j, :) = (/ jitter(j), 0.d0 /)
        enddo

        do j = 1, npl + 1
            j_mass(j) = mass(j) / 1.26686534d17
        enddo

        bestpar_4 = (/ a(7 * npl + ndset + 1), &
                dsqrt(covar(7 * npl + ndset + 1, 7 * npl + ndset + 1)), &
                a(7 * npl + ndset + 2), &
                dsqrt(covar(7 * npl + ndset + 2, 7 * npl + ndset + 2)), &
                dble(ndata), dble(mfit), rms, &
                chisq / dble(ndata - mfit), t0, (j_mass(i + 1), i = 1, npl), &
                (ap(i) / 1.49597892d11, i = 1, npl) /)
    endif

    if(writeflag_fit>0) then
        dt = ((t_max - t0) + model_max + model_min) / dble(nt - 1)
        do i = 1, nt
            x(i) = (i - 1) * dt * 8.64d4
        enddo
        call RVKEP_ewcop_fin (x, a, ymod, ymod_pl, dyda, ma, nt, epsil, deltat,&
         hkl, coplar_inc)
        do i = 1, nt
!            write(*,*) i,ymod(i), a(7 * npl * ndset + 1), & 
!                                   a(7 * npl * ndset + 2),(ymod_pl(j,i),j=1,npl)                 

            fit_array(i, :) = (/ t0 + x(i) / 8.64d4, &
                    ymod(i) + a(7 * npl + ndset + 1) * (x(i) / 8.64d4)&
                            + a(7 * npl + ndset + 2) * (x(i) / 8.64d4)**2, &
                            (ymod_pl(j,i), j=1,npl)/)

        enddo
    endif
    return
end

subroutine RVKEP_keplm (x, a, y, y_pl, dyda, ma, ts, hkl)
    implicit none
    real(8) :: PI, TWOPI
    parameter (PI = 3.14159265358979d0)
    parameter (TWOPI = 2.0d0 * PI)
    integer :: npl, ndset, idset, ma, i, j, NDSMAX, ts
    integer :: hkl, gr_flag
    parameter (NDSMAX = 20)
    integer :: idsmax(NDSMAX)
    real(8) :: x, y, a(ma), dyda(ma), mass(10), ap(10), y_pl(10)
    real(8) :: cosw, sinw, capm, cape, cose, sine, cosf
    real(8) :: sinf, fac1, fac2, fac3
    real(8) :: orbel_ehybrid, omega(10), capmm(10), ecc(10)
    real(8) :: wm, sinwm, coswm, sin2wm, cos2wm, sin3wm, cos3wm, omegad(10)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    y = 0.d0

    if (hkl.eq.0) then

        do i = 1, npl
            j = 6 * (i - 1)

            if (a(j + 2)<0.d0) then !if P<0, set P>0
                a(j + 2) = abs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a(j + 4) = a(j + 4) + PI
                a(j + 1) = abs(a(j + 1))
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
            endif
            if (a(j + 3)<0.d0) then !if e<0, set e>0 and w=w+PI, M0=M0-PI
                a(j + 3) = abs(a(j + 3))
                a(j + 4) = a(j + 4) + PI
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
                a(j + 5) = a(j + 5) - PI
                if (a(j + 5)<0.d0) a(j + 5) = a(j + 5) + 2.d0 * PI
            endif
            if (a(j + 4)<0.d0) a(j + 4) = dmod(a(j + 4) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            !if (a(j+6)<0.d0) a(j+6) = dmod(a(j+6)+2.d0*PI, 2.d0*PI)

            if (a(j + 4)>2.d0 * PI) a(j + 4) = dmod(a(j + 4), 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            !if (a(j+6)>2.d0*PI) a(j+6) = dmod(a(j+6), 2.d0*PI)

            ecc(i) = a(j + 3)
            omega(i) = a(j + 4)
            capmm(i) = a(j + 5)

            if(gr_flag.ne.0) call MA_J_keplm (a, ma, npl, 1.0d0, &
                    mass, ap, hkl, gr_flag)

            omegad(i) = a(j + 6)
        enddo
    else
        do i = 1, npl
            j = 6 * (i - 1)
            if (a(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a(j + 4) = -1.d0 * a(j + 4) !which is h = -h, k = -k
                a(j + 3) = -1.d0 * a(j + 3)
                a(j + 1) = abs(a(j + 1))
            endif

            ecc(i) = dsqrt(a(j + 3)**2 + a(j + 4)**2)
            omega(i) = atan2(a(j + 3), a(j + 4))

            if(omega(i)<0.d0)omega(i) = dmod(omega(i) + 2.d0 * PI, 2.d0 * PI)
            if(omega(i)>0.d0)omega(i) = dmod(omega(i), 2.d0 * PI)
            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)

            capmm(i) = a(j + 5) - omega(i)

            if(capmm(i)<0.d0)capmm(i) = dmod(capmm(i) + 2.d0 * PI, 2.d0 * PI)
            if(capmm(i)>0.d0)capmm(i) = dmod(capmm(i), 2.d0 * PI)
        enddo
    endif
    if (hkl.eq.0) then
        do j = 1, npl

            i = 6 * (j - 1)
            cosw = dcos(omega(j) + omegad(j) * x / 365.25d0)
            sinw = dsin(omega(j) + omegad(j) * x / 365.25d0)
            capm = TWOPI * x / a(2 + i) + capmm(j)
            capm = dmod(capm, 2.d0 * PI)

            cape = ORBEL_EHYBRID (ecc(j), capm)

            cose = dcos(cape)
            sine = dsin(cape)

            cosf = (cose - ecc(j)) / (1.d0 - ecc(j) * cose)
            sinf = (dsqrt(1.d0 - ecc(j)**2) * sine) / (1.d0 - ecc(j) * cose)

            fac1 = cosw * cosf - sinw * sinf + ecc(j) * cosw

            fac2 = (cosw * sinf + sinw * cosf) / (1.d0 - ecc(j) * cose)**2
            fac3 = -a(1 + i) * dsqrt(1.d0 - ecc(j)**2) * fac2

            y_pl(j) = a(1 + i) * fac1
            y = y + a(1 + i) * fac1

            dyda(1 + i) = fac1
            dyda(2 + i) = -TWOPI * fac3 * x / a(2 + i)**2
            dyda(3 + i) = -a(1 + i) * sine * (2.d0 - ecc(j)**2 - ecc(j) * cose) * &
                    fac2 / dsqrt(1.d0 - ecc(j)**2)
            dyda(4 + i) = -a(1 + i) * (sinw * cosf + cosw * sinf + ecc(j) * sinw)
            dyda(5 + i) = fac3
            dyda(6 + i) = -a(1 + i) * (sinw * cosf + cosw * sinf&
                    + ecc(j) * sinw) * x / 365.25d0
        enddo
    else
        do j = 1, npl
            i = 6 * (j - 1)

            if (ecc(j)>1.d-2) then
                cosw = dcos(omega(j))
                sinw = dsin(omega(j))

                capm = TWOPI * x / a(2 + i) + capmm(j)
                capm = dmod(capm, 2.d0 * PI)
                !             write(*,*) capm
                cape = ORBEL_EHYBRID (ecc(j), capm)
                cose = dcos(cape)
                sine = dsin(cape)
                cosf = (cose - ecc(j)) / (1.d0 - ecc(j) * cose)
                sinf = (dsqrt(1.d0 - ecc(j)**2) * sine) / (1.d0 - ecc(j) * cose)

                fac1 = cosw * cosf - sinw * sinf + ecc(j) * cosw
                fac2 = cosw * sinf + sinw * cosf
                fac3 = -a(1 + i) * dsqrt(1.d0 - ecc(j)**2) * fac2 / &
                        (1.d0 - ecc(j) * cose)**2

                y_pl(j) = a(1 + i) * fac1
                y = y + a(1 + i) * fac1
                dyda(1 + i) = fac1
                dyda(2 + i) = -TWOPI * fac3 * x / a(2 + i)**2
                dyda(3 + i) = -a(1 + i) * fac2 * ((2.d0 - ecc(j)**2 - ecc(j) * cose) * &
                        sinw * sine / dsqrt(1.d0 - ecc(j)**2) - &
                        dsqrt(1.d0 - ecc(j) * 2) * cosw / ecc(j)) / &
                        (1.d0 - ecc(j) * cose)**2 - &
                        a(1 + i) * fac2 * cosw / ecc(j)
                dyda(4 + i) = -a(1 + i) * fac2 * ((2.d0 - ecc(j)**2 - ecc(j) * cose) * &
                        cosw * sine / dsqrt(1.d0 - ecc(j)**2) + &
                        dsqrt(1.d0 - ecc(j) * 2) * sinw / ecc(j)) / &
                        (1.d0 - ecc(j) * cose)**2 + &
                        a(1 + i) * fac2 * sinw / ecc(j)
                dyda(5 + i) = fac3
                dyda(6 + i) = dyda(4 + i) * x / 365.25d0
            else
                wm = TWOPI * x / a(2 + i) + a(5 + i)
                wm = dmod(wm, 2.d0 * PI)

                coswm = dcos(wm)
                sinwm = dsin(wm)
                cos2wm = dcos(2.d0 * wm)
                sin2wm = dsin(2.d0 * wm)
                cos3wm = dcos(3.d0 * wm)
                sin3wm = dsin(3.d0 * wm)

                fac1 = coswm + a(3 + i) * sin2wm - a(4 + i) * (1.d0 - cos2wm) - &
                        a(3 + i)**2 * (0.875d0 * coswm + 1.125d0 * cos3wm) - &
                        a(3 + i) * a(4 + i) * (0.25d0 * sinwm - 2.25d0 * sin3wm) - &
                        a(4 + i)**2 * 1.125d0 * (coswm - cos3wm)
                fac3 = -sinwm + 2.d0 * a(3 * i) * cos2wm - 2.d0 * a(4 + i) * sin2wm + &
                        a(3 + i)**2 * (0.875d0 * sinwm - 3.375d0 * sin3wm) - &
                        a(3 + i) * a(4 + i) * (0.25d0 * coswm - 6.75d0 * cos3wm) + &
                        a(4 + i)**2 * (1.125d0 * coswm - 3.375 * sin3wm)

                y_pl(j) = a(1 + i) * fac1
                y = y + a(1 + i) * fac1
                dyda(1 + i) = fac1
                dyda(2 + i) = -a(1 + i) * TWOPI * fac3 * x / a(2 + i)**2
                dyda(3 + i) = a(1 + i) * (sin2wm - &
                        a(3 + i) * (1.75d0 * coswm + 2.25d0 * cos3wm) - &
                        a(4 + i) * (0.25d0 * sinwm - 2.25d0 * sin3wm))
                dyda(4 + i) = a(1 + i) * (-1.d0 + a(4 + i) * cos2wm - &
                        a(3 + i)**(0.25d0 * sinwm - 2.25d0 * sin3wm) - &
                        a(4 + i) * 2.25d0 * (coswm - cos3wm))
                dyda(5 + i) = a(1 + i) * fac3
                dyda(6 + i) = dyda(4 + i) * x / 365.25d0
            endif
        enddo
    endif

    y = y + a(6 * npl + ts)
    dyda(6 * npl + ts) = 1.d0

    y = y + a(6 * npl + ndset + 1) * x + a(6 * npl + ndset + 2) * x**2

    dyda(6 * npl + ndset + 1) = x
    dyda(6 * npl + ndset + 2) = x**2

    do i = ts + 1, ndset
        dyda(6 * npl + i) = 0.d0
    enddo
    return
end

subroutine RVKEP_dynamoplus (x, a, y, y_pl, dyda, ma, ts, k)
    implicit none
    real(8) :: PI, TWOPI
    parameter (PI = 3.14159265358979d0)
    parameter (TWOPI = 2.0d0 * PI)
    integer :: npl, ndset, idset, ma, i, j, NDSMAX, ts, k
    parameter (NDSMAX = 20)
    integer :: idsmax(NDSMAX), gr_flag
    real(8) :: x, y, a(ma), dyda(ma), y_pl(npl)
    real(8) :: cosw, sinw, capm, cape, cose, sine, cosf, sinf, fac1, fac2, fac3
    real(8) :: orbel_ehybrid

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    do j = 1, k
        i = 7 * (j - 1)
        cosw = dcos(a(4 + i))
        sinw = dsin(a(4 + i))

        capm = TWOPI * x / (a(2 + i) * 86400.d0) + a(5 + i)
        capm = dmod(capm, 2.d0 * PI)

        cape = ORBEL_EHYBRID (a(3 + i), capm)
        cose = dcos(cape)
        sine = dsin(cape)

        cosf = (cose - a(3 + i)) / (1.d0 - a(3 + i) * cose)
        sinf = (dsqrt(1.d0 - a(3 + i)**2) * sine) / (1.d0 - a(3 + i) * cose)

        fac1 = cosw * cosf - sinw * sinf + a(3 + i) * cosw

        fac2 = (cosw * sinf + sinw * cosf) / (1.d0 - a(3 + i) * cose)**2
        fac3 = -a(1 + i) * dsqrt(1.d0 - a(3 + i)**2) * fac2

        y_pl(k) = a(1 + i) * fac1
        y = y + a(1 + i) * fac1
        dyda(1 + i) = fac1
        dyda(2 + i) = -TWOPI * fac3 * x / (a(2 + i) * 86400.d0)**2
        dyda(3 + i) = -a(1 + i) * sine * (2.d0 - a(3 + i)**2 - &
                a(3 + i) * cose) * fac2 / dsqrt(1.d0 - a(3 + i)**2)
        dyda(4 + i) = -a(1 + i) * (sinw * cosf + cosw * sinf)
        dyda(5 + i) = fac3
    enddo

    do i = 1, ts - 1
        dyda(7 * npl + i) = 0.d0
    enddo

!    write(*,*) ts, npl
    dyda(7 * npl + ts) = 1.d0

    dyda(7 * npl + ndset + 1) = x
    dyda(7 * npl + ndset + 2) = x * x

    do i = ts + 1, ndset
        dyda(7 * npl + i) = 0.d0
    enddo
    return
end

subroutine split_parameters(a, a_kep, a_dyn, dynamical_planets, k, d)
    implicit none
    integer :: MMAX, NPLMAX, NDSMAX
    parameter (MMAX = 200, NPLMAX = 10, NDSMAX = 20)
    real(8) :: a(MMAX), a_kep(MMAX), a_dyn(MMAX), PI
    parameter (PI = 3.14159265358979d0)
    integer :: dynamical_planets(NPLMAX), i, j, k, d
    integer :: npl, ndset, idset, gr_flag
    integer :: idsmax(NDSMAX)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    k = 0
    d = 0

    do i = 1, npl
        if (dynamical_planets(i).eq.1.d0) then
            do j = 1, 7
                a_dyn(7 * d + j) = a(7 * (i - 1) + j)
            enddo
            d = d + 1
        else
            do j = 1, 7
                a_kep(7 * k + j) = a(7 * (i - 1) + j)
            enddo
            k = k + 1
        endif
    enddo

    do i = 1, k
        j = 7 * (i - 1)
        a_kep(j + 2) = a_kep(j + 2) / 8.64d4
    enddo
end

subroutine RVKEP_dynamo (t, a, ymod, ymod_pl, ma, ndata, epsil, dt, hkl, &
        dynamical_planets, ts, coplar_inc)
    implicit none
    real(8) :: PI, TWOPI, eps, dt
    parameter (PI = 3.14159265358979d0, eps = 1.d-6)
    parameter (TWOPI = 2.0d0 * PI)
    integer :: npl, ma, i, j, NPLMAX, na, ndset, NDSMAX, ndata
    parameter (NPLMAX = 10, NDSMAX = 20)
    real(8) :: t(ndata), ymod(ndata), ymod_pl(npl,ndata), a(ma)
    real(8) :: mstar, ap(NPLMAX), mass(NPLMAX), epsil, a2(ma)
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX) 
    real(8) :: vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)
    real(8) :: xj(NPLMAX), yj(NPLMAX), zj(NPLMAX) 
    real(8) :: vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
    real(8) :: rpl(NPLMAX), rhill(NPLMAX)
    real(8) :: sini(NPLMAX)
    integer :: hkl, idset, nbod, gr_flag, coplar_inc
    integer :: dynamical_planets(npl), k, d, ts(ndata)
    real(8) :: a_kep(ma), a_dyn(ma), ymod_kep(ndata), dyda_kep(ma)

    integer :: idsmax(NDSMAX)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    common mstar, sini

    do i = 1, ma
        a2(i) = a(i)
    enddo

    if (hkl.eq.0) then

        do i = 1, npl
            j = 7 * (i - 1)
            !a2(j+2) = 2.d0*PI/(a2(j+2)*8.64d4)
            if (a2(j + 2)<0.d0) then !if P<0, set P>0
                a2(j + 2) = abs(a2(j + 2))
            endif

            if (a2(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a2(j + 4) = a2(j + 4) + PI
                a2(j + 1) = abs(a2(j + 1))
                if (a2(j + 4)>2.d0 * PI) a2(j + 4) = a2(j + 4) - 2.d0 * PI
            endif
            if (a2(j + 3)<0.d0) then !if e<0, set e>0 and w=w+PI, M0=M0-PI
                a2(j + 3) = abs(a2(j + 3))
                a2(j + 4) = a2(j + 4) + PI
                if (a2(j + 4)>2.d0 * PI) a2(j + 4) = a2(j + 4) - 2.d0 * PI
                a2(j + 5) = a2(j + 5) - PI
                if (a2(j + 5)<0.d0) a2(j + 5) = a2(j + 5) + 2.d0 * PI
            endif

            if(a2(j + 3)>=1.d0) then !if e>=1 set it to 0.99 to prevent errors
                a2(j + 3) = 0.99d0
            endif
            if (a2(j + 4)<0.d0) a2(j + 4) = dmod(a2(j + 4) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 5)<0.d0) a2(j + 5) = dmod(a2(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 4)>2.d0 * PI) a2(j + 4) = dmod(a2(j + 4), 2.d0 * PI)
            if (a2(j + 5)>2.d0 * PI) a2(j + 5) = dmod(a2(j + 5), 2.d0 * PI)
            if (a2(j + 6)<0.d0) a2(j + 6) = dmod(a2(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 7)<0.d0) a2(j + 7) = dmod(a2(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 6)>2.d0 * PI) a2(j + 6) = dmod(a2(j + 6), 2.d0 * PI)
            if (a2(j + 7)>2.d0 * PI) a2(j + 7) = dmod(a2(j + 7), 2.d0 * PI)
        enddo
    else
        do i = 1, npl
            j = 7 * (i - 1)
            !a2(j+2) = 2.d0*PI/(a2(j+2)*8.64d4)
            if (a2(j + 2)<0.d0) then  ! if P<0, set P>0
                a2(j + 2) = abs(a2(j + 2))
            endif

            if (a2(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a2(j + 4) = -1.d0 * a2(j + 4) !which is h = -h, k = -k
                a2(j + 3) = -1.d0 * a2(j + 3)
                a2(j + 1) = abs(a2(j + 1))
            endif

            if (a2(j + 5)<0.d0) a2(j + 5) = dmod(a2(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 6)<0.d0) a2(j + 6) = dmod(a2(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 7)<0.d0) a2(j + 7) = dmod(a2(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a2(j + 5)>2.d0 * PI) a2(j + 5) = dmod(a2(j + 5), 2.d0 * PI)
            if (a2(j + 6)>2.d0 * PI) a2(j + 6) = dmod(a2(j + 6), 2.d0 * PI)
            if (a2(j + 7)>2.d0 * PI) a2(j + 7) = dmod(a2(j + 7), 2.d0 * PI)
        enddo
    endif

    if (coplar_inc.ne.0) then
        do i = 2, npl
            j = 7 * (i - 1)
            a2(j + 6) = a2(6)
        enddo
    endif


    ymod(:) = 0.d0
    ymod_kep(:) = 0.d0

    if (npl.ne.0) then
        call split_parameters(a2, a_kep, a_dyn, dynamical_planets, k, d)

        nbod = d + 1
        na = 7 * d

        do i = 1, ndata !initialize ymod
            call RVKEP_dynamoplus(t(i), a_kep, ymod_kep(i), ymod_pl(:,i), &
                dyda_kep, ma - na, ts(i), k)
        enddo


        !get ymod first
        if(d.ne.0) then
            call MA_J_cop_fin (a_dyn, na, d, mstar, mass, ap, hkl)
            call GENINIT_J3_ewcop (nbod, ap, a_dyn, &
                    mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)
            call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
                    xh, yh, zh, vxh, vyh, vzh)
            call integrate_cop_fin(ymod, ymod_pl, t, nbod, ndata, mass, &
                    xh, yh, zh, vxh, vyh, vzh, epsil, dt)
        endif
    endif

    do i = 1, ndata
        ymod(i) = ymod(i) + ymod_kep(i)
    enddo
    return
end

subroutine RVKEP_ewcop_fin (t, a, ymod, ymod_pl, dyda, ma, ndata, epsil, &
        deltat, hkl, coplar_inc)
    implicit none
    real(8) :: PI, TWOPI, eps, epsil, deltat
    parameter (PI = 3.14159265358979d0, eps = 1.d-6)
    parameter (TWOPI = 2.0d0 * PI)
    integer :: npl, ma, i, j, NPLMAX, MMAX, na, ndset, NDSMAX, idset, ndata, nbod
    parameter (NPLMAX = 10, NDSMAX = 20, MMAX = 200)
    real(8) :: t(ndata), ymod(ndata), ymod_pl(npl,ndata)
    real(8) :: a(ma), dyda(20000, MMAX)
    real(8) :: mstar, ap(NPLMAX), mass(NPLMAX)
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) :: vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)
    real(8) :: xj(NPLMAX), yj(NPLMAX), zj(NPLMAX) 
    real(8) :: vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
    real(8) :: rpl(NPLMAX), rhill(NPLMAX)
    real(8) :: ah(ma), ahh(ma), ymodhb(ndata), ymodha(ndata)
    real(8) :: ymodha_pl(npl,ndata), ymodhb_pl(npl,ndata)
    real(8) :: sini, factor
    integer :: idsmax(NDSMAX), hkl, gr_flag, coplar_inc

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    common mstar, sini

    nbod = npl + 1
    na = 7 * npl

    if (hkl.eq.0) then

        do i = 1, npl
            j = 7 * (i - 1)
            !a2(j+2) = 2.d0*PI/(a2(j+2)*8.64d4)
            if (a(j + 2)<0.d0) then  ! if P<0, set P>0
                a(j + 2) = abs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a(j + 4) = a(j + 4) + PI
                a(j + 1) = abs(a(j + 1))
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
            endif
            if (a(j + 3)<0.d0) then !if e<0, set e>0 and w=w+PI, M0=M0-PI
                a(j + 3) = abs(a(j + 3))
                a(j + 4) = a(j + 4) + PI
                if (a(j + 4)>2.d0 * PI) a(j + 4) = a(j + 4) - 2.d0 * PI
                a(j + 5) = a(j + 5) - PI
                if (a(j + 5)<0.d0) a(j + 5) = a(j + 5) + 2.d0 * PI
            endif

            if(a(j + 3)>=1.d0) then !if e>=1 set it to 0.99 to prevent errors
                a(j + 3) = 0.99d0
            endif
            if (a(j + 4)<0.d0) a(j + 4) = dmod(a(j + 4) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 4)>2.d0 * PI) a(j + 4) = dmod(a(j + 4), 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            if (a(j + 6)<0.d0) a(j + 6) = dmod(a(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 7)<0.d0) a(j + 7) = dmod(a(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 6)>2.d0 * PI) a(j + 6) = dmod(a(j + 6), 2.d0 * PI)
            if (a(j + 7)>2.d0 * PI) a(j + 7) = dmod(a(j + 7), 2.d0 * PI)
        enddo
    else
        do i = 1, npl
            j = 7 * (i - 1)
            !a(j+2) = 2.d0*PI/(a(j+2)*8.64d4)
            if (a(j + 2)<0.d0) then  ! if P<0, set P>0
                a(j + 2) = abs(a(j + 2))
            endif

            if (a(j + 1)<0.d0) then !if K<0, set K>0 and w = w+PI
                a(j + 4) = -1.d0 * a(j + 4) !which is h = -h, k = -k
                a(j + 3) = -1.d0 * a(j + 3)
                a(j + 1) = abs(a(j + 1))
            endif

            if (a(j + 5)<0.d0) a(j + 5) = dmod(a(j + 5) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 6)<0.d0) a(j + 6) = dmod(a(j + 6) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 7)<0.d0) a(j + 7) = dmod(a(j + 7) + 2.d0 * PI, 2.d0 * PI)
            if (a(j + 5)>2.d0 * PI) a(j + 5) = dmod(a(j + 5), 2.d0 * PI)
            if (a(j + 6)>2.d0 * PI) a(j + 6) = dmod(a(j + 6), 2.d0 * PI)
            if (a(j + 7)>2.d0 * PI) a(j + 7) = dmod(a(j + 7), 2.d0 * PI)
        enddo
    endif

    if (coplar_inc.ne.0) then
        do i = 2, npl
             j = 7 * (i - 1)
             a(j + 6) = a(6)
        enddo
    endif
     

    if (npl.ne.0) then
        !get ymod first
        call MA_J_cop_fin (a, ma, npl, mstar, mass, ap, hkl)

        call GENINIT_J3_ewcop (nbod, ap, a, &
                mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)

        call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
                xh, yh, zh, vxh, vyh, vzh)

        do i = 1, ndata !initialize ymod
            ymod(i) = 0.0
        enddo

        call integrate_cop_fin(ymod, ymod_pl, t, nbod, ndata, mass, &
                xh, yh, zh, vxh, vyh, vzh, epsil, deltat)

        do i = 1, ma
            do j = 1, ndata
                dyda(j, i) = 0.0
            enddo
        enddo

        factor = 1.0d0

        do i = 0, npl - 1
            ahh(7 * i + 1) = 0.05d0 * factor
            ahh(7 * i + 2) = 8640d0
            ahh(7 * i + 3) = 0.0001d0 * factor
            ahh(7 * i + 4) = 0.01d0 * factor
            ahh(7 * i + 5) = 0.01d0 * factor
            ahh(7 * i + 6) = 0.01d0 * factor
            ahh(7 * i + 7) = 0.01d0 * factor
            !ahh(7*npl + ndset + 1) = 0.02d0*factor
        enddo
        !calculate dyda
        do i = 1, na

            do j = 1, na
                ah(j) = a(j)
            enddo

            ! get ah of back
            if (mod(i, 7).ne.44) then
                ah(i) = a(i) - ahh(i)

                call MA_J_cop_fin (ah, ma, npl, mstar, mass, ap, hkl)
                call GENINIT_J3_ewcop (nbod, ap, ah, &
                        mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)
                call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
                        xh, yh, zh, vxh, vyh, vzh)
                call integrate_cop_fin(ymodhb, ymodhb_pl,t, nbod, ndata, mass, &
                        xh, yh, zh, vxh, vyh, vzh, epsil, deltat)
                !get ah of ahead
                ah(i) = a(i) + ahh(i)

                call MA_J_cop_fin (ah, ma, npl, mstar, mass, ap, hkl)
                call GENINIT_J3_ewcop (nbod, ap, ah, &
                        mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)
                call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
                        xh, yh, zh, vxh, vyh, vzh)
                call integrate_cop_fin(ymodha, ymodha_pl,t, nbod, ndata, mass, &
                        xh, yh, zh, vxh, vyh, vzh, epsil, deltat)
                !calculate the ith dyda
                do j = 1, ndata
                    dyda(j, i) = (ymodha(j) - ymodhb(j)) / (2.d0 * ahh(i))
                enddo
            else
                ah(i) = a(i) - ahh(i)
                !ah(i) = dmod(a(i) - ahh(i),  0.0d0 )
                !if (ah(j+6)<=0.d0) ah(j+6) = dmod(ah(j+6) + PI,  PI )
                !if (ah(j+6)>=PI)   ah(j+6) = dmod(ah(j+6) +0.00001,  PI)

                call MA_J_cop_fin (ah, ma, npl, mstar, mass, ap, hkl)
                call GENINIT_J3_ewcop (nbod, ap, ah, &
                        mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)
                call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
                        xh, yh, zh, vxh, vyh, vzh)
                call integrate_cop_fin(ymodhb, ymodhb_pl,t, nbod, ndata, mass, &
                        xh, yh, zh, vxh, vyh, vzh, epsil, deltat)
                !get ah of ahead
                ah(i) = a(i) + ahh(i)
                !ah(i) = dmod(a(i) + ahh(i),  PI )
                call MA_J_cop_fin (ah, ma, npl, mstar, mass, ap, hkl)
                call GENINIT_J3_ewcop (nbod, ap, ah, &
                        mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)
                call coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
                        xh, yh, zh, vxh, vyh, vzh)
                call integrate_cop_fin(ymodha,ymodha_pl, t, nbod, ndata, mass, &
                        xh, yh, zh, vxh, vyh, vzh, epsil, deltat)
                !calculate the ith dyda
                do j = 1, ndata
                    dyda(j, i) = (ymodha(j) - ymodhb(j)) / (2.d0 * ahh(i))
                enddo
            endif
        enddo
    endif
    return
end

subroutine MA_J_keplm (a, ma, npl, m0, mass, ap, hkl, gr_flag)
    implicit none
    real(8) :: m0, PI, TWOPI, THIRD, GMSUN, dm, MSUN
    integer :: npl, ma, i, j, NPLMAX, hkl, gr_flag
    parameter (NPLMAX = 10)
    real(8) :: mm(NPLMAX), ecc, corr
    real(8) :: a(ma), mass(NPLMAX), ap(NPLMAX), mpold(NPLMAX), mtotal
    parameter (THIRD = 1.d0 / 3.d0)
    parameter (PI = 3.14159265358979d0, TWOPI = 2.d0 * PI)
    parameter (GMSUN = 1.32712440018d20, MSUN = 1.32712440018d20)

    !*******G is set to be unit, and s, m, kg as unit of time, length and mass expectively.

    do j = 1, npl
        i = 6 * (j - 1)
        mm(j) = 2.d0 * PI / (a(i + 2) * 8.64d4)
    enddo

    do i = 0, npl - 1
        if (hkl.eq.0) then
            ecc = a(6 * i + 3)
        else
            ecc = dsqrt(a(6 * i + 3)**2 + a(6 * i + 4)**2) !only for h, k
        endif

        mass(1) = m0
        mpold(i + 1) = 0.d0
        1010      continue
        if (i.eq.0) then
            mtotal = m0
            mass(i + 2) = a(6 * i + 1) * (TWOPI / mm(i + 1) * (m0 + mpold(i + 1))**2 / &
                    (TWOPI * GMSUN))**THIRD * &
                    dsqrt(1.d0 - ecc**2)
        else
            mtotal = m0
            do j = 0, i - 1
                mtotal = mtotal + mass(j + 2)
            enddo
            mass(i + 2) = a(6 * i + 1) * (TWOPI / mm(i + 1) * (mtotal&
                    + mpold(i + 1))**2 / (TWOPI * GMSUN))**THIRD * &
                    dsqrt(1.d0 - ecc**2)
        endif

        dm = dabs(mass(i + 2) - mpold(i + 1)) / mass(i + 2)
        mpold(i + 1) = mass(i + 2)
        if (dm>0) goto 1010

        ap(i + 1) = (GMSUN * (mtotal + mass(i + 2)) * (1.d0 / mm(i + 1))&
                **2)**THIRD
    enddo

    do i = 1, npl + 1
        mass(i) = mass(i) * MSUN
    enddo

    if(gr_flag.ne.0) then
        do i = 1, npl
            j = 6 * (i - 1)
            call gr_corr(ap(i), a(j + 3), mass(1), corr)
            a(j + 6) = corr * 365.25
        enddo
    endif
    return
end

! MA_J calculates the actual masses and Jacobi semimajor axes of a two-planet
! system for assumed sin(i) using the parameters P, K and e from a
! two-Kepler fit.
subroutine MA_J_cop_fin (a, ma, npl, m0, mass, ap, hkl)
    implicit none
    real(8) :: m0, PI, TWOPI, THIRD, GMSUN, dm, MSUN, ecc
    integer :: npl, ma, i, j, NPLMAX, hkl
    parameter (NPLMAX = 10)
    real(8) :: a(ma), mass(NPLMAX), ap(NPLMAX), mpold(NPLMAX), mtotal
    parameter (THIRD = 1.d0 / 3.d0)
    parameter (PI = 3.14159265358979d0, TWOPI = 2.d0 * PI)
    parameter (GMSUN = 1.32712440018d20, MSUN = 1.32712440018d20)

    !*******G is set to be unit, and s, m, kg as unit of time, length and mass expectively.

    do i = 0, npl - 1
        if (hkl.eq.0) then
            ecc = a(7 * i + 3)
        else
            ecc = dsqrt(a(7 * i + 3)**2 + a(7 * i + 4)**2) !only for h, k
        endif

        mass(1) = m0
        mpold(i + 1) = 0.d0
        101        continue
        if (i.eq.0) then
            mtotal = m0
            mass(i + 2) = dabs(a(7 * i + 1)) * (a(7 * i + 2) * (m0 + mpold(i + 1))**2 / &
                    (TWOPI * GMSUN))**THIRD * dsqrt(1.d0 - ecc**2)&
                    / dabs(dsin(a(7 * i + 6)))
        else
            mtotal = m0
            do j = 0, i - 1
                mtotal = mtotal + mass(j + 2)
            enddo
            mass(i + 2) = dabs(a(7 * i + 1)) * (a(7 * i + 2) * (mtotal&
                    + mpold(i + 1))**2 / (TWOPI * GMSUN))**THIRD&
                    * dsqrt(1.d0 - ecc**2) / dabs(dsin(a(7 * i + 6)))
        endif

        dm = dabs(mass(i + 2) - mpold(i + 1)) / mass(i + 2)
        mpold(i + 1) = mass(i + 2)
        if (dm>0) goto 101

        ap(i + 1) = (GMSUN * (mtotal + mass(i + 2)) * &
                (a(7 * i + 2) / TWOPI)**2)**THIRD
    enddo

    do i = 1, npl + 1
        mass(i) = mass(i) * MSUN
    enddo
    return
end

! GENINIT_J3 reads Jacobi orbital elements of nbod planets and generates
! initial position and velocity in Jacobi coords.
! This version outputs rpl and rhill.
! Last modified by Man Hoi Lee, Aug 16, 2003.
subroutine GENINIT_J3_ewcop (nbod, ap, a, &
        mass, xj, yj, zj, vxj, vyj, vzj, rpl, rhill, hkl)
    implicit none                
    real(8) :: SMASSYR, MSUN, PI, eps, THIRD
    parameter (PI = 3.14159265358979d0, eps = 1.d-7)
    parameter (SMASSYR = 4.d0 * PI * PI)
    parameter (MSUN = 1.32712440018d20)
    parameter (THIRD = 1.d0 / 3.d0)

    integer ::  nbod, NPLMAX, i, j, hkl
    parameter (NPLMAX = 10)
    real(8) :: mass(NPLMAX), ecc
    real(8) :: xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
    real(8) :: vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
    real(8) :: rpl(NPLMAX), rhill(NPLMAX)
    real(8) :: frho3, ap(NPLMAX), a(NPLMAX)
    real(8) :: gm, inc, capom, omega, capm
    integer ::  ialpha

    ! SET F/RHO^(1/3) FOR RADIUS (RHO IN G/CM^3) TO 1.D0 FOR NOW.
    frho3 = 1.d0

    xj(1) = 0.d0
    yj(1) = 0.d0
    zj(1) = 0.d0
    vxj(1) = 0.d0
    vyj(1) = 0.d0
    vzj(1) = 0.d0

    ialpha = -1
    inc = 0.d0
    capom = 0.d0
    gm = mass(1)

    do i = 2, nbod
        j = 7 * (i - 2)
        if (hkl.eq.0) then
            ecc = a(j + 3)
            omega = a(j + 4)
            capm = a(j + 5)
        else
            ecc = dsqrt(a(j + 3)**2 + a(j + 4)**2)
            omega = datan2(a(j + 3), a(j + 4))
            if (omega<0.d0) omega = dmod(omega + 2.d0 * PI, 2.d0 * PI)
            if (omega>2.d0 * PI) omega = dmod(omega, 2.d0 * PI)

            capm = a(j + 5) - omega

            if (capm<0.d0) capm = dmod(capm + 2.d0 * PI, 2.d0 * PI)
            if (capm>2.d0 * PI) capm = dmod(capm, 2.d0 * PI)
        endif

        gm = gm + mass(i)
        rpl(i) = frho3 * (1.5d0 * mass(i) / 2.d0 * PI)**THIRD
        rhill(i) = ap(i - 1) * (mass(i) / (3.d0 * mass(1)))**THIRD

        call ORBEL_EL2XV (gm, ialpha, ap(i - 1), ecc, a(j + 6), a(j + 7), &
                omega, capm, xj(i), yj(i), zj(i), vxj(i), vyj(i), vzj(i))
    enddo
    return
end

subroutine integrate_cop_fin(ymod, ymod_pl, t, nbod, ndata, mass, &
        xh, yh, zh, vxh, vyh, vzh, eps, dt)
 
    include 'rvmod.inc'

    integer ::  IO_NBITS
    parameter(IO_NBITS = 6)
    real(8) :: mass(NPLMAX), j2rp2, j4rp4
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) :: vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)

    integer ::  i1st
    integer ::  nbod, ntp, nleft
    integer ::  iub, iuj, iud, iue

    real(8) :: tstop, dt
    real(8) :: time

    integer ::  ndata, i, j, nd, flag
    real(8) :: ymod(ndata), ymod_pl(nbod-1,ndata), t(ndata)
    real(8) :: h, eps, mtotal

    ntp = 0
    ! Prompt and read name of planet data file
    j2rp2 = 0.0d0
    j4rp4 = 0.0d0

    ! Initialize initial time and times for first output and first dump
    time = 0.d0
    tstop = t(ndata)

    iub = 20
    iuj = 30
    iud = 40
    iue = 60

    nleft = ntp
    i1st = 0

    do i = 1, ndata !initialize ymod
        ymod(i) = 0.0
!        write(*,*) nbod
!        do j = 1, nbod-1
!           ymod_pl(j,i) = 0.0
!        enddo
    enddo
    !---------here is the big loop--------------------------
    nd = 1
    !------output the first ymod if time of dataset begins from 0
    if (abs(t(1))<1.d-10) then
        mtotal = 0.d0
        do i = 1, nbod
            mtotal = mass(i) + mtotal
        enddo
        do i = 2, nbod
            ymod_pl(i-1,nd) = mass(i) / mtotal * vzh(i)            
            ymod(nd) = ymod(nd) + mass(i) / mtotal * vzh(i)
        enddo
        nd = nd + 1
    endif

    !-------loop-----
    do while(time<=tstop)
        h = dt
        flag = 0 !flag for controling output
        do i = 1, ndata
            if ((time<t(i)).and.((t(i) - time)<=dt)) then
                flag = 1
                h = t(i) - time
                exit
            endif
        enddo

        call bs_step2(nbod, mass, j2rp2, j4rp4, &
                xh, yh, zh, vxh, vyh, vzh, h, eps)

        if (flag.eq.1) then
            mtotal = 0.d0
            do i = 1, nbod
                mtotal = mtotal + mass(i)
            enddo
            do i = 2, nbod
                j = 7 * (i - 2)
                ymod_pl(i-1,nd) = mass(i) / mtotal * vzh(i)                           
                ymod(nd) = ymod(nd) + mass(i) / mtotal * vzh(i)
            enddo
            nd = nd + 1
        endif

        time = time + h
    enddo
    return
end

subroutine bs_step2(nbod, mass, j2rp2, j4rp4, &
        xh, yh, zh, vxh, vyh, vzh, dt, eps)
    include 'rvmod.inc'
    include 'bs.inc'

    !...  Inputs Only:
    integer ::  nbod
    real(8) :: mass(nbod), dt, j2rp2, j4rp4

    !...  Inputs and Outputs:
    real(8) :: xh(nbod), yh(nbod), zh(nbod)
    real(8) :: vxh(nbod), vyh(nbod), vzh(nbod)

    !...  Internals
    integer ::  i, ntpi
    real(8) :: xb(NPLMAX), yb(NPLMAX), zb(NPLMAX), eps
    real(8) :: vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)
    real(8) :: ybs(6, (NTPMAX + NPLMAX)), tfake, dttmp, msys

    !...  Executable code

    !...  set things up if this is the initial call
    !...  Convert to barycentric coords

    call coord_h2b(nbod, mass, xh, yh, zh, vxh, vyh, vzh, &
            xb, yb, zb, vxb, vyb, vzb, msys)

    !call coord_h2b_tp(ntp,xht,yht,zht,vxht,vyht,vzht, &
    !     xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1), &
    !     xbt,ybt,zbt,vxbt,vybt,vzbt)

    !...  copy to the big array
    do i = 1, nbod
        ybs(1, i) = xb(i)
        ybs(2, i) = yb(i)
        ybs(3, i) = zb(i)
        ybs(4, i) = vxb(i)
        ybs(5, i) = vyb(i)
        ybs(6, i) = vzb(i)
    enddo

    ntpi = 0

    tfake = 0.0d0
    dttmp = dt

    do while((abs(tfake - dt) / dt) > 1.0e-7) !just to be real safe
        call bs_int_pl(nbod, ntpi, mass, j2rp2, j4rp4, &
                tfake, dttmp, ybs, eps)
        dttmp = dt - tfake
    enddo

    !...  put things back
    do i = 1, nbod
        xb(i) = ybs(1, i)
        yb(i) = ybs(2, i)
        zb(i) = ybs(3, i)
        vxb(i) = ybs(4, i)
        vyb(i) = ybs(5, i)
        vzb(i) = ybs(6, i)
    enddo

    !...  Convert back to helio. coords at the end of the step
    call coord_b2h(nbod, xb, yb, zb, vxb, vyb, vzb, &
            xh, yh, zh, vxh, vyh, vzh)
    return
end

! MRQMIN attempts to reduce the chi^2 of a fit by the Levenberg-Marquardt
! method. It uses COVSRT, GAUSSJ, and MRQCOF.
! From Numerical Recipes.
subroutine MRQMIN_dynamo (x, y, sig, ndata, a, ia, ma, ts, covar, alpha, &
        nca, chisq, funcs, alamda, loglik, jitt, hkl)
    implicit none
    integer ::  ma, nca, ndata, ia(ma), MMAX, NDSMAX 
    integer ::  ts(ndata)    
    real(8) :: alamda, chisq, a(ma), alpha(nca, nca), covar(nca, nca), &
            sig(ndata), x(ndata), y(ndata), loglik
    external funcs
    parameter (MMAX = 200, NDSMAX = 20)
    integer ::  j, k, l, mfit, hkl
    real(8) :: ochisq, atry(MMAX), beta(MMAX), da(MMAX), jitt(NDSMAX)
    save ochisq, atry, beta, da, mfit

    ! Initialization.
    if (alamda<0.d0) then
        mfit = 0
        do j = 1, ma
            if (ia(j).ne.0) mfit = mfit + 1
        enddo
        alamda = 0.001d0

        call MRQCOF_dynamo (x, y, sig, ndata, a, ia, ma, ts, &
                alpha, beta, nca, chisq, funcs, loglik, jitt, hkl)

        ochisq = chisq
        do j = 1, ma
            atry(j) = a(j)
        enddo
    endif

    ! Alter linearized fitting matrix by augmenting diagonal elements.
    do j = 1, mfit
        do k = 1, mfit
            covar(j, k) = alpha(j, k)
        enddo
        covar(j, j) = alpha(j, j) * (1.d0 + alamda)
        da(j) = beta(j)
    enddo

    ! Matrix solution.
    call GAUSSJ_dynamo (covar, mfit, nca, da, 1, 1)

    ! Evaluate covariance matrix once converged.
    if (alamda.eq.0.d0) then
        call COVSRT (covar, nca, ma, ia, mfit)
        call COVSRT (alpha, nca, ma, ia, mfit)
        return
    endif

    j = 0
    do l = 1, ma
        if (ia(l).ne.0) then
            j = j + 1
            atry(l) = a(l) + da(j)
        endif
    enddo

    call MRQCOF_dynamo (x, y, sig, ndata, atry, ia, ma, &
            ts, covar, da, nca, chisq, funcs, loglik, jitt, hkl)

    if (chisq<ochisq) then
        !Accept new solution.
        alamda = 0.1d0 * alamda
        ochisq = chisq
        do j = 1, mfit
            do k = 1, mfit
                alpha(j, k) = covar(j, k)
            enddo
            beta(j) = da(j)
        enddo
        do l = 1, ma
            a(l) = atry(l)
        enddo
    else
        !Increase alamda and return.
        alamda = 10.d0 * alamda
        chisq = ochisq
    endif
    return
end

! MRQMIN attempts to reduce the chi^2 of a fit by the Levenberg-Marquardt
! method. It uses COVSRT, GAUSSJ, and MRQCOF.
! From Numerical Recipes.
subroutine MRQMIN_dynlm (x, ts, y, sig, ndata, a, ia, ma, covar, alpha, &
        nca, chisq, funcs, alamda, loglik, jitt, epsil, deltat, hkl, coplar_inc)
    implicit none
    integer ::  ma, nca, ndata, ia(ma), MMAX, NDSMAX, ts(ndata)
    parameter (MMAX = 200, NDSMAX = 20)
    integer ::  npl, ndset, idset, idsmax(NDSMAX), hkl, gr_flag
    real(8) :: alamda, chisq, a(ma), alpha(nca, nca), covar(nca, nca), &
            sig(ndata), x(ndata), y(ndata), loglik, jitt(NDSMAX)
    external funcs

    integer ::  j, k, l, mfit, coplar_inc
    real(8) :: ochisq, atry(MMAX), beta(MMAX), da(MMAX), epsil, deltat
    save ochisq, atry, beta, da, mfit

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    ! Initialization.
    if (alamda<0.d0) then
        mfit = 0
        do j = 1, ma
            if (ia(j).ne.0) mfit = mfit + 1
        enddo
        alamda = 0.001d0

        call MRQCOF_dynlm (x, ts, y, sig, ndata, a, ia, ma, alpha, beta, &
                nca, chisq, funcs, loglik, jitt, epsil, deltat, hkl, coplar_inc)

        ochisq = chisq

        do j = 1, ma
            atry(j) = a(j)

        enddo
    endif
    ! Alter linearized fitting matrix by augmenting diagonal elements.
    do j = 1, mfit
        do k = 1, mfit
            covar(j, k) = alpha(j, k)
        enddo
        covar(j, j) = alpha(j, j) * (1.d0 + alamda)
        da(j) = beta(j)
    enddo

    ! Matrix solution.
    call GAUSSJ_dynlm (covar, mfit, nca, da, 1, 1)

    ! Evaluate covariance matrix once converged.
    if (alamda.eq.0.d0) then
        call COVSRT (covar, nca, ma, ia, mfit)
        call COVSRT (alpha, nca, ma, ia, mfit)
        return
    endif

    j = 0
    do l = 1, ma
        if (ia(l).ne.0) then
            j = j + 1
            atry(l) = a(l) + da(j)

        endif
    enddo

    call MRQCOF_dynlm (x, ts, y, sig, ndata, atry, ia, ma, covar, da, &
            nca, chisq, funcs, loglik, jitt, epsil, deltat, hkl, coplar_inc)

    if (chisq<ochisq) then
        !Accept new solution.
        alamda = 0.1d0 * alamda
        ochisq = chisq
        do j = 1, mfit
            do k = 1, mfit
                alpha(j, k) = covar(j, k)
            enddo
            beta(j) = da(j)
        enddo
        do l = 1, ma
            a(l) = atry(l)
        enddo
    else
        !Increase alamda and return.
        alamda = 10.d0 * alamda
        chisq = ochisq
    endif
    return
end

! MRQCOF evaluates the linearized fitting matrix alpha, the vector beta,
! and chi^2.
! From Numerical Recipes.
subroutine MRQCOF_dynamo (x, y, sig, ndata, a, ia, ma, ts, alpha, &
        beta, nalp, chisq, funcs, loglik, jitt, hkl)
    implicit none
    integer ::  npl, ndset, idset, ma, nalp, ndata, ia(ma), NDSMAX, MMAX
    parameter (NDSMAX = 20, MMAX = 200)
    integer ::  idsmax(NDSMAX), ts(ndata), hkl, gr_flag
    real(8) :: chisq, a(ma), alpha(nalp, nalp), beta(ma), sig(ndata), &
            x(ndata), y(ndata), loglik, TWOPI, jitt(NDSMAX)
    parameter (TWOPI = 2.d0 * 3.14159265358979d0)
    external funcs
    integer ::  mfit, i, j, k, l, m
    real(8) :: dy, sig2i, wt, ymod,ymod_pl(npl), dyda(MMAX)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag

    mfit = 0
    do j = 1, ma
        if (ia(j).ne.0) mfit = mfit + 1
    enddo

    ! Initialize (symmetric) alpha and beta.
    do j = 1, mfit
        do k = 1, j
            alpha(j, k) = 0.d0
        enddo
        beta(j) = 0.d0
    enddo
    chisq = 0.d0
    loglik = 0.d0

    ! Loop over all data.
    idset = 1
    do i = 1, ndata
        idset = ts(i)
        call FUNCS (x(i), a, ymod, ymod_pl, dyda, ma, idset, hkl)
        sig2i = 1.d0 / (sig(i)**2 + jitt(idset)**2)
        dy = y(i) - ymod
        j = 0

        do l = 1, ma
            if (ia(l).ne.0) then
                j = j + 1
                wt = dyda(l) * sig2i
                k = 0
                do m = 1, l
                    if (ia(m).ne.0) then
                        k = k + 1
                        alpha(j, k) = alpha(j, k) + wt * dyda(m)
                    endif
                enddo
                beta(j) = beta(j) + dy * wt
            endif
        enddo

        !Compute chi^2 and loglik.
        chisq = chisq + dy * dy * sig2i
        loglik = loglik - 0.5 * dy * dy * sig2i - &
                0.5 * dlog(TWOPI * (sig(i)**2 + &
                        jitt(idset)**2))
    enddo
    ! Fill in the symmetric side.
    do j = 2, mfit
        do k = 1, j - 1
            alpha(k, j) = alpha(j, k)
        enddo
    enddo
    return
end

! MRQCOF evaluates the linearized fitting matrix alpha, the vector beta,
! and chi^2.
! From Numerical Recipes.
subroutine MRQCOF_dynlm (x, ts, y, sig, ndata, a, ia, ma, alpha, beta, nalp, &
        chisq, funcs, loglik, jitt, epsil, deltat, hkl, coplar_inc)
    implicit none
    integer ::  npl, ndset, idset, ma, nalp, ndata, ia(ma), NDSMAX, MMAX
    parameter (NDSMAX = 20, MMAX = 200)
    integer ::  idsmax(NDSMAX), ts(ndata), hkl
    real(8) :: chisq, a(ma), alpha(nalp, nalp), beta(ma), sig(ndata), &
            x(ndata), y(ndata), loglik, TWOPI, epsil, deltat
    parameter (TWOPI = 2.d0 * 3.14159265358979d0)
    external funcs
    integer ::  mfit, i, j, k, l, m, gr_flag, coplar_inc
    real(8) :: dy, sig2i, wt, ymod(ndata)
    real(8) :: ymod_pl(npl,ndata), dyda(20000, MMAX), jitt(NDSMAX)

    common /DSBLK/ npl, ndset, idsmax, idset, gr_flag
    save dyda

    mfit = 0
    do j = 1, ma
        if (ia(j).ne.0) mfit = mfit + 1
    enddo

    ! Initialize (symmetric) alpha and beta.
    do j = 1, mfit
        do k = 1, j
            alpha(j, k) = 0.d0
        enddo
        beta(j) = 0.d0
    enddo
    chisq = 0.d0
    loglik = 0.d0

    ! Loop over all data.
    call FUNCS (x, a, ymod, ymod_pl, dyda, ma, ndata, epsil, deltat, hkl, coplar_inc)

    do i = 1, ndata
        idset = ts(i)
        ymod(i) = ymod(i) + a(7 * npl + idset)
        dyda(i, 7 * npl + idset) = 1.d0

        !lin. trend:
        ymod(i) = ymod(i) + a(7 * npl + ndset + 1) * (x(i) / 86400.d0)&
                + a(7 * npl + ndset + 2) * (x(i) / 86400.d0)**2
        dyda(i, 7 * npl + ndset + 1) = (x(i) / 86400.d0)
        dyda(i, 7 * npl + ndset + 2) = (x(i) / 86400.d0)**2

        do j = idset + 1, ndset
            dyda(i, 7 * npl + j) = 0.d0
        enddo
        sig2i = 1.d0 / (sig(i)**2 + jitt(idset)**2)
        dy = y(i) - ymod(i)
        j = 0
        do l = 1, ma
            if (ia(l).ne.0) then
                j = j + 1
                wt = dyda(i, l) * sig2i
                k = 0
                do m = 1, l
                    if (ia(m).ne.0) then
                        k = k + 1
                        alpha(j, k) = alpha(j, k) + wt * dyda(i, m)
                    endif
                enddo
                beta(j) = beta(j) + dy * wt
            endif
        enddo

        chisq = chisq + dy * dy * sig2i
        loglik = loglik - 0.5 * dy * dy * sig2i - &
                0.5 * dlog(TWOPI * (sig(i)**2 + &
                        jitt(idset)**2))
    enddo

    j = 0
    ! Fill in the symmetric side.
    do j = 2, mfit
        do k = 1, j - 1
            alpha(k, j) = alpha(j, k)
        enddo
    enddo
    j = 0
    return
end

! GAUSSJ solves linear equation by Gauss-Jordan elimination.
! From Numerical Recipes.
subroutine GAUSSJ_dynamo (a, n, np, b, m, mp)
    implicit none
    integer ::  m, mp, n, np, NMAX
    real(8) :: a(np, np), b(np, mp)
    parameter (NMAX = 50)
    integer ::  i, icol, irow, j, k, l, ll, kkk, &
            indxc(NMAX), indxr(NMAX), ipiv(NMAX)
    real(8) :: big, dum, pivinv

    icol = 0
    irow = 0

    do j = 1, n
        ipiv(j) = 0
    enddo

    ! Main loop over the columns to be reduced.
    do i = 1, n
        !Search for a pivot element.
        big = 0.d0
        do j = 1, n
            if (ipiv(j).ne.1) then
                do k = 1, n
                    if (ipiv(k).eq.0) then
                        if (dabs(a(j, k))>=big) then
                            big = dabs(a(j, k))
                            irow = j
                            icol = k
                        endif
                    endif
                enddo
            endif
        enddo
        ipiv(icol) = ipiv(icol) + 1

        !Interchange rows to put the pivot element on the diagonal.
        if (irow.ne.icol) then
            do l = 1, n
                dum = a(irow, l)
                a(irow, l) = a(icol, l)
                a(icol, l) = dum
            enddo
        endif
        indxr(i) = irow
        indxc(i) = icol

        !Divide pivot row by the pivot element.
        pivinv = 1.d0 / a(icol, icol)
        a(icol, icol) = 1.d0
        do l = 1, n
            a(icol, l) = a(icol, l) * pivinv
        enddo
        do l = 1, m
            b(icol, l) = b(icol, l) * pivinv
        enddo

        !Reduce the rows, except for the pivot one.
        do ll = 1, n
            if (ll.ne.icol) then
                dum = a(ll, icol)
                a(ll, icol) = 0.d0
                do l = 1, n
                    a(ll, l) = a(ll, l) - a(icol, l) * dum
                enddo
                do l = 1, m
                    b(ll, l) = b(ll, l) - b(icol, l) * dum
                enddo
            endif
        enddo
    enddo
    ! Unscramble the solution in view of the column interchanges.
    do l = n, 1, -1
        if (indxr(l).ne.indxc(l)) then
            do k = 1, n
                dum = a(k, indxr(l))
                a(k, indxr(l)) = a(k, indxc(l))
                a(k, indxc(l)) = dum
            enddo
            do kkk = 1, m
                dum = b(indxr(l), kkk)
                b(indxr(l), kkk) = b(indxc(l), kkk)
                b(indxc(l), kkk) = dum
            enddo
        endif
    enddo
    return
end

! GAUSSJ solves linear equation by Gauss-Jordan elimination.
! From Numerical Recipes.
subroutine GAUSSJ_dynlm (a, n, np, b, m, mp)
    implicit none
    integer ::  m, mp, n, np, NMAX
    real(8) :: a(np, np), b(np, mp)
    parameter (NMAX = 51)
    integer ::  i, icol, irow, j, k, l, ll, indxc(NMAX), indxr(NMAX), ipiv(NMAX)
    real(8) :: big, dum, pivinv

    icol = 0
    irow = 0

    do j = 1, n
        ipiv(j) = 0
    enddo
    ! Main loop over the columns to be reduced.
    do i = 1, n
        !Search for a pivot element.
        big = 0.d0
        do j = 1, n
            if (ipiv(j).ne.1) then
                do k = 1, n
                    if (ipiv(k).eq.0) then
                        if (dabs(a(j, k))>=big) then
                            big = dabs(a(j, k))
                            irow = j
                            icol = k
                        endif
                    endif
                enddo
            endif
        enddo
        ipiv(icol) = ipiv(icol) + 1

        !Interchange rows to put the pivot element on the diagonal.
        if (irow.ne.icol) then
            do l = 1, n
                dum = a(irow, l)
                a(irow, l) = a(icol, l)
                a(icol, l) = dum
            enddo
            do l = 1, m
                dum = b(irow, l)
                b(irow, l) = b(icol, l)
                b(icol, l) = dum
            enddo
        endif

        indxr(i) = irow
        indxc(i) = icol

        if (a(icol, icol).ne.0.d0) then
            !Divide pivot row by the pivot element.
            pivinv = 1.d0 / a(icol, icol)

            a(icol, icol) = 1.d0
            do l = 1, n
                a(icol, l) = a(icol, l) * pivinv
            enddo
            do l = 1, m
                b(icol, l) = b(icol, l) * pivinv
            enddo
        endif

        !Reduce the rows, except for the pivot one.
        do ll = 1, n
            if (ll.ne.icol) then
                dum = a(ll, icol)
                a(ll, icol) = 0.d0
                do l = 1, n
                    a(ll, l) = a(ll, l) - a(icol, l) * dum
                enddo
                do l = 1, m
                    b(ll, l) = b(ll, l) - b(icol, l) * dum
                enddo
            endif
        enddo
    enddo
    ! Unscramble the solution in view of the column interchanges.
    do l = n, 1, -1
        if (indxr(l).ne.indxc(l)) then
            do k = 1, n
                dum = a(k, indxr(l))
                a(k, indxr(l)) = a(k, indxc(l))
                a(k, indxc(l)) = dum
            enddo
        endif
    enddo
    return
end

! COVSRT expands in storage the covariance matrix covar, so as to take
! into account parameters that are being held fixed. (For the latter,
! return zero covariances.)
! From Numerical Recipes.
subroutine COVSRT (covar, npc, ma, ia, mfit)
    implicit none
    integer ::  ma, mfit, npc, ia(ma)
    real(8) :: covar(npc, npc)
    integer ::  i, j, k
    real(8) :: swap

    do i = mfit + 1, ma
        do j = 1, i
            covar(i, j) = 0.d0
            covar(j, i) = 0.d0
        enddo
    enddo

    k = mfit
    do j = ma, 1, -1
        if (ia(j).ne.0) then
            do i = 1, ma
                swap = covar(i, k)
                covar(i, k) = covar(i, j)
                covar(i, j) = swap
            enddo
            do i = 1, ma
                swap = covar(k, i)
                covar(k, i) = covar(j, i)
                covar(j, i) = swap
            enddo
            k = k - 1
        endif
    enddo
    return
end

!*****************************************************************************
!*                          ORBEL_EL2XV.F
!*****************************************************************************
!     PURPOSE: To compute cartesian positions and velocities given
!               central mass, ialpha ( = +1 for hyp., 0 for para. and
!               -1 for ellipse), and orbital elements.
!       input:
!            gm       ==> G times central mass (real scalar)
!         ialpha   ==> conic section type ( see PURPOSE, integer scalar)
!         a        ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            inc      ==> inclination  (real scalar)
!            capom    ==> longitude of ascending node (real scalar)
!         omega    ==> argument of perihelion (real scalar)
!         capm     ==> mean anomoly(real scalar)
!       Output:
!            x,y,z    ==>  position of object (real scalars)
!            vx,vy,vz ==>  velocity of object (real scalars)
!
!     ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
!     REMARKS: All angles are in RADIANS
!
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  May 11, 1992.
!     REVISIONS: May 26 - now use better Kepler solver for ellipses
!                 and hyperbolae called EHYBRID.F and FHYBRID.F
!**********************************************************************
subroutine orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm, &
        x, y, z, vx, vy, vz)
    include 'rvmod.inc'

    !...  Inputs Only:
    integer ::  ialpha
    real(8) :: gm, a, e, inc, capom, omega, capm

    !...  Outputs:
    real(8) :: x, y, z, vx, vy, vz

    !...  Internals:
    real(8) :: cape, capf, zpara, em1
    real(8) :: sp, cp, so, co, si, ci
    real(8) :: d11, d12, d13, d21, d22, d23
    real(8) :: scap, ccap, shcap, chcap
    real(8) :: sqe, sqgma, xfac1, xfac2, ri, vfac1, vfac2
    real(8) :: orbel_ehybrid, orbel_fhybrid, orbel_zget

    vfac1 = 0
    vfac2 = 0
    xfac1 = 0
    xfac2 = 0

    !...  Executable code
    if(e<0.0) then
        !write(*,*) ' ERROR in orbel_el2xv: e<0, setting e=0!!!'
        e = 0.0
        !ialpha = -1
    endif

    !...    check for inconsistencies between ialpha and e
    em1 = e - 1.d0
    if(&
            ((ialpha.eq.0) .and. (abs(em1)>TINY))  .or.&
                    ((ialpha<0) .and. (e>1.0d0))  .or.&
                    ((ialpha>0) .and. (e<1.0d0)))  then
        e = 0.0
        ialpha = -1

        !write(*,*) 'ERROR in orbel_el2xv: ialpha and e inconsistent'
        !write(*,*) '                       ialpha = ',ialpha
        !write(*,*) '                            e = ',e
    endif

    ! Generate rotation matrices (on p. 42 of Fitzpatrick)
    call orbel_scget(omega, sp, cp)
    call orbel_scget(capom, so, co)
    call orbel_scget(inc, si, ci)
    d11 = cp * co - sp * so * ci
    d12 = cp * so + sp * co * ci
    d13 = sp * si
    d21 = -sp * co - cp * so * ci
    d22 = -sp * so + cp * co * ci
    d23 = cp * si

    ! Get the other quantities depending on orbit type ( i.e. IALPHA)
    if (ialpha .eq. -1) then
        cape = orbel_ehybrid(e, capm)
        call orbel_scget(cape, scap, ccap)
        sqe = dsqrt(1.d0 - e * e)
        sqgma = dsqrt(gm * a)
        xfac1 = a * (ccap - e)
        xfac2 = a * sqe * scap
        ri = 1.d0 / (a * (1.d0 - e * ccap))
        vfac1 = -ri * sqgma * scap
        vfac2 = ri * sqgma * sqe * ccap
    endif

    if (ialpha .eq. +1) then
        capf = orbel_fhybrid(e, capm)
        call orbel_schget(capf, shcap, chcap)
        sqe = dsqrt(e * e - 1.d0)
        sqgma = dsqrt(gm * a)
        xfac1 = a * (e - chcap)
        xfac2 = a * sqe * shcap
        ri = 1.d0 / (a * (e * chcap - 1.d0))
        vfac1 = -ri * sqgma * shcap
        vfac2 = ri * sqgma * sqe * chcap
    endif

    if (ialpha .eq. 0) then
        zpara = orbel_zget(capm)
        sqgma = dsqrt(2.d0 * gm * a)
        xfac1 = a * (1.d0 - zpara * zpara)
        xfac2 = 2.d0 * a * zpara
        ri = 1.d0 / (a * (1.d0 + zpara * zpara))
        vfac1 = -ri * sqgma * zpara
        vfac2 = ri * sqgma
    endif

    x = d11 * xfac1 + d21 * xfac2
    y = d12 * xfac1 + d22 * xfac2
    z = d13 * xfac1 + d23 * xfac2
    vx = d11 * vfac1 + d21 * vfac2
    vy = d12 * vfac1 + d22 * vfac2
    vz = d13 * vfac1 + d23 * vfac2
    return
end

!**********************************************************************
!                    ORBEL_EHYBRID.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: For e < 0.18 uses fast routine ESOLMD
!                 For larger e but less than 0.8, uses EGET
!                 For e > 0.8 uses EHIE
!     REMARKS: Only EHIE brings M and E into range (0,TWOPI)
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 25,1992.
!     REVISIONS: 2/26/93 hfl
!**********************************************************************
real(8) function orbel_ehybrid(e, m)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) e, m

    !...  Internals:
    real(8) orbel_esolmd, orbel_eget, orbel_ehie

    !...  Executable code
    if(e < 0.18d0) then
        orbel_ehybrid = orbel_esolmd(e, m)
    else
        if(e <= 0.8d0) then
            orbel_ehybrid = orbel_eget(e, m)
        else
            orbel_ehybrid = orbel_ehie(e, m)
        endif
    endif

    return
end

!**********************************************************************
!                    ORBEL_EHIE.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Use Danby's quartic for 3 iterations.
!                Eqn. is f(x) = x - e*sin(x+M). Note  that
!                 E = x + M. First guess is very good for e near 1.
!                 Need to first get M between 0. and PI and use
!                 symmetry to return right answer if M between PI and 2PI
!     REMARKS: Modifies M so that both E and M are in range (0,TWOPI)
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 25,1992.
!     REVISIONS:
!**********************************************************************
real(8) function orbel_ehie(e, m)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: e, m

    !...  Internals:
    integer ::  iflag, nper, niter, NMAX
    real(8) :: dx, x, sa, ca, esa, eca, f, fp

    parameter (NMAX = 3)
    !...  Executable code
    ! In this section, bring M into the range (0,TWOPI) and if
    ! the result is greater than PI, solve for (TWOPI - M).
    iflag = 0
    nper = int(m / TWOPI)
    m = m - nper * TWOPI
    if (m < 0.d0) m = m + TWOPI

    if (m>PI) then
        m = TWOPI - m
        iflag = 1
    endif

    ! Make a first guess that works well for e near 1.
    x = (6.d0 * m)**(1.d0 / 3.d0) - m
    niter = 0

    ! Iteration loop
    do niter = 1, NMAX
        call orbel_scget(x + m, sa, ca)
        esa = e * sa
        eca = e * ca
        f = x - esa
        fp = 1.d0 - eca
        dx = -f / fp
        dx = -f / (fp + 0.5d0 * dx * esa)
        dx = -f / (fp + 0.5d0 * dx * (esa + 0.3333333333333333d0 * eca * dx))
        x = x + dx
    enddo

    orbel_ehie = m + x

    if (iflag.eq.1) then
        orbel_ehie = TWOPI - orbel_ehie
        m = TWOPI - m
    endif
    return
end

!------------------------------------------------------------------
!**********************************************************************
!                    ORBEL_FHYBRID.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           n ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON
!               For larger N, uses FGET
!     REMARKS:
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 26,1992.
!     REVISIONS:
!     REVISIONS: 2/26/93 hfl
!**********************************************************************
real(8) function orbel_fhybrid(e, n)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: e, n
    !...  Internals:
    real(8) :: abn
    real(8) :: orbel_flon, orbel_fget

    !...  Executable code
    abn = n
    if(n<0.d0) abn = -abn
    if(abn < 0.636d0 * e - 0.6d0) then
        orbel_fhybrid = orbel_flon(e, n)
    else
        orbel_fhybrid = orbel_fget(e, n)
    endif
    return
end

!**********************************************************************
!                          ORBEL_SCGET.F
!**********************************************************************
!     PURPOSE:  Given an angle, efficiently compute sin and cos.
!
!        Input:
!             angle ==> angle in radians (real scalar)
!
!        Output:
!             sx    ==>  sin(angle)  (real scalar)
!             cx    ==>  cos(angle)  (real scalar)
!
!     ALGORITHM: Obvious from the code
!     REMARKS: The HP 700 series won't return correct answers for sin
!       and cos if the angle is bigger than 3e7. We first reduce it
!       to the range [0,2pi) and use the sqrt rather than cos (it's faster)
!       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  May 6, 1992.
!     REVISIONS:
!**********************************************************************
subroutine orbel_scget(angle, sx, cx)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: angle

    !...  Output:
    real(8) :: sx, cx

    !... Internals:
    integer ::  nper
    real(8) :: x
    real(8) :: PI3BY2
    parameter(PI3BY2 = 1.5d0 * PI)

    !...  Executable code
    nper = int(angle / TWOPI)
    x = angle - nper * TWOPI
    if(x<0.d0) then
        x = x + TWOPI
    endif
    sx = sin(x)
    cx = sqrt(1.d0 - sx * sx)
    if((x > PIBY2) .and. (x <PI3BY2)) then
        cx = -cx
    endif
    return
end

!**********************************************************************
!                    ORBEL_EGET.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!                  orbel_eget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Quartic convergence from Danby
!     REMARKS: For results very near roundoff, give it M between
!           0 and 2*pi. One can condition M before calling EGET
!           by calling my double precision function MOD2PI(M).
!           This is not done within the routine to speed it up
!           and because it works fine even for large M.
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 7, 1992.
!     REVISIONS: May 21, 1992.  Now have it go through EXACTLY two iterations
!                with the premise that it will only be called if
!                 we have an ellipse with e between 0.15 and 0.8
!**********************************************************************
real(8) function orbel_eget(e, m)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: e, m

    !...  Internals:
    real(8) :: x, sm, cm, sx, cx
    real(8) :: es, ec, f, fp, fpp, fppp, dx

    !...  Executable code
    ! Function to solve Kepler's eqn for E (here called
    ! x) for given e and M. returns value of x.
    ! MAY 21 : FOR e < 0.18 use ESOLMD for speed and sufficient accuracy
    ! MAY 21 : FOR e > 0.8 use EHIE - this one may not converge fast enough.

    call orbel_scget(m, sm, cm)

    !  begin with a guess accurate to order ecc**3
    x = m + e * sm * (1.d0 + e * (cm + e * (1.d0 - 1.5d0 * sm * sm)))

    !  Go through one iteration for improved estimate
    call orbel_scget(x, sx, cx)
    es = e * sx
    ec = e * cx
    f = x - es - m
    fp = 1.d0 - ec
    fpp = es
    fppp = ec
    dx = -f / fp
    dx = -f / (fp + dx * fpp / 2.d0)
    dx = -f / (fp + dx * fpp / 2.d0 + dx * dx * fppp / 6.d0)
    orbel_eget = x + dx

    ! Do another iteration.
    ! For m between 0 and 2*pi this seems to be enough to
    ! get near roundoff error for eccentricities between 0 and 0.8
    x = orbel_eget
    call orbel_scget(x, sx, cx)
    es = e * sx
    ec = e * cx
    f = x - es - m
    fp = 1.d0 - ec
    fpp = es
    fppp = ec
    dx = -f / fp
    dx = -f / (fp + dx * fpp / 2.d0)
    dx = -f / (fp + dx * fpp / 2.d0 + dx * dx * fppp / 6.d0)

    orbel_eget = x + dx
    return
end

!**********************************************************************
!                    ORBEL_ESOLMD.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           m ==> mean anomaly. (real scalar)
!             Returns:
!                orbel_esolmd ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Some sort of quartic convergence from Wisdom.
!     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
!         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
!                ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI
!     INCLUDES: needs SCGET.F
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 7, 1992.
!     REVISIONS: 2/26/93 hfl
!**********************************************************************
real(8) function orbel_esolmd(e, m)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: e, m

    !...  Internals:
    real(8) :: x, sm, cm, sx, cx
    real(8) :: es, ec, f, fp, fpp, fppp, dx

    !...  Executable code
    !...    Function to solve Kepler's eqn for E (here called
    !...    x) for given e and M. returns value of x.
    call orbel_scget(m, sm, cm)
    x = m + e * sm * (1.d0 + e * (cm + e * (1.d0 - 1.5d0 * sm * sm)))

    call orbel_scget(x, sx, cx)
    es = e * sx
    ec = e * cx
    f = x - es - m
    fp = 1.d0 - ec
    fpp = es
    fppp = ec
    dx = -f / fp
    dx = -f / (fp + dx * fpp / 2.d0)
    dx = -f / (fp + dx * fpp / 2.d0 + dx * dx * fppp / 6.d0)

    orbel_esolmd = x + dx
    return
end

!***************************************************************************
!                  TU4_GETACCB.F
!*************************************************************************
! GETACCB returns the bary. acc. on each of n mutually
! interacting objects by direct pairwise summation
!
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of planets (real array)
!                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xb,yb,zb      ==>  position of planets in beri coord
!                                    (real arrays)
!             Output:
!               axb,ayb,azb   ==>  accel in beri coord (real arrays)
!
! Remarks:  Based on Martin's NB4M routines
! Authors:  Martin Duncan
! Date:    3/8/93
! Last revision: 4/5/95
subroutine tu4_getaccb(nbod, mass, j2rp2, j4rp4, xb, yb, zb, axb, ayb, azb)
    include 'rvmod.inc'

    !...  Inputs Only:
    integer ::  nbod
    real(8) :: mass(nbod), j2rp2, j4rp4
    real(8) :: xb(nbod), yb(nbod), zb(nbod)

    !...  Output
    real(8) :: axb(nbod), ayb(nbod), azb(nbod)

    !...  Internals
    real(8) :: xx, yy, zz, rr2, fac, mi, mj
    real(8) :: axx, ayy, azz, fac1
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX), irh(NPLMAX)
    real(8) :: aoblx(NPLMAX), aobly(NPLMAX), aoblz(NPLMAX)
    integer ::  i, j

    !...  executable code
    xh(1) = 0.0d0
    yh(1) = 0.0d0
    zh(1) = 0.0d0

    !...  do sun part first
    axb(1) = 0.0d0
    ayb(1) = 0.0d0
    azb(1) = 0.0d0
    i = 1
    do j = i + 1, nbod
        mi = mass(i)
        mj = mass(j)
        xx = xb(i) - xb(j)
        yy = yb(i) - yb(j)
        zz = zb(i) - zb(j)
        rr2 = xx**2 + yy**2 + zz**2

        fac1 = 1.d0 / dsqrt(rr2)
        fac = fac1 / rr2

        !..      save for the J2 and J4 calculations
        xh(j) = -xx
        yh(j) = -yy
        zh(j) = -zz
        irh(j) = fac1

        axx = xx * fac
        ayy = yy * fac
        azz = zz * fac
        axb(i) = axb(i) - axx * mj
        ayb(i) = ayb(i) - ayy * mj
        azb(i) = azb(i) - azz * mj
        axb(j) = axx * mi
        ayb(j) = ayy * mi
        azb(j) = azz * mi
    enddo

    do i = 2, nbod - 1
        do j = i + 1, nbod
            mi = mass(i)
            mj = mass(j)
            xx = xb(i) - xb(j)
            yy = yb(i) - yb(j)
            zz = zb(i) - zb(j)
            rr2 = xx**2 + yy**2 + zz**2
            fac = 1.d0 / (rr2 * dsqrt(rr2))
            axx = xx * fac
            ayy = yy * fac
            azz = zz * fac
            axb(i) = axb(i) - axx * mj
            ayb(i) = ayb(i) - ayy * mj
            azb(i) = azb(i) - azz * mj
            axb(j) = axb(j) + axx * mi
            ayb(j) = ayb(j) + ayy * mi
            azb(j) = azb(j) + azz * mi
        enddo
    enddo

    if(j2rp2.ne.0.0d0) then
        call obl_acc(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, &
                aoblx, aobly, aoblz)
        do i = 1, nbod
            axb(i) = axb(i) + aoblx(i)
            ayb(i) = ayb(i) + aobly(i)
            azb(i) = azb(i) + aoblz(i)
        enddo
    endif
    return
end

!                          COORD_H2B.F
!***********************************************************************
!     PURPOSE: Converts from Heliocentric to Barycentric coords.
!     ARGUMENTS:  Input is
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!                   mass(*) ==>  masses (real array)
!                 xh(*),yh(*),zh(*) ==> heliocentric particle coords
!                                          (real array)
!                 vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
!                                             (real array)
!                 Returned are
!                    xb(*),yb(*),zb(*) ==> bary. particle positions
!                                          (real array)
!                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
!                                            (real array)
!                    msys              ==>  Total mass of of system
!                                            (real scalar)
!     Authors:  Martin Duncan
!     ALGORITHM: Obvious
!     WRITTEN:  Jan 27/93
!     REVISIONS: 2/22/94  HFL
subroutine coord_h2b(nbod, mass, xh, yh, zh, vxh, vyh, vzh, &
        xb, yb, zb, vxb, vyb, vzb, msys)
    include 'rvmod.inc'

    !...  Inputs:
    integer ::  nbod
    real(8) :: mass(NPLMAX)
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) :: vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)

    !...  Outputs:
    real(8) :: xb(NPLMAX), yb(NPLMAX), zb(NPLMAX)
    real(8) :: vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)

    !...  Internals:
    real(8) :: msys, xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp
    integer ::  n

    !...  Executable code
    msys = mass(1)
    xtmp = 0.d0
    ytmp = 0.d0
    ztmp = 0.d0
    vxtmp = 0.d0
    vytmp = 0.d0
    vztmp = 0.d0

    do n = 2, nbod
        msys = msys + mass(n)
        xtmp = xtmp + mass(n) * xh(n)
        ytmp = ytmp + mass(n) * yh(n)
        ztmp = ztmp + mass(n) * zh(n)
        vxtmp = vxtmp + mass(n) * vxh(n)
        vytmp = vytmp + mass(n) * vyh(n)
        vztmp = vztmp + mass(n) * vzh(n)
    enddo

    xb(1) = -xtmp / msys
    yb(1) = -ytmp / msys
    zb(1) = -ztmp / msys
    vxb(1) = -vxtmp / msys
    vyb(1) = -vytmp / msys
    vzb(1) = -vztmp / msys

    do n = 2, nbod
        xb(n) = xh(n) + xb(1)
        yb(n) = yh(n) + yb(1)
        zb(n) = zh(n) + zb(1)
        vxb(n) = vxh(n) + vxb(1)
        vyb(n) = vyh(n) + vyb(1)
        vzb(n) = vzh(n) + vzb(1)
    enddo
    return
end

!***********************************************************************
!                          COORD_B2H.F
!***********************************************************************
!     PURPOSE: Converts from Barycentric to Helio coords.
!     ARGUMENTS:  Input is
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!                   mass(*) ==>  masses (real array)
!                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
!                                 SYMMETRY IN SUBROUTINE CALLS
!                 xb(*),yb(*),zb(*) ==> Barycentric particle coords
!                                          (real array)
!                 vxb(*),vyb(*),vzb(*) ==> Barycentric particle velocities
!                                             (real array)
!                 Returned are
!                    xh(*),yh(*),zh(*) ==> Helio particle positions
!                                          (real array)
!                    vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
!                                            (real array)
!
!     ALGORITHM: Obvious
!     REMARKS:  Can of course use this to get coords. relative to any body.
!              by changing the one subtracted off.
!
!     Authors:  Martin Duncan
!     WRITTEN:  Jan 27/93
!     REVISIONS: 2/17/95  HFL
subroutine coord_b2h(nbod, xb, yb, zb, vxb, vyb, vzb, &
        xh, yh, zh, vxh, vyh, vzh)
    include 'rvmod.inc'

    !...  Inputs:
    integer ::  nbod
    real(8) :: xb(NPLMAX), yb(NPLMAX), zb(NPLMAX)
    real(8) :: vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)

    !...  Outputs:
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) :: vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)
    
   

    !...  Internals:
    integer ::  n

    !...  Executable code
    do n = 1, nbod
        xh(n) = xb(n) - xb(1)
        yh(n) = yb(n) - yb(1)
        zh(n) = zb(n) - zb(1)
        vxh(n) = vxb(n) - vxb(1)
        vyh(n) = vyb(n) - vyb(1)
        vzh(n) = vzb(n) - vzb(1)
    enddo
    return
end

!***************************************************************************
!                  OBL_ACC.F
!*************************************************************************
! OBL_ACC returns the BARYCENTRIC x,y,z components of the acc. on NBOD
! particles due to the oblateness of mass(1) using
! the values of J2RP2 and J4RP4 passed into the routine.
! (J2RP2 for example is the product of
! J_2 times the square of the central body's radius)
! Here we return the net acc. produced
! only by the J2 and J4 terms (i.e. including
! neither the monopole nor higher order terms).
!
!
!             Input:
!                 nbod     ==>  number of massive bodies (incl. central one)
!                 mass(*)  ==>  masses of particles (real(8) array)
!                 j2rp2    ==>  scaled value of j2 moment (real(8) scalar)
!                 j4rp4    ==>  scaled value of j4 moment (real(8) scalar)
!                                    (real(8) vectors)
!                 xh(*),yh(*),zh(*)   ==>  HELIO. positions of particles
!                 irh(*)   ==> 1./ magnitude of radius vector (real(8) vector)
!                                (passed in to save calcs.)
!             Output:
!               aoblx(*),aobly(*),aoblz(*)  ==>  BARY. components of accel
!                                        (real(8) vectors)
!
! Remarks:  aoblx(1) (for example) contains x-component of
!           bary. acc. of central body
! Authors:  Martin Duncan
! Date:    3/4/94
! Last revision:
subroutine obl_acc(nbod, mass, j2rp2, j4rp4, xh, yh, zh, irh, &
        aoblx, aobly, aoblz)
    include 'rvmod.inc'

    !...  Inputs Only:
    integer ::  nbod
    real(8) :: j2rp2, j4rp4
    real(8) :: mass(NPLMAX)
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX), irh(NPLMAX)

    !...  Output
    real(8) :: aoblx(NPLMAX), aobly(NPLMAX), aoblz(NPLMAX)

    !...  Internals
    integer ::  n
    real(8) :: rinv2, t0, t1, t2, t3
    real(8) :: fac1, fac2

    !...  executable code
    ! First get the bary acc. of each "planet" due to central oblate "sun"
    do n = 2, nbod
        ! Note that here we assume we know inverse of radius rather than calc. it
        ! from (x,y,z) to save the ddsqrt.
        rinv2 = irh(n)**2
        t0 = -mass(1) * rinv2 * rinv2 * irh(n)
        t1 = 1.5d0 * j2rp2
        t2 = zh(n) * zh(n) * rinv2
        t3 = 1.875d0 * j4rp4 * rinv2

        fac1 = t0 * (t1 - t3 - (5.d0 * t1 - (14.d0 - 21.d0 * t2) * t3) * t2)
        fac2 = 2.d0 * t0 * (t1 - (2.d0 - (14.d0 * t2 / 3.d0)) * t3)

        aoblx(n) = fac1 * xh(n)
        aobly(n) = fac1 * yh(n)
        aoblz(n) = (fac1 + fac2) * zh(n)
    enddo
    ! Now compute the bary. acc. of Sun due to all the planets
    aoblx(1) = 0.d0
    aobly(1) = 0.d0
    aoblz(1) = 0.d0
    do n = 2, nbod
        aoblx(1) = aoblx(1) - mass(n) * aoblx(n) / mass(1)
        aobly(1) = aobly(1) - mass(n) * aobly(n) / mass(1)
        aoblz(1) = aoblz(1) - mass(n) * aoblz(n) / mass(1)
    enddo
    return
end

!**********************************************************************
!                    ORBEL_FGET.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_fget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
!           Cel. Mech. ".  Quartic convergence from Danby's book.
!     REMARKS:
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: 2/26/93 hfl
!**********************************************************************
real(8) function orbel_fget(e, capn)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: e, capn

    !...  Internals:
    integer ::  i, IMAX
    real(8) :: tmp, x, shx, chx
    real(8) :: esh, ech, f, fp, fpp, fppp, dx
    PARAMETER (IMAX = 10)

    !...  Executable code
    ! Function to solve "Kepler's eqn" for F (here called
    ! x) for given e and CAPN.
    !  begin with a guess proposed by Danby
    if(capn < 0.d0) then
        tmp = -2.d0 * capn / e + 1.8d0
        x = -log(tmp)
    else
        tmp = +2.d0 * capn / e + 1.8d0
        x = log(tmp)
    endif

    orbel_fget = x

    do i = 1, IMAX
        call orbel_schget(x, shx, chx)
        esh = e * shx
        ech = e * chx
        f = esh - x - capn
        !        write(6,*) 'i,x,f : ',i,x,f
        fp = ech - 1.d0
        fpp = esh
        fppp = ech
        dx = -f / fp
        dx = -f / (fp + dx * fpp / 2.d0)
        dx = -f / (fp + dx * fpp / 2.d0 + dx * dx * fppp / 6.d0)
        orbel_fget = x + dx
        !   If we have converged here there's no point in going on
        if(abs(dx) <= TINY) RETURN
        x = orbel_fget
    enddo
    write(6, *) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE'
    return
end

!**********************************************************************
!                    ORBEL_ZGET.F
!**********************************************************************
!     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola
!          given Q (Fitz. notation.)
!
!             Input:
!                           q ==>  parabola mean anomaly. (real scalar)
!             Returns:
!                  orbel_zget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
!     REMARKS: For a parabola we can solve analytically.
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: May 27 - corrected it for negative Q and use power
!            series for small Q.
!**********************************************************************
real(8) function orbel_zget(q)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: q

    !...  Internals:
    integer ::  iflag
    real(8) :: x, tmp

    !...  Executable code
    iflag = 0
    if(q<0.d0) then
        iflag = 1
        q = -q
    endif

    if (q<1.d-3) then
        orbel_zget = q * (1.d0 - (q * q / 3.d0) * (1.d0 - q * q))
    else
        x = 0.5d0 * (3.d0 * q + dsqrt(9.d0 * (q**2) + 4.d0))
        tmp = x**(1.d0 / 3.d0)
        orbel_zget = tmp - 1.d0 / tmp
    endif

    if(iflag .eq.1) then
        orbel_zget = -orbel_zget
        q = -q
    endif
    return
end

!**********************************************************************
!                        ORBEL_SCHGET.F
!**********************************************************************
!     PURPOSE:  Given an angle, efficiently compute sinh and cosh.
!
!        Input:
!             angle ==> angle in radians (real scalar)
!
!        Output:
!             shx    ==>  sinh(angle)  (real scalar)
!             chx    ==>  cosh(angle)  (real scalar)
!
!     ALGORITHM: Obvious from the code
!     REMARKS: Based on the routine SCGET for sine's and cosine's.
!       We use the dsqrt rather than cosh (it's faster)
!       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
!       OR OVERFLOWS WILL OCCUR!
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  May 6, 1992.
!     REVISIONS:
!**********************************************************************
subroutine orbel_schget(angle, shx, chx)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: angle

    !...  Output:
    real(8) :: shx, chx

    !----
    !...  Executable code

    shx = sinh(angle)
    chx = dsqrt(1.d0 + shx * shx)
    return
end

!**********************************************************************
!                    ORBEL_FLON.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_flon ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Uses power series for N in terms of F and Newton,s method
!     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
!     AUTHOR: M. Duncan
!     DATE WRITTEN: May 26, 1992.
!     REVISIONS:
!**********************************************************************
real(8) function orbel_flon(e, capn)
    include 'rvmod.inc'

    !...  Inputs Only:
    real(8) :: e, capn

    !...  Internals:
    integer ::  iflag, i, IMAX
    real(8) :: a, b, sq, biga, bigb
    real(8) :: x, x2
    real(8) :: f, fp, dx
    real(8) :: diff
    real(8) :: a0, a1, a3, a5, a7, a9, a11
    real(8) :: b1, b3, b5, b7, b9, b11
    PARAMETER (IMAX = 10)
    PARAMETER (a11 = 156.d0, a9 = 17160.d0, a7 = 1235520.d0)
    PARAMETER (a5 = 51891840.d0, a3 = 1037836800.d0)
    PARAMETER (b11 = 11.d0 * a11, b9 = 9.d0 * a9, b7 = 7.d0 * a7)
    PARAMETER (b5 = 5.d0 * a5, b3 = 3.d0 * a3)

    !...  Executable code
    ! Function to solve "Kepler's eqn" for F (here called
    ! x) for given e and CAPN. Only good for smallish CAPN
    iflag = 0
    if(capn < 0.d0) then
        iflag = 1
        capn = -capn
    endif

    a1 = 6227020800.d0 * (1.d0 - 1.d0 / e)
    a0 = -6227020800.d0 * capn / e
    b1 = a1

    !  Set iflag nonzero if capn < 0., in which case solve for -capn
    !  and change the sign of the final answer for F.
    !  Begin with a reasonable guess based on solving the cubic for small F
    a = 6.d0 * (e - 1.d0) / e
    b = -6.d0 * capn / e
    sq = dsqrt(0.25 * b * b + a * a * a / 27.d0)
    biga = (-0.5 * b + sq)**0.3333333333333333d0
    bigb = -(+0.5 * b + sq)**0.3333333333333333d0
    x = biga + bigb
    orbel_flon = x
    ! If capn is tiny (or zero) no need to go further than cubic even for
    ! e =1.
    if(capn < TINY) go to 100

    do i = 1, IMAX
        x2 = x * x
        f = a0 + x * (a1 + x2 * (a3 + x2 * (a5 + x2 * (a7 + x2 * (a9 + x2 * (a11 + x2))))))
        fp = b1 + x2 * (b3 + x2 * (b5 + x2 * (b7 + x2 * (b9 + x2 * (b11 + 13.d0 * x2)))))
        dx = -f / fp
        !        write(6,*) 'i,dx,x,f : '
        !        write(6,432) i,dx,x,f
        !432        format(1x,i3,3(2x,1p1e22.15))
        orbel_flon = x + dx
        !   If we have converged here there's no point in going on
        if(abs(dx) <= TINY) go to 100
        x = orbel_flon
    enddo

    ! Abnormal return here - we've gone thru the loop
    ! IMAX times without convergence
    if(iflag .eq. 1) then
        orbel_flon = -orbel_flon
        capn = -capn
    endif
    write(6, *) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE'
    diff = e * sinh(orbel_flon) - orbel_flon - capn
    write(6, *) 'N, F, ecc*sinh(F) - F - N : '
    write(6, *) capn, orbel_flon, diff
    return

    !  Normal return here, but check if capn was originally negative
100 if(iflag .eq. 1) then
        orbel_flon = -orbel_flon
        capn = -capn
    endif
    return
end

!***********************************************************************
!                          COORD_J2H.F
!***********************************************************************
!     PURPOSE: Converts from Jacobi to Helio coords.
!     ARGUMENTS:  Input is
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!                   mass(*) ==>  masses (real array)
!                 xj(*),yj(*),zj(*) ==> Jacobi particle coords
!                                          (real array)
!                 vxj(*),vyj(*),vzj(*) ==> Jacobi particle velocities
!                                             (real array)
!                 Returned are
!                    xh(*),yh(*),zh(*) ==> Helio particle positions
!                                          (real array)
!                    vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
!                                            (real array)
!
!     ALGORITHM: See my notes on Nov 21.
!
!     Authors:  Martin Duncan
!     WRITTEN:  Jan 27/93
!     REVISIONS: 2/17/95  HFL
subroutine coord_j2h(nbod, mass, xj, yj, zj, vxj, vyj, vzj, &
        xh, yh, zh, vxh, vyh, vzh)
    include 'rvmod.inc'

    !...  Inputs:
    integer ::  nbod
    real(8) :: mass(NPLMAX)
    real(8) :: xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
    real(8) :: vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)

    !...  Outputs:
    real(8) :: xh(NPLMAX), yh(NPLMAX), zh(NPLMAX)
    real(8) :: vxh(NPLMAX), vyh(NPLMAX), vzh(NPLMAX)

    !...  Internals:
    integer ::  n
    real(8) :: sumx, sumy, sumz, sumvx, sumvy, sumvz
    real(8) :: eta(NPLMAX)

    !...  Executable code
    ! First calc. the array eta(*) then convert to jacobi coords
    eta(1) = mass(1)
    do n = 2, nbod
        eta(n) = eta(n - 1) + mass(n)
    enddo

    xh(1) = 0.d0
    yh(1) = 0.d0
    zh(1) = 0.d0
    vxh(1) = 0.d0
    vyh(1) = 0.d0
    vzh(1) = 0.d0

    xh(2) = xj(2)
    yh(2) = yj(2)
    zh(2) = zj(2)
    vxh(2) = vxj(2)
    vyh(2) = vyj(2)
    vzh(2) = vzj(2)

    sumx = mass(2) * xj(2) / eta(2)
    sumy = mass(2) * yj(2) / eta(2)
    sumz = mass(2) * zj(2) / eta(2)
    sumvx = mass(2) * vxj(2) / eta(2)
    sumvy = mass(2) * vyj(2) / eta(2)
    sumvz = mass(2) * vzj(2) / eta(2)

    do n = 3, nbod
        xh(n) = xj(n) + sumx
        yh(n) = yj(n) + sumy
        zh(n) = zj(n) + sumz
        vxh(n) = vxj(n) + sumvx
        vyh(n) = vyj(n) + sumvy
        vzh(n) = vzj(n) + sumvz

        if(n<nbod) then
            sumx = sumx + mass(n) * xj(n) / eta(n)
            sumy = sumy + mass(n) * yj(n) / eta(n)
            sumz = sumz + mass(n) * zj(n) / eta(n)
            sumvx = sumvx + mass(n) * vxj(n) / eta(n)
            sumvy = sumvy + mass(n) * vyj(n) / eta(n)
            sumvz = sumvz + mass(n) * vzj(n) / eta(n)
        endif
    enddo
    return
end

!**************************************************************************
!                              BS_INT
!**************************************************************************
! This is the subroutine that does the the bs integration step.
!
!             Input:
!              nbod          ==> number of planets  (int scalar)
!              ntp           ==> number of test particles  (int scalar)
!              mass          ==>  mass of bodies (real array)
!              j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!              istat         ==>  status of the test paricles
!                                       (2d integer array)
!                                   istat(i,1) = 0 ==> active:  = 1 not
!                                   istat(i,2) = -1 ==> Danby did not work
!               x            ==> initial value independent variable
!                                        (real scalar)
!               h0           ==> stepsize  (real scalar)
!               y            ==> initial value dependent variables
!                                     (real array)
!               eps          ==> local truncation error tolerance
!
!             Output:
!                 y  ==> final value dependent variables  (real array)
!                 x  ==> final value independent variable (real scalar)
!                 h0 ==> recommended stepsize for the next call (real scalar)
!
! Remarks:  Based on Renu's code: mass,j2rp2,j4rp4 and istat are
!           just passed on to bs_der
! Authors:  Hal Levison
! Date:    5/17/93
! Last revision: 2/24/94
subroutine bs_int_pl(nbod, ntp, mass, j2rp2, j4rp4, x, h0, y, eps)
    include 'rvmod.inc'
    include 'bs.inc'

    !...  Inputs Only:
    integer ::  nbod, ntp
    real(8) :: mass(nbod), h0, eps, j2rp2, j4rp4

    !...  Input & Output
    real(8) :: x, y(N6DBS)

    !...  Internals
    real(8) :: tp(NTEMP), dy(N6DBS), d(6), alt(10), lt(10)
    integer ::  idiv, i, ii, m, l, m1, k, mmk, i1, i1max, ik, n
    real(8) :: xa, xb, varm, fl, h, hd, flt, eta2, dta, yb, c, b1
    real(8) :: den, dtn, b, var, varma
    logical lbreak

    data lt/1, 2, 3, 4, 6, 8, 12, 16, 24, 32/
    data alt/1.d0, 2.d0, 3.d0, 4.d0, 6.d0, 8.d0, 12.d0, 16.d0, 24.d0, 32.d0/

    save lt, alt
    d(:) = 0.d0
    
    varma = 0.d0

    !...  Executable code
    n = 6 * (nbod + ntp)
    xa = x

    call bs_der_pl(nbod, mass, j2rp2, j4rp4, y, dy)
    do i = 1, n
        ii = 12 * i
        tp(ii - 1) = dabs(y(i))
        !
        if(tp(ii - 1)<eps) then
            tp(ii - 1) = eps
        endif
        !
        tp(ii - 4) = dy(i)
        tp(ii) = y(i)
    enddo

    do idiv = 0, NTRYS
        xb = h0 + xa
        !successive extrapolations
        m = 1
        lbreak = .true.
        do while((m<=10) .and. lbreak)
            l = int(lt(m))
            fl = alt(m)
            varm = 0.d0
            m1 = min0(m - 1, 6)
            !calculation of d(k)=(h(m-k)/h(m))**2 for equation (6)
            if(m1.ne.0) then
                do k = 1, m1
                    mmk = m - k
                    flt = alt(mmk)
                    d(k) = (fl / flt)**2
                enddo
            endif
            h = h0 / fl
            hd = 0.5d0 * h
            !integration
            do i = 1, n
                ii = 12 * i
                tp(ii - 3) = tp(ii)
                y(i) = tp(ii) + hd * tp(ii - 4) !equation (3b)
            enddo
            i1max = 2 * l - 1
            x = xa

            do i1 = 1, i1max
                x = x + hd
                call bs_der_pl(nbod, mass, j2rp2, j4rp4, y, dy)
                do i = 1, n
                    ii = 12 * i
                    tp(ii - 1) = dmax1(tp(ii - 1), dabs(y(i)))
                    eta2 = tp(ii - 3) + h * dy(i)
                    tp(ii - 3) = y(i)
                    y(i) = eta2
                enddo
            enddo

            call bs_der_pl(nbod, mass, j2rp2, j4rp4, y, dy)
            do i = 1, n
                ii = 12 * i
                dta = tp(ii - 11)
                yb = (tp(ii - 3) + y(i) + hd * dy(i)) / 2.d0 !equation (3d)
                !extrapolated values
                c = yb !equation (6b)
                tp(ii - 11) = yb !equation (6a)

                if(m1.ne.0) then
                    do k = 1, m1
                        b1 = d(k) * dta
                        den = b1 - c
                        dtn = dta
                        if(den.ne.0.d0) then
                            b = (c - dta) / den
                            dtn = c * b !equation (6c)
                            c = b1 * b !equation (6d)
                        endif
                        ik = ii - 11 + k
                        dta = tp(ik)
                        tp(ik) = dtn !equation (6c)
                        yb = yb + dtn !equation (6f)
                    enddo
                    var = dabs(tp(ii - 2) - yb) / tp(ii - 1)
                    varm = dmax1(varm, var)
                endif
                tp(ii - 2) = yb
            enddo

            if(m>3) then
                if(varm<=eps) then !return results to calling program
                    x = xb
                    do i = 1, n
                        ii = 12 * i
                        y(i) = tp(ii - 2)
                    enddo
                    h0 = h0 * 1.5d0 * 0.6d0**(m - 1 - m1) !recommend a new step size
                    return
                endif

                if(varm>=varma) then !calculation did not converge
                    !start again with half the step size
                    lbreak = .false.
                endif
            endif
            varma = varm
            m = m + 1
        enddo
        h0 = h0 / 2.d0
    enddo
    write(*, *) ' ERROR (b_int): lack of convergence !!!! '
end

!**************************************************************************
!                              BS_DER
!**************************************************************************
! This is the subroutine that calculates the derivatives of the independant var
!
!             Input:
!              nbod  ==> number of planets  (int scalar)
!              ntp   ==> number of test particles  (int scalar)
!              mass  ==>  mass of bodies (real array)
!      j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!             istat  ==>  status of the test paricles
!                           (2d integer array)
!                            istat(i,1) = 0 ==> active:  = 1 not
!                            istat(i,2) = -1 ==> Danby did not work
!               ybs  ==> values dependent variables  (real array)
!
!             Output:
!                dy  ==> derivatives of the independant var (real array)
!
! Remarks:  This used TU4 routines !!
! Authors:  Hal Levison
! Date:    5/17/93
! Last revision: 2/24/94
subroutine bs_der_pl(nbod, mass, j2rp2, j4rp4, ybs, dy)
    include 'rvmod.inc'
    include 'bs.inc'

    !...  Inputs Only:
    integer ::  nbod
    real(8) :: mass(nbod), j2rp2, j4rp4
    real(8) :: ybs(6, (NTPMAX + NPLMAX))

    !...  Output
    real(8) :: dy(6, (NTPMAX + NPLMAX))

    !...  Internals
    integer ::  i
    real(8) :: xb(NPLMAX), yb(NPLMAX), zb(NPLMAX)
    real(8) :: vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)
    real(8) :: axb(NPLMAX), ayb(NPLMAX), azb(NPLMAX)

    !...  Executable code
    !...  move things so that I can deal with it
    do i = 1, nbod
        xb(i) = ybs(1, i)
        yb(i) = ybs(2, i)
        zb(i) = ybs(3, i)
        vxb(i) = ybs(4, i)
        vyb(i) = ybs(5, i)
        vzb(i) = ybs(6, i)
    enddo

    call tu4_getaccb(nbod, mass, j2rp2, j4rp4, xb, yb, zb, axb, ayb, azb)
    !.... moves things intp dy array
    do i = 1, nbod
        dy(1, i) = ybs(4, i)
        dy(2, i) = ybs(5, i)
        dy(3, i) = ybs(6, i)
        dy(4, i) = axb(i)
        dy(5, i) = ayb(i)
        dy(6, i) = azb(i)
    enddo
    return
end
