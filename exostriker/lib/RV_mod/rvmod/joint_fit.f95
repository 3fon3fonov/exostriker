subroutine joint_timing_likelihood(elements,planet_mass,stellar_mass,epoch,start_time,end_time,model_mode,input_axes, &
        step_days,tolerance,max_steps,jacobi_mode,samples_per_orbit,max_samples,max_events, &
        dataset,planet,number,observed,uncertainty,event_count,event_model, &
        data_model,residual,point_loglik,by_dataset,loglik,status,nobs,npl)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: model_mode
    real(8),intent(in) :: input_axes(npl)
    integer,intent(in) :: nobs,npl,max_steps,jacobi_mode,samples_per_orbit,max_samples,max_events
    integer,intent(in) :: dataset(nobs),planet(nobs),number(nobs)
    real(8),intent(in) :: elements(npl,6),planet_mass(npl),stellar_mass,epoch,start_time,end_time
    real(8),intent(in) :: step_days,tolerance,observed(nobs),uncertainty(nobs)
    integer,intent(out) :: event_count(npl),status
    real(8),intent(out) :: event_model(max_events,npl,4),data_model(nobs),residual(nobs)
    real(8),intent(out) :: point_loglik(nobs),by_dataset(10),loglik
    real(8),parameter :: logtwopi=1.837877066409345483560659472811d0
    ! A normalization scale only: timing must not assume a 1-Rsun star
    ! and reject valid ultra-short-period planets around small stars.
    real(8),parameter :: geometry_scale=1.d-6
    real(8),parameter :: radius_au=geometry_scale*6.957d8/1.49597892d11
    real(8),allocatable :: grid(:),geometry(:,:,:)
    real(8) :: axes(npl),cadence,scale,lo,hi,left,right,middle,h,g0,g1,g,root,clock_time
    real(8) :: r(3),v(3),compensation(10),corrected,total,span,ngrid_real
    integer :: ngrid,i,j,k,iteration,d,n,allocation_status
    status=7;event_count=0;event_model=0.d0;data_model=0.d0;residual=0.d0
    point_loglik=0.d0;by_dataset=0.d0;loglik=-huge(1.d0)
    if(npl<1 .or. npl>9 .or. nobs<0 .or. max_events<1 .or. max_samples<4) return
    if(samples_per_orbit<32 .or. max_steps<1 .or. jacobi_mode<0 .or. jacobi_mode>1) return
    if(.not.all(ieee_is_finite(elements)) .or. .not.all(ieee_is_finite(planet_mass))) return
    if(.not.all(ieee_is_finite([stellar_mass,epoch,start_time,end_time,step_days,tolerance]))) return
    if(stellar_mass<=0.d0 .or. step_days<=0.d0 .or. end_time<=start_time) return
    if(tolerance<1.d-15 .or. tolerance>1.d-3 .or. any(planet_mass<0.d0)) return
    if(any(elements(:,1)<=0.d0) .or. any(elements(:,2)<0.d0) .or. any(elements(:,2)>=1.d0)) return
    if(any(elements(:,5)<0.d0) .or. any(elements(:,5)>180.d0)) return
    if(any(dataset<1) .or. any(dataset>10) .or. any(planet<1) .or. any(planet>npl)) return
    if(any(number<1) .or. .not.all(ieee_is_finite(observed))) return
    if(.not.all(ieee_is_finite(uncertainty)) .or. any(uncertainty<=0.d0)) return
    scale=minval(elements(:,1)*(1.d0-elements(:,2))**1.5d0/sqrt(1.d0+elements(:,2)))
    cadence=scale/dble(samples_per_orbit)
    lo=start_time-epoch;hi=end_time-epoch;span=hi-lo
    ngrid_real=span/cadence+5.d0
    status=9
    if(.not.ieee_is_finite(ngrid_real)) return
    if(ngrid_real>max_samples .or. ngrid_real*npl*6.d0>64000000.d0) return
    ngrid=ceiling(span/cadence)+5
    allocate(grid(ngrid),geometry(ngrid,npl,6),stat=allocation_status)
    if(allocation_status/=0) return
    ! A sample on either side permits minima exactly at either bound.
    do i=1,ngrid
        grid(i)=lo-cadence+(span+2.d0*cadence)*dble(i-1)/dble(ngrid-1)
    enddo
    axes=input_axes/geometry_scale
    call rv_transit_geometry_core(grid,elements,planet_mass,axes,stellar_mass,geometry_scale, &
        model_mode,step_days,tolerance,max_steps,jacobi_mode,geometry,status,ngrid,npl)
    if(status/=0) return
    do j=1,npl
        do i=1,ngrid-1
            g0=dot_product(geometry(i,j,1:2),geometry(i,j,4:5))
            g1=dot_product(geometry(i+1,j,1:2),geometry(i+1,j,4:5))
            if(g0>0.d0 .or. g1<=0.d0) cycle
            h=grid(i+1)-grid(i);left=0.d0;right=1.d0
            do iteration=1,64
                middle=0.5d0*(left+right)
                call interpolate(middle,r,v)
                g=dot_product(r(1:2),v(1:2))
                if(g<=0.d0) then
                    left=middle
                else
                    right=middle
                endif
                if((right-left)*h<1.d-12) exit
            enddo
            root=grid(i)+0.5d0*(left+right)*h
            call interpolate(0.5d0*(left+right),r,v)
            if(r(3)<=0.d0 .or. root<lo-1.d-10 .or. root>hi+1.d-10) cycle
            root=max(lo,min(hi,root));n=event_count(j)+1
            status=10
            if(n>max_events) return
            event_count(j)=n;clock_time=epoch+root
            event_model(n,j,:)=[dble(n),clock_time, &
                sqrt(dot_product(r(1:2),r(1:2)))*radius_au, &
                sqrt(dot_product(v(1:2),v(1:2)))*radius_au]
        enddo
    enddo
    compensation=0.d0;status=8
    do k=1,nobs
        j=planet(k);n=number(k);d=dataset(k)
        if(n>event_count(j)) return
        data_model(k)=event_model(n,j,2)
        residual(k)=observed(k)-data_model(k)
        point_loglik(k)=-0.5d0*((residual(k)/uncertainty(k))**2+2.d0*log(uncertainty(k))+logtwopi)
        if(.not.ieee_is_finite(point_loglik(k))) return
        corrected=point_loglik(k)-compensation(d);total=by_dataset(d)+corrected
        compensation(d)=(total-by_dataset(d))-corrected;by_dataset(d)=total
    enddo
    loglik=sum(by_dataset)
    if(.not.ieee_is_finite(loglik)) then
        status=7;return
    endif
    status=0
contains
    subroutine interpolate(u,position,velocity)
        real(8),intent(in) :: u
        real(8),intent(out) :: position(3),velocity(3)
        real(8) :: u2,u3,p0(3),p1(3),v0(3),v1(3)
        u2=u*u;u3=u2*u
        p0=geometry(i,j,1:3);p1=geometry(i+1,j,1:3)
        v0=geometry(i,j,4:6);v1=geometry(i+1,j,4:6)
        position=(2*u3-3*u2+1)*p0+(u3-2*u2+u)*h*v0 &
            +(-2*u3+3*u2)*p1+(u3-u2)*h*v1
        velocity=((6*u2-6*u)*p0+(-6*u2+6*u)*p1)/h &
            +(3*u2-4*u+1)*v0+(3*u2-2*u)*v1
    end subroutine interpolate
end subroutine joint_timing_likelihood

! Native bounded Gaussian-likelihood optimization. No interpreter callbacks.
! Packed inputs: 10 orbital entries per planet (K,P,e,w,M,i,node,t0,aR,RpR),
! stellar mass/radius, RV trends, RV offset/jitter pairs, LC offset/jitter/LD4/
! linear/quadratic/dilution, and ten TTV rho/RpR pairs. Indices are one-based.
subroutine native_joint_fit(initial,free_index,lower,upper,scale,packed,map, &
        meta,data,sample_time,sample_owner,ld_laws,config,settings,method,maxiter,event_capacity, &
        parameters,models,variances,loglik,components,covariance,errors,diagnostics, &
        final_event_count,final_event_model,final_rv,final_flux,final_parts,final_ma,npar,nfree,npack,nobs,nsample,npl, &
        optimizer_options,prior_options)
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer,intent(in) :: npar,nfree,npack,nobs,nsample,npl,method,maxiter,event_capacity
    integer,intent(in) :: free_index(nfree),map(npack),meta(nobs,4),sample_owner(nsample),ld_laws(60),config(13)
    real(8),optional,intent(in) :: optimizer_options(9),prior_options(npar,6)
    real(8) :: controls(9),ftol,xtol,gtol,priors(npar,6)
    logical :: optimizing
    real(8),intent(in) :: initial(npar),lower(nfree),upper(nfree),scale(nfree),packed(npack)
    real(8),intent(in) :: data(nobs,4),sample_time(nsample),settings(8)
    real(8),intent(out) :: parameters(npar),models(nobs),variances(nobs),loglik,components(5)
    real(8),intent(out) :: covariance(nfree,nfree),errors(nfree)
    integer(8),intent(out) :: diagnostics(9)
    integer,intent(out) :: final_event_count(npl)
    real(8),intent(out) :: final_event_model(event_capacity,npl,4)
    real(8),intent(out) :: final_rv(nobs),final_flux(nobs),final_parts(nobs,npl),final_ma(npl)
    real(8) :: q(nfree),ql(nfree),qu(nfree),best(nfree),mu(nobs),var(nobs),f,bf,lambda, &
        grad(nfree),info(nfree,nfree),a(nfree,nfree),delta(nfree),candidate(nfree),cf, &
        cm(nobs),cv(nobs),dm(nobs,nfree),dv(nobs,nfree),h(nfree),qp(nfree),qm(nfree), &
        fp,fm,mp(nobs),mm(nobs),vp(nobs),vm(nobs),step,ratio,spread,diameter, &
        vertices(nfree,nfree+1),values(nfree+1),centre(nfree),xr(nfree),xe(nfree),fr,fe, &
        cp(5),nanval,threshold,inv(nfree,nfree),unit(nfree),z(nfree),origin(nfree)
    integer :: i,j,k,it,evals,code,bestcode,covcode,order(nfree+1),tmp,accepted,stagnant,light_calls,timing_calls
    integer :: shared_calls,trial_steps
    integer(8) :: shared_steps
    real(8) :: weight(nobs)
    real(8) :: cached_parts(nobs,npl)
    real(8) :: cached_rv(nobs),cached_flux(nobs),cached_ma(9),physical_key(10*npl+2), &
        light_key(60,4),timing_key(10*npl+2)
    real(8),allocatable :: cached_events(:,:,:)
    integer :: cached_counts(npl)
    logical :: ok,converged,light_valid,timing_valid,blocked(nfree)
    controls=[0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.001d0,1.d0,1.d-6]
    if(present(optimizer_options)) then
        where(optimizer_options>0.d0) controls=optimizer_options
    endif
    ftol=settings(7);if(controls(3)>0.d0) ftol=controls(3)
    xtol=settings(7);if(controls(4)>0.d0) xtol=controls(4)
    gtol=settings(7);if(controls(5)>0.d0) gtol=controls(5)
    priors=0.d0
    if(present(prior_options)) priors=prior_options
    optimizing=.true.
    nanval=ieee_value(0.d0,ieee_quiet_nan)
    parameters=initial;models=nanval;variances=nanval;loglik=-huge(1.d0);components=0.d0
    final_event_count=0;final_event_model=0.d0;final_rv=0.d0;final_flux=1.d0;final_parts=0.d0;final_ma=0.d0
    covariance=nanval;errors=nanval;diagnostics=[integer(8) :: 2,0,0,1,0,0,0,0,0]
    evals=0;light_calls=0;timing_calls=0;shared_calls=0;shared_steps=0
    if(nobs<1 .or. npl<0 .or. npl>9 .or. npar<1 .or. nfree<0) return
    if(npack/=10*npl+684 .or. method<0 .or. method>2 .or. maxiter<0) return
    if(any(free_index<1) .or. any(free_index>npar) .or. any(map<0) .or. any(map>npar)) return
    if(any(config(10:13)<0) .or. any(config(10:13)>1)) return
    if(any(config(1:7)<0) .or. any(config(1:7)>1)) return
    if(event_capacity/=config(8)) return
    if(config(8)<1 .or. config(8)>200000 .or. config(9)<1 .or. config(9)>60) return
    if(any(meta(:,1)<1) .or. any(meta(:,1)>5)) return
    if(any(meta(:,2)<1) .or. any(meta(:,2)>60)) return
    if(any(meta(:,3)<1) .or. any(meta(:,4)<1)) return
    if(any(meta(:,1)>2 .and. meta(:,3)>npl)) return
    if(any(meta(:,1)>2 .and. meta(:,2)>10)) return
    if(.not.all(ieee_is_finite(data)) .or. any(data(:,3)<=0.d0)) return
    if(.not.all(ieee_is_finite(settings)) .or. .not.all(ieee_is_finite(sample_time))) return
    if(settings(2)<=0.d0 .or. settings(3)<1.d-15 .or. settings(3)>1.d-3 .or. settings(7)<=0.d0) return
    if(any(meta(:,1)>2) .and. settings(5)<=settings(4)) return
    if(any(ieee_is_nan(lower)) .or. any(ieee_is_nan(upper))) return
    if(any(scale<=0.d0) .or. any(lower>upper)) return
    if(.not.all(ieee_is_finite(initial)) .or. .not.all(ieee_is_finite(scale))) return
    do i=1,nfree
        if(any(free_index(i+1:)==free_index(i))) return
    enddo
    if(.not.all(ieee_is_finite(priors))) return
    do i=1,npar
        if(priors(i,3)/=0.d0 .and. priors(i,2)<=0.d0) return
        if(priors(i,6)/=0.d0) then
            if(priors(i,4)<=0.d0 .or. priors(i,5)<=priors(i,4) .or. initial(i)<=0.d0) return
        endif
    enddo
    origin=initial(free_index);ql=(lower-origin)/scale;qu=(upper-origin)/scale
    if(any(ql>0.d0) .or. any(qu<0.d0)) return
    weight=1.d0
    do i=1,nobs
        if(meta(i,1)>=3) weight(i)=dble(config(meta(i,1)+7))
    enddo
    if(sum(weight)==0.d0) return
    light_valid=.false.;timing_valid=.false.;physical_key=0.d0;timing_key=0.d0;light_key=0.d0
    cached_rv=0.d0;cached_flux=1.d0;cached_parts=0.d0;cached_ma=0.d0
    if(any(meta(:,1)>2)) allocate(cached_events(config(8),npl,4))
    q=0.d0;best=q;bf=huge(1.d0)/100.d0;bestcode=0
    call objective(q,f,mu,var,code)
    if(code/=0) then
        diagnostics=[integer(8) :: 3,0,evals,1,code,light_calls,timing_calls,shared_calls,shared_steps];return
    endif
    converged=method==0 .or. nfree==0;lambda=controls(7);it=0;stagnant=0
    if(method==1 .and. nfree>0) then
        do it=1,maxiter
            if(controls(2)>0.d0 .and. evals>=int(controls(2))) exit
            call derivatives(q,mu,var,grad,info,ok)
            if(.not.ok) exit
            ! At an active bound, an outward descent direction is infeasible.
            ! Solve for the remaining coordinates rather than clipping a
            ! coupled unconstrained step repeatedly. Reconsider each iteration
            ! so a coordinate can leave its bound when the gradient reverses.
            blocked=(q<=ql+1.d-10 .and. grad>0.d0) .or. &
                    (q>=qu-1.d-10 .and. grad<0.d0)
            unit=-grad
            where(blocked) grad=0.d0
            if(maxval(abs(grad)/sqrt(max(1.d-30,abs([(info(i,i),i=1,nfree)]))))< &
                max(gtol,sqrt(epsilon(1.d0)*max(1.d0,abs(f))))) then
                converged=.true.;exit
            endif
            a=info
            do i=1,nfree
                a(i,i)=a(i,i)+lambda*max(abs(info(i,i)),1.d-12)
            enddo
            call bounded_step(a,unit,delta,ok)
            if(.not.ok) then
                lambda=lambda*10.d0
                if(lambda>1.d18) exit
                cycle
            endif
            candidate=max(ql,min(qu,q+delta));delta=candidate-q
            call objective(candidate,cf,cm,cv,code)
            if(code==0 .and. cf<f) then
                ratio=f-cf;q=candidate;f=cf;mu=cm;var=cv;lambda=max(lambda/3.d0,1.d-12)
                if(maxval(abs(delta))<xtol .and. ratio<ftol) then
                    converged=.true.;exit
                endif
            else
                lambda=lambda*10.d0
                if(lambda>1.d18) exit
            endif
        enddo
    else if(method==2 .and. nfree>0) then
        vertices=0.d0;values(1)=f
        do i=1,nfree
            vertices(:,i+1)=q
            step=min(controls(6),qu(i))
            if(step<1.d-8) step=max(-controls(6),ql(i))
            vertices(i,i+1)=step
            call objective(vertices(:,i+1),values(i+1),cm,cv,code)
        enddo
        do it=1,maxiter
            if(controls(2)>0.d0 .and. evals>=int(controls(2))) exit
            order=[(i,i=1,nfree+1)]
            do i=2,nfree+1
                tmp=order(i);j=i-1
                do while(j>=1)
                    if(values(order(j))<=values(tmp)) exit
                    order(j+1)=order(j);j=j-1
                enddo
                order(j+1)=tmp
            enddo
            vertices=vertices(:,order);values=values(order)
            spread=maxval(abs(values-values(1)))
            diameter=maxval(abs(vertices-spread_shape(vertices(:,1))))
            if(controls(4)==0.d0) xtol=sqrt(settings(7))
            if(spread<ftol .and. diameter<xtol) then
                converged=.true.;exit
            endif
            centre=sum(vertices(:,1:nfree),dim=2)/nfree
            xr=max(ql,min(qu,2.d0*centre-vertices(:,nfree+1)))
            call objective(xr,fr,cm,cv,code)
            if(fr<values(1)) then
                xe=max(ql,min(qu,centre+2.d0*(xr-centre)))
                call objective(xe,fe,cm,cv,code)
                if(fe<fr) then
                    vertices(:,nfree+1)=xe;values(nfree+1)=fe
                else
                    vertices(:,nfree+1)=xr;values(nfree+1)=fr
                endif
            else if(fr<values(nfree)) then
                vertices(:,nfree+1)=xr;values(nfree+1)=fr
            else
                if(fr<values(nfree+1)) then
                    xe=centre+0.5d0*(xr-centre);threshold=fr
                else
                    xe=centre+0.5d0*(vertices(:,nfree+1)-centre);threshold=values(nfree+1)
                endif
                xe=max(ql,min(qu,xe));call objective(xe,fe,cm,cv,code)
                if(fe<threshold) then
                    vertices(:,nfree+1)=xe;values(nfree+1)=fe
                else
                    do j=2,nfree+1
                        vertices(:,j)=0.5d0*(vertices(:,1)+vertices(:,j))
                        call objective(vertices(:,j),values(j),cm,cv,code)
                    enddo
                endif
            endif
        enddo
    endif
    optimizing=.false.
    q=best;call objective(q,f,mu,var,code)
    if(allocated(cached_events)) then
        final_event_count=cached_counts;final_event_model=cached_events
    endif
    j=0
    do i=1,nobs
        if(meta(i,1)>2) cycle
        j=j+1;final_rv(i)=cached_rv(j);final_flux(i)=cached_flux(j);final_parts(i,:)=cached_parts(j,:)
    enddo
    final_ma=cached_ma(1:npl)
    parameters(free_index)=origin+scale*q;models=mu;variances=var;loglik=-f
    components=0.d0
    do i=1,nobs
        j=meta(i,1)
        if(j>=3) then
            if(config(j+7)==0) cycle
        endif
        components(j)=components(j)-0.5d0*((data(i,2)-mu(i))**2/var(i)+log(var(i))+ &
            1.83787706640934548356d0)
    enddo
    ! Publish the data likelihood independently of the posterior objective,
    ! avoiding cancellation for very strong priors.
    if(any(priors(:,3)/=0.d0) .or. any(priors(:,6)/=0.d0)) loglik=sum(components)
    covcode=merge(0,1,converged)
    if(nfree>0 .and. method/=0 .and. converged) then
        call derivatives(q,mu,var,grad,info,ok)
        if(.not.ok) covcode=2
        if(method==2 .and. config(1)==0 .and. ok) then
            ! Observed negative-log-likelihood Hessian for Keplerian Simplex,
            ! Fisher/Gauss-Newton for interacting systems, as for legacy RV.
            do i=1,nfree
                h(i)=max(2.d-3,abs(q(i))*2.d-4)
                if(q(i)-h(i)<=ql(i) .or. q(i)+h(i)>=qu(i)) then
                    ok=.false.;covcode=3
                endif
            enddo
            if(ok) then
                do i=1,nfree
                    qp=q;qm=q;qp(i)=qp(i)+h(i);qm(i)=qm(i)-h(i)
                    call objective(qp,fp,mp,vp,code);if(code/=0) ok=.false.
                    call objective(qm,fm,mm,vm,code);if(code/=0) ok=.false.
                    info(i,i)=(fp-2.d0*f+fm)/h(i)**2
                    do j=1,i-1
                        qp=q;qp(i)=q(i)+h(i);qp(j)=q(j)+h(j)
                        call objective(qp,fp,mp,vp,code);if(code/=0) ok=.false.
                        qp(j)=q(j)-h(j);call objective(qp,fm,mp,vp,code);if(code/=0) ok=.false.
                        ratio=fp-fm;qp(i)=q(i)-h(i);qp(j)=q(j)+h(j)
                        call objective(qp,fp,mp,vp,code);if(code/=0) ok=.false.
                        qp(j)=q(j)-h(j);call objective(qp,fm,mp,vp,code);if(code/=0) ok=.false.
                        info(i,j)=(ratio-fp+fm)/(4.d0*h(i)*h(j));info(j,i)=info(i,j)
                    enddo
                enddo
            endif
        endif
        if(any(q<=ql+1.d-7) .or. any(q>=qu-1.d-7)) then
            ok=.false.;covcode=4
        endif
        if(ok) then
            do i=1,nfree
                unit=0.d0;unit(i)=1.d0;call solve(info,unit,z,ok)
                if(.not.ok) then
                    covcode=5;exit
                endif
                inv(:,i)=z
            enddo
        endif
        if(ok) then
            do i=1,nfree
                do j=1,nfree
                    covariance(i,j)=inv(i,j)*scale(i)*scale(j)
                enddo
            enddo
            errors=sqrt(max(0.d0,[(covariance(i,i),i=1,nfree)]))
        else
            if(covcode==0) covcode=2
        endif
    endif
    diagnostics=[integer(8) :: merge(0,1,converged),min(it,maxiter),evals,covcode,0, &
        light_calls,timing_calls,shared_calls,shared_steps]
contains
    function spread_shape(v) result(out)
        real(8),intent(in) :: v(nfree)
        real(8) :: out(nfree,nfree+1)
        integer :: l
        do l=1,nfree+1
            out(:,l)=v
        enddo
    end function
    subroutine solve(matrix,rhs,solution,success)
        real(8),intent(in) :: matrix(nfree,nfree),rhs(nfree)
        real(8),intent(out) :: solution(nfree)
        logical,intent(out) :: success
        real(8) :: l(nfree,nfree),y(nfree),v
        integer :: ii,jj
        l=0.d0;success=.false.
        do ii=1,nfree
            do jj=1,ii
                v=matrix(ii,jj)-dot_product(l(ii,1:jj-1),l(jj,1:jj-1))
                if(ii==jj) then
                    if(v<=1.d-14*max(1.d-30,abs(matrix(ii,ii)))) return
                    l(ii,jj)=sqrt(v)
                else
                    l(ii,jj)=v/l(jj,jj)
                endif
            enddo
        enddo
        do ii=1,nfree
            y(ii)=(rhs(ii)-dot_product(l(ii,1:ii-1),y(1:ii-1)))/l(ii,ii)
        enddo
        do ii=nfree,1,-1
            solution(ii)=(y(ii)-dot_product(l(ii+1:nfree,ii),solution(ii+1:nfree)))/l(ii,ii)
        enddo
        success=all(ieee_is_finite(solution))
    end subroutine
    subroutine bounded_step(matrix,rhs,solution,success)
        ! Minimize the damped quadratic inside the box. When a trial direction
        ! reaches a bound, retain that displacement and recompute the other
        ! coordinates, including its cross terms. Clipping a full unconstrained
        ! step loses those terms and can repeatedly reject correlated steps.
        real(8),intent(in) :: matrix(nfree,nfree),rhs(nfree)
        real(8),intent(out) :: solution(nfree)
        logical,intent(out) :: success
        real(8) :: reduced(nfree,nfree),residual(nfree),direction(nfree), &
            low(nfree),high(nfree),alpha,fraction,violation,worst
        integer :: active(nfree),ii,iteration,hit,release
        low=ql-q;high=qu-q;solution=0.d0;active=0;success=.false.
        where(blocked .and. q<=ql+1.d-10) active=-1
        where(blocked .and. q>=qu-1.d-10) active=1
        do iteration=1,10*nfree+10
            reduced=matrix;residual=rhs-matmul(matrix,solution)
            do ii=1,nfree
                if(active(ii)==0) cycle
                reduced(ii,:)=0.d0;reduced(:,ii)=0.d0
                reduced(ii,ii)=1.d0;residual(ii)=0.d0
            enddo
            call solve(reduced,residual,direction,success)
            if(.not.success) return
            alpha=1.d0;hit=0
            do ii=1,nfree
                if(active(ii)/=0 .or. direction(ii)==0.d0) cycle
                if(direction(ii)>0.d0) then
                    fraction=(high(ii)-solution(ii))/direction(ii)
                else
                    fraction=(low(ii)-solution(ii))/direction(ii)
                endif
                if(fraction>=alpha) cycle
                alpha=max(0.d0,fraction);hit=ii
            enddo
            solution=max(low,min(high,solution+alpha*direction))
            if(hit/=0) then
                active(hit)=merge(1,-1,direction(hit)>0.d0)
                cycle
            endif
            residual=rhs-matmul(matrix,solution);worst=1.d-10;release=0
            do ii=1,nfree
                if(active(ii)==0) cycle
                violation=-active(ii)*residual(ii)/sqrt(max(abs(matrix(ii,ii)),1.d-30))
                if(violation<=worst) cycle
                worst=violation;release=ii
            enddo
            if(release==0) return
            active(release)=0
        enddo
        success=.false.
    end subroutine
    subroutine derivatives(x,m,v,g,information,success)
        real(8),intent(in) :: x(nfree),m(nobs),v(nobs)
        real(8),intent(out) :: g(nfree),information(nfree,nfree)
        logical,intent(out) :: success
        real(8) :: xp(nfree),xm(nfree),mplus(nobs),mminus(nobs),vplus(nobs),vminus(nobs), &
            ff,hp,hm,res(nobs),mg(nobs,nfree),vg(nobs,nfree)
        integer :: ii,jj,sp,sm
        success=.false.
        do ii=1,nfree
            step=controls(8)*max(1.d-4,abs(x(ii))*1.d-5)
            xp=x;xm=x;xp(ii)=min(qu(ii),x(ii)+step);xm(ii)=max(ql(ii),x(ii)-step)
            hp=xp(ii)-x(ii);hm=x(ii)-xm(ii)
            if(hp+hm<=0.d0) return
            call objective(xp,ff,mplus,vplus,sp);call objective(xm,ff,mminus,vminus,sm)
            if(sp/=0 .or. sm/=0) return
            mg(:,ii)=(mplus-mminus)/(hp+hm);vg(:,ii)=(vplus-vminus)/(hp+hm)
        enddo
        res=m-data(:,2)
        do ii=1,nfree
            g(ii)=sum(weight*(res/v*mg(:,ii)+0.5d0*(1.d0/v-res**2/v**2)*vg(:,ii)))
            do jj=1,nfree
                information(ii,jj)=sum(weight*(mg(:,ii)*mg(:,jj)/v+0.5d0*vg(:,ii)*vg(:,jj)/v**2))
            enddo
        enddo
        do ii=1,nfree
            jj=free_index(ii)
            ff=origin(ii)+scale(ii)*x(ii)
            if(priors(jj,3)/=0.d0) then
                g(ii)=g(ii)+(ff-priors(jj,1))*scale(ii)/priors(jj,2)**2
                information(ii,ii)=information(ii,ii)+(scale(ii)/priors(jj,2))**2
            endif
            if(priors(jj,6)/=0.d0) then
                if(ff<=0.d0) return
                g(ii)=g(ii)+scale(ii)/ff
                information(ii,ii)=information(ii,ii)-(scale(ii)/ff)**2
            endif
        enddo
        success=all(ieee_is_finite(g)) .and. all(ieee_is_finite(information))
    end subroutine
    subroutine prior_penalty(p,penalty,istat)
        real(8),intent(in) :: p(npar)
        real(8),intent(out) :: penalty
        integer,intent(out) :: istat
        integer :: ii
        penalty=0.d0;istat=1
        do ii=1,npar
            if(priors(ii,3)/=0.d0) then
                penalty=penalty+0.5d0*((p(ii)-priors(ii,1))/priors(ii,2))**2+ &
                    log(priors(ii,2))+0.91893853320467274178d0
            endif
            if(priors(ii,6)/=0.d0) then
                if(p(ii)<=0.d0) return
                penalty=penalty+log(p(ii))+log(log(priors(ii,5))-log(priors(ii,4)))
            endif
        enddo
        if(ieee_is_finite(penalty)) istat=0
    end subroutine
    subroutine objective(x,value,model,variance,istat)
        real(8),intent(in) :: x(nfree)
        real(8),intent(out) :: value,model(nobs),variance(nobs)
        integer,intent(out) :: istat
        real(8) :: p(npar),v(npack),planets(npl,9),raw(npl,10),sorted(npl,9),masses(9), &
            sm,rs,rho,pi,inner,kbar,coeff,ml,mh,mid,axis(npl),tzero(9),tset(60,6),rset(60,2), &
            base(60,3),trend(2),rv(nobs),flux(nobs),res(nobs),plog(nobs),byset(60,2), &
            comp(2),ll,chi(2),radii(9),parts(nobs,npl),ma(9),te(9),el(npl,6), &
            ev(max(1,config(8)),npl,4),tm(nobs),tr(nobs),tl(nobs),tb(10), &
            rsky,vsky,b,rr,rstar,depth,duration,ov1,ov2,triangle,rvtime(nobs),obs(nobs),err(nobs),cov(nobs), &
            sorted_axis(npl),dummy(nobs),penalty
        real(8) :: query_rv(nsample),query_flux(nsample),query_parts(nsample,npl)
        integer :: query_kind(nsample),query_set(nsample),sample_count(nobs)
        integer :: ip,jp,kk,idx,offset,nr,nt,inds(nobs),td(nobs),tp(nobs),tn(nobs),cnt(npl), &
            permutation(npl),inverse(npl),tt,sc,n
        pi=acos(-1.d0);istat=1;value=huge(1.d0)/100.d0;model=0.d0;variance=1.d0
        if(optimizing .and. controls(2)>0.d0 .and. evals>=int(controls(2))) return
        evals=evals+1
        p=initial;p(free_index)=origin+scale*x;v=packed
        call prior_penalty(p,penalty,istat)
        if(istat/=0) return
        istat=1
        do ip=1,npack
            if(map(ip)>0) v(ip)=p(map(ip))
        enddo
        if(.not.all(ieee_is_finite(v))) return
        do ip=1,npl
            raw(ip,:)=v(10*(ip-1)+1:10*ip)
        enddo
        offset=10*npl;sm=v(offset+1);rs=v(offset+2);trend=v(offset+3:offset+4)
        if(sm<=0.d0 .or. rs<=0.d0) return
        do ip=1,60
            rset(ip,:)=v(offset+4+2*(ip-1)+1:offset+4+2*ip)
            tset(ip,:)=v(offset+124+9*(ip-1)+1:offset+124+9*(ip-1)+6)
            base(ip,:)=v(offset+124+9*(ip-1)+7:offset+124+9*ip)
        enddo
        if(config(5)==1) rset(:,2)=rset(config(9),2)
        if(any(rset(:,2)<0.d0) .or. any(tset(:,2)<0.d0) .or. any(base(:,3)<0.d0)) return
        if(npl==0) then
            ! Flat physical model: retain the same offsets, trends and variances.
            physical_key=v(1:2);light_key=tset(:,3:6);light_valid=.true.
            goto 710
        endif
        planets(:,1:7)=raw(:,1:7);planets(:,8:9)=raw(:,9:10)
        if(config(3)==1) then
            planets(:,3)=sqrt(raw(:,3)**2+raw(:,4)**2)
            planets(:,4)=atan2(raw(:,3),raw(:,4))*180.d0/pi
            planets(:,5)=raw(:,5)-planets(:,4)
        endif
        if(config(4)==1 .and. config(1)==1) planets(:,6)=planets(1,6)
        if(any(planets(:,1)<0.d0) .or. any(planets(:,2)<=0.d0)) return
        if(any(planets(:,3)<0.d0) .or. any(planets(:,3)>=1.d0)) return
        if(any(planets(:,6)<=0.d0) .or. any(planets(:,6)>=180.d0)) return
        if(config(6)==1) then
            rho=raw(1,9)
            call rv_transit_density(planets(:,2),rho,sm,axis,rs,sc,npl)
            if(sc/=0) return
            planets(:,8)=axis
        endif
        permutation=[(ip,ip=1,npl)]
        do ip=2,npl
            kk=permutation(ip);jp=ip-1
            do while(jp>=1)
                if(planets(permutation(jp),2)<=planets(kk,2)) exit
                permutation(jp+1)=permutation(jp);jp=jp-1
            enddo
            permutation(jp+1)=kk
        enddo
        masses=0.d0;inner=sm
        do ip=1,npl
            kk=permutation(ip);inverse(kk)=ip
            kbar=planets(kk,1)*sqrt(1.d0-planets(kk,3)**2)/abs(sin(planets(kk,6)*pi/180.d0))
            coeff=(2.d0*pi*1.32712440018d20/(planets(kk,2)*86400.d0))**(1.d0/3.d0)
            ml=0.d0;mh=max(2.d0*kbar*inner**(2.d0/3.d0)/coeff,1.d-12)
            do tt=1,100
                if(coeff*mh/(inner+mh)**(2.d0/3.d0)>=kbar) exit
                mh=mh*2.d0
            enddo
            do tt=1,80
                mid=(ml+mh)*0.5d0
                if(coeff*mid/(inner+mid)**(2.d0/3.d0)>kbar) then
                    mh=mid
                else
                    ml=mid
                endif
            enddo
            masses(ip)=(ml+mh)*0.5d0
            if(kbar==0.d0) masses(ip)=0.d0
            inner=inner+masses(ip)
        enddo
        sorted=planets(permutation,:);tzero=0.d0;tzero(1:npl)=raw(permutation,8)
        if(.not.any(meta(:,1)==2) .and. config(6)==0) then
            inner=sm
            do ip=1,npl
                sorted(ip,8)=(1.32712440018d20*(inner+masses(ip))* &
                    (sorted(ip,2)*86400.d0/(2.d0*pi))**2)**(1.d0/3.d0)/(rs*6.957d8)
                inner=inner+masses(ip)
            enddo
        endif
710     nr=count(meta(:,1)<=2);nt=nobs-nr
        if(config(1)==1 .and. config(13)==1 .and. nt>0) then
            if(.not.light_valid .or. .not.timing_valid .or. any(physical_key/=v(1:10*npl+2)) &
                .or. any(light_key/=tset(:,3:6))) then
                query_kind=1;query_set=1
                if(nr>0) then
                    inds(1:nr)=pack([(ip,ip=1,nobs)],meta(:,1)<=2)
                    do ip=1,nsample
                        query_kind(ip)=meta(inds(sample_owner(ip)),1)
                        query_set(ip)=meta(inds(sample_owner(ip)),2)
                    enddo
                endif
                shared_calls=shared_calls+1
                call joint_shared_nbody(sample_time,query_kind,query_set,sorted(:,2:7),masses(1:npl), &
                    sorted(:,9),sm,rs,settings(1),settings(4),settings(5),ld_laws,tset,settings(2), &
                    min(settings(3),1.d-14),1000000,config(2),config(8),query_rv,query_flux, &
                    cached_counts,cached_events,trial_steps,query_parts,sc,nsample,npl)
                shared_steps=shared_steps+trial_steps
                if(sc/=0) then
                    istat=sc;return
                endif
                if(nr>0) then
                    cached_rv(1:nr)=0.d0;cached_flux(1:nr)=0.d0;cached_parts(1:nr,:)=0.d0;sample_count=0
                    do ip=1,nsample
                        kk=sample_owner(ip);sample_count(kk)=sample_count(kk)+1
                        cached_rv(kk)=cached_rv(kk)+query_rv(ip)
                        cached_flux(kk)=cached_flux(kk)+query_flux(ip)
                        cached_parts(kk,:)=cached_parts(kk,:)+query_parts(ip,:)
                    enddo
                    if(any(sample_count(1:nr)==0)) return
                    cached_rv(1:nr)=cached_rv(1:nr)/sample_count(1:nr)
                    cached_flux(1:nr)=cached_flux(1:nr)/sample_count(1:nr)
                    do ip=1,nr
                        cached_parts(ip,:)=cached_parts(ip,:)/sample_count(ip)
                    enddo
                endif
                cached_ma=0.d0;cached_ma(1:npl)=modulo(sorted(:,5),360.d0)
                physical_key=v(1:10*npl+2);timing_key=physical_key;light_key=tset(:,3:6)
                light_valid=.true.;timing_valid=.true.
            endif
        endif
        if(nr>0) then
            inds(1:nr)=pack([(ip,ip=1,nobs)],meta(:,1)<=2)
            do ip=1,nr
                kk=inds(ip);rvtime(ip)=data(kk,1);obs(ip)=data(kk,2);err(ip)=data(kk,3);cov(ip)=data(kk,4)
            enddo
            sc=0
            if(.not.light_valid .or. any(physical_key/=v(1:10*npl+2)) .or. any(light_key/=tset(:,3:6))) then
            light_calls=light_calls+1
            call rv_transit_likelihood(sample_time,sample_owner,rvtime(1:nr),meta(inds(1:nr),2), &
                meta(inds(1:nr),1),obs(1:nr),err(1:nr),cov(1:nr),sorted,masses,sm,rs,settings(1), &
                config(1),config(2),config(7),tzero,ld_laws,tset,rset,base,trend,64,settings(2), &
                settings(3),1000000,rv(1:nr),flux(1:nr),tm(1:nr),res(1:nr),tr(1:nr), &
                plog(1:nr),byset,comp,ll,chi,radii,parts(1:nr,:),ma,te,sc,nsample,nr,npl)
            if(sc/=0) then
                istat=sc;return
            endif
            cached_rv(1:nr)=rv(1:nr);cached_flux(1:nr)=flux(1:nr);cached_ma=ma;cached_parts(1:nr,:)=parts(1:nr,:)
            physical_key=v(1:10*npl+2);light_key=tset(:,3:6);light_valid=.true.
            endif
            do ip=1,nr
                kk=inds(ip);idx=meta(kk,2)
                if(meta(kk,1)==1) then
                    model(kk)=cached_rv(ip)+trend(1)*data(kk,1)+trend(2)*data(kk,1)**2+rset(idx,1)
                    variance(kk)=data(kk,3)**2+rset(idx,2)**2
                else
                    model(kk)=(cached_flux(ip)+base(idx,1)*data(kk,4)+ &
                        base(idx,2)*data(kk,4)**2)*base(idx,3)+1.d0-base(idx,3)-tset(idx,1)
                    variance(kk)=data(kk,3)**2+tset(idx,2)**2
                endif
            enddo
            if(config(7)==1) sorted(:,5)=cached_ma(1:npl)
        endif
        if(nt>0) then
            inds(1:nt)=pack([(ip,ip=1,nobs)],meta(:,1)>2)
            td(1:nt)=meta(inds(1:nt),2);tp(1:nt)=inverse(meta(inds(1:nt),3));tn(1:nt)=meta(inds(1:nt),4)
            el=sorted(:,2:7);sorted_axis=sorted(:,8)*rs
            dummy=1.d0
            if(.not.timing_valid .or. any(timing_key/=v(1:10*npl+2))) then
            timing_calls=timing_calls+1
            call joint_timing_likelihood(el,masses(1:npl),sm,settings(1),settings(4),settings(5), &
                config(1),sorted_axis,settings(2),min(settings(3),1.d-14),1000000,config(2),512,1000000, &
                config(8),td(1:nt),tp(1:nt),tn(1:nt),data(inds(1:nt),1),dummy(1:nt),cnt,ev,tm(1:nt), &
                tr(1:nt),tl(1:nt),tb,ll,sc,nt,npl)
            if(sc/=0) then
                istat=sc;return
            endif
            cached_events=ev;cached_counts=cnt;timing_key=v(1:10*npl+2);timing_valid=.true.
            else
                ev=cached_events;cnt=cached_counts
            endif
            do ip=1,nt
                kk=inds(ip);jp=tp(ip);n=tn(ip);idx=meta(kk,2)
                if(n>cnt(jp) .or. n>config(8)) then
                    istat=8;return
                endif
                select case(meta(kk,1))
                case(3)
                    model(kk)=ev(n,jp,2)
                case(4,5)
                    rho=v(offset+664+2*(idx-1)+1);rr=v(offset+664+2*idx)
                    ! With LC data, use the physical photodynamical radius and Rp/R.
                    if(nr>0 .and. any(meta(:,1)==2)) then
                        rstar=rs*6.957d8/1.49597892d11;rr=sorted(jp,9)
                    else
                        if(rho<=0.d0) return
                        rstar=(3.d0*sm*1.9884098709677423d30/(4.d0*pi*rho))**(1.d0/3.d0)/1.49597892d11
                    endif
                    if(rr<0.d0) return
                    rsky=ev(n,jp,3);vsky=ev(n,jp,4);b=rsky/rstar;duration=0.d0;depth=0.d0
                    if(b<1.d0+rr .and. vsky>0.d0) then
                        duration=2.d0*rstar*sqrt((1.d0+rr)**2-b*b)/vsky
                        if(b<=abs(1.d0-rr)) then
                            depth=min(1.d0,rr*rr)
                        else
                            ov1=acos(max(-1.d0,min(1.d0,(b*b+1.d0-rr*rr)/(2.d0*b))))
                            ov2=acos(max(-1.d0,min(1.d0,(b*b+rr*rr-1.d0)/(2.d0*b*rr))))
                            triangle=sqrt(max(0.d0,(-b+1.d0+rr)*(b+1.d0-rr)*(b-1.d0+rr)*(b+1.d0+rr)))
                            depth=(ov1+rr*rr*ov2-0.5d0*triangle)/pi
                        endif
                    endif
                    if(meta(kk,1)==4) model(kk)=duration
                    if(meta(kk,1)==5) model(kk)=1.d0-depth
                end select
                variance(kk)=data(kk,3)**2
            enddo
        endif
        if(any(variance<=0.d0) .or. .not.all(ieee_is_finite(model))) return
        value=0.5d0*sum(weight*((data(:,2)-model)**2/variance+log(variance)+1.83787706640934548356d0))
        value=value+penalty
        if(.not.ieee_is_finite(value)) then
            value=huge(1.d0)/100.d0;return
        endif
        istat=0
        if(value<bf) then
            bf=value;best=x
        endif
    end subroutine
end subroutine native_joint_fit

! Sampler entry: keep model arrays inside Fortran and return only lnL/status.
subroutine native_joint_loglik(initial,packed,map,meta,data,sample_time,sample_owner, &
        ld_laws,config,settings,event_capacity,loglik,status,npar,npack,nobs,nsample,npl)
    implicit none
    interface
subroutine native_joint_fit(initial,free_index,lower,upper,scale,packed,map, &
        meta,data,sample_time,sample_owner,ld_laws,config,settings,method,maxiter,event_capacity, &
        parameters,models,variances,loglik,components,covariance,errors,diagnostics, &
        final_event_count,final_event_model,final_rv,final_flux,final_parts,final_ma,npar,nfree,npack,nobs,nsample,npl, &
        optimizer_options,prior_options)
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer,intent(in) :: npar,nfree,npack,nobs,nsample,npl,method,maxiter,event_capacity
    integer,intent(in) :: free_index(nfree),map(npack),meta(nobs,4),sample_owner(nsample),ld_laws(60),config(13)
    real(8),optional,intent(in) :: optimizer_options(9),prior_options(npar,6)
    real(8) :: controls(9),ftol,xtol,gtol
    logical :: optimizing
    real(8),intent(in) :: initial(npar),lower(nfree),upper(nfree),scale(nfree),packed(npack)
    real(8),intent(in) :: data(nobs,4),sample_time(nsample),settings(8)
    real(8),intent(out) :: parameters(npar),models(nobs),variances(nobs),loglik,components(5)
    real(8),intent(out) :: covariance(nfree,nfree),errors(nfree)
    integer(8),intent(out) :: diagnostics(9)
    integer,intent(out) :: final_event_count(npl)
    real(8),intent(out) :: final_event_model(event_capacity,npl,4)
    real(8),intent(out) :: final_rv(nobs),final_flux(nobs),final_parts(nobs,npl),final_ma(npl)
end subroutine native_joint_fit
    end interface
    integer,intent(in) :: npar,npack,nobs,nsample,npl,event_capacity
    integer,intent(in) :: map(npack),meta(nobs,4),sample_owner(nsample),ld_laws(60),config(13)
    real(8),intent(in) :: initial(npar),packed(npack),data(nobs,4),sample_time(nsample),settings(8)
    real(8),intent(out) :: loglik
    integer,intent(out) :: status
    real(8) :: parameters(npar),models(nobs),variances(nobs),components(5),covariance(1,1),errors(1)
    real(8) :: events(event_capacity,npl,4),rv(nobs),flux(nobs),parts(nobs,npl),ma(npl)
    integer :: counts(npl)
    integer(8) :: diagnostics(9)
    call native_joint_fit(initial,[npar],[initial(npar)],[initial(npar)],[1.d0],packed,map, &
        meta,data,sample_time,sample_owner,ld_laws,config,settings,0,0,event_capacity, &
        parameters,models,variances,loglik,components,covariance,errors,diagnostics, &
        counts,events,rv,flux,parts,ma,npar,1,npack,nobs,nsample,npl)
    status=int(diagnostics(1))
end subroutine native_joint_loglik
