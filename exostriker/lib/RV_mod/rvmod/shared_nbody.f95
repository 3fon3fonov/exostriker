! One outward trajectory per side of the element epoch. Observation samples and
! foreground conjunctions use degree-seven Hermite dense output from positions,
! velocities, Newtonian accelerations and jerks at accepted step endpoints.
! Unlike a separate dense transit-search grid, observations do not change the
! integration mesh. This routine never interpolates between osculating elements.
subroutine joint_shared_nbody(times,kind,dataset,elements,planet_mass,radius_ratio,stellar_mass,stellar_radius, &
        epoch,start_time,end_time,ld_laws,set_parameters,step_days,tolerance,max_steps,jacobi_mode,max_events, &
        rv,flux,event_count,event_model,steps,rv_parts,status,nobs,npl)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: nobs,npl,max_steps,jacobi_mode,max_events,kind(nobs),dataset(nobs),ld_laws(60)
    real(8),intent(in) :: times(nobs),elements(npl,6),planet_mass(npl),radius_ratio(npl)
    real(8),intent(in) :: stellar_mass,stellar_radius,epoch,start_time,end_time,step_days,tolerance,set_parameters(60,6)
    real(8),intent(out) :: rv(nobs),flux(nobs),event_model(max_events,npl,4),rv_parts(nobs,npl)
    integer,intent(out) :: event_count(npl),steps,status
    real(8),parameter :: pi=3.1415926535897932384626433832795d0
    real(8),parameter :: day=86400.d0,gsun=1.32712440018d20,rsun=6.957d8,au=1.49597892d11
    real(8) :: gm(10),inner(9),semiaxis(9),initial(10,6),state(10,6),jacobi(6,10),ybs(6,4106),seed(6,10)
    real(8) :: xb(10),yb(10),zb(10),vx(10),vy(10),vz(10),msys,msum,rstar
    real(8) :: gx(64),gw(64),norms(60),ld_coeff(4),distances(9),sizes(9),angles(9)
    real(8) :: acc0(3,10),jerk0(3,10),acc1(3,10),jerk1(3,10),old(6,10),polynomial(3,10,0:7)
    real(8) :: a(3,10),b(3,10),c(3,10),deltaj(3,10),eval(6,10),r(3),v(3)
    real(8) :: clock,previous,h,target,step,leftmost,rightmost,lo,hi,u,lft,rgt,middle,g0,g1,g,root,tmp(4)
    integer :: i,j,k,nbod,branch,direction,next_sample,first_positive,index,iterations,n,d,m
    logical :: valid,need_dense,transit_cross(npl),search_events
    status=1;rv=0.d0;flux=1.d0;rv_parts=0.d0;event_count=0;event_model=0.d0;steps=0
    if(npl<1 .or. npl>9 .or. nobs<1 .or. max_steps<1 .or. max_events<1) return
    if(.not.all(ieee_is_finite(times)) .or. .not.all(ieee_is_finite(elements))) return
    if(.not.all(ieee_is_finite(planet_mass)) .or. .not.all(ieee_is_finite(radius_ratio))) return
    if(.not.all(ieee_is_finite([stellar_mass,stellar_radius,epoch,start_time,end_time,step_days,tolerance]))) return
    if(.not.all(ieee_is_finite(set_parameters))) return
    if(stellar_mass<=0.d0 .or. stellar_radius<=0.d0 .or. step_days<=0.d0 .or. end_time<start_time) return
    if(tolerance<1.d-15 .or. tolerance>1.d-3 .or. jacobi_mode<0 .or. jacobi_mode>1) return
    if(any(planet_mass<0.d0) .or. any(radius_ratio<0.d0)) return
    if(any(elements(:,1)<=0.d0) .or. any(elements(:,2)<0.d0) .or. any(elements(:,2)>=1.d0)) return
    if(any(elements(:,5)<0.d0) .or. any(elements(:,5)>180.d0)) return
    if(any(kind<1) .or. any(kind>2) .or. any(dataset<1) .or. any(dataset>60)) return
    do i=2,nobs
        if(times(i)<times(i-1)) return
    enddo
    do j=2,npl
        if(elements(j,1)<elements(j-1,1)) return
    enddo
    call rv_transit_quadrature(64,gx,gw)
    do d=1,60
        call rv_transit_profile(ld_laws(d),set_parameters(d,3:6),norms(d),status)
        if(status/=0) return
    enddo
    rstar=stellar_radius*rsun;gm=0.d0;gm(1)=stellar_mass*gsun;msum=gm(1);nbod=npl+1
    step=step_days*day
    do j=1,npl
        gm(j+1)=planet_mass(j)*gsun;msum=msum+gm(j+1);inner(j)=msum
        if(jacobi_mode==1) inner(j)=gm(1)*msum/(msum-gm(j+1))
        semiaxis(j)=(inner(j)*(elements(j,1)*day/(2.d0*pi))**2)**(1.d0/3.d0)
        if(any(kind==2) .and. semiaxis(j)*(1.d0-elements(j,2))<=rstar) then
            status=1;return
        endif
        step=min(step,elements(j,1)*day*(1.d0-elements(j,2))**1.5d0 &
            /sqrt(1.d0+elements(j,2))/32.d0)
    enddo
    call kepler_state(0.d0,jacobi)
    if(status==2) return
    initial=transpose(jacobi)
    call coord_j2h(nbod,gm,initial(:,1),initial(:,2),initial(:,3),initial(:,4),initial(:,5),initial(:,6), &
        state(:,1),state(:,2),state(:,3),state(:,4),state(:,5),state(:,6))
    call coord_h2b(nbod,gm,state(:,1),state(:,2),state(:,3),state(:,4),state(:,5),state(:,6), &
        xb,yb,zb,vx,vy,vz,msys)
    seed=0.d0
    do j=1,nbod
        seed(:,j)=[xb(j),yb(j),zb(j),vx(j),vy(j),vz(j)]
    enddo
    lo=(start_time-epoch)*day;hi=(end_time-epoch)*day
    ! Equal search bounds request observation models only.
    search_events=end_time>start_time
    leftmost=min(0.d0,times(1)*day);rightmost=max(0.d0,times(nobs)*day)
    if(search_events) then
        ! Pad the bounds, so an endpoint conjunction has a genuine bracket.
        leftmost=min(leftmost,lo-step);rightmost=max(rightmost,hi+step)
    endif
    first_positive=nobs+1
    do i=1,nobs
        if(times(i)>=0.d0) then
            first_positive=i;exit
        endif
    enddo
    do i=1,nobs
        if(times(i)==0.d0) call observation(i,seed)
    enddo
    do branch=1,2
        direction=2*branch-3
        target=leftmost
        if(direction==1) target=rightmost
        if(target==0.d0) cycle
        ybs=0.d0;ybs(:,1:10)=seed;clock=0.d0
        next_sample=first_positive-1
        if(direction==1) then
            next_sample=first_positive
            do while(next_sample<=nobs)
                if(times(next_sample)>0.d0) exit
                next_sample=next_sample+1
            enddo
        endif
        do while(abs(target-clock)>max(1.d-8,abs(target)*2.d-15))
            steps=steps+1
            if(steps>max_steps) then
                status=4;return
            endif
            old=ybs(:,1:10);previous=clock
            h=sign(min(step,abs(target-clock)),target-clock)
            call bs_int_pl(nbod,0,gm,0.d0,0.d0,clock,h,ybs,tolerance)
            if(.not.all(ieee_is_finite(ybs(:,1:nbod))) .or. .not.ieee_is_finite(clock)) then
                status=5;return
            endif
            if((clock-previous)*direction<=0.d0) then
                status=5;return
            endif
            h=clock-previous
            need_dense=.false.;transit_cross=.false.
            if(next_sample>=1 .and. next_sample<=nobs) then
                need_dense=(times(next_sample)*day-clock)*direction<=max(1.d-8,abs(clock)*2.d-15)
            endif
            if(search_events .and. max(previous,clock)>=lo .and. min(previous,clock)<=hi) then
                do j=1,npl
                    if(old(3,j+1)<=old(3,1) .and. ybs(3,j+1)<=ybs(3,1)) cycle
                    r=old(1:3,j+1)-old(1:3,1);v=old(4:6,j+1)-old(4:6,1)
                    g0=dot_product(r(1:2),v(1:2))
                    r=ybs(1:3,j+1)-ybs(1:3,1);v=ybs(4:6,j+1)-ybs(4:6,1)
                    g1=dot_product(r(1:2),v(1:2))
                    if(direction==1) then
                        transit_cross(j)=g0<=0.d0 .and. g1>0.d0
                    else
                        transit_cross(j)=g0>0.d0 .and. g1<=0.d0
                    endif
                enddo
            endif
            if(.not.need_dense .and. .not.any(transit_cross)) cycle
            ! Build dense output only for steps containing observations/events.
            call acceleration(old,acc0,jerk0,valid)
            if(.not.valid) then
                status=5;return
            endif
            call acceleration(ybs(:,1:10),acc1,jerk1,valid)
            if(.not.valid) then
                status=5;return
            endif
            polynomial(:,:,0)=old(1:3,:);polynomial(:,:,1)=h*old(4:6,:)
            polynomial(:,:,2)=0.5d0*h*h*acc0;polynomial(:,:,3)=h**3*jerk0/6.d0
            a=ybs(1:3,1:10)-sum(polynomial(:,:,0:3),dim=3)
            b=h*ybs(4:6,1:10)-polynomial(:,:,1)-2*polynomial(:,:,2)-3*polynomial(:,:,3)
            c=h*h*acc1-2*polynomial(:,:,2)-6*polynomial(:,:,3)
            deltaj=h**3*jerk1-6*polynomial(:,:,3)
            polynomial(:,:,4)=35*a-15*b+2.5d0*c-deltaj/6.d0
            polynomial(:,:,5)=-84*a+39*b-7*c+deltaj/2.d0
            polynomial(:,:,6)=70*a-34*b+6.5d0*c-deltaj/2.d0
            polynomial(:,:,7)=-20*a+10*b-2*c+deltaj/6.d0
            do while(next_sample>=1 .and. next_sample<=nobs)
                if((times(next_sample)*day-clock)*direction>max(1.d-8,abs(clock)*2.d-15)) exit
                u=max(0.d0,min(1.d0,(times(next_sample)*day-previous)/h))
                call interpolate(u,eval)
                call observation(next_sample,eval)
                next_sample=next_sample+direction
            enddo
            do j=1,npl
                if(.not.transit_cross(j)) cycle
                lft=0.d0;rgt=1.d0
                do iterations=1,64
                    middle=(lft+rgt)*0.5d0;call interpolate(middle,eval)
                    r=eval(1:3,j+1)-eval(1:3,1);v=eval(4:6,j+1)-eval(4:6,1)
                    g=dot_product(r(1:2),v(1:2))
                    if(g*direction<=0.d0) then
                        lft=middle
                    else
                        rgt=middle
                    endif
                    if((rgt-lft)*abs(h)<day*1.d-12) exit
                enddo
                u=0.5d0*(lft+rgt);root=previous+u*h
                call interpolate(u,eval)
                r=eval(1:3,j+1)-eval(1:3,1);v=eval(4:6,j+1)-eval(4:6,1)
                if(r(3)<=0.d0 .or. root<lo-day*1.d-10 .or. root>hi+day*1.d-10) cycle
                root=max(lo,min(hi,root));n=event_count(j)+1
                if(n>max_events) then
                    status=10;return
                endif
                event_count(j)=n
                event_model(n,j,:)=[0.d0,epoch+root/day,sqrt(dot_product(r(1:2),r(1:2)))/au, &
                    sqrt(dot_product(v(1:2),v(1:2)))*day/au]
            enddo
        enddo
    enddo
    ! Backward events were encountered in reverse time order. Restore the
    ! chronological, one-based numbering independently for every planet.
    do j=1,npl
        do i=2,event_count(j)
            tmp=event_model(i,j,:);k=i-1
            do while(k>=1)
                if(event_model(k,j,2)<=tmp(2)) exit
                event_model(k+1,j,:)=event_model(k,j,:);k=k-1
            enddo
            event_model(k+1,j,:)=tmp
        enddo
        do i=1,event_count(j)
            event_model(i,j,1)=dble(i)
        enddo
    enddo
    status=0
contains
    subroutine acceleration(s,acc,jerk,ok)
        real(8),intent(in) :: s(6,10)
        real(8),intent(out) :: acc(3,10),jerk(3,10)
        logical,intent(out) :: ok
        real(8) :: dx(3),dv(3),r2,r3,radial,force(3),deriv(3)
        integer :: ii,jj
        acc=0.d0;jerk=0.d0;ok=.false.
        do ii=1,nbod-1
            do jj=ii+1,nbod
                dx=s(1:3,jj)-s(1:3,ii);dv=s(4:6,jj)-s(4:6,ii);r2=dot_product(dx,dx)
                if(r2<=0.d0) return
                r3=1.d0/(r2*sqrt(r2));radial=3.d0*dot_product(dx,dv)/r2
                force=dx*r3;deriv=(dv-radial*dx)*r3
                acc(:,ii)=acc(:,ii)+gm(jj)*force;acc(:,jj)=acc(:,jj)-gm(ii)*force
                jerk(:,ii)=jerk(:,ii)+gm(jj)*deriv;jerk(:,jj)=jerk(:,jj)-gm(ii)*deriv
            enddo
        enddo
        ok=all(ieee_is_finite(acc)) .and. all(ieee_is_finite(jerk))
    end subroutine
    subroutine interpolate(f,s)
        real(8),intent(in) :: f
        real(8),intent(out) :: s(6,10)
        integer :: power
        s(1:3,:)=polynomial(:,:,7);s(4:6,:)=7*polynomial(:,:,7)
        do power=6,0,-1
            s(1:3,:)=s(1:3,:)*f+polynomial(:,:,power)
            if(power>=1) s(4:6,:)=s(4:6,:)*f+power*polynomial(:,:,power)
        enddo
        s(4:6,:)=s(4:6,:)/h
    end subroutine
    subroutine observation(ix,s)
        integer,intent(in) :: ix
        real(8),intent(in) :: s(6,10)
        real(8) :: xx,yy,separation
        integer :: jp,nn,ds
        rv(ix)=-s(6,1)
        do jp=1,npl
            rv_parts(ix,jp)=gm(jp+1)*(s(6,jp+1)-s(6,1))/msum
        enddo
        ds=dataset(ix);nn=0
        do jp=1,npl
            if(s(3,jp+1)<=s(3,1) .or. radius_ratio(jp)==0.d0) cycle
            xx=(s(1,jp+1)-s(1,1))/rstar;yy=(s(2,jp+1)-s(2,1))/rstar
            separation=sqrt(xx*xx+yy*yy)
            if(separation>=1.d0+radius_ratio(jp)) cycle
            nn=nn+1;distances(nn)=separation;sizes(nn)=radius_ratio(jp);angles(nn)=atan2(yy,xx)
        enddo
        ld_coeff=set_parameters(ds,3:6)
        call rv_transit_disk(nn,distances(1:nn),sizes(1:nn),angles(1:nn), &
            ld_laws(ds),ld_coeff,64,gx,gw,norms(ds),flux(ix))
    end subroutine
    subroutine kepler_state(t,s)
        real(8), intent(in) :: t
        real(8), intent(out) :: s(6,10)
        real(8) :: eccentricity, inc, node, omega, anomaly, axis, gm_inner
        real(8) :: E,lo,hi,f,trial,den,x,y,vx,vy,sq,cp,sp,co,so,ci,si
        real(8) :: d1(3),d2(3)
        integer :: p,k
        s=0.d0
        do p=1,npl
            eccentricity=elements(p,2); omega=elements(p,3)*pi/180.d0
            anomaly=modulo(elements(p,4)*pi/180.d0+2.d0*pi*t/(elements(p,1)*day),2.d0*pi)
            inc=elements(p,5)*pi/180.d0; node=elements(p,6)*pi/180.d0
            axis=semiaxis(p); gm_inner=inner(p)
            anomaly=modulo(anomaly+pi,2.d0*pi)-pi
            lo=-pi; hi=pi; E=anomaly
            do k=1,100
                f=E-eccentricity*sin(E)-anomaly
                if(abs(f)<2.d-15) exit
                if(f<0.d0) then
                    lo=E
                else
                    hi=E
                endif
                trial=E-f/(1.d0-eccentricity*cos(E))
                if(trial<=lo .or. trial>=hi) trial=0.5d0*(lo+hi)
                E=trial
            enddo
            if(k>100) then
                status=2; return
            endif
            ! Direct cos is essential at conjunction: recovering cos from
            ! sqrt(1-sin**2), as in the legacy helper, loses timing precision.
            cp=cos(omega); sp=sin(omega); co=cos(node); so=sin(node)
            ci=cos(inc); si=sin(inc)
            d1=[cp*co-sp*so*ci,cp*so+sp*co*ci,sp*si]
            d2=[-sp*co-cp*so*ci,-sp*so+cp*co*ci,cp*si]
            sq=sqrt(1.d0-eccentricity**2); den=1.d0-eccentricity*cos(E)
            x=axis*(cos(E)-eccentricity); y=axis*sq*sin(E)
            vx=-sqrt(gm_inner/axis)*sin(E)/den
            vy=sqrt(gm_inner/axis)*sq*cos(E)/den
            s(1:3,p+1)=d1*x+d2*y
            s(4:6,p+1)=d1*vx+d2*vy
        enddo
    end subroutine

end subroutine joint_shared_nbody
