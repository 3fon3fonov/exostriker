! Astrometric predictions in ExTRA's east/north convention. No observations,
! uncertainties, likelihoods or optimizer state are used by this entry point.
subroutine astrometry_model(times,orbits,standard,observer,epoch,scan_angle, &
        east,north,along_scan,parts,status,nobs,npl)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: nobs,npl
    ! Orbits: P [day], e, omega/i/node [rad], periastron JD, stellar axis [mas].
    ! Standard: RA/Dec [deg], parallax [mas], mu_alpha_star/mu_delta [mas/yr].
    real(8),intent(in) :: times(nobs),orbits(npl,7),standard(5),observer(nobs,3)
    real(8),intent(in) :: epoch,scan_angle(nobs)
    real(8),intent(out) :: east(nobs),north(nobs),along_scan(nobs),parts(nobs,npl,2)
    integer,intent(out) :: status
    real(8),parameter :: pi=3.1415926535897932384626433832795d0
    real(8) :: ra,dec,p,ecc,w,inc,node,axis,ma,ea,lo,hi,f,trial,x,y
    real(8) :: aa,bb,ff,gg,pa,pd
    integer :: i,j,k
    status=1
    east=0.d0
    north=0.d0
    along_scan=0.d0
    parts=0.d0
    if(nobs<1 .or. npl<0 .or. npl>9) return
    if(.not.all(ieee_is_finite(times)) .or. .not.all(ieee_is_finite(orbits))) return
    if(.not.all(ieee_is_finite(standard)) .or. .not.all(ieee_is_finite(observer))) return
    if(.not.ieee_is_finite(epoch) .or. .not.all(ieee_is_finite(scan_angle))) return
    if(abs(standard(2))>=90.d0 .or. standard(3)<0.d0) return
    if(any(orbits(:,1)<=0.d0) .or. any(orbits(:,2)<0.d0) .or. any(orbits(:,2)>=1.d0)) return
    if(any(orbits(:,7)<0.d0)) return
    ra=standard(1)*pi/180.d0
    dec=standard(2)*pi/180.d0
    do i=1,nobs
        ! Preserve ExTRA.standard_model(tangential=True), including its RA factor.
        pa=(-observer(i,1)*sin(ra)+observer(i,2)*cos(ra))/cos(dec)
        pd=(-observer(i,1)*cos(ra)-observer(i,2)*sin(ra))*sin(dec)+observer(i,3)*cos(dec)
        east(i)=standard(4)*(times(i)-epoch)/365.25d0+standard(3)*pa
        north(i)=standard(5)*(times(i)-epoch)/365.25d0+standard(3)*pd
    enddo
    do j=1,npl
        p=orbits(j,1)
        ecc=orbits(j,2)
        w=orbits(j,3)
        inc=orbits(j,4)
        node=orbits(j,5)
        axis=orbits(j,7)
        aa=axis*(cos(w)*cos(node)-sin(w)*sin(node)*cos(inc))
        bb=axis*(cos(w)*sin(node)+sin(w)*cos(node)*cos(inc))
        ff=axis*(-sin(w)*cos(node)-cos(w)*sin(node)*cos(inc))
        gg=axis*(-sin(w)*sin(node)+cos(w)*cos(node)*cos(inc))
        do i=1,nobs
            ma=modulo(2.d0*pi*modulo(times(i)-orbits(j,6),p)/p+pi,2.d0*pi)-pi
            ea=ma
            lo=-pi
            hi=pi
            do k=1,100
                f=ea-ecc*sin(ea)-ma
                if(abs(f)<2.d-15) exit
                if(f<0.d0) then
                    lo=ea
                else
                    hi=ea
                endif
                trial=ea-f/(1.d0-ecc*cos(ea))
                if(trial<=lo .or. trial>=hi) trial=0.5d0*(lo+hi)
                ea=trial
            enddo
            if(k>100) then
                status=2
                return
            endif
            x=cos(ea)-ecc
            y=sqrt(1.d0-ecc**2)*sin(ea)
            parts(i,j,:)=[bb*x+gg*y,aa*x+ff*y]
            east(i)=east(i)+parts(i,j,1)
            north(i)=north(i)+parts(i,j,2)
        enddo
    enddo
    ! Gaia position angle measured from north towards east, in radians.
    along_scan=east*sin(scan_angle)+north*cos(scan_angle)
    status=0
end subroutine astrometry_model

! Search foreground conjunctions independently of the RV/flux sampling grid.
! Empty observation arrays stay inside Fortran for f2py compatibility.
subroutine timing_model(elements,masses,stellar_mass,epoch,start_time,end_time,mode,axes, &
        step,tolerance,max_steps,jacobi_mode,samples_per_orbit,max_samples,max_events,counts,events,status,npl)
    implicit none
    integer,intent(in) :: npl,mode,max_steps,jacobi_mode,samples_per_orbit,max_samples,max_events
    real(8),intent(in) :: elements(npl,6),masses(npl),stellar_mass,epoch,start_time,end_time,axes(npl)
    real(8),intent(in) :: step,tolerance
    integer,intent(out) :: counts(npl),status
    real(8),intent(out) :: events(max_events,npl,4)
    integer :: empty_int(0)
    real(8) :: empty_real(0),prediction(0),residual(0),point_loglik(0),by_dataset(10),loglik
    call joint_timing_likelihood(elements,masses,stellar_mass,epoch,start_time,end_time,mode,axes, &
        step,tolerance,max_steps,jacobi_mode,samples_per_orbit,max_samples,max_events, &
        empty_int,empty_int,empty_int,empty_real,empty_real,counts,events, &
        prediction,residual,point_loglik,by_dataset,loglik,status,0,npl)
end subroutine timing_model

subroutine astrometry_nbody(times,elements,masses,axes,stellar_mass,standard,observer,epoch, &
        scan_angle,step,tolerance,max_steps,jacobi_mode,east,north,along_scan,parts,status,nobs,npl)
    implicit none
    integer,intent(in) :: nobs,npl,max_steps,jacobi_mode
    real(8),intent(in) :: times(nobs),elements(npl,6),masses(npl),axes(npl),stellar_mass
    real(8),intent(in) :: standard(5),observer(nobs,3),epoch,scan_angle(nobs),step,tolerance
    real(8),intent(out) :: east(nobs),north(nobs),along_scan(nobs),parts(nobs,npl,2)
    integer,intent(out) :: status
    real(8) :: geometry(nobs,npl,6),empty_orbits(0,7),empty_parts(nobs,0,2),factor
    integer :: j
    parts=0.d0
    call astrometry_model(times,empty_orbits,standard,observer,epoch,scan_angle, &
        east,north,along_scan,empty_parts,status,nobs,0)
    if(status/=0) return
    ! Tiny normalization radius avoids imposing a solar-radius collision limit.
    call rv_transit_geometry_core(times-epoch,elements,masses,axes,stellar_mass,1.d-6, &
        1,step,tolerance,max_steps,jacobi_mode,geometry,status,nobs,npl)
    if(status/=0) return
    do j=1,npl
        ! Heliocentric positions -> barycentric stellar reflex; R_sun -> AU.
        factor=-masses(j)/(stellar_mass+sum(masses))*1.d-6*6.957d8/1.49597892d11*standard(3)
        parts(:,j,1)=geometry(:,j,2)*factor
        parts(:,j,2)=geometry(:,j,1)*factor
        east=east+parts(:,j,1)
        north=north+parts(:,j,2)
    enddo
    along_scan=east*sin(scan_angle)+north*cos(scan_angle)
end subroutine astrometry_nbody

! Complete geometric transit diagnostics, independent of the plotting grid.
! The output time columns remain relative to epoch to retain duration precision.
subroutine transit_event_model(elements,masses,axes,radii,stellar_mass,stellar_radius, &
        epoch,start_time,end_time,mode,step,tolerance,max_steps,jacobi_mode, &
        samples_per_orbit,max_samples,max_events,law,coeff,nquad,counts,events,status,npl)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    integer,intent(in) :: npl,mode,max_steps,jacobi_mode,samples_per_orbit,max_samples,max_events,law,nquad
    real(8),intent(in) :: elements(npl,6),masses(npl),axes(npl),radii(npl),stellar_mass,stellar_radius
    real(8),intent(in) :: epoch,start_time,end_time,step,tolerance,coeff(4)
    integer,intent(out) :: counts(npl),status
    ! midpoint, contacts 1..4, T14, T23, impact, sky speed, depth, mid-depth, depth time.
    real(8),intent(out) :: events(max_events,npl,12)
    real(8),allocatable :: grid(:),geometry(:,:,:)
    real(8) :: cadence,padding,lo,hi,span,required,h,left,right,middle,g0,g1,g,root
    real(8) :: r(3),v(3),contact_times(4),impact,nan,minimum,peak,probe,value,aa,bb,c,d,fc,fd
    real(8) :: norm,mid_flux
    integer :: ngrid,i,j,k,n,iteration,allocation_status
    counts=0
    events=0.d0
    status=1
    if(npl<1 .or. npl>9 .or. max_events<1 .or. max_samples<4 .or. samples_per_orbit<32) return
    if(.not.all(ieee_is_finite(elements)) .or. .not.all(ieee_is_finite(radii))) return
    if(.not.all(ieee_is_finite([epoch,start_time,end_time,stellar_radius]))) return
    if(end_time<=start_time .or. stellar_radius<=0.d0 .or. any(radii<0.d0)) return
    if(any(elements(:,1)<=0.d0) .or. any(elements(:,2)<0.d0) .or. any(elements(:,2)>=1.d0)) return
    if(nquad<16 .or. nquad>1024) return
    call rv_transit_profile(law,coeff,norm,status)
    if(status/=0) return
    cadence=minval(elements(:,1)*(1.d0-elements(:,2))**1.5d0/sqrt(1.d0+elements(:,2))) &
        /dble(samples_per_orbit)
    padding=maxval(elements(:,1))/2.d0
    lo=start_time-epoch
    hi=end_time-epoch
    span=hi-lo+2.d0*padding
    required=span/cadence+2.d0
    status=9
    if(.not.ieee_is_finite(required) .or. required>max_samples .or. required*npl*6.d0>64000000.d0) return
    ngrid=ceiling(span/cadence)+1
    allocate(grid(ngrid),geometry(ngrid,npl,6),stat=allocation_status)
    if(allocation_status/=0) return
    do i=1,ngrid
        grid(i)=lo-padding+span*dble(i-1)/dble(ngrid-1)
    enddo
    call rv_transit_geometry_core(grid,elements,masses,axes,stellar_mass,stellar_radius, &
        mode,step,tolerance,max_steps,jacobi_mode,geometry,status,ngrid,npl)
    if(status/=0) return
    nan=ieee_value(0.d0,ieee_quiet_nan)
    do j=1,npl
        if(radii(j)==0.d0) cycle
        do i=1,ngrid-1
            g0=dot_product(geometry(i,j,1:2),geometry(i,j,4:5))
            g1=dot_product(geometry(i+1,j,1:2),geometry(i+1,j,4:5))
            if(g0>0.d0 .or. g1<=0.d0) cycle
            h=grid(i+1)-grid(i)
            left=grid(i)
            right=grid(i+1)
            do iteration=1,64
                middle=0.5d0*(left+right)
                call interpolate(middle,r,v)
                g=dot_product(r(1:2),v(1:2))
                if(g<=0.d0) then
                    left=middle
                else
                    right=middle
                endif
                if(right-left<1.d-12) exit
            enddo
            root=0.5d0*(left+right)
            call interpolate(root,r,v)
            impact=sqrt(dot_product(r(1:2),r(1:2)))
            if(root<lo-1.d-10 .or. root>hi+1.d-10 .or. r(3)<=0.d0 .or. impact>=1.d0+radii(j)) cycle
            n=counts(j)+1
            status=10
            if(n>max_events) return
            call contact(1.d0+radii(j),-1,contact_times(1))
            if(status/=0) return
            call contact(1.d0+radii(j),1,contact_times(4))
            if(status/=0) return
            contact_times(2:3)=nan
            if(impact<abs(1.d0-radii(j))) then
                call contact(abs(1.d0-radii(j)),-1,contact_times(2))
                if(status/=0) return
                call contact(abs(1.d0-radii(j)),1,contact_times(3))
                if(status/=0) return
            endif
            mid_flux=flux_at(root)
            minimum=mid_flux
            peak=root
            ! Search the full event, including limb-brightened profiles.
            k=0
            do iteration=0,64
                probe=contact_times(1)+(contact_times(4)-contact_times(1))*dble(iteration)/64.d0
                value=flux_at(probe)
                if(value<minimum) then
                    minimum=value
                    peak=probe
                    k=iteration
                endif
            enddo
            aa=contact_times(1)+(contact_times(4)-contact_times(1))*dble(max(0,k-1))/64.d0
            bb=contact_times(1)+(contact_times(4)-contact_times(1))*dble(min(64,k+1))/64.d0
            c=bb-(bb-aa)*0.6180339887498949d0
            d=aa+(bb-aa)*0.6180339887498949d0
            fc=flux_at(c)
            fd=flux_at(d)
            do iteration=1,64
                if(bb-aa<1.d-12) exit
                if(fc<fd) then
                    bb=d
                    d=c
                    fd=fc
                    c=bb-(bb-aa)*0.6180339887498949d0
                    fc=flux_at(c)
                else
                    aa=c
                    c=d
                    fc=fd
                    d=aa+(bb-aa)*0.6180339887498949d0
                    fd=flux_at(d)
                endif
            enddo
            if(fc<minimum) then
                minimum=fc
                peak=c
            endif
            if(fd<minimum) then
                minimum=fd
                peak=d
            endif
            if(status/=0) return
            call interpolate(root,r,v)
            counts(j)=n
            events(n,j,:)=[root,contact_times,contact_times(4)-contact_times(1), &
                contact_times(3)-contact_times(2),impact,sqrt(dot_product(v(1:2),v(1:2))), &
                1.d0-minimum,1.d0-mid_flux,peak]
        enddo
    enddo
    status=0
contains
    subroutine interpolate(t,position,velocity)
        real(8),intent(in) :: t
        real(8),intent(out) :: position(3),velocity(3)
        real(8) :: u,u2,u3,dt,p0(3),p1(3),v0(3),v1(3)
        integer :: index
        index=max(1,min(ngrid-1,int((t-grid(1))/(grid(ngrid)-grid(1))*dble(ngrid-1))+1))
        dt=grid(index+1)-grid(index)
        u=(t-grid(index))/dt
        u2=u*u
        u3=u2*u
        p0=geometry(index,j,1:3)
        p1=geometry(index+1,j,1:3)
        v0=geometry(index,j,4:6)
        v1=geometry(index+1,j,4:6)
        position=(2*u3-3*u2+1)*p0+(u3-2*u2+u)*dt*v0+(-2*u3+3*u2)*p1+(u3-u2)*dt*v1
        velocity=((6*u2-6*u)*p0+(-6*u2+6*u)*p1)/dt+(3*u2-4*u+1)*v0+(3*u2-2*u)*v1
    end subroutine interpolate

    subroutine contact(limit,direction,tcontact)
        real(8),intent(in) :: limit
        integer,intent(in) :: direction
        real(8),intent(out) :: tcontact
        real(8) :: a,b,mid,position(3),velocity(3),edge
        integer :: index,iter
        status=11
        tcontact=nan
        index=i
        if(direction>0) index=i+1
        do while(index>=1 .and. index<=ngrid)
            edge=grid(index)
            if((edge-root)*direction>0.d0) then
                if(sqrt(sum(geometry(index,j,1:2)**2))>=limit) then
                    a=min(root,edge)
                    b=max(root,edge)
                    do iter=1,64
                        mid=0.5d0*(a+b)
                        call interpolate(mid,position,velocity)
                        if((sqrt(sum(position(1:2)**2))-limit)*direction>0.d0) then
                            b=mid
                        else
                            a=mid
                        endif
                        if(b-a<1.d-12) exit
                    enddo
                    tcontact=0.5d0*(a+b)
                    status=0
                    return
                endif
            endif
            index=index+direction
        enddo
    end subroutine contact

    function flux_at(t) result(value)
        real(8),intent(in) :: t
        real(8) :: value,position(3),velocity(3),distance(1),flux(1)
        integer :: flux_status
        call interpolate(t,position,velocity)
        distance(1)=sqrt(sum(position(1:2)**2))
        call rv_transit_event_flux(distance,radii(j),law,coeff,nquad,flux,flux_status,1)
        if(flux_status/=0) status=flux_status
        value=flux(1)
    end function flux_at
end subroutine transit_event_model
