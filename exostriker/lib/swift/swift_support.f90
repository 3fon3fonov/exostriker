! Shared control for the serialized in-process API. No I/O or process exits.
module swift_memory_support
    use iso_fortran_env, only: int64
    implicit none
    integer :: swift_status = 0, ntracked = 0
    integer :: tracked_slot(50) = 0
    integer(int64) :: clock_start = 0, clock_rate = 0, clock_max = 0
    real(8) :: time_limit = 0.d0
contains
    subroutine memory_reset(n, timeout)
        integer, intent(in) :: n
        real(8), intent(in) :: timeout
        integer j
        swift_status = 0
        ntracked = n
        tracked_slot = 0
        do j=1,n
            tracked_slot(j) = j+1
        enddo
        time_limit = timeout
        call system_clock(clock_start, clock_rate, clock_max)
        if(clock_rate <= 0) swift_status = 4
    end subroutine

    subroutine swift_poll()
        integer(int64) now, elapsed
        if(swift_status /= 0) return
        call system_clock(now)
        elapsed = now-clock_start
        if(now < clock_start) elapsed = (clock_max-clock_start)+now+1_int64
        if(real(elapsed,8)/real(clock_rate,8) >= time_limit) swift_status = 2
    end subroutine

    subroutine swift_remove(ip)
        integer, intent(in) :: ip
        where(tracked_slot == ip) tracked_slot = 0
        where(tracked_slot > ip) tracked_slot = tracked_slot-1
    end subroutine

    subroutine swift_merge(survivor, removed)
        integer, intent(in) :: survivor, removed
        where(tracked_slot == removed) tracked_slot = survivor
    end subroutine
end module

! The numerical merger code calls these hooks. Bookkeeping is in memory;
! array compaction is reported separately by discard_mass_reorder5.
subroutine io_discard_merge(time,ip1,ip2,m1,r1,x1,y1,z1,vx1,vy1, &
        vz1,m2,r2,x2,y2,z2,vx2,vy2,vz2,mn,rn,xn,yn,zn,vxn,vyn,vzn)
    use swift_memory_support, only: swift_merge
    implicit none
    integer ip1,ip2
    real(8) time,m1,r1,x1,y1,z1,vx1,vy1,vz1,m2,r2,x2,y2,z2,vx2,vy2,vz2
    real(8) mn,rn,xn,yn,zn,vxn,vyn,vzn
    call swift_merge(ip1,ip2)
end subroutine

subroutine io_discard_mass(iflg,time,id,m1,r1,x1,y1,z1,vx1,vy1,vz1,iu,iwhy,fopenstat)
    implicit none
    integer iflg,id,iu,iwhy
    real(8) time,m1,r1,x1,y1,z1,vx1,vy1,vz1
    character*(*) fopenstat
    ! No file is opened. swift_remove is called when the arrays are compacted.
end subroutine
