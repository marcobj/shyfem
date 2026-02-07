module m_set_random_seed2
  !! High‑quality seed initialization suitable for double‑precision RNG use
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  implicit none
contains

  pure function mix64(x) result(y)
    !! SplitMix64 mixing function (fast, good diffusion)
    integer(int64), intent(inout) :: x
    integer(int64) :: y

    x = x + int(z'9E3779B97F4A7C15', int64)
    y = x
    y = ieor(y, ishft(y, -30)); y = y * int(z'BF58476D1CE4E5B9', int64)
    y = ieor(y, ishft(y, -27)); y = y * int(z'94D049BB133111EB', int64)
    y = ieor(y, ishft(y, -31))
  end function mix64


  subroutine set_random_seed2()
    !! Initialize the Fortran RNG with a strong seed suitable for double‑precision runs.
    integer :: sze, i
    integer, allocatable :: seed(:)
    integer :: vals(8)
    integer :: cnt, rate, cntmax
    integer(int64) :: state, acc

    ! Query required seed size
    call random_seed(size = sze)
    if (sze < 1) sze = 1      ! Robust fallback (rare)
    allocate(seed(sze))

    ! --- Collect entropy ---
    call date_and_time(values = vals)          ! vals: year,month,day,tz,h,m,s,ms
    call system_clock(cnt, rate, cntmax)

    ! --- Build a 64-bit entropy accumulator ---
    acc = int(vals(1), int64)
    acc = ieor(acc, ishft(int(vals(2),int64),  6))   ! month
    acc = ieor(acc, ishft(int(vals(3),int64), 12))   ! day
    acc = ieor(acc, ishft(int(vals(5),int64), 18))   ! hour
    acc = ieor(acc, ishft(int(vals(6),int64), 24))   ! minute
    acc = ieor(acc, ishft(int(vals(7),int64), 30))   ! sec
    acc = ieor(acc, ishft(int(vals(8),int64), 36))   ! msec
    acc = ieor(acc, ishft(int(cnt   ,int64),  1))    ! cycle counter
    acc = ieor(acc, ishft(int(rate  ,int64), 43))    ! counter rate
    acc = ieor(acc, ishft(int(cntmax,int64), 51))    ! counter range

    if (acc == 0_int64) acc = 1_int64   ! avoid degenerate state

    ! --- Produce as many 32‑bit seeds as required ---
    state = acc
    do i = 1, sze
       state = mix64(state)
       seed(i) = int( iand(state, int(huge(0_int32), int64)), int32 )
       seed(i) = abs(seed(i)) + 1       ! avoid zeros
    end do

    call random_seed(put = seed)
    deallocate(seed)
  end subroutine set_random_seed2

end module m_set_random_seed2
