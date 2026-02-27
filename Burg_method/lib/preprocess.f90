module preprocess
    use my_prec
    implicit none
contains

!=========================================================
! Centering: subtract mean
!=========================================================
subroutine center_series(x)

    real(mp), intent(inout) :: x(:)
    real(mp)                :: mean

    mean = sum(x) / real(size(x), mp)
    x = x - mean

end subroutine center_series


!=========================================================
! Remove linear trend using least squares
!
! x(i) -> x(i) - (a*i + b)
!=========================================================
subroutine detrend_linear_time(t, x)
    use my_prec
    real(mp), intent(in)    :: t(:)   ! MJD
    real(mp), intent(inout) :: x(:)

    integer :: n, i
    real(mp) :: st, sx, stt, stx
    real(mp) :: a, b, denom

    n = size(x)

    st  = 0.0_mp
    sx  = 0.0_mp
    stt = 0.0_mp
    stx = 0.0_mp

    do i = 1, n
        st  = st  + t(i)
        sx  = sx  + x(i)
        stt = stt + t(i)*t(i)
        stx = stx + t(i)*x(i)
    end do

    denom = real(n,mp)*stt - st*st
    if (abs(denom) < 1.0e-14_mp) return

    a = (real(n,mp)*stx - st*sx) / denom
    b = (sx - a*st) / real(n,mp)

    do i = 1, n
        x(i) = x(i) - (a*t(i) + b)
    end do
end subroutine detrend_linear_time


!=========================================================
! Make complex series z = x - i*y
!=========================================================
subroutine make_complex_series(x, y, z)

    real(mp), intent(in)                    :: x(:), y(:)
    complex(mp), allocatable, intent(out)   :: z(:)

    integer :: n

    n = size(x)
    if (size(y) /= n) error stop "make_complex_series: size mismatch"

    allocate(z(n))
    z = cmplx(x, -1.0_mp*y, kind=mp)
end subroutine make_complex_series

end module preprocess
