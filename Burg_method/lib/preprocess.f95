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
subroutine detrend_linear(x)

    real(mp), intent(inout) :: x(:)

    integer                 :: n, i
    real(mp)                :: sx, si, sii, six
    real(mp)                :: a, b
    real(mp)                :: denom

    n = size(x)

    sx  = 0.0_mp
    si  = 0.0_mp
    sii = 0.0_mp
    six = 0.0_mp

    do i = 1, n
        sx  = sx  + x(i)
        si  = si  + real(i, mp)
        sii = sii + real(i*i, mp)
        six = six + real(i, mp) * x(i)
    end do

    denom = real(n,mp)*sii - si*si

    if (abs(denom) < 1.0e-12_mp) return

    a = (real(n,mp)*six - si*sx) / denom
    b = (sx - a*si) / real(n,mp)

    do i = 1, n
        x(i) = x(i) - (a*real(i,mp) + b)
    end do

end subroutine detrend_linear


!=========================================================
! Make complex series z = x - i*y
!=========================================================
subroutine make_complex_series(x, y, z)

    real(mp), intent(in)  :: x(:), y(:)
    complex(mp), allocatable, intent(out) :: z(:)

    integer :: n

    n = size(x)
    if (size(y) /= n) error stop "make_complex_series: size mismatch"

    allocate(z(n))
    z = cmplx(x, -1_mp*y, kind=mp)

end subroutine make_complex_series

end module preprocess
