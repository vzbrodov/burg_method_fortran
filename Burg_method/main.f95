program main
    use my_prec
    use preprocess
    implicit none
    real(mp), allocatable :: x(:), y(:)
    integer :: n, i

    n = 10
    allocate(x(n), y(n))

    do i = 1, n
        x(i) = 2.0_mp*i + 5.0_mp
        y(i) = -1.0_mp*i + 3.0_mp
    end do

    print *, "Original x:"
    print *, x
    print *, "Original y:"
    print *, y

    call center_series(x)
    call center_series(y)

    print *, "After centering x:"
    print *, x
    print *, "After centering y:"
    print *, y

    ! --- линейный детренд
    call detrend_linear(x)
    call detrend_linear(y)

    print *, "After detrend x:"
    print *, x
    print *, "After detrend y:"
    print *, y

end program main
