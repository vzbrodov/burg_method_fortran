program test_burg_method
  use burg_method
  implicit none

  complex(sp), allocatable :: samples_sp(:), coefficients_sp(:)
  complex(dp), allocatable :: samples_dp(:), coefficients_dp(:)
  complex(sp) :: reference(15)
  real(sp) :: variance_sp
  real(dp) :: variance_dp
  real(sp) :: frequency_sp(5), power_sp(5)
  real(dp) :: frequency_dp(5), power_dp(5)
  integer :: failures, stat

  failures = 0
  call read_complex_series('example/data/marple_test64.dat', samples_sp)
  call read_reference('test/reference/marple_ar15_coefficients.dat', &
    reference)

  call burg_fit(samples_sp, 15, coefficients_sp, variance_sp)
  call check(size(coefficients_sp) == 15, 'single-precision order')
  call check(maxval(abs(coefficients_sp-reference)) < 2.0e-5_sp, &
    'Marple single-precision coefficients')
  call check(abs(variance_sp-5.408684e-3_sp) < 1.0e-8_sp, &
    'Marple residual variance')

  frequency_sp = [-0.5_sp, -0.25_sp, 0.0_sp, 0.25_sp, 0.499_sp]
  call burg_spectrum(coefficients_sp, variance_sp, frequency_sp, power_sp)
  call check(all(power_sp > 0.0_sp), 'positive single-precision spectrum')

  allocate(samples_dp(size(samples_sp)))
  samples_dp = cmplx(real(samples_sp, dp), aimag(samples_sp), kind=dp)
  call burg_fit(samples_dp, 15, coefficients_dp, variance_dp)
  frequency_dp = real(frequency_sp, dp)
  call burg_spectrum(coefficients_dp, variance_dp, frequency_dp, power_dp)
  call check(size(coefficients_dp) == 15, 'double-precision order')
  call check(variance_dp > 0.0_dp, 'positive double-precision variance')
  call check(all(power_dp > 0.0_dp), 'positive double-precision spectrum')

  call burg_fit(samples_sp, size(samples_sp), coefficients_sp, variance_sp, &
    stat)
  call check(stat == BURG_INVALID_ORDER, 'invalid order status')

  call burg_fit([cmplx(0.0_sp, 0.0_sp, kind=sp), &
    cmplx(0.0_sp, 0.0_sp, kind=sp)], 1, coefficients_sp, variance_sp, stat)
  call check(stat == BURG_DEGENERATE_SERIES, 'degenerate series status')

  call burg_spectrum(reference, -1.0_sp, frequency_sp, power_sp, stat)
  call check(stat == BURG_INVALID_VARIANCE, 'invalid variance status')

  if (failures > 0) then
    print '(A,I0)', 'FAILED tests: ', failures
    error stop 1
  end if
  print '(A)', 'All burg_method tests passed.'

contains

  subroutine check(condition, name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name

    if (.not. condition) then
      failures = failures + 1
      print '(A)', 'FAIL: ' // name
    end if
  end subroutine check


  subroutine read_complex_series(filename, values)
    character(len=*), intent(in) :: filename
    complex(sp), allocatable, intent(out) :: values(:)

    character(len=512) :: line
    complex(sp) :: value
    integer :: unit, io_status, count, i

    open(newunit=unit, file=filename, status='old', action='read', &
      iostat=io_status)
    if (io_status /= 0) error stop 'Cannot open Marple test data'

    count = 0
    do
      read(unit, '(A)', iostat=io_status) line
      if (io_status /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=io_status) value
      if (io_status /= 0) error stop 'Invalid Marple test datum'
      count = count + 1
    end do

    rewind(unit)
    allocate(values(count))
    i = 0
    do
      read(unit, '(A)', iostat=io_status) line
      if (io_status /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=io_status) value
      if (io_status /= 0) error stop 'Invalid Marple test datum'
      i = i + 1
      values(i) = value
    end do
    close(unit)
  end subroutine read_complex_series


  subroutine read_reference(filename, values)
    character(len=*), intent(in) :: filename
    complex(sp), intent(out) :: values(:)

    character(len=512) :: line
    real(sp) :: real_part, imaginary_part
    integer :: unit, io_status, lag

    open(newunit=unit, file=filename, status='old', action='read', &
      iostat=io_status)
    if (io_status /= 0) error stop 'Cannot open reference coefficients'

    values = cmplx(0.0_sp, 0.0_sp, kind=sp)
    do
      read(unit, '(A)', iostat=io_status) line
      if (io_status /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=io_status) lag, real_part, imaginary_part
      if (io_status /= 0) error stop 'Invalid reference coefficient'
      if (lag >= 1 .and. lag <= size(values)) then
        values(lag) = cmplx(real_part, imaginary_part, kind=sp)
      end if
    end do
    close(unit)
  end subroutine read_reference

end program test_burg_method
