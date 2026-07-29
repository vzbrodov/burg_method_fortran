module burg_method
  use, intrinsic :: iso_fortran_env, only : real32, real64
  implicit none
  private

  integer, parameter, public :: sp = real32
  integer, parameter, public :: dp = real64

  integer, parameter, public :: BURG_SUCCESS = 0
  integer, parameter, public :: BURG_INVALID_SIZE = 1
  integer, parameter, public :: BURG_INVALID_ORDER = 2
  integer, parameter, public :: BURG_INVALID_VARIANCE = 3
  integer, parameter, public :: BURG_DEGENERATE_SERIES = 4
  integer, parameter, public :: BURG_UNSTABLE_MODEL = 5

  public :: burg_fit
  public :: burg_spectrum

  interface burg_fit
    module procedure burg_fit_sp
    module procedure burg_fit_dp
  end interface burg_fit

  interface burg_spectrum
    module procedure burg_spectrum_sp
    module procedure burg_spectrum_dp
  end interface burg_spectrum

contains

  ! Fit a complex autoregressive model with Burg's recursion.
  !
  ! The returned coefficients define
  !
  !   x(n) + a(1) x(n-1) + ... + a(order) x(n-order) = e(n).
  !
  ! If STAT is omitted, invalid input or a numerical failure terminates the
  ! program with ERROR STOP. If STAT is present, the routine returns one of
  ! the public BURG_* status constants instead.
  subroutine burg_fit_sp(samples, order, coefficients, variance, stat)
    complex(sp), intent(in) :: samples(:)
    integer, intent(in) :: order
    complex(sp), allocatable, intent(out) :: coefficients(:)
    real(sp), intent(out) :: variance
    integer, optional, intent(out) :: stat

    complex(sp), allocatable :: forward_error(:), backward_error(:)
    complex(sp), allocatable :: ar(:)
    complex(sp) :: reflection, old_left, new_forward, new_backward
    real(sp) :: denominator
    integer :: i, k, j, status

    status = validate_fit_input(size(samples), order)
    if (present(stat)) stat = status
    variance = 0.0_sp
    if (status /= BURG_SUCCESS) then
      if (present(stat)) return
      if (status == BURG_INVALID_SIZE) then
        error stop 'burg_fit: at least one sample is required'
      else
        error stop 'burg_fit: order must satisfy 0 <= order < size(samples)'
      end if
    end if

    allocate(forward_error(size(samples)), backward_error(size(samples)))
    allocate(ar(0:order))
    forward_error = samples
    backward_error = samples
    ar = cmplx(0.0_sp, 0.0_sp, kind=sp)
    ar(0) = cmplx(1.0_sp, 0.0_sp, kind=sp)

    variance = sum(abs(samples)**2) / real(size(samples), sp)
    if (order == 0) then
      allocate(coefficients(0))
      return
    end if

    denominator = 0.0_sp
    do i = 1, size(samples)
      denominator = denominator + 2.0_sp * real( &
        forward_error(i) * conjg(forward_error(i)), sp)
    end do
    denominator = denominator &
      - real(forward_error(1) * conjg(forward_error(1)), sp) &
      - real(backward_error(size(samples)) &
        * conjg(backward_error(size(samples))), sp)

    do k = 0, order - 1
      if (denominator <= tiny(denominator)) then
        if (present(stat)) then
          stat = BURG_DEGENERATE_SERIES
          variance = 0.0_sp
          return
        end if
        error stop 'burg_fit: the prediction-error denominator is zero'
      end if

      reflection = cmplx(0.0_sp, 0.0_sp, kind=sp)
      do i = 1, size(samples) - k - 1
        reflection = reflection + forward_error(i+k+1) &
          * conjg(backward_error(i))
      end do
      reflection = -2.0_sp * reflection / denominator

      if (abs(reflection) >= 1.0_sp) then
        if (present(stat)) then
          stat = BURG_UNSTABLE_MODEL
          variance = 0.0_sp
          return
        end if
        error stop 'burg_fit: unstable reflection coefficient'
      end if

      variance = (1.0_sp - abs(reflection)**2) * variance

      ! Update each symmetric coefficient pair exactly once. This ordering
      ! follows Marple's reference implementation.
      do j = 1, (k + 1) / 2
        old_left = ar(j)
        ar(j) = old_left + reflection * conjg(ar(k+1-j))
        if (j /= k + 1 - j) then
          ar(k+1-j) = ar(k+1-j) + reflection * conjg(old_left)
        end if
      end do
      ar(k+1) = reflection

      do i = 1, size(samples) - k - 1
        new_forward = forward_error(i+k+1) &
          + reflection * backward_error(i)
        new_backward = backward_error(i) &
          + conjg(reflection) * forward_error(i+k+1)
        forward_error(i+k+1) = new_forward
        backward_error(i) = new_backward
      end do

      denominator = (1.0_sp - real(reflection*conjg(reflection), sp)) &
        * denominator &
        - real(forward_error(k+2) * conjg(forward_error(k+2)), sp) &
        - real(backward_error(size(samples)-k-1) &
          * conjg(backward_error(size(samples)-k-1)), sp)
    end do

    allocate(coefficients(order))
    coefficients = ar(1:order)
  end subroutine burg_fit_sp


  ! Double-precision implementation of BURG_FIT.
  subroutine burg_fit_dp(samples, order, coefficients, variance, stat)
    complex(dp), intent(in) :: samples(:)
    integer, intent(in) :: order
    complex(dp), allocatable, intent(out) :: coefficients(:)
    real(dp), intent(out) :: variance
    integer, optional, intent(out) :: stat

    complex(dp), allocatable :: forward_error(:), backward_error(:)
    complex(dp), allocatable :: ar(:)
    complex(dp) :: reflection, old_left, new_forward, new_backward
    real(dp) :: denominator
    integer :: i, k, j, status

    status = validate_fit_input(size(samples), order)
    if (present(stat)) stat = status
    variance = 0.0_dp
    if (status /= BURG_SUCCESS) then
      if (present(stat)) return
      if (status == BURG_INVALID_SIZE) then
        error stop 'burg_fit: at least one sample is required'
      else
        error stop 'burg_fit: order must satisfy 0 <= order < size(samples)'
      end if
    end if

    allocate(forward_error(size(samples)), backward_error(size(samples)))
    allocate(ar(0:order))
    forward_error = samples
    backward_error = samples
    ar = cmplx(0.0_dp, 0.0_dp, kind=dp)
    ar(0) = cmplx(1.0_dp, 0.0_dp, kind=dp)

    variance = sum(abs(samples)**2) / real(size(samples), dp)
    if (order == 0) then
      allocate(coefficients(0))
      return
    end if

    denominator = 0.0_dp
    do i = 1, size(samples)
      denominator = denominator + 2.0_dp * real( &
        forward_error(i) * conjg(forward_error(i)), dp)
    end do
    denominator = denominator &
      - real(forward_error(1) * conjg(forward_error(1)), dp) &
      - real(backward_error(size(samples)) &
        * conjg(backward_error(size(samples))), dp)

    do k = 0, order - 1
      if (denominator <= tiny(denominator)) then
        if (present(stat)) then
          stat = BURG_DEGENERATE_SERIES
          variance = 0.0_dp
          return
        end if
        error stop 'burg_fit: the prediction-error denominator is zero'
      end if

      reflection = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do i = 1, size(samples) - k - 1
        reflection = reflection + forward_error(i+k+1) &
          * conjg(backward_error(i))
      end do
      reflection = -2.0_dp * reflection / denominator

      if (abs(reflection) >= 1.0_dp) then
        if (present(stat)) then
          stat = BURG_UNSTABLE_MODEL
          variance = 0.0_dp
          return
        end if
        error stop 'burg_fit: unstable reflection coefficient'
      end if

      variance = (1.0_dp - abs(reflection)**2) * variance

      do j = 1, (k + 1) / 2
        old_left = ar(j)
        ar(j) = old_left + reflection * conjg(ar(k+1-j))
        if (j /= k + 1 - j) then
          ar(k+1-j) = ar(k+1-j) + reflection * conjg(old_left)
        end if
      end do
      ar(k+1) = reflection

      do i = 1, size(samples) - k - 1
        new_forward = forward_error(i+k+1) &
          + reflection * backward_error(i)
        new_backward = backward_error(i) &
          + conjg(reflection) * forward_error(i+k+1)
        forward_error(i+k+1) = new_forward
        backward_error(i) = new_backward
      end do

      denominator = (1.0_dp - real(reflection*conjg(reflection), dp)) &
        * denominator &
        - real(forward_error(k+2) * conjg(forward_error(k+2)), dp) &
        - real(backward_error(size(samples)-k-1) &
          * conjg(backward_error(size(samples)-k-1)), dp)
    end do

    allocate(coefficients(order))
    coefficients = ar(1:order)
  end subroutine burg_fit_dp


  ! Evaluate the MEM power spectral density at normalized frequencies.
  !
  ! FREQUENCY is measured in cycles per sample. For complex data it is
  ! usually appropriate to evaluate the complete interval [-0.5, 0.5).
  subroutine burg_spectrum_sp(coefficients, variance, frequency, power, stat)
    complex(sp), intent(in) :: coefficients(:)
    real(sp), intent(in) :: variance
    real(sp), intent(in) :: frequency(:)
    real(sp), intent(out) :: power(:)
    integer, optional, intent(out) :: stat

    complex(sp) :: denominator
    real(sp) :: two_pi
    integer :: i, k

    if (present(stat)) stat = BURG_SUCCESS
    power = 0.0_sp
    if (size(power) /= size(frequency)) then
      if (present(stat)) then
        stat = BURG_INVALID_SIZE
        return
      end if
      error stop 'burg_spectrum: frequency and power sizes differ'
    end if
    if (variance < 0.0_sp) then
      if (present(stat)) then
        stat = BURG_INVALID_VARIANCE
        power = 0.0_sp
        return
      end if
      error stop 'burg_spectrum: variance must be non-negative'
    end if

    two_pi = 2.0_sp * acos(-1.0_sp)
    do i = 1, size(frequency)
      denominator = cmplx(1.0_sp, 0.0_sp, kind=sp)
      do k = 1, size(coefficients)
        denominator = denominator + coefficients(k) * exp(cmplx( &
          0.0_sp, -two_pi * frequency(i) * real(k, sp), kind=sp))
      end do
      power(i) = variance / abs(denominator)**2
    end do
  end subroutine burg_spectrum_sp


  ! Double-precision implementation of BURG_SPECTRUM.
  subroutine burg_spectrum_dp(coefficients, variance, frequency, power, stat)
    complex(dp), intent(in) :: coefficients(:)
    real(dp), intent(in) :: variance
    real(dp), intent(in) :: frequency(:)
    real(dp), intent(out) :: power(:)
    integer, optional, intent(out) :: stat

    complex(dp) :: denominator
    real(dp) :: two_pi
    integer :: i, k

    if (present(stat)) stat = BURG_SUCCESS
    power = 0.0_dp
    if (size(power) /= size(frequency)) then
      if (present(stat)) then
        stat = BURG_INVALID_SIZE
        return
      end if
      error stop 'burg_spectrum: frequency and power sizes differ'
    end if
    if (variance < 0.0_dp) then
      if (present(stat)) then
        stat = BURG_INVALID_VARIANCE
        power = 0.0_dp
        return
      end if
      error stop 'burg_spectrum: variance must be non-negative'
    end if

    two_pi = 2.0_dp * acos(-1.0_dp)
    do i = 1, size(frequency)
      denominator = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do k = 1, size(coefficients)
        denominator = denominator + coefficients(k) * exp(cmplx( &
          0.0_dp, -two_pi * frequency(i) * real(k, dp), kind=dp))
      end do
      power(i) = variance / abs(denominator)**2
    end do
  end subroutine burg_spectrum_dp


  integer function validate_fit_input(sample_count, order) result(status)
    integer, intent(in) :: sample_count, order

    if (sample_count < 1) then
      status = BURG_INVALID_SIZE
    else if (order < 0 .or. order >= sample_count) then
      status = BURG_INVALID_ORDER
    else
      status = BURG_SUCCESS
    end if
  end function validate_fit_input

end module burg_method
