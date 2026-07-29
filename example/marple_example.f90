program marple_example
  use burg_method, only : sp, burg_fit, burg_spectrum
  implicit none

  integer, parameter :: ar_order = 15
  character(len=512) :: input_path
  complex(sp), allocatable :: samples(:), coefficients(:)
  real(sp), allocatable :: frequency(:), power(:)
  real(sp) :: variance
  integer :: i, unit

  input_path = 'example/data/marple_test64.dat'
  if (command_argument_count() >= 1) call get_command_argument(1, input_path)

  call read_complex_series(trim(input_path), samples)
  if (size(samples) <= ar_order) then
    error stop 'marple_example: the input must contain more than 15 samples'
  end if

  call burg_fit(samples, ar_order, coefficients, variance)

  ! A complex-valued series generally has an asymmetric spectrum, so the
  ! example evaluates the complete normalized-frequency interval.
  allocate(frequency(size(samples)), power(size(samples)))
  do i = 1, size(samples)
    frequency(i) = -0.5_sp + real(i-1, sp) / real(size(samples), sp)
  end do
  call burg_spectrum(coefficients, variance, frequency, power)

  open(newunit=unit, file='marple_coefficients.txt', status='replace', &
    action='write')
  write(unit, '(A)') '# lag  real(a)  imag(a)'
  write(unit, '(I4,2(1X,ES24.16))') 0, 1.0_sp, 0.0_sp
  do i = 1, size(coefficients)
    write(unit, '(I4,2(1X,ES24.16))') i, real(coefficients(i)), &
      aimag(coefficients(i))
  end do
  close(unit)

  open(newunit=unit, file='marple_spectrum.txt', status='replace', &
    action='write')
  write(unit, '(A)') '# normalized_frequency  power'
  do i = 1, size(frequency)
    write(unit, '(2(1X,ES24.16))') frequency(i), power(i)
  end do
  close(unit)

  print '(A,I0)', 'Samples: ', size(samples)
  print '(A,I0)', 'AR order: ', ar_order
  print '(A,ES14.6)', 'Residual variance: ', variance
  print '(A)', 'Coefficients: marple_coefficients.txt'
  print '(A)', 'Spectrum: marple_spectrum.txt'

contains

  subroutine read_complex_series(filename, values)
    character(len=*), intent(in) :: filename
    complex(sp), allocatable, intent(out) :: values(:)

    character(len=512) :: line
    complex(sp) :: value
    integer :: unit, io_status, count, i

    open(newunit=unit, file=filename, status='old', action='read', &
      iostat=io_status)
    if (io_status /= 0) error stop 'Cannot open the example data file'

    count = 0
    do
      read(unit, '(A)', iostat=io_status) line
      if (io_status /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=io_status) value
      if (io_status /= 0) error stop 'Invalid complex value in example data'
      count = count + 1
    end do
    if (count == 0) error stop 'The example data file is empty'

    rewind(unit)
    allocate(values(count))
    i = 0
    do
      read(unit, '(A)', iostat=io_status) line
      if (io_status /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=io_status) value
      if (io_status /= 0) error stop 'Invalid complex value in example data'
      i = i + 1
      values(i) = value
    end do
    close(unit)
  end subroutine read_complex_series

end program marple_example
