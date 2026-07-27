program test64_burg
  use my_prec
  use eop_io
  use mem_lib
  implicit none

  integer                       :: n, i, unit
  integer, parameter            :: ar_order = 15
  real(mp), allocatable         :: power_work(:), power(:), frequency(:)
  complex(mp), allocatable      :: z(:), a_coeffs(:)
  real(mp)                      :: var

  call read_complex_series('test64', z)
  n = size(z)

  if (n /= 64) error stop "test64_burg: expected exactly 64 samples"

  allocate(power_work(n/2), power(n), frequency(n), a_coeffs(ar_order))

  call burg_mem(z, n, ar_order, power_work, var, a_coeffs)

  ! A full normalized-frequency grid is important for a complex series:
  ! its spectrum is generally not symmetric about zero.
  do i = 1, n
     frequency(i) = -0.5_mp + real(i-1, mp)/real(n, mp)
  end do
  call mem_spectrum_grid(a_coeffs, ar_order, var, frequency, power)

  open(newunit=unit, file='burg_coefficients.txt', status='replace', &
       action='write')
  write(unit, '(A)') '# lag  real(a)  imag(a)'
  write(unit, '(I4,2(1X,ES24.16))') 0, 1.0_mp, 0.0_mp
  do i = 1, ar_order
     write(unit, '(I4,2(1X,ES24.16))') i, real(a_coeffs(i)), &
                                      aimag(a_coeffs(i))
  end do
  close(unit)

  open(newunit=unit, file='burg_spectrum.txt', status='replace', &
       action='write')
  write(unit, '(A)') '# normalized_frequency  power'
  do i = 1, n
     write(unit, '(2(1X,ES24.16))') frequency(i), power(i)
  end do
  close(unit)

  print '(A,I0)', 'Samples: ', n
  print '(A,I0)', 'AR order: ', ar_order
  print '(A,ES14.6)', 'Residual variance: ', var
  print '(A)', 'Coefficients: burg_coefficients.txt'
  print '(A)', 'Spectrum: burg_spectrum.txt'

  deallocate(z, power_work, power, frequency, a_coeffs)
end program test64_burg
