program test
  use my_prec
  use preprocess
  use eop_io
  use mem_lib
  implicit none

  integer :: n, i
  integer, parameter            :: m_test = 2500
  real(mp), allocatable         :: power(:), mjd(:), pmx(:), pmy(:)
  complex(mp), allocatable      :: z(:), a_coeffs(:)
  real(mp)                      :: var

  call read_eop_pm('code_p.eop.txt', mjd, pmx, pmy, n)
  print *, 'n=', n

  call center_series(pmx)
  call center_series(pmy)
  call detrend_linear_time(mjd, pmx)
  call detrend_linear_time(mjd, pmy)

  allocate(z(n), power(n/2), a_coeffs(m_test))
  call make_complex_series(pmx, pmy, z)


  ! --- 1. Burg coefficients ---
  call burg_mem(z, n, m_test, power, var, a_coeffs)


  ! --- 2. MEM spectrum ---
  call mem_spectrum(a_coeffs, m_test, var, power)


  ! --- 3. Output ---
  open(10,file='mem.txt',status='replace')
  do i = 1, n/2
     write(10,'(F10.6,ES14.6)') real(i-1,mp)/real(n,mp), power(i)
  end do
  close(10)

  deallocate(z, power, mjd, pmx, pmy, a_coeffs)
end program test
