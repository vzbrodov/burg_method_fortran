module mem_lib
  use my_prec
  implicit none

contains

subroutine burg_mem(x, N, m, power, var_out, a_coeffs)
  integer, intent(in)                   :: N, m
  complex(mp), intent(in)               :: x(N)
  real(mp), intent(out)                 :: power(N/2)
  real(mp), intent(out)                 :: var_out
  complex(mp), intent(out), optional    :: a_coeffs(m)

  integer                               :: k, i, j
  complex(mp)                           :: mu, t1, t2
  complex(mp), allocatable              :: f(:), b(:), Ak(:)
  real(mp)                              :: Dk, rho

  if (N < 2) error stop "burg_mem: at least two samples are required"
  if (m < 1 .or. m >= N) error stop "burg_mem: order must satisfy 1 <= m < N"
  if (size(power) /= N/2) error stop "burg_mem: power must have size N/2"

  allocate(f(N), b(N), Ak(0:m))

  ! --- initialization ---
  Ak = (0.0_mp, 0.0_mp)
  Ak(0) = (1.0_mp, 0.0_mp)

  f = x
  b = x
  rho = sum(abs(x)**2) / real(N, mp)

  ! --- initial Dk ---
  Dk = 0.0_mp
  do i = 1, N
     Dk = Dk + 2.0_mp * real( f(i) * conjg(f(i)) )
  end do
  Dk = Dk - real( f(1) * conjg(f(1)) ) &
           - real( b(N) * conjg(b(N)) )

  ! ======================
  ! === Burg recursion ===
  ! ======================
  do k = 0, m-1

     mu = (0.0_mp, 0.0_mp)
     do i = 1, N-k-1
        mu = mu + f(i+k+1) * conjg(b(i))
     end do
     if (Dk <= tiny(Dk)) error stop "burg_mem: non-positive denominator"
     mu = -2.0_mp * mu / Dk
     if (abs(mu) >= 1.0_mp) error stop "burg_mem: unstable reflection coefficient"
     rho = (1.0_mp - abs(mu)**2) * rho

     ! Levinson recursion.  Update each symmetric pair exactly once;
     ! this ordering also follows Marple's reference implementation.
     do j = 1, (k+1)/2
        t1 = Ak(j)
        Ak(j) = t1 + mu * conjg(Ak(k+1-j))
        if (j /= k+1-j) then
           Ak(k+1-j) = Ak(k+1-j) + mu * conjg(t1)
        end if
     end do
     Ak(k+1) = mu

     do i = 1, N-k-1
        t1 = f(i+k+1) + mu * b(i)
        t2 = b(i)     + conjg(mu) * f(i+k+1)
        f(i+k+1) = t1
        b(i)     = t2
     end do

     Dk = (1.0_mp - real(mu*conjg(mu))) * Dk &
          - real( f(k+2) * conjg(f(k+2)) ) &
          - real( b(N-k-1) * conjg(b(N-k-1)) )

  end do

  if (present(a_coeffs)) then
     do i = 1, m
        a_coeffs(i) = Ak(i)
     end do
  end if

  var_out = rho
  power = 0.0_mp

  deallocate(f, b, Ak)

end subroutine burg_mem

subroutine mem_spectrum(a, m, var, power)
  use my_prec
  implicit none

  integer, intent(in)       :: m
  complex(mp), intent(in)   :: a(m)
  real(mp), intent(in)      :: var
  real(mp), intent(out)     :: power(:)

  integer                   :: i, nfreq
  real(mp), allocatable     :: frequency(:)

  nfreq = size(power)
  allocate(frequency(nfreq))

  do i = 1, nfreq
     frequency(i) = real(i-1, mp) / real(2*nfreq, mp)
  end do

  call mem_spectrum_grid(a, m, var, frequency, power)
  deallocate(frequency)

end subroutine mem_spectrum

subroutine mem_spectrum_grid(a, m, var, frequency, power)
  use my_prec
  implicit none

  integer, intent(in)       :: m
  complex(mp), intent(in)   :: a(m)
  real(mp), intent(in)      :: var
  real(mp), intent(in)      :: frequency(:)
  real(mp), intent(out)     :: power(:)

  integer                   :: i, k
  real(mp)                  :: twopi
  complex(mp)               :: denom, ex

  if (size(power) /= size(frequency)) then
     error stop "mem_spectrum_grid: frequency/power size mismatch"
  end if

  twopi = 2.0_mp * acos(-1.0_mp)

  do i = 1, size(frequency)
     denom = (1.0_mp, 0.0_mp)
     do k = 1, m
        ex = exp(cmplx(0.0_mp, -twopi*frequency(i)*real(k,mp), mp))
        denom = denom + a(k) * ex
     end do

     power(i) = var / abs(denom)**2
  end do

end subroutine mem_spectrum_grid
end module mem_lib
