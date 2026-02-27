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

  integer                               :: k, i, j, maxn
  complex(mp)                           :: mu, t1, t2
  complex(mp), allocatable              :: f(:), b(:), Ak(:)
  real(mp)                              :: Dk

  allocate(f(N), b(N), Ak(0:m))

  ! --- initialization ---
  Ak = (0.0_mp, 0.0_mp)
  Ak(0) = (1.0_mp, 0.0_mp)

  f = x
  b = x

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
     mu = -2.0_mp * mu / Dk

     maxn = (k+1)/2
     do j = 0, maxn
        t1 = Ak(j)       + mu * conjg(Ak(k+1-j))
        t2 = Ak(k+1-j)   + mu * conjg(Ak(j))
        Ak(j)       = t1
        Ak(k+1-j)   = t2
     end do

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

  var_out = Dk / real(N, mp)
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

  integer                   :: i, k, nfreq
  real(mp)                  :: f, twopi
  complex(mp)               :: denom, ex

  twopi = 2.0_mp * acos(-1.0_mp)
  nfreq = size(power)

  do i = 1, nfreq
     f = real(i-1, mp) / real(2*nfreq, mp)

     denom = (1.0_mp, 0.0_mp)
     do k = 1, m
        ex = exp( cmplx(0.0_mp, -twopi * f * real(k,mp), mp) )
        denom = denom + a(k) * ex
     end do

     power(i) = var / real( denom * conjg(denom) )
  end do

end subroutine mem_spectrum
end module mem_lib
