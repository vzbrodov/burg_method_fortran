# Burg’s Method (MEM) for Complex Time Series in Fortran

This project provides a **Fortran implementation of Burg’s Method** for estimating
**autoregressive (AR) models** and **Maximum Entropy Method (MEM) spectra**
from **complex-valued time series**.

The code is designed for scientific and geophysical applications where signals
are naturally complex (for example, two coupled real-valued time series).



## Features

- Burg’s algorithm for **complex-valued signals**
- Computes:
  - AR coefficients
  - MEM power spectrum
- Modular and easy to integrate


## File Structure


 `my_prec.f90` --> Floating-point precision definition (`mp`) 
 `mem_lib.f90` --> Burg algorithm and MEM spectrum routines 
 `preprocess.f90` --> Data preprocessing (centering, detrending)
 `eop_io.f90` --> Example input routines
 `main.f90` --> Example test program


## Usage

### 1. Prepare Complex Input Data

The input signal must be a **complex array**:

```fortran
complex(mp), allocatable :: z(:)
allocate(z(n))

! Fill z(:) with complex values
integer, parameter :: m = 20
complex(mp) :: a_coeffs(m)
real(mp) :: var
real(mp) :: power(n/2)

call fast_burg_mem(z, n, m, power, var, a_coeffs)
```

After this call:

- `a_coeffs(1:m)` contains the AR coefficients

- `var` contains the residual variance

### 2. Compute MEM Spectrum

```fortran
call mem_spectrum(a_coeffs, m, var, power)
f(i) = (i-1) / n      ! frequency range [0, 0.5)

```

### 3. Example Output

```fortran
print *, 'Residual variance:', var
print *, 'First AR coefficient:', a_coeffs(1)
print *, 'First 5 spectral values:', power(1:5)

```
