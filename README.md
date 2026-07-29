# burg-method-fortran

[![CI](https://github.com/vzbrodov/burg_method_fortran/actions/workflows/ci.yml/badge.svg)](https://github.com/vzbrodov/burg_method_fortran/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/1163219836.svg)](https://doi.org/10.5281/zenodo.21683153)

A small Fortran library for Burg autoregressive modelling and Maximum
Entropy Method (MEM) spectral estimation of **complex-valued** time series.

[English](#english) · [Русский](#русский)

## English

### Features

- Burg estimation of complex AR models;
- MEM power spectra on arbitrary normalized-frequency grids;
- generic single- and double-precision APIs (`real32` and `real64`);
- explicit status codes for applications that must handle errors themselves;
- no external runtime dependencies;
- an `fpm` package, a tested reference example, and continuous integration.

The library contains signal-processing routines only. It makes no assumptions
about the physical origin of the samples.

### Quick start

Add the repository to an
[Fortran Package Manager](https://fpm.fortran-lang.org/) project:

```toml
[dependencies]
burg-method = { git = "https://github.com/vzbrodov/burg_method_fortran.git" }
```

The public API is provided by the `burg_method` module:

```fortran
use burg_method, only : dp, burg_fit, burg_spectrum

complex(dp)              :: samples(128)
complex(dp), allocatable :: coefficients(:)
real(dp)                 :: variance
real(dp)                 :: frequency(1024), power(1024)
integer                  :: i

! Fill SAMPLES with a uniformly sampled complex-valued series.

call burg_fit(samples, order=20, coefficients=coefficients, &
  variance=variance)

do i = 1, size(frequency)
  frequency(i) = -0.5_dp + real(i-1, dp) / real(size(frequency), dp)
end do
call burg_spectrum(coefficients, variance, frequency, power)
```

Frequencies are normalized in cycles per sample. For complex data, evaluate
the whole interval `[-0.5, 0.5)` because the spectrum need not be symmetric.
Multiply normalized frequencies by the sampling rate to obtain physical
frequency units.

`burg_fit` returns coefficients for the convention

```text
x(n) + a(1)x(n-1) + ... + a(m)x(n-m) = e(n),
```

and `variance` is the final prediction-error variance. `burg_spectrum`
evaluates

```text
P(f) = variance / |1 + Σ a(k) exp(-i 2π f k)|².
```

The generic procedures accept either `complex(sp)`/`real(sp)` or
`complex(dp)`/`real(dp)` arguments. All arguments in one call must have the
same precision.

### Error handling

Both public procedures accept an optional integer `stat` argument. When it is
present, they return one of:

| Constant | Meaning |
| --- | --- |
| `BURG_SUCCESS` | Successful calculation |
| `BURG_INVALID_SIZE` | Incompatible or empty arrays |
| `BURG_INVALID_ORDER` | AR order is outside `0 <= order < N` |
| `BURG_INVALID_VARIANCE` | Negative prediction-error variance |
| `BURG_DEGENERATE_SERIES` | Zero prediction-error denominator |
| `BURG_UNSTABLE_MODEL` | Reflection coefficient is not inside the unit circle |

Without `stat`, an invalid call terminates with a descriptive `error stop`.

### Marple reference example

The repository retains the classic 64-sample complex test sequence from
S. Lawrence Marple's
[*Digital Spectral Analysis: With Applications*](https://books.google.com/books?id=D-1QAAAAMAAJ)
(Prentice-Hall, 1987). The example fits an AR(15) model in single precision,
matching the arithmetic used for the published reference coefficients.

```sh
fpm run --example marple
```

The program writes:

- `marple_coefficients.txt` — AR coefficients, including `a(0) = 1`;
- `marple_spectrum.txt` — a 64-point MEM spectrum on `[-0.5, 0.5)`.

You may pass another file containing one Fortran-form complex value per line:

```sh
fpm run --example marple -- path/to/complex_series.dat
```

The bundled data and coefficients are also used as a regression test.

### Citation

If this library contributes to your research, please cite the archived
software release: [burg-method-fortran v0.1.0](https://doi.org/10.5281/zenodo.21683155).
GitHub also provides APA and BibTeX entries through the
**Cite this repository** menu.

### Build and test

With `fpm`:

```sh
fpm test
fpm run --example marple
```

With a Fortran 2008 compiler directly:

```sh
mkdir -p build
gfortran -std=f2008 -Wall -Wextra -fcheck=all -J build \
  src/burg_method.f90 test/test_burg_method.f90 \
  -o build/test_burg_method
./build/test_burg_method
```

### Repository layout

```text
src/                         Library source and public burg_method module
example/marple_example.f90   Complete command-line example
example/data/                Marple's 64-sample complex sequence
test/                        Regression and API tests
fpm.toml                     Fortran Package Manager manifest
```

## Русский

`burg-method-fortran` — компактная библиотека на Fortran для построения
авторегрессионных моделей методом Бёрга и оценки спектра методом максимальной
энтропии (MEM) для **комплекснозначных** временных рядов.

### Возможности

- оценка комплексных AR-моделей методом Бёрга;
- MEM-спектр на произвольной сетке нормированных частот;
- единый API для одинарной и двойной точности (`real32` и `real64`);
- явные коды состояния для обработки ошибок вызывающей программой;
- отсутствие внешних зависимостей времени выполнения;
- пакет `fpm`, проверяемый эталонный пример и непрерывное тестирование.

Библиотека содержит только алгоритмы обработки сигналов и ничего не
предполагает о физическом происхождении данных.

### Быстрый старт

Добавьте репозиторий как зависимость
[Fortran Package Manager](https://fpm.fortran-lang.org/):

```toml
[dependencies]
burg-method = { git = "https://github.com/vzbrodov/burg_method_fortran.git" }
```

Публичный API находится в модуле `burg_method`:

```fortran
use burg_method, only : dp, burg_fit, burg_spectrum

complex(dp)              :: samples(128)
complex(dp), allocatable :: coefficients(:)
real(dp)                 :: variance
real(dp)                 :: frequency(1024), power(1024)
integer                  :: i

! Заполните SAMPLES отсчётами равномерно дискретизированного комплексного ряда.

call burg_fit(samples, order=20, coefficients=coefficients, &
  variance=variance)

do i = 1, size(frequency)
  frequency(i) = -0.5_dp + real(i-1, dp) / real(size(frequency), dp)
end do
call burg_spectrum(coefficients, variance, frequency, power)
```

Частота задаётся в циклах на отсчёт. Для комплексного ряда следует
рассчитывать весь интервал `[-0.5, 0.5)`, поскольку его спектр в общем случае
несимметричен. Для перехода к физическим единицам умножьте нормированную
частоту на частоту дискретизации.

`burg_fit` использует соглашение

```text
x(n) + a(1)x(n-1) + ... + a(m)x(n-m) = e(n),
```

а `variance` содержит итоговую дисперсию ошибки предсказания.
`burg_spectrum` вычисляет

```text
P(f) = variance / |1 + Σ a(k) exp(-i 2π f k)|².
```

Обобщённые процедуры принимают аргументы `complex(sp)`/`real(sp)` или
`complex(dp)`/`real(dp)`. В одном вызове точность всех аргументов должна
совпадать.

### Обработка ошибок

Обе публичные процедуры принимают необязательный целочисленный аргумент
`stat`. Если он передан, процедура возвращает один из кодов:

| Константа | Значение |
| --- | --- |
| `BURG_SUCCESS` | Вычисление выполнено |
| `BURG_INVALID_SIZE` | Пустые или несовместимые массивы |
| `BURG_INVALID_ORDER` | Порядок AR не удовлетворяет `0 <= order < N` |
| `BURG_INVALID_VARIANCE` | Отрицательная дисперсия ошибки |
| `BURG_DEGENERATE_SERIES` | Нулевой знаменатель ошибки предсказания |
| `BURG_UNSTABLE_MODEL` | Коэффициент отражения вне единичного круга |

Если `stat` не передан, некорректный вызов завершается понятным сообщением
`error stop`.

### Эталонный пример Марпла

В репозитории сохранена классическая комплексная последовательность из
64 отсчётов из книги S. Lawrence Marple
[*Digital Spectral Analysis: With Applications*](https://books.google.com/books?id=D-1QAAAAMAAJ)
(Prentice-Hall, 1987). Пример строит модель AR(15) в одинарной точности,
чтобы соответствовать арифметике опубликованных эталонных коэффициентов.

```sh
fpm run --example marple
```

Программа создаёт:

- `marple_coefficients.txt` — коэффициенты AR, включая `a(0) = 1`;
- `marple_spectrum.txt` — MEM-спектр из 64 точек на интервале `[-0.5, 0.5)`.

Можно передать другой файл с одним комплексным числом в формате Fortran на
строку:

```sh
fpm run --example marple -- path/to/complex_series.dat
```

Встроенные данные и коэффициенты также используются в регрессионном тесте.

### Цитирование

Если библиотека использовалась в исследовании, пожалуйста, процитируйте
архивированный релиз:
[burg-method-fortran v0.1.0](https://doi.org/10.5281/zenodo.21683155).
Готовые записи в форматах APA и BibTeX доступны через меню GitHub
**Cite this repository**.

### Сборка и тестирование

Через `fpm`:

```sh
fpm test
fpm run --example marple
```

Напрямую компилятором с поддержкой Fortran 2008:

```sh
mkdir -p build
gfortran -std=f2008 -Wall -Wextra -fcheck=all -J build \
  src/burg_method.f90 test/test_burg_method.f90 \
  -o build/test_burg_method
./build/test_burg_method
```

### Структура репозитория

```text
src/                         Исходный код библиотеки и модуль burg_method
example/marple_example.f90   Полный пример командной строки
example/data/                Комплексная последовательность Марпла
test/                        Регрессионные тесты и тесты API
fpm.toml                     Манифест Fortran Package Manager
```
