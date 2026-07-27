# Burg Method for Complex Time Series in Fortran

## English

This repository contains a Fortran implementation of Burg's method for complex-valued time series and Maximum Entropy Method (MEM) spectral estimation.

The current example is focused on Earth orientation / polar motion analysis:

- reading `PM-X` and `PM-Y` from an EOP text file;
- centering and detrending the series;
- removing the annual harmonic;
- suppressing the Chandler component;
- building a complex signal;
- estimating AR coefficients with Burg's method;
- computing MEM spectra for positive and negative frequencies.

### Main files

- `Burg_method/main.f90` - example program that runs the full workflow
- `Burg_method/lib/my_prec.f90` - floating-point precision definition
- `Burg_method/lib/mem_lib.f90` - Burg AR estimation and MEM spectrum routines
- `Burg_method/lib/preprocess.f90` - centering, detrending, harmonic removal, and filtering
- `Burg_method/lib/eop_io.f90` - input/output for EOP polar motion data
- `Burg_method/lib/rawdata_io.f90` - additional text data readers


### Algorithm summary

The implementation follows this idea:

1. Convert two real-valued components into one complex series.
2. Estimate autoregressive coefficients using Burg's recursion.
3. Use the AR model and residual variance to evaluate the MEM spectrum.
4. Compute spectra separately for positive and negative frequencies.

### Project status

This repository looks like a research / working project rather than a packaged library. The README therefore describes the current structure and workflow of the example program, not a fully generalized API.

### Running the `test64` example

From the repository root:

```sh
gfortran -std=f2008 -Wall -Wextra -fcheck=all \
  Burg_method/lib/my_prec.f90 \
  Burg_method/lib/eop_io.f90 \
  Burg_method/lib/mem_lib.f90 \
  Burg_method/main.f90 \
  -o burg_test64
./burg_test64
```

The program reads 64 complex samples from `test64`, fits an AR(15)
model, and writes:

- `burg_coefficients.txt` — lag, real part, and imaginary part of each
  coefficient (including the fixed coefficient `a(0) = 1`);
- `burg_spectrum.txt` — the 64-point MEM spectrum on the normalized
  frequency grid `[-0.5, 0.5)` cycles per sample.

The `test64` values match Marple's published test sequence.  The
coefficients can differ slightly from the printed table because that table
was produced with single-precision complex arithmetic, while this project
uses double precision (`mp = 8`).

---

## Русский

Этот репозиторий содержит реализацию метода Бёрга на Fortran для комплекснозначных временных рядов и оценивания спектра методом максимальной энтропии (MEM).

Текущий пример ориентирован на анализ движения полюса / параметров ориентации Земли:

- чтение `PM-X` и `PM-Y` из текстового файла EOP;
- центрирование и удаление линейного тренда;
- удаление годовой гармоники;
- подавление чандлеровской компоненты;
- формирование комплексного ряда;
- оценивание коэффициентов AR методом Бёрга;
- вычисление MEM-спектров для положительных и отрицательных частот.

### Основные файлы

- `Burg_method/main.f90` - пример программы, выполняющий полный цикл обработки
- `Burg_method/lib/my_prec.f90` - задание рабочей точности
- `Burg_method/lib/mem_lib.f90` - оценка AR-коэффициентов и расчёт MEM-спектра
- `Burg_method/lib/preprocess.f90` - центрирование, детрендинг, удаление гармоник и фильтрация
- `Burg_method/lib/eop_io.f90` - ввод/вывод данных полюсного движения
- `Burg_method/lib/rawdata_io.f90` - дополнительные процедуры чтения текстовых данных

### Кратко об алгоритме

Идея реализации такая:

1. Две вещественные компоненты объединяются в один комплексный ряд.
2. Коэффициенты авторегрессии оцениваются рекурсией Бёрга.
3. По AR-модели и остаточной дисперсии вычисляется MEM-спектр.
4. Спектры считаются отдельно для положительных и отрицательных частот.

### Состояние проекта

По структуре это скорее исследовательский / рабочий проект, чем оформленная библиотека. Поэтому README описывает текущую организацию кода и сценарий запуска примера, а не универсальный публичный API.

### Запуск примера `test64`

Из корня репозитория:

```sh
gfortran -std=f2008 -Wall -Wextra -fcheck=all \
  Burg_method/lib/my_prec.f90 \
  Burg_method/lib/eop_io.f90 \
  Burg_method/lib/mem_lib.f90 \
  Burg_method/main.f90 \
  -o burg_test64
./burg_test64
```

Программа читает 64 комплексные точки из `test64`, строит модель AR(15)
и создаёт:

- `burg_coefficients.txt` — номер лага, действительная и мнимая части
  каждого коэффициента (включая фиксированный коэффициент `a(0) = 1`);
- `burg_spectrum.txt` — MEM-спектр из 64 точек на нормированной сетке
  частот `[-0.5, 0.5)` цикла на отсчёт.

Значения в `test64` соответствуют опубликованной тестовой
последовательности Марпла. Коэффициенты могут немного отличаться от
печатной таблицы: она была рассчитана с одинарной точностью комплексной
арифметики, а проект использует двойную точность (`mp = 8`).
