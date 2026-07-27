module eop_io
    use my_prec
    implicit none
contains

!=========================================================
! Read one complex value per line in Fortran (re, im) form
!=========================================================
subroutine read_complex_series(filename, z)

    character(len=*), intent(in)         :: filename
    complex(mp), allocatable, intent(out):: z(:)

    integer            :: unit, ios, count, i
    character(len=512) :: line
    complex(mp)        :: value

    count = 0
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop "read_complex_series: cannot open input file"

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        read(line, *, iostat=ios) value
        if (ios /= 0) error stop "read_complex_series: invalid complex value"
        count = count + 1
    end do

    if (count == 0) error stop "read_complex_series: no data"
    rewind(unit)
    allocate(z(count))

    i = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        read(line, *, iostat=ios) value
        if (ios /= 0) error stop "read_complex_series: invalid complex value"
        i = i + 1
        z(i) = value
    end do
    close(unit)

end subroutine read_complex_series


!=========================================================
! Read PM-X and PM-Y from IERS EOP file
!
! Input:
!   filename  - path to code_p.eop.txt
!
! Output:
!   mjd(:)
!   pmx(:)
!   pmy(:)
!   n         - number of data rows
!=========================================================
subroutine      read_eop_pm(filename, mjd, pmx, pmy, n)

    character(len=*), intent(in)        :: filename
    real(mp), allocatable, intent(out)  :: mjd(:), pmx(:), pmy(:)
    integer, intent(out)                :: n

    integer :: unit, ios, count
    character(len=512) :: line
    real(mp) :: tmp_mjd, tmp_pmx, tmp_pmy

    ! --- First pass: count valid data lines
    count = 0
    open(newunit=unit, file=filename, status='old', action='read')

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        if (line(1:1) /= '#') then
            if (len_trim(line) > 0) count = count + 1
        end if
    end do

    rewind(unit)

    allocate(mjd(count), pmx(count), pmy(count))

    ! --- Second pass: read data
    n = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        if (line(1:1) == '#') cycle
        if (len_trim(line) == 0) cycle

        read(line, *) tmp_mjd, tmp_pmx, tmp_pmy

        n = n + 1
        mjd(n) = tmp_mjd
        pmx(n) = tmp_pmx
        pmy(n) = tmp_pmy
    end do

    close(unit)

end subroutine read_eop_pm


!=========================================================
! Save two real arrays into text file
!
! Output format:
!   col1   col2
!=========================================================
subroutine write_two_columns(filename, x, y)

    character(len=*), intent(in)    :: filename
    real(mp), intent(in)            :: x(:), y(:)
    integer                         :: unit, i, n

    n = size(x)
    if (size(y) /= n) error stop "write_two_columns: size mismatch"

    open(newunit=unit, file=filename, status='replace', action='write')

    do i = 1, n
        write(unit,'(2ES24.16)') x(i), y(i)
    end do

    close(unit)

end subroutine write_two_columns


endmodule
