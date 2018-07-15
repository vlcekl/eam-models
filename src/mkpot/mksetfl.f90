!
!  File name: combine.f90
!  Date:      2010/10/04 00:09
!  Author:    Lukas Vlcek
! 

program mksetfl
    implicit none
    real*4 :: etime, elapsed(2)
    character*80 :: filname
    character*12 :: fftype, alat
    integer*8 :: i, nel, m, n, nrho, nr, na, j
    real*8 :: fac, a, c, e, drho, dr, cutoff, lc, ma, rho, r, rmin
    real*8 :: d, c0, c1, c2
    character*2, dimension(:), allocatable :: el
    real*8, dimension(:), allocatable :: fembd, edens, fpair

    fac = 1.60217662e-19*6.022140857e23/1000.0 ! eV to kJ/mol
    rmin = 1.8

    call getarg(1, filname)

    ! read parameters
    open(1, file=filname, status='old')

    read(1, *) nel
    allocate(el(nel))
    read(1, *) (el(i), i = 1, nel)
    read(1, *) nrho, drho, nr, dr, cutoff
    read(1, *) na, ma, lc, alat

    read(1, *) fftype
    if (fftype .eq. 'SC') then
        read(1, *) a, c, e, m, n
    else if (fftype .eq. 'FS') then
        d = 4.400224
        c = 3.25
        cutoff = d
        read(1, *) a, c0, c1, c2

!    A = 1.896373
!    c0 = 47.1346499
!    c1 = -33.7665655
!    c2 = 6.2541999

    else
        print *, 'Uknown potential ', fftype
        stop
    end if

    close(1, status='keep')

    !dr = cutoff/float(nr-1)

    ! generate potential functions

    allocate(fembd(nrho), edens(nr), fpair(nr))

    do i = 1, nrho
        rho = real(i-1)*drho
        fembd(i) = -a*sqrt(rho)
    end do

    do i = 1, nr
        r = real(i-1)*dr
        !if (r < rmin) then
        !    r = rmin
        !end if
        if (r <= d) then
            edens(i) = (r-d)**2
        else
            edens(i) = 0.0
        end if
    end do

    do i = 1, nr
        r = real(i-1)*dr
        !if (r < rmin) then
        !    r = rmin
        !end if
        if (r <= c) then
            fpair(i) = r*((r-c)**2*(c0 + c1*r + c2*r**2))
        else
            fpair(i) = 0.0
        end if
    end do

    ! write the potential file
    write(*, *) 'Comment 1'
    write(*, *) 'Comment 1'
    write(*, *) 'Comment 1'

    write(*, '(I5,1X,3A2)') nel, (el(i), i = 1, nel)

    write(*, '(I5,E24.16,I5,E24.16,E24.16)') nrho, drho, nr, dr, cutoff

    write(*, '(I5,1X,F14.4,1X,F14.4,1X,A10)') na, ma, lc, alat

    do i = 1, nrho/5
        write(*,'(5E24.16)') (fembd((i-1)*5+j), j = 1,5)
    end do
    do i = 1, nr/5
        write(*,'(5E24.16)') (edens((i-1)*5+j), j = 1,5)
    end do
    do i = 1, nr/5
        write(*,'(5E24.16)') (fpair((i-1)*5+j), j = 1,5)
    end do

    deallocate(el, fembd, edens, fpair)

    stop
end program mksetfl

!  end of combine.f90
