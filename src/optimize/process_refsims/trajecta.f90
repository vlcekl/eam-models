program trajecta
    implicit none
    character*80 :: saux, filname
    integer*8 :: it, ip, i, j, ti, itmin, nt, np
    integer*8 :: nsm, nsp, nfm, nfp, npx
    real*8 :: r, ri, rr, dx, dy, dz, x, y, z, a
    real*8 :: lx, ly, lz, xy, xz, yz
    real*8 :: lxmin, lxmax, lymin, lymax, lzmin, lzmax, rcin, rcout, rc
    real*8 :: rfs_c, rfs_d
    real*8, dimension(:,:), allocatable :: scm, scp, fsm, fsp
    real*8, dimension(:), allocatable :: scmx, fsmx
    real*8, dimension(:), allocatable :: xi, yi, zi

    itmin = 1
    nt = 10001
    !itmin = 0
    !nt = 1
    np = 1024

    a = 3.165
    rcin = a*sqrt(2.0)
    rfs_c = 3.25
    rfs_d = 4.11476
    rcout = a*2.0

    allocate(xi(np), yi(np), zi(np))

    nsm = 7
    nsp = 7
    nfm = 4
    nfp = 4
    allocate(scm(nt, nsm), scp(nt, nsp), fsm(nt, nfm), fsp(nt, nfp))
    scm = 0.0 ; scp = 0.0 ; fsm = 0.0 ; fsp = 0.0
    allocate(scmx(nsm), fsmx(nfm))
    scmx = 0.0 ; fsmx = 0.0

    call getarg(1, filname)
    open(1, file=filname, status='old')
    !open(2, file='data',status='unknown')

    do it = 1, nt
        read(1, *) saux
        read(1, *) saux
        read(1, *) saux
        read(1, *) npx
        read(1, *) saux, saux, saux, saux
        !print *, 'len ', saux
        if (saux == 'pp') then
            read(1, *) lxmin, lxmax
            read(1, *) lymin, lymax
            read(1, *) lzmin, lzmax
            xy = 0.0
            xz = 0.0
            yz = 0.0
        else
            read(1, *) lxmin, lxmax, xy
            read(1, *) lymin, lymax, xz
            read(1, *) lzmin, lzmax, yz
        end if
        read(1, *) saux
        lx = lxmax - lxmin
        ly = lymax - lymin
        lz = lzmax - lzmin
        do ip = 1, npx
            read(1, *) i, ti, x, y, z
            xi(i) = x
            yi(i) = y
            zi(i) = z
        end do

        if (it <= itmin) cycle

        do i = 1, npx
            scmx = 0.0
            fsmx = 0.0
            do j = 1, npx
                if (j == i) cycle

                dx = xi(j) - xi(i)
                dy = yi(j) - yi(i)
                dz = zi(j) - zi(i)

                if (dz > 0.5*lz) then
                    dz = dz - lz
                    dx = dx - xz
                    dy = dy - yz
                else if (dz < -0.5*lz) then
                    dz = dz + lz
                    dx = dx + xz
                    dy = dy + yz
                end if
                if (dy > 0.5*ly) then
                    dy = dy - ly
                    dx = dx - xy
                else if (dy < -0.5*ly) then
                    dy = dy + ly
                    dx = dx + xy
                end if
                if (dx > 0.5*lx) then
                    dx = dx - lx
                else if (dx < -0.5*lx) then
                    dx = dx + lx
                end if

                rr = dx*dx + dy*dy + dz*dz

                if (rr < rcout*rcout) then
                    r = sqrt(rr)
                    ri = 1.0/r
                    scp(it,1) = scp(it,1) + ri**8
                    scp(it,2) = scp(it,2) + ri**9
                    scp(it,3) = scp(it,3) + ri**10
                    scp(it,4) = scp(it,4) + ri**11
                    scp(it,5) = scp(it,5) + ri**12
                    scp(it,6) = scp(it,6) + ri**13
                    scp(it,7) = scp(it,7) + ri**14
                    scmx(1) = scmx(1) + ri**4
                    scmx(2) = scmx(2) + ri**5
                    scmx(3) = scmx(3) + ri**6
                    scmx(4) = scmx(4) + ri**7
                    scmx(5) = scmx(5) + ri**8
                    scmx(6) = scmx(6) + ri**9
                    scmx(7) = scmx(7) + ri**10
                    if (r < rfs_d) then
                        rc = abs(r - rfs_d)
                        fsmx(1) = fsmx(1) + rc
                        fsmx(2) = fsmx(2) + rc*rc
                        fsmx(3) = fsmx(3) + rc*rc*rc
                        fsmx(4) = fsmx(4) + rc*rc*rc*rc
                        if (r < rfs_c) then
                            rc = abs(r - rfs_c)
                            fsp(it,1) = fsp(it,1) + rc*rc
                            fsp(it,2) = fsp(it,2) + rc*rc*r
                            fsp(it,3) = fsp(it,3) + rc*rc*r*r
                            fsp(it,4) = fsp(it,4) + rc*rc*r*r*r
                        end if
                    end if
                end if
            end do
            scm(it,:) = scm(it,:) + sqrt(scmx(:))
            fsm(it,:) = fsm(it,:) + sqrt(fsmx(:))
        end do

        write(*,*) it-1
        write(*,*) fsp(it,:)
        write(*,*) fsm(it,:)
        write(*,*) scp(it,:)
        write(*,*) scm(it,:)

    end do

    !close(2, status='keep')
    close(1, status='keep')

    stop
end program trajecta

