program trajecta
    implicit none
    character*80 :: saux, filname
    integer*8 :: it, ip, i, j, ti, l, lmax, mi, ns, nc, iaux, ii
    integer*8 :: itmin, nmol, im, jm, is, js, np, nt, io, ih1, ih2 
    real*8 :: r, rr, r6, dx, dy, dz, qi, x, y, z, dr, dri, rcut, erf, ferf, M_PI
    real*8 :: fx, f1, f2, f3, aaa, ccc, qqq, s6, s12, sqq, utot, fac, faux, pq
    real*8 :: siref, epref, qref, uq, pwa, pwc, pi, ulj, ucoul, press
    real*8 :: sumq, ri, lx, ly, lz, xlo, xhi, ylo, yhi, zlo, zhi, vol, v2p, viriw
    integer*8, dimension(:,:,:), allocatable :: hhh
    integer*8, dimension(2) :: ntyp
    real*8, dimension(4) :: qr
    real*8, dimension(:), allocatable :: xi, yi, zi, xc, yc, zc
    real*8, dimension(:), allocatable :: ut, uw, uqs, uql, pt, tt
    real*8, dimension(:,:), allocatable :: xs, ys, zs

    M_PI = 3.14159265359
    fac = (1.602176487d-19)**2*6.02214179d23*1d10/(4.*3.14159265359*8.854187817d-12)/1000.0/4.184

    !print *, 'fac', fac

    qr(1) = 0.0
    qr(2) = 0.5
    qr(3) = 0.5
    qr(4) = -1.0

!    aaa = 4.0*0.14980054611663479923*3.21591921335**12
!    ccc = 4.0*0.14980054611663479923*3.21591921335**6
!    qqq = 1.07554
!    qqq = 0.0
!    siref = 3.21591921335
!    epref = 0.14980054611663479923
!    qref = 1.07554
!    siref = 3.1790081
!    epref = 0.15497033221797323135 
!    qref =  1.08181752
    siref = 3.12217297613
    epref = 0.19639051593331739961
    qref = 1.082914
    aaa = 4.0*epref*siref**12
    ccc = 4.0*epref*siref**6
    qqq = qref

    itmin = 10000
    nt = 200000
    !itmin = 0
    !nt = 10
    !nt = 10310*2

    dr = 0.05
    dri = 1.0/dr

!    lx = 2.0*9.852275
!    ly = 2.0*9.852275
!    lz = 2.0*9.852275
!    rcut = 0.5*lx
!    lx = 2.0*8.510424
    lx = 2.0*7.389207
    !lx = 2.0*9.309817
    ly = lx
    lz = lx
    !rcut = 8.5
    rcut = 0.5*lx
    erf = 78.0
    ferf = (erf - 1.0)/((2.0*erf + 1.0)*rcut**3)

    ntyp(2) = 108 !256
    !ntyp(2) = 216 !256
    ntyp(1) = 2*ntyp(2) !512
    nmol = ntyp(2)
    np = sum(ntyp)
    ns = 4

    allocate(xi(np), yi(np), zi(np), xc(nmol), yc(nmol), zc(nmol))
    allocate(xs(nmol, ns), ys(nmol, ns), zs(nmol, ns))

    ii = nt - itmin
    allocate(ut(ii), uw(ii), uqs(ii), uql(ii), pt(ii), tt(ii))
    ut = 0.0 ; uw = 0.0 ; uqs = 0.0 ; uql = 0.0 ; pt = 0.0 ; tt = 0.0

    open(1, file='enerun', status='old')
    do it = 1, nt
        if (it <= itmin) then
            read(1, *) saux
        else 
            ii = it - itmin
            read(1, *) iaux, faux, ut(ii), faux, faux, uw(ii), uqs(ii), uql(ii), faux, faux, faux, pt(ii), tt(ii)
        end if
    end do

    close(1, status='keep')

    lmax = int(rcut*dri) + 1 + 1
    allocate (hhh(3, nt-itmin, lmax))
    hhh = 0

    open(1, file='tip4p.dump', status='old')
    open(2, file='data',status='unknown')
    write(2,*) '#refpar', siref, epref, qref, 1.0

    do it = 1, nt
        if (mod(it,100) == 0) print *, it
        read(1, *) saux
        read(1, *) saux

        read(1, *) saux
        read(1, *) saux

        read(1, *) saux
        read(1, *) xlo, xhi
        read(1, *) ylo, yhi
        read(1, *) zlo, zhi
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        vol = lx*ly*lz

        read(1, *) saux

        do ip = 1, np
            read(1, *) i, mi, ti, qi, x, y, z
            xi(i) = x
            yi(i) = y
            zi(i) = z
        end do

        do im = 1, nmol
            io  = 3*(im-1) + 1
            ih1 = 3*(im-1) + 2
            ih2 = 3*(im-1) + 3

            dx = xi(ih1) - xi(io)
            dy = yi(ih1) - yi(io)
            dz = zi(ih1) - zi(io)

            if (dx > 0.5*lx) then
                xi(ih1) = xi(ih1) - lx
            else if (dx < -0.5*lx) then
                xi(ih1) = xi(ih1) + lx
            end if
            if (dy > 0.5*ly) then
                yi(ih1) = yi(ih1) - ly
            else if (dy < -0.5*ly) then
                yi(ih1) = yi(ih1) + ly
            end if
            if (dz > 0.5*lz) then
                zi(ih1) = zi(ih1) - lz
            else if (dz < -0.5*lz) then
                zi(ih1) = zi(ih1) + lz
            end if

            dx = xi(ih2) - xi(io)
            dy = yi(ih2) - yi(io)
            dz = zi(ih2) - zi(io)

            if (dx > 0.5*lx) then
                xi(ih2) = xi(ih2) - lx
            else if (dx < -0.5*lx) then
                xi(ih2) = xi(ih2) + lx
            end if
            if (dy > 0.5*ly) then
                yi(ih2) = yi(ih2) - ly
            else if (dy < -0.5*ly) then
                yi(ih2) = yi(ih2) + ly
            end if
            if (dz > 0.5*lz) then
                zi(ih2) = zi(ih2) - lz
            else if (dz < -0.5*lz) then
                zi(ih2) = zi(ih2) + lz
            end if

            !xc(im) = (15.9994*xi(io) + 1.0079*(xi(ih1) + xi(ih2)))/(15.9994+2.0*1.0079)
            !yc(im) = (15.9994*yi(io) + 1.0079*(yi(ih1) + yi(ih2)))/(15.9994+2.0*1.0079)
            !zc(im) = (15.9994*zi(io) + 1.0079*(zi(ih1) + zi(ih2)))/(15.9994+2.0*1.0079)
            xc(im) = xi(io)
            yc(im) = yi(io)
            zc(im) = zi(io)

!            if (xc(im) > 0.5*lx) then
!                xc(im) = xc(im) - lx
!                xi(io) = xi(io) - lx
!                xi(ih1) = xi(ih1) - lx
!                xi(ih2) = xi(ih2) - lx
!            else if (xc(im) < 0.5*lx) then
!                xc(im) = xc(im) + lx
!                xi(io) = xi(io) + lx
!                xi(ih1) = xi(ih1) + lx
!                xi(ih2) = xi(ih2) + lx
!            end if
!            if (yc(im) > 0.5*ly) then
!                yc(im) = yc(im) - ly
!                yi(io) = yi(io) - ly
!                yi(ih1) = yi(ih1) - ly
!                yi(ih2) = yi(ih2) - ly
!            else if (yc(im) < 0.5*ly) then
!                yc(im) = yc(im) + ly
!                yi(io) = yi(io) + ly
!                yi(ih1) = yi(ih1) + ly
!                yi(ih2) = yi(ih2) + ly
!            end if
!            if (zc(im) > 0.5*lz) then
!                zc(im) = zc(im) - lz
!                zi(io) = zi(io) - lz
!                zi(ih1) = zi(ih1) - lz
!                zi(ih2) = zi(ih2) - lz
!            else if (zc(im) < 0.5*lz) then
!                zc(im) = zc(im) + lz
!                zi(io) = zi(io) + lz
!                zi(ih1) = zi(ih1) + lz
!                zi(ih2) = zi(ih2) + lz
!            end if

            xs(im, 1) = xi(io)
            ys(im, 1) = yi(io)
            zs(im, 1) = zi(io)
            xs(im, 2) = xi(ih1)
            ys(im, 2) = yi(ih1)
            zs(im, 2) = zi(ih1)
            xs(im, 3) = xi(ih2)
            ys(im, 3) = yi(ih2)
            zs(im, 3) = zi(ih2)
            r =     (xi(ih1) + xi(ih2) - 2.0*xi(io))**2
            r = r + (yi(ih1) + yi(ih2) - 2.0*yi(io))**2
            r = r + (zi(ih1) + zi(ih2) - 2.0*zi(io))**2
            ri = 1.0/sqrt(r)
            xs(im, 4) = (xi(ih1) + xi(ih2) - 2.0*xi(io))*ri*0.1546 + xi(io)
            ys(im, 4) = (yi(ih1) + yi(ih2) - 2.0*yi(io))*ri*0.1546 + yi(io)
            zs(im, 4) = (zi(ih1) + zi(ih2) - 2.0*zi(io))*ri*0.1546 + zi(io)
            !print *, 'C', xc(im), yc(im), zc(im)
            !do i = 1, 4
            !    print *, i, xs(im,i), ys(im,i), zs(im,i)
            !end do
        end do

        if (it <= itmin) cycle
        ii = it - itmin

        s6 = 0.0
        s12 = 0.0
        sqq = 0.0
        sumq = 0.0
        nc = 0
        do im = 1, nmol-1
            !print *, ((xs(im,1)-xs(im,4))**2+(ys(im,1)-ys(im,4))**2+(zs(im,1)-zs(im,4))**2)**0.5
            do jm = im+1, nmol
                dx = xc(jm) - xc(im)
                dy = yc(jm) - yc(im)
                dz = zc(jm) - zc(im)

                if (dx > 0.5*lx) then
                    dx = dx - lx
                else if (dx < -0.5*lx) then
                    dx = dx + lx
                end if
                if (dy > 0.5*ly) then
                    dy = dy - ly
                else if (dy < -0.5*ly) then
                    dy = dy + ly
                end if
                if (dz > 0.5*lz) then
                    dz = dz - lz
                else if (dz < -0.5*lz) then
                    dz = dz + lz
                end if
                rr = dx*dx + dy*dy + dz*dz
                if (rr > rcut*rcut) then
                    hhh(1,ii,lmax) = hhh(1,ii,lmax) + 1
                    hhh(2,ii,lmax) = hhh(2,ii,lmax) + 4
                    hhh(3,ii,lmax) = hhh(3,ii,lmax) + 4
                    cycle
                end if

                do is = 1, ns
                    do js = 1, ns
                        dx = xs(jm, js) - xs(im, is)
                        dy = ys(jm, js) - ys(im, is)
                        dz = zs(jm, js) - zs(im, is)
                        
                        if (dz > 0.5*lz) then
                            dz = dz - lz
                        else if (dz < -0.5*lz) then
                            dz = dz + lz
                        end if
                        if (dy > 0.5*ly) then
                            dy = dy - ly
                        else if (dy < -0.5*ly) then
                            dy = dy + ly
                        end if
                        if (dx > 0.5*lx) then
                            dx = dx - lx
                        else if (dx < -0.5*lx) then
                            dx = dx + lx
                        end if
                        
                        rr = dx*dx + dy*dy + dz*dz
                        r = sqrt(rr)

                        if (is == 1 .and. js == 1) then
                            r6 = 1.0/rr**3
                            s6 = s6 + r6
                            s12 = s12 + r6*r6
                        else if (is > 1 .and. js > 1) then
                            !sqq = sqq + qr(is)*qr(js)*(1./r + rr*ferf)
                            !sqq = sqq + (1./r + rr*ferf)
                            sqq = sqq + qr(is)*qr(js)*(1./r + rr*ferf)
                            if (im == 10 .and. jm == 20) then
                                !print *, it, is, js, r, qr(is)*qr(js)*(1./r + rr*ferf), sqq
                            end if
                            sumq = sumq + qr(is)*qr(js)
                            !print *, im, jm, is, js, qr(is), qr(js), qr(is)*qr(js)/r, qr(is)*qr(js)*(1.0/r-rr*ferf),r
                        end if

                        l = int(r*dri) + 1
                        if (l > lmax) then
                            if (is == 1 .and. js == 1) then
                                hhh(1,ii,lmax) = hhh(1,ii,lmax) + 1
                            else if ((is==2.or.is==3).and.(js==2.or.js==3)) then
                                hhh(3,ii,lmax) = hhh(3,ii,lmax) + 1
                            else if (is < 4 .and. js < 4) then
                                hhh(2,ii,lmax) = hhh(2,ii,lmax) + 1
                            end if
                        else
                            if (is == 1 .and. js == 1) then
                                hhh(1,ii,l) = hhh(1,ii,l) + 1
                            else if ((is==2.or.is==3).and.(js==2.or.js==3)) then
                                hhh(3,ii,l) = hhh(3,ii,l) + 1
                            else if (is < 4 .and. js < 4) then
                                hhh(2,ii,l) = hhh(2,ii,l) + 1
                            end if
                        end if
                    end do
                end do
                nc = nc + 1
                !print *, 'sqq', nc, fac*sqq*qqq*qqq
            end do
        end do

        utot = aaa*s12 - ccc*s6  + qqq*qqq*fac*sqq
        !utot = qqq*qqq*fac*sqq
        ulj = aaa*s12 - ccc*s6 !  + qqq*qqq*fac*sqq
        v2p = 1000.0/6.02214179e23/(1e-30)/101325 * 4.184
        !viriw = (12*aaa*s12 - 6*ccc*s6)*v2p/(3.0*vol) + 4.559390401983079*298.15
        !viriw = (12*aaa*s12 - 6*ccc*s6)*v2p/(3.0*vol) + 4.559390401983079*283.03617
        viriw = (12*aaa*s12 - 6*ccc*s6)*v2p/(3.0*vol) + 4.8372*tt(ii) ! i.g. + vdw
        pwa = (12*s12)*v2p/(3.0*vol)
        pwc = (6*s6)*v2p/(3.0*vol)
        pi = 4.8372*tt(ii)  ! ideal gas
        pq = (pt(ii) - (pwa*aaa - pwc*ccc + pi))/qqq**2
        !pq = pt(ii) - viriw
        ucoul = uqs(ii)+uql(ii)
        uq = ucoul/qqq**2
        press = pi + pwa*aaa - pwc*ccc + pq*qqq**2
        !print *, 'pp', pt(ii), press, pi, pwa, aaa, pwc, ccc, pq, qqq**2

        write(2,*) '#', ii
        write(2,*) ut(ii), ulj, ucoul, s12, s6, sqq*fac, uq, pt(ii), pwa, pwc, pi, pq, tt(ii)
        do l = 1, lmax
            write(2,'(3(1X,I6))') hhh(1,ii,l), hhh(2,ii,l), hhh(3,ii,l)
        end do
    end do

    close(2, status='keep')
    close(1, status='keep')

    open(1, file='href',status='unknown')
    do l = 1, lmax
        fx = 1.0/real(nt-itmin)
        f1 = fx*float(sum(hhh(1,1:nt-itmin,l)))
        f2 = fx*float(sum(hhh(2,1:nt-itmin,l)))
        f3 = fx*float(sum(hhh(3,1:nt-itmin,l)))
        write(1,*) l*dr,f1,f2,f3
    end do
    close(1,status='keep')

    open(1, file='histog',status='unknown')
    do l = 1, lmax
        r = real(l-1)*dr
        fx = 1.0/(4.0/3.0*M_PI*((r+0.5*dr)**3 - (r-0.5*dr)**3))
        fx = fx/real(nt-itmin)
        fx = lx*ly*lz*fx
        f1 = fx*float(sum(hhh(1,1:nt-itmin,l)))/float((ntyp(2)-1)*ntyp(2)/2)
        f2 = fx*float(sum(hhh(2,1:nt-itmin,l)))/float(ntyp(1)*ntyp(2))
        f3 = fx*float(sum(hhh(3,1:nt-itmin,l)))/float((ntyp(1)-1)*ntyp(1)/2-ntyp(2))
        write(1,*) l*dr,f1,f2,f3
    end do
    close(1,status='keep')

    stop
end program trajecta

