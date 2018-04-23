!
!  File name: combine.f90
!  Date:      2010/10/04 00:09
!  Author:    Lukas Vlcek
! 

program combine
    implicit none
    real*4 :: etime, elapsed(2)
    character*80 :: filname
    character*12 :: sname
    integer*8 :: i, nmax, ncom, ip, npa, np0, np1, np2, iaux, it, np
    integer*8 :: ia, i0, i1, i2, nfp, nfm, nsp, nsm
    real*8 :: amin, amax, da, fac, fni, bmin
    real*8 :: temp, beta, finmax, shee, kt, uave
    real*8, dimension(:), allocatable :: aaa, cc0, cc1, cc2
    real*8, dimension(:), allocatable :: ge, he, se, db
    real*8, dimension(:), allocatable :: hee, eee, uuu
    real*8, dimension(:,:), allocatable :: scm, scp, fsm, fsp

    fac = 1.60217662e-19*6.022140857e23/1000.0 ! eV to kJ/mol
    np = 128
    np = 111
    np = 10
    fni = 1.0/float(np)
    nmax = 1001

    allocate(aaa(100), cc0(100), cc1(100), cc2(100))

    call getarg(1, filname)

    open(1, file=filname, status='old')

    read(1, *) sname, npa
    read(1, *) amin, amax, da
    do i=1,npa
        aaa(i) = amin + float(i-1)*da
    end do
    aaa = aaa*fac

    read(1, *) sname, np0
    read(1, *) amin, amax, da
    do i=1,np0
        cc0(i) = amin + float(i-1)*da
    end do

    read(1, *) sname, np1
    read(1, *) amin, amax, da
    do i=1,np1
        cc1(i) = amin + float(i-1)*da
    end do

    read(1, *) sname, np2
    read(1, *) amin, amax, da
    do i=1,np2
        cc2(i) = amin + float(i-1)*da
    end do

    close(1, status='keep')

    ncom = npa*np0*np1*np2
    print *, ncom
    allocate(ge(ncom), he(ncom), se(ncom), db(ncom))
    ge = 0.0 ; he = 0.0 ; se = 0.0 ; db = 0.0

    ! target energies from simulation
    allocate(hee(nmax))
    allocate(eee(nmax),uuu(nmax))
    finmax = 1.0/real(nmax)
    call getarg(2, filname)
    open(1, file=filname, status='old')
    read(1,*) temp
    do it = 1, nmax
        read(1,*) iaux, hee(it)
    end do
    close(1, status='keep')

    !temp = 300.0
    beta = 1.0/(8.314472*temp/1000.0) ! (kJ/mol)^-1
    kt = 1./beta

    hee = hee*fac 
    shee = sum(hee)*finmax

    ! reference cfg properties
    nsm = 7 ; nsp = 7 ; nfm = 4 ; nfp = 4
    allocate(scm(nmax, nsm), scp(nmax, nsp), fsm(nmax, nfm), fsp(nmax, nfp))
    scm = 0.0 ; scp = 0.0 ; fsm = 0.0 ; fsp = 0.0

    call getarg(3, filname)
    open(1, file=filname, status='old')
    do it = 1, nmax
        read(1,*)iaux
        read(1,*)fsp(it,:)
        read(1,*)fsm(it,:)
        read(1,*)scp(it,:)
        read(1,*)scm(it,:)
    end do
    close(1, status='keep')

    bmin = 10000.0
    ip = 0
    do ia = 1, npa
        do i0 = 1, np0
            do i1 = 1, np1
                do i2 = 1, np2
                    ip = ip + 1

                    uuu(:) = beta*(cc0(i0)*fsp(:,1) + cc1(i1)*fsp(:,2) + cc2(i2)*fsp(:,3)- aaa(ia)*fsm(:,1) - hee(:))
    
                    uave = sum(uuu)*finmax
                    uuu = uuu - uave
    
                    eee = exp(-uuu)
                    ge(ip) = sum(eee)
    
                    he(ip) = (sum(eee*(beta*hee + uuu + uave))/ge(ip) - beta*shee)*kt
                    he(ip) = he(ip)*fni

                    ge(ip) = -log(ge(ip)*finmax)
    
                    eee = exp(-0.5*(uuu - ge(ip)))
                    !eee = exp(-0.5*(uuu))
                    db(ip) = -2.0*log(sum(eee)*finmax)
                    
                    ge(ip) = (ge(ip) + uave)*kt * fni
                    se(ip) = he(ip) - ge(ip)

                    if (db(ip) < bmin) then
                        bmin = db(ip)
                        print *,aaa(ia)/fac,cc0(i0),cc1(i1),cc2(i2),ge(ip),he(ip),se(ip),db(ip)
                    end if
                end do
            end do
        end do
    end do

    ! print results
    open(1, file='optout', status='unknown')
    ip = 0
    do ia = 1, npa
        do i0 = 1, np0
            do i1 = 1, np1
                do i2 = 1, np2
                    ip = ip + 1
                    write(1,10)aaa(ia)/fac,cc0(i0),cc1(i1),cc2(i2),ge(ip),he(ip),se(ip),db(ip)
                end do
            end do
        end do
    end do
    close(1, status='keep')

    stop
10  format(4(f8.2,1x),4(f14.6,1x))
end program combine

!  end of combine.f90
