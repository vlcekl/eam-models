C-----------------------------------------------------------------------
C G. Bonny, June 2016
C Implementation of the WRe potential version 20160414
C Stiffened to ZBL function
C All units are eV and Angstrom
C Unpublished so far
C W is EAM2 from Marinica 2013 but modified as for WHHe (Bonny 2014)
C Re is version 20160313
C WRe is version 20160414
C-----------------------------------------------------------------------
C W = 1, Re = 2
C-----------------------------------------------------------------------

      program potGen
      implicit double precision (a-h,o-z)
      
      parameter (ngrid = 6000 )
      
      common/SCALE/dscale(3)
      common/GAUGE/C(3),S(3)

      character*8 ltype
      dimension work(ngrid)
      
      call initialise
C Some tests for the pure species before transformation
      a0BCC = 2.85d0
      a0FCC = 3.6d0
      da = 1.d-1
      eps = 1.d-10
      ! Pure #W
      a = a0BCC; ip = 1; ltype = 'BCC'
      call minimize(eps,a,da,ip,ltype)
      dens0 = eqDens(a,ip,ltype)
      E0 = energy(a,ip,ltype)
      write(*,*)'BCC W:'
      write(*,*)'LATTICE PARAMETER   = ',a
      write(*,*)'EQUILIBRIUM DENSITY = ',dens0
      write(*,*)'ENERGY              = ',E0
      a = a0FCC; ip = 1; ltype = 'FCC'
      call minimize(eps,a,da,ip,ltype)
      dens0 = eqDens(a,ip,ltype)
      E0 = energy(a,ip,ltype)
      write(*,*)'FCC W:'
      write(*,*)'LATTICE PARAMETER   = ',a
      write(*,*)'EQUILIBRIUM DENSITY = ',dens0
      write(*,*)'ENERGY              = ',E0
C Write Spline coefficents for Yuri Osetsky's code
C Minimum cut-off is essential!

      Npoints = 5000

      ! pair
      open(unit=10,file='iapp_WW.dat')
      rmin = 0.01d0
      rmax = 5.4604375d0
      r2min = rmin**2
      r2max = rmax**2
      dr2 = (r2max-r2min)/dble(Npoints)
      do i=1,Npoints+1
        r2 = r2min + dr2*dble(i-1)
        rc = dsqrt(r2)
        work(i) = vpair(rc,1,1)
      enddo
      write(10,*)rmax,r2max,r2min,1.d0/dr2,Npoints+1
      do i=1,Npoints+1
        r2 = r2min + dr2*dble(i-1)
        if (i.ne.1                         ) fm1 = work(i-1)
        if (i.ne.Npoints+1                 ) fp1 = work(i+1)
        if (i.ne.Npoints+1.and.i.ne.Npoints) fp2 = work(i+2)
        if (i.eq.1                         ) fm1 = 2.d0*work(1)-work(2)
        if (i.eq.Npoints+1                 ) fp1 = 0.d0
        if (i.eq.Npoints+1                 ) fp2 = 0.d0
        if (i.eq.Npoints                   ) fp2 = work(Npoints)

        f0  =  work(i)
        call abcd(fm1,f0,fp1,fp2,dr2,a,b,cc,d)
        write(10,*)a,b,cc,d,r2
      enddo
      close(10)

      ! edens
      open(unit=10,file='iapn_W.dat')
      rmin = 0.01d0
      rmax = 5.4604375d0
      r2min = rmin**2
      r2max = rmax**2
      dr2 = (r2max-r2min)/dble(Npoints)
      do i=1,Npoints+1
        r2 = r2min + dr2*dble(i-1)
        rc = dsqrt(r2)
        work(i) = elepot(rc,1)
      enddo
      write(10,*)rmax,r2max,r2min,1.d0/dr2,Npoints+1
      do i=1,Npoints+1
        r2 = r2min + dr2*dble(i-1)
        if (i.ne.1                         ) fm1 = work(i-1)
        if (i.ne.Npoints+1                 ) fp1 = work(i+1)
        if (i.ne.Npoints+1.and.i.ne.Npoints) fp2 = work(i+2)
        if (i.eq.1                         ) fm1 = 2.d0*work(1)-work(2)
        if (i.eq.Npoints+1                 ) fp1 = 0.d0
        if (i.eq.Npoints+1                 ) fp2 = 0.d0
        if (i.eq.Npoints                   ) fp2 = work(Npoints)

        f0  =  work(i)
        call abcd(fm1,f0,fp1,fp2,dr2,a,b,cc,d)
        write(10,*)a,b,cc,d,r2
      enddo
      close(10)

      ! embedding
      open(unit=10,file='iape_W.dat')
      rmin = 0.d0
      rmax = 300.d0
      dr2 = (rmax-rmin)/dble(Npoints)
      do i=1,Npoints+3
        r2 = rmin + dr2*dble(i-1)
        work(i) = demb(r2,1)
      enddo
      write(10,*)rmax,dr2,Npoints+1
      do i=1,Npoints+1
        if (i.ne.1) fm1 = work(i-1)
        if (i.eq.1) fm1 = 2.d0*work(1)-work(2)

        f0  = work(i)
        fp1 = work(i+1)
        fp2 = work(i+2)
        call abcd(fm1,f0,fp1,fp2,dr2,a,b,cc,d)
        write(10,*)a,b,cc,d
      enddo
      close(10)
      write(*,*)'***Tables format Osetsky generated***'

C LAMMPS format
      rmin = 0.d0
      rmax = 5.4604375d0
      rhomin = 0.d0
      rhomax = 10.d0
      !alattW = 3.140d0
      alattW = 3.1880d0
      amW = 183.840d0
      NW = 74
      dr = (rmax-rmin)/dble(Npoints)
      drho = (rhomax-rhomin)/dble(Npoints)
      open(unit=10,file='WW.eam.fs')
      write(10,*)'COMMENT'
      write(10,*)'COMMENT'
      write(10,*)'COMMENT'
      write(10,*)'1 W'
      write(10,*)Npoints,drho,Npoints,dr,rmax
      write(10,*)NW,amW,alattW,' bcc'
      write(10,200)(demb(rhomin+dble(i-1)*drho,1),i=1,Npoints)
      write(10,200)(elepot(dble(i-1)*dr,1),i=1,Npoints)
      write(10,200)(dble(i-1)*dr*vpair(dble(i-1)*dr,1,1),i=1,Npoints)
      close(10)
      write(*,*)'***LAMMPS table generated***'

  100 format(I10,3F19.11)
  200 format(4E20.12)
  300 format(2E20.12)
  101 format('POTENTIAL PAIR ',I2,I2,I7,F8.4,F8.4)
  102 format('POTENTIAL DENS ',   I2,I7,F8.4,F8.4)
  103 format('POTENTIAL EMBED',   I2,I7,F8.4,F10.4)
  104 format('POTENTIAL SDENS',I2,I2,I7,F8.4,F8.4)
  105 format('POTENTIAL SEMBED',  I2,I7,F8.4,F10.4)
      
      end

***************Initialisation of all variables and tables***************
      subroutine initialise
      implicit double precision(a-h,o-z)

      parameter (ngrid = 6000 )

      common/LATTIBCC/rnuBCC(15),nnuBCC(15)
      common/LATTIFCC/rnuFCC(15),nnuFCC(15)
      common/GAUGE/C(3),S(3)
      common/SCALE/dscale(3)
      common /interact/ frho(ngrid,3),rhor(ngrid,3),zr(ngrid,3),
     $  rini,rfin,drar,rhoini,rhofin,drhoar,rcutsq,npoints,npoints1
C Gbonny Feb '10
      common/ZERO/nzero

      data C/0.d0,0.d0,0.d0/
      data S/1.d0,1.d0,1.d0/
      data dscale/1.d0,1.d0,1.d0/

C BCC
      data rnuBCC/1.,1.15470054,1.63299316,1.91485422,2.,2.30940108,
     $ 2.51661148,2.58198890,2.82842712,3.,3.,3.26598632,
     $ 3.41565026,3.46410162,3.46410162/
      data nnuBCC/8,6,12,24,8,6,4*24,8,12,48,6,24/
C FCC
      data rnuFCC/1.,1.41421356,1.73205081,2.,2.23606798,2.44948974,
     $ 2.64575131,2.82842712,3.,3.,3.16227766,3.31662479,
     $ 3.46410161,3.60555128,3.60555128/
      data nnuFCC/12,6,24,12,24,8,48,6,24,12,3*24,48,24/
      
      return
      end
      
***************Energy/lattice parameter optimisation routines***********
      function eqDens(a,ip,ltype)
      implicit double precision (a-h,o-z)
      
      common/LATTIBCC/rnuBCC(15),nnuBCC(15)
      common/LATTIFCC/rnuFCC(15),nnuFCC(15)

      dimension rnu(15),nnu(15)
      character*8 ltype

      if(ltype.eq.'BCC') then
        rnn = dsqrt(3.d0)/2.d0 * a
        do i=1,15; rnu(i) = rnuBCC(i); nnu(i) = nnuBCC(i); enddo
      elseif(ltype.eq.'FCC') then
        rnn = dsqrt(2.d0)/2.d0 * a
        do i=1,15; rnu(i) = rnuFCC(i); nnu(i) = nnuFCC(i); enddo
      else
        stop 'effgauge : Wrong lattice type'
      endif

      nmax = 15

      rho = 0.d0
      do i=1,nmax
        r = rnu(i)*rnn
        rmult = dble(nnu(i))
        rho = rho + elepot(r,ip) * rmult
      enddo
      eqDens = rho
      
      return
      end
      
      function energy(a,ip,ltype)
      implicit double precision (a-h,o-z)

      common/LATTIBCC/rnuBCC(15),nnuBCC(15)
      common/LATTIFCC/rnuFCC(15),nnuFCC(15)

      dimension rnu(15),nnu(15)
      character*8 ltype

      if(ltype.eq.'BCC') then
        rnn = dsqrt(3.d0)/2.d0 * a
        do i=1,15; rnu(i) = rnuBCC(i); nnu(i) = nnuBCC(i); enddo
      elseif(ltype.eq.'FCC') then
        rnn = dsqrt(2.d0)/2.d0 * a
        do i=1,15; rnu(i) = rnuFCC(i); nnu(i) = nnuFCC(i); enddo
      else
        stop 'effgauge : Wrong lattice type'
      endif

      nmax = 15

      rho = 0.d0
      epair = 0.d0
      do i=1,nmax
        r = rnu(i)*rnn
        rmult = dble(nnu(i))
        rho = rho + elepot(r,ip) * rmult
        epair = epair + 0.5d0 * rmult * vpair(r,ip,ip)
      enddo
      energy = epair + demb(rho,ip)

      return
      end
      
      subroutine minimize(eps,a,da,ip,ltype)
      implicit double precision (a-h,o-z)

      character*8 ltype

      maxfun=500   ! Maximum calls to ENERGY
      ic = 0

      step = da
      ea = energy(a,ip,ltype)
  100 continue
      b = a + step
      eb = energy(b,ip,ltype)
      if(eb.gt.ea) then
        step = -step/2.d0
      endif

      a = b
      ea = eb

      ic = ic + 1
      if(ic.gt.maxfun) then
        print*,'Warning : maxfun reached in MINIMIZE '
        return
      endif

      if(step**2.gt.eps**2) goto 100
      a = a + step

      return
      end

***************Gauge transformations************************************
      function vpair(r,ip1,ip2)
      implicit double precision (a-h,o-z)

      common/GAUGE/C(3),S(3)

      call fpair(pair,r,ip1,ip2)
      vpair = pair
      
      return
      end
      
      function elepot(r,ip)
      implicit double precision (a-h,o-z)

      common/GAUGE/C(3),S(3)
      common/SCALE/dscale(3)

      call fdens(dens,r,ip)
      elepot = dens

      return
      end
      
      function demb(rho,ip)
      implicit double precision (a-h,o-z)

      common/GAUGE/C(3),S(3)
      common/SCALE/dscale(3)

      call fembed(F,rho,ip,0)
      demb = F

      return
      end
      
!     subroutine effgauge(a,ip,ltype)
!     implicit double precision (a-h,o-z)
!
!     common/LATTIBCC/rnuBCC(15),nnuBCC(15)
!     common/LATTIFCC/rnuFCC(15),nnuFCC(15)
!     common/GAUGE/C(3),S(3)
!
!     dimension rnu(15),nnu(15)
!     character*8 ltype
!
!     if(ltype.eq.'BCC') then
!       rnn = dsqrt(3.d0)/2.d0 * a
!       do i=1,15; rnu(i) = rnuBCC(i); nnu(i) = nnuBCC(i); enddo
!     elseif(ltype.eq.'FCC') then
!       rnn = dsqrt(2.d0)/2.d0 * a
!       do i=1,15; rnu(i) = rnuFCC(i); nnu(i) = nnuFCC(i); enddo
!     else
!       stop 'effgauge : Wrong lattice type'
!     endif
!
!     nmax = 15
!
!     rho = 0.d0
!     do i=1,nmax
!       r = rnu(i)*rnn
!       rmult = dble(nnu(i))
!       call fdens(dens,r,ip)
!       rho = rho + dens*rmult
!     enddo
!     call fembed(F,rho,ip,1)
!     C(ip) = - F
!     S(ip) = 1.d0/rho
!
!     return
!     end
*************************************** Connection fucntion*************
C----------------------------------------------------------------------
C Weighting function for continuous match up to the 2nd derivative
C between x1 and x2, of fun1 and fun2 : fun = W*fun1 + (1-W)*fun2.
C fun1 holds for x < x1 and fun2 for x > x2
C key = 0,1,2 ---> W, W', W".
C
      function connect(x,x1,x2,key)
      implicit double precision (a-h,o-z)

      parameter(two_3rd=2.d0/3.d0)

      z =(x1+x2-2.d0*x)/(x2-x1)

      if(key.eq.0) then
        if(x .lt. x1)then
          connect= 1.d0
        else if(x .gt. x2)then
          connect= 0.d0
        else
          connect= (0.2d0*z*z*z*z*z- two_3rd*z*z*z+ z)*0.9375d0+ 0.5d0
        end if
      elseif(key.eq.1) then
        if(x .lt. x1 .OR. x .gt. x2) then
          connect= 0.d0
        else
          dz_x = -2.d0/(x2-x1)
          connect=  (z*z*z*z - 2.d0*z*z + 1.d0)*0.9375d0*dz_x
        end if
      elseif(key.eq.2) then
        if(x .lt. x1 .OR. x .gt. x2) then
          connect= 0.d0
        else
          dz_x = -2.d0/(x2-x1)
          connect= (z*z*z - z)*3.75d0*dz_x*dz_x
        end if
      else
        stop '*** WEIGHT: undefined key ***'
      endif

      return
      end
***************Parameterisation of the interactions*********************
      subroutine fpair(epair,r,it1,it2)
      implicit double precision (a-h,o-z)
      
      parameter (ngrid = 6000 )
      
      common /interact/ frho(ngrid,3),rhor(ngrid,3),zr(ngrid,3),
     $  rini,rfin,drar,rhoini,rhofin,drhoar,rcutsq,npoints,npoints1

      if( (it1.eq.1).and.(it2.eq.1) ) then ! W-W EAM2
      rb=0.529177210818181818d0
      r2 = 2.d0
      r1 = 1.d0
      Z = 74.d0
      !print *, 'r', r
      if(r.ge.r2) then ! Cubic spline part
        ! Make effective gauge already here
        ! Now spline is made between ZBL and Veff
        epair = 0.d0
        if(r.ge.5.4604375d0) return
        epair = -5.95486537d2*(2.22-r)**3*H(2.22-r)
     $          +3.42401535d2*(2.46-r)**3*H(2.46-r)
     $          -4.15497479d2*(2.5648975-r)**3*H(2.5648975-r)
     $          +3.71316746d2*(2.629795-r)**3*H(2.629795-r)
     $          -1.07274610d2*(2.6946925-r)**3*H(2.6946925-r)
     $          -2.93282898d0*(2.8663175-r)**3*H(2.8663175-r)
     $          +1.51298572d1*(2.973045-r)**3*H(2.973045-r)
     $          -1.37454674d1*(3.0797725-r)**3*H(3.0797725-r)
     $          +8.16724633d0*(3.5164725-r)**3*H(3.5164725-r)
     $          -3.16850352d0*(3.846445-r)**3*H(3.846445-r)
     $          +8.74607042d-1*(4.1764175-r)**3*H(4.1764175-r)
     $          -1.07878768d0*(4.700845-r)**3*H(4.700845-r)
     $          +1.97939820d0*(4.8953-r)**3*H(4.8953-r)
     $          -8.26149646d-1*(5.089755-r)**3*H(5.089755-r)
     $          -1.62437297d-1*(5.3429525-r)**3*H(5.3429525-r)
     $          -1.45694910d0*(5.401695-r)**3*H(5.401695-r)
     $          +1.42153527d0*(5.4604375-r)**3*H(5.4604375-r)
        !call fdens(dens,r,1)
        !print *, 'epair1', r, epair, epair*r
        !epair = epair! - 2.d0*1.848055990d0*dens/2.232322602d-1
        !print *, 'epair2', r, epair, epair*r
        return
      endif
      if(r.le.r1) then
        epair = 0.d0
        rs = 0.88534d0*rb/((2.d0**0.5d0)*(Z**(1.d0/3.d0)))
        x = r/rs
        qe_sq = 14.3992d0  ! *0.99834058980680 ![eV*A] charge of electron
        if(r.eq.0.d0)return
        epair = Z**2*qe_sq/r * ( 0.1818d0*dexp(-3.2d0*x)
     $                    +0.5099d0*dexp(-0.9423d0*x)
     $                    +0.2802d0*dexp(-0.4029d0*x)
     $                    +0.02817d0*dexp(-0.2016d0*x) )
        return
      endif
      if((r.lt.r2).and.(r.gt.r1))then
        rs = 0.88534d0*rb/((2.d0**0.5d0)*(Z**(1.d0/3.d0)))
        x = r/rs
        qe_sq = 14.3992d0  ! *0.99834058980680 ![eV*A] charge of electron
        fun1 = Z**2*qe_sq/r * ( 0.1818d0*dexp(-3.2d0*x)
     $                    +0.5099d0*dexp(-0.9423d0*x)
     $                    +0.2802d0*dexp(-0.4029d0*x)
     $                    +0.02817d0*dexp(-0.2016d0*x) )
        fun2 =-5.95486537d2*(2.22-r)**3*H(2.22-r)
     $        +3.42401535d2*(2.46-r)**3*H(2.46-r)
     $        -4.15497479d2*(2.5648975-r)**3*H(2.5648975-r)
     $        +3.71316746d2*(2.629795-r)**3*H(2.629795-r)
     $        -1.07274610d2*(2.6946925-r)**3*H(2.6946925-r)
     $        -2.93282898d0*(2.8663175-r)**3*H(2.8663175-r)
     $        +1.51298572d1*(2.973045-r)**3*H(2.973045-r)
     $        -1.37454674d1*(3.0797725-r)**3*H(3.0797725-r)
     $        +8.16724633d0*(3.5164725-r)**3*H(3.5164725-r)
     $        -3.16850352d0*(3.846445-r)**3*H(3.846445-r)
     $        +8.74607042d-1*(4.1764175-r)**3*H(4.1764175-r)
     $        -1.07878768d0*(4.700845-r)**3*H(4.700845-r)
     $        +1.97939820d0*(4.8953-r)**3*H(4.8953-r)
     $        -8.26149646d-1*(5.089755-r)**3*H(5.089755-r)
     $        -1.62437297d-1*(5.3429525-r)**3*H(5.3429525-r)
     $        -1.45694910d0*(5.401695-r)**3*H(5.401695-r)
     $        +1.42153527d0*(5.4604375-r)**3*H(5.4604375-r)
        call fdens(dens,r,1)
        !fun2 = fun2 !- 2.d0*1.848055990d0*dens/2.232322602d-1
        epair = connect(r,r1,r2,0)*fun1 +
     $        (1.d0-connect(r,r1,r2,0))*fun2
        return
      endif
      else
      stop '***SPECIES UNDEFINED***'
      endif
      
      return
      end
      
      subroutine fdens(dens,r,it)
      implicit double precision (a-h,o-z)
      
      parameter (ngrid = 6000 )
      
      common /interact/ frho(ngrid,3),rhor(ngrid,3),zr(ngrid,3),
     $  rini,rfin,drar,rhoini,rhofin,drhoar,rcutsq,npoints,npoints1
     
      dens = 0.d0
      if(it.eq.1) then ! W EAM2
        ! Becomes constant as described in Marinica 2013
        ! Thus here no modification to stiffen necessary
        ! Already make effective gauge here: rescale with S
        if(r.ge.4.9d0) return
        if(r.le.2.0) then
          rc = 2.0
          dens = +1.29203135d+2*(2.5-rc)**3*H(2.5-rc)
     $           -9.25284676d+0*(2.8663175-rc)**3*H(2.8663175-rc)
     $           +2.63967417d+0*(3.1-rc)**3*H(3.1-rc)
     $           +2.06410900d+0*(3.5-rc)**3*H(3.5-rc)
     $           +6.03258240d-2*(4.9-rc)**3*H(4.9-rc)
          return
        endif
        dens = +1.29203135d+2*(2.5-r)**3*H(2.5-r)
     $         -9.25284676d+0*(2.8663175-r)**3*H(2.8663175-r)
     $         +2.63967417d+0*(3.1-r)**3*H(3.1-r)
     $         +2.06410900d+0*(3.5-r)**3*H(3.5-r)
     $         +6.03258240d-2*(4.9-r)**3*H(4.9-r)
        return
      else
      stop '***SPECIES UNDEFINED***'
      endif
      
      return
      end
      
      subroutine fembed(embed,rho,it,id)
      implicit double precision (a-h,o-z)
      
      parameter (ngrid = 6000 )
      
      common /interact/ frho(ngrid,3),rhor(ngrid,3),zr(ngrid,3),
     $  rini,rfin,drar,rhoini,rhofin,drhoar,rcutsq,npoints,npoints1
C Gbonny Feb '10
      common/ZERO/nzero

      if(id.ne.0) stop 'fembed: derivative not implemented'

      if(it.eq.1) then ! W EAM2
        ! transform into effective gauge already here
        !ro = rho/2.232322602d-1
        ro = rho
        if(rho.le.1.359141225d0) then
          embed = -3.27882058d0*ro**(0.5d0)
     $       -9.7768186d-4*ro**2 ! + ro*1.848055990d0
!          embed = -4.18353755d+00*ro**(0.5d0)
!     $       -9.63668936d-03*ro**2 ! + ro*1.848055990d0
        else
          embed = -3.27882058d0*rho**(0.5d0)
     $       -9.7768186d-4*rho**2 ! + ro*1.848055990d0
!          embed = -4.18353755d+00*ro**(0.5d0)
!     $       -9.63668936d-03*ro**2 ! + ro*1.848055990d0
!          embed = -5.524855802d00 + 2.317313103d-1 * rho
!     $            -3.665345949d-2 * rho**2 + 8.989367404d-3 * rho**3
        endif
        return
      else
        stop '***SPECIES UNDEFINED***'
      endif

      return
      end
      
      function aMorse(r,r0,alpha)  ! Morse function
      implicit double precision (a-h,o-z)
      
      aMorse = dexp(-2.d0*alpha*(r-r0)) - 2.d0*dexp(-alpha*(r-r0))
      
      return
      end
      
      function Psi(x)  ! Cut off function
      implicit double precision (a-h,o-z)
      
      if(x.ge.0.d0) then
        Psi = 0.d0
      else
        Psi = x**4/(1.d0 + x**4)
      endif
      
      return
      end
      
      function H(x) ! Step function
      implicit double precision (a-h,o-z)
      
      if (x.gt.0.d0) then
        H = 1.d0
      else
        H = 0.D0
      endif
      
      return
      end
      
      subroutine abcd ( fm1,f,fp1,fp2,dx,a,b,c,d )
      implicit double precision (a-h,o-z)

      a  =  f
      b  = (fp1-fm1)/2.d0
      b1 = (fp2-a)/2.d0
      d  = (b1 + b + 2.d0*(a-fp1))
      c  = (b1-b-3.d0*d)/2.d0
      
      b = b/dx
      c = c/dx**2
      d = d/dx**3

      return
      end
