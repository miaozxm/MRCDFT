!==============================================================================!
! MODULE EM                                                                    !
!                                                                              !
! This module calculates the electromagnetic multipole operators               !
!                                                                              !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!

Module EM
    use Constants,only: r64
    implicit none
    
contains

    subroutine calculate_Qlm(l,iphi,it,Qlm,pQlm,cQlm,pcQlm)
        !-------------------------------------------------------------------------------
        !    calculate <q1| Q_lm R|q2>/<q1|R|q2>  = \sum_{m1 m2}<m1|Q_lm|m2> rho_{m2m1}
        !    calculate <q1| Q^\dagger_lm R|q2>/<q1|R|q2>  = \sum_{m1 m2}<m1|Q_lm|m2> rho_{m1m2}
        !  where R = R(alpha,beta,gamma, phi_n,phi_p) 
        !    = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
        !--------------------------------------------------------------------------------
        use Constants,only: itx
        use Globals, only: BS, mix
        integer :: l,iphi,it
        integer :: m,ifg,ndsp,i0sp,m1,nr1,nl1,nj1,nm1,m2,nr2,nl2,nj2,nm2
        complex(r64), dimension(-l:l,itx):: Qlm,pQlm,cQlm,pcQlm
        real(r64) :: Qlm_m1m2
        do m = -l,l
            Qlm(m,it) = (0.d0,0.d0)
            pQlm(m,it) = (0.d0,0.d0)
            cQlm(m,it) = (0.d0,0.d0)
            pcQlm(m,it) = (0.d0,0.d0)
            do ifg = 1,2
                ndsp = BS%HO_sph%idsp(1,ifg)
                i0sp = BS%HO_sph%iasp(1,ifg)
                do m1 = 1,ndsp
                    nr1= BS%HO_sph%nljm(i0sp+m1,1) ! n_r
                    nl1= BS%HO_sph%nljm(i0sp+m1,2) ! l
                    nj1= BS%HO_sph%nljm(i0sp+m1,3) ! j +1/2
                    nm1= BS%HO_sph%nljm(i0sp+m1,4) ! m_j + 1/2
                    do m2 = 1,ndsp
                        nr2= BS%HO_sph%nljm(i0sp+m2,1) ! n_r
                        nl2= BS%HO_sph%nljm(i0sp+m2,2) ! l
                        nj2= BS%HO_sph%nljm(i0sp+m2,3) ! j +1/2
                        nm2= BS%HO_sph%nljm(i0sp+m2,4) ! m_j + 1/2
                        ! <m1|Q_lm|m2>
                        call multipole_matrix_elements(nr1,nl1,nj1,nm1,l,m,nr2,nl2,nj2,nm2,Qlm_m1m2)
                        ! <m1|Q_lm|m2> rho_{m2m1}
                        Qlm(m,it) = Qlm(m,it) + Qlm_m1m2*mix%rho_mm(m2,m1,ifg,iphi,it)
                        pQlm(m,it) = pQlm(m,it) + Qlm_m1m2*mix%prho_mm(m2,m1,ifg,iphi,it)
                        ! <m1|Q_lm|m2> rho_{m1m2}
                        cQlm(m,it) = cQlm(m,it) + Qlm_m1m2*mix%rho_mm(m1,m2,ifg,iphi,it)
                        pcQlm(m,it) = pcQlm(m,it) + Qlm_m1m2*mix%prho_mm(m1,m2,ifg,iphi,it)
                    end do
                end do 
            end do
        end do
    end subroutine

    subroutine calculate_r2(iphi,it,r2,pr2)
        !--------------------------------------------------
        !  calculate <q1| r^2 R|q2>/<q1|R|q2>
        !  = \sum_{m1 m2} <m1|r^2|m2> rho_{m2m1}.
        !  where <m1|r^2|m2> = <n1 l1|r^2|n2 l2> when l1==l2 and j1==j2 and mj1==mj2
        !  and R = R(alpha,beta,gamma, phi_n,phi_p) = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
        !--------------------------------------------------
        use Constants,only: itx
        use Globals, only: BS, mix
        integer :: l,iphi,it
        integer :: ifg,ndsp,i0sp,m1,nr1,nl1,nj1,nm1,m2,nr2,nl2,nj2,nm2
        complex(r64) :: r2,pr2
        r2 = 0.0d0
        pr2 = 0.0d0
        do ifg = 1,2
            ndsp = BS%HO_sph%idsp(1,ifg)
            i0sp = BS%HO_sph%iasp(1,ifg)
            do m1 = 1,ndsp
                nr1= BS%HO_sph%nljm(i0sp+m1,1) ! n_r
                nl1= BS%HO_sph%nljm(i0sp+m1,2) ! l
                nj1= BS%HO_sph%nljm(i0sp+m1,3) ! j +1/2
                nm1= BS%HO_sph%nljm(i0sp+m1,4) ! m_j + 1/2
                do m2 = 1,ndsp
                    nr2= BS%HO_sph%nljm(i0sp+m2,1) ! n_r
                    nl2= BS%HO_sph%nljm(i0sp+m2,2) ! l
                    nj2= BS%HO_sph%nljm(i0sp+m2,3) ! j +1/2
                    nm2= BS%HO_sph%nljm(i0sp+m2,4) ! m_j + 1/2
                    if(nl1==nl2 .and. nj1==nj2 .and. nm1==nm2) then 
                        r2 = r2 + rl_nl(nr1,nl1,2,nr2,nl2)*mix%rho_mm(m2,m1,ifg,iphi,it)
                        pr2 = pr2 + rl_nl(nr1,nl1,2,nr2,nl2)*mix%prho_mm(m2,m1,ifg,iphi,it)
                    end if 
                end do
            end do 
        end do
    end subroutine

    subroutine multipole_matrix_elements(n1,l1,j1,m1,l,m,n2,l2,j2,m2,mpm) 
        !-----------------------------------------------------------------------
        !      Calculates <n1,l1,J1,M1|Q_(lm)|n2,l2,J2,M2> 
        !                =<n1,l1| r^l |n2,l2>*<l1,J1,M1| Y_(lm) |l2,J2,M2> 
        !      where, Q_(lm)= r^l Y_(lm), x=r/b 
        !      J1 = j1 -1/2, M1 = m1 - 1/2 and also for J2,M2 (half-intervalues)
        !-----------------------------------------------------------------------
        integer :: n1,l1,j1,m1,l,m,n2,l2,j2,m2
        real(r64) :: mpm
        mpm=rl_nl(n1,l1,l,n2,l2)*ylm_ljm(l1,j1,m1,l,m,l2,j2,m2) 
    end subroutine
    
    double precision function rl_nl(n1,l1,l,n2,l2)
        !-----------------------------------------------------------------------
        !               calculate <n1,l1| r^l |n2,l2>
        !
        !     calculates the radial functions for the spherical oscillator
        !     the wave function R_nl(r) of the spherical oscillator are: 
        !     phi(r,Omega) = b^(-3/2) * R_nl(r) * Y_ljm(Omega) 
        !     R_nl(r) = N_nl * r**l * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
        !     N_nl    = sqrt(2 * (n-1)!/(n+l-1/2)!)     and    x=r/b
        !     n=1,2,3,...    
        !     R_nl is normalized in such way that the norm integral reads
        !     \int dr r**2 R_nl(r)^2 = 1 
        !     However, in the 
        !     Ref: S.G. Nilsson, Mat-fys Medd 29 (1955)16,1-69
        !     the radial quantum number n=n-1, i.e.,
        !      R_nl(r) = N_nl * r**l * L^(l+1/2)_n(x*x) * exp(-x*x/2), n=0,1,2,...
        !      one should pay attention to
        !-----------------------------------------------------------------------
        ! Note: gm2(n) =  gamma(n+1/2)
        !       fak(n) =  n!
        !       iv(n)  =  (-1)**n
        use Constants, only: igfv
        use Globals, only: gfv, BS
        integer :: n1,l1,l,n2,l2,nn1,nn2,mu1,mu2,tm,smin,smax,s
        real(r64) :: b0, fac1,facs
        b0 = BS%HO_sph%b0
        rl_nl = 0.d0
        ! n = 0, 1, 2    
        nn1 = n1 - 1
        nn2 = n2 - 1
        !---------
        ! l2+l >= l1 >= l2-l
        if(l1.gt.(l2+l).or.l1.lt.l2-l) return
        ! N2 + l >= N1 >= N2-l;  where N = 2n+l
        if((2*nn1+l1).gt.(2*nn2+l2+l).or.(2*nn1+l1).lt.(2*nn2+l2-l)) return
        ! l1-l2+l and l2-l1+l should be even number
        if(l1-l2+l > igfv .or. l2-l1+l>igfv) stop 'ifgv too small!'
        if(gfv%iv(l1-l2+l).ne.1.or.gfv%iv(l2-l1+l).ne.1) return
        !----------
        mu2  = (l1-l2+l)/2   !nu
        mu1  = (l2-l1+l)/2   !nu'
        tm   = (l2+l1+l)/2   !t-1/2
        if(nn1>igfv.or.nn2>igfv.or.nn2+l2+1>igfv.or.nn1+l1+1>igfv.or.mu1>igfv.or.mu2>igfv) stop 'ifgv too small!'
        fac1  = gfv%fak(nn1)*gfv%fak(nn2)/(gfv%gm2(nn2+l2+1)*gfv%gm2(nn1+l1+1))
        facs  = sqrt(fac1)*gfv%fak(mu1)*gfv%fak(mu2)
        smin = max(0,nn1-mu1,nn2-mu2)
        smax = min(nn1,nn2)
        do s = smin,smax
            if(s.lt.0) cycle
            if((nn1-s).lt.0) cycle
            if((nn2-s).lt.0) cycle
            if((s+mu2-nn2).lt.0) cycle
            if((s+mu1-nn1).lt.0) cycle
            if(tm+s+1>igfv.or.s>igfv.or.nn2-s>igfv.or.nn1-s>igfv.or.s+mu2-nn2>igfv.or.s+mu1-nn1>igfv) stop 'ifgv too small!'
            rl_nl = rl_nl + facs*gfv%gm2(tm+s+1)/(gfv%fak(s)*gfv%fak(nn2-s)*gfv%fak(nn1-s)*gfv%fak(s+mu2-nn2)*gfv%fak(s+mu1-nn1))
        enddo
        rl_nl = b0**l*gfv%iv(nn1+nn2)*rl_nl
    end function

    double precision function ylm_ljm(l1,j1,m1,l,m,l2,j2,m2)
        !-----------------------------------------------------------------------
        !     Calculates <l1,J1,M1|Y_(lm)|l2,J2,M2>
        !     J1 = j1 -1/2, M1 = m1 - 1/2 and also for J2,M2 (half-intervalues)
        !-----------------------------------------------------------------------
        use Constants, only: igfv
        use Globals, only: gfv
        integer :: l1,j1,m1,l,m,l2,j2,m2,m1r
        ylm_ljm = 0.d0
        m1r = 1 - m1 ! -M1 + 1/2
        if(j1-m1>igfv) stop 'igfv too small'
        ylm_ljm = gfv%iv(j1-m1)*wigner3j(j1,l,j2,m1r,m,m2,IS=1)*yl_lj(l1,j1,l,l2,j2)
    end function

    double precision function wigner3j(j1,l,j2,m1,m,m2,IS)
        !---------------------------------------------------------------------
        ! Calculates Wigner-3j-coefficient 
        !       ( J1  l  J2 )
        !       ( M1  m  M2 )
        ! icase = 0: (integer values for J1, J2, M1, M2)
        ! where J1 = j1, J2 = j2
        !       M1 = m1, M2 = m2 
        ! icase = 1: (half integer valus for J1, J2, M1, M2)
        ! where J1 = j1 - 1/2, J2 = j2 - 1/2
        !       M1 = m1 - 1/2, M2 = m2 - 1/2
        !-------------------------------------------------------------------
        use Constants, only: igfv
        use Globals,only: gfv
        integer :: j1,l,j2,m1,m,m2,IS
        integer :: i0,i1,i2,i3,i4,i5,n2,n1,n
        wigner3j = 0.d0

        if(IS == 0) then
            if (m1+m+m2.ne.0) return 
            i0 = j1+l+j2+1
            if (igfv.lt.i0) stop 'in wignei: igfv too small'
            i1 = j1+l-j2
            i2 = j1-m1
            i3 = l+m
            i4 = j2-l+m1
            i5 = j2-j1-m
            n2 = min0(i1,i2,i3) 
            n1 = max0(0,-i4,-i5)
            if (n1.gt.n2) return
            do n = n1,n2
                if(n>igfv.or.i1-n>igfv.or.i2-n>igfv.or.i3-n>igfv.or.i4+n>igfv.or.i5+n>igfv) stop 'igfv too small'
                wigner3j = wigner3j + gfv%iv(n)*gfv%fi(n)*gfv%fi(i1-n)*gfv%fi(i2-n)*gfv%fi(i3-n)*gfv%fi(i4+n)*gfv%fi(i5+n)
            end do 
            wigner3j = gfv%iv(i2+i3)*wigner3j*gfv%wfi(i0)*gfv%wf(i1)*gfv%wf(i2+i4)*gfv%wf(i3+i5)*&
                gfv%wf(i2)*gfv%wf(j1+m1)*gfv%wf(l-m)*gfv%wf(i3)*gfv%wf(j2-m2)*gfv%wf(j2+m2)
        else if(IS == 1) then 
            if (m1+m+m2.ne.1) return
            i0 = j1+l+j2
            i1 = j1+l-j2
            i2 = j1-m1
            i3 = l+m
            i4 = j2-l+m1-1
            i5 = j2-j1-m
            n2 = min0(i1,i2,i3) 
            n1 = max0(0,-i4,-i5)
            if (n1.gt.n2) return
            do n = n1,n2
                if(n>igfv.or.i1-n>igfv.or.i2-n>igfv.or.i3-n>igfv.or.i4+n>igfv.or.i5+n>igfv) stop 'igfv too small'
                wigner3j = wigner3j + gfv%iv(n)*gfv%fi(n)*gfv%fi(i1-n)*gfv%fi(i2-n)*gfv%fi(i3-n)*gfv%fi(i4+n)*gfv%fi(i5+n)
            end do
            if(i2+i3>igfv.or.i0>igfv.or.i1>igfv.or.i2+i4>igfv.or.i3+i5>igfv.or.j1+m1-1>igfv.or.&
                i2>igfv.or.i3>igfv.or.l-m>igfv.or.j2+m2-1>igfv.or.j2-m2>igfv) stop 'igfv too small' 
            wigner3j = -gfv%iv(i2+i3)*wigner3j*gfv%wfi(i0)*gfv%wf(i1)*gfv%wf(i2+i4)*gfv%wf(i3+i5)*gfv%wf(j1+m1-1)*&
                gfv%wf(i2)*gfv%wf(i3)*gfv%wf(l-m)*gfv%wf(j2+m2-1)*gfv%wf(j2-m2)
        else 
            stop 'not implement!'
        end if 

        return 
    end function

    double precision function yl_lj(l1,j1,l,l2,j2)
        !----------------------------------------------------------------------
        !     Reduced matrix element of Y_l in jj-coupling                 
        !       < l1,J1 || Y_l || l2,J2 >
        !     where J1 = j1 - 1/2, J2 = j2 - 1/2
        !     parity is not checked
        !     ATTENTION: In our code, ls-coupling (not sl-coupling) is used.
        !     In this case, we need an additional factor (-1)^(j2-j1+l)
        !----------------------------------------------------------------------
        use Constants, only: igfv
        use Globals, only: gfv
        integer :: l1,j1,l,l2,j2
        real(r64) :: a = 0.564189583547756d0 ! a=2/sqrt(4*pi) 
        yl_lj = 0.d0
        if((-1)**(l1+l2+l).ne.1) return
        if(j1>igfv.or.j2>igfv.or.2*l+1>igfv.or.j2-j1+l>igfv) stop 'igfv too small' 
        yl_lj = -a*gfv%sq(j1)*gfv%sq(j2)*gfv%sq(2*l+1)*gfv%iv(j1) &
                *wigner3j(j1,l,j2,0,0,1,IS=1)*gfv%iv(j2-j1+l) ! the last factor comes from ls-coupling
    end function

    subroutine reduced_multipole_matrix_elements(n1,l1,j1,l,n2,l2,j2,rmpm)
        !-----------------------------------------------------------------------
        !      Calculates <n1,l1,J1 || Q_l||n2,l2,J2> 
        !                =<n1,l1| r^l |n2,l2>*<l1,J1|| Y_l ||l2,J2> 
        !      where, Q_(lm)= r^l Y_(lm) 
        !      J1 = j1 -1/2  and also for J2 (half-intervalues)
        !-----------------------------------------------------------------------
        integer :: n1,l1,j1,l,n2,l2,j2
        real(r64) :: rmpm
        rmpm = rl_nl(n1,l1,l,n2,l2)*yl_lj(l1,j1,l,l2,j2)
    end subroutine

    subroutine reduced_monopole_matrix_elements(n1,l1,j1,n2,l2,j2,rdpm)
        !-----------------------------------------------------------------------
        !      Calculates \sqrt{4pi} <n1,l1,J1 || r^2 Y_0||n2,l2,J2> 
        !                = \sqrt{4pi} <n1,l1| r^2 |n2,l2>*<l1,J1|| Y_0 ||l2,J2> 
        !      where, Q_(lm)= r^l Y_(lm)
        !      J1 = j1 -1/2  and also for J2 (half-intervalues)
        !-----------------------------------------------------------------------
        use Constants, only: pi
        integer :: n1,l1,j1,n2,l2,j2
        real(r64) :: rdpm
        
        ! Way1
        rdpm = sqrt(4.d0*pi)*rl_nl(n1,l1,2,n2,l2)*yl_lj(l1,j1,0,l2,j2)

        ! Way2
        ! if (n1==n2 .and. j1==j2) then
        !     rdpm = rl_nl(n1,l1,2,n2,l2)*sqrt(2.d0*j1)
        ! else 
        !     rdpm = 0.d0
        ! end if
    end subroutine
end Module EM