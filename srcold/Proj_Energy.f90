Module Energy 
    use Constants, only: r64,ngr,ntheta,nphi,hbc
    implicit none
    contains
    
    subroutine calculate_sigma_nabla_Spherical
        !--------------------------------------------------------------------------------------------------------------------------
        ! calculates matrix elements:
        !   <n|\sigma\nabla|m> =  \int dr \Psi_{n}^\dagger(r,\sigma) \boldsymbol{\sigma}\cdot\boldsymbol{\nabla} \Psi_{m}(r,\sigma)
        ! for Fermions in the Spherical oscillator wave functions. 
        !
        !The wave functions for the spherical oscillator is 
        !    |m> = |n_r l j m_j> 
        !        = \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s} R_{n_r l}(r) Y_{l m_l}(\theta,\phi) \chi_{m_{s}}^{\sigma}
        !    where C^{j m_j}_{l m_l 1/2 m_s} is Clebsch–Gordan coefficients.
        !
        ! Note:
        !     m is the wave functions basis of large component.
        !     n is the wave functions basis of small component.
        ! Output:
        !     ME%spds(m1,m2) = - <m1|\sigma\nabla|m2> 
        !
        !---------------------------------------------------------------------------------------------------------------------------
        use Constants, only: zero,r64,ngr,nr_max_sph,l_max_sph
        use Globals, only: BS,spatial_grid,ME
        integer :: ib,nf,ng,i0f,i0g,nr2,nl2,nj2,nm2,nl2p,nr1,nl1,nj1,nm1,ix,m1,m2
        real(r64) :: radint,s,x
        ib = 1
        nf = BS%HO_sph%idsp(ib,1)
        ng = BS%HO_sph%idsp(ib,2)
        i0f =BS%HO_sph%iasp(ib,1)
        i0g =BS%HO_sph%iasp(ib,2)
        allocate(ME%spds(ng,nf))
        do m2 = 1,nf ! large component
            nr2 = BS%HO_sph%nljm(i0f+m2,1)
            nl2 = BS%HO_sph%nljm(i0f+m2,2)
            nj2 = BS%HO_sph%nljm(i0f+m2,3)
            nm2 = BS%HO_sph%nljm(i0f+m2,4)      
            nl2p = 2*nj2 - nl2 -1 
            do m1 = 1,ng ! small component
                nr1 = BS%HO_sph%nljm(i0g+m1,1)
                nl1 = BS%HO_sph%nljm(i0g+m1,2)
                nj1 = BS%HO_sph%nljm(i0g+m1,3)
                nm1 = BS%HO_sph%nljm(i0g+m1,4) 
                ! D. A. Varshalovich, 1988, page 200 
                radint = zero
                s = zero
                if(nm1.ne.nm2.or.nj1.ne.nj2) goto 10
                if(nl1.ne.nl2p) goto 10   
                if(nr1.gt.nr_max_sph.or.nl1.gt.l_max_sph) stop 'nrx or nlx is too small'
                s = -1.d0
                do ix = 1,ngr
                    x = spatial_grid%x(ix)
                    radint = radint + x**2* BS%HO_sph%Rnl(ix,nr1,nl1)*BS%HO_sph%Rnl1(ix,nr2,nl2) &
                            -((nj2-0.5d0)*(nj2+0.5d0)-nl2*(nl2+1.d0)-0.75)*x &
                            *BS%HO_sph%Rnl(ix,nr1,nl1)*BS%HO_sph%Rnl(ix,nr2,nl2)
                enddo !ir  
            10  ME%spds(m1,m2) = -s* radint 
            enddo !m1
        enddo !m2
        return
    end subroutine

    subroutine Kinetic_term(iphi_n,iphi_p,E_kin,pE_kin)
        !--------------------------------------------------------------------------------------
        !  calculate kinetic energy term:
        !-------------------------------------------------------------------------------------
        use Globals,only: BS,ME,force,mix
        integer,intent(in) :: iphi_n,iphi_p
        complex(r64), intent(out) :: E_kin,pE_kin
        integer :: nf,ng,m1,m2,i
        complex(r64) :: E_kin_n,E_kin_p,pE_kin_n,pE_kin_p
        E_kin_n = (0.0d0,0.0d0)
        E_kin_p = (0.0d0,0.0d0)
        pE_kin_n = (0.0d0,0.0d0)
        pE_kin_p = (0.0d0,0.0d0)
        nf = BS%HO_sph%idsp(1,1)
        ng = BS%HO_sph%idsp(1,2)
        ! -\sum_{mn} {rho_mn^{+-} \int [ \Psi_n^\dagger sigma nabla \Psi_{m} ]dr}  + \sum_{mn}{rho_nm^{-+} \int [ \Psi_m^\dagger sigma nabla \Psi_{n} ]dr}
        ! Note:
        ! 1) -\int [\Psi_n^\dagger sigma nabla \Psi_{m} ]dr have been calculated and stored in ME%spds
        ! 2)  \int [ \Psi_m^\dagger sigma nabla \Psi_{n} ]dr = -\int [\Psi_n^\dagger sigma nabla \Psi_{m} ]dr
        ! 3)  Because the nabla operator of `ME%spds` is actually a partial derivative with respect to x, 
        !     and since \partial r = b \partial x, it needs to be divided by b.
        do m1 = 1, nf
            do m2 = 1, ng
                E_kin_n = E_kin_n + (mix%rho_mm(m1,m2,3,iphi_n,1) + mix%rho_mm(m2,m1,4,iphi_n,1))*ME%spds(m2,m1)/BS%HO_sph%b0
                E_kin_p = E_kin_p + (mix%rho_mm(m1,m2,3,iphi_p,2) + mix%rho_mm(m2,m1,4,iphi_p,2))*ME%spds(m2,m1)/BS%HO_sph%b0
                pE_kin_n = pE_kin_n + (mix%prho_mm(m1,m2,3,iphi_n,1) + mix%prho_mm(m2,m1,4,iphi_n,1))*ME%spds(m2,m1)/BS%HO_sph%b0
                pE_kin_p = pE_kin_p + (mix%prho_mm(m1,m2,3,iphi_p,2) + mix%prho_mm(m2,m1,4,iphi_p,2))*ME%spds(m2,m1)/BS%HO_sph%b0
            end do 
        end do
        ! \int  ( m*\tilde{rho}_S -m*\tilde{rho}_V ) dr
        do i = 1, ngr*ntheta*nphi
            E_kin_n = E_kin_n - force%masses%amu*mix%rho_V_it(i,iphi_n,1) + force%masses%amu*mix%rho_S_it(i,iphi_n,1)
            E_kin_p = E_kin_p - force%masses%amu*mix%rho_V_it(i,iphi_p,2) + force%masses%amu*mix%rho_S_it(i,iphi_p,2)
            pE_kin_n = pE_kin_n - force%masses%amu*mix%prho_V_it(i,iphi_n,1) + force%masses%amu*mix%prho_S_it(i,iphi_n,1)
            pE_kin_p = pE_kin_p - force%masses%amu*mix%prho_V_it(i,iphi_p,2) + force%masses%amu*mix%prho_S_it(i,iphi_p,2)
        end do
        E_kin = (E_kin_n + E_kin_p)*hbc
        pE_kin = (pE_kin_n + pE_kin_p)*hbc
    end subroutine
    subroutine EDF_terms(iphi_n,iphi_p,E_EDF,pE_EDF)
        !--------------------------------------------------------------------------------------
        !  calculate energy density functional terms:
        !
        !  
        ! -------------------------------------------------------------------------------------------------------------
        !  For a function F, the integral in coordinate space is:
        !  \int F(\boldsymbol{r}) d\bodysymbol{r}
        !       = \int F(r,theta,phi) r^2 sin(theta) dr dtheta dphi
        !       = \sum F(r,theta,phi) r^2 sin(theta) w(r) w(theta) w(phi)
        !       = \sum F(r,theta,phi)/b^3*x^2*sin(theta)*w(x)*w(theta)*w(phi)
        !  Note the mixed density, current and tensor are multiplied by x^2 sin(theta wx wtheta wphi) and divided by b^(-3).
        !  
        !  \int F^2(\boldsymbol{r}) d\bodysymbol{r}
        !       = \sum F^2(r,theta,phi)/b^3*x^2*sin(theta)*w(x)*w(theta)*w(phi)
        !       = \sum (F(r,theta,phi)/b^3*x^2 sin(theta) w(x) w(theta) w(phi))^2 / (r^2 sin(theta)w(r) w(theta) w(phi))
        !--------------------------------------------------------------------------------------------------------------
        use Globals, only: force,spatial_grid,mix
        integer,intent(in) :: iphi_n,iphi_p
        complex(r64),intent(out) :: E_EDF,pE_EDF
        integer :: i
        complex(r64) :: E_S,E_V,E_TS,E_TV,E_cur,E_nl,E_dS,E_dV,E_dTS,E_dTV,&
                        rho_S,rho_V,rho_TS,rho_TV,j_V1,j_V2,j_V3,j_V_square,j_TV1,j_TV2,j_TV3,j_TV_square,&
                        d2rho_S,d2rho_V,d2rho_TS,d2rho_TV,&
                        pE_S,pE_V,pE_TS,pE_TV,pE_cur,pE_nl,pE_dS,pE_dV,pE_dTS,pE_dTV,&
                        prho_S,prho_V,prho_TS,prho_TV,pj_V1,pj_V2,pj_V3,pj_V_square,pj_TV1,pj_TV2,pj_TV3,pj_TV_square,&
                        pd2rho_S,pd2rho_V,pd2rho_TS,pd2rho_TV
        E_S = (0.0d0,0.0d0)
        E_V = (0.0d0,0.0d0)
        E_TS = (0.0d0,0.0d0)
        E_TV = (0.0d0,0.0d0)
        E_cur = (0.0d0,0.0d0)
        E_nl = (0.0d0,0.0d0)
        E_dS = (0.0d0,0.0d0)
        E_dV = (0.0d0,0.0d0)
        E_dTS = (0.0d0,0.0d0)
        E_dTV = (0.0d0,0.0d0)
        ! 
        pE_S = (0.0d0,0.0d0)
        pE_V = (0.0d0,0.0d0)
        pE_TS = (0.0d0,0.0d0)
        pE_TV = (0.0d0,0.0d0)
        pE_cur = (0.0d0,0.0d0)
        pE_nl = (0.0d0,0.0d0)
        pE_dS = (0.0d0,0.0d0)
        pE_dV = (0.0d0,0.0d0)
        pE_dTS = (0.0d0,0.0d0)
        pE_dTV = (0.0d0,0.0d0)
        do i = 1,ngr*ntheta*nphi
            ! Scalar, vector, tensor scalar, and tensor vector densities
            rho_S = mix%rho_S_it(i,iphi_n,1) + mix%rho_S_it(i,iphi_p,2)
            rho_V = mix%rho_V_it(i,iphi_n,1) + mix%rho_V_it(i,iphi_p,2)
            rho_TS = mix%rho_S_it(i,iphi_n,1) - mix%rho_S_it(i,iphi_p,2)
            rho_TV = mix%rho_V_it(i,iphi_n,1) - mix%rho_V_it(i,iphi_p,2)
            !
            prho_S = mix%prho_S_it(i,iphi_n,1) + mix%prho_S_it(i,iphi_p,2)
            prho_V = mix%prho_V_it(i,iphi_n,1) + mix%prho_V_it(i,iphi_p,2)
            prho_TS = mix%prho_S_it(i,iphi_n,1) - mix%prho_S_it(i,iphi_p,2)
            prho_TV = mix%prho_V_it(i,iphi_n,1) - mix%prho_V_it(i,iphi_p,2)
            ! current
            j_V1 = mix%j_V1_it(i,iphi_n,1) + mix%j_V1_it(i,iphi_p,2)
            j_V2 = mix%j_V2_it(i,iphi_n,1) + mix%j_V2_it(i,iphi_p,2)
            j_V3 = mix%j_V3_it(i,iphi_n,1) + mix%j_V3_it(i,iphi_p,2)
            j_V_square = j_V1**2 + j_V2**2 + j_V3**2
            j_TV1 = mix%j_V1_it(i,iphi_n,1) - mix%j_V1_it(i,iphi_p,2)
            j_TV2 = mix%j_V2_it(i,iphi_n,1) - mix%j_V2_it(i,iphi_p,2)
            j_TV3 = mix%j_V3_it(i,iphi_n,1) - mix%j_V3_it(i,iphi_p,2)
            j_TV_square = j_TV1**2 + j_TV2**2 + j_TV3**2
            !
            pj_V1 = mix%pj_V1_it(i,iphi_n,1) + mix%pj_V1_it(i,iphi_p,2)
            pj_V2 = mix%pj_V2_it(i,iphi_n,1) + mix%pj_V2_it(i,iphi_p,2)
            pj_V3 = mix%pj_V3_it(i,iphi_n,1) + mix%pj_V3_it(i,iphi_p,2)
            pj_V_square = pj_V1**2 + pj_V2**2 + pj_V3**2
            pj_TV1 = mix%pj_V1_it(i,iphi_n,1) - mix%pj_V1_it(i,iphi_p,2)
            pj_TV2 = mix%pj_V2_it(i,iphi_n,1) - mix%pj_V2_it(i,iphi_p,2)
            pj_TV3 = mix%pj_V3_it(i,iphi_n,1) - mix%pj_V3_it(i,iphi_p,2)
            pj_TV_square = pj_TV1**2 + pj_TV2**2 + pj_TV3**2
            ! second derivative
            d2rho_S = mix%d2rho_S_it(i,iphi_n,1) + mix%d2rho_S_it(i,iphi_p,2)
            d2rho_V = mix%d2rho_V_it(i,iphi_n,1) + mix%d2rho_V_it(i,iphi_p,2)
            d2rho_TS = mix%d2rho_S_it(i,iphi_n,1) - mix%d2rho_S_it(i,iphi_p,2)
            d2rho_TV = mix%d2rho_V_it(i,iphi_n,1) - mix%d2rho_V_it(i,iphi_p,2)
            !
            pd2rho_S = mix%pd2rho_S_it(i,iphi_n,1) + mix%pd2rho_S_it(i,iphi_p,2)
            pd2rho_V = mix%pd2rho_V_it(i,iphi_n,1) + mix%pd2rho_V_it(i,iphi_p,2)
            pd2rho_TS = mix%pd2rho_S_it(i,iphi_n,1) - mix%pd2rho_S_it(i,iphi_p,2)
            pd2rho_TV = mix%pd2rho_V_it(i,iphi_n,1) - mix%pd2rho_V_it(i,iphi_p,2)


            E_S = E_S + 1.d0/2.d0*force%couplg%ggsig*rho_S**2/spatial_grid%wwsp(i)
            E_V = E_V + 1.d0/2.d0*force%couplg%ggome*rho_V**2/spatial_grid%wwsp(i)
            E_TS = E_TS + 1.d0/2.d0*force%couplg%ggdel*rho_TS**2/spatial_grid%wwsp(i)
            E_TV = E_TV + 1.d0/2.d0*force%couplg%ggrho*rho_TV**2/spatial_grid%wwsp(i)
            E_cur = E_cur + 1.d0/2.d0*force%couplg%ggome*j_V_square/spatial_grid%wwsp(i) &
                        + 1.d0/2.d0*force%couplg%ggrho*j_TV_square/spatial_grid%wwsp(i) &
                        - 1.d0/4.d0*force%coupnl%gggamv*j_V_square**2/spatial_grid%wwsp(i)**3 &
                        + 1.d0/2.d0*force%coupnl%gggamv*rho_V**2*j_V_square/spatial_grid%wwsp(i)**3
            E_nl = E_nl + 1.d0/3.d0*force%coupnl%ggbet*rho_S**3/spatial_grid%wwsp(i)**2 &
                        + 1.d0/4.d0*force%coupnl%gggams*rho_S**4/spatial_grid%wwsp(i)**3 &
                        + 1.d0/4.d0*force%coupnl%gggamv*rho_V**4/spatial_grid%wwsp(i)**3
            E_dS = E_dS + 1.d0/2.d0*force%coupld%ddsig*rho_S*d2rho_S/spatial_grid%wwsp(i)
            E_dV = E_dV + 1.d0/2.d0*force%coupld%ddome*rho_V*d2rho_V/spatial_grid%wwsp(i)
            E_dTS = E_dTS + 1.d0/2.d0*force%coupld%dddel*rho_TS*d2rho_TS/spatial_grid%wwsp(i)
            E_dTV = E_dTV + 1.d0/2.d0*force%coupld%ddrho*rho_TV*d2rho_TV/spatial_grid%wwsp(i)
            !
            pE_S = pE_S + 1.d0/2.d0*force%couplg%ggsig*prho_S**2/spatial_grid%wwsp(i)
            pE_V = pE_V + 1.d0/2.d0*force%couplg%ggome*prho_V**2/spatial_grid%wwsp(i)
            pE_TS = pE_TS + 1.d0/2.d0*force%couplg%ggdel*prho_TS**2/spatial_grid%wwsp(i)
            pE_TV = pE_TV + 1.d0/2.d0*force%couplg%ggrho*prho_TV**2/spatial_grid%wwsp(i)
            pE_cur = pE_cur + 1.d0/2.d0*force%couplg%ggome*pj_V_square/spatial_grid%wwsp(i) &
                        + 1.d0/2.d0*force%couplg%ggrho*pj_TV_square/spatial_grid%wwsp(i) &
                        - 1.d0/4.d0*force%coupnl%gggamv*pj_V_square**2/spatial_grid%wwsp(i)**3 &
                        + 1.d0/2.d0*force%coupnl%gggamv*prho_V**2*pj_V_square/spatial_grid%wwsp(i)**3
            pE_nl = pE_nl + 1.d0/3.d0*force%coupnl%ggbet*prho_S**3/spatial_grid%wwsp(i)**2 &
                        + 1.d0/4.d0*force%coupnl%gggams*prho_S**4/spatial_grid%wwsp(i)**3 &
                        + 1.d0/4.d0*force%coupnl%gggamv*prho_V**4/spatial_grid%wwsp(i)**3
            pE_dS = pE_dS + 1.d0/2.d0*force%coupld%ddsig*prho_S*pd2rho_S/spatial_grid%wwsp(i)
            pE_dV = pE_dV + 1.d0/2.d0*force%coupld%ddome*prho_V*pd2rho_V/spatial_grid%wwsp(i)
            pE_dTS = pE_dTS + 1.d0/2.d0*force%coupld%dddel*prho_TS*pd2rho_TS/spatial_grid%wwsp(i)
            pE_dTV = pE_dTV + 1.d0/2.d0*force%coupld%ddrho*prho_TV*pd2rho_TV/spatial_grid%wwsp(i)
        end do
        E_EDF = (E_S + E_V + E_TS + E_TV - E_cur + E_nl + E_dS + E_dV + E_dTS + E_dTV )*hbc
        pE_EDF = (pE_S + pE_V + pE_TS + pE_TV - pE_cur + pE_nl + pE_dS + pE_dV + pE_dTS + pE_dTV )*hbc
    end subroutine
    subroutine Coulomb_term(iphi_p,E_cou,pE_cou)
        !--------------------------------------------------------------------------------------
        !  calculate Coulomb energy term:
        !
        !  E_cou = 1/2*e\int V_cou(\boldsymbol{r}) \rho^p_V(\boldsymbol{r}) d\boldsymbol{r}
        !
        !--------------------------------------------------------------------------------------
        use Constants, only: alphi
        use Globals, only: BS,spatial_grid,mix
        integer :: iphi_p
        complex(r64) :: E_cou,pE_cou
        integer :: ix1,itheta1,iphi1,i1,ix2,itheta2,iphi2,i2,i
        real(r64) :: b0,x1,r1,theta1,phi1,x2,r2,theta2,phi2,cosOmega12,r12_square,r12
        complex(r64), dimension(ngr*ntheta*nphi) :: V_cou,pV_cou
        E_cou = (0.0d0,0.0d0)
        pE_cou = (0.0d0,0.0d0)
        ! calculate scalar potential V_cou
        b0 = BS%HO_sph%b0
        do ix1 = 1, ngr
            x1 = spatial_grid%x(ix1) ! r/b
            r1 = x1*b0
            do itheta1 = 1, ntheta
                theta1 = spatial_grid%theta(itheta1) ! theta
                do iphi1 = 1, nphi 
                    phi1 = spatial_grid%phi(iphi1) ! phi
                    i1 = ix1 + (itheta1-1)*ngr + (iphi1-1)*ngr*ntheta
                    V_cou(i1) = (0.0d0,0.0d0)
                    pV_cou(i1) = (0.0d0,0.0d0)
                    do ix2 = 1, ngr
                        x2 = spatial_grid%x(ix2) ! r/b
                        r2 = x2*b0
                        do itheta2 = 1, ntheta
                            theta2 = spatial_grid%theta(itheta2) ! theta
                            do iphi2 = 1, nphi 
                                phi2 = spatial_grid%phi(iphi2) ! phi
                                i2 = ix2 + (itheta2-1)*ngr + (iphi2-1)*ngr*ntheta
                                ! cos(\Omega_12) = cos(theta1)cos(theta2) + sin(theta1)sin(theta2)cos(phi1-phi2)
                                cosOmega12 = dcos(theta1)*dcos(theta2) + dsin(theta1)*dsin(theta2)*dcos(phi1-phi2)
                                ! |r1 -r2|^2 = r1^2 + r2^2 - 2r1r2cos(\Omega_12)
                                r12_square = r1**2 + r2**2 - 2.d0*r1*r2*cosOmega12
                                if(r12_square < 0.0d0) r12_square = 0.01d0
                                r12 = dsqrt(r12_square)
                                V_cou(i1) = V_cou(i1) + 1.d0/(2.d0*alphi)*r12*mix%d2rho_V_it(i2,iphi_p,2)
                                pV_cou(i1) = pV_cou(i1) + 1.d0/(2.d0*alphi)*r12*mix%pd2rho_V_it(i2,iphi_p,2)
                            end do 
                        end do 
                    end do 
                end do 
            end do
        end do 
        do i = 1, ngr*ntheta*nphi
            E_cou = E_cou + 1.d0/2.d0*V_cou(i)*mix%rho_V_it(i,iphi_p,2)
            pE_cou = pE_cou + 1.d0/2.d0*pV_cou(i)*mix%prho_V_it(i,iphi_p,2)
        end do
        E_cou = E_cou*hbc
        pE_cou = pE_cou*hbc
    end subroutine
    subroutine Pairing_term(iphi_n,iphi_p,E_pair,pE_pair)
        !----------------------------------------------------------------------------------------------------------
        !  calculate pairing energy term:
        !  E_pair^\tau(q1, q2, alpha, beta, gamma, phi_\tau) 
        !              = 1/4 \int \kappa*_\tau \kappa_\tau V_\tau d\boldsymbol{r}
        !  where \kappa = \kappa(\boldsymbol{r}, q1, q2, alpha, beta, gamma, phi_\tau)
        !         \tau = n is neutron and \tau = p is proton
        !  E_pair = E_pair^n + E_pair^p
        !-----------------------------------------------------------------------------------------------------------
        use Globals, only: pairing,spatial_grid,mix,option
        use BCS, only: initial_pairing_field
        integer,intent(in) :: iphi_n,iphi_p
        complex(r64),intent(out) :: E_pair,pE_pair
        integer :: i
        complex(r64) :: E_pair_n,E_pair_p,pE_pair_n,pE_pair_p
        if(option%CDFTType ==0) call initial_pairing_field(.False.)
        E_pair_n = (0.0d0,0.0d0)
        E_pair_p = (0.0d0,0.0d0)
        pE_pair_n = (0.0d0,0.0d0)
        pE_pair_p = (0.0d0,0.0d0)
        do i = 1, ngr*ntheta*nphi
            E_pair_n = E_pair_n + 1.d0/4.d0*mix%kappac_it(i,iphi_n,1)*mix%kappa_it(i,iphi_n,1)*(-pairing%vpair(1))/spatial_grid%wwsp(i)
            E_pair_p = E_pair_p + 1.d0/4.d0*mix%kappac_it(i,iphi_p,2)*mix%kappa_it(i,iphi_p,2)*(-pairing%vpair(2))/spatial_grid%wwsp(i)
            pE_pair_n = pE_pair_n + 1.d0/4.d0*mix%pkappac_it(i,iphi_n,1)*mix%pkappa_it(i,iphi_n,1)*(-pairing%vpair(1))/spatial_grid%wwsp(i)
            pE_pair_p = pE_pair_p + 1.d0/4.d0*mix%pkappac_it(i,iphi_p,2)*mix%pkappa_it(i,iphi_p,2)*(-pairing%vpair(2))/spatial_grid%wwsp(i)
        end do
        E_pair = E_pair_n + E_pair_p
        pE_pair = pE_pair_n + pE_pair_p
    end subroutine
    subroutine Center_of_Mass_Correction_term(E_cm)
        !--------------------------------------------------------------------------------------
        !  calculate center of mass correction energy term:
        !--------------------------------------------------------------------------------------
        use Globals,only: Proj_option,wf1,wf2,BS
        complex(r64),intent(out) :: E_cm
        ! zero
        if(Proj_option%icm==0) then
            E_cm = (0.0d0,0.0d0)
        !q-dependent. average 
        else if ( Proj_option%icm==1 ) then
            E_cm = (wf1%ecm + wf2%ecm)/2.d0
        else if (Proj_option%icm==2) then 
            E_cm = -0.75d0*BS%HO_sph%hom
        else 
            write(*,*) "[Center_of_Mass_Correction_term] Worning : Invalid value for icm "
        end if
    end subroutine
End Module