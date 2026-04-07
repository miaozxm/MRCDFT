Module Eccentricity
    use Constants,only: r64
    implicit none
    integer, parameter :: max_n = 2 ! max n of f_n
    real(r64), dimension(:,:,:,:,:), allocatable:: f_n_array ! f_n_array(ifg1, m1, ifg2, m2, n)
    real(r64), dimension(:,:,:,:,:,:), allocatable:: Q_array ! Q_array(ifg1,m1,ifg2,m2,n,mu)
    logical, dimension(-max_n:max_n) :: precomputed_n = .False.
contains

    subroutine calculate_Eccentri_n(n,iphi,it,Eccentri,pEccentri,Q_mu,pQ_mu)
        !--------------------------------------------------
        !  calculate <q1| E_n R|q2>/<q1|R|q2>
        !  where R = R(alpha,beta,gamma, phi_n,phi_p) = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
        !--------------------------------------------------
        use Constants,only: itx
        use Globals, only: BS, mix
        integer :: n,iphi,it
        integer :: ifg1,ifg2,ifg3,ifg4,m1,m2,mu,m3,m4,iter,total_iter
        complex(r64) :: Eccentri(2),pEccentri(2),Q_mu(-n:n),pQ_mu(-n:n)
        real(r64) :: e_1B,e_2B 
        Eccentri = (0.d0,0.d0)
        pEccentri = (0.d0,0.d0)
        Q_mu = (0.d0,0.d0)
        pQ_mu = (0.d0,0.d0)

        if(.not. precomputed_n(n)) then 
            if(n>max_n) then 
                write(*,*) 'Warning: max_n is small than input n ! f_n will not be precomputed.'
            else 
                call precompute_f_n(n)
                call precompute_f_n(-n)
                call precompute_Q(n)
                precomputed_n(n) = .True.
            end if
        end if

        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
        !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2,mu,ifg3,m3,ifg4,m4,e_1B,e_2B) &
        !$omp reduction(+:Eccentri,pEccentri,Q_mu,pQ_mu)
        !$OMP DO COLLAPSE(1) SCHEDULE(static)
        do iter = 1,total_iter
            if(iter <= BS%HO_sph%idsp(1,1)) then
                ifg1 = 1
                m1 = iter
            else
                ifg1 = 2
                m1 = iter - BS%HO_sph%idsp(1,1)
            end if
            do ifg2 = 1,2
            do m2 = 1,BS%HO_sph%idsp(1,ifg2)
                ! 1 Body
                if(ifg1==ifg2) then
                    call eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,n,e_1B)
                    Eccentri(1) = Eccentri(1) + e_1B*mix%rho_mm(m2,m1,ifg1,iphi,it)
                    pEccentri(1) = pEccentri(1) + e_1B*mix%prho_mm(m2,m1,ifg1,iphi,it)
                end if 
                ! to calculate nn-pp and pp-nn term
                if(ifg1==ifg2) then
                    do mu = -n, n 
                        Q_mu(mu) = Q_mu(mu) +  Q(ifg1,m1,ifg2,m2,n,mu)*mix%rho_mm(m2,m1,ifg1,iphi,it)
                        pQ_mu(mu) = pQ_mu(mu) +  Q(ifg1,m1,ifg2,m2,n,mu)*mix%rho_mm(m2,m1,ifg1,iphi,it)
                    end do
                end if 
                do ifg3 = 1,2
                do m3 = 1,BS%HO_sph%idsp(1,ifg3)
                    do ifg4 = 1,2
                    do m4 = 1,BS%HO_sph%idsp(1,ifg4)
                        ! 2 Body
                        if(ifg1==ifg3 .and. ifg2==ifg4) then
                            call eccentricity_matrix_element_two_body(ifg1,m1,ifg2,m2,ifg4,m4,ifg3,m3,n,e_2B)
                            !  e_{m1 m2 m4 m3}*<c^+_{m1}c^+_{m2}c_{m3}c_{m4}> = e_{m1 m2 m4 m3}*rho_4321
                            Eccentri(2) = Eccentri(2)+ e_2B*( & 
                                          mix%rho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%rho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                        - mix%rho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%rho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                        + mix%kappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%kappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it))
                            pEccentri(2) = pEccentri(2) + e_2B*( &
                                          mix%prho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%prho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                        - mix%prho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%prho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                        + mix%pkappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%pkappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it))
                        end if 

                    end do
                    end do
                end do 
                end do 
            end do 
            end do 
        end do 
        !$OMP END PARALLEL
    end subroutine

    subroutine eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,n,EME1B)
        !-----------------------------------------------------------
        !   e^{n}_{m1 m2} = \sum_m f^n_{m1 m} f^{-n}_{m m2} 
        ! 
        !   e^{n}_{p q}  = \sum_{m p'} (-1)^m (Q_nm)_{pp'}(Q_n-m)_{p'q}  
        !-----------------------------------------------------------
        use Globals, only: BS
        integer,intent(in) :: ifg1,m1,ifg2,m2,n
        real(r64) :: EME1B
        integer :: ndsp,ifg,m,mu
        EME1B = 0.d0
        ifg = ifg1
        ndsp = BS%HO_sph%idsp(1,ifg)

        ! do m = 1, ndsp
        !     EME1B = EME1B + f_n(ifg1,m1,ifg,m,n)*f_n(ifg,m,ifg2,m2,-n)
        ! end do
        
        EME1B = 0.d0
        do mu = -n,n
            do m = 1, ndsp
                EME1B = EME1B + (-1)**mu*Q(ifg1,m1,ifg,m,n,mu)*Q(ifg,m,ifg2,m2,n,-mu)
            end do 
        end do 
    end subroutine

    subroutine eccentricity_matrix_element_two_body(ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,n,EME2B)
        !-------------------------------------------------------------------------------------
        !    
        !    e^{n}_{m1 m2 m3 m4} = - f^n_{m1 m4}f^{-n}_{m2 m3}
        ! or in antisymmetrized form:
        !    e^{n}_{m1 m2 m3 m4} = 1/4 * 2(f^n_{m1 m3}f^{-n}_{m2 m4} - f^{-n}_{m1 m4}f^n_{m2 m3}
        !
        !    e^{n}_{1 2 3 4} = - \sum_{m} (-1)^m (Q_nm)_{1 4} (Q_n-m)_{2 3} 
        !-------------------------------------------------------------------------------------
        use Globals, only: BS
        integer :: ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,n,mu
        real(r64) :: EME2B
        ! EME2B = - f_n(ifg1,m1,ifg4,m4,n)*f_n(ifg2,m2,ifg3,m3,-n)
        ! EME2B = 1.d0/4.d0*2.d0*(f_n(ifg1,m1,ifg3,m3,n)*f_n(ifg2,m2,ifg4,m4,-n)-f_n(ifg1,m1,ifg4,m4,-n)*f_n(ifg2,m2,ifg3,m3,n)) 
        EME2B = 0.d0
        do mu = -n,n
            EME2B = EME2B - (-1)**mu * Q(ifg1,m1,ifg4,m4,n,mu)*Q(ifg2,m2,ifg3,m3,n,-mu) 
        end do
    end subroutine

    double precision function f_n(ifg1,m1,ifg2,m2,n)
        integer,intent(in) :: ifg1,ifg2, m1, m2, n
        if(precomputed_n(n)) then 
            f_n = f_n_array(ifg1,m1,ifg2,m2,n)
        else
            f_n = compute_Q(ifg1,m1,ifg2,m2,abs(n),n)
        end if 
    end

    double precision function Q(ifg1,m1,ifg2,m2,n,mu)
        integer,intent(in) :: ifg1,ifg2, m1, m2, n, mu
        if(precomputed_n(n)) then 
            Q = Q_array(ifg1,m1,ifg2,m2,n,mu)
        else
            Q = compute_Q(ifg1,m1,ifg2,m2,n,mu)
        end if 
    end

    subroutine precompute_f_n(n)
        !----------------------------------------------------
        ! calculate <nr1 nl1 nj1 nm1| F_n | nr1 nl2 nj2 nm2> 
        ! and store to f_n_array
        !----------------------------------------------------
        use Globals, only: BS
        integer,intent(in) :: n
        integer :: ifg1,m1,ifg2,m2,dim_m_max,iter,total_iter
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        if(.not. allocated(f_n_array)) allocate(f_n_array(2,dim_m_max,2,dim_m_max,-max_n:max_n),source=0.d0)

        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
        !$omp parallel do collapse(1) default(shared) &
        !$omp private(iter,ifg1,m1,ifg2,m2) schedule(static)
        do iter = 1,total_iter
            if(iter <= BS%HO_sph%idsp(1,1)) then
                ifg1 = 1
                m1 = iter
            else
                ifg1 = 2
                m1 = iter - BS%HO_sph%idsp(1,1)
            end if
                do ifg2 = 1,2
                    do m2 = 1, BS%HO_sph%idsp(1,ifg2)
                        f_n_array(ifg1,m1,ifg2,m2,n) =  compute_Q(ifg1,m1,ifg2,m2,abs(n),n)
                    end do
                end do
        end do
        !$omp end parallel do
    end subroutine


    subroutine precompute_Q(n)
        use Globals, only: BS
        integer,intent(in) :: n
        integer :: ifg1,m1,ifg2,m2,mu,dim_m_max,iter,total_iter
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        if(.not. allocated(Q_array)) allocate(Q_array(2,dim_m_max,2,dim_m_max,max_n,-max_n:max_n),source=0.d0)

        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
        !$omp parallel do collapse(1) default(shared) &
        !$omp private(iter,ifg1,m1,ifg2,m2,mu) schedule(static)
        do iter = 1,total_iter
            if(iter <= BS%HO_sph%idsp(1,1)) then
                ifg1 = 1
                m1 = iter
            else
                ifg1 = 2
                m1 = iter - BS%HO_sph%idsp(1,1)
            end if
                do ifg2 = 1,2
                    do m2 = 1, BS%HO_sph%idsp(1,ifg2)
                        do mu = -n, n
                            Q_array(ifg1,m1,ifg2,m2,n,mu) = compute_Q(ifg1,m1,ifg2,m2,n,mu)
                        end do
                    end do
                end do
        end do
        !$omp end parallel do
    end subroutine

    double precision function compute_Q(ifg1,m1,ifg2,m2,n,mu)
        !----------------------------------------------------
        !   calculate <nr1 nl1 nj1 nm1| r^n Y_{n mu} | nr1 nl2 nj2 nm2>
        !---------------------------------------------------
        use Globals, only: BS
        use EM, only: multipole_matrix_elements
        integer,intent(in) :: ifg1,ifg2, m1, m2, n, mu
        integer :: i0sp1,nr1,nl1,nj1,nm1,i0sp2,nr2,nl2,nj2,nm2
        real(r64) :: mpme

        i0sp1 = BS%HO_sph%iasp(1,ifg1)
        nr1 = BS%HO_sph%nljm(i0sp1+m1,1)
        nl1 = BS%HO_sph%nljm(i0sp1+m1,2)
        nj1 = BS%HO_sph%nljm(i0sp1+m1,3)
        nm1 = BS%HO_sph%nljm(i0sp1+m1,4)

        i0sp2 = BS%HO_sph%iasp(1,ifg2)
        nr2 = BS%HO_sph%nljm(i0sp2+m2,1)
        nl2 = BS%HO_sph%nljm(i0sp2+m2,2)
        nj2 = BS%HO_sph%nljm(i0sp2+m2,3)
        nm2 = BS%HO_sph%nljm(i0sp2+m2,4)

        call multipole_matrix_elements(nr1,nl1,nj1,nm1,n,mu,nr2,nl2,nj2,nm2,mpme)
        compute_Q =  mpme
    end function

    integer function indexfg(ifg1,ifg2)
        integer, intent(in) :: ifg1, ifg2
        if(ifg1==1.and.ifg2==1) then
            indexfg = 1 ! ++
        else if(ifg1==2.and.ifg2==2) then
            indexfg = 2 ! --
        else if(ifg1==1.and.ifg2==2) then
            indexfg = 3 ! +-
        else if(ifg1==2.and.ifg2==1) then
            indexfg = 4 ! -+
        else 
            stop 'wrong ifg1 and ifg2'
        end if 
    end function

    subroutine  calculate_Eccentricity_kernel_by_density_matrix_element
        !-------------------------------------------------------------------------------------------
        !   <J   K_1 q_1 Pi  | E_n |J   K_2 q_2 Pi>
        ! = \sum_{m1 m2} e_{m1 m2} \rho_{m2 m1}
        !  + \sum_{m1 m2 m3 m4} e_{m1 m2 m3 m4} \rho_{m4 m3 m2 m1}
        ! where e_{m1 m2} = \sum_m f^n_{m1 m} f^{-n}_{m m2}
        !       e_{m1 m2 m3 m4} = - f^n_{m1 m3}f^{-n}_{m2 m4}
        !       \rho_{m2 m1} = <q1|c^+_{m2}c_{m1} P^{J}_{K1 K2} P^N P^Z P^{Pi}|q2>
        !      \rho_{m4 m3 m2 m1} = <q1|c^+_{m1}c^+_{m2}c_{m3}c_{m4} P^{J}_{K1 K2}P^N P^Z P^{Pi}|q2>
        !-------------------------------------------------------------------------------------------
        use Globals, only: gcm_space,Proj_option,BS,kernels
        use Proj_Density, only: calculate_one_body_density_matrix_element,calculate_two_body_density_matrix_element
        integer :: J,K1_start,K1_end,K2_start,K2_end,K1,K2,iParity,Parity,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,total_iter,iter,mu
        complex(r64) :: ME1B(2),ME2B(2),Eccentricity_1B(2),Eccentricity_2B(2),Q_nm(-2:2,2),Eccentricity_NP
        real(r64) :: e_1B,e_2B 
        write(*,'(5x,A)') 'calculate_Eccentricity_kernel_by_density_matrix_element ... '
        do J = 0,0
            if(Proj_option%AMPtype==0 .or. Proj_option%AMPtype==1) then
                K1_start = 0
                K1_end = 0
                K2_start = 0
                K2_end = 0
            else
                K1_start = -J
                K1_end = J
                K2_start = -J
                K2_end = J
            end if
            do K1 = K1_start, K1_end
                do K2 = K2_start, K2_end
                    do iParity = 1, 2
                        Parity = (-1)**(iParity+1) ! 1: +1 , 2: -1
                        if( Parity /= (-1)**J) cycle
                        Eccentricity_1B = (0.d0,0.d0)
                        Eccentricity_2B = (0.d0,0.d0)
                        Q_nm = (0.d0,0.d0)
                        Eccentricity_NP = (0.d0,0.d0)
                        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
                        !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,mu,e_1B,ME1B,e_2B,ME2B) &
                        !$omp reduction(+:Eccentricity_1B,Eccentricity_2B,Q_nm)
                        !$OMP DO COLLAPSE(1) SCHEDULE(static)
                        do iter = 1,total_iter
                            if(iter <= BS%HO_sph%idsp(1,1)) then
                                ifg1 = 1
                                m1 = iter
                            else
                                ifg1 = 2
                                m1 = iter - BS%HO_sph%idsp(1,1)
                            end if
                            do ifg2 = 1, 2
                            do m2 = 1, BS%HO_sph%idsp(1,ifg2)
                                ! 1 Body
                                if(ifg1==ifg2) then
                                    call eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,2,e_1B) ! e_{m1 m2}

                                    call calculate_one_body_density_matrix_element(J,K1,K2,Parity,ifg2,m2,ifg1,m1,ME1B) ! rho_{m2 m1}
                                    Eccentricity_1B(1) = Eccentricity_1B(1) +  e_1B*ME1B(1)
                                    Eccentricity_1B(2) = Eccentricity_1B(2) +  e_1B*ME1B(2)
                                    do mu = -2, 2 
                                        ! neutron part 
                                        Q_nm(mu,1) = Q_nm(mu,1) +  Q(ifg1,m1,ifg2,m2,2,mu)*ME1B(1)
                                        Q_nm(mu,2) = Q_nm(mu,2) +  Q(ifg1,m1,ifg2,m2,2,mu)*ME1B(2)
                                    end do 
                                end if  
                                do ifg3 = 1, 2
                                do m3 = 1, BS%HO_sph%idsp(1,ifg3)
                                    do ifg4 = 1, 2
                                    do m4 = 1, BS%HO_sph%idsp(1,ifg4)
                                        ! 2 Body
                                        if(ifg1==ifg3 .and. ifg2==ifg4) then
                                            call eccentricity_matrix_element_two_body(ifg1,m1,ifg2,m2,ifg4,m4,ifg3,m3,2,e_2B) ! e_{m1 m2 m4 m3}: - f^n_{m1 m3}f^{-n}_{m2 m4}
                                            if(abs(e_2B) > 1.E-8) then ! It can save a lot of time !
                                                call calculate_two_body_density_matrix_element(J,K1,K2,Parity,ifg4,m4,ifg3,m3,ifg2,m2,ifg1,m1,ME2B) !  rho_{m4 m3 m2 m1}: <c^+_{m1}c^+_{m2}c_{m3}c_{m4}>
                                                Eccentricity_2B(1) = Eccentricity_2B(1) +  e_2B*ME2B(1) ! e_{m1 m2 m4 m3} * <c+_m1 c+_m2 c_m3 c_m4> = e_{m1 m2 m3 m4} * rho_{m4 m3 m2 m1}
                                                Eccentricity_2B(2) = Eccentricity_2B(2) +  e_2B*ME2B(2)
                                            end if 
                                        end if 
                                    end do 
                                    end do
                                end do
                                end do 
                            end do
                            end do 
                        end do

                        !$OMP END PARALLEL
                        kernels%Eccentricity_KK_byDensity(J,K1,K2,1,iParity,1) = Eccentricity_1B(1) 
                        kernels%Eccentricity_KK_byDensity(J,K1,K2,1,iParity,2) = Eccentricity_2B(1)
                        kernels%Eccentricity_KK_byDensity(J,K1,K2,2,iParity,1) = Eccentricity_1B(2) 
                        kernels%Eccentricity_KK_byDensity(J,K1,K2,2,iParity,2) = Eccentricity_2B(2)
                        do mu = -2, 2
                             Eccentricity_NP = Eccentricity_NP +  (-1)**mu*Q_nm(mu,1)*Q_nm(-mu,2) + (-1)**mu*Q_nm(mu,2)*Q_nm(-mu,1)
                        end do
                        kernels%Eccentricity_KK_byDensity(J,K1,K2,1,iParity,3) = Eccentricity_NP
                    end do 
                end do 
            end do 
        end do
    end subroutine

end Module Eccentricity