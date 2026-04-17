Module Eccentricity
    use Constants,only: r64
    implicit none
    integer, parameter :: max_n = 2 ! max n of f_n
    real(r64), dimension(:,:,:,:,:), allocatable:: f_n_array ! f_n_array(ifg1, m1, ifg2, m2, n)
    real(r64), dimension(:,:,:,:,:,:), allocatable:: Q_array ! Q_array(ifg1,m1,ifg2,m2,n,mu)
    real(r64), dimension(:,:,:,:,:), allocatable :: Ecc_ME1B_array ! Ecc_ME1B_array(ifg1,m1,ifg2,m2,n)
    real(r64), dimension(:,:,:,:,:,:,:), allocatable :: Ecc_ME2B_array ! Ecc_ME2B_array(ifg1,ifg2,m1,m2,m3,m4,n)
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
        logical :: all_zero
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
                call precompute_eccentricity_matrix_element_two_body(n)
                precomputed_n(n) = .True.
            end if
        end if

        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
        !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2,mu,ifg3,m3,ifg4,m4,all_zero,e_1B,e_2B) &
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
            if(ifg1 /= ifg2) cycle
            do m2 = 1,BS%HO_sph%idsp(1,ifg2)
                ! 1 Body
                call eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,n,e_1B)
                Eccentri(1) = Eccentri(1) + e_1B*mix%rho_mm(m2,m1,ifg1,iphi,it)
                pEccentri(1) = pEccentri(1) + e_1B*mix%prho_mm(m2,m1,ifg1,iphi,it)

                ! check if all Q_m1m2(n,mu) are zero
                all_zero = .true.
                do mu = -n, n
                    if (Q(ifg1,m1,ifg2,m2,n,mu) /= 0.0d0) then
                        all_zero = .false.
                        exit
                    end if
                end do
                if (all_zero) cycle

                ! to calculate nn-pp and pp-nn term
                do mu = -n, n 
                    Q_mu(mu) = Q_mu(mu) +  Q(ifg1,m1,ifg2,m2,n,mu)*mix%rho_mm(m2,m1,ifg1,iphi,it)
                    pQ_mu(mu) = pQ_mu(mu) +  Q(ifg1,m1,ifg2,m2,n,mu)*mix%prho_mm(m2,m1,ifg1,iphi,it)
                end do

                do ifg3 = 1,2
                do m3 = 1,BS%HO_sph%idsp(1,ifg3)
                    do ifg4 = 1,2
                    if(ifg3 /= ifg4) cycle
                    do m4 = 1,BS%HO_sph%idsp(1,ifg4)
                        ! 2 Body
                        ! call eccentricity_matrix_element_two_body(ifg1,m1,ifg3,m3,ifg4,m4,ifg2,m2,n,e_2B)
                        e_2B = Ecc_ME2B_array(ifg1,ifg3,m1,m3,m4,m2,n)
                        if(e_2B == 0.d0) cycle
                        !  e_{m1 m3 m4 m2}*<c^+_{m1}c^+_{m3}c_{m2}c_{m4}> = e_{m1 m3 m4 m2}*rho_4231
                        Eccentri(2) = Eccentri(2)+ e_2B*( & 
                                        mix%rho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%rho_mm(m2,m3,indexfg(ifg2,ifg3),iphi,it) &
                                    - mix%rho_mm(m2,m1,indexfg(ifg2,ifg1),iphi,it)*mix%rho_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it) &
                                    + mix%kappa01c_mm(m1,m3,indexfg(ifg1,ifg3),iphi,it)*mix%kappa10_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it))
                        pEccentri(2) = pEccentri(2) + e_2B*( &
                                        mix%prho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%prho_mm(m2,m3,indexfg(ifg2,ifg3),iphi,it) &
                                    - mix%prho_mm(m2,m1,indexfg(ifg2,ifg1),iphi,it)*mix%prho_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it) &
                                    + mix%pkappa01c_mm(m1,m3,indexfg(ifg1,ifg3),iphi,it)*mix%pkappa10_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it))
                    end do
                    end do
                end do 
                end do 
            end do 
            end do 
        end do 
        !$OMP END PARALLEL
    end subroutine

    subroutine calculate_Eccentri_n_MM(n,iphi,it,Eccentri,pEccentri,Q_mu,pQ_mu,Each_2B_Term,pEach_2B_Term)
        !--------------------------------------------------
        !  calculate <q1| E_n R|q2>/<q1|R|q2>
        !  where R = R(alpha,beta,gamma, phi_n,phi_p) = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
        !--------------------------------------------------
        use Constants,only: itx
        use Globals, only: BS, mix
        use MathMethods, only: zGEMM,zTrace,zGEMM_Trace
        integer :: n,iphi,it
        complex(r64) :: Eccentri(2),pEccentri(2),Q_mu(-n:n),pQ_mu(-n:n),Each_2B_Term(3),pEach_2B_Term(3)
        integer :: ifg1,ifg2,mu,dim1,dim2,dim_m_max
        complex(r64), allocatable,dimension(:,:) :: QQ12,Q12,Q34,Matrix1,Matrix2,Matrix3, &
                rho21,rho43,rho41,rho23,kappa42,kappac13,prho21,prho43,prho41,prho23,pkappa42,pkappac13
        complex(r64) :: zero,one 

        zero = (0.0D0, 0.0D0)
        one = (1.0D0, 0.0D0)
        Eccentri = (0.d0,0.d0)
        pEccentri = (0.d0,0.d0)
        Q_mu = (0.d0,0.d0)
        pQ_mu = (0.d0,0.d0)
        Each_2B_Term = (0.d0,0.d0)
        pEach_2B_Term = (0.d0,0.d0)
        
        if(.not. precomputed_n(n)) then 
            if(n>max_n) then 
                write(*,*) 'Warning: max_n is small than input n ! f_n will not be precomputed.'
            else 
                call precompute_f_n(n)
                call precompute_f_n(-n)
                call precompute_Q(n)
                call precompute_eccentricity_matrix_element_one_body(n)
                ! call precompute_eccentricity_matrix_element_two_body(n)
                precomputed_n(n) = .True.
            end if
        end if

        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        allocate(QQ12(dim_m_max,dim_m_max),Q12(dim_m_max,dim_m_max),Q34(dim_m_max,dim_m_max), &
                 Matrix1(dim_m_max,dim_m_max),Matrix2(dim_m_max,dim_m_max),Matrix3(dim_m_max,dim_m_max), &
                 rho21(dim_m_max,dim_m_max),rho43(dim_m_max,dim_m_max),rho41(dim_m_max,dim_m_max), &
                 rho23(dim_m_max,dim_m_max),kappa42(dim_m_max,dim_m_max),kappac13(dim_m_max,dim_m_max), &
                 prho21(dim_m_max,dim_m_max),prho43(dim_m_max,dim_m_max),prho41(dim_m_max,dim_m_max), &
                 prho23(dim_m_max,dim_m_max),pkappa42(dim_m_max,dim_m_max),pkappac13(dim_m_max,dim_m_max), &
                 source=(0.d0,0.d0))
        do ifg1 = 1,2
            dim1 = BS%HO_sph%idsp(1,ifg1)
            ! 1 Body
            QQ12(1:dim1,1:dim1) = Ecc_ME1B_array(ifg1,1:dim1,ifg1,1:dim1,n)
            rho21(1:dim1,1:dim1) = mix%rho_mm(1:dim1,1:dim1,indexfg(ifg1,ifg1),iphi,it)
            prho21(1:dim1,1:dim1) = mix%prho_mm(1:dim1,1:dim1,indexfg(ifg1,ifg1),iphi,it)
            call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,QQ12,dim_m_max,rho21,dim_m_max,zero,Matrix1,dim_m_max)
            Eccentri(1) = Eccentri(1) + zTrace(Matrix1,dim1)
            call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,QQ12,dim_m_max,prho21,dim_m_max,zero,Matrix1,dim_m_max)
            pEccentri(1) = pEccentri(1) + zTrace(Matrix1,dim1)
            do ifg2 = 1,2
                dim2 = BS%HO_sph%idsp(1,ifg2)
                rho43(1:dim2,1:dim2) = mix%rho_mm(1:dim2,1:dim2,indexfg(ifg2,ifg2),iphi,it)
                rho41(1:dim2,1:dim1) = mix%rho_mm(1:dim2,1:dim1,indexfg(ifg2,ifg1),iphi,it)
                rho23(1:dim1,1:dim2) = mix%rho_mm(1:dim1,1:dim2,indexfg(ifg1,ifg2),iphi,it)
                kappa42(1:dim2,1:dim1) = mix%kappa10_mm(1:dim2,1:dim1,indexfg(ifg2,ifg1),iphi,it)
                kappac13(1:dim1,1:dim2) = mix%kappa01c_mm(1:dim1,1:dim2,indexfg(ifg1,ifg2),iphi,it)
                prho43(1:dim2,1:dim2) = mix%prho_mm(1:dim2,1:dim2,indexfg(ifg2,ifg2),iphi,it)
                prho41(1:dim2,1:dim1) = mix%prho_mm(1:dim2,1:dim1,indexfg(ifg2,ifg1),iphi,it)
                prho23(1:dim1,1:dim2) = mix%prho_mm(1:dim1,1:dim2,indexfg(ifg1,ifg2),iphi,it)
                pkappa42(1:dim2,1:dim1) = mix%pkappa10_mm(1:dim2,1:dim1,indexfg(ifg2,ifg1),iphi,it)
                pkappac13(1:dim1,1:dim2) = mix%pkappa01c_mm(1:dim1,1:dim2,indexfg(ifg1,ifg2),iphi,it)
                ! 2B
                do mu = -n, n
                    ! direct term(nn or pp)
                    Q12(1:dim1,1:dim1) =  Q_array(ifg1,1:dim1,ifg1,1:dim1,n,mu)
                    Q34(1:dim2,1:dim2) =  Q_array(ifg2,1:dim2,ifg2,1:dim2,n,-mu)
                    call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,Q12,dim_m_max,rho21,dim_m_max,zero,Matrix1,dim_m_max)
                    call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim2,one,Q34,dim_m_max,rho43,dim_m_max,zero,Matrix2,dim_m_max)
                    if(ifg2==1) Q_mu(mu) = Q_mu(mu) + zTrace(Matrix1,dim1)
                    Eccentri(2) = Eccentri(2) - (-1)**mu*zTrace(Matrix1,dim1)*zTrace(Matrix2,dim2)
                    Each_2B_Term(1) = Each_2B_Term(1) - (-1)**mu*zTrace(Matrix1,dim1)*zTrace(Matrix2,dim2)
                    call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,Q12,dim_m_max,prho21,dim_m_max,zero,Matrix1,dim_m_max)
                    call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim2,one,Q34,dim_m_max,prho43,dim_m_max,zero,Matrix2,dim_m_max)
                    if(ifg2==1) pQ_mu(mu) = pQ_mu(mu) + zTrace(Matrix1,dim1)
                    pEccentri(2) = pEccentri(2) - (-1)**mu*zTrace(Matrix1,dim1)*zTrace(Matrix2,dim2)
                    pEach_2B_Term(1) = pEach_2B_Term(1) - (-1)**mu*zTrace(Matrix1,dim1)*zTrace(Matrix2,dim2)

                    ! exchange term
                    call zGEMM('N' ,'N' ,dim1,dim2,dim1,one,Q12,dim_m_max,rho23,dim_m_max,zero,Matrix1,dim_m_max)
                    call zGEMM('N' ,'N' ,dim2,dim1,dim2,one,Q34,dim_m_max,rho41,dim_m_max,zero,Matrix2,dim_m_max)
                    call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim2,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                    Eccentri(2) = Eccentri(2) + (-1)**mu*zTrace(Matrix3,dim1)
                    Each_2B_Term(2) = Each_2B_Term(2) + (-1)**mu*zTrace(Matrix3,dim1)
                    call zGEMM('N' ,'N' ,dim1,dim2,dim1,one,Q12,dim_m_max,prho23,dim_m_max,zero,Matrix1,dim_m_max)
                    call zGEMM('N' ,'N' ,dim2,dim1,dim2,one,Q34,dim_m_max,prho41,dim_m_max,zero,Matrix2,dim_m_max)
                    call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim2,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                    pEccentri(2) = pEccentri(2) + (-1)**mu*zTrace(Matrix3,dim1)
                    pEach_2B_Term(2) = pEach_2B_Term(2) + (-1)**mu*zTrace(Matrix3,dim1)

                    ! kappa term
                    call zGEMM('N' ,'T' ,dim2,dim1,dim1,one,kappa42,dim_m_max,Q12,dim_m_max,zero,Matrix1,dim_m_max)
                    call zGEMM('N' ,'N' ,dim1,dim2,dim2,one,kappac13,dim_m_max,Q34,dim_m_max,zero,Matrix2,dim_m_max)
                    call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim1,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                    Eccentri(2) = Eccentri(2) + (-1)**mu*zTrace(Matrix3,dim2)
                    Each_2B_Term(3) = Each_2B_Term(3) + (-1)**mu*zTrace(Matrix3,dim2)
                    call zGEMM('N' ,'T' ,dim2,dim1,dim1,one,pkappa42,dim_m_max,Q12,dim_m_max,zero,Matrix1,dim_m_max)
                    call zGEMM('N' ,'N' ,dim1,dim2,dim2,one,pkappac13,dim_m_max,Q34,dim_m_max,zero,Matrix2,dim_m_max)
                    call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim1,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                    pEccentri(2) = pEccentri(2) + (-1)**mu*zTrace(Matrix3,dim2)
                    pEach_2B_Term(3) = pEach_2B_Term(3) + (-1)**mu*zTrace(Matrix3,dim2)
                end do 
            end do 
        end do
        deallocate(QQ12,Q12,Q34,Matrix1,Matrix2,Matrix3,rho21,rho43,rho41,rho23,kappa42,kappac13,prho21,prho43,prho41,prho23,pkappa42,pkappac13)
    end subroutine

    ! --------      -----------
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

    subroutine eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,n,EME1B)
        !-----------------------------------------------------------
        !   e^{n}_{m1 m2} = \sum_m f^n_{m1 m} f^{-n}_{m m2} 
        ! 
        !   e^{n}_{p q}  = \sum_{m p'} (-1)^m (Q_nm)_{pp'}(Q_n-m)_{p'q}  
        !-----------------------------------------------------------
        use Globals, only: BS
        integer,intent(in) :: ifg1,m1,ifg2,m2,n
        real(r64) :: EME1B
        if(precomputed_n(n)) then 
            EME1B = Ecc_ME1B_array(ifg1,m1,ifg2,m2,n)
        else
            EME1B = compute_eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,n)
        end if 
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
        integer :: ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,n
        real(r64) :: EME2B
        if(precomputed_n(n)) then 
            EME2B = Ecc_ME2B_array(ifg1,ifg2,m1,m2,m3,m4,n)
        else
            EME2B = compute_eccentricity_matrix_element_two_body(ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,n)
        end if
    end subroutine

    ! ---- precalculate and store ------
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

    subroutine precompute_eccentricity_matrix_element_one_body(n)
        use Globals, only: BS
        integer,intent(in) :: n
        integer :: ifg1,m1,ifg2,m2,mu,dim_m_max,iter,total_iter
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        if(.not. allocated(Ecc_ME1B_array)) allocate(Ecc_ME1B_array(2,dim_m_max,2,dim_m_max,max_n),source=0.d0)

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
                        Ecc_ME1B_array(ifg1,m1,ifg2,m2,n) = compute_eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,n)
                    end do
                end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine precompute_eccentricity_matrix_element_two_body(n)
        use Globals, only: BS
        integer,intent(in) :: n
        integer :: ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,mu,dim_m_max,iter,total_iter
        logical :: all_zero
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        if(.not. allocated(Ecc_ME2B_array)) allocate(Ecc_ME2B_array(2,2,dim_m_max,dim_m_max,dim_m_max,dim_m_max,max_n),source=0.d0)

        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
        !$omp parallel do collapse(1) default(shared) &
        !$omp private(iter,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,mu,all_zero) schedule(static)
        do iter = 1,total_iter
            if(iter <= BS%HO_sph%idsp(1,1)) then
                ifg1 = 1
                m1 = iter
            else
                ifg1 = 2
                m1 = iter - BS%HO_sph%idsp(1,1)
            end if
            do ifg4 = 1,2
            if(ifg1/=ifg4) cycle
            do m4 = 1, BS%HO_sph%idsp(1,ifg4)
                    all_zero = .true.
                    do mu = -n, n
                        if (Q(ifg1, m1, ifg4, m4, n, mu) /= 0.0d0) then
                            all_zero = .false.
                            exit
                        end if
                    end do
                    if (all_zero) cycle

                    do ifg2 = 1,2
                    do m2 = 1,BS%HO_sph%idsp(1,ifg2)
                        do ifg3 = 1,2
                        if(ifg2/=ifg3) cycle
                        do m3 = 1,BS%HO_sph%idsp(1,ifg3)
                            Ecc_ME2B_array(ifg1,ifg2,m1,m2,m3,m4,n) = compute_eccentricity_matrix_element_two_body(ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,n)
                        end do
                        end do
                    end do
                    end do
            end do
            end do
        end do
        !$omp end parallel do
    end subroutine

    ! ------ function to calculate -----
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

    double precision function compute_eccentricity_matrix_element_one_body(ifg1,m1,ifg2,m2,n)
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
        compute_eccentricity_matrix_element_one_body = EME1B
    end function

    double precision function compute_eccentricity_matrix_element_two_body(ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,n)
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
        compute_eccentricity_matrix_element_two_body = EME2B
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


    ! ------------------------------------------------------------
    !--------------------------------------------------------------
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
        complex(r64) :: ME1B(2),ME2B(4),Eccentricity_1B(2),Eccentricity_2B(2),Eccentricity_NP(2)
        real(r64) :: e_1B,e_2B,QQ
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
                        Eccentricity_NP = (0.d0,0.d0)
                        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
                        !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,mu,e_1B,ME1B,e_2B,QQ,ME2B) &
                        !$omp reduction(+:Eccentricity_1B,Eccentricity_2B,Eccentricity_NP)
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
                                    if(abs(e_1B) > 1.E-8) then
                                        call calculate_one_body_density_matrix_element(J,K1,K2,Parity,ifg2,m2,ifg1,m1,ME1B) ! rho_{m2 m1}
                                        Eccentricity_1B(1) = Eccentricity_1B(1) +  e_1B*ME1B(1)
                                        Eccentricity_1B(2) = Eccentricity_1B(2) +  e_1B*ME1B(2)
                                    end if 
                                end if  
                                do ifg3 = 1, 2
                                do m3 = 1, BS%HO_sph%idsp(1,ifg3)
                                    do ifg4 = 1, 2
                                    do m4 = 1, BS%HO_sph%idsp(1,ifg4)
                                        ! 2 Body
                                        if(ifg1==ifg3 .and. ifg2==ifg4) then
                                            call eccentricity_matrix_element_two_body(ifg1,m1,ifg2,m2,ifg4,m4,ifg3,m3,2,e_2B) ! e_{m1 m2 m4 m3}: - f^n_{m1 m3}f^{-n}_{m2 m4}
                                            QQ = (0.d0,0.d0)
                                            do mu = -2, 2 
                                                QQ = QQ - (-1)**mu*Q(ifg1,m1,ifg3,m3,2,mu)*Q(ifg2,m2,ifg4,m4,2,-mu) ! np and pn part
                                            end do
                                            if(abs(e_2B)> 1.E-8 .or. abs(QQ) > 1.E-8) then ! It can save a lot of time !
                                                call calculate_two_body_density_matrix_element(J,K1,K2,Parity,ifg4,m4,ifg3,m3,ifg2,m2,ifg1,m1,ME2B) !  rho_{m4 m3 m2 m1}: <c^+_{m1}c^+_{m2}c_{m3}c_{m4}>
                                                Eccentricity_2B(1) = Eccentricity_2B(1) +  e_2B*ME2B(1) ! e_{m1 m2 m4 m3} * <c+_m1 c+_m2 c_m3 c_m4> = e_{m1 m2 m3 m4} * rho_{m4 m3 m2 m1}
                                                Eccentricity_2B(2) = Eccentricity_2B(2) +  e_2B*ME2B(2)
                                                Eccentricity_NP(1) = Eccentricity_NP(1) +  QQ*ME2B(3) ! np
                                                Eccentricity_NP(2) = Eccentricity_NP(2) +  QQ*ME2B(4) ! pn
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
                        kernels%Eccentricity_KK_byDensity(J,K1,K2,3,iParity,1) = Eccentricity_NP(1)
                        kernels%Eccentricity_KK_byDensity(J,K1,K2,3,iParity,2) = Eccentricity_NP(2)
                        
                    end do 
                end do 
            end do 
        end do
    end subroutine

end Module Eccentricity