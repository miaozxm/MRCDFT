!==============================================================================!
! MODULE JNsquare                                                              !
!                                                                              !
! This module calculates the ...                                               !
!                                                                              !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
Module JNsquare
    use Constants,only: r64
    implicit none
    complex(r64), dimension(:,:,:,:,:), allocatable:: angular_matrix_element_array ! angular_matrix_element_array(Jx/Jy/Jz,ifg1, m1, ifg2, m2)
    complex(r64), dimension(:,:,:,:), allocatable :: JJ_ME1B_array ! JJ_ME1B_array(ifg1,m1,ifg2,m2) = \sum_i <n1 l1 j1 m1|J_i^2|n2 l2 j2 m2>
    logical :: precomputed_angular_matrix_element = .False.
contains

! ------------------------------------------------------------------------
!                                J^2
! -------------------------------------------------------------------------
subroutine calculate_Jsquare_and_J(iphi,it,J2,pJ2,J_i,pJ_i)
    !---------------------------------------------------------------------------
    ! 
    !  calculate  
    !  1)   <q1| J^2 R |q2>/<q1|R|q2>  = \sum_i <q1| J_i^2 R |q2>/<q1|R|q2> 
    !  2)   <q1| J_i R |q2>/<q1|R|q2> 
    !  for neutron and porton where i= x,y,z ;  
    !  and R = R(alpha,beta,gamma, phi_n,phi_p) 
    !        = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
    ! 
    ! -------------------------------------------------------------------------------
    !   <q1| J_i^2 R |q2>/<q1|R|q2>  
    !    =  \sum_{m1 m2 m3 m4}  (J_i)_{m1 m2} (J_i)_{m2 m4} rho_{m4 m1} \delta_{m2 m3}
    !    -  \sum_{m1 m2 m3 m4}  (J_i)_{m1 m3} (J_i)_{m2 m4} rho_{m4 m3 m2 m1}
    !   where rho_{m3 m1} = <q1|c^+_m3 c_m1 R|q2>/<q1|R|q2>
    !         rho_{m4 m3 m2 m1} = <q1|c^+_m1 c^+_m2 c_m3 c_m4 R|q2>/<q1|R|q2>
    !   and rho_{m4 m3 m2 m1} = rho_{m4 m1}rho_{m3 m2} - rho_{m3 m1}rho_{m4 m2} + kappa01*(m1 m2)kappa10(m4 m3)
    !---------------------------------------------------------------------------
    use Globals, only: BS, mix
    integer :: iphi,it
    complex(r64) :: J2,pJ2,J12_J23,rho_4321,prho_4321,J13_J24,J_i(3),pJ_i(3)
    integer :: ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,total_iter,iter

    J2 = (0.0d0,0.0d0)
    pJ2 = (0.0d0,0.0d0)
    J_i = (0.0d0,0.0d0)
    pJ_i = (0.0d0,0.0d0)

    if(.not. precomputed_angular_matrix_element) then
        call precompute_angular_matrix_element 
    end if

    total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,J12_J23,rho_4321,prho_4321,J13_J24) &
    !$omp reduction(+:J2,pJ2,J_i,pJ_i)
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
        do m2 = 1, BS%HO_sph%idsp(1,ifg2)

            ! calculate <q1| J_i R |q2>/<q1|R|q2> 
            if(ifg1==ifg2) then
                J_i(1) = J_i(1) + angular_matrix_element(1,ifg1,m1,ifg2,m2)*mix%rho_mm(m2,m1,ifg1,iphi,it)
                J_i(2) = J_i(2) + angular_matrix_element(2,ifg1,m1,ifg2,m2)*mix%rho_mm(m2,m1,ifg1,iphi,it)
                J_i(3) = J_i(3) + angular_matrix_element(3,ifg1,m1,ifg2,m2)*mix%rho_mm(m2,m1,ifg1,iphi,it)
                pJ_i(1) = pJ_i(1) + angular_matrix_element(1,ifg1,m1,ifg2,m2)*mix%prho_mm(m2,m1,ifg1,iphi,it)
                pJ_i(2) = pJ_i(2) + angular_matrix_element(2,ifg1,m1,ifg2,m2)*mix%prho_mm(m2,m1,ifg1,iphi,it)
                pJ_i(3) = pJ_i(3) + angular_matrix_element(3,ifg1,m1,ifg2,m2)*mix%prho_mm(m2,m1,ifg1,iphi,it)
            end if 

            ! calculate \sum_i <q1| J_i^2 R |q2>/<q1|R|q2>
            do ifg3 = 1,2
            do m3 = 1, BS%HO_sph%idsp(1,ifg3)
                ! 1B
                if(ifg1==ifg3.and. ifg1 == ifg2) then
                    J12_J23 = angular_matrix_element(1,ifg1,m1,ifg2,m2)*angular_matrix_element(1,ifg2,m2,ifg3,m3) &
                            + angular_matrix_element(2,ifg1,m1,ifg2,m2)*angular_matrix_element(2,ifg2,m2,ifg3,m3) &
                            + angular_matrix_element(3,ifg1,m1,ifg2,m2)*angular_matrix_element(3,ifg2,m2,ifg3,m3)
                    J2 = J2 + J12_J23 * mix%rho_mm(m3,m1,ifg1,iphi,it)
                    pJ2 = pJ2 + J12_J23 * mix%prho_mm(m3,m1,ifg1,iphi,it)
                end if

                do ifg4 = 1,2
                do m4 = 1, BS%HO_sph%idsp(1,ifg4) 
                    if(ifg1==ifg3 .and. ifg2==ifg4) then
                        ! nonzero (J_i)_{m2 m4} 
                        if( angular_matrix_element(1,ifg2,m2,ifg4,m4)==(0.d0,0.d0) .and. angular_matrix_element(2,ifg2,m2,ifg4,m4)==(0.d0,0.d0) & 
                            .and. angular_matrix_element(3,ifg2,m2,ifg4,m4)==(0.d0,0.d0) ) cycle 
                        ! 2B
                        J13_J24 = angular_matrix_element(1,ifg1,m1,ifg3,m3)*angular_matrix_element(1,ifg2,m2,ifg4,m4)  &
                                    + angular_matrix_element(2,ifg1,m1,ifg3,m3)*angular_matrix_element(2,ifg2,m2,ifg4,m4)  &
                                    + angular_matrix_element(3,ifg1,m1,ifg3,m3)*angular_matrix_element(3,ifg2,m2,ifg4,m4)
                        rho_4321  = mix%rho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%rho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                    - mix%rho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%rho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                    + mix%kappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%kappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it)
                        prho_4321 =  mix%prho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%prho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                    - mix%prho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%prho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                    + mix%pkappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%pkappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it)
                        J2 = J2 - J13_J24*rho_4321
                        pJ2 = pJ2 - J13_J24*prho_4321
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

function angular_matrix_element(case,ifg1,m1,ifg2,m2)
    complex(r64) :: angular_matrix_element
    integer :: case,ifg1,m1,ifg2,m2
    if(precomputed_angular_matrix_element) then
        angular_matrix_element = angular_matrix_element_array(case,ifg1,m1,ifg2,m2)
    else 
        angular_matrix_element = compute_angular_matrix_element(case,ifg1,m1,ifg2,m2)
    end if 
end

subroutine calculate_Jsquare_and_J_MM(iphi,it,J2,pJ2,J_i,pJ_i)
    !---------------------------------------------------------------------------
    ! 
    !  calculate  
    !  1)   <q1| J^2 R |q2>/<q1|R|q2>  = \sum_i <q1| J_i^2 R |q2>/<q1|R|q2> 
    !  2)   <q1| J_i R |q2>/<q1|R|q2> 
    !  for neutron and porton where i= x,y,z ;  
    !  and R = R(alpha,beta,gamma, phi_n,phi_p) 
    !        = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
    ! 
    ! -------------------------------------------------------------------------------
    !   <q1| J_i^2 R |q2>/<q1|R|q2>  
    !    =  \sum_{m1 m2 m3 m4}  (J_i)_{m1 m2} (J_i)_{m2 m4} rho_{m4 m1} \delta_{m2 m3}
    !    -  \sum_{m1 m2 m3 m4}  (J_i)_{m1 m2} (J_i)_{m3 m4} rho_{m4 m2 m3 m1}
    !   where rho_{m4 m1} = <q1|c^+_m4 c_m1 R|q2>/<q1|R|q2>
    !         rho_{m4 m2 m3 m1} = <q1|c^+_m1 c^+_m3 c_m2 c_m4 R|q2>/<q1|R|q2>
    !   and rho_{m4 m2 m3 m1} = rho_{m4 m1}rho_{m2 m3} - rho_{m2 m1}rho_{m4 m3} + kappa01*(m1 m3)kappa10(m4 m2)
    !---------------------------------------------------------------------------
    use Globals, only: BS, mix
    use MathMethods, only: zGEMM,zTrace,zGEMM_Trace
    integer :: iphi,it
    complex(r64) :: J2,pJ2,J_i(3),pJ_i(3)
    integer :: ifg1,ifg2,i,dim1,dim2,dim_m_max
    complex(r64), allocatable,dimension(:,:) :: JJ12,J12,J34,Matrix1,Matrix2,Matrix3, &
            rho21,rho43,rho41,rho23,kappa42,kappac13,prho21,prho43,prho41,prho23,pkappa42,pkappac13
    complex(r64) :: zero,one 
    zero = (0.0D0, 0.0D0)
    one = (1.0D0, 0.0D0)
    J2 = zero
    pJ2 = zero
    J_i = zero
    pJ_i = zero
    if(.not. precomputed_angular_matrix_element) then
        call precompute_angular_matrix_element 
        call precompute_JJ_matrix_element_one_body
    end if

    dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))

    allocate(JJ12(dim_m_max,dim_m_max),J12(dim_m_max,dim_m_max),J34(dim_m_max,dim_m_max), &
             Matrix1(dim_m_max,dim_m_max),Matrix2(dim_m_max,dim_m_max),Matrix3(dim_m_max,dim_m_max), &
             rho21(dim_m_max,dim_m_max),rho43(dim_m_max,dim_m_max),rho41(dim_m_max,dim_m_max), &
             rho23(dim_m_max,dim_m_max),kappa42(dim_m_max,dim_m_max),kappac13(dim_m_max,dim_m_max), &
             prho21(dim_m_max,dim_m_max),prho43(dim_m_max,dim_m_max),prho41(dim_m_max,dim_m_max), &
             prho23(dim_m_max,dim_m_max),pkappa42(dim_m_max,dim_m_max),pkappac13(dim_m_max,dim_m_max), &
             source=(0.d0,0.d0))
    do ifg1 = 1,2
        dim1 = BS%HO_sph%idsp(1,ifg1)
        ! 1 Body
        JJ12(1:dim1,1:dim1) = JJ_ME1B_array(ifg1,1:dim1,ifg1,1:dim1)
        rho21(1:dim1,1:dim1) = mix%rho_mm(1:dim1,1:dim1,indexfg(ifg1,ifg1),iphi,it)
        prho21(1:dim1,1:dim1) = mix%prho_mm(1:dim1,1:dim1,indexfg(ifg1,ifg1),iphi,it)
        call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,JJ12,dim_m_max,rho21,dim_m_max,zero,Matrix1,dim_m_max)
        J2 = J2 + zTrace(Matrix1,dim1)
        call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,JJ12,dim_m_max,prho21,dim_m_max,zero,Matrix1,dim_m_max)
        pJ2 = pJ2 + zTrace(Matrix1,dim1)

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
            do i = 1, 3
                ! direct term(nn or pp)
                J12(1:dim1,1:dim1) =  angular_matrix_element_array(i,ifg1,1:dim1,ifg1,1:dim1)
                J34(1:dim2,1:dim2) =  angular_matrix_element_array(i,ifg2,1:dim2,ifg2,1:dim2)
                call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,J12,dim_m_max,rho21,dim_m_max,zero,Matrix1,dim_m_max)
                call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim2,one,J34,dim_m_max,rho43,dim_m_max,zero,Matrix2,dim_m_max)
                if(ifg2==1) J_i(i) = J_i(i) + zTrace(Matrix1,dim1)
                J2 = J2 + zTrace(Matrix1,dim1)*zTrace(Matrix2,dim2)
                call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim1,one,J12,dim_m_max,prho21,dim_m_max,zero,Matrix1,dim_m_max)
                call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim2,one,J34,dim_m_max,prho43,dim_m_max,zero,Matrix2,dim_m_max)
                if(ifg2==1) pJ_i(i) = pJ_i(i) + zTrace(Matrix1,dim1)
                pJ2 = pJ2 + zTrace(Matrix1,dim1)*zTrace(Matrix2,dim2)

                ! exchange term
                call zGEMM('N' ,'N' ,dim1,dim2,dim1,one,J12,dim_m_max,rho23,dim_m_max,zero,Matrix1,dim_m_max)
                call zGEMM('N' ,'N' ,dim2,dim1,dim2,one,J34,dim_m_max,rho41,dim_m_max,zero,Matrix2,dim_m_max)
                call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim2,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                J2 = J2 - zTrace(Matrix3,dim1)
                call zGEMM('N' ,'N' ,dim1,dim2,dim1,one,J12,dim_m_max,prho23,dim_m_max,zero,Matrix1,dim_m_max)
                call zGEMM('N' ,'N' ,dim2,dim1,dim2,one,J34,dim_m_max,prho41,dim_m_max,zero,Matrix2,dim_m_max)
                call zGEMM_Trace('N' ,'N' ,dim1,dim1,dim2,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                pJ2 = pJ2 - zTrace(Matrix3,dim1)

                ! kappa term
                call zGEMM('N' ,'T' ,dim2,dim1,dim1,one,kappa42,dim_m_max,J12,dim_m_max,zero,Matrix1,dim_m_max)
                call zGEMM('N' ,'N' ,dim1,dim2,dim2,one,kappac13,dim_m_max,J34,dim_m_max,zero,Matrix2,dim_m_max)
                call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim1,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                J2 = J2 - zTrace(Matrix3,dim2)
                call zGEMM('N' ,'T' ,dim2,dim1,dim1,one,pkappa42,dim_m_max,J12,dim_m_max,zero,Matrix1,dim_m_max)
                call zGEMM('N' ,'N' ,dim1,dim2,dim2,one,pkappac13,dim_m_max,J34,dim_m_max,zero,Matrix2,dim_m_max)
                call zGEMM_Trace('N' ,'N' ,dim2,dim2,dim1,one,Matrix1,dim_m_max,Matrix2,dim_m_max,zero,Matrix3,dim_m_max)
                pJ2 = pJ2 - zTrace(Matrix3,dim2)
            end do 
        end do 
    end do
    deallocate(JJ12,J12,J34,Matrix1,Matrix2,Matrix3,rho21,rho43,rho41,rho23,kappa42,kappac13,prho21,prho43,prho41,prho23,pkappa42,pkappac13)
end subroutine

subroutine precompute_angular_matrix_element
    !--------------------------------------------------
    ! store the calculated <n1 l1 j1 m1|J_i|n2 l2 j2 m2> 
    !  to angular_matrix_element_array
    !--------------------------------------------------
    use Globals, only: BS
    integer :: ifg1,m1,ifg2,m2,dim_m_max,iter,total_iter

    dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
    if(.not. allocated(angular_matrix_element_array)) allocate(angular_matrix_element_array(3,2,dim_m_max,2,dim_m_max),source=(0.d0,0.d0))

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
                    angular_matrix_element_array(1,ifg1,m1,ifg2,m2) = compute_angular_matrix_element(1,ifg1,m1,ifg2,m2)
                    angular_matrix_element_array(2,ifg1,m1,ifg2,m2) = compute_angular_matrix_element(2,ifg1,m1,ifg2,m2)
                    angular_matrix_element_array(3,ifg1,m1,ifg2,m2) = compute_angular_matrix_element(3,ifg1,m1,ifg2,m2)
                end do
            end do
    end do
    !$omp end parallel do
    precomputed_angular_matrix_element = .True.
end subroutine

subroutine precompute_JJ_matrix_element_one_body
    !--------------------------------------------------
    ! store the calculated \sum_i <n1 l1 j1 m1|J_i^2|n2 l2 j2 m2> 
    !--------------------------------------------------
    use Globals, only: BS
    integer :: ifg1,m1,ifg2,m2,dim_m_max,iter,total_iter
    dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
    if(.not. allocated(JJ_ME1B_array)) allocate(JJ_ME1B_array(2,dim_m_max,2,dim_m_max),source=(0.d0,0.d0))

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
                    JJ_ME1B_array(ifg1,m1,ifg2,m2) = compute_JJ_matrix_element_one_body(ifg1,m1,ifg2,m2)
                end do
            end do
    end do
    !$omp end parallel do

end subroutine

function compute_angular_matrix_element(case,ifg1,m1,ifg2,m2)
    !--------------------------------------------------
    !   calculate <n1 l1 j1 m1|J_i|n2 l2 j2 m2>
    !       where i = 1,2,3 corresponds to x,y,z
    ! ------------------------------------------------
    use Globals, only: BS
    complex(r64) :: compute_angular_matrix_element
    integer :: case,ifg1,m1,ifg2,m2
    integer ::  i0sp1, nr1, nl1, nj1, nm1
    integer ::  i0sp2, nr2, nl2, nj2, nm2
    real(r64) :: nj1_half, nm1_half
    compute_angular_matrix_element = (0.0d0,0.0d0)
    
    i0sp1 = BS%HO_sph%iasp(1,ifg1)
    nr1 = BS%HO_sph%nljm(i0sp1+m1,1)
    nl1 = BS%HO_sph%nljm(i0sp1+m1,2)
    nj1 = BS%HO_sph%nljm(i0sp1+m1,3)
    nm1 = BS%HO_sph%nljm(i0sp1+m1,4)
    nj1_half = nj1 - 0.5d0 ! j
    nm1_half = nm1 - 0.5d0 ! m_j

    i0sp2 = BS%HO_sph%iasp(1,ifg2)
    nr2 = BS%HO_sph%nljm(i0sp2+m2,1)
    nl2 = BS%HO_sph%nljm(i0sp2+m2,2)
    nj2 = BS%HO_sph%nljm(i0sp2+m2,3)
    nm2 = BS%HO_sph%nljm(i0sp2+m2,4)

    if(nr1/=nr2 .or. nl1 /= nl2 .or. nj1/= nj2) return 
    if(case==1) then
        ! <n1 l1 j1 m1|Jx|n2 l2 j2 m2>
        if (nm2==(nm1-1)) compute_angular_matrix_element = 0.5d0*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half-1.0d0))
        if (nm2==(nm1+1)) compute_angular_matrix_element = 0.5d0*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half+1.0d0))
    else if(case==2) then
        ! <n1 l1 j1 m1|Jy|n2 l2 j2 m2>
        if (nm2==(nm1-1)) compute_angular_matrix_element =  -(0.0, 0.5d0)*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half-1.0d0))
        if (nm2==(nm1+1)) compute_angular_matrix_element =   (0.0, 0.5d0)*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half+1.0d0)) 
    else if(case==3) then
        ! <n1 l1 j1 m1|Jz|n2 l2 j2 m2>
        if(nm1==nm2) compute_angular_matrix_element = nm1_half
    else 
        stop 'wrong case in compute_angular_matrix_element'
    end if
end 

double precision function compute_JJ_matrix_element_one_body(ifg1,m1,ifg2,m2)
    !-----------------------------------------------------------
    !  <n1 l1 j1 m1|J^2|n2 l2 j2 m2> = \sum_i <n1 l1 j1 m1|J_i^2|n2 l2 j2 m2>
    !   = \sum_i \sum_{m} <n1 l1 j1 m1|J_i|n l j m><n l j m|J_i|n2 l2 j2 m2>
    !-----------------------------------------------------------
    use Globals, only: BS
    integer,intent(in) :: ifg1,m1,ifg2,m2
    real(r64) :: EME1B
    integer :: ndsp,ifg,m,i
    EME1B = 0.d0
    ifg = ifg1
    ndsp = BS%HO_sph%idsp(1,ifg)
    EME1B = 0.d0
    do i = 1, 3
        do m = 1, ndsp
            EME1B = EME1B + angular_matrix_element_array(i,ifg1,m1,ifg,m)*angular_matrix_element_array(i,ifg,m,ifg2,m2)
        end do 
    end do 
    compute_JJ_matrix_element_one_body = EME1B
end function

! ------------------------------------------------------------------------
!                                N^2
! -------------------------------------------------------------------------
subroutine calculate_N2(iphi,it,N2,pN2)
    !-------------------------------------------------------------------------------------------------
    !
    !   <q1| N^2 R |q2>/<q1|R|q2>  =  \sum_{m m'} <q1|c^+_m c_m  c^+_m' c_m' R|q2>/<q1|R|q2>
    !     =   \sum_{m} rho_{m m} - \sum_{m m'} rho_{m' m m' m}
    !
    !  where             R  = R(alpha,beta,gamma, phi_n,phi_p) 
    !                       = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N} 
    !             rho_{m m} = <q1|c^+_m c_m R|q2>/<q1|R|q2>
    !       rho_{m' m m' m} = <q1|c^+_m c^+_m' c_m c_m' R|q2>/<q1|R|q2>
    !                       = rho_{m' m}rho_{m m'} - rho_{m m}rho_{m' m'} + kappa01*(m m')kappa10(m' m)
    !--------------------------------------------------------------------------------------------------
    use Globals, only: BS, mix
    integer :: iphi,it
    complex(r64) :: N2,pN2
    integer :: ifg1,m1,ifg2,m2,total_iter,iter

    N2 = (0.0d0,0.0d0)
    pN2 = (0.0d0,0.0d0)

    total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2) &
    !$omp reduction(+:N2,pN2)
    !$OMP DO COLLAPSE(1) SCHEDULE(static)

    do iter = 1,total_iter
        if(iter <= BS%HO_sph%idsp(1,1)) then
            ifg1 = 1
            m1 = iter
        else
            ifg1 = 2
            m1 = iter - BS%HO_sph%idsp(1,1)
        end if
        N2 = N2 + mix%rho_mm(m1,m1,ifg1,iphi,it)
        pN2 = pN2 + mix%prho_mm(m1,m1,ifg1,iphi,it)
        do ifg2 = 1,2
        do m2 = 1, BS%HO_sph%idsp(1,ifg2)
            N2 = N2 - (  mix%rho_mm(m2,m1,indexfg(ifg2,ifg1),iphi,it)*mix%rho_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it) &
                       - mix%rho_mm(m1,m1,indexfg(ifg1,ifg1),iphi,it)*mix%rho_mm(m2,m2,indexfg(ifg2,ifg2),iphi,it) &
                       + mix%kappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%kappa10_mm(m2,m1,indexfg(ifg2,ifg1),iphi,it))
            pN2 = pN2 - (  mix%prho_mm(m2,m1,indexfg(ifg2,ifg1),iphi,it)*mix%prho_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it) &
                         - mix%prho_mm(m1,m1,indexfg(ifg1,ifg1),iphi,it)*mix%prho_mm(m2,m2,indexfg(ifg2,ifg2),iphi,it) &
                         + mix%pkappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%pkappa10_mm(m2,m1,indexfg(ifg2,ifg1),iphi,it))
        end do
        end do
    end do
    !$OMP END PARALLEL
end subroutine


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
end Module