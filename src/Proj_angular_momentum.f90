!==============================================================================!
! MODULE AM                                                                    !
!                                                                              !
! This module calculates the angular momentum operators                        !
!                                                                              !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
Module AM
    use Constants,only: r64
    implicit none
    complex(r64), dimension(:,:,:,:,:), allocatable:: angular_matrix_element_array ! angular_matrix_element_array(Jx/Jy/Jz,ifg1, m1, ifg2, m2)
    logical :: precomputed_angular_matrix_element = .False.
contains


subroutine calculate_Jsquare(iphi,it,J2,pJ2)
    !---------------------------------------------------------------------------
    ! 
    !  calculate calculate <q1| J^2 R |q2>/<q1|R|q2>  = <q1| J_x^2 R |q2>/<q1|R|q2> 
    !                   + <q1| J_y^2 R |q2>/<q1|R|q2> + <q1| J_z^2 R |q2>/<q1|R|q2>
    ! where R = R(alpha,beta,gamma, phi_n,phi_p) 
    !   = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
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
    complex(r64) :: J2,pJ2,J12_J24,rho_4321,prho_4321,J13_J24
    integer :: ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,total_iter,iter

    J2 = (0.0d0,0.0d0)
    pJ2 = (0.0d0,0.0d0)

    if(.not. precomputed_angular_matrix_element) then
        call precompute_angular_matrix_element
        write(*,*) 'precompute_angular_matrix_element' 
    end if

    total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,J12_J24,rho_4321,prho_4321,J13_J24) &
    !$omp reduction(+:J2,pJ2)
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
            do ifg3 = 1,2
            do m3 = 1, BS%HO_sph%idsp(1,ifg3)
                do ifg4 = 1,2
                do m4 = 1, BS%HO_sph%idsp(1,ifg4)
                    if(ifg1==ifg3 .and. ifg2==ifg4) then
                        ! nonzero (J_i)_{m2 m4} 
                        if( angular_matrix_element(1,ifg2,m2,ifg4,m4)==(0.d0,0.d0) .and.  angular_matrix_element(2,ifg2,m2,ifg4,m4)==(0.d0,0.d0) & 
                            .or. angular_matrix_element(3,ifg2,m2,ifg4,m4)==(0.d0,0.d0) ) cycle
                        ! 1B
                        if(ifg2==ifg3 .and. m2==m3) then
                            if (abs(mix%rho_mm(m4,m1,ifg1,iphi,it)) > 1.E-10 .or. abs(mix%prho_mm(m4,m1,ifg1,iphi,it)) > 1.E-10) then 
                                J12_J24 = angular_matrix_element(1,ifg1,m1,ifg2,m2)*angular_matrix_element(1,ifg2,m2,ifg4,m4) &
                                        + angular_matrix_element(2,ifg1,m1,ifg2,m2)*angular_matrix_element(2,ifg2,m2,ifg4,m4) &
                                        + angular_matrix_element(3,ifg1,m1,ifg2,m2)*angular_matrix_element(3,ifg2,m2,ifg4,m4)
                                J2 = J2 + J12_J24 * mix%rho_mm(m4,m1,ifg1,iphi,it)
                                pJ2 = pJ2 + J12_J24 * mix%prho_mm(m4,m1,ifg1,iphi,it)
                            end if 
                        end if 
                        ! 2B
                        rho_4321  = mix%rho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%rho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                    - mix%rho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%rho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                    + mix%kappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%kappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it)
                        prho_4321 =  mix%prho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%prho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                    - mix%prho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%prho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                    + mix%pkappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%pkappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it)
                        J13_J24 = angular_matrix_element(1,ifg1,m1,ifg3,m3)*angular_matrix_element(1,ifg2,m2,ifg4,m4)  &
                                    + angular_matrix_element(2,ifg1,m1,ifg3,m3)*angular_matrix_element(2,ifg2,m2,ifg4,m4)  &
                                    + angular_matrix_element(3,ifg1,m1,ifg3,m3)*angular_matrix_element(3,ifg2,m2,ifg4,m4)
                                    J2 = J2 + J13_J24*rho_4321
                                    pJ2 = pJ2 + J13_J24*prho_4321
                        J2 = J2 + J13_J24*rho_4321
                        pJ2 = pJ2 + J13_J24*prho_4321
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
    nj1_half = nj1 - 0.5 ! j
    nm1_half = nm1 - 0.5 ! m_j

    i0sp2 = BS%HO_sph%iasp(1,ifg2)
    nr2 = BS%HO_sph%nljm(i0sp2+m2,1)
    nl2 = BS%HO_sph%nljm(i0sp2+m2,2)
    nj2 = BS%HO_sph%nljm(i0sp2+m2,3)
    nm2 = BS%HO_sph%nljm(i0sp2+m2,4)

    if(nr1/=nr2 .or. nl1 /= nl2 .or. nj1== nj2) return 

    if(case==1) then
        ! <n1 l1 j1 m1|Jx|n2 l2 j2 m2>
        if (m2==m1+1) compute_angular_matrix_element = 0.5d0*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half+1.0d0))
        if (m2==m1-1) compute_angular_matrix_element = 0.5d0*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half-1.0d0))
    else if(case==2) then
        ! <n1 l1 j1 m1|Jy|n2 l2 j2 m2>
        if (m2==m1+1) compute_angular_matrix_element = -(0.0, 0.5d0)*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half+1.0d0)) 
        if (m2==m1-1) compute_angular_matrix_element =  (0.0, 0.5d0)*sqrt(nj1_half*(nj1_half+1.0d0)-nm1_half*(nm1_half-1.0d0))
    else if(case==3) then
        ! <n1 l1 j1 m1|Jz|n2 l2 j2 m2>
        if(nm1==nm2) compute_angular_matrix_element = nm1_half
    else 
        stop 'wrong case in compute_angular_matrix_element'
    end if
end 

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

end Module AM