!==============================================================================!
! MODULE Mixed                                                                 !
!                                                                              !
! This module calculates the                                                   !
!      norm overlap                                                            !
!      mixed ... matrix elements                                               !
!      mixed ... in coordinate space                                           !
!                                                                              !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
Module Mixed
    use Constants, only: r64,nh2x,nkx
    use Globals, only: mix
    implicit none
    complex(r64),dimension(nh2x,-nkx:nkx,2) :: fg_rotated,pfg_rotated
    complex(r64),dimension(:,:,:),allocatable :: R, pR, RT_inv, pRT_inv, D, pD, D_inv, pD_inv
    complex(r64),dimension(2) :: RT_Deter, pRT_Deter, D_Deter, pD_Deter
contains

subroutine determine_truncated_dimension
    !----------------------------------------------------------------------------------
    !     determination of cutoff zeta in Dirac space zeta1 and zeta2 for different 
    ! mesh points in q-space.
    !     zeta1 and zeta2 are chosen as different values to make sure that the truncated 
    ! matrices T_22 or D, to be square !! 
    !
    ! Input/Output:
    !     * truncated_dim: Final truncated dimension
    !     * eps_occupation: truncation threshold
    !     * truncated_k: Mapping from original to truncated indices
    !-----------------------------------------------------------------------------------
    use Constants, only: r64,nkx,eps_occupation
    use Globals, only: wf1,wf2
    use MathMethods, only: sort
    integer :: it,k,truncated_dim_1,truncated_dim_2,truncated_index
    real(r64), dimension(nkx) :: sorted_v2d_1,sorted_v2d_2
    do it = 1, 2
        ! count states when 2*v^2 large than occupation epsilon
        truncated_dim_1 = 0
        do k = 1, wf1%nk(it)
            sorted_v2d_1(k) = wf1%v2d(k,it)
            if(wf1%v2d(k,it) > eps_occupation) truncated_dim_1 = truncated_dim_1 +1
        end do
        truncated_dim_2 = 0
        do k = 1, wf2%nk(it) 
            sorted_v2d_2(k) = wf2%v2d(k,it)
            if(wf2%v2d(k,it) > eps_occupation) truncated_dim_2 = truncated_dim_2 +1
        end do
        ! order the occupation probabilities
        call sort(wf1%nk(it), sorted_v2d_1, descending=.True.)
        call sort(wf2%nk(it), sorted_v2d_2, descending=.True.)
        ! choose the minimal dimension and set the truncation epsilon 
        if(truncated_dim_1 < truncated_dim_2) then
            wf1%truncated_dim(it)= truncated_dim_1
            wf2%truncated_dim(it)= truncated_dim_1
            wf1%eps_occupation(it) = eps_occupation
            wf2%eps_occupation(it) = sorted_v2d_2(truncated_dim_1+1)
        else if(truncated_dim_1 > truncated_dim_2) then
            wf1%truncated_dim(it)= truncated_dim_2
            wf2%truncated_dim(it)= truncated_dim_2
            wf1%eps_occupation(it) = sorted_v2d_1(truncated_dim_2+1)
            wf2%eps_occupation(it) = eps_occupation
        else
            wf1%truncated_dim(it)= truncated_dim_1
            wf2%truncated_dim(it)= truncated_dim_2
            wf1%eps_occupation(it) = eps_occupation
            wf2%eps_occupation(it) = eps_occupation
        endif
        ! set index
        truncated_index = 0
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) > wf1%eps_occupation(it)) then 
                truncated_index = truncated_index +1
                wf1%truncated_index(k,it) = truncated_index
            endif
        end do
        truncated_index = 0
        do k = 1, wf2%nk(it)
            if(wf2%v2d(k,it) > wf2%eps_occupation(it)) then 
                truncated_index = truncated_index +1
                wf2%truncated_index(k,it) = truncated_index
            endif
        end do
    end do
    return
end subroutine


subroutine calculate_mixed_DensCurrTens_and_norm_overlap(alpha, beta, gamma)
    use Globals, only: Proj_option,projection_mesh
    real(r64), intent(in) :: alpha, beta, gamma
    integer :: it, iphi
    real(r64) :: phi
    call initialize_variables
    do it = 1,2
        call calculate_Rotated_Wavefunction(alpha, beta, gamma, it)
        call calculate_Rotation_Matrix(alpha, beta, gamma, it)
        do iphi = 1, projection_mesh%nphi(it)
            phi = iphi*projection_mesh%dphi(it)
            call calculate_D_Matrix(phi,it)
            ! rho_mm', kappa_mm', kappa*_mm'
            call calculate_mixed_density_tensor_matrix_elements(iphi,phi,it)
            ! rho_S_it, rho_V_it, ...
            call calculate_mixed_density_current_tensor_in_coordinate_space(iphi,it)
            ! <q1|R(alpha,beta,gamma)e^{i phi_it N}|q2>
            call calculate_norm_overlap(iphi,phi,it)
        end do
    end do
    call deallocate_variables
end subroutine

subroutine initialize_variables
    use Constants, only: ngr,ntheta,nphi,itx
    use Globals, only: gcm_space,kernels,wf1,wf2,BS,projection_mesh,mix
    integer :: Jmax,truncated_dim_max_1,truncated_dim_max_2,dim_m_max,nphi_max

    truncated_dim_max_1 = max(wf1%truncated_dim(1),wf1%truncated_dim(2))
    truncated_dim_max_2 = max(wf2%truncated_dim(1),wf2%truncated_dim(2))
    allocate(R(2*truncated_dim_max_1, 2*truncated_dim_max_2, 2))
    allocate(pR(2*truncated_dim_max_1, 2*truncated_dim_max_2, 2))
    allocate(RT_inv(2*truncated_dim_max_1,2*truncated_dim_max_2,2))
    allocate(pRT_inv(2*truncated_dim_max_1,2*truncated_dim_max_2,2))
    allocate(D(2*truncated_dim_max_1, 2*truncated_dim_max_2, 2))
    allocate(pD(2*truncated_dim_max_1, 2*truncated_dim_max_2, 2))
    allocate(D_inv(2*truncated_dim_max_1, 2*truncated_dim_max_2, 2))
    allocate(pD_inv(2*truncated_dim_max_1, 2*truncated_dim_max_2, 2))

    ! --- Global mix ---
    nphi_max = max(projection_mesh%nphi(1),projection_mesh%nphi(2))
    if(.not. allocated(mix%norm)) allocate(mix%norm(nphi_max,2))
    if(.not. allocated(mix%pnorm)) allocate(mix%pnorm(nphi_max,2))

    dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
    ! rho_mm^{++/--/+-/-+}(alpha,beta,gamma,phi,it) with fixed alpha, beta,gamma
    if(.not. allocated(mix%rho_mm)) allocate(mix%rho_mm(dim_m_max,dim_m_max,4,nphi_max,itx))      ! 4: \rho^++,\rho^--,\rho^+-,\rho^-+; itx: nutron, proton  
    if(.not. allocated(mix%prho_mm)) allocate(mix%prho_mm(dim_m_max,dim_m_max,4,nphi_max,itx))
    if(.not. allocated(mix%kappa_mm)) allocate(mix%kappa_mm(dim_m_max,dim_m_max,2,nphi_max,itx))  ! 2: \kappa^++,\kappa^--; itx: nutron, proton  
    if(.not. allocated(mix%pkappa_mm)) allocate(mix%pkappa_mm(dim_m_max,dim_m_max,2,nphi_max,itx))
    if(.not. allocated(mix%kappac_mm)) allocate(mix%kappac_mm(dim_m_max,dim_m_max,2,nphi_max,itx))  ! 2: \kappa^++,\kappa^--; itx: nutron, proton
    if(.not. allocated(mix%pkappac_mm)) allocate(mix%pkappac_mm(dim_m_max,dim_m_max,2,nphi_max,itx))

    if(.not. allocated(mix%kappa10_mm)) allocate(mix%kappa10_mm(dim_m_max,dim_m_max,4,nphi_max,itx))  ! 2: \kappa10^++,\kappa10^--, \kappa10^+-,\kappa10^-+; itx: nutron, proton  
    if(.not. allocated(mix%pkappa10_mm)) allocate(mix%pkappa10_mm(dim_m_max,dim_m_max,4,nphi_max,itx))
    if(.not. allocated(mix%kappa01c_mm)) allocate(mix%kappa01c_mm(dim_m_max,dim_m_max,4,nphi_max,itx))  ! 2: \kappa01*^++,\kappa01*^--,\kappa01*^+-,\kappa01*^-+; itx: nutron, proton
    if(.not. allocated(mix%pkappa01c_mm)) allocate(mix%pkappa01c_mm(dim_m_max,dim_m_max,4,nphi_max,itx))

    ! rho_S_it(alpha,beta,gamma,phi,it) with fixed alpha, beta,gamma
    if(.not. allocated(mix%rho_S_it)) allocate(mix%rho_S_it(ngr*ntheta*nphi,nphi_max,itx)) ! ((x,theta,phi),phi_it,it)
    if(.not. allocated(mix%rho_V_it)) allocate(mix%rho_V_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%d2rho_S_it)) allocate(mix%d2rho_S_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%d2rho_V_it)) allocate(mix%d2rho_V_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%prho_S_it)) allocate(mix%prho_S_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%prho_V_it)) allocate(mix%prho_V_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%pd2rho_S_it)) allocate(mix%pd2rho_S_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%pd2rho_V_it)) allocate(mix%pd2rho_V_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%j_V1_it)) allocate(mix%j_V1_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%j_V2_it)) allocate(mix%j_V2_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%j_V3_it)) allocate(mix%j_V3_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%pj_V1_it)) allocate(mix%pj_V1_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%pj_V2_it)) allocate(mix%pj_V2_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%pj_V3_it)) allocate(mix%pj_V3_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%kappa_it)) allocate(mix%kappa_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%kappac_it)) allocate(mix%kappac_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%pkappa_it)) allocate(mix%pkappa_it(ngr*ntheta*nphi,nphi_max,itx))
    if(.not. allocated(mix%pkappac_it)) allocate(mix%pkappac_it(ngr*ntheta*nphi,nphi_max,itx))
end subroutine

subroutine deallocate_variables
    deallocate(R,RT_inv,D,D_inv,pR,pRT_inv,pD,pD_inv)
end subroutine  

subroutine calculate_Rotated_Wavefunction(alpha, beta, gamma,it)
    !---------------------------------------------------------------------------------------------------------------------------
    !  Rotate the wave function coefficients (q_2) for large (F) and small (G) components:
    !   F(Omega)_mk = \sum_{m'} R(Omega)_mm' F_{m'k}
    !   G(Omega)_nk = \sum_{n'} R(Omega)_nn' G_{n'k}
    !   ....
    ! 
    ! Omega(alpha,beta,gamma): Rotation parameters
    !      m,       m'       : Indices spherical HO basis of the large component, and |m> = |n_r l j m_j>
    !      n,       n'       : Indices spherical HO basis of the small component
    !          k             : Indices of canonical basis (or single particle basis)
    ! R(Omega)_mm' = <m|R(alpha,beta,gamma)|m'> 
    !              = \delta_{n_r n'_r} \delta_{ll'} \delta_{jj'}e^{i m_j alpha}d^j_{m_j m'_j}e^{i m_j gamma}
    ! ------------
    ! Time Reverse：
    !   T|m> = (-1)^{l+j-m_j}|n_r l j -m_j> = (-1)^{l+j-m_j}|-m>， where |-m> = |n_r l j -m_j>
    !                             (      \sum_m F_{m Tk} |m>      )       (      \sum_m F_{mk} |m>     )
    !   T|k> = |Tk> = |\bar{k}> = (                               )  = T  (                            )
    !                             (     i\sum_n G_{n Tk} |n>      )       (      i\sum_n G_{nk} |n>    )
    !
    !         (   \sum_m  F_{mk} (-1)^{l+j-m_j}  |-m> )   (    \sum_m F_{-mk} (-1)^{l+j+m_j} |m> )   
    !       = (                                       ) = (                                      )
    !         ( -i\sum_n  G_{nk} (-1)^{l+j-m_j}  |-n> )   (  -i\sum_n G_{-nk} (-1)^{l+j+m_j} |n> )
    !
    ! Parity：  
    !    P|m> = (-1)^l |m>
    !                   (   \sum_m F_{m Pk} |m>    )      (   \sum_m F_{mk} |m> )
    !    P|k> = | Pk> = (                          ) =  P (                     )
    !                   (  i\sum_n G_{n Pk} |n>    )      (  i\sum_n G_{nk} |n> )
    !
    !           (   \sum_m F_{mk} P|m> )   (   \sum_m (-1)^l F_{mk} |m> )
    !       =   (                      ) = (                            )
    !           ( -i\sum_n G_{nk} P|n> )   ( -i\sum_n (-1)^l G_{nk} |n> )
    ! ----------
    ! Output:
    !   F(Omega)_mk     = \sum_{m'} R(Omega)_mm' F_{m'  k}
    !   F(Omega)_m{Tk}  = \sum_{m'} R(Omega)_mm' F_{m' Tk} 
    !                   = \sum_{m'} R(Omega)_mm' F_{-m'k} (-1)^{l'+j'+m'_j} 
    !                   = \sum_{m'} R(Omega)_m-m' F_{m'k} (-1)^{l'+j'-m'_j}
    !   F(Omega)_m{Pk}  = \sum_{m'} R(Omega)_mm' F_{m' Pk} 
    !                   = \sum_{m'} R(Omega)_mm' F_{m'k} (-1)^l'
    !   F(Omega)_m{TPk} = \sum_{m'} R(Omega)_m-m' F_{m'k} (-1)^{l'+j'-m'_j+l'}
    !   G(Omega)_nk     = \sum_{n'} R(Omega)_nn' G_{n'k}
    !   G(Omega)_n{Tk}  = \sum_{n'} R(Omega)_nn' G_{n' Tk}
    !                   = \sum_{n'} R(Omega)_nn' G_{-nk} (-1)^{l'+j'+m'_j+1}
    !                   = \sum_{n'} R(Omega)_n-n' G_{nk} (-1)^{l'+j'-m'_j+1}
    !   G(Omega)_n{Pk}  = \sum_{n'} R(Omega)_nn' G_{n' Pk}
    !                   = \sum_{n'} R(Omega)_nn' G_{nk} (-1)^{l'+1}
    !   G(Omega)_n{TPk} = \sum_{n'} R(Omega)_n-n' G_{nk} (-1)^{l'+j'-m'_j+1+l'+1}
    !------------------------------------------------------------------------------------------------------------------------
    use Globals, only: BS,wf2
    use Basis, only: djmk
    integer ::  it,ib,ifg,i0sp,ndsp,nfgsp,k,m2,nr2,nl2,nj2,nm2,m3,nr3,nl3,nj3,nm3,&
                m3_reversed,nr3_reversed,nl3_reversed,nj3_reversed,nm3_reversed
    real(r64) :: alpha,beta,gamma,cos_beta,nj2_half,nm2_half,nj3_half,nm3_half,nm3_reversed_half
    complex(r64) :: calpha, cgamma,fgrot,fgrotr,pfgrot,pfgrotr,Rm2m3_11,Rm2m3_12
    
    calpha = DCMPLX(0.d0,alpha)
    cgamma = DCMPLX(0.d0,gamma)
    cos_beta = dcos(beta)
    ib = 1
    do ifg = 1, 2 ! loop over large and small components
        i0sp = BS%HO_sph%iasp(ib,ifg)
        ndsp = BS%HO_sph%idsp(ib,ifg)
        nfgsp = (ifg-1)*BS%HO_sph%idsp(ib,1)
        do k = 1, wf2%nk(it) ! loop over single particle basis
            if(wf2%v2d(k,it) <= wf2%eps_occupation(it)) cycle
            do m2 = 1, ndsp ! loop over sph. HO basis
                nr2= BS%HO_sph%nljm(i0sp+m2,1) ! n_r
                nl2= BS%HO_sph%nljm(i0sp+m2,2) ! l
                nj2= BS%HO_sph%nljm(i0sp+m2,3) ! j +1/2
                nm2= BS%HO_sph%nljm(i0sp+m2,4) ! m_j + 1/2
                nj2_half = nj2 - 0.5 ! j
                nm2_half = nm2 - 0.5 ! m_j
                fgrot = (0.d0,0.d0)
                fgrotr = (0.d0,0.d0)
                pfgrot = (0.d0,0.d0)
                pfgrotr = (0.d0,0.d0)
                do m3 = 1, ndsp
                    nr3= BS%HO_sph%nljm(i0sp+m3,1) ! n_r
                    nl3= BS%HO_sph%nljm(i0sp+m3,2) ! l
                    nj3= BS%HO_sph%nljm(i0sp+m3,3) ! j +1/2
                    nm3= BS%HO_sph%nljm(i0sp+m3,4) ! m_j + 1/2
                    nj3_half = nj3 - 0.5 ! j
                    nm3_half = nm3 - 0.5 ! m_j

                    ! time reversed |-m> = |n_r l j -m_j>
                    ! m3_reversed = m3 - 2*nm3 + 1
                    ! nr3_reversed= BS%HO_sph%nljm(i0sp+m3_reversed,1) ! n_r
                    ! nl3_reversed= BS%HO_sph%nljm(i0sp+m3_reversed,2) ! l
                    ! nj3_reversed= BS%HO_sph%nljm(i0sp+m3_reversed,3) ! j +1/2
                    ! nm3_reversed= BS%HO_sph%nljm(i0sp+m3_reversed,4) ! - m_j+1/2
                    ! if(nr3.ne.nr3_reversed .and. nl3 .ne. nl3_reversed .and. nj3.ne.nj3_reversed) stop 'Wrong time reversed!'
                    ! if(nm3 + nm3_reversed .ne. 1) stop 'wrong nm3_reversed'
                    nm3_reversed = 1 - nm3 ! -m_j + 1/2
                    nm3_reversed_half = nm3_reversed - 0.5 ! -m_j
                    if(nr3.ne.nr2 .or. nl3.ne.nl2 .or. nj3.ne.nj2) cycle
                    Rm2m3_11=CDEXP(nm2_half*calpha+nm3_half*cgamma)*djmk(nj2,nm2,nm3,cos_beta,1) ! <m|R(alpha,beta,gamma)|m'>
                    Rm2m3_12=CDEXP(nm2_half*calpha+nm3_reversed_half*cgamma)*djmk(nj2,nm2,nm3_reversed,cos_beta,1) ! <m|R(alpha,beta,gamma)|-m'>
                    fgrot   = fgrot   + Rm2m3_11*wf2%fg(nfgsp+m3,k,it)
                    fgrotr  = fgrotr  + Rm2m3_12*wf2%fg(nfgsp+m3,k,it)*(-1)**(nl3+nj3-nm3+ifg-1)
                    pfgrot  = pfgrot  + Rm2m3_11*wf2%fg(nfgsp+m3,k,it)*(-1)**(nl3+ifg-1)
                    pfgrotr = pfgrotr + Rm2m3_12*wf2%fg(nfgsp+m3,k,it)*(-1)**(nl3+nj3-nm3+ifg-1+nl3+ifg-1)
                end do
                if((nfgsp+m2)> nh2x .or. k> nkx ) stop "nh2x or nkx too small!"
                fg_rotated(nfgsp+m2, k,it)  = fgrot
                fg_rotated(nfgsp+m2,-k,it)  = fgrotr
                pfg_rotated(nfgsp+m2,k,it)  = pfgrot
                pfg_rotated(nfgsp+m2,-k,it) = pfgrotr
            end do
        end do
    end do 
end subroutine

subroutine calculate_Rotation_Matrix(alpha,beta,gamma,it)
    !---------------------------------------------------------------------------------------------------------------------
    ! calculate the rotation matrix R (in the canonical basis), the (R^T)^{-1} and the determinant of R^T
    !   R_{kl} = <k| R(alpha,beta,gamma) |l> = <k| R(Omega) |l> 
    !          = \sum_{mm'} F^*_{mk}(q_1) F_{m'l}(q_2) R(Omega)_mm'
    !          + \sum_{nn'} G^*_{nk}(q_1) G_{n'l}(q_2) R(Omega)_nn'
    !   R_{k\bar{l}} = <k| R(alpha,beta,gamma) |\bar{l}> = <k| R(Omega) |\bar{l}> 
    !          = \sum_{mm'} F^*_{mk}(q_1) F_{m'l}(q_2) (-1)^{l'+j'-m'_j}  R(Omega)_m-m'
    !          + \sum_{nn'} G^*_{nk}(q_1) G_{n'l}(q_2) (-1)^{l'+j'-m'_j+1}R(Omega)_n-n'
    !   R_\bar{k}l = - R*_k\bar{l}
    !   R_\bar{k}\bar{l} = R*_kl
    !
    !   where |k> is the canonical basis state and |\bar{k}> = T|k>  is the time-reversed state of |k>.
    !   
    ! Returns:
    !   R, RT_inverse, RT_Deter
    !   pR, pRT_inverse, pRT_Deter
    !---------------------------------------------------------------------------------------------------------------------
    use Globals, only: wf1,wf2,BS
    use Basis, only: djmk
    use MathMethods, only: clingd
    integer,intent(in)::  it
    real(r64),intent(in) :: alpha,beta,gamma
    integer :: l,k,ifg,i0sp,ndsp,nfgsp,m2,m1,nr2,nl2,nj2,nm2,nr1,nl1,nj1,nm1,nm2_reversed,nm1_reversed,&
                truncated_k,truncated_l,truncated_dim,ifl
    real(r64) :: cos_beta,nj2_half,nm2_half,nm2_reversed_half,nj1_half,nm1_half,nm1_reversed_half
    complex(r64) :: calpha, cgamma,Rm1m2_11,Rm1m2_12,Rkl_11,Rkl_12,pRkl_11,pRkl_12,tmp_RT_Deter,tmp_pRT_Deter
    complex(r64),dimension(:,:), allocatable :: tmp_RT,tmp_RT_inv,tmp_pRT,tmp_pRT_inv
    if (wf1%truncated_dim(it) /=wf2%truncated_dim(it) ) stop "Wrong truncated dimension!"
    truncated_dim = wf1%truncated_dim(it)
    ! calculate R 
    calpha = DCMPLX(0.d0,alpha)
    cgamma = DCMPLX(0.d0,gamma)
    cos_beta = dcos(beta)
    do l = 1, wf2%nk(it)
        if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle 
            Rkl_11 = (0.d0,0.d0)
            Rkl_12 = (0.d0,0.d0)
            pRkl_11 = (0.d0,0.d0)
            pRkl_12 = (0.d0,0.d0)
            do ifg = 1, 2
                i0sp = BS%HO_sph%iasp(1,ifg)
                ndsp = BS%HO_sph%idsp(1,ifg)
                nfgsp = (ifg-1)*BS%HO_sph%idsp(1,1)
                do m2 = 1, ndsp
                    nr2= BS%HO_sph%nljm(i0sp+m2,1) ! n'_r
                    nl2= BS%HO_sph%nljm(i0sp+m2,2) ! l'
                    nj2= BS%HO_sph%nljm(i0sp+m2,3) ! j' +1/2
                    nm2= BS%HO_sph%nljm(i0sp+m2,4) ! m'_j + 1/2
                    nj2_half = nj2 - 0.5 ! j'
                    nm2_half = nm2 - 0.5 ! m'_j
                    nm2_reversed = 1 - nm2 ! -m'_j + 1/2
                    nm2_reversed_half = nm2_reversed - 0.5 ! -m'_j
                    do m1 = 1, ndsp
                        nr1= BS%HO_sph%nljm(i0sp+m1,1) ! n_r
                        nl1= BS%HO_sph%nljm(i0sp+m1,2) ! l
                        nj1= BS%HO_sph%nljm(i0sp+m1,3) ! j +1/2
                        nm1= BS%HO_sph%nljm(i0sp+m1,4) ! m_j + 1/2
                        nj1_half = nj1 - 0.5 ! j
                        nm1_half = nm1 - 0.5 ! m_j
                        nm1_reversed = 1 - nm1 ! -m_j + 1/2
                        nm1_reversed_half = nm1_reversed - 0.5 ! -m_j
                        if(nr1.ne.nr2 .or. nl1.ne.nl2 .or. nj1.ne.nj2) cycle
                        Rm1m2_11=CDEXP(nm1_half*calpha+nm2_half*cgamma)*djmk(nj1,nm1,nm2,cos_beta,1) ! <m|R(alpha,beta,gamma)|m'>
                        Rm1m2_12=CDEXP(nm1_half*calpha+nm2_reversed_half*cgamma)*djmk(nj1,nm1,nm2_reversed,cos_beta,1) ! <m|R(alpha,beta,gamma)|-m'>
                        Rkl_11  = Rkl_11  + wf1%fg(nfgsp+m1,k,it)*wf2%fg(nfgsp+m2,l,it)*Rm1m2_11
                        Rkl_12  = Rkl_12  + wf1%fg(nfgsp+m1,k,it)*wf2%fg(nfgsp+m2,l,it)*Rm1m2_12*(-1)**(nl2+nj2-nm2+ifg-1)
                        pRkl_11 = pRkl_11 + wf1%fg(nfgsp+m1,k,it)*wf2%fg(nfgsp+m2,l,it)*Rm1m2_11*(-1)**(nl2+ifg-1)
                        pRkl_12 = pRkl_12 + wf1%fg(nfgsp+m1,k,it)*wf2%fg(nfgsp+m2,l,it)*Rm1m2_12*(-1)**(nl2+nj2-nm2+ifg-1+nl2+ifg-1)
                    end do
                end do
            end do
            truncated_k = wf1%truncated_index(k,it)
            truncated_l = wf2%truncated_index(l,it)     
            ! store R: 2truncated_dim x 2truncated_dim
            R(truncated_k,truncated_l,it) = Rkl_11                                      ! R_kl
            R(truncated_k,truncated_dim+truncated_l,it) = Rkl_12                        ! R_k\bar{l}
            R(truncated_dim+truncated_k,truncated_l,it) = -DCONJG(Rkl_12)               ! R_\bar{k}l
            R(truncated_dim+truncated_k,truncated_dim+truncated_l,it) = DCONJG(Rkl_11)  ! R_\bar{k}\bar{l}
            pR(truncated_k,truncated_l,it) = pRkl_11
            pR(truncated_k,truncated_dim+truncated_l,it) = pRkl_12
            pR(truncated_dim+truncated_k,truncated_l,it) = -DCONJG(pRkl_12)
            pR(truncated_dim+truncated_k,truncated_dim+truncated_l,it) = DCONJG(pRkl_11)
        end do
    end do
    if(truncated_k/= wf1%truncated_dim(it)) stop "truncated_k = wf1%truncated_dim(it)"
    if(truncated_l/= wf2%truncated_dim(it)) stop "truncated_l = wf2%truncated_dim(it)"
    ! calculate the inverse of R^T ((R^T)^{-1}) and its determinant 
    allocate(tmp_RT(2*truncated_dim, 2*truncated_dim))
    allocate(tmp_RT_inv(2*truncated_dim, 2*truncated_dim))
    allocate(tmp_pRT(2*truncated_dim, 2*truncated_dim))
    allocate(tmp_pRT_inv(2*truncated_dim, 2*truncated_dim))
    do l = 1, 2*truncated_dim
        do k = 1, 2*truncated_dim
            tmp_RT(k,l) = R(l,k,it) !R^T
            tmp_RT_inv(k,l) = 0.d0
            tmp_pRT(k,l) = pR(l,k,it)
            tmp_pRT_inv(k,l) = 0.d0
        end do
        tmp_RT_inv(l,l) = CMPLX(1.d0,0.d0)
        tmp_pRT_inv(l,l) = CMPLX(1.d0,0.d0)
    end do
    call clingd(2*truncated_dim,2*truncated_dim,2*truncated_dim,2*truncated_dim,tmp_RT,tmp_RT_inv,tmp_RT_Deter,ifl)
    if(ifl==-1) stop 'WARNING!!!! Matrix R has zero determinant !!!! Redefine the Hilbert Space by changing eps_occupation in Constants Module'
    call clingd(2*truncated_dim,2*truncated_dim,2*truncated_dim,2*truncated_dim,tmp_pRT,tmp_pRT_inv,tmp_pRT_Deter,ifl)
    if(ifl==-1) stop 'WARNING!!!! Matrix pR has zero determinant !!!! Redefine the Hilbert Space by changing eps_occupation in Constants Module'
    ! store the (R^T)^{-1}
    do l = 1, 2*truncated_dim
        do k = 1, 2*truncated_dim
            RT_inv(k,l,it) = tmp_RT_inv(k,l)
            pRT_inv(k,l,it) = tmp_pRT_inv(k,l)
        end do
    end do
    ! store the determinant (R^T)^{-1}
    RT_Deter(it) = tmp_RT_Deter
    pRT_Deter(it) = tmp_pRT_Deter
end subroutine

subroutine calculate_D_Matrix(phi,it)
    !--------------------------------------------------------------------------------------
    !  calculate D matrix,  inverse of D matrix D^{-1} and the determinant of D
    !          D = U^T(q_1) (R^T)^{-1} U*(q_2) + V^T(q_1) R^T V*(q_2)
    ! 
    ! Note:
    !   R_\bar{k}l = - R*_k\bar{l}
    !   R_\bar{k}\bar{l} = R*_kl
    !
    ! Returns:
    !   D, D_inverse, D_Deter
    !   pD, pD_inverse, pD_Deter
    !--------------------------------------------------------------------------------------
    use Globals, only: wf1,wf2
    use MathMethods, only: clingd
    integer,intent(in)::  it
    real(r64),intent(in) :: phi
    complex(r64) :: cphi,vl,ei2phi,tmp_D_Deter,tmp_pD_Deter
    integer :: truncated_dim,k,l,truncated_k,truncated_l,ifl
    real(r64) :: uk,vk,ul
    complex(r64),dimension(:,:), allocatable :: tmp_D, tmp_D_inv, tmp_pD, tmp_pD_inv
    cphi   = cmplx(0.d0,phi)
    ei2phi = cdexp(2*cphi) ! e^{2i\phi}
    if (wf1%truncated_dim(it) /=wf2%truncated_dim(it) ) stop "Wrong truncated dimension!"
    truncated_dim = wf1%truncated_dim(it)
    ! calculate D
    do k = 1, wf1%nk(it)
        if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
        vk = dsqrt(wf1%v2(k,it))
        uk = dsqrt(1-wf1%v2(k,it))
        do l = 1 , wf2%nk(it)
            if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
            vl = dsqrt(wf2%v2(l,it))*ei2phi
            ul = dsqrt(1-wf2%v2(l,it))
            truncated_k = wf1%truncated_index(k,it)
            truncated_l = wf2%truncated_index(l,it)
            D(truncated_k,truncated_l,it) = uk*ul*RT_inv(truncated_k,truncated_l,it) &                             ! D_kl
                                        + vk*vl*DCONJG(R(truncated_k,truncated_l,it))
            D(truncated_k,truncated_dim+truncated_l,it) = uk*ul*RT_inv(truncated_k,truncated_dim+truncated_l,it) & ! D_k\bar{l}
                                        + vk*vl*DCONJG(R(truncated_k,truncated_dim+truncated_l,it))
            D(truncated_dim+truncated_k,truncated_l,it) = uk*ul*RT_inv(truncated_dim+truncated_k,truncated_l,it) & ! D_k\bar{l}
                                        + vk*vl*DCONJG(R(truncated_dim+truncated_k,truncated_l,it))
            D(truncated_dim+truncated_k,truncated_dim+truncated_l,it)  &                                               ! D_\bar{k}\bar{l}
                                        = uk*ul*RT_inv(truncated_dim+truncated_k,truncated_dim+truncated_l,it) &
                                        + vk*vl*DCONJG(R(truncated_dim+truncated_k,truncated_dim+truncated_l,it))
            
            ! check D_{\bar{k} \bar{l}} = D_{k l}^* and D_{\bar{k} l} = - D_{k \bar{l}}^* 
            ! but may not satisfied after PNP 
            ! if(D(truncated_dim+truncated_k,truncated_l,it) /= - DCONJG(D(truncated_k,truncated_dim+truncated_l,it))) stop "Wrong D matrix!"
            ! if(D(truncated_dim+truncated_k,truncated_dim+truncated_l,it) /= DCONJG(D(truncated_k,truncated_l,it)))   stop "Wrong D matrix!"
            
            ! parity part
            pD(truncated_k,truncated_l,it) = uk*ul*pRT_inv(truncated_k,truncated_l,it) &                             ! D_kl
                                        + vk*vl*DCONJG(pR(truncated_k,truncated_l,it))
            pD(truncated_k,truncated_dim+truncated_l,it) = uk*ul*pRT_inv(truncated_k,truncated_dim+truncated_l,it) & ! D_k\bar{l}
                                        + vk*vl*DCONJG(pR(truncated_k,truncated_dim+truncated_l,it))
            pD(truncated_dim+truncated_k,truncated_l,it) = uk*ul*pRT_inv(truncated_dim+truncated_k,truncated_l,it) & ! D_k\bar{l}
                                        + vk*vl*DCONJG(pR(truncated_dim+truncated_k,truncated_l,it))
            pD(truncated_dim+truncated_k,truncated_dim+truncated_l,it)  &                                               ! D_\bar{k}\bar{l}
                                        = uk*ul*pRT_inv(truncated_dim+truncated_k,truncated_dim+truncated_l,it) &
                                        + vk*vl*DCONJG(pR(truncated_dim+truncated_k,truncated_dim+truncated_l,it))
            ! but may not satisfied after PNP 
            ! if(pD(truncated_dim+truncated_k,truncated_l,it) /= - DCONJG(pD(truncated_k,truncated_dim+truncated_l,it))) stop "Wrong pD matrix!"
            ! if(pD(truncated_dim+truncated_k,truncated_dim+truncated_l,it) /= DCONJG(pD(truncated_k,truncated_l,it)))   stop "Wrong pD matrix!"
        end do
    end do
    ! calculate the inverse of D matrix and its determinant
    allocate(tmp_D(2*truncated_dim, 2*truncated_dim))
    allocate(tmp_D_inv(2*truncated_dim, 2*truncated_dim))
    allocate(tmp_pD(2*truncated_dim, 2*truncated_dim))
    allocate(tmp_pD_inv(2*truncated_dim, 2*truncated_dim))
    do l = 1, 2*truncated_dim
        do k = 1, 2*truncated_dim
            tmp_D(k,l) = D(k,l,it)
            tmp_D_inv(k,l) = (0.d0,0.d0)
            tmp_pD(k,l) = pD(k,l,it)
            tmp_pD_inv(k,l) = (0.d0,0.d0)
        end do
        tmp_D_inv(l,l) = CMPLX(1.d0,0.d0)
        tmp_pD_inv(l,l) = CMPLX(1.d0,0.d0)
    end do
    call clingd(2*truncated_dim,2*truncated_dim,2*truncated_dim,2*truncated_dim,tmp_D,tmp_D_inv,tmp_D_Deter,ifl)
    if(ifl==-1) stop 'WARNING!!!! Matrix D has zero determinant !!!! '
    call clingd(2*truncated_dim,2*truncated_dim,2*truncated_dim,2*truncated_dim,tmp_pD,tmp_pD_inv,tmp_pD_Deter,ifl)
    if(ifl==-1) stop 'WARNING!!!! Matrix pD has zero determinant !!!! '
    ! store the D^{-1}
    do l = 1, 2*truncated_dim
        do k = 1, 2*truncated_dim
            D_inv(k,l,it) = tmp_D_inv(k,l)
            pD_inv(k,l,it) = tmp_pD_inv(k,l)
        end do
    end do
    ! store the determinant of D: det(D)
    D_Deter(it) = tmp_D_Deter
    pD_Deter(it) = tmp_pD_Deter
    deallocate(tmp_D, tmp_D_inv, tmp_pD, tmp_pD_inv)
end subroutine

subroutine calculate_mixed_density_tensor_matrix_elements(iphi,phi,it)
    !--------------------------------------------------------------------------------------
    !   calculate mixed density matrix elements: rho_mm' = <q1|c^+_m' c_m R|q2>/<q1|R|q2>.
    !   calculate mixed tensor matrix elements: kappa_mm' = <q1|c_m' c_m R|q2>/<q1|R|q2>.
    !   calculate mixed tensor matrix elements: kappa*_mm' = <q1|c+_m' c+_m R|q2>/<q1|R|q2>.
    !   
    !   calculate mixed tensor matrix elements: kappa10_mm' = <q1|c_m' c_m R|q2>/<q1|R|q2>.
    !   calculate mixed tensor matrix elements: kappa01*_mm' = <q1|c+_m' c+_m R|q2>/<q1|R|q2>.
    !--------------------------------------------------------------------------------------
    use Globals, only: wf1,wf2,BS
    use MathMethods, only : cmatmul_ABC
    integer,intent(in)::  iphi,it
    real(r64), intent(in) :: phi
    complex(r64) :: cphi,ei2phi
    cphi   = DCMPLX(0.d0,phi)
    ei2phi = cdexp(2*cphi)
    call mixed_density_matrix_elments
    call mixed_tensor_matrix_elements
    call mixed_kappa_matrix_elements(with_f=.False.)
    contains
    subroutine mixed_density_matrix_elments
        !--------------------------------------------------------------------------------------
        !  calculate mixed density matrix elements:
        !   rho_mm' = <q1|c^+_m' c_m R|q2>/<q1|R|q2>.
        !
        !  \rho_mm'++ = R_mm'' (F_2)_m''q (V_2*)_ql (D^-1)_lk (V_1^T)_kp (F_1^T)_pm'
        !  \rho_nn'-- = R_nn'' (G_2)_n''q (V_2*)_ql (D^-1)_lk (V_1^T)_kp (G_1^T)_pn'
        !  \rho_mn'+- = R_mm'' (F_2)_m''q (V_2*)_ql (D^-1)_lk (V_1^T)_kp (G_1^T)_pn'
        !  \rho_nm'-+ = R_nn'' (G_2)_n''q (V_2*)_ql (D^-1)_lk (V_1^T)_kp (F_1^T)_pm'
        !  where m and m' are the spherical harmonic basis indices of the large component,
        !        n and n' are the spherical harmonic basis indices of the small component.
        !
        ! Note:
        !  1) calculate_Rotated_Wavefunction has been called to calculate:
        !       F(Omega)_mk = \sum_{m'} R(Omega)_mm' F_{m'k}
        !       G(Omega)_nk = \sum_{n'} R(Omega)_nn' G_{n'k}
        !  2) V_1 and V_2 are in the canonical basis, 
        !     i.e.
        !         ( ...             )
        !    V_1 =(     0    v_p    )
        !         (    -v_p   0     )
        !         (             ... )
        ! Returns:
        !  * rho_mm(:,:,ifg,iphi,it): ifg=1: \rho_mm'++, ifg=2: \rho_nn'--, ifg=3: \rho_mn'+-, ifg=4: \rho_nm'-+
        !  * prho_mm(:,:,ifg,iphi,it)
        !--------------------------------------------------------------------------------------
        integer :: dim_m_max,truncated_dim,ifg,i0sp,ndsp,nfgsp,l,truncated_l,m1,k,truncated_k,nl1,nj1,nm1,nm1_reversed, &
                    m1_reversed, i0f,i0g,nf,ng,m2,nl2,nj2,nm2,nm2_reversed,m2_reversed
        complex(r64) :: vl
        real(r64) :: ul,uk,vk,nj1_half,nm1_half,nm1_reversed_half,nj2_half,nm2_half,nm2_reversed_half
        complex(r64), dimension(:,:),allocatable :: A,pA,B,pB,C,pC,tmp,F,pF,RF,pRF,G,pG,RG,pRG
        if (wf1%truncated_dim(it) /=wf2%truncated_dim(it) ) stop "Wrong truncated dimension!"
        truncated_dim = wf1%truncated_dim(it)
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        allocate(A(dim_m_max,2*truncated_dim))
        allocate(pA(dim_m_max,2*truncated_dim))
        allocate(B(2*truncated_dim,2*truncated_dim))
        allocate(pB(2*truncated_dim,2*truncated_dim))
        allocate(C(2*truncated_dim,dim_m_max))
        allocate(pC(2*truncated_dim,dim_m_max))
        allocate(tmp(dim_m_max,dim_m_max))
        allocate(F(2*truncated_dim,BS%HO_sph%idsp(1,1)))
        allocate(pF(2*truncated_dim,BS%HO_sph%idsp(1,1)))
        allocate(G(2*truncated_dim,BS%HO_sph%idsp(1,2)))
        allocate(pG(2*truncated_dim,BS%HO_sph%idsp(1,2)))
        allocate(RF(BS%HO_sph%idsp(1,1),2*truncated_dim))
        allocate(pRF(BS%HO_sph%idsp(1,1),2*truncated_dim))
        allocate(RG(BS%HO_sph%idsp(1,2),2*truncated_dim))
        allocate(pRG(BS%HO_sph%idsp(1,2),2*truncated_dim))
        ! rho^++ and rho^--
        do ifg = 1, 2
            i0sp = BS%HO_sph%iasp(1,ifg)
            ndsp = BS%HO_sph%idsp(1,ifg)
            nfgsp = (ifg-1)*BS%HO_sph%idsp(1,1)

            do l = 1, wf2%nk(it) ! loop over single particle basis
                if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
                vl = dsqrt(wf2%v2(l,it))*ei2phi
                ! ul = dsqrt(1-wf2%v2(l,it))
                truncated_l = wf2%truncated_index(l,it)
                do m1 = 1,ndsp
                    A(m1,truncated_l) = fg_rotated(nfgsp+m1,l, it)*vl
                    A(m1,truncated_dim+truncated_l) = fg_rotated(nfgsp+m1,-l, it)*vl
                    pA(m1,truncated_l) = pfg_rotated(nfgsp+m1,l, it)*vl
                    pA(m1,truncated_dim+truncated_l) = pfg_rotated(nfgsp+m1,-l, it)*vl  
                end do
                do k = 1, wf1%nk(it)
                    if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                    ! vk = dsqrt(wf1%v2(k,it))
                    ! uk = dsqrt(1-wf1%v2(k,it))
                    truncated_k = wf1%truncated_index(k,it)
                    B(truncated_l,truncated_k) = D_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    B(truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_dim+truncated_l,truncated_k,it)
                    B(truncated_dim+truncated_l,truncated_k) = -D_inv(truncated_l,truncated_dim+truncated_k,it)
                    B(truncated_dim+truncated_l,truncated_dim+truncated_k) = D_inv(truncated_l,truncated_k,it)
                    pB(truncated_l,truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    pB(truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_dim+truncated_l,truncated_k,it)
                    pB(truncated_dim+truncated_l,truncated_k) = -pD_inv(truncated_l,truncated_dim+truncated_k,it)
                    pB(truncated_dim+truncated_l,truncated_dim+truncated_k) = pD_inv(truncated_l,truncated_k,it)
                end do
            end do
            do k = 1, wf1%nk(it)
                if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                vk = dsqrt(wf1%v2(k,it))
                ! uk = dsqrt(1-wf1%v2(k,it))
                truncated_k = wf1%truncated_index(k,it)
                do m1 = 1,ndsp !|m1>
                    ! nr1= BS%HO_sph%nljm(i0sp+m1,1) ! n_r
                    nl1= BS%HO_sph%nljm(i0sp+m1,2) ! l
                    nj1= BS%HO_sph%nljm(i0sp+m1,3) ! j +1/2
                    nm1= BS%HO_sph%nljm(i0sp+m1,4) ! m_j + 1/2
                    nj1_half = nj1 - 0.5 ! j
                    nm1_half = nm1 - 0.5 ! m_j
                    ! nm1_reversed = 1 - nm1 ! -m_j + 1/2
                    ! nm1_reversed_half = nm1_reversed - 0.5 ! -m_j
                    m1_reversed = m1 - 2*nm1 + 1 ! location of |-m1>
                    C(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*vk
                    C(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*vk*(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                    pC(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*vk
                    pC(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*vk*(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                end do
            end do

            call cmatmul_ABC(A(1:ndsp,1:2*truncated_dim), B(1:2*truncated_dim,1:2*truncated_dim), C(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp) ! tmp = ABC
            mix%rho_mm(:,:,ifg,iphi,it) = tmp ! ifg = 1 is rho^++ and ifg = 2 is rho^--
            call cmatmul_ABC(pA(1:ndsp,1:2*truncated_dim), pB(1:2*truncated_dim,1:2*truncated_dim), pC(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp) ! tmp = ABC
            mix%prho_mm(:,:,ifg,iphi,it) = tmp
        end do
        ! rho^+- and rho^-+
        i0f = BS%HO_sph%iasp(1,1)
        i0g = BS%HO_sph%iasp(1,2)
        nf = BS%HO_sph%idsp(1,1)
        ng = BS%HO_sph%idsp(1,2)
        do l = 1, wf2%nk(it)
            if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
            vl = dsqrt(wf2%v2(l,it))*ei2phi
            ! ul = dsqrt(1-wf2%v2(l,it))
            truncated_l = wf2%truncated_index(l,it)
            do m1 = 1, nf
                RF(m1,truncated_l) = fg_rotated(m1,l, it)*vl
                RF(m1,truncated_dim+truncated_l) = fg_rotated(m1,-l, it)*vl
                pRF(m1,truncated_l) = pfg_rotated(m1,l, it)*vl
                pRF(m1,truncated_dim+truncated_l) = pfg_rotated(m1,-l, it)*vl
            end do
            do m2 = 1, ng
                RG(m2,truncated_l) = fg_rotated(nf+m2,l, it)*vl
                RG(m2,truncated_dim+truncated_l) = fg_rotated(nf+m2,-l, it)*vl
                pRG(m2,truncated_l) = pfg_rotated(nf+m2,l, it)*vl
                pRG(m2,truncated_dim+truncated_l) = pfg_rotated(nf+m2,-l, it)*vl
            end do
            do k = 1, wf1%nk(it)
                if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                ! vk = dsqrt(wf1%v2(k,it))
                ! uk = dsqrt(1-wf1%v2(k,it))
                truncated_k = wf1%truncated_index(k,it)
                B(truncated_l,truncated_k) = D_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                B(truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_dim+truncated_l,truncated_k,it)
                B(truncated_dim+truncated_l,truncated_k) = -D_inv(truncated_l,truncated_dim+truncated_k,it)
                B(truncated_dim+truncated_l,truncated_dim+truncated_k) = D_inv(truncated_l,truncated_k,it)
                pB(truncated_l,truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                pB(truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_dim+truncated_l,truncated_k,it)
                pB(truncated_dim+truncated_l,truncated_k) = -pD_inv(truncated_l,truncated_dim+truncated_k,it)
                pB(truncated_dim+truncated_l,truncated_dim+truncated_k) = pD_inv(truncated_l,truncated_k,it)
            end do
        end do 
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
            vk = dsqrt(wf1%v2(k,it))
            ! uk = dsqrt(1-wf1%v2(k,it))
            truncated_k = wf1%truncated_index(k,it)
            do m1 = 1, nf
                ! nr1= BS%HO_sph%nljm(i0sp+m1,1) ! n_r
                nl1= BS%HO_sph%nljm(i0f+m1,2) ! l
                nj1= BS%HO_sph%nljm(i0f+m1,3) ! j +1/2
                nm1= BS%HO_sph%nljm(i0f+m1,4) ! m_j + 1/2
                nj1_half = nj1 - 0.5 ! j
                nm1_half = nm1 - 0.5 ! m_j
                ! nm1_reversed = 1 - nm1 ! -m_j + 1/2
                ! nm1_reversed_half = nm1_reversed - 0.5 ! -m_j
                m1_reversed = m1 - 2*nm1 + 1 ! location of |-m1>
                F(truncated_k,m1) = wf1%fg(m1,k,it)*vk
                F(truncated_dim+truncated_k,m1) = wf1%fg(m1_reversed,k,it)*vk*(-1)**(Int(nl1+nj1_half+nm1_half))
                pF(truncated_k,m1) = wf1%fg(+m1,k,it)*vk
                pF(truncated_dim+truncated_k,m1) = wf1%fg(m1_reversed,k,it)*vk*(-1)**(Int(nl1+nj1_half+nm1_half))
            end do
            do m2 = 1, ng
                nl2= BS%HO_sph%nljm(i0g+m2,2) ! l
                nj2= BS%HO_sph%nljm(i0g+m2,3) ! j +1/2
                nm2= BS%HO_sph%nljm(i0g+m2,4) ! m_j + 1/2
                nj2_half = nj2 - 0.5 ! j
                nm2_half = nm2 - 0.5 ! m_j
                m2_reversed = m2 - 2*nm2 + 1 ! location of |-m2>
                G(truncated_k,m2) = wf1%fg(nf+m2,k,it)*vk
                G(truncated_dim+truncated_k,m2)  = - wf1%fg(nf+m2_reversed,k,it)*vk*(-1)**(Int(nl2+nj2_half+nm2_half))
                pG(truncated_k,m2) = wf1%fg(nf+m2,k,it)*vk
                pG(truncated_dim+truncated_k,m2) = - wf1%fg(nf+m2_reversed,k,it)*vk*(-1)**(Int(nl2+nj2_half+nm2_half))
            end do
        end do 
        call cmatmul_ABC(RF(1:nf,1:2*truncated_dim), B(1:2*truncated_dim,1:2*truncated_dim), G(1:2*truncated_dim,1:ng), & 
                        tmp(1:nf,1:ng), nf, 2*truncated_dim, 2*truncated_dim, ng)
        mix%rho_mm(:,:,3,iphi,it) = tmp ! ifg = 3 is rho^+-
        call cmatmul_ABC(pRF(1:nf,1:2*truncated_dim), pB(1:2*truncated_dim,1:2*truncated_dim), pG(1:2*truncated_dim,1:ng),&
                        tmp(1:nf,1:ng), nf, 2*truncated_dim, 2*truncated_dim, ng)
        mix%prho_mm(:,:,3,iphi,it) = tmp
        call cmatmul_ABC(RG(1:ng,1:2*truncated_dim), B(1:2*truncated_dim,1:2*truncated_dim), F(1:2*truncated_dim,1:nf), &
                        tmp(1:ng,1:nf), ng, 2*truncated_dim, 2*truncated_dim, nf)
        mix%rho_mm(:,:,4,iphi,it) = tmp ! ifg = 4 is rho^-+
        call cmatmul_ABC(pRG(1:ng,1:2*truncated_dim), pB(1:2*truncated_dim,1:2*truncated_dim), pF(1:2*truncated_dim,1:nf), &
                        tmp(1:ng,1:nf), ng, 2*truncated_dim, 2*truncated_dim, nf)
        mix%prho_mm(:,:,4,iphi,it) = tmp
    end subroutine
    subroutine mixed_tensor_matrix_elements
        !--------------------------------------------------------------------------------------
        !  calculate mixed tensor matrix elements:
        !   kappa_mm'  = <q1|c_m' c_m R|q2>/<q1|R|q2>.
        !   kappa_mm'* = <q1|c+_m' c+_m R|q2>/<q1|R|q2>.
        !
        !  \kappa_mm'++  = \sum_{k,l>0; k,l<0} \sqrt{f2_l} \sqrt{f1_k} v_l u_k 
        !                [R_mm'' (F_2)_m''l (D^-1)_\bar{l}\bar{k} (F_1^T)_km' - R_mm'' (F_2)_m''\bar{l} (D^-1)_l\bar{k} (F_1^T)_km']
        !  \kappa_nn'--  = \sum_{k,l>0; k,l<0} \sqrt{f2_l} \sqrt{f1_k} v_l u_k 
        !                [R_nn'' (G_2)_n''l (D^-1)_\bar{l}\bar{k} (G_1^T)_kn' - R_nn'' (G_2)_n''\bar{l} (D^-1)_l\bar{k} (G_1^T)_kn']
        !  \kappa*_mm'++ = \sum_{k,l>0; k,l<0} \sqrt{f2_l} \sqrt{f1_k} u_l v_k 
        !                [R*_mm'' (F_2)_m''l (D^-1)_\bar{l}\bar{k} (F_1^T)_km' - R*_mm'' (F_2)_m''\bar{l} (D^-1)_l\bar{k} (F_1^T)_km']
        !  \kappa*_nn'-- = \sum_{k,l>0; k,l<0} \sqrt{f2_l} \sqrt{f1_k} u_l v_k 
        !                [R*_nn'' (G_2)_n''l (D^-1)_\bar{l}\bar{k} (G_1^T)_kn' - R*_nn'' (G_2)_n''\bar{l} (D^-1)_l\bar{k} (G_1^T)_kn']
        !
        !  where m and m' are the spherical harmonic basis indices of the large component,
        !        n and n' are the spherical harmonic basis indices of the small component.
        !
        !
        ! Returns:
        !  * kappa_mm(:,:,ifg,iphi,it)  : ifg=1: \kappa_mm'++ ,  ifg=2: \kappa_nn'--
        !  * kappac_mm(:,:,ifg,iphi,it) : ifg=1: \kappa*_mm'++,  ifg=2: \kappa*_nn'--
        !  * pkappa_mm(:,:,ifg,iphi,it)
        !  * pkappac_mm(:,:,ifg,iphi,it)
        !--------------------------------------------------------------------------------------
        integer :: truncated_dim,dim_m_max,ifg,i0sp,ndsp,nfgsp,l,truncated_l,m1,k,truncated_k,nl1,nj1,nm1,nm1_reversed, &
                    m1_reversed
        complex(r64) :: vl
        real(r64) :: f_l,ul,uk,vk,f_k,nj1_half,nm1_half
        complex(r64), dimension(:,:),allocatable :: Av,pAv,Au,pAu,B,pB,Cu,pCu,Cv,pCv,tmp
        if (wf1%truncated_dim(it) /=wf2%truncated_dim(it) ) stop "Wrong truncated dimension!"
        truncated_dim = wf1%truncated_dim(it)
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        allocate(Av(dim_m_max,2*truncated_dim))
        allocate(pAv(dim_m_max,2*truncated_dim))
        allocate(Au(dim_m_max,2*truncated_dim))
        allocate(pAu(dim_m_max,2*truncated_dim))
        allocate(B(2*truncated_dim,2*truncated_dim))
        allocate(pB(2*truncated_dim,2*truncated_dim))
        allocate(Cu(2*truncated_dim,dim_m_max))
        allocate(pCu(2*truncated_dim,dim_m_max))
        allocate(Cv(2*truncated_dim,dim_m_max))
        allocate(pCv(2*truncated_dim,dim_m_max))
        allocate(tmp(dim_m_max,dim_m_max))
        do ifg = 1, 2
            i0sp = BS%HO_sph%iasp(1,ifg)
            ndsp = BS%HO_sph%idsp(1,ifg)
            nfgsp = (ifg-1)*BS%HO_sph%idsp(1,1)
            do l = 1, wf2%nk(it) ! loop over single particle basis
                if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
                vl = sqrt(wf2%v2(l,it))*ei2phi
                ul = dsqrt(1-wf2%v2(l,it))
                f_l = wf2%skk(l,it)
                truncated_l = wf2%truncated_index(l,it)
                do m1 = 1, ndsp
                    Av(m1,truncated_l) = fg_rotated(nfgsp+m1,l,it)*vl*dsqrt(f_l)
                    Av(m1,truncated_dim+truncated_l) = fg_rotated(nfgsp+m1,-l, it)*vl*dsqrt(f_l)
                    pAv(m1,truncated_l) = pfg_rotated(nfgsp+m1,l,it)*vl*dsqrt(f_l)
                    pAv(m1,truncated_dim+truncated_l) = pfg_rotated(nfgsp+m1,-l, it)*vl*dsqrt(f_l)
                    Au(m1,truncated_l) = DCONJG(fg_rotated(nfgsp+m1,l,it))*ul*dsqrt(f_l)
                    Au(m1,truncated_dim+truncated_l) = DCONJG(fg_rotated(nfgsp+m1,-l,it))*ul*dsqrt(f_l)
                    pAu(m1,truncated_l) = DCONJG(pfg_rotated(nfgsp+m1,l,it))*ul*dsqrt(f_l)
                    pAu(m1,truncated_dim+truncated_l) = DCONJG(pfg_rotated(nfgsp+m1,-l, it))*ul*dsqrt(f_l)
                end do 
                do k = 1, wf1%nk(it)
                    if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                    truncated_k = wf1%truncated_index(k,it)
                    B(truncated_l,truncated_k) = D_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    B(truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_dim+truncated_l,truncated_k,it)
                    B(truncated_dim+truncated_l,truncated_k) = -D_inv(truncated_l,truncated_dim+truncated_k,it)
                    B(truncated_dim+truncated_l,truncated_dim+truncated_k) = D_inv(truncated_l,truncated_k,it)
                    pB(truncated_l,truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    pB(truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_dim+truncated_l,truncated_k,it)
                    pB(truncated_dim+truncated_l,truncated_k) = -pD_inv(truncated_l,truncated_dim+truncated_k,it)
                    pB(truncated_dim+truncated_l,truncated_dim+truncated_k) = pD_inv(truncated_l,truncated_k,it)
                end do
            end do
            do k = 1, wf1%nk(it)
                if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                vk = dsqrt(wf1%v2(k,it))
                uk = dsqrt(1-wf1%v2(k,it))
                f_k = wf1%skk(k,it)
                truncated_k = wf1%truncated_index(k,it)
                do m1 = 1,ndsp !|m1>
                    ! nr1= BS%HO_sph%nljm(i0sp+m1,1) ! n_r
                    nl1= BS%HO_sph%nljm(i0sp+m1,2) ! l
                    nj1= BS%HO_sph%nljm(i0sp+m1,3) ! j +1/2
                    nm1= BS%HO_sph%nljm(i0sp+m1,4) ! m_j + 1/2
                    nj1_half = nj1 - 0.5 ! j
                    nm1_half = nm1 - 0.5 ! m_j
                    ! nm1_reversed = 1 - nm1 ! -m_j + 1/2
                    ! nm1_reversed_half = nm1_reversed - 0.5 ! -m_j
                    m1_reversed = m1 - 2*nm1 + 1 ! location of |-m1>
                    Cu(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*uk*dsqrt(f_k)
                    Cu(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*uk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1)) ! F_{m \bar{k}} =  F_{-mk} (-1)^{l+j+m_j}; G_{n \bar{k}}= -(-1)^{l+j+m_j} 
                    pCu(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*uk*dsqrt(f_k)
                    pCu(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*uk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                    Cv(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*vk*dsqrt(f_k)
                    Cv(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*vk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                    pCv(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*vk*dsqrt(f_k)
                    pCv(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*vk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                end do
            end do
            call cmatmul_ABC(Av(1:ndsp,1:2*truncated_dim), B(1:2*truncated_dim,1:2*truncated_dim), Cu(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            mix%kappa_mm(:,:,ifg,iphi,it) = tmp ! ifg = 1 is kappa^++ and ifg = 2 is kappa^--
            call cmatmul_ABC(pAv(1:ndsp,1:2*truncated_dim), pB(1:2*truncated_dim,1:2*truncated_dim), pCu(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            mix%pkappa_mm(:,:,ifg,iphi,it) = tmp

            call cmatmul_ABC(Au(1:ndsp,1:2*truncated_dim), B(1:2*truncated_dim,1:2*truncated_dim), Cv(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            mix%kappac_mm(:,:,ifg,iphi,it) = tmp ! ifg = 1 is kappa*^++ and ifg = 2 is kappa*^--
            call cmatmul_ABC(pAu(1:ndsp,1:2*truncated_dim), pB(1:2*truncated_dim,1:2*truncated_dim), pCv(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            mix%pkappac_mm(:,:,ifg,iphi,it) = tmp
        end do
    end subroutine
    subroutine mixed_kappa_matrix_elements(with_f)
        !--------------------------------------------------------------------------------------
        !  calculate mixed tensor matrix elements:
        !   kappa10_mm'  = <q1|c_m' c_m R|q2>/<q1|R|q2>.
        !   kappa01_m'm* = <q1|c+_m' c+_m R|q2>/<q1|R|q2>.
        !
        !  \kappa10_mm'++  = \sqrt{f_1}\sqrt{f_2} R_mm'' (F_2)_m''q (V_2*)_ql (D^-1)_lk (U_1^T)_kp (F_1^T)_pm'
        !  \kappa10_nn'--  = \sqrt{f_1}\sqrt{f_2} R_nn'' (G_2)_n''q (V_2*)_ql (D^-1)_lk (U_1^T)_kp (G_1^T)_pn'
        !  \kappa10_mn+-  = \sqrt{f_1}\sqrt{f_2} R_mm'' (F_2)_m''q (V_2*)_ql (D^-1)_lk (U_1^T)_kp (G_1^T)_pn
        !  \kappa10_nm-+  = \sqrt{f_1}\sqrt{f_2} R_nn'' (G_2)_n''q (V_2*)_ql (D^-1)_lk (U_1^T)_kp (F_1^T)_pm
        !  \kappa01_m'm++* = \sqrt{f_1}\sqrt{f_2} R*_mm'' (F_2)_m''q (U_2*)_ql (D^-1)_lk (V_1^T)_kp (F_1^T)_pm'
        !  \kappa01_n'n--* = \sqrt{f_1}\sqrt{f_2} R*_nn'' (G_2)_n''q (U_2*)_ql (D^-1)_lk (V_1^T)_kp (G_1^T)_pn'
        !  \kappa01_mn+-* = \sqrt{f_1}\sqrt{f_2} R*_nn'' (G_2)_n''q (U_2*)_ql (D^-1)_lk (V_1^T)_kp (F_1^T)_pm
        !  \kappa01_nm-+* = \sqrt{f_1}\sqrt{f_2} R*_mm'' (F_2)_m''q (U_2*)_ql (D^-1)_lk (V_1^T)_kp (G_1^T)_pn
        !  where m and m' are the spherical harmonic basis indices of the large component,
        !        n and n' are the spherical harmonic basis indices of the small component.
        !
        ! Note:
        !  1) calculate_Rotated_Wavefunction has been called to calculate:
        !       F(Omega)_mk = \sum_{m'} R(Omega)_mm' F_{m'k}
        !       G(Omega)_nk = \sum_{n'} R(Omega)_nn' G_{n'k}
        !  2) U_1, U_2 and V_1, V_2 are in the canonical basis, 
        !     i.e.
        !         ( ...             )        ( ...             )
        !    U_1 =(     u_p    0    )   V_1 =(     0    v_p    )
        !         (      0    u_p   )        (    -v_p   0     )
        !         (             ... )        (             ... )
        ! Returns:
        !  * kappa10_mm(:,:,ifg,iphi,it): ifg=1: \kappa10_mm'++, ifg=2: \kappa10_nn'--, ifg=3: \kappa10_mn^+-, ifg=4: \kappa10_nm^-+
        !  * kappa01c_mm(:,:,ifg,iphi,it): ifg=1: \kappa01_m'm++*, ifg=2: \kappa01_n'n--*, ifg=3: \kappa01_mn^+-*, ifg=4: \kappa01_nm^-+*
        !  * pkappa10_mm(:,:,ifg,iphi,it)
        !  * pkappa01c_mm(:,:,ifg,iphi,it)
        !--------------------------------------------------------------------------------------
        integer :: truncated_dim,dim_m_max,ifg,i0sp,ndsp,nfgsp,l,truncated_l,m1,k,truncated_k,nl1,nj1,nm1,nm1_reversed, &
                    m1_reversed,m2,i0f,i0g,nf,ng,nl2,nj2,nm2,nm2_reversed,m2_reversed
        complex(r64) :: vl
        real(r64) :: f_l,ul,uk,vk,f_k,nj1_half,nm1_half,nj2_half,nm2_half
        complex(r64), dimension(:,:),allocatable :: Av,pAv,Au,pAu,B10,pB10,B01,pB01,Cu,pCu,Cv,pCv,tmp,RFv,pRFv,RFu,pRFu,RGv, &
                                                    pRGv,RGu,pRGu,Fu,pFu,Fv,pFv,Gu,pGu,Gv,pGv
        logical :: with_f
        if (wf1%truncated_dim(it) /=wf2%truncated_dim(it) ) stop "Wrong truncated dimension!"
        truncated_dim = wf1%truncated_dim(it)
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        allocate(Av(dim_m_max,2*truncated_dim))
        allocate(pAv(dim_m_max,2*truncated_dim))
        allocate(Au(dim_m_max,2*truncated_dim))
        allocate(pAu(dim_m_max,2*truncated_dim))
        allocate(B10(2*truncated_dim,2*truncated_dim))
        allocate(pB10(2*truncated_dim,2*truncated_dim))
        allocate(B01(2*truncated_dim,2*truncated_dim))
        allocate(pB01(2*truncated_dim,2*truncated_dim))
        allocate(Cu(2*truncated_dim,dim_m_max))
        allocate(pCu(2*truncated_dim,dim_m_max))
        allocate(Cv(2*truncated_dim,dim_m_max))
        allocate(pCv(2*truncated_dim,dim_m_max))
        allocate(tmp(dim_m_max,dim_m_max))
        allocate(RFv(BS%HO_sph%idsp(1,1),2*truncated_dim))
        allocate(pRFv(BS%HO_sph%idsp(1,1),2*truncated_dim))
        allocate(RFu(BS%HO_sph%idsp(1,1),2*truncated_dim))
        allocate(pRFu(BS%HO_sph%idsp(1,1),2*truncated_dim))
        allocate(RGv(BS%HO_sph%idsp(1,2),2*truncated_dim))
        allocate(pRGv(BS%HO_sph%idsp(1,2),2*truncated_dim))
        allocate(RGu(BS%HO_sph%idsp(1,2),2*truncated_dim))
        allocate(pRGu(BS%HO_sph%idsp(1,2),2*truncated_dim))
        allocate(Fu(2*truncated_dim,BS%HO_sph%idsp(1,1)))
        allocate(pFu(2*truncated_dim,BS%HO_sph%idsp(1,1)))
        allocate(Fv(2*truncated_dim,BS%HO_sph%idsp(1,1)))
        allocate(pFv(2*truncated_dim,BS%HO_sph%idsp(1,1)))
        allocate(Gu(2*truncated_dim,BS%HO_sph%idsp(1,2)))
        allocate(pGu(2*truncated_dim,BS%HO_sph%idsp(1,2)))
        allocate(Gv(2*truncated_dim,BS%HO_sph%idsp(1,2)))
        allocate(pGv(2*truncated_dim,BS%HO_sph%idsp(1,2)))
        ! \kappa10_mm'^++ and \kappa10_nn'^--; \kappa01_m'm^++* and \kappa01_n'n^--*
        do ifg = 1, 2
            i0sp = BS%HO_sph%iasp(1,ifg)
            ndsp = BS%HO_sph%idsp(1,ifg)
            nfgsp = (ifg-1)*BS%HO_sph%idsp(1,1)
            do l = 1, wf2%nk(it) ! loop over single particle basis
                if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
                vl = sqrt(wf2%v2(l,it))*ei2phi
                ul = dsqrt(1-wf2%v2(l,it))
                if(with_f) then 
                    f_l = wf2%skk(l,it)
                else 
                    f_l = 1.d0
                end if 
                truncated_l = wf2%truncated_index(l,it)
                do m1 = 1, ndsp
                    Av(m1,truncated_l) = fg_rotated(nfgsp+m1,l,it)*vl*dsqrt(f_l)
                    Av(m1,truncated_dim+truncated_l) = fg_rotated(nfgsp+m1,-l, it)*vl*dsqrt(f_l)
                    pAv(m1,truncated_l) = pfg_rotated(nfgsp+m1,l,it)*vl*dsqrt(f_l)
                    pAv(m1,truncated_dim+truncated_l) = pfg_rotated(nfgsp+m1,-l, it)*vl*dsqrt(f_l)
                    Au(m1,truncated_l) = DCONJG(fg_rotated(nfgsp+m1,l,it))*ul*dsqrt(f_l)
                    Au(m1,truncated_dim+truncated_l) = DCONJG(fg_rotated(nfgsp+m1,-l,it))*ul*dsqrt(f_l)
                    pAu(m1,truncated_l) = DCONJG(pfg_rotated(nfgsp+m1,l,it))*ul*dsqrt(f_l)
                    pAu(m1,truncated_dim+truncated_l) = DCONJG(pfg_rotated(nfgsp+m1,-l, it))*ul*dsqrt(f_l)
                end do 
                do k = 1, wf1%nk(it)
                    if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                    truncated_k = wf1%truncated_index(k,it)
                    B10(truncated_l,truncated_k) = D_inv(truncated_dim+truncated_l,truncated_k,it)
                    B10(truncated_l,truncated_dim+truncated_k) = D_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    B10(truncated_dim+truncated_l,truncated_k) = -D_inv(truncated_l,truncated_k,it)
                    B10(truncated_dim+truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_l,truncated_dim+truncated_k,it)
                    pB10(truncated_l,truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_k,it)
                    pB10(truncated_l,truncated_dim+truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    pB10(truncated_dim+truncated_l,truncated_k) = -pD_inv(truncated_l,truncated_k,it)
                    pB10(truncated_dim+truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_l,truncated_dim+truncated_k,it)
                    !
                    B01(truncated_l,truncated_k) = D_inv(truncated_l,truncated_dim+truncated_k,it)
                    B01(truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_l,truncated_k,it)
                    B01(truncated_dim+truncated_l,truncated_k) = D_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    B01(truncated_dim+truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_dim+truncated_l,truncated_k,it)
                    pB01(truncated_l,truncated_k) = pD_inv(truncated_l,truncated_dim+truncated_k,it)
                    pB01(truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_l,truncated_k,it)
                    pB01(truncated_dim+truncated_l,truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                    pB01(truncated_dim+truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_dim+truncated_l,truncated_k,it)
                end do
            end do
            do k = 1, wf1%nk(it)
                if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                vk = dsqrt(wf1%v2(k,it))
                uk = dsqrt(1-wf1%v2(k,it))
                if(with_f) then 
                    f_k = wf1%skk(k,it)
                else 
                    f_k = 1.d0
                end if 
                truncated_k = wf1%truncated_index(k,it)
                do m1 = 1,ndsp !|m1>
                    ! nr1= BS%HO_sph%nljm(i0sp+m1,1) ! n_r
                    nl1= BS%HO_sph%nljm(i0sp+m1,2) ! l
                    nj1= BS%HO_sph%nljm(i0sp+m1,3) ! j +1/2
                    nm1= BS%HO_sph%nljm(i0sp+m1,4) ! m_j + 1/2
                    nj1_half = nj1 - 0.5 ! j
                    nm1_half = nm1 - 0.5 ! m_j
                    ! nm1_reversed = 1 - nm1 ! -m_j + 1/2
                    ! nm1_reversed_half = nm1_reversed - 0.5 ! -m_j
                    m1_reversed = m1 - 2*nm1 + 1 ! location of |-m1>
                    Cu(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*uk*dsqrt(f_k)
                    Cu(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*uk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                    pCu(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*uk*dsqrt(f_k)
                    pCu(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*uk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                    Cv(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*vk*dsqrt(f_k)
                    Cv(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*vk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                    pCv(truncated_k,m1) = wf1%fg(nfgsp+m1,k,it)*vk*dsqrt(f_k)
                    pCv(truncated_dim+truncated_k,m1) = wf1%fg(nfgsp+m1_reversed,k,it)*vk*dsqrt(f_k) &
                                                        *(-1)**(Int(nl1+nj1_half+nm1_half+ifg-1))
                end do
            end do
            ! kappa10_mm
            call cmatmul_ABC(Av(1:ndsp,1:2*truncated_dim), B10(1:2*truncated_dim,1:2*truncated_dim), Cu(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            mix%kappa10_mm(:,:,ifg,iphi,it) = tmp ! ifg = 1 is kappa10^++ and ifg = 2 is kappa10^--
            call cmatmul_ABC(pAv(1:ndsp,1:2*truncated_dim), pB10(1:2*truncated_dim,1:2*truncated_dim), pCu(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            mix%pkappa10_mm(:,:,ifg,iphi,it) = tmp
            
            ! kappa01c_mm
            call cmatmul_ABC(Au(1:ndsp,1:2*truncated_dim), B01(1:2*truncated_dim,1:2*truncated_dim), Cv(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            do m1 = 1, ndsp
                do m2 = 1,ndsp
                    mix%kappa01c_mm(m1,m2,ifg,iphi,it) = tmp(m2,m1) ! ifg = 1 is kappa01^++* and ifg = 2 is kappa01^--*
                end do
            end do
            call cmatmul_ABC(pAu(1:ndsp,1:2*truncated_dim), pB01(1:2*truncated_dim,1:2*truncated_dim), pCv(1:2*truncated_dim,1:ndsp), tmp(1:ndsp,1:ndsp), &
                                ndsp, 2*truncated_dim, 2*truncated_dim, ndsp)
            do m1 = 1, ndsp
                do m2 = 1,ndsp
                    mix%pkappa01c_mm(m1,m2,ifg,iphi,it) = tmp(m2,m1)
                end do
            end do
        end do
        ! \kappa10_mn^+- and \kappa10_nm^-+; \kappa01_mn^+-* and \kappa01_nm^-+*
        i0f = BS%HO_sph%iasp(1,1)
        i0g = BS%HO_sph%iasp(1,2)
        nf = BS%HO_sph%idsp(1,1)
        ng = BS%HO_sph%idsp(1,2)
        do l = 1, wf2%nk(it) ! loop over single particle basis
            if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
            vl = sqrt(wf2%v2(l,it))*ei2phi
            ul = dsqrt(1-wf2%v2(l,it))
            if(with_f) then 
                f_l = wf2%skk(l,it)
            else 
                f_l = 1.d0
            end if 
            truncated_l = wf2%truncated_index(l,it)
            do m1 = 1, nf
                RFv(m1,truncated_l) = fg_rotated(m1,l,it)*vl*dsqrt(f_l)
                RFv(m1,truncated_dim+truncated_l) = fg_rotated(m1,-l, it)*vl*dsqrt(f_l)
                pRFv(m1,truncated_l) = pfg_rotated(m1,l,it)*vl*dsqrt(f_l)
                pRFv(m1,truncated_dim+truncated_l) = pfg_rotated(m1,-l, it)*vl*dsqrt(f_l)
                RFu(m1,truncated_l) = DCONJG(fg_rotated(m1,l,it))*ul*dsqrt(f_l)
                RFu(m1,truncated_dim+truncated_l) = DCONJG(fg_rotated(m1,-l,it))*ul*dsqrt(f_l)
                pRFu(m1,truncated_l) = DCONJG(pfg_rotated(m1,l,it))*ul*dsqrt(f_l)
                pRFu(m1,truncated_dim+truncated_l) = DCONJG(pfg_rotated(m1,-l, it))*ul*dsqrt(f_l)
            end do 
            do m2 = 1, ng
                RGv(m2,truncated_l) = fg_rotated(nf+m2,l,it)*vl*dsqrt(f_l)
                RGv(m2,truncated_dim+truncated_l) = fg_rotated(nf+m2,-l, it)*vl*dsqrt(f_l)
                pRGv(m2,truncated_l) = pfg_rotated(nf+m2,l,it)*vl*dsqrt(f_l)
                pRGv(m2,truncated_dim+truncated_l) = pfg_rotated(nf+m2,-l, it)*vl*dsqrt(f_l)
                RGu(m2,truncated_l) = DCONJG(fg_rotated(nf+m2,l,it))*ul*dsqrt(f_l)
                RGu(m2,truncated_dim+truncated_l) = DCONJG(fg_rotated(nf+m2,-l,it))*ul*dsqrt(f_l)
                pRGu(m2,truncated_l) = DCONJG(pfg_rotated(nf+m2,l,it))*ul*dsqrt(f_l)
                pRGu(m2,truncated_dim+truncated_l) = DCONJG(pfg_rotated(nf+m2,-l, it))*ul*dsqrt(f_l)
            end do 

            do k = 1, wf1%nk(it)
                if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
                truncated_k = wf1%truncated_index(k,it)
                B10(truncated_l,truncated_k) = D_inv(truncated_dim+truncated_l,truncated_k,it)
                B10(truncated_l,truncated_dim+truncated_k) = D_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                B10(truncated_dim+truncated_l,truncated_k) = -D_inv(truncated_l,truncated_k,it)
                B10(truncated_dim+truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_l,truncated_dim+truncated_k,it)
                pB10(truncated_l,truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_k,it)
                pB10(truncated_l,truncated_dim+truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                pB10(truncated_dim+truncated_l,truncated_k) = -pD_inv(truncated_l,truncated_k,it)
                pB10(truncated_dim+truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_l,truncated_dim+truncated_k,it)
                !
                B01(truncated_l,truncated_k) = D_inv(truncated_l,truncated_dim+truncated_k,it)
                B01(truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_l,truncated_k,it)
                B01(truncated_dim+truncated_l,truncated_k) = D_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                B01(truncated_dim+truncated_l,truncated_dim+truncated_k) = -D_inv(truncated_dim+truncated_l,truncated_k,it)
                pB01(truncated_l,truncated_k) = pD_inv(truncated_l,truncated_dim+truncated_k,it)
                pB01(truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_l,truncated_k,it)
                pB01(truncated_dim+truncated_l,truncated_k) = pD_inv(truncated_dim+truncated_l,truncated_dim+truncated_k,it)
                pB01(truncated_dim+truncated_l,truncated_dim+truncated_k) = -pD_inv(truncated_dim+truncated_l,truncated_k,it)
            end do
        end do
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
            vk = dsqrt(wf1%v2(k,it))
            uk = dsqrt(1-wf1%v2(k,it))
            if(with_f) then 
                f_k = wf1%skk(k,it)
            else 
                f_k = 1.d0
            end if 
            truncated_k = wf1%truncated_index(k,it)
            do m1 = 1,nf
                ! nr1= BS%HO_sph%nljm(i0f+m1,1) ! n_r
                nl1= BS%HO_sph%nljm(i0f+m1,2) ! l
                nj1= BS%HO_sph%nljm(i0f+m1,3) ! j +1/2
                nm1= BS%HO_sph%nljm(i0f+m1,4) ! m_j + 1/2
                nj1_half = nj1 - 0.5 ! j
                nm1_half = nm1 - 0.5 ! m_j
                ! nm1_reversed = 1 - nm1 ! -m_j + 1/2
                ! nm1_reversed_half = nm1_reversed - 0.5 ! -m_j
                m1_reversed = m1 - 2*nm1 + 1 ! location of |-m1>
                Fu(truncated_k,m1) = wf1%fg(m1,k,it)*uk*dsqrt(f_k)
                Fu(truncated_dim+truncated_k,m1) = wf1%fg(m1_reversed,k,it)*uk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl1+nj1_half+nm1_half))
                pFu(truncated_k,m1) = wf1%fg(m1,k,it)*uk*dsqrt(f_k)
                pFu(truncated_dim+truncated_k,m1) = wf1%fg(m1_reversed,k,it)*uk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl1+nj1_half+nm1_half))
                Fv(truncated_k,m1) = wf1%fg(m1,k,it)*vk*dsqrt(f_k)
                Fv(truncated_dim+truncated_k,m1) = wf1%fg(m1_reversed,k,it)*vk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl1+nj1_half+nm1_half))
                pFv(truncated_k,m1) = wf1%fg(m1,k,it)*vk*dsqrt(f_k)
                pFv(truncated_dim+truncated_k,m1) = wf1%fg(m1_reversed,k,it)*vk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl1+nj1_half+nm1_half))
            end do
            do m2 = 1,ng
                ! nr2= BS%HO_sph%nljm(i0g+m2,1) ! n_r
                nl2= BS%HO_sph%nljm(i0g+m2,2) ! l
                nj2= BS%HO_sph%nljm(i0g+m2,3) ! j +1/2
                nm2= BS%HO_sph%nljm(i0g+m2,4) ! m_j + 1/2
                nj2_half = nj2 - 0.5 ! j
                nm2_half = nm2 - 0.5 ! m_j
                ! nm2_reversed = 1 - nm2 ! -m_j + 1/2
                ! nm2_reversed_half = nm2_reversed - 0.5 ! -m_j
                m2_reversed = m2 - 2*nm2 + 1 ! location of |-m2>
                Gu(truncated_k,m2) = wf1%fg(nf+m2,k,it)*uk*dsqrt(f_k)
                Gu(truncated_dim+truncated_k,m2) = -wf1%fg(nf+m2_reversed,k,it)*uk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl2+nj2_half+nm2_half))
                pGu(truncated_k,m2) = wf1%fg(nf+m2,k,it)*uk*dsqrt(f_k)
                pGu(truncated_dim+truncated_k,m2) = -wf1%fg(nf+m2_reversed,k,it)*uk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl2+nj2_half+nm2_half))
                Gv(truncated_k,m2) = wf1%fg(nf+m2,k,it)*vk*dsqrt(f_k)
                Gv(truncated_dim+truncated_k,m2) = -wf1%fg(nf+m2_reversed,k,it)*vk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl2+nj2_half+nm2_half))
                pGv(truncated_k,m2) = wf1%fg(nf+m2,k,it)*vk*dsqrt(f_k)
                pGv(truncated_dim+truncated_k,m2) = -wf1%fg(nf+m2_reversed,k,it)*vk*dsqrt(f_k) &
                                                    *(-1)**(Int(nl2+nj2_half+nm2_half))
            end do
        end do 
        ! kappa10    
        call cmatmul_ABC(RFv(1:nf,1:2*truncated_dim), B10(1:2*truncated_dim,1:2*truncated_dim), Gu(1:2*truncated_dim,1:ng), & 
                        tmp(1:nf,1:ng), nf, 2*truncated_dim, 2*truncated_dim, ng)
        mix%kappa10_mm(:,:,3,iphi,it) = tmp ! ifg = 3 is kappa10^+-
        call cmatmul_ABC(pRFv(1:nf,1:2*truncated_dim), pB10(1:2*truncated_dim,1:2*truncated_dim), pGu(1:2*truncated_dim,1:ng), & 
                        tmp(1:nf,1:ng), nf, 2*truncated_dim, 2*truncated_dim, ng)
        mix%pkappa10_mm(:,:,3,iphi,it) = tmp ! ifg = 3 is pkappa10^+-
        call cmatmul_ABC(RGv(1:ng,1:2*truncated_dim), B10(1:2*truncated_dim,1:2*truncated_dim), Fu(1:2*truncated_dim,1:nf), & 
                        tmp(1:ng,1:nf), ng, 2*truncated_dim, 2*truncated_dim, nf)
        mix%kappa10_mm(:,:,4,iphi,it) = tmp ! ifg = 4 is kappa10^-+
        call cmatmul_ABC(pRGv(1:ng,1:2*truncated_dim), pB10(1:2*truncated_dim,1:2*truncated_dim), pFu(1:2*truncated_dim,1:nf), & 
                        tmp(1:ng,1:nf), ng, 2*truncated_dim, 2*truncated_dim, nf)
        mix%pkappa10_mm(:,:,4,iphi,it) = tmp ! ifg = 4 is pkappa10^-+
        ! kappa01*
        call cmatmul_ABC(RFu(1:nf,1:2*truncated_dim), B01(1:2*truncated_dim,1:2*truncated_dim), Gv(1:2*truncated_dim,1:ng), &
                        tmp(1:nf,1:ng), nf, 2*truncated_dim, 2*truncated_dim, ng)
        do m1 = 1, ng
            do m2 = 1, nf
                mix%kappa01c_mm(m1,m2,4,iphi,it) = tmp(m2,m1) ! ! ifg = 4 is kappa01^-+*
            end do
        end do
        call cmatmul_ABC(pRFu(1:nf,1:2*truncated_dim), pB01(1:2*truncated_dim,1:2*truncated_dim), pGv(1:2*truncated_dim,1:ng), &
                        tmp(1:nf,1:ng), nf, 2*truncated_dim, 2*truncated_dim, ng)
        do m1 = 1, ng
            do m2 = 1, nf
                mix%pkappa01c_mm(m1,m2,4,iphi,it) = tmp(m2,m1) ! ! ifg = 4 is pkappa01^-+*
            end do
        end do
        call cmatmul_ABC(RGu(1:ng,1:2*truncated_dim), B01(1:2*truncated_dim,1:2*truncated_dim), Fv(1:2*truncated_dim,1:nf), &
                        tmp(1:ng,1:nf), ng, 2*truncated_dim, 2*truncated_dim, nf)
        do m1 = 1, nf
            do m2 = 1, ng
                mix%kappa01c_mm(m1,m2,3,iphi,it) = tmp(m2,m1) ! ! ifg = 3 is kappa01^+-*
            end do
        end do
        call cmatmul_ABC(pRGu(1:ng,1:2*truncated_dim), pB01(1:2*truncated_dim,1:2*truncated_dim), pFv(1:2*truncated_dim,1:nf), &
                        tmp(1:ng,1:nf), ng, 2*truncated_dim, 2*truncated_dim, nf)
        do m1 = 1, nf
            do m2 = 1, ng
                mix%pkappa01c_mm(m1,m2,3,iphi,it) = tmp(m2,m1) ! ! ifg = 3 is pkappa01^+-*
            end do
        end do
    end subroutine
end subroutine

subroutine calculate_mixed_density_current_tensor_in_coordinate_space(iphi,it)
    !---------------------------------------------------------------------------------------------------------------
    !  calculate mixed density, current and tensor in coordinate space
    !   \rho_S
    !   \rho_V
    !   \nabla^2\rho_S
    !   \nabla^2\rho_V
    !   j_V1, j_V2, j_V3
    !   \kappa
    !   \kappa*
    ! 
    !  Note
    !  Because Am, Dm, Y1, Y2, Y3 are multiplied by sqrt(x^2 sin(theta) wx wtheta wphi) and divided by b^(-3/2), so the
    !  the mixed density, current and tensor are multiplied by (x^2 sin(theta) wx wtheta wphi) and divided by b^(-3).
    !
    !  The true density,current and tensor originally depends on (r, \theta, \phi); when divided by b^(-3), they are 
    !  now expressed as a function of (x, \theta, \phi).
    !  The true density, current, and tensor can be obtained as  
    !      F(r,theta,phi) = b^(-3)*f(x,theta,phi)/(x^2 sin(theta) wx wtheta wphi), 
    !  where x = r/b
    !----------------------------------------------------------------------------------------------------------------

    integer,intent(in)::  iphi,it
    call mixed_density
    call mixed_current
    call mixed_tensor
    contains
    subroutine mixed_density
        !---------------------------------------------------------------------------------------------------------------
        !  calculate mixed density in coordinate space
        !
        !  \rho_S = \sum_{mm'\sigma}\rho_mm'++\Phi^\dagger_m'\Phi_m - \sum_{nn'\sigma}\rho__nn'--\Phi^\dagger_n'\Phi_n
        !  \rho_V = \sum_{mm'\sigma}\rho_mm'++\Phi^\dagger_m'\Phi_m + \sum_{nn'\sigma}\rho__nn'--\Phi^\dagger_n'\Phi_n
        !  \nabla^2\rho_S = \sum_{mm'\sigma}\rho_mm'++ \nabla^2(\Phi^\dagger_m'\Phi_m) 
        !                   - \sum_{nn'\sigma}\rho__nn'--\nabla^2(\Phi^\dagger_n'\Phi_n)
        !  \nabla^2\rho_V = \sum_{mm'\sigma}\rho_mm'++ \nabla^2(\Phi^\dagger_m'\Phi_m)
        !                   + \sum_{nn'\sigma}\rho__nn'--\nabla^2(\Phi^\dagger_n'\Phi_n)
        !
        !  where \Phi_m = \Phi_{m}(r, \theta, \phi, \sigma)
        !     m is the large component spherical harmonic basis index
        !     n is the small component spherical harmonic basis index
        !-----------------------------------------------------------------------------------------------------------------
        use Constants, only: ngr,ntheta,nphi 
        use Globals, only: BS
        use omp_lib
        integer :: i,ifg,basis_dim,m1,m2,ms
        complex(r64) :: Am,Bm,Dm,Em,Y1,Y2,Y3,Z1,Z2,Z3
        complex(r64), dimension(2) :: tmp_rho1, tmp_rho2,tmp_prho1,tmp_prho2
        !$omp parallel do default(none) & 
        !$omp shared(BS,mix,iphi,it) &
        !$omp private(i,ifg,basis_dim,m1,m2,ms,Am,Bm,Dm,Em,Y1,Y2,Y3,Z1,Z2,Z3,tmp_rho1,tmp_rho2,tmp_prho1,tmp_prho2)
        
        do i = 1, ngr*ntheta*nphi
            do ifg = 1,2 
                tmp_rho1(ifg) = (0.d0,0.d0)
                tmp_rho2(ifg) = (0.d0,0.d0)
                tmp_prho1(ifg) = (0.d0,0.d0)
                tmp_prho2(ifg) = (0.d0,0.d0)
                basis_dim = BS%HO_sph%idsp(1,ifg)
                do m1 = 1, basis_dim
                    do m2 = 1, basis_dim
                        do ms = 0,1
                            Am = BS%HO_sph%Am(m1,i,ms,ifg)
                            Bm = conjg(BS%HO_sph%Am(m2,i,ms,ifg))
                            Dm = BS%HO_sph%Dm(m1,i,ms,ifg)
                            Em = conjg(BS%HO_sph%Dm(m2,i,ms,ifg))
                            Y1 = BS%HO_sph%Y1(m1,i,ms,ifg)
                            Y2 = BS%HO_sph%Y2(m1,i,ms,ifg)
                            Y3 = BS%HO_sph%Y3(m1,i,ms,ifg)
                            Z1 = conjg(BS%HO_sph%Y1(m2,i,ms,ifg))
                            Z2 = conjg(BS%HO_sph%Y2(m2,i,ms,ifg))
                            Z3 = conjg(BS%HO_sph%Y3(m2,i,ms,ifg))
                            tmp_rho1(ifg) = tmp_rho1(ifg) + mix%rho_mm(m1,m2,ifg,iphi,it)*Am*Bm
                            tmp_rho2(ifg) = tmp_rho2(ifg) + mix%rho_mm(m1,m2,ifg,iphi,it)*(Dm*Bm+Am*Em+2*Y1*Z1+2*Y2*Z2+2*Y3*Z3)
                            tmp_prho1(ifg) = tmp_prho1(ifg) + mix%prho_mm(m1,m2,ifg,iphi,it)*Am*Bm
                            tmp_prho2(ifg) = tmp_prho2(ifg) + mix%prho_mm(m1,m2,ifg,iphi,it)*(Dm*Bm+Am*Em+2*Y1*Z1+2*Y2*Z2+2*Y3*Z3)
                        end do
                    end do
                end do
            end do
            mix%rho_S_it(i,iphi,it) =  tmp_rho1(1) - tmp_rho1(2)
            mix%rho_V_it(i,iphi,it) =  tmp_rho1(1) + tmp_rho1(2)
            mix%d2rho_S_it(i,iphi,it) = tmp_rho2(1) - tmp_rho2(2)
            mix%d2rho_V_it(i,iphi,it) = tmp_rho2(1) + tmp_rho2(2)
            !---
            mix%prho_S_it(i,iphi,it) =  tmp_prho1(1) - tmp_prho1(2)
            mix%prho_V_it(i,iphi,it) =  tmp_prho1(1) + tmp_prho1(2)
            mix%pd2rho_S_it(i,iphi,it) = tmp_prho2(1) - tmp_prho2(2)
            mix%pd2rho_V_it(i,iphi,it) = tmp_prho2(1) + tmp_prho2(2)
        end do 
        !$omp end parallel do
    end subroutine
    subroutine mixed_current
        !--------------------------------------------------------------------------------------------------------------------------------
        !  calculate mixed current in coordinate space
        !
        !  j_V1 = -i\sum_{mn'\sigma}\rho_mn'+-\Phi^\dagger_n' \sigma_1 \Phi_m + i\sum_{nm'\sigma}\rho_nm'-+\Phi^\dagger_m' \sigma_1 \Phi_n
        !  j_V2 = -i\sum_{mn'\sigma}\rho_mn'+-\Phi^\dagger_n' \sigma_2 \Phi_m + i\sum_{nm'\sigma}\rho_nm'-+\Phi^\dagger_m' \sigma_2 \Phi_n
        !  j_V3 = -i\sum_{mn'\sigma}\rho_mn'+-\Phi^\dagger_n' \sigma_3 \Phi_m + i\sum_{nm'\sigma}\rho_nm'-+\Phi^\dagger_m' \sigma_3 \Phi_n
        !  where \Phi_m = \Phi_{m}(r, \theta, \phi, \sigma)
        !     m is the large component spherical harmonic basis index
        !     n is the small component spherical harmonic basis index
        !---------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: ngr,ntheta,nphi
        use Globals, only: BS
        integer :: i,ifg,m1,m2,ms
        complex(r64) :: ci,Am,Bm,Cm
        complex(r64), dimension(2) :: tmp_j1,tmp_j2,tmp_j3,tmp_pj1,tmp_pj2,tmp_pj3
        ci = dcmplx(0,1) 
        !$omp parallel do default(none) &
        !$omp shared(ci,BS, mix, iphi, it) &
        !$omp private(i,ifg,m1,m2,ms,Am,Bm,Cm,tmp_j1,tmp_j2,tmp_j3,tmp_pj1,tmp_pj2,tmp_pj3)
        do i = 1, ngr*ntheta*nphi
            do ifg = 1,2
                tmp_j1(ifg) = (0.d0,0.d0)
                tmp_j2(ifg) = (0.d0,0.d0)
                tmp_j3(ifg) = (0.d0,0.d0)
                tmp_pj1(ifg) = (0.d0,0.d0)
                tmp_pj2(ifg) = (0.d0,0.d0)
                tmp_pj3(ifg) = (0.d0,0.d0)
                do m1 = 1, BS%HO_sph%idsp(1,ifg)
                    do m2 = 1, BS%HO_sph%idsp(1,3-ifg)
                        do ms = 0,1 
                            Am = BS%HO_sph%Am(m1,i,ms,ifg)
                            Bm = conjg(BS%HO_sph%Am(m2,i,1-ms,3-ifg))
                            Cm = conjg(BS%HO_sph%Am(m2,i,ms,3-ifg))
                            tmp_j1(ifg) = tmp_j1(ifg) + mix%rho_mm(m1,m2,2+ifg,iphi,it)*Am*Bm
                            tmp_j2(ifg) = tmp_j2(ifg) + mix%rho_mm(m1,m2,2+ifg,iphi,it)*Am*Bm*(-1)**ms
                            tmp_j3(ifg) = tmp_j3(ifg) + mix%rho_mm(m1,m2,2+ifg,iphi,it)*Am*Cm*(-1)**(ms-1)
                            tmp_pj1(ifg) = tmp_pj1(ifg) + mix%prho_mm(m1,m2,2+ifg,iphi,it)*Am*Bm
                            tmp_pj2(ifg) = tmp_pj2(ifg) + mix%prho_mm(m1,m2,2+ifg,iphi,it)*Am*Bm*(-1)**ms
                            tmp_pj3(ifg) = tmp_pj3(ifg) + mix%prho_mm(m1,m2,2+ifg,iphi,it)*Am*Cm*(-1)**(ms-1)
                        end do
                    end do
                end do 
            end do
            mix%j_V1_it(i,iphi,it) =  (- tmp_j1(1) + tmp_j1(2))*ci
            mix%j_V2_it(i,iphi,it) =   - tmp_j2(1) + tmp_j2(2)
            mix%j_V3_it(i,iphi,it) =  (- tmp_j3(1) + tmp_j3(2))*ci
            mix%pj_V1_it(i,iphi,it) =  (- tmp_pj1(1) + tmp_pj1(2))*ci
            mix%pj_V2_it(i,iphi,it) =   - tmp_pj2(1) + tmp_pj2(2)
            mix%pj_V3_it(i,iphi,it) =  (- tmp_pj3(1) + tmp_pj3(2))*ci
        end do 
        !$omp end parallel do
    end subroutine
    subroutine mixed_tensor
        !--------------------------------------------------------------------------------------------------------------------
        !  calculate mixed tensor matrix elements in coordinate space
        !
        !  \kappa  = \sum_{mm'\sigma} \kappa_mm'++\Phi*_m'\Phi_m + \sum_{nn'\sigma} \kappa_nn'--\Phi*_n'\Phi_n
        !  \kappa* = \sum_{mm'\sigma}\kappa*_mm'++\Phi*_m'\Phi_m + \sum_{nn'\sigma}\kappa*_nn'--\Phi*_n'\Phi_n
        !  where \Phi_m = \Phi_{m}(r, \theta, \phi, \sigma)
        !     m is the large component spherical harmonic basis index
        !     n is the small component spherical harmonic basis index
        !----------------------------------------------------------------------------------------------------------------------
        use Constants, only: ngr,ntheta,nphi
        use Globals, only: BS
        integer :: i,ifg,basis_dim,m1,m2,ms
        complex(r64) :: Am,Bm
        complex(r64), dimension(2) ::tmp_kappa,tmp_kappac,tmp_pkappa,tmp_pkappac
        !$omp parallel do default(none) &
        !$omp shared(BS, mix, iphi, it) &
        !$omp private(i,ifg,basis_dim,m1,m2,ms,Am,Bm,tmp_kappa,tmp_kappac,tmp_pkappa,tmp_pkappac)
        do i = 1, ngr*ntheta*nphi
            do ifg = 1,2 
                tmp_kappa(ifg) = (0.d0,0.d0)
                tmp_kappac(ifg) = (0.d0,0.d0)
                tmp_pkappa(ifg) = (0.d0,0.d0)
                tmp_pkappac(ifg) = (0.d0,0.d0)
                basis_dim = BS%HO_sph%idsp(1,ifg)
                do m1 = 1, basis_dim
                    do m2 = 1, basis_dim
                        do ms = 0,1
                            Am = BS%HO_sph%Am(m1,i,ms,ifg)
                            Bm = conjg(BS%HO_sph%Am(m2,i,ms,ifg))
                            tmp_kappa(ifg) = tmp_kappa(ifg) + mix%kappa_mm(m1,m2,ifg,iphi,it)*Am*Bm
                            tmp_kappac(ifg) = tmp_kappac(ifg) + mix%kappac_mm(m1,m2,ifg,iphi,it)*Am*Bm
                            tmp_pkappa(ifg) = tmp_pkappa(ifg) + mix%pkappa_mm(m1,m2,ifg,iphi,it)*Am*Bm
                            tmp_pkappac(ifg) = tmp_pkappac(ifg) + mix%pkappac_mm(m1,m2,ifg,iphi,it)*Am*Bm
                        end do
                    end do
                end do
            end do
            mix%kappa_it(i,iphi,it) =  tmp_kappa(1) + tmp_kappa(2)
            mix%kappac_it(i,iphi,it) = tmp_kappac(1) + tmp_kappac(2)
            mix%pkappa_it(i,iphi,it) = tmp_pkappa(1) + tmp_pkappa(2)
            mix%pkappac_it(i,iphi,it) = tmp_pkappac(1) + tmp_pkappac(2)
        end do
        !$omp end parallel do
    end subroutine
end subroutine

subroutine calculate_norm_overlap(iphi,phi,it)
    !--------------------------------------------------------------------------------------
    !  calculate norm overlap
    ! 
    !  1) n(q_1, q_2, Theta) =  n(q_1, q_2, alpha, beta, gamma, phi)
    !       = \sqrt{det(D) det(R)}
    !  2) Pfaffian: Robledo, PRC79,021302(R)(2009)
    !  3) Pfaffian: Bertsch and Robledo, PRL108, 042505 (2012)
    !--------------------------------------------------------------------------------------
    use Globals, only: Proj_option
    integer,intent(in)::  iphi,it
    real(r64) :: phi
    complex(r64) :: norm_overlap, pnorm_overlap

    if(Proj_option%ihf==1) then
        ! way 1 
        mix%norm(iphi,it) = cdsqrt(D_Deter(it) * RT_Deter(it))
        mix%pnorm(iphi,it) = cdsqrt(pD_Deter(it) * pRT_Deter(it))
        write(222,*) 'it', it,'iphi', iphi,  'norm', mix%norm(iphi,it), 'pnorm',mix%pnorm(iphi,it)
    else if(Proj_option%ihf==2) then
        ! way 2 
        call Pfaffian_Robledo()
    else if(Proj_option%ihf==3) then
        ! way 3
        call Pfaffian_Bertsch(phi,it,norm_overlap,pnorm_overlap)
        mix%norm(iphi,it) = norm_overlap
        mix%pnorm(iphi,it) = pnorm_overlap
        ! write(222,*) 'it', it,'iphi', iphi,  'norm', mix%norm(iphi,it), 'pnorm',mix%pnorm(iphi,it)
        ! Pfaffian_Bertsch1 is same as Pfaffian_Bertsch
        ! call Pfaffian_Bertsch1(phi,it,norm_overlap,pnorm_overlap)
        ! mix%norm(iphi,it) = norm_overlap
        ! mix%pnorm(iphi,it) = pnorm_overlap
        ! write(222,*) 'it', it,'iphi', iphi,  'norm', mix%norm(iphi,it), 'pnorm',mix%pnorm(iphi,it)

    end if 
    contains

    subroutine Pfaffian_Robledo()
        stop 'Not implement yet!'
    end subroutine
    subroutine Pfaffian_Bertsch(phi,it,norm_overlap,pnorm_overlap)
        !---------------------------------------------------------------
        !   Bertsch and Robledo, PRL108, 042505 (2012).
        !
        !    1. construct the skew-symmetric matrix M (4n * 4n)
        !          (  V(a)^T U^(a) ,  V(a)^T R V^{(b)*}  )
        !      M = (   ------------   ------------------ )
        !          (               ,  U^{(b)+} V{(b)^+}  )
        !
        !          (      M^(a) ,           M^(c)        )
        !        = (   ------------   ------------------ )
        !          (               ,        M^(b)        )
        !
        !    2. calculate the Pfaffian of M
        !
        !    3. calculate the overlap <qa|qb> = (-1)^n * Pf(M)/ (N^(a)*N^(b))
        !       where normalization factor N^(a) = Prod_i v^(a)_i.
        !
        !      Pf(M): Pfaffian of a skew-symmetric matrix M
        !-----------------------------------------------------------------
        use Globals, only: wf1, wf2
        use MathMethods,only: cmatmul_ABC, ZPfaffianF
        real(r64), intent(in):: phi
        integer, intent(in) :: it
        complex(r64), intent(out):: norm_overlap,pnorm_overlap 
        complex(r64),dimension(:,:),allocatable :: M, pM, V1T_matrix, V2_matrix,R_matrix,pR_matrix, tmp_D, tmp_pD
        integer, dimension(:,:),allocatable:: IPIV
        integer :: dim1,dim2, k,l,truncated_k,truncated_l,truncated_k2,truncated_l2
        real(r64) :: vk,uk,ul
        complex(r64) :: vl,N1,N2,cphi,eiphi,ei2phi,Pf_M,Pf_pM
        cphi   = cmplx(0.d0,phi)
        eiphi = cdexp(cphi)
        ei2phi = cdexp(2*cphi)
        dim1 = wf1%truncated_dim(it)
        dim2 = wf2%truncated_dim(it)
        allocate(M((dim1+dim2)*2,(dim1+dim2)*2),pM((dim1+dim2)*2,(dim1+dim2)*2),source=(0.d0,0.d0))
        ! V(a)^T U^(a)
        N1= (1.d0,0.d0)
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
            vk = dsqrt(wf1%v2(k,it))
            uk = dsqrt(1-wf1%v2(k,it))
            truncated_k = wf1%truncated_index(k,it)
            truncated_k2 = truncated_k*2
            M(truncated_k2-1,truncated_k2) = -vk*uk
            M(truncated_k2,truncated_k2-1) = vk*uk
            pM(truncated_k2-1,truncated_k2) = -vk*uk
            pM(truncated_k2,truncated_k2-1) = vk*uk
            N1 = N1*vk
        end do
    
        ! U^{(b)+} V{(b)^*}
        N2= (1.d0,0.d0)
        do l = 1, wf2%nk(it)
            if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
            vl = dsqrt(wf2%v2(l,it)) *ei2phi
            ul = dsqrt(1-wf2%v2(l,it))
            truncated_l = wf2%truncated_index(l,it)
            truncated_l2 = truncated_l*2 + dim1*2
            M(truncated_l2-1,truncated_l2) = vl*ul
            M(truncated_l2,truncated_l2-1) = -vl*ul
            pM(truncated_l2-1,truncated_l2) = vl*ul
            pM(truncated_l2,truncated_l2-1) = -vl*ul
            N2 = N2*vl
        end do
    
        ! V(a)^T R V^{(b)*}
        allocate(V1T_matrix(dim1*2, dim1*2),V2_matrix(dim2*2,dim2*2),R_matrix(dim1*2, dim2*2),pR_matrix(dim1*2, dim2*2),source=(0.d0,0.d0))
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
            vk = dsqrt(wf1%v2(k,it))
            ! uk = dsqrt(1-wf1%v2(k,it))
            truncated_k = wf1%truncated_index(k,it)
            truncated_k2 = truncated_k*2
            V1T_matrix(truncated_k2-1,truncated_k2) = -vk
            V1T_matrix(truncated_k2,truncated_k2-1) = vk
        end do 
        do l = 1, wf2%nk(it)
            if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
            vl = dsqrt(wf2%v2(l,it)) * ei2phi
            ! ul = dsqrt(1-wf2%v2(l,it))
            truncated_l = wf2%truncated_index(l,it)
            truncated_l2 = truncated_l*2
            V2_matrix(truncated_l2-1,truncated_l2) = vl
            V2_matrix(truncated_l2,truncated_l2-1) = -vl
        end do
        do truncated_k = 1, dim1
            truncated_k2 = truncated_k*2
            do truncated_l = 1, dim2
                truncated_l2 = truncated_l*2
                R_matrix(truncated_k2-1,truncated_l2-1) = R(truncated_k,truncated_l,it)
                R_matrix(truncated_k2-1,truncated_l2)   = R(truncated_k,dim2+truncated_l,it)
                R_matrix(truncated_k2,truncated_l2-1)   = R(dim1+truncated_k,truncated_l,it)
                R_matrix(truncated_k2,truncated_l2)     = R(dim1+truncated_k,dim2+truncated_l,it)
                pR_matrix(truncated_k2-1,truncated_l2-1) = pR(truncated_k,truncated_l,it)
                pR_matrix(truncated_k2-1,truncated_l2)   = pR(truncated_k,dim2+truncated_l,it)
                pR_matrix(truncated_k2,truncated_l2-1)   = pR(dim1+truncated_k,truncated_l,it)
                pR_matrix(truncated_k2,truncated_l2)     = pR(dim1+truncated_k,dim2+truncated_l,it)
            end do 
        end do
        allocate(tmp_D(dim1*2,dim2*2),tmp_pD(dim1*2,dim2*2) )
        call cmatmul_ABC(V1T_matrix,R_matrix,V2_matrix,tmp_D ,dim1*2,dim1*2,dim2*2,dim2*2)
        call cmatmul_ABC(V1T_matrix,pR_matrix,V2_matrix,tmp_pD,dim1*2,dim1*2,dim2*2,dim2*2)
        ! 
        do k = 1, dim1*2
            do l = 1, dim2*2
                M(k,dim1*2+l) = tmp_D(k,l)
                M(dim1*2+l,k) = - tmp_D(k,l)
                pM(k,dim1*2+l) = tmp_pD(k,l)
                pM(dim1*2+l,k) = - tmp_pD(k,l)
            end do
        end do 
    
        allocate(IPIV((dim1+dim2)*2,(dim1+dim2)*2)) 
        call ZPfaffianF(M,(dim1+dim2)*2,(dim1+dim2)*2,IPIV,Pf_M)
        call ZPfaffianF(pM,(dim1+dim2)*2,(dim1+dim2)*2,IPIV,Pf_pM)
        norm_overlap = (-1)**(dim1)*Pf_M/(N1*N2)
        pnorm_overlap = (-1)**(dim1)*Pf_pM/(N1*N2)
    end subroutine
    subroutine Pfaffian_Bertsch1(phi,it,norm_overlap,pnorm_overlap)
        !---------------------------------------------------------------
        !   Bertsch and Robledo, PRL108, 042505 (2012).
        !
        !    1. construct the skew-symmetric matrix M (4n * 4n)
        !          (  V(a)^T U^(a) ,  V(a)^T R V^{(b)*}  )
        !      M = (   ------------   ------------------ )
        !          (               ,  U^{(b)+} V{(b)^+}  )
        !
        !          (      M^(a) ,           M^(c)        )
        !        = (   ------------   ------------------ )
        !          (               ,        M^(b)        )
        !
        !    2. calculate the Pfaffian of M
        !
        !    3. calculate the overlap <qa|qb> = (-1)^n * Pf(M)/ (N^(a)*N^(b))
        !       where normalization factor N^(a) = Prod_i v^(a)_i.
        !
        !      Pf(M): Pfaffian of a skew-symmetric matrix M
        !-----------------------------------------------------------------
        use Globals, only: wf1, wf2
        use MathMethods,only: cmatmul_ABC, ZPfaffianF
        real(r64), intent(in):: phi
        integer, intent(in) :: it
        complex(r64), intent(out):: norm_overlap,pnorm_overlap 
        complex(r64),dimension(:,:),allocatable :: M, pM, V1T_matrix, V2_matrix,R_matrix,pR_matrix, tmp_D, tmp_pD
        integer, dimension(:,:),allocatable:: IPIV
        integer :: dim1,dim2, k,l,truncated_k,truncated_l,truncated_k2,truncated_l2
        real(r64) :: vk,uk,ul
        complex(r64) :: vl,N1,N2,cphi,eiphi,ei2phi,Pf_M,Pf_pM
        cphi   = cmplx(0.d0,phi)
        eiphi = cdexp(cphi)
        ei2phi = cdexp(2*cphi)
        dim1 = wf1%truncated_dim(it)
        dim2 = wf2%truncated_dim(it)
        allocate(M((dim1+dim2)*2,(dim1+dim2)*2),pM((dim1+dim2)*2,(dim1+dim2)*2),source=(0.d0,0.d0))
        ! V(a)^T U^(a)
        N1= (1.d0,0.d0)
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
            vk = dsqrt(wf1%v2(k,it))
            uk = dsqrt(1-wf1%v2(k,it))
            truncated_k = wf1%truncated_index(k,it)
            truncated_k2 = truncated_k*2
            M(truncated_k2-1,truncated_k2) = -vk*uk
            M(truncated_k2,truncated_k2-1) = vk*uk
            pM(truncated_k2-1,truncated_k2) = -vk*uk
            pM(truncated_k2,truncated_k2-1) = vk*uk
            N1 = N1*vk
        end do
    
        ! U^{(b)+} V{(b)^*}
        N2= (1.d0,0.d0)
        do l = 1, wf2%nk(it)
            if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
            vl = dsqrt(wf2%v2(l,it)) !*ei2phi
            ul = dsqrt(1-wf2%v2(l,it))
            truncated_l = wf2%truncated_index(l,it)
            truncated_l2 = truncated_l*2 + dim1*2
            M(truncated_l2-1,truncated_l2) = vl*ul
            M(truncated_l2,truncated_l2-1) = -vl*ul
            pM(truncated_l2-1,truncated_l2) = vl*ul
            pM(truncated_l2,truncated_l2-1) = -vl*ul
            N2 = N2*vl
        end do
    
        ! V(a)^T R V^{(b)*}
        allocate(V1T_matrix(dim1*2, dim1*2),V2_matrix(dim2*2,dim2*2),R_matrix(dim1*2, dim2*2),pR_matrix(dim1*2, dim2*2),source=(0.d0,0.d0))
        do k = 1, wf1%nk(it)
            if(wf1%v2d(k,it) <= wf1%eps_occupation(it)) cycle
            vk = dsqrt(wf1%v2(k,it))
            ! uk = dsqrt(1-wf1%v2(k,it))
            truncated_k = wf1%truncated_index(k,it)
            truncated_k2 = truncated_k*2
            V1T_matrix(truncated_k2-1,truncated_k2) = -vk
            V1T_matrix(truncated_k2,truncated_k2-1) = vk
        end do 
        do l = 1, wf2%nk(it)
            if(wf2%v2d(l,it) <= wf2%eps_occupation(it)) cycle
            vl = dsqrt(wf2%v2(l,it)) * eiphi
            ! ul = dsqrt(1-wf2%v2(l,it))
            truncated_l = wf2%truncated_index(l,it)
            truncated_l2 = truncated_l*2
            V2_matrix(truncated_l2-1,truncated_l2) = vl
            V2_matrix(truncated_l2,truncated_l2-1) = -vl
        end do
        do truncated_k = 1, dim1
            truncated_k2 = truncated_k*2
            do truncated_l = 1, dim2
                truncated_l2 = truncated_l*2
                R_matrix(truncated_k2-1,truncated_l2-1) = R(truncated_k,truncated_l,it)
                R_matrix(truncated_k2-1,truncated_l2)   = R(truncated_k,dim2+truncated_l,it)
                R_matrix(truncated_k2,truncated_l2-1)   = R(dim1+truncated_k,truncated_l,it)
                R_matrix(truncated_k2,truncated_l2)     = R(dim1+truncated_k,dim2+truncated_l,it)
                pR_matrix(truncated_k2-1,truncated_l2-1) = pR(truncated_k,truncated_l,it)
                pR_matrix(truncated_k2-1,truncated_l2)   = pR(truncated_k,dim2+truncated_l,it)
                pR_matrix(truncated_k2,truncated_l2-1)   = pR(dim1+truncated_k,truncated_l,it)
                pR_matrix(truncated_k2,truncated_l2)     = pR(dim1+truncated_k,dim2+truncated_l,it)
            end do 
        end do
        allocate(tmp_D(dim1*2,dim2*2),tmp_pD(dim1*2,dim2*2) )
        call cmatmul_ABC(V1T_matrix,R_matrix,V2_matrix,tmp_D ,dim1*2,dim1*2,dim2*2,dim2*2)
        call cmatmul_ABC(V1T_matrix,pR_matrix,V2_matrix,tmp_pD,dim1*2,dim1*2,dim2*2,dim2*2)
        ! 
        do k = 1, dim1*2
            do l = 1, dim2*2
                M(k,dim1*2+l) = tmp_D(k,l)
                M(dim1*2+l,k) = - tmp_D(k,l)
                pM(k,dim1*2+l) = tmp_pD(k,l)
                pM(dim1*2+l,k) = - tmp_pD(k,l)
            end do
        end do 
    
        allocate(IPIV((dim1+dim2)*2,(dim1+dim2)*2)) 
        call ZPfaffianF(M,(dim1+dim2)*2,(dim1+dim2)*2,IPIV,Pf_M)
        call ZPfaffianF(pM,(dim1+dim2)*2,(dim1+dim2)*2,IPIV,Pf_pM)
        norm_overlap = (-1)**(dim1)*Pf_M/(N1*N2)
        pnorm_overlap = (-1)**(dim1)*Pf_pM/(N1*N2)
    end subroutine
end subroutine

End Module Mixed