Module TD
    !------------------------------------------------------------------------------
    ! This module calculates transition density matrix element
    ! TODO: 
    !       1) implement the symmetry of rho_mm
    !-------------------------------------------------------------------------------
    use Constants, only: r64
    use Globals, only: TDs
    implicit none
    contains
    subroutine calculate_reduced_transition_density_matrix_element(q1,q2)
        !-----------------------------------------------------------------
        !
        !   calculate and store the reduced transition denstiy
        ! 
        ! 
        !  1B:
        !  TDs%reduced_TDME1B   =  <Jf Kf Parity_f qf|| M_lambda || Ji Ki Parity_i qi>
        !  TDs%reduced_TDME1B_c =  <Jf Kf Parity_f qi|| M_lambda || Ji Ki Parity_i qf>
        !      where M_lambda = [c_a^\dagger \tilde{c}_b]_\lambda
        !
        !-----------------------------------------------------------------
        use Globals, only: gcm_space,Proj_option
        integer, intent(in) :: q1,q2
        integer :: J,Ji,Jf,lambda,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,ifg,a,b,Parity_i,Parity_f,iPi,iPf
        complex(r64) :: reduced_TDME1B(2),reduced_TDME1B_c(2)
        logical :: q2_q1_Symmetry
        write(*,'(5x,A)') 'calculate_reduced_transition_density_matrix_element ...'
        call set_nlj_mapping
        if(q1/=q2 .and. Proj_option%Kernel_Symmetry==1) then
            q2_q1_Symmetry = .True.
        else
            q2_q1_Symmetry = .False.
        end if 
        ! reduced 1B transition density
        if(Proj_option%AMPtype==0 .or. Proj_option%AMPtype==1) then
            if(.not. allocated(TDs%reduced_TDME1B)) allocate(TDs%reduced_TDME1B(gcm_space%Jmin:gcm_space%Jmax+TDs%lambda_max, 0:0, 2,&        ! Jf,Kf,Pf
                                                                            0:TDs%lambda_max, &                                               ! lambda    
                                                                            gcm_space%Jmin:gcm_space%Jmax, 0:0, 2,&                           ! Ji,Ki,Pi
                                                                            2, TDs%nlj_length(2), TDs%nlj_length(2), 2 ),source=(0.d0,0.d0))  ! ifg,a,b,it  
            if(.not. allocated(TDs%reduced_TDME1B_c)) allocate(TDs%reduced_TDME1B_c(gcm_space%Jmin:gcm_space%Jmax+TDs%lambda_max, 0:0, 2,&        ! Jf,Kf,Pf
                                                                                0:TDs%lambda_max,&                                                ! lambda
                                                                                gcm_space%Jmin:gcm_space%Jmax, 0:0, 2,&                           ! Ji,Ki,Pi
                                                                                2, TDs%nlj_length(2), TDs%nlj_length(2), 2 ),source=(0.d0,0.d0))  ! ifg,a,b,it
        else
            if(.not. allocated(TDs%reduced_TDME1B)) allocate(TDs%reduced_TDME1B(gcm_space%Jmin:gcm_space%Jmax+TDs%lambda_max, -gcm_space%Jmax-TDs%lambda_max:gcm_space%Jmax+TDs%lambda_max, 2, &
                                                                            0:TDs%lambda_max, &
                                                                            gcm_space%Jmin:gcm_space%Jmax, -gcm_space%Jmax:gcm_space%Jmax, 2, &
                                                                            2, TDs%nlj_length(2), TDs%nlj_length(2), 2 ),source=(0.d0,0.d0))
            if(.not. allocated(TDs%reduced_TDME1B_c)) allocate(TDs%reduced_TDME1B_c(gcm_space%Jmin:gcm_space%Jmax+TDs%lambda_max, -gcm_space%Jmax-TDs%lambda_max:gcm_space%Jmax+TDs%lambda_max, 2, &
                                                                            0:TDs%lambda_max, &
                                                                            gcm_space%Jmin:gcm_space%Jmax, -gcm_space%Jmax:gcm_space%Jmax, 2, &
                                                                            2, TDs%nlj_length(2), TDs%nlj_length(2), 2 ),source=(0.d0,0.d0))
        end if 
        do Ji = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            do Jf = Ji, Ji+TDs%lambda_max
                if(Proj_option%AMPtype==0 .or. Proj_option%AMPtype==1) then
                    Ki_start = 0
                    Ki_end = 0
                    Kf_start = 0
                    Kf_end = 0
                else
                    Ki_start = -Ji
                    Ki_end = Ji
                    Kf_start = -Jf
                    Kf_end = Jf
                end if 
                do lambda = abs(Jf-Ji), min(TDs%lambda_max,Jf+Ji)
                    do Kf = Kf_start, Kf_end
                        do Ki = Ki_start, Ki_end 
                            do ifg = 1, 2
                                !$OMP PARALLEL DEFAULT(shared) PRIVATE(a,b,Parity_i,iPi,iPf,Parity_f,reduced_TDME1B,reduced_TDME1B_c)
                                !$OMP DO COLLAPSE(2) SCHEDULE(static)
                                do a = 1, TDs%nlj_length(ifg)
                                    do b = 1, TDs%nlj_length(ifg)
                                        Parity_i = (-1)**Ji ! In the axially symmetric case, the kernel is non-zero only when Parity_i * (-1)^J_i = 1
                                        iPi = (3-Parity_i)/2 ! +1: 1, -1: 2
                                        do iPf = 1,2
                                            Parity_f = (-1)**(iPf+1) ! 1: +1 , 2: -1
                                            call calculate_reduced_one_body_transition_density_matrix_element(Jf,Kf,Parity_f,lambda,Ji,Ki,Parity_i,ifg,a,b,reduced_TDME1B,reduced_TDME1B_c,q2_q1_Symmetry)
                                            ! q1-q2
                                            ! neutron
                                            TDs%reduced_TDME1B(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,1) =  reduced_TDME1B(1)
                                            ! proton
                                            TDs%reduced_TDME1B(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,2) =  reduced_TDME1B(2)
                                            if(q2_q1_Symmetry) then 
                                                ! q2-q1
                                                ! neutron
                                                TDs%reduced_TDME1B_c(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,1) =  reduced_TDME1B_c(1)
                                                ! proton
                                                TDs%reduced_TDME1B_c(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,2) =  reduced_TDME1B_c(2)
                                            end if
                                        end do 
                                    end do 
                                end do 
                                !$OMP END PARALLEL
                            end do
                        end do 
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine set_nlj_mapping
        !------------------------------------------------------------------------------------------------------------------
        !
        !  Extract the unique (n_r,l,j) from the spherical harmonic oscillator basis (n_r l j m_j).
        !  ----------
        !  TDs%nlj_index(index,ifg) : Stores the starting index in BS%HO_sph%nljm for the index-th unique (n_r,l,j).
        !  TDs%nlj_length(ifg)      : Total number of unique (n_r,l,j) (ifg = 1 for large component, 2 for small component).
        !------------------------------------------------------------------------------------------------------------------
        use Globals, only: BS,TDs
        integer :: ifg,ndsp,i0sp,index,m,nr,nl,nj,nm,a,b,a_index,nra,nla,nja,b_index,nrb,nlb,njb
        do ifg = 1,2
            ndsp = BS%HO_sph%idsp(1,ifg)
            i0sp = BS%HO_sph%iasp(1,ifg)
            index = 1
            TDs%nlj_index(index,ifg) = 1
            do m = 2,ndsp
                nr= BS%HO_sph%nljm(i0sp+m,1) ! n_r
                nl= BS%HO_sph%nljm(i0sp+m,2) ! l
                nj= BS%HO_sph%nljm(i0sp+m,3) ! j +1/2
                nm= BS%HO_sph%nljm(i0sp+m,4) ! m_j + 1/2
                if(nr/=BS%HO_sph%nljm(i0sp+m-1,1) .or. nl/=BS%HO_sph%nljm(i0sp+m-1,2) .or. nj/=BS%HO_sph%nljm(i0sp+m-1,3)) then
                    index = index + 1
                    TDs%nlj_index(index,ifg) = m
                end if
            end do
            TDs%nlj_length(ifg) = index
        end do
        ! test nlj 
        do ifg = 1, 2
            i0sp = BS%HO_sph%iasp(1,ifg)
            do a = 1, TDs%nlj_length(ifg)
                do b = a+1, TDs%nlj_length(ifg)
                    a_index = TDs%nlj_index(a,ifg)
                    nra= BS%HO_sph%nljm(i0sp+a_index,1) ! n_r
                    nla= BS%HO_sph%nljm(i0sp+a_index,2) ! l
                    nja= BS%HO_sph%nljm(i0sp+a_index,3) ! j +1/2
                    b_index = TDs%nlj_index(b,ifg)
                    nrb= BS%HO_sph%nljm(i0sp+b_index,1) ! n_r
                    nlb= BS%HO_sph%nljm(i0sp+b_index,2) ! l
                    njb= BS%HO_sph%nljm(i0sp+b_index,3) ! j +1/2
                    if(nra==nrb .and. nla==nlb .and. nja==njb) then
                        stop 'Worong nlj!'
                    end if 
                end do 
            end do
        end do 
    end subroutine
    
    subroutine calculate_reduced_one_body_transition_density_matrix_element(Jf,Kf,Parity_f,lambda,Ji,Ki,Parity_i,ifg,a,b,reduced_TDME1B,reduced_TDME1B_c,q2_q1_Symmetry)
        !-----------------------------------------------------------------------------------------------------------------------------
        !
        !   calculate reduced 1 body transition density
        !           <Jf Kf Parity_f qf|| M_lambda || Ji Ki Parity_i qi>
        ! 
        !   where M_lambda = [c_a^\dagger \tilde{c}_b]_\lambda
        !
        !  Relation between Clebsch–Gordan coefficients and Wigner 3j symbols:
        !  C^{lambda mu}_{ja ma jb -mb} = -(-1)**(lambda+mu) * dsqrt(2*lambda+1) * wigner3j(nja,lambda,njb,nma,-mu,1-nmb,IS=1)
        !  here j_i,J_f are hale integers.
        !
        ! q_1 , q_2 exchange :
        !           <Jf Kf Parity_f qi|| M_lambda || Ji Ki Parity_i qf>
        ! -----------------------------------------------
        !  Note: 
        !  1)Using axial symmetry, we do not calculate the mixed densities ( Proj_densities%rho_mm  or Proj_densities%prho_mm) for beta in [pi/2, pi] .
        !   However, it can be shown that the contribution from [pi/2, pi ] is the same as that from ([0, pi/2]),
        !   differing only by a factor of Parity_i*(-1)^{J_i}. 
        !   Here we take Parity_i = (-1)^{J_i}), so we can simply multiply by 2 without integrating over [pi/2, pi ].
        !------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: BS,projection_mesh,Proj_option,nucleus_attributes,mix,Proj_densities
        use Basis, only: djmk
        use EM, only: wigner3j
        integer, intent(in) :: Jf,Ji,lambda,Kf,Ki,ifg,a,b,Parity_f,Parity_i
        complex(r64),intent(out) :: reduced_TDME1B(2),reduced_TDME1B_c(2)
        integer :: i0sp,a_index,nra,nla,nja,b_index,nrb,nlb,njb,mu,K1,ialpha,ibeta,igamma,L_n,L_p,phi_n_index,phi_p_index,m1,m2,nma,nmb,it,nu,fac_Parity
        real(r64) :: alpha,beta,gamma,w,phi_n,phi_p
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac_AMP,emiNphi,emiZphi,fac_PNP,fac,norm,pnorm,cfac1,cfac_AMP
        logical :: q2_q1_Symmetry
        i0sp = BS%HO_sph%iasp(1,ifg)
        a_index = TDs%nlj_index(a,ifg)
        nra= BS%HO_sph%nljm(i0sp+a_index,1) ! n_r
        nla= BS%HO_sph%nljm(i0sp+a_index,2) ! l
        nja= BS%HO_sph%nljm(i0sp+a_index,3) ! j +1/2
        b_index = TDs%nlj_index(b,ifg)
        nrb= BS%HO_sph%nljm(i0sp+b_index,1) ! n_r
        nlb= BS%HO_sph%nljm(i0sp+b_index,2) ! l
        njb= BS%HO_sph%nljm(i0sp+b_index,3) ! j +1/2
        reduced_TDME1B = (0.d0, 0.d0) 
        reduced_TDME1B_c = (0.d0, 0.d0) 

        fac_Parity = (1 + Parity_i*Parity_f*(-1)**(nla+nlb))/2 ! 0 or 1
        if (fac_Parity==0) return
        if (fac_Parity/=1) stop 'wrong Parity_i or Parity_f !'

        do ialpha = 1, projection_mesh%nalpha
            do ibeta = 1, projection_mesh%nbeta
                ! If we have implemented the symmetry of rho_mm, then these lines are not needed.
                if(Proj_option%Euler_Symmetry == 2) then
                    stop '[calculate_reduced_one_body_transition_density]: Euler_Symmetry=2, Not yet implemented! You should set the Symmetry of Euler angles as 0!'
                end if
                if(Proj_option%Euler_Symmetry==1 .and. ibeta>(projection_mesh%nbeta+1)/2) then
                    cycle
                    ! Using the tensor symmetry, it can be proven that when Parity_i = (-1)**Ji
                    !  the contribution from (pi/2, pi]) is the same as that from (0, pi/2).
                end if 
                do igamma = 1, projection_mesh%ngamma
                    do mu = -lambda,lambda
                        K1 = Kf - mu
                        alpha = projection_mesh%alpha(ialpha)
                        calpha = DCMPLX(0.d0,alpha)
                        beta = projection_mesh%beta(ibeta)
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        if(Proj_option%AMPtype==0) then
                            fac_AMP = 1
                            cfac_AMP = 1
                        else if (Proj_option%AMPtype==1) then
                            ! In fact, when checking the AMP type (equal to 1) , 
                            ! the values of (alpha) and (gamma) have already been set to zero.
                            ! The following four assignments are unnecessary.
                            alpha = 0.d0
                            calpha = DCMPLX(0.d0,alpha)
                            gamma = 0.d0
                            cgamma = DCMPLX(0.d0,gamma)

                            w = projection_mesh%wbeta(ibeta)
                            fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,K1,Ki,dcos(beta),0)
                            fac_AMP = fac1*w
                            ! factor of qi-qf exchange
                            cfac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,Ki,K1,dcos(beta),0)
                            cfac_AMP = cfac1*w
                        else
                            cpi = DCMPLX(0.d0,pi) ! i*pi
                            w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                            fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,K1,Ki,dcos(beta),0)*CDEXP(-K1*calpha-Ki*cgamma)
                            fac2 = 1.0d0 + (-1)**mu*CDEXP(-K1*cpi) + CDEXP(-Ki*cpi) + (-1)**mu*CDEXP(-K1*cpi-Ki*cpi)! D2 symmetry is required, with alpha, gamma in [0, pi].
                            fac_AMP = fac1*fac2*w
                            ! factor of qi-qf exchange
                            cfac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,Ki,K1,dcos(beta),0)*CDEXP(Ki*calpha+K1*cgamma)
                            cfac_AMP = cfac1*fac2*w
                        end if
                        L_n = projection_mesh%nphi(1)
                        L_p = projection_mesh%nphi(2)
                        do phi_n_index = 1, L_n
                            phi_n =  phi_n_index*projection_mesh%dphi(1)
                            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
                            do phi_p_index = 1, L_p
                                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                                fac_PNP = 1.d0/(L_n*L_p)*emiNphi*emiZphi
                                do m1 = a_index, a_index+2*nja-1
                                    nma = BS%HO_sph%nljm(i0sp+m1,4) ! m_j + 1/2  
                                    do m2 = b_index, b_index+2*njb-1
                                        nmb = BS%HO_sph%nljm(i0sp+m2,4) ! m_j + 1/2
                                        norm = Proj_densities%norm(phi_n_index,1,ialpha,ibeta,igamma)*Proj_densities%norm(phi_p_index,2,ialpha,ibeta,igamma)
                                        pnorm = Proj_densities%pnorm(phi_n_index,1,ialpha,ibeta,igamma)*Proj_densities%pnorm(phi_p_index,2,ialpha,ibeta,igamma)
                                        ! ----------- q1-q2 (qf-qi) ---------------- 
                                        fac = fac_AMP*fac_PNP*(2*Jf+1)*(-1)**(Jf-Kf)*wigner3j(Jf,lambda,Ji,-Kf,mu,K1,IS=0)* &
                                            (-1)**(njb-nmb)*(-1)**(lambda+mu+1)*dsqrt(2*lambda+1.d0)*wigner3j(nja,lambda,njb,nma,-mu,1-nmb,IS=1)! (-1)^{jb-mb}C^{lambda mu}_{ja ma jb -mb}
                                        ! neutron part
                                        it = 1 
                                        reduced_TDME1B(it) = reduced_TDME1B(it) + fac*(norm*Proj_densities%rho_mm(m2,m1,ifg,phi_n_index,it,ialpha,ibeta,igamma)+ &
                                                        Parity_i*pnorm*Proj_densities%prho_mm(m2,m1,ifg,phi_n_index,it,ialpha,ibeta,igamma))/2.0d0
                                        ! proton part
                                        it = 2
                                        reduced_TDME1B(it) = reduced_TDME1B(it) + fac*(norm*Proj_densities%rho_mm(m2,m1,ifg,phi_p_index,it,ialpha,ibeta,igamma)+ &
                                                        Parity_i*pnorm*Proj_densities%prho_mm(m2,m1,ifg,phi_p_index,it,ialpha,ibeta,igamma))/2.0d0  
                                        
                                        !------------- q2-q1(qi-qf) -----------------
                                        if(q2_q1_Symmetry) then
                                            do nu = -lambda,lambda
                                                fac = cfac_AMP*(2*Jf+1)*(-1)**(Jf-Kf)*wigner3j(Jf,lambda,Ji,-Kf,mu,K1,IS=0)* &
                                                    (-1)**(mu-nu)*djmk(lambda,-nu,-mu,dcos(beta),IS=0)*CDEXP(nu*calpha+mu*cgamma)*&
                                                    (-1)**(njb-nmb)*(-1)**(lambda+nu+1)*dsqrt(2*lambda+1.d0)*wigner3j(nja,lambda,njb,nma,-nu,1-nmb,IS=1)
                                                ! neutron part
                                                it = 1 
                                                reduced_TDME1B_c(it) = reduced_TDME1B_c(it) + fac*DCONJG(fac_PNP*norm*Proj_densities%rho_mm(m1,m2,ifg,phi_n_index,it,ialpha,ibeta,igamma)+ &
                                                                Parity_f*fac_PNP*pnorm*Proj_densities%prho_mm(m1,m2,ifg,phi_n_index,it,ialpha,ibeta,igamma))/2.0d0
                                                ! proton part
                                                it = 2
                                                reduced_TDME1B_c(it) = reduced_TDME1B_c(it) + fac*DCONJG(fac_PNP*norm*Proj_densities%rho_mm(m1,m2,ifg,phi_p_index,it,ialpha,ibeta,igamma)+ &
                                                                Parity_f*fac_PNP*pnorm*Proj_densities%prho_mm(m1,m2,ifg,phi_p_index,it,ialpha,ibeta,igamma))/2.0d0
                                            end do 
                                        end if 
                                    end do 
                                end do
                            end do
                        end do
                    end do 
                end do 
            end do 
        end do       
        
        if(Proj_option%Euler_Symmetry==1 .and. Proj_option%AMPtype==1) then
            ! If we have implemented the symmetry of rho_mm, No need to multiply by 2.
            reduced_TDME1B(1) = fac_Parity*reduced_TDME1B(1)*2.d0
            reduced_TDME1B(2) = fac_Parity*reduced_TDME1B(2)*2.d0  
            reduced_TDME1B_c(1) = fac_Parity*reduced_TDME1B_c(1)*2.d0
            reduced_TDME1B_c(2) = fac_Parity*reduced_TDME1B_c(2)*2.d0  
        else if(Proj_option%Euler_Symmetry==0) then
            reduced_TDME1B(1) = fac_Parity*reduced_TDME1B(1)
            reduced_TDME1B(2) = fac_Parity*reduced_TDME1B(2)
            reduced_TDME1B_c(1) = fac_Parity*reduced_TDME1B_c(1)
            reduced_TDME1B_c(2) = fac_Parity*reduced_TDME1B_c(2)
        else 
            write(*,*) 'AMPtype=',Proj_option%AMPtype
            write(*,*) 'Euler_Symmetry=',Proj_option%Euler_Symmetry
            stop "Wrong AMPtype or Euler_Symmetry ! "
        end if 
    end subroutine 

end Module TD