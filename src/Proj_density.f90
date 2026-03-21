Module Proj_Density
    !------------------------------------------------------------------------------
    ! This module calculates the density matrix element after projection
    !-------------------------------------------------------------------------------
    use Constants, only: r64
    use Globals, only: Proj_densities
    implicit none
contains
    subroutine store_mix_density_matrix_elements(ialpha,ibeta,igamma)
        !---------------------------------------------------
        ! store mix density matrix elements
        !-------------------------------------------------
        use Constants, only: itx
        use Globals, only: BS,mix,projection_mesh,Proj_option,Proj_densities
        integer,intent(in) :: ialpha,ibeta,igamma
        integer :: dim_m_max, nphi_max, nalpha, nbeta, ngamma
        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        nphi_max = max(projection_mesh%nphi(1),projection_mesh%nphi(2))
        nalpha = projection_mesh%nalpha
        nbeta = projection_mesh%nbeta
        ngamma = projection_mesh%ngamma
        if(.not. allocated(Proj_densities%norm)) allocate(Proj_densities%norm(nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0))
        if(.not. allocated(Proj_densities%pnorm)) allocate(Proj_densities%pnorm(nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0))
        if(.not. allocated(Proj_densities%rho_mm)) allocate(Proj_densities%rho_mm(dim_m_max,dim_m_max,4,nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0)) !(m,m',++/--,phi_it,it,alpha,beta,gamma)
        if(.not. allocated(Proj_densities%prho_mm)) allocate(Proj_densities%prho_mm(dim_m_max,dim_m_max,4,nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0))
        if(.not. allocated(Proj_densities%kappa01c_mm)) allocate(Proj_densities%kappa01c_mm(dim_m_max,dim_m_max,4,nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0))
        if(.not. allocated(Proj_densities%pkappa01c_mm)) allocate(Proj_densities%pkappa01c_mm(dim_m_max,dim_m_max,4,nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0))
        if(.not. allocated(Proj_densities%kappa10_mm)) allocate(Proj_densities%kappa10_mm(dim_m_max,dim_m_max,4,nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0))
        if(.not. allocated(Proj_densities%pkappa10_mm)) allocate(Proj_densities%pkappa10_mm(dim_m_max,dim_m_max,4,nphi_max,itx,nalpha,nbeta,ngamma),source=(0.d0,0.d0))

        ! we have to implement the symmetry of rho_mm
        if(Proj_option%Euler_Symmetry>=1 .and. ialpha > (nalpha+1)/2 ) then
            write(*,*) '[alpha symmetry of rho_mm] Not yet implemented.... '
            return 
        end if
        if(Proj_option%Euler_Symmetry>=1 .and. ibeta > (nbeta+1)/2 ) then
            if(Proj_option%Euler_Symmetry == 1) then
                ! write(*,*) '[beta symmetry of rho_mm] Not yet implemented  .... ' 
                ! we can simply multiply by 2 without integrating over [pi/2, pi ]
            else if(Proj_option%Euler_Symmetry == 2) then
                write(*,*) '[beta symmetry of rho_mm] Not yet implemented  .... ' 
            end if 
            return
        end if 

        Proj_densities%norm(:,:,ialpha,ibeta,igamma) = mix%norm(:,:)
        Proj_densities%pnorm(:,:,ialpha,ibeta,igamma) = mix%pnorm(:,:)
        Proj_densities%rho_mm(:,:,1:4,:,:,ialpha,ibeta,igamma) = mix%rho_mm(:,:,1:4,:,:)
        Proj_densities%prho_mm(:,:,1:4,:,:,ialpha,ibeta,igamma) = mix%prho_mm(:,:,1:4,:,:)
        Proj_densities%kappa01c_mm(:,:,1:4,:,:,ialpha,ibeta,igamma) = mix%kappa01c_mm(:,:,1:4,:,:)
        Proj_densities%pkappa01c_mm(:,:,1:4,:,:,ialpha,ibeta,igamma) = mix%pkappa01c_mm(:,:,1:4,:,:)
        Proj_densities%kappa10_mm(:,:,1:4,:,:,ialpha,ibeta,igamma) = mix%kappa10_mm(:,:,1:4,:,:)
        Proj_densities%pkappa10_mm(:,:,1:4,:,:,ialpha,ibeta,igamma) = mix%pkappa10_mm(:,:,1:4,:,:)
    end subroutine

    subroutine calculate_density_matrix_element(q1,q2)
        !-----------------------------------------------------------------
        !
        !   calculate and store the denstiy matrix elements
        !
        !  Note:
        !      Due to the extremely large memory usage, the two-body density
        !      is not stored in Proj_densities%ME2B.
        !-----------------------------------------------------------------
        use Globals, only: gcm_space,Proj_option,BS,kernels
        integer, intent(in) :: q1,q2
        integer :: dim_m_max,J,K1_start,K1_end,K2_start,K2_end,K1,K2,iParity,Parity,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,total_iter,iter
        complex(r64) :: ME1B(2),ME2B(2)
        complex(r64) :: N(2), N2(2)
        ! logical :: q2_q1_Symmetry
        write(*,'(5x,A)') 'calculate_density_matrix_element ...'
        ! if(q1/=q2 .and. Proj_option%Kernel_Symmetry==1) then
        !     q2_q1_Symmetry = .True.
        ! else
        !     q2_q1_Symmetry = .False.
        ! end if

        dim_m_max = max(BS%HO_sph%idsp(1,1), BS%HO_sph%idsp(1,2))
        if(Proj_option%AMPtype==0 .or. Proj_option%AMPtype==1) then
            if(.not. allocated(Proj_densities%ME1B)) allocate(Proj_densities%ME1B(gcm_space%Jmin:gcm_space%Jmax, 0:0, 0:0, 2,& ! Jf, K, K', Pi(+/-)
                                                             2, 2, dim_m_max,dim_m_max)) ! it, ifg(++/--), m, m'
            ! if(.not. allocated(Proj_densities%ME2B)) allocate(Proj_densities%ME2B(gcm_space%Jmin:gcm_space%Jmax, 0:0, 0:0, 2,& ! Jf, K, K', Pi(+/-)
            !                                                  2, 2, dim_m_max,dim_m_max,dim_m_max,dim_m_max)) ! it, ifg(++/--), m1, m2, m3, m4
        else
            if(.not. allocated(Proj_densities%ME1B)) allocate(Proj_densities%ME1B(gcm_space%Jmin:gcm_space%Jmax, -gcm_space%Jmax:gcm_space%Jmax, -gcm_space%Jmax:gcm_space%Jmax, 2,& ! Jf, K, K', Pi(+/-)
                                                             2, 2, dim_m_max,dim_m_max)) ! it, ifg(++/--), m, m'
            ! if(.not. allocated(Proj_densities%ME2B)) allocate(Proj_densities%ME2B(gcm_space%Jmin:gcm_space%Jmax, -gcm_space%Jmax:gcm_space%Jmax, -gcm_space%Jmax:gcm_space%Jmax, 2,& ! Jf, K, K', Pi(+/-)
            !                                                  2, 2, dim_m_max,dim_m_max,dim_m_max,dim_m_max)) ! it, ifg(++/--), m1, m2, m3, m4
        end if
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
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
                        N = 0.d0
                        N2 = 0.d0 
                        total_iter = BS%HO_sph%idsp(1,1) + BS%HO_sph%idsp(1,2)
                        !$OMP PARALLEL DEFAULT(shared) PRIVATE(iter,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,ME1B,ME2B) &
                        !$omp reduction(+:N,N2)
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
                                if(Proj_option%DsType == 1 .or. Proj_option%DsType == 3) then 
                                    ! 1 Body
                                    if(ifg1==ifg2) then
                                        call calculate_one_body_density_matrix_element(J,K1,K2,Parity,ifg1,m1,ifg2,m2,ME1B)
                                        ! neutron
                                        Proj_densities%ME1B(J,K1,K2,iParity,1,ifg1,m1,m2) = ME1B(1)
                                        ! proton
                                        Proj_densities%ME1B(J,K1,K2,iParity,2,ifg1,m1,m2) = ME1B(2)
                                        ! check particle number
                                        if (ifg1==ifg2 .and. m1==m2) then
                                            N(1) = N(1) + Proj_densities%ME1B(J,K1,K2,iParity,1,ifg1,m1,m2)
                                            N(2) = N(2) + Proj_densities%ME1B(J,K1,K2,iParity,2,ifg1,m1,m2)
                                        end if
                                    end if
                                end if
                                if(Proj_option%DsType == 2 .or. Proj_option%DsType == 3) then
                                    do ifg3 = 1, 2
                                    do m3 = 1, BS%HO_sph%idsp(1,ifg3)
                                        do ifg4 = 1, 2
                                        do m4 = 1, BS%HO_sph%idsp(1,ifg4)
                                            ! 2 Body
                                            ! call calculate_two_body_density_matrix_element(J,K1,K2,Parity,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,ME2B)
                                            ! ! neutron
                                            ! Proj_densities%ME2B(J,K1,K2,iParity,1,ifg1,m1,m2,m3,m4) = Real(ME2B(1))
                                            ! ! proton
                                            ! Proj_densities%ME2B(J,K1,K2,iParity,2,ifg1,m1,m2,m3,m4) = Real(ME2B(2))
                                            ! check particle number
                                            if ((ifg1==ifg3 .and. m1==m3) .and. (ifg2==ifg4 .and. m2==m4))  then
                                                call calculate_two_body_density_matrix_element(J,K1,K2,Parity,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,ME2B)
                                                N2(1) = N2(1) + ME2B(1)
                                                N2(2) = N2(2) + ME2B(2)
                                            end if 
                                        end do 
                                        end do
                                    end do
                                    end do 
                                end if 
                            end do
                            end do
                        end do
                        !$OMP END PARALLEL
                        ! write(*,"(8x,4(a,i3,4x),a,f15.10)") "J=",J,"K=",K1,"K'=",K2,'Parity=',Parity,'    N     =',Real(N(1)/(kernels%N_KK(J,K1,K2,iParity)+1.0d-30))
                        ! write(*,"(8x,4(a,i3,4x),a,f15.10)") "J=",J,"K=",K1,"K'=",K2,'Parity=',Parity,'    Z     =',Real(N(2)/(kernels%N_KK(J,K1,K2,iParity)+1.0d-30))
                        ! write(*,"(8x,4(a,i3,4x),a,f15.10)") "J=",J,"K=",K1,"K'=",K2,'Parity=',Parity,'    N(1-N)=',Real(N2(1)/(kernels%N_KK(J,K1,K2,iParity)+1.0d-30))
                        ! write(*,"(8x,4(a,i3,4x),a,f15.10)") "J=",J,"K=",K1,"K'=",K2,'Parity=',Parity,'    Z(1-Z)=',Real(N2(2)/(kernels%N_KK(J,K1,K2,iParity)+1.0d-30))
                    end do 
                end do 
            end do 
        end do
    end subroutine

    subroutine calculate_one_body_density_matrix_element(J,K1,K2,Parity,ifg1,m1,ifg2,m2,ME1B)
        !--------------------------------------------------------------------------------------------------
        !
        !   calculate reduced 1 body density matrix element
        !        rho_{m1 m2}^{J K1 K2,NZ, Parity} = <q1|c^+_{m2}c_{m1} P^{J}_{K1 K2} P^N P^Z P^{Parity}|q2>
        ! 
        !--------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: projection_mesh,Proj_option,nucleus_attributes
        use Basis, only: djmk
        integer,intent(in) :: J,K1,K2,Parity,ifg1,m1,ifg2,m2
        complex(r64), intent(out) :: ME1B(2)
        integer :: ialpha,ibeta,igamma,L_n,L_p,phi_n_index,phi_p_index,it
        real(r64) :: alpha,beta,gamma,w,phi_n,phi_p
        complex(r64) :: calpha,cgamma,fac_AMP,cpi,fac1,fac2,emiNphi,emiZphi,fac_PNP,norm,pnorm
        ME1B = (0.d0, 0.d0) 
        do ialpha = 1, projection_mesh%nalpha
            do ibeta = 1, projection_mesh%nbeta
                do igamma = 1, projection_mesh%ngamma
                    ! If we have implemented the symmetry of rho_mm, then these lines are not needed.
                    if(Proj_option%Euler_Symmetry == 2) then
                        stop '[calculate_one_body_density_matrix_element]: Euler_Symmetry=2, Not yet implemented! You should set the Symmetry of Euler angles as 0!'
                    end if
                    if(Proj_option%Euler_Symmetry==1 .and. ibeta>(projection_mesh%nbeta+1)/2) then
                        cycle
                        ! Using the tensor symmetry, it can be proven that when Parity = (-1)**J
                        !  the contribution from (pi/2, pi]) is the same as that from (0, pi/2).
                        ! stop '[calculate_one_body_density_matrix_element]: Euler_Symmetry=1, Not yet implemented! You should set the Symmetry of Euler angles as 0!'
                    end if

                    alpha = projection_mesh%alpha(ialpha)
                    calpha = DCMPLX(0.d0,alpha)
                    beta = projection_mesh%beta(ibeta)
                    gamma = projection_mesh%gamma(igamma)
                    cgamma = DCMPLX(0.d0,gamma)
                    if(Proj_option%AMPtype==0) then
                        fac_AMP = 1
                    else if (Proj_option%AMPtype==1) then
                        w = projection_mesh%wbeta(ibeta)
                        fac1 = (2*J+1)/(2.0d0)*dsin(beta)*djmk(J,K1,K2,dcos(beta),0)
                        fac_AMP = fac1*w
                    else
                        cpi = DCMPLX(0.d0,pi) ! i*pi
                        w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                        fac1 = (2*J+1)/(8.0d0*pi**2)*dsin(beta)*djmk(J,K1,K2,dcos(beta),0)*CDEXP(-K1*calpha-K2*cgamma)
                        fac2 = 1.0d0 + CDEXP(-K1*cpi) + CDEXP(-K2*cpi) + CDEXP(-K1*cpi-K2*cpi) ! D2 symmetry is required, with alpha, gamma in [0, pi].
                        fac_AMP = fac1*fac2*w
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
                            norm = Proj_densities%norm(phi_n_index,1,ialpha,ibeta,igamma)*Proj_densities%norm(phi_p_index,2,ialpha,ibeta,igamma)
                            pnorm = Proj_densities%pnorm(phi_n_index,1,ialpha,ibeta,igamma)*Proj_densities%pnorm(phi_p_index,2,ialpha,ibeta,igamma)
                            ! neutron part
                            it = 1 
                            ME1B(it) = ME1B(it) + fac_AMP*fac_PNP*(norm*Proj_densities%rho_mm(m1,m2,indexfg(ifg1,ifg2),phi_n_index,it,ialpha,ibeta,igamma)+ &
                                                    Parity*pnorm*Proj_densities%prho_mm(m1,m2,indexfg(ifg1,ifg2),phi_n_index,it,ialpha,ibeta,igamma))/2.0d0
                            ! proton part
                            it = 2
                            ME1B(it) = ME1B(it) + fac_AMP*fac_PNP*(norm*Proj_densities%rho_mm(m1,m2,indexfg(ifg1,ifg2),phi_p_index,it,ialpha,ibeta,igamma)+ &
                                                    Parity*pnorm*Proj_densities%prho_mm(m1,m2,indexfg(ifg1,ifg2),phi_p_index,it,ialpha,ibeta,igamma))/2.0d0
                        end do
                    end do
                end do
            end do
        end do 

        if(Proj_option%Euler_Symmetry==1 .and. Proj_option%AMPtype==1 ) then
            ! If we have implemented the symmetry of rho_mm, No need to multiply by 2.
            ME1B(1) = ME1B(1)*2.d0
            ME1B(2) = ME1B(2)*2.d0   
        else if(Proj_option%Euler_Symmetry==0 .or. Proj_option%AMPtype==0) then
            ME1B(1) = ME1B(1)
            ME1B(2) = ME1B(2)
        else 
            write(*,*) 'AMPtype=',Proj_option%AMPtype
            write(*,*) 'Euler_Symmetry=',Proj_option%Euler_Symmetry
            stop "Wrong AMPtype or Euler_Symmetry ! "
        end if 
    end subroutine

    subroutine calculate_two_body_density_matrix_element(J,K1,K2,Parity,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4,ME2B)
        !---------------------------------------------------------------------------------------------------------------
        !
        !   calculate reduced 2 body density matrix element
        !      rho_{m1 m2 m3 m4}^{J K1 K2,NZ, Parity} = <q1|c^+_{m4}c^+_{m3}c_{m2}c_{m1}P^{J}_{K1 K2}P^NP^ZP^{Parity}|q2>
        ! 
        !----------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: projection_mesh,Proj_option,nucleus_attributes
        use Basis, only: djmk
        integer,intent(in) :: J,K1,K2,Parity,ifg1,m1,ifg2,m2,ifg3,m3,ifg4,m4
        complex(r64), intent(out) :: ME2B(2)
        integer :: ialpha,ibeta,igamma,L_n,L_p,phi_n_index,phi_p_index,it
        real(r64) :: alpha,beta,gamma,w,phi_n,phi_p
        complex(r64) :: calpha,cgamma,fac_AMP,cpi,fac1,fac2,emiNphi,emiZphi,fac_PNP,norm,pnorm
        ME2B = (0.d0, 0.d0)
        do ialpha = 1, projection_mesh%nalpha
            do ibeta = 1, projection_mesh%nbeta
                do igamma = 1, projection_mesh%ngamma
                    ! If we have implemented the symmetry of rho_mm, then these lines are not needed.
                    if(Proj_option%Euler_Symmetry == 2) then
                        stop '[calculate_two_body_density_matrix_element]: Euler_Symmetry=2, Not yet implemented! You should set the Symmetry of Euler angles as 0!'
                    end if
                    if(Proj_option%Euler_Symmetry==1 .and. ibeta>(projection_mesh%nbeta+1)/2) then
                        cycle
                        ! Using the tensor symmetry, it can be proven that when Parity = (-1)**J
                        !  the contribution from (pi/2, pi]) is the same as that from (0, pi/2).
                        ! stop '[calculate_two_body_density_matrix_element]: Euler_Symmetry=1, Not yet implemented! You should set the Symmetry of Euler angles as 0!'
                    end if

                    alpha = projection_mesh%alpha(ialpha)
                    calpha = DCMPLX(0.d0,alpha)
                    beta = projection_mesh%beta(ibeta)
                    gamma = projection_mesh%gamma(igamma)
                    cgamma = DCMPLX(0.d0,gamma)
                    if(Proj_option%AMPtype==0) then
                        fac_AMP = 1
                    else if (Proj_option%AMPtype==1) then
                        w = projection_mesh%wbeta(ibeta)
                        fac1 = (2*J+1)/(2.0d0)*dsin(beta)*djmk(J,K1,K2,dcos(beta),0)
                        fac_AMP = fac1*w
                    else
                        cpi = DCMPLX(0.d0,pi) ! i*pi
                        w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                        fac1 = (2*J+1)/(8.0d0*pi**2)*dsin(beta)*djmk(J,K1,K2,dcos(beta),0)*CDEXP(-K1*calpha-K2*cgamma)
                        fac2 = 1.0d0 + CDEXP(-K1*cpi) + CDEXP(-K2*cpi) + CDEXP(-K1*cpi-K2*cpi) ! D2 symmetry is required, with alpha, gamma in [0, pi].
                        fac_AMP = fac1*fac2*w
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
                            norm = Proj_densities%norm(phi_n_index,1,ialpha,ibeta,igamma)*Proj_densities%norm(phi_p_index,2,ialpha,ibeta,igamma)
                            pnorm = Proj_densities%pnorm(phi_n_index,1,ialpha,ibeta,igamma)*Proj_densities%pnorm(phi_p_index,2,ialpha,ibeta,igamma)
                            ! neutron part
                            it = 1
                            ME2B(it) = ME2B(it) + fac_AMP*fac_PNP* &
                                            (norm* &
                                            (Proj_densities%rho_mm(m1,m4,indexfg(ifg1,ifg4),phi_n_index,it,ialpha,ibeta,igamma)*Proj_densities%rho_mm(m2,m3,indexfg(ifg2,ifg3),phi_n_index,it,ialpha,ibeta,igamma) &
                                            -Proj_densities%rho_mm(m2,m4,indexfg(ifg2,ifg4),phi_n_index,it,ialpha,ibeta,igamma)*Proj_densities%rho_mm(m1,m3,indexfg(ifg1,ifg3),phi_n_index,it,ialpha,ibeta,igamma) &
                                            +Proj_densities%kappa01c_mm(m4,m3,indexfg(ifg4,ifg3),phi_n_index,it,ialpha,ibeta,igamma)*Proj_densities%kappa10_mm(m1,m2,indexfg(ifg1,ifg2),phi_n_index,it,ialpha,ibeta,igamma)) &
                                            +Parity*pnorm* &
                                            (Proj_densities%prho_mm(m1,m4,indexfg(ifg1,ifg4),phi_n_index,it,ialpha,ibeta,igamma)*Proj_densities%prho_mm(m2,m3,indexfg(ifg2,ifg3),phi_n_index,it,ialpha,ibeta,igamma) &
                                            -Proj_densities%prho_mm(m2,m4,indexfg(ifg2,ifg4),phi_n_index,it,ialpha,ibeta,igamma)*Proj_densities%prho_mm(m1,m3,indexfg(ifg1,ifg3),phi_n_index,it,ialpha,ibeta,igamma) &
                                            +Proj_densities%pkappa01c_mm(m4,m3,indexfg(ifg4,ifg3),phi_n_index,it,ialpha,ibeta,igamma)*Proj_densities%pkappa10_mm(m1,m2,indexfg(ifg1,ifg2),phi_n_index,it,ialpha,ibeta,igamma)) &
                                            )/2.0d0
                            ! proton part
                            it = 2
                            ME2B(it) = ME2B(it) + fac_AMP*fac_PNP* &
                                            (norm* &
                                            (Proj_densities%rho_mm(m1,m4,indexfg(ifg1,ifg4),phi_p_index,it,ialpha,ibeta,igamma)*Proj_densities%rho_mm(m2,m3,indexfg(ifg2,ifg3),phi_p_index,it,ialpha,ibeta,igamma) &
                                            -Proj_densities%rho_mm(m2,m4,indexfg(ifg2,ifg4),phi_p_index,it,ialpha,ibeta,igamma)*Proj_densities%rho_mm(m1,m3,indexfg(ifg1,ifg3),phi_p_index,it,ialpha,ibeta,igamma) &
                                            +Proj_densities%kappa01c_mm(m4,m3,indexfg(ifg4,ifg3),phi_p_index,it,ialpha,ibeta,igamma)*Proj_densities%kappa10_mm(m1,m2,indexfg(ifg1,ifg2),phi_p_index,it,ialpha,ibeta,igamma)) &
                                            +Parity*pnorm* &
                                            (Proj_densities%prho_mm(m1,m4,indexfg(ifg1,ifg4),phi_p_index,it,ialpha,ibeta,igamma)*Proj_densities%prho_mm(m2,m3,indexfg(ifg2,ifg3),phi_p_index,it,ialpha,ibeta,igamma) &
                                            -Proj_densities%prho_mm(m2,m4,indexfg(ifg2,ifg4),phi_p_index,it,ialpha,ibeta,igamma)*Proj_densities%prho_mm(m1,m3,indexfg(ifg1,ifg3),phi_p_index,it,ialpha,ibeta,igamma) &
                                            +Proj_densities%pkappa01c_mm(m4,m3,indexfg(ifg4,ifg3),phi_p_index,it,ialpha,ibeta,igamma)*Proj_densities%pkappa10_mm(m1,m2,indexfg(ifg1,ifg2),phi_p_index,it,ialpha,ibeta,igamma)) &
                                            )/2.0d0
                        end do
                    end do
                end do
            end do
        end do

        if(Proj_option%Euler_Symmetry==1 .and. Proj_option%AMPtype==1) then
            ! If we have implemented the symmetry of rho_mm, No need to multiply by 2.
            ME2B(1) = ME2B(1)*2.d0
            ME2B(2) = ME2B(2)*2.d0   
        else if(Proj_option%Euler_Symmetry==0 .or. Proj_option%AMPtype==0) then
            ME2B(1) = ME2B(1)
            ME2B(2) = ME2B(2)
        else 
            write(*,*) 'AMPtype=',Proj_option%AMPtype
            write(*,*) 'Euler_Symmetry=',Proj_option%Euler_Symmetry
            stop "Wrong AMPtype or Euler_Symmetry ! "
        end if 
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

End Module Proj_Density