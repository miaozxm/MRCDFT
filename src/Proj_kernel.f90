!==============================================================================!
! MODULE Kernel                                                                !
!                                                                              !
! This module calculates  kernels                                              !
!          <J K_1 q_1 Pi | O | J K_2 q_2 Pi >   where O = 1 or N or H          !
!     <J_f K_f q_1 Pi_f ||Q_2|| J_i K_i q_2 Pi_i >                             !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
Module Kernel
    use Constants, only: r64
    implicit none
    complex(r64),allocatable :: Norm_PNP_AMParray(:,:,:), pNorm_PNP_AMParray(:,:,:), &
                                Etot_PNP_AMParray(:,:,:), pEtot_PNP_AMParray(:,:,:), &
                                Particle_PNP_AMParray(:,:,:,:), pParticle_PNP_AMParray(:,:,:,:), &
                                N2_PNP_AMParray(:,:,:,:), pN2_PNP_AMParray(:,:,:,:), &
                                J2_PNP_AMParray(:,:,:), pJ2_PNP_AMParray(:,:,:), &
                                Q2m_PNP_AMParray(:,:,:,:,:),pQ2m_PNP_AMParray(:,:,:,:,:), &
                                cQ2m_PNP_AMParray(:,:,:,:,:),pcQ2m_PNP_AMParray(:,:,:,:,:), &
                                r2_PNP_AMParray(:,:,:,:), pr2_PNP_AMParray(:,:,:,:),&
                                Eccentri_PNP_AMParray(:,:,:,:), pEccentri_PNP_AMParray(:,:,:,:),&
                                r2_2b_PNP_AMParray(:,:,:,:,:), pr2_2b_PNP_AMParray(:,:,:,:,:)


    contains
    
    subroutine set_projection_mesh_points()
        use Constants, only: pi
        use Globals, only: input_par,Proj_option,projection_mesh,gcm_space
        use MathMethods, only: GaussLegendre_X1toX2
        integer :: it
        ! Set the parameters for the PNP mesh (gauge angles)
        ! for even-even nuclei, N and Z are even number, so we can reduce [0, 2pi] to [0,pi]
        projection_mesh%nphi(1) = input_par%nphi
        projection_mesh%nphi(2) = input_par%nphi ! same as neutron
        projection_mesh%dphi(1) = pi/projection_mesh%nphi(1)
        projection_mesh%dphi(2) = pi/projection_mesh%nphi(2)
        if(Proj_option%PNPtype==0) then 
            ! no PNP
            do it = 1,2
                projection_mesh%nphi(it) = 1
                projection_mesh%dphi(it) = 0.0d0
            end do 
            ! write(*,*) 'no PNP'
        end if

        ! Set the parameters for the AMP mesh (Euler angles)
        projection_mesh%nalpha = input_par%nalpha
        projection_mesh%nbeta = input_par%nbeta
        projection_mesh%ngamma = input_par%ngamma
        if(.not. allocated(projection_mesh%alpha)) allocate(projection_mesh%alpha(projection_mesh%nalpha))
        if(.not. allocated(projection_mesh%walpha)) allocate(projection_mesh%walpha(projection_mesh%nalpha))
        if(.not. allocated(projection_mesh%beta)) allocate(projection_mesh%beta(projection_mesh%nbeta))
        if(.not. allocated(projection_mesh%wbeta)) allocate(projection_mesh%wbeta(projection_mesh%nbeta))
        if(.not. allocated(projection_mesh%gamma)) allocate(projection_mesh%gamma(projection_mesh%ngamma))
        if(.not. allocated(projection_mesh%wgamma)) allocate(projection_mesh%wgamma(projection_mesh%ngamma))
        call GaussLegendre_X1toX2(0.d0,pi,projection_mesh%alpha,projection_mesh%walpha,projection_mesh%nalpha) ! reduce [0, 2pi] to [0,pi] ! D2 symmetry is required.
        call GaussLegendre_X1toX2(0.d0,pi,projection_mesh%beta,projection_mesh%wbeta,projection_mesh%nbeta) ! [0,pi]
        call GaussLegendre_X1toX2(0.d0,pi,projection_mesh%gamma,projection_mesh%wgamma,projection_mesh%ngamma)  ! reduce [0, 2pi] to [0,pi] ! D2 symmetry is required.
        if(Proj_option%AMPtype==0) then
            ! alpha
            projection_mesh%nalpha = 1
            projection_mesh%alpha(1) = 0.d0
            projection_mesh%walpha(1) = pi
            ! beta
            projection_mesh%nbeta = 1
            projection_mesh%beta(1) = 0.d0
            projection_mesh%wbeta(1) = pi
            ! gamma
            projection_mesh%ngamma = 1
            projection_mesh%gamma(1) = 0.d0
            projection_mesh%wgamma(1) = pi
            ! J
            gcm_space%Jmin = 0
            gcm_space%Jmax = 0
            ! write(*,*) 'no AMP'
        else if(Proj_option%AMPtype==1) then! 1DAMP
            ! alpha
            projection_mesh%nalpha = 1
            projection_mesh%alpha(1) = 0.d0
            projection_mesh%walpha(1) = pi
            ! 
            ! gamma
            projection_mesh%ngamma = 1
            projection_mesh%gamma(1) = 0.d0
            projection_mesh%wgamma(1) = pi
            ! write(*,*) '1DAMP'
            
        else 
            ! write(*,*) '3DAMP'
        end if

    end subroutine
    
    !----------------------------------------
    !   AMP Integration
    !----------------------------------------
    subroutine calculate_Kernel
        !----------------------------------------------------------------------
        ! 1D-AMP+GCM for axial-quadrupole deformed states  
        !      <q1| O R(beta)|q2> and <q1|O R(beta) P|q2>
        ! where R is rotation operator and P is parity operator
        !----------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: projection_mesh,Proj_option
        integer :: nalpha,nbeta,ngamma

        nalpha = projection_mesh%nalpha
        nbeta = projection_mesh%nbeta
        ngamma = projection_mesh%ngamma
        allocate(Norm_PNP_AMParray(nalpha,nbeta,ngamma), pNorm_PNP_AMParray(nalpha,nbeta,ngamma), &
                 Etot_PNP_AMParray(nalpha,nbeta,ngamma), pEtot_PNP_AMParray(nalpha,nbeta,ngamma), &
                 Particle_PNP_AMParray(nalpha,nbeta,ngamma,2), pParticle_PNP_AMParray(nalpha,nbeta,ngamma,2), &
                 N2_PNP_AMParray(nalpha,nbeta,ngamma,2), pN2_PNP_AMParray(nalpha,nbeta,ngamma,2), &
                 J2_PNP_AMParray(nalpha,nbeta,ngamma), pJ2_PNP_AMParray(nalpha,nbeta,ngamma), &
                 Q2m_PNP_AMParray(nalpha,nbeta,ngamma,-2:2,2),pQ2m_PNP_AMParray(nalpha,nbeta,ngamma,-2:2,2),&
                 cQ2m_PNP_AMParray(nalpha,nbeta,ngamma,-2:2,2),pcQ2m_PNP_AMParray(nalpha,nbeta,ngamma,-2:2,2),&
                 r2_PNP_AMParray(nalpha,nbeta,ngamma,2),pr2_PNP_AMParray(nalpha,nbeta,ngamma,2),&
                 Eccentri_PNP_AMParray(nalpha,nbeta,ngamma,5),pEccentri_PNP_AMParray(nalpha,nbeta,ngamma,5),&
                 r2_2b_PNP_AMParray(nalpha,nbeta,ngamma,5,3),pr2_2b_PNP_AMParray(nalpha,nbeta,ngamma,5,3),source=(0.0d0,0.0d0))

        call calculate_overlaps_arrays

        if(Proj_option%PPtype==0) then 
            ! no PP
            call reset_parityOverlap
        end if
        

        ! Integration over AMP mesh points(alpha beta gamma)
        call calcualate_Norm_Hamiltonian_ParticleNumber_kernels  ! calculate <J K_1 q_1 Pi | O | J K_2 q_2 Pi >
        if( Proj_option%checkN2J2 == 1) then 
            call calculate_N2_kernels ! <J   K_f q_1 Pi  | N^2 |J   K_i q_2 Pi>
            call calculate_J2_kernels ! <J   K_f q_1 Pi  | J^2 |J   K_i q_2 Pi>
        end if 
        call calculate_EM_kernels ! calcualate <J_f K_f q_1 Pi_f ||T_lambda|| J_i K_i q_2 Pi_i >
        call calculate_E0_kernel ! <J   K_f q_1 Pi  | r^2 |J   K_i q_2 Pi>
        call calculate_r2_2body_kernel ! <J   K_f q_1 Pi  | r^2_2b |J   K_i q_2 Pi>
        if( Proj_option%EccentriType == 1 .or. Proj_option%EccentriType == 3) then
            call calculate_Eccentricity_kernel ! <J   K_f q_1 Pi  | E_n |J   K_i q_2 Pi>
        end if 
        
        deallocate( Norm_PNP_AMParray,pNorm_PNP_AMParray, &
                    Etot_PNP_AMParray,pEtot_PNP_AMParray, &
                    Particle_PNP_AMParray,pParticle_PNP_AMParray, &
                    N2_PNP_AMParray, pN2_PNP_AMParray, &
                    J2_PNP_AMParray, pJ2_PNP_AMParray, &
                    Q2m_PNP_AMParray,pQ2m_PNP_AMParray,&
                    cQ2m_PNP_AMParray,pcQ2m_PNP_AMParray,&
                    r2_PNP_AMParray,pr2_PNP_AMParray,&
                    Eccentri_PNP_AMParray,pEccentri_PNP_AMParray,&
                    r2_2b_PNP_AMParray, pr2_2b_PNP_AMParray)
    end subroutine

    subroutine calculate_overlaps_arrays
        use Globals, only: projection_mesh,Proj_option,MPI_Infor
        use Proj_Density, only: store_mix_density_matrix_elements
        use Tools, only: adjust_left
        integer :: nalpha,nbeta,ngamma,ialpha,ibeta,igamma,mu
        real(r64) :: alpha, beta, gamma
        character(len=*),parameter :: format1 = "(5x,'alpha:',i3,'/',a,'beta:',i3,'/',a,'gamma:',i3,'/',a)"
        nalpha = projection_mesh%nalpha
        nbeta = projection_mesh%nbeta
        ngamma = projection_mesh%ngamma
        do ialpha = 1, nalpha ! loop of alpha
            alpha = projection_mesh%alpha(ialpha)
            do ibeta = 1, nbeta  ! loop of beta [0,pi]
                beta = projection_mesh%beta(ibeta)
                do igamma = 1, ngamma ! loop of gamma
                    gamma = projection_mesh%gamma(igamma)
                    if (MPI_Infor%rank == 0) write(*,format1,advance='no') ialpha,adjust_left(nalpha,3),ibeta,adjust_left(nbeta,3),igamma,adjust_left(ngamma,3)
                    !##########################################################
                    !#    using the symmetry(D2 and Axial+Parity) of the Euler angles
                    !##########################################################
                    ! D2
                    if(Proj_option%AMPtype==2 .and. Proj_option%Euler_Symmetry==2 .and. ialpha > (nalpha+1)/2 ) then
                        if (MPI_Infor%rank == 0) write(*,'(A)') '(alpha) symmetry(D2).'
                        ! because <O R(alpha,beta,gamma)> = <O R(pi-alpha,beta,pi-gamma)>, 
                        Norm_PNP_AMParray(ialpha,ibeta,igamma)  = Norm_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma)
                        pNorm_PNP_AMParray(ialpha,ibeta,igamma) = pNorm_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma)
                        Etot_PNP_AMParray(ialpha,ibeta,igamma) = Etot_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma)
                        pEtot_PNP_AMParray(ialpha,ibeta,igamma) = pEtot_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma)
                        Particle_PNP_AMParray(ialpha,ibeta,igamma,1) = Particle_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,1)
                        Particle_PNP_AMParray(ialpha,ibeta,igamma,2) = Particle_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,2) 
                        pParticle_PNP_AMParray(ialpha,ibeta,igamma,1) = pParticle_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,1)
                        pParticle_PNP_AMParray(ialpha,ibeta,igamma,2) = pParticle_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,2)
                        ! N^2
                        N2_PNP_AMParray(ialpha,ibeta,igamma,:) = N2_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:)
                        pN2_PNP_AMParray(ialpha,ibeta,igamma,:) = pN2_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:)
                        ! J^2
                        J2_PNP_AMParray(ialpha,ibeta,igamma) = J2_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma)
                        pJ2_PNP_AMParray(ialpha,ibeta,igamma) = pJ2_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma)
                        ! because <T_{lm} R(alpha,beta,gamma)> =(-1)^l*<T_{l-m} R(pi-alpha,beta,pi-gamma)>, 
                        ! Q2m
                        do mu =-2,2
                            Q2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**2*Q2m_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,-mu,:)
                            pQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**2*pQ2m_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,-mu,:)
                            cQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**2*cQ2m_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,-mu,:)
                            pcQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**2*pcQ2m_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,-mu,:)
                        end do 
                        ! r2
                        r2_PNP_AMParray(ialpha,ibeta,igamma,:) = r2_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:)
                        pr2_PNP_AMParray(ialpha,ibeta,igamma,:) = pr2_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:)
                        ! r2_2b
                        r2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:) = r2_2b_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:,:)
                        pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:) = pr2_2b_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:,:)
                        ! 
                        Eccentri_PNP_AMParray(ialpha,ibeta,igamma,:) = Eccentri_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:)
                        pEccentri_PNP_AMParray(ialpha,ibeta,igamma,:) = pEccentri_PNP_AMParray(nalpha+1-ialpha,ibeta,ngamma+1-igamma,:)
                        ! store density matrix elements
                        call store_mix_density_matrix_elements(ialpha,ibeta,igamma)
                        cycle
                    end if                   
                    if(Proj_option%AMPtype==2 .and. Proj_option%Euler_Symmetry==2 .and. ibeta > (nbeta+1)/2 ) then
                        if (MPI_Infor%rank == 0) write(*,'(A)') '(beta) symmetry(D2).'
                        ! because <O R(alpha,beta,gamma)> = <O R(alpha,pi-beta,pi-gamma)>, 
                        Norm_PNP_AMParray(ialpha,ibeta,igamma)  = Norm_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        pNorm_PNP_AMParray(ialpha,ibeta,igamma) = pNorm_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        Etot_PNP_AMParray(ialpha,ibeta,igamma) = Etot_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        pEtot_PNP_AMParray(ialpha,ibeta,igamma) = pEtot_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        Particle_PNP_AMParray(ialpha,ibeta,igamma,1) = Particle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,1)
                        Particle_PNP_AMParray(ialpha,ibeta,igamma,2) = Particle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,2) 
                        pParticle_PNP_AMParray(ialpha,ibeta,igamma,1) = pParticle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,1)
                        pParticle_PNP_AMParray(ialpha,ibeta,igamma,2) = pParticle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,2)
                        ! N^2
                        N2_PNP_AMParray(ialpha,ibeta,igamma,:) = N2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:)
                        pN2_PNP_AMParray(ialpha,ibeta,igamma,:) = pN2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:)
                        ! J^2
                        J2_PNP_AMParray(ialpha,ibeta,igamma) = J2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        pJ2_PNP_AMParray(ialpha,ibeta,igamma) = pJ2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        ! because <T_{lm} R(alpha,beta,gamma)> = (-1)^m*<T_{lm} R(alpha,pi-beta,pi-gamma)>
                        ! Q2m
                        do mu =-2,2
                            Q2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**mu*Q2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:)
                            pQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**mu*pQ2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:)
                            cQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**mu*cQ2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:)
                            pcQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**mu*pcQ2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:)
                        end do 
                        ! r2
                        r2_PNP_AMParray(ialpha,ibeta,igamma,:) = r2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:) 
                        pr2_PNP_AMParray(ialpha,ibeta,igamma,:) = pr2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:) 
                        ! r2_2b
                        r2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:) = r2_2b_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:,:)
                        pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:) = pr2_2b_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:,:)
                        ! 
                        Eccentri_PNP_AMParray(ialpha,ibeta,igamma,:) = Eccentri_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:) 
                        pEccentri_PNP_AMParray(ialpha,ibeta,igamma,:) = pEccentri_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:) 
                        ! store density matrix elements
                        call store_mix_density_matrix_elements(ialpha,ibeta,igamma)
                        cycle
                    end if
                    ! Axial+Parity
                    if(Proj_option%AMPtype==1 .and. Proj_option%Euler_Symmetry==1 .and. ibeta > (nbeta+1)/2) then 
                        if (MPI_Infor%rank == 0) write(*,'(A)') '(beta) symmetry (Axially+Parity).'
                        ! because <O R(0,beta,0)> = <O R(0,pi-beta,0)P> 
                        ! because <O R(0,beta,0)P> = <O R(0,pi-beta,0)>
                        ! In fact, in this case  ialpha = igamma = 1 
                        Norm_PNP_AMParray(ialpha,ibeta,igamma)  = pNorm_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        pNorm_PNP_AMParray(ialpha,ibeta,igamma) = Norm_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        Etot_PNP_AMParray(ialpha,ibeta,igamma) = pEtot_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma) 
                        pEtot_PNP_AMParray(ialpha,ibeta,igamma) = Etot_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        Particle_PNP_AMParray(ialpha,ibeta,igamma,1) = pParticle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,1)
                        Particle_PNP_AMParray(ialpha,ibeta,igamma,2) = pParticle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,2)
                        pParticle_PNP_AMParray(ialpha,ibeta,igamma,1) = Particle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,1)
                        pParticle_PNP_AMParray(ialpha,ibeta,igamma,2) = Particle_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,2)
                        ! N^2
                        N2_PNP_AMParray(ialpha,ibeta,igamma,:) = pN2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:)
                        pN2_PNP_AMParray(ialpha,ibeta,igamma,:) = N2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:) 
                        ! J^2
                        J2_PNP_AMParray(ialpha,ibeta,igamma) = pJ2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)
                        pJ2_PNP_AMParray(ialpha,ibeta,igamma) = J2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma)

                        ! because <T_{lm} R(0,beta,0)> = (-1)^m*<T_{lm} R(0,pi-beta,0)P>
                        ! because <T_{lm} R(0,beta,0)P> = (-1)^m*<T_{lm} R(0,pi-beta,0)>
                        ! Q2m
                        do mu =-2,2
                            Q2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**mu*pQ2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:)
                            pQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**mu*Q2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:)
                            cQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) =  (-1)**mu*pcQ2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:) 
                            pcQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,:) = (-1)**mu*cQ2m_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,mu,:)
                        end do
                        ! r2
                        r2_PNP_AMParray(ialpha,ibeta,igamma,:) = pr2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:)
                        pr2_PNP_AMParray(ialpha,ibeta,igamma,:) = r2_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:)
                        ! r2_2b
                        r2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:) = pr2_2b_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:,:)
                        pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:) = r2_2b_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:,:)
                        !
                        Eccentri_PNP_AMParray(ialpha,ibeta,igamma,:) = pEccentri_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:) 
                        pEccentri_PNP_AMParray(ialpha,ibeta,igamma,:) = Eccentri_PNP_AMParray(ialpha,nbeta+1-ibeta,ngamma+1-igamma,:) 
                        ! store density matrix elements 
                        call store_mix_density_matrix_elements(ialpha,ibeta,igamma)
                        cycle 
                    end if 

                    !############################################################
                    !#   calculate the overlap at given Euler angles
                    !#############################################################
                    if (MPI_Infor%rank == 0) write(*,'(A)') 'calculate_overlaps_after_PNP_at_Euler_angles ...'
                    call calculate_overlaps_after_PNP_at_Euler_angles(alpha,beta,gamma,         &
                        Norm_PNP_AMParray(ialpha,ibeta,igamma),         & ! <q_1|   R(alpha,beta,gamma)  |q_2 >
                        pNorm_PNP_AMParray(ialpha,ibeta,igamma),        & ! <q_1|   R(alpha,beta,gamma) P|q_2 >
                        Etot_PNP_AMParray(ialpha,ibeta,igamma),         & ! <q_1| H R(alpha,beta,gamma)  |q_2 >
                        pEtot_PNP_AMParray(ialpha,ibeta,igamma),        & ! <q_1| H R(alpha,beta,gamma) P|q_2 >
                        Particle_PNP_AMParray(ialpha,ibeta,igamma,:),   & ! <q_1| N R(alpha,beta,gamma)  |q_2 > for neutron  and proton
                        pParticle_PNP_AMParray(ialpha,ibeta,igamma,:),  & ! <q_1| N R(alpha,beta,gamma) P|q_2 > for neutron and proton
                        N2_PNP_AMParray(ialpha,ibeta,igamma,:),           & ! <q_1| N^2 R(alpha,beta,gamma)    |q_2 >
                        pN2_PNP_AMParray(ialpha,ibeta,igamma,:),          & ! <q_1| N^2 R(alpha,beta,gamma)  P |q_2 >
                        J2_PNP_AMParray(ialpha,ibeta,igamma),           & ! <q_1| J^2 R(alpha,beta,gamma)    |q_2 >
                        pJ2_PNP_AMParray(ialpha,ibeta,igamma),          & ! <q_1| J^2 R(alpha,beta,gamma)  P |q_2 >
                        Q2m_PNP_AMParray(ialpha,ibeta,igamma,:,:),      & ! <q_1| Q_{2 mu} R(alpha,beta,gamma)  |q_2 >
                        pQ2m_PNP_AMParray(ialpha,ibeta,igamma,:,:),     & ! <q_1| Q_{2 mu} R(alpha,beta,gamma) P|q_2 >
                        cQ2m_PNP_AMParray(ialpha,ibeta,igamma,:,:),     & ! <q_1| Q^{\dagger}_{2 mu} R(alpha,beta,gamma)   |q_2 >
                        pcQ2m_PNP_AMParray(ialpha,ibeta,igamma,:,:),    & ! <q_1| Q^{\dagger}_{2 mu} R(alpha,beta,gamma) P |q_2 >
                        r2_PNP_AMParray(ialpha,ibeta,igamma,:),         & ! <q_1| r^2 R(alpha,beta,gamma)    |q_2 >
                        pr2_PNP_AMParray(ialpha,ibeta,igamma,:),        & ! <q_1| r^2 R(alpha,beta,gamma)  P |q_2 >
                        r2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:),      & ! <q_1| r^2_2b R(alpha,beta,gamma)  |q_2 >
                        pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,:,:),    &
                        Eccentri_PNP_AMParray(ialpha,ibeta,igamma,:),   & ! <q_1| E_n R(alpha,beta,gamma)    |q_2 >
                        pEccentri_PNP_AMParray(ialpha,ibeta,igamma,:))    ! <q_1| E_n R(alpha,beta,gamma)  P |q_2 >            
                    ! store density matrix elements
                    call store_mix_density_matrix_elements(ialpha,ibeta,igamma)
                end do
            end do 
        end do
    end subroutine

    subroutine reset_parityOverlap

        pNorm_PNP_AMParray(:,:,:) = Norm_PNP_AMParray(:,:,:)
        pEtot_PNP_AMParray(:,:,:) = Etot_PNP_AMParray(:,:,:)
        pParticle_PNP_AMParray(:,:,:,:) = Particle_PNP_AMParray(:,:,:,:)
        pN2_PNP_AMParray(:,:,:,:) = N2_PNP_AMParray(:,:,:,:)
        pJ2_PNP_AMParray(:,:,:) = J2_PNP_AMParray(:,:,:)
        pQ2m_PNP_AMParray(:,:,:,:,:) = Q2m_PNP_AMParray(:,:,:,:,:)
        pcQ2m_PNP_AMParray(:,:,:,:,:) = cQ2m_PNP_AMParray(:,:,:,:,:)
        pr2_PNP_AMParray(:,:,:,:) = r2_PNP_AMParray(:,:,:,:)
        pEccentri_PNP_AMParray(:,:,:,:) = Eccentri_PNP_AMParray(:,:,:,:)
        pr2_2b_PNP_AMParray(:,:,:,:,:) = r2_2b_PNP_AMParray(:,:,:,:,:)
    end subroutine

    subroutine calcualate_Norm_Hamiltonian_ParticleNumber_kernels
        !-------------------------------------------------------------------------------------------------------------------------------
        !  
        !       <J K_1 q_1 Pi | O | J K_2 q_2 Pi >
        !   =   <q_1| O P^{J}_{K_1 K_2} P^{Pi} |q_2 >
        !   =  (2J_i +1)/(8*pi^2) \int d alpha d beta d gamma  D^{J*}_{K_1 K_2}(alpha,beta,gamma) 
        !                          * <q_1| O R(alpha,beta,gamma) P^{Pi} |q_2 >
        !  where O = 1 or N or H
        !  Note:
        !  1) fac2 = 1.0d0 + CDEXP(-K1*cpi) + CDEXP(-K2*cpi) + CDEXP(-K1*cpi-K2*cpi)
        !     This factor allows alpha and gamma to be reduced to [0, pi], but D2 symmetry is required.
        !-------------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: gcm_space,projection_mesh,kernels,Proj_option
        use Basis, only: djmk
        integer :: ialpha,ibeta,igamma,J,K1_start,K1_end,K2_start,K2_end,K1,K2,it
        real(r64) :: alpha, beta, gamma, w
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac,Norm_PNP_Euler,pNorm_PNP_Euler,Etot_PNP_Euler,pEtot_PNP_Euler,Particle_PNP_Euler(2),pParticle_PNP_Euler(2)
        kernels%N_KK = (0.d0,0.d0)
        kernels%H_KK = (0.d0,0.d0)
        kernels%X_KK = (0.d0,0.d0)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            ! H, N kernels
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
            do ialpha = 1, projection_mesh%nalpha
                alpha = projection_mesh%alpha(ialpha)
                calpha = DCMPLX(0.d0,alpha)
                do ibeta = 1, projection_mesh%nbeta
                    beta = projection_mesh%beta(ibeta)
                    do igamma = 1, projection_mesh%ngamma
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        do K1 = K1_start, K1_end
                            do K2 = K2_start, K2_end
                                if(Proj_option%AMPtype==0) then
                                    fac = 1.d0
                                else if (Proj_option%AMPtype==1) then
                                    w = projection_mesh%wbeta(ibeta)
                                    fac1 = (2*J+1)/(2.0d0)*dsin(beta)*djmk(J,K1,K2,dcos(beta),0)
                                    fac = fac1*w
                                else
                                    cpi = DCMPLX(0.d0,pi) ! i*pi
                                    w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                    fac1 = (2*J+1)/(8.0d0*pi**2)*dsin(beta)*djmk(J,K1,K2,dcos(beta),0)*CDEXP(-K1*calpha-K2*cgamma)
                                    fac2 = 1.0d0 + CDEXP(-K1*cpi) + CDEXP(-K2*cpi) + CDEXP(-K1*cpi-K2*cpi) ! D2 symmetry is required, with alpha, gamma in [0, pi].
                                    fac = fac1*fac2*w
                                end if
                                Norm_PNP_Euler = Norm_PNP_AMParray(ialpha,ibeta,igamma)
                                pNorm_PNP_Euler = pNorm_PNP_AMParray(ialpha,ibeta,igamma)
                                Etot_PNP_Euler = Etot_PNP_AMParray(ialpha,ibeta,igamma)
                                pEtot_PNP_Euler = pEtot_PNP_AMParray(ialpha,ibeta,igamma)
                                Particle_PNP_Euler(1) = Particle_PNP_AMParray(ialpha,ibeta,igamma,1)
                                Particle_PNP_Euler(2) = Particle_PNP_AMParray(ialpha,ibeta,igamma,2)
                                pParticle_PNP_Euler(1) = pParticle_PNP_AMParray(ialpha,ibeta,igamma,1)
                                pParticle_PNP_Euler(2) = pParticle_PNP_AMParray(ialpha,ibeta,igamma,2)
                                !<J K_1 q_1 Pi | H | J K_2 q_2 Pi >
                                ! Pi = +
                                kernels%N_KK(J,K1,K2,1) = kernels%N_KK(J,K1,K2,1) + fac*(Norm_PNP_Euler + pNorm_PNP_Euler)/2.0d0
                                ! Pi = -
                                kernels%N_KK(J,K1,K2,2) = kernels%N_KK(J,K1,K2,2) + fac*(Norm_PNP_Euler - pNorm_PNP_Euler)/2.0d0

                                !<J K_1 q_1 Pi | J K_2 q_2 Pi >
                                ! Pi = +
                                kernels%H_KK(J,K1,K2,1) = kernels%H_KK(J,K1,K2,1) + fac*(Etot_PNP_Euler + pEtot_PNP_Euler)/2.0d0
                                ! Pi = -
                                kernels%H_KK(J,K1,K2,2) = kernels%H_KK(J,K1,K2,2) + fac*(Etot_PNP_Euler - pEtot_PNP_Euler)/2.0d0

                                do it =1,2
                                     !<J K_1 q_1 Pi | N | J K_2 q_2 Pi >
                                    ! Pi = +
                                    kernels%X_KK(J,K1,K2,it,1) = kernels%X_KK(J,K1,K2,it,1) + fac*(Particle_PNP_Euler(it) + pParticle_PNP_Euler(it))/2.0d0
                                    ! Pi = -
                                    kernels%X_KK(J,K1,K2,it,2) = kernels%X_KK(J,K1,K2,it,2) + fac*(Particle_PNP_Euler(it) - pParticle_PNP_Euler(it))/2.0d0
                                end do

                            end do 
                        end do 
                    end do 
                end do 
            end do
        end do
    end subroutine

    subroutine calculate_N2_kernels
        !-------------------------------------------------------------------------------------------------------------------------------
        !  
        !      <J_f K_f q_1 Pi_i| N^2 |J_i K_i q_2 Pi_i> 
        !   =  <J   K_f q_1 Pi  | N^2 |J   K_i q_2 Pi> 
        !   =  <q_1| J^2 P^{J}_{K_f K_i} P^{Pi} |q_2 >
        !   =  (2J_i +1)/(pi^2) \int d alpha d beta d gamma  D^{J*}_{K_f K_i}(alpha,beta,gamma) 
        !                          * <q_1| N^2 R(alpha,beta,gamma) P^{Pi} |q_2 >
        !-------------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: gcm_space,projection_mesh,kernels,Proj_option
        use Basis, only: djmk
        integer :: ialpha,ibeta,igamma,J,Ji,Jf,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,it
        real(r64) :: alpha, beta, gamma, w
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac
        kernels%N2_KK = (0.d0,0.d0)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            Ji = J
            Jf = J
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
            do ialpha = 1, projection_mesh%nalpha
                alpha = projection_mesh%alpha(ialpha)
                calpha = DCMPLX(0.d0,alpha)
                do ibeta = 1, projection_mesh%nbeta
                    beta = projection_mesh%beta(ibeta)
                    do igamma = 1, projection_mesh%ngamma
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        do Kf = Kf_start, Kf_end
                            do Ki = Ki_start, Ki_end            
                                if(Proj_option%AMPtype==0) then
                                    fac = 1
                                else if (Proj_option%AMPtype==1) then
                                    w = projection_mesh%wbeta(ibeta)
                                    fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)
                                    fac = fac1*w
                                else
                                    cpi = DCMPLX(0.d0,pi) ! i*pi
                                    w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                    fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)*CDEXP(-Kf*calpha-Ki*cgamma)
                                    fac2 = 1.0d0 + CDEXP(-Kf*cpi) + CDEXP(-Ki*cpi) + CDEXP(-Kf*cpi-Ki*cpi)
                                    fac = fac1*fac2*w
                                end if

                                do it =1,2
                                   ! <Ji Kf q1 Pi|N^2|Ji Ki q2 Pi> 
                                   ! Pi = +
                                   kernels%N2_KK(Ji,Kf,Ki,it,1) = kernels%N2_KK(Ji,Kf,Ki,it,1) + fac*(N2_PNP_AMParray(ialpha,ibeta,igamma,it) + pN2_PNP_AMParray(ialpha,ibeta,igamma,it))/2.d0
                                   ! Pi = -
                                   kernels%N2_KK(Ji,Kf,Ki,it,2) = kernels%N2_KK(Ji,Kf,Ki,it,2) + fac*(N2_PNP_AMParray(ialpha,ibeta,igamma,it) - pN2_PNP_AMParray(ialpha,ibeta,igamma,it))/2.d0
                               end do

                            end do 
                        end do 
                    end do 
                end do 
            end do
        end do
    end subroutine

    subroutine calculate_J2_kernels
        !-------------------------------------------------------------------------------------------------------------------------------
        !  
        !      <J_f K_f q_1 Pi_i| J^2 |J_i K_i q_2 Pi_i> 
        !   =  <J   K_f q_1 Pi  | J^2 |J   K_i q_2 Pi> 
        !   =  <q_1| J^2 P^{J}_{K_f K_i} P^{Pi} |q_2 >
        !   =  (2J_i +1)/(pi^2) \int d alpha d beta d gamma  D^{J*}_{K_f K_i}(alpha,beta,gamma) 
        !                          * <q_1| J^2 R(alpha,beta,gamma) P^{Pi} |q_2 >
        !-------------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: gcm_space,projection_mesh,kernels,Proj_option
        use Basis, only: djmk
        integer :: ialpha,ibeta,igamma,J,Ji,Jf,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,it
        real(r64) :: alpha, beta, gamma, w
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac
        kernels%J2_KK = (0.d0,0.d0)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            Ji = J
            Jf = J
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
            do ialpha = 1, projection_mesh%nalpha
                alpha = projection_mesh%alpha(ialpha)
                calpha = DCMPLX(0.d0,alpha)
                do ibeta = 1, projection_mesh%nbeta
                    beta = projection_mesh%beta(ibeta)
                    do igamma = 1, projection_mesh%ngamma
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        do Kf = Kf_start, Kf_end
                            do Ki = Ki_start, Ki_end            
                                if(Proj_option%AMPtype==0) then
                                    fac = 1
                                else if (Proj_option%AMPtype==1) then
                                    w = projection_mesh%wbeta(ibeta)
                                    fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)
                                    fac = fac1*w
                                else
                                    cpi = DCMPLX(0.d0,pi) ! i*pi
                                    w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                    fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)*CDEXP(-Kf*calpha-Ki*cgamma)
                                    fac2 = 1.0d0 + CDEXP(-Kf*cpi) + CDEXP(-Ki*cpi) + CDEXP(-Kf*cpi-Ki*cpi)
                                    fac = fac1*fac2*w
                                end if

                                ! <Jf Kf q1 Pi|J^2|Ji Ki q2 Pi> 
                                ! Pi = +
                                kernels%J2_KK(Ji,Kf,Ki,1) = kernels%J2_KK(Ji,Kf,Ki,1) + fac*(J2_PNP_AMParray(ialpha,ibeta,igamma) + pJ2_PNP_AMParray(ialpha,ibeta,igamma))/2.d0
                                ! Pi = -
                                kernels%J2_KK(Ji,Kf,Ki,2) = kernels%J2_KK(Ji,Kf,Ki,2) + fac*(J2_PNP_AMParray(ialpha,ibeta,igamma) - pJ2_PNP_AMParray(ialpha,ibeta,igamma))/2.d0
                            end do 
                        end do 
                    end do 
                end do 
            end do
        end do
    end subroutine

    subroutine calculate_EM_kernels
        !-------------------------------------------------------------------------------------------------------------------------------
        !  
        !       <J_f K_f q_1 Pi_f ||T_lambda|| J_i K_i q_2 Pi_i >
        !   = fac_parity * sqrt(2J_f + 1) * \sum_{K_1 mu} C^{J_f K_f}_{J_i K_1 lambda mu}
        !                           *<q_1| T_{lambda mu} P^{J_i}_{K_1 K_i} P^{Pi_i} |q_2 >
        !  where fac_parity =  (1 + Pi_i*Pi_f*(-1)^lambda )/2 while T is electric multipole operator
        !        fac_parity =  (1 - Pi_i*Pi_f*(-1)^lambda )/2 while T is magnetic multipole operator  
        ! 
        !       <q_1| T_{lambda mu} P^{J_i}_{K_1 K_i} P^{Pi_i} |q_2 >
        !     =  (2J_i +1)/(8*pi^2) \int d alpha d beta d gamma  D^{J_i*}_{K_1 K_i}(alpha,beta,gamma) 
        !                          * <q_1| T_{lambda mu} R(alpha,beta,gamma) P^{Pi_i} |q_2 >
        !                       --------------------------------
        !
        !  Relation between Clebsch–Gordan coefficients and Wigner 3j symbols:
        !    C^{J_f K_f}_{J_i K_1 lambda mu} = (-1)^{J_f - K_f} * sqrt(2J_f + 1) * wigner3j(J_f lambda J_i -K_f mu K_1)
        !  Here, J_i,J_f are integers.
        !  =================================================================================================================
        !   q_1 , q_2 exchange :
        ! 
        !       <J_f K_f q_2 Pi_f ||T_lambda|| J_i K_i q_1 Pi_i >
        !   = fac_parity * sqrt(2J_f + 1) * \sum_{K_1 mu} C^{J_f K_f}_{J_i K_1 lambda mu}
        !                           *<q_2| T_{lambda mu} P^{J_i}_{K_1 K_i} P^{Pi_i} |q_1 >
        !
        !
       !       <q_2| T_{lambda mu} P^{J_i}_{K_1 K_i} P^{Pi_i} |q_1 >
        !     =  (2J_i +1)/(8*pi^2) \int d alpha d beta d gamma  D^{J_i*}_{K_1 K_i}(alpha,beta,gamma) 
        !                          * <q_2| T_{lambda mu} R(alpha,beta,gamma) P^{Pi_i} |q_1 >
        !     = (2J_i +1)/(8*pi^2) \int d alpha d beta d gamma  D^{J_i}_{K_i K_1}(alpha,beta,gamma) 
        !                           *\sum_{nu} (-1)^{mu-nu}  D^{lambda*}_{-nu-mu}(alpha,beta,gamma)  
        !                           *(<q_1| M^{\dagger}_{lambda mu} R(alpha,beta,gamma) P^{Pi_f} |q_2 >)^*
        !
        !  ===================================================================================================================
        !  
        !  Note:
        !    1) fac2 = 1.0d0 + (-1)**mu*CDEXP(-K1*cpi) + CDEXP(-Ki*cpi) + (-1)**mu*CDEXP(-K1*cpi-Ki*cpi)
        !       This factor allows alpha and gamma to be reduced to [0, pi], but D2 symmetry is required.
        !-------------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: gcm_space,projection_mesh,kernels,Proj_option
        use Basis, only: djmk
        use EM, only: wigner3j
        integer :: ialpha,ibeta,igamma,J,Ji,Jf,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,mu,K1,it,nu
        real(r64) :: alpha, beta, gamma, w, Pi_i, Pi_f
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac
        kernels%Q2_KK_12 = (0.d0,0.d0)
        kernels%Q2_KK_21 = (0.d0,0.d0)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            ! Q2 Kernels
            Ji = J
            Jf = J + 2
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
            do ialpha = 1, projection_mesh%nalpha
                alpha = projection_mesh%alpha(ialpha)
                calpha = DCMPLX(0.d0,alpha)
                do ibeta = 1, projection_mesh%nbeta
                    beta = projection_mesh%beta(ibeta)
                    do igamma = 1, projection_mesh%ngamma
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        do Kf = Kf_start, Kf_end
                            do Ki = Ki_start, Ki_end            
                                do mu = -2, 2 ! lambda = 2
                                    K1 = Kf - mu
                                    if(Proj_option%AMPtype==0) then
                                        fac = 1.d0
                                    else if (Proj_option%AMPtype==1) then
                                        w = projection_mesh%wbeta(ibeta)
                                        fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,K1,Ki,dcos(beta),0) ! Note: djmk(Ji,K1,Ki,dcos(beta),0) = d^{Ji}_{K1 Ki}(beta)
                                        fac = fac1*w
                                    else
                                        cpi = DCMPLX(0.d0,pi) ! cpi = i*pi
                                        w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                        fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,K1,Ki,dcos(beta),0)*CDEXP(-K1*calpha-Ki*cgamma) ! Note: CDEXP(-Ki*calpha-K1*cgamma) = e^(-i*K1*alpha-i*Ki*gamma)
                                        fac2 = 1.0d0 + (-1)**mu*CDEXP(-K1*cpi) + CDEXP(-Ki*cpi) + (-1)**mu*CDEXP(-K1*cpi-Ki*cpi) ! D2 symmetry is required, with alpha, gamma in [0, pi].
                                        fac = fac1*fac2*w
                                    end if
                                    do it =1,2

                                        ! The kernel is non-zero when Pi_i * Pi_f * (-1)^lambda = 1.
                                        ! when lambda = 2, Pi_i is same as Pi_f (so Pi_f is not explicitly labeled here).

                                        ! <Ji+2 Kf q1 Pi_f||Q2||Ji Ki q2 Pi_i> 
                                        ! Pi_i = + 
                                        Pi_i =  1.d0
                                        kernels%Q2_KK_12(Ji,Kf,Ki,it,1) = kernels%Q2_KK_12(Ji,Kf,Ki,it,1) + fac*(2*Jf+1)*(-1)**(Jf-Kf)*wigner3j(Jf,2,Ji,-Kf,mu,K1,IS=0)*&
                                            (Q2m_PNP_AMParray(ialpha,ibeta,igamma,mu,it) + Pi_i*pQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,it))/2.d0
                                        ! Pi_i = - 
                                        Pi_i =  -1.d0
                                        kernels%Q2_KK_12(Ji,Kf,Ki,it,2) = kernels%Q2_KK_12(Ji,Kf,Ki,it,2) + fac*(2*Jf+1)*(-1)**(Jf-Kf)*wigner3j(Jf,2,Ji,-Kf,mu,K1,IS=0)*&
                                            (Q2m_PNP_AMParray(ialpha,ibeta,igamma,mu,it) + Pi_i*pQ2m_PNP_AMParray(ialpha,ibeta,igamma,mu,it))/2.d0
                                    end do 
                                end do
                                ! q_i, q_f exchange
                                do mu = -2, 2 ! lambda = 2
                                    K1 = Kf - mu
                                    if(Proj_option%AMPtype==0) then
                                        fac = 1
                                    else if (Proj_option%AMPtype==1) then
                                        w = projection_mesh%wbeta(ibeta)
                                        fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,Ki,K1,dcos(beta),0) ! Note: djmk(Ji,Ki,K1,dcos(beta),0) = d^{Ji}_{Ki K1}(beta)
                                        fac = fac1*w
                                    else
                                        cpi = DCMPLX(0.d0,pi) ! cpi = i*pi
                                        w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                        fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,Ki,K1,dcos(beta),0)*CDEXP(Ki*calpha+K1*cgamma) ! Note: CDEXP(Ki*calpha+K1*cgamma) =  e^(i*Ki*alpha+i*K1*gamma)
                                        fac2 = 1.0d0 + (-1)**mu*CDEXP(-K1*cpi) + CDEXP(-Ki*cpi) + (-1)**mu*CDEXP(-K1*cpi-Ki*cpi)
                                        fac = fac1*fac2*w
                                    end if
                                    do nu = -2, 2 ! lambda = 2
                                        do it =1,2
                                            ! <Ji+2 Kf q2 Pi_f ||Q2||Ji Ki q1 Pi_i> 
                                            ! Pi_i = +
                                            Pi_f =  1.d0 ! because Pi_i*Pi_f*(-1)^2=1, so Pi_f = + 
                                            kernels%Q2_KK_21(Ji,Kf,Ki,it,1) = kernels%Q2_KK_21(Ji,Kf,Ki,it,1) + fac*(2*Jf+1)*(-1)**(Jf-Kf)*wigner3j(Jf,2,Ji,-Kf,mu,K1,IS=0)*&
                                                (-1)**(mu-nu)*djmk(2,-nu,-mu,dcos(beta),0)*CDEXP(nu*calpha+mu*cgamma)*&
                                                DCONJG(cQ2m_PNP_AMParray(ialpha,ibeta,igamma,nu,it) + Pi_f*pcQ2m_PNP_AMParray(ialpha,ibeta,igamma,nu,it))/2.d0
                                            ! Pi_i = -
                                            Pi_f =  -1.d0 ! because Pi_i*Pi_f*(-1)^2=1, so Pi_f = - 
                                            kernels%Q2_KK_21(Ji,Kf,Ki,it,2) = kernels%Q2_KK_21(Ji,Kf,Ki,it,2) + fac*(2*Jf+1)*(-1)**(Jf-Kf)*wigner3j(Jf,2,Ji,-Kf,mu,K1,IS=0)*&
                                                (-1)**(mu-nu)*djmk(2,-nu,-mu,dcos(beta),0)*CDEXP(nu*calpha+mu*cgamma)*&
                                                DCONJG(cQ2m_PNP_AMParray(ialpha,ibeta,igamma,nu,it) + Pi_f*pcQ2m_PNP_AMParray(ialpha,ibeta,igamma,nu,it))/2.d0
                                        end do 
                                    end do 
                                end do
                            end do 
                        end do 
                    end do 
                end do 
            end do
        end do
    end subroutine

    subroutine calculate_E0_kernel
        !-------------------------------------------------------------------------------------------------------------------------------
        !  
        !      <J_f K_f q_1 Pi_i| r2 |J_i K_i q_2 Pi_i> 
        !   =  <J   K_f q_1 Pi  | r2 |J   K_i q_2 Pi> 
        !   =  <q_1| r2 P^{J}_{K_f K_i} P^{Pi} |q_2 >
        !   =  (2J_i +1)/(pi^2) \int d alpha d beta d gamma  D^{J*}_{K_f K_i}(alpha,beta,gamma) 
        !                          * <q_1| r2 R(alpha,beta,gamma) P^{Pi} |q_2 >
        !-------------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: gcm_space,projection_mesh,kernels,Proj_option
        use Basis, only: djmk
        integer :: ialpha,ibeta,igamma,J,Ji,Jf,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,it
        real(r64) :: alpha, beta, gamma, w
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac
        kernels%E0_KK = (0.d0,0.d0)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            Ji = J
            Jf = J
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
            do ialpha = 1, projection_mesh%nalpha
                alpha = projection_mesh%alpha(ialpha)
                calpha = DCMPLX(0.d0,alpha)
                do ibeta = 1, projection_mesh%nbeta
                    beta = projection_mesh%beta(ibeta)
                    do igamma = 1, projection_mesh%ngamma
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        do Kf = Kf_start, Kf_end
                            do Ki = Ki_start, Ki_end            
                                if(Proj_option%AMPtype==0) then
                                    fac = 1
                                else if (Proj_option%AMPtype==1) then
                                    w = projection_mesh%wbeta(ibeta)
                                    fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)
                                    fac = fac1*w
                                else
                                    cpi = DCMPLX(0.d0,pi) ! i*pi
                                    w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                    fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)*CDEXP(-Kf*calpha-Ki*cgamma)
                                    fac2 = 1.0d0 + CDEXP(-Kf*cpi) + CDEXP(-Ki*cpi) + CDEXP(-Kf*cpi-Ki*cpi)
                                    fac = fac1*fac2*w
                                end if
                                do it =1,2
                                    ! <Jf Kf q1 Pi|r2|Ji Ki q2 Pi> 
                                    ! Pi = +
                                    kernels%E0_KK(Ji,Kf,Ki,it,1) = kernels%E0_KK(Ji,Kf,Ki,it,1) + fac*(r2_PNP_AMParray(ialpha,ibeta,igamma,it) + pr2_PNP_AMParray(ialpha,ibeta,igamma,it))/2.d0
                                    ! Pi = -
                                    kernels%E0_KK(Ji,Kf,Ki,it,2) = kernels%E0_KK(Ji,Kf,Ki,it,2) + fac*(r2_PNP_AMParray(ialpha,ibeta,igamma,it) - pr2_PNP_AMParray(ialpha,ibeta,igamma,it))/2.d0
                                end do 
                            end do 
                        end do 
                    end do 
                end do 
            end do
        end do
    end subroutine

    subroutine calculate_Eccentricity_kernel
        !-------------------------------------------------------------------------------------------------------------------------------
        !  
        !      <J_f K_f q_1 Pi_i| E_n |J_i K_i q_2 Pi_i> 
        !   =  <J   K_f q_1 Pi  | E_n |J   K_i q_2 Pi> 
        !   =  <q_1| r2 P^{J}_{K_f K_i} P^{Pi} |q_2 >
        !   =  (2J_i +1)/(pi^2) \int d alpha d beta d gamma  D^{J*}_{K_f K_i}(alpha,beta,gamma) 
        !                          * <q_1| E_n R(alpha,beta,gamma) P^{Pi} |q_2 >
        !-------------------------------------------------------------------------------------------------------------------------------------
        use Constants, only: pi
        use Globals, only: gcm_space,projection_mesh,kernels,Proj_option
        use Basis, only: djmk
        integer :: ialpha,ibeta,igamma,J,Ji,Jf,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,it
        real(r64) :: alpha, beta, gamma, w
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac
        kernels%Eccentricity_KK = (0.d0,0.d0)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            Ji = J
            Jf = J
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
            do ialpha = 1, projection_mesh%nalpha
                alpha = projection_mesh%alpha(ialpha)
                calpha = DCMPLX(0.d0,alpha)
                do ibeta = 1, projection_mesh%nbeta
                    beta = projection_mesh%beta(ibeta)
                    do igamma = 1, projection_mesh%ngamma
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        do Kf = Kf_start, Kf_end
                            do Ki = Ki_start, Ki_end            
                                if(Proj_option%AMPtype==0) then
                                    fac = 1
                                else if (Proj_option%AMPtype==1) then
                                    w = projection_mesh%wbeta(ibeta)
                                    fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)
                                    fac = fac1*w
                                else
                                    cpi = DCMPLX(0.d0,pi) ! i*pi
                                    w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                    fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)*CDEXP(-Kf*calpha-Ki*cgamma)
                                    fac2 = 1.0d0 + CDEXP(-Kf*cpi) + CDEXP(-Ki*cpi) + CDEXP(-Kf*cpi-Ki*cpi)
                                    fac = fac1*fac2*w
                                end if
                                ! <Jf Kf q1 Pi|E_n|Ji Ki q2 Pi> 
                                ! Pi = +
                                kernels%Eccentricity_KK(Ji,Kf,Ki,1,1) = kernels%Eccentricity_KK(Ji,Kf,Ki,1,1) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,1) + pEccentri_PNP_AMParray(ialpha,ibeta,igamma,1))/2.d0 ! 1B
                                kernels%Eccentricity_KK(Ji,Kf,Ki,1,2) = kernels%Eccentricity_KK(Ji,Kf,Ki,1,2) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,2) + pEccentri_PNP_AMParray(ialpha,ibeta,igamma,2))/2.d0 ! 2B
                                kernels%Eccentricity_KK(Ji,Kf,Ki,1,3) = kernels%Eccentricity_KK(Ji,Kf,Ki,1,3) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,3) + pEccentri_PNP_AMParray(ialpha,ibeta,igamma,3))/2.d0 ! direct term
                                kernels%Eccentricity_KK(Ji,Kf,Ki,1,4) = kernels%Eccentricity_KK(Ji,Kf,Ki,1,4) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,4) + pEccentri_PNP_AMParray(ialpha,ibeta,igamma,4))/2.d0 ! exchange term
                                kernels%Eccentricity_KK(Ji,Kf,Ki,1,5) = kernels%Eccentricity_KK(Ji,Kf,Ki,1,5) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,5) + pEccentri_PNP_AMParray(ialpha,ibeta,igamma,5))/2.d0 ! kappa term
                                ! Pi = -
                                kernels%Eccentricity_KK(Ji,Kf,Ki,2,1) = kernels%Eccentricity_KK(Ji,Kf,Ki,2,1) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,1) - pEccentri_PNP_AMParray(ialpha,ibeta,igamma,1))/2.d0 ! 1B
                                kernels%Eccentricity_KK(Ji,Kf,Ki,2,2) = kernels%Eccentricity_KK(Ji,Kf,Ki,2,2) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,2) - pEccentri_PNP_AMParray(ialpha,ibeta,igamma,2))/2.d0 ! 1B 
                                kernels%Eccentricity_KK(Ji,Kf,Ki,2,3) = kernels%Eccentricity_KK(Ji,Kf,Ki,2,3) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,3) - pEccentri_PNP_AMParray(ialpha,ibeta,igamma,3))/2.d0 ! direct term
                                kernels%Eccentricity_KK(Ji,Kf,Ki,2,4) = kernels%Eccentricity_KK(Ji,Kf,Ki,2,4) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,4) - pEccentri_PNP_AMParray(ialpha,ibeta,igamma,4))/2.d0 ! exchange term
                                kernels%Eccentricity_KK(Ji,Kf,Ki,2,5) = kernels%Eccentricity_KK(Ji,Kf,Ki,2,5) + fac*(Eccentri_PNP_AMParray(ialpha,ibeta,igamma,5) - pEccentri_PNP_AMParray(ialpha,ibeta,igamma,5))/2.d0 ! kappa term
                            end do 
                        end do 
                    end do 
                end do 
            end do
        end do
    end subroutine

    subroutine calculate_r2_2body_kernel
        ! 计算 <J Kf q1 Pi| r^2_2B |J Ki q2 Pi>
        use Constants, only: pi
        use Globals, only: gcm_space,projection_mesh,kernels,Proj_option
        use Basis, only: djmk
        integer :: ialpha,ibeta,igamma,J,Ji,Jf,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,it,it_chan
        real(r64) :: alpha, beta, gamma, w
        complex(r64) :: calpha,cgamma,cpi,fac1,fac2,fac

        kernels%r2_2b_KK = (0.d0, 0.d0)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            Ji = J
            Jf = J
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
            do ialpha = 1, projection_mesh%nalpha
                alpha = projection_mesh%alpha(ialpha)
                calpha = DCMPLX(0.d0,alpha)
                do ibeta = 1, projection_mesh%nbeta
                    beta = projection_mesh%beta(ibeta)
                    do igamma = 1, projection_mesh%ngamma
                        gamma = projection_mesh%gamma(igamma)
                        cgamma = DCMPLX(0.d0,gamma)
                        do Kf = Kf_start, Kf_end
                            do Ki = Ki_start, Ki_end
                                if(Proj_option%AMPtype==0) then
                                    fac = 1
                                else if (Proj_option%AMPtype==1) then
                                    w = projection_mesh%wbeta(ibeta)
                                    fac1 = (2*Ji+1)/(2.0d0)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)
                                    fac = fac1*w
                                else
                                    cpi = DCMPLX(0.d0,pi)
                                    w = projection_mesh%walpha(ialpha)*projection_mesh%wbeta(ibeta)*projection_mesh%wgamma(igamma)
                                    fac1 = (2*Ji+1)/(8.0d0*pi**2)*dsin(beta)*djmk(Ji,Kf,Ki,dcos(beta),0)*CDEXP(-Kf*calpha-Ki*cgamma)
                                    fac2 = 1.0d0 + CDEXP(-Kf*cpi) + CDEXP(-Ki*cpi) + CDEXP(-Kf*cpi-Ki*cpi)
                                    fac = fac1*fac2*w
                                end if
                                do it_chan = 1,3
                                    ! Pi = +
                                    kernels%r2_2b_KK(Ji,Kf,Ki,1,1,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,1,1,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,1,it_chan) + pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,1,it_chan))/2.d0 ! 1B
                                    kernels%r2_2b_KK(Ji,Kf,Ki,1,2,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,1,2,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,2,it_chan) + pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,2,it_chan))/2.d0 ! 2B
                                    kernels%r2_2b_KK(Ji,Kf,Ki,1,3,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,1,3,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,3,it_chan) + pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,3,it_chan))/2.d0 ! direct term
                                    kernels%r2_2b_KK(Ji,Kf,Ki,1,4,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,1,4,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,4,it_chan) + pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,4,it_chan))/2.d0 ! exchange term
                                    kernels%r2_2b_KK(Ji,Kf,Ki,1,5,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,1,5,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,5,it_chan) + pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,5,it_chan))/2.d0 ! kappa term
                                    ! Pi = -
                                    kernels%r2_2b_KK(Ji,Kf,Ki,2,1,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,2,1,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,1,it_chan) - pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,1,it_chan))/2.d0 ! 1B
                                    kernels%r2_2b_KK(Ji,Kf,Ki,2,2,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,2,2,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,2,it_chan) - pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,2,it_chan))/2.d0 ! 2B
                                    kernels%r2_2b_KK(Ji,Kf,Ki,2,3,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,2,3,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,3,it_chan) - pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,3,it_chan))/2.d0 ! direct term
                                    kernels%r2_2b_KK(Ji,Kf,Ki,2,4,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,2,4,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,4,it_chan) - pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,4,it_chan))/2.d0 ! exchange term
                                    kernels%r2_2b_KK(Ji,Kf,Ki,2,5,it_chan) = kernels%r2_2b_KK(Ji,Kf,Ki,2,5,it_chan) + fac*(r2_2b_PNP_AMParray(ialpha,ibeta,igamma,5,it_chan) - pr2_2b_PNP_AMParray(ialpha,ibeta,igamma,5,it_chan))/2.d0 ! kappa term
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine
    !------------------------------------------------
    !   PNP Integration
    !------------------------------------------------
    subroutine calculate_overlaps_after_PNP_at_Euler_angles(alpha, beta, gamma, Norm_PNP, pNorm_PNP, Etot_PNP, pEtot_PNP, &
        Particle_PNP, pParticle_PNP,N2_PNP,pN2_PNP,J2_PNP, pJ2_PNP, Q2m_PNP, pQ2m_PNP,cQ2m_PNP,pcQ2m_PNP,r2_PNP, pr2_PNP, &
        r2_2b_PNP, pr2_2b_PNP,Eccentri_PNP,pEccentri_PNP)

        use Globals, only: Proj_option
        use Mixed, only: calculate_mixed_DensCurrTens_and_norm_overlap
        real(r64), intent(in) :: alpha, beta, gamma
        complex(r64),intent(out) :: Norm_PNP, pNorm_PNP, Etot_PNP, pEtot_PNP, Particle_PNP(2), pParticle_PNP(2),N2_PNP(2),pN2_PNP(2),J2_PNP, pJ2_PNP, &
                                    r2_PNP(2), pr2_PNP(2),r2_2b_PNP(5,3), pr2_2b_PNP(5,3),Eccentri_PNP(5),pEccentri_PNP(5)
        complex(r64), dimension(-2:2,2),intent(out) :: Q2m_PNP, pQ2m_PNP,cQ2m_PNP,pcQ2m_PNP

        ! 1) calcualate mixed ... matrix elements
        ! 2) calcualate mixed ... in coordinate space 
        ! 3) calcualate norm overlap
        call calculate_mixed_DensCurrTens_and_norm_overlap(alpha, beta, gamma)

        ! Integration over PNP mesh points(phi_it) and spatial coordinates (r,theta,phi)  
        call calculate_norm_overlap_and_particle_number_after_PNP(Norm_PNP,pNorm_PNP,Particle_PNP,pParticle_PNP)
        call calculate_Rotated_Energy_after_PNP(Etot_PNP,pEtot_PNP)
        if(Proj_option%checkN2J2 == 1) then
            call calculate_N2_after_PNP(N2_PNP,pN2_PNP)
            call calculate_J2_after_PNP(J2_PNP, pJ2_PNP)
        end if 
        call calculate_Qlm_after_PNP(Q2m_PNP,pQ2m_PNP,cQ2m_PNP,pcQ2m_PNP)
        call calculate_r2_after_PNP(r2_PNP, pr2_PNP)
        ! 直接调用r² 2B 
        call calculate_r2_2body_after_PNP(r2_2b_PNP, pr2_2b_PNP)
        if( Proj_option%EccentriType == 1 .or. Proj_option%EccentriType == 3) then
            call calculate_Eccentricity_after_PNP(Eccentri_PNP,pEccentri_PNP)
        end if 
    end subroutine

    subroutine calculate_norm_overlap_and_particle_number_after_PNP(Norm_PNP,pNorm_PNP,Particle_PNP,pParticle_PNP)
        !---------------------------------------------------------------------
        !      Norm_PNP: <q_1|   R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >  
        !     pNorm_PNP: <q_1|   R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !  Particle_PNP: <q_1| N R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >  
        ! pParticle_PNP: <q_1| N R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !---------------------------------------------------------------------
        use Constants, only: ngr,ntheta,nphi
        use Globals, only: projection_mesh,nucleus_attributes,mix
        complex(r64),intent(out) :: Norm_PNP,pNorm_PNP,Particle_PNP(2),pParticle_PNP(2)
        integer :: L_n,L_p,phi_n_index, phi_p_index,i
        real(r64) :: phi_n,phi_p
        complex(r64) :: emiNphi,emiZphi,fac,pfac
        Norm_PNP = (0.0d0,0.0d0)
        pNorm_PNP = (0.0d0,0.0d0)
        Particle_PNP(1) = (0.0d0,0.0d0)
        Particle_PNP(2) = (0.0d0,0.0d0)
        pParticle_PNP(1) = (0.0d0,0.0d0)
        pParticle_PNP(2) = (0.0d0,0.0d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)
        do  phi_n_index = 1, L_n
            phi_n =  phi_n_index*projection_mesh%dphi(1)
            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
            do  phi_p_index = 1, L_p
                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                fac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2)
                pfac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)
                Norm_PNP = Norm_PNP + fac
                pNorm_PNP = pNorm_PNP + pfac
                do i = 1, ngr*ntheta*nphi
                    Particle_PNP(1) = Particle_PNP(1) + fac*mix%rho_V_it(i,phi_n_index,1)
                    Particle_PNP(2) = Particle_PNP(2) + fac*mix%rho_V_it(i,phi_p_index,2)
                    pParticle_PNP(1) = pParticle_PNP(1) + pfac*mix%prho_V_it(i,phi_n_index,1)
                    pParticle_PNP(2) = pParticle_PNP(2) + pfac*mix%prho_V_it(i,phi_p_index,2)
                end do
            end do
        end do
        ! write(*,*) "N:", Particle_PNP(1)/Norm_PNP, "Z:", Particle_PNP(2)/Norm_PNP
        ! write(*,*) "pN:", pParticle_PNP(1)/pNorm_PNP, "pZ:", pParticle_PNP(2)/pNorm_PNP
    end subroutine

    subroutine calculate_Rotated_Energy_after_PNP(Etot_PNP,pEtot_PNP)
        !-----------------------------------------------------------------
        !    Etot_PNP: <q_1| H R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >  
        !   pEtot_PNP: <q_1| H R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !------------------------------------------------------------------
        !    E = E_kin + E_EDF + E_pair + E_cou + E_cm
        !------------------------------------------------------------------
        use Globals, only: projection_mesh,nucleus_attributes,mix
        use Energy, only: Kinetic_term,EDF_terms,Coulomb_term,Pairing_term,Center_of_Mass_Correction_term
        complex(r64),intent(out) :: Etot_PNP,pEtot_PNP
        integer :: L_n,L_p,phi_n_index, phi_p_index
        real(r64) :: phi_n,phi_p
        complex(r64) :: E_cou,pE_cou,E_kin,pE_kin,E_EDF,pE_EDF,E_pair,pE_pair,E_cm,emiNphi,emiZphi
        Etot_PNP = (0.0d0,0.0d0)
        pEtot_PNP= (0.0d0,0.0d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)
        do  phi_p_index = 1, L_p
            phi_p =  phi_p_index*projection_mesh%dphi(2) 
            emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
            call Coulomb_term(phi_p_index,E_cou,pE_cou)
            do  phi_n_index = 1, L_n
                phi_n =  phi_n_index*projection_mesh%dphi(1)
                emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
                call Kinetic_term(phi_n_index,phi_p_index,E_kin,pE_kin)
                call EDF_terms(phi_n_index,phi_p_index,E_EDF,pE_EDF)
                call Pairing_term(phi_n_index,phi_p_index,E_pair,pE_pair)
                call Center_of_Mass_Correction_term(E_cm)
                Etot_PNP = Etot_PNP+ 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2)*&
                    (E_kin+E_EDF+E_pair+E_cou+E_cm)
                pEtot_PNP= pEtot_PNP+ 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)*&
                    (pE_kin+pE_EDF+pE_pair+pE_cou+E_cm)
            end do
        end do
        contains
    end subroutine

    subroutine calculate_N2_after_PNP(N2_PNP, pN2_PNP)
        !----------------------------------------------------------------
        !  N2_PNP: <q_1| N^2 R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >
        ! pN2_PNP: <q_1| N^2 R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !---------------------------------------------------------------
        use Globals, only: projection_mesh,nucleus_attributes,mix
        use JNsquare, only: calculate_N2
        complex(r64),intent(out) :: N2_PNP(2), pN2_PNP(2)
        integer :: L_n,L_p,phi_n_index,phi_p_index,it
        real(r64) :: phi_n,phi_p
        complex(r64) :: N2,pN2,emiNphi,emiZphi,fac,pfac
        complex(r64), allocatable :: N2_arry(:,:), pN2_arry(:,:)

        N2_PNP = (0.d0, 0.d0)
        pN2_PNP = (0.d0, 0.d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)

        ! calculate <q1| N^2 R |q2>/<q1|R|q2> for proton and neutron separately.
        ! where R = R(alpha,beta,gamma, phi_n,phi_p) = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
        allocate(N2_arry(max(L_n,L_p),2),pN2_arry(max(L_n,L_p),2))
        do phi_n_index = 1, L_n
            it = 1
            call calculate_N2(phi_n_index,it,N2,pN2)
            N2_arry(phi_n_index,it) = N2
            pN2_arry(phi_n_index,it) = pN2
        end do 
        do phi_p_index = 1, L_p
            it = 2
            call calculate_N2(phi_p_index,it,N2,pN2)
            N2_arry(phi_p_index,it) = N2
            pN2_arry(phi_p_index,it) = pN2
        end do 

        do  phi_n_index = 1, L_n
            phi_n =  phi_n_index*projection_mesh%dphi(1)
            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
            do  phi_p_index = 1, L_p
                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                fac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2)
                pfac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)
                it = 1
                N2_PNP(it) = N2_PNP(it) + fac*N2_arry(phi_n_index,it)
                pN2_PNP(it) = pN2_PNP(it) + pfac*pN2_arry(phi_n_index,it)
                it= 2
                N2_PNP(it) = N2_PNP(it) + fac*N2_arry(phi_p_index,it)
                pN2_PNP(it) = pN2_PNP(it) + pfac*pN2_arry(phi_p_index,it)
            end do
        end do
        deallocate(N2_arry,pN2_arry)
    end subroutine

    subroutine calculate_J2_after_PNP(J2_PNP, pJ2_PNP)
        !----------------------------------------------------------------
        !  J2_PNP: <q_1| J^2 R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >
        ! pJ2_PNP: <q_1| J^2 R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !---------------------------------------------------------------
        use Globals, only: projection_mesh,nucleus_attributes,mix
        use JNsquare, only: calculate_Jsquare_and_J_MM
        complex(r64),intent(out) :: J2_PNP, pJ2_PNP
        integer :: L_n,L_p,phi_n_index,phi_p_index,it
        real(r64) :: phi_n,phi_p
        complex(r64) :: J2,pJ2,emiNphi,emiZphi,fac,pfac,J_i(3),pJ_i(3)
        complex(r64), allocatable :: Jsquare_arry(:,:), pJsquare_arry(:,:),J_array(:,:,:),pJ_array(:,:,:)

        J2_PNP = (0.d0, 0.d0)
        pJ2_PNP = (0.d0, 0.d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)

        ! calculate <q1| J^2 R |q2>/<q1|R|q2> and  <q1| J_i R |q2>/<q1|R|q2> for proton and neutron separately.
        ! where R = R(alpha,beta,gamma, phi_n,phi_p) = e^{i alpha J_z} e^{i beta J_y} e^{i gamma J_z} e^{i phi_n N} e^{i phi_p N}
        allocate(Jsquare_arry(max(L_n,L_p),2),pJsquare_arry(max(L_n,L_p),2),J_array(max(L_n,L_p),2,3),pJ_array(max(L_n,L_p),2,3))
        do phi_n_index = 1, L_n
            it = 1
            call calculate_Jsquare_and_J_MM(phi_n_index,it,J2,pJ2,J_i,pJ_i)
            Jsquare_arry(phi_n_index,it) = J2
            pJsquare_arry(phi_n_index,it) = pJ2
            J_array(phi_n_index,it,:) = J_i(:)
            pJ_array(phi_n_index,it,:) = pJ_i(:)
        end do 
        do phi_p_index = 1, L_p
            it = 2
            call calculate_Jsquare_and_J_MM(phi_p_index,it,J2,pJ2,J_i,pJ_i)
            Jsquare_arry(phi_p_index,it) = J2
            pJsquare_arry(phi_p_index,it) = pJ2
            J_array(phi_p_index,it,:) = J_i(:)
            pJ_array(phi_p_index,it,:) = pJ_i(:)
        end do 

        do  phi_n_index = 1, L_n
            phi_n =  phi_n_index*projection_mesh%dphi(1)
            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
            do  phi_p_index = 1, L_p
                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                fac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2)
                pfac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)
                ! <J^2_n> = <Jx^2_n> + <Jy^2_n> + <Jz^2_n>
                it = 1
                J2_PNP = J2_PNP + fac*Jsquare_arry(phi_n_index,it)
                pJ2_PNP = pJ2_PNP + pfac*pJsquare_arry(phi_n_index,it)
                ! <J^2_p> = <Jx^2_p> + <Jy^2_p> + <Jz^2_p>
                it= 2
                J2_PNP = J2_PNP + fac*Jsquare_arry(phi_p_index,it)
                pJ2_PNP = pJ2_PNP + pfac*pJsquare_arry(phi_p_index,it)
                ! 2<J_p J_n> = 2<Jx>_n<Jx>_p + 2<Jy>_n<Jy>_p + 2<Jz>_n<Jz>_p
                J2_PNP = J2_PNP + fac*2.d0*J_array(phi_n_index,1,1)*J_array(phi_p_index,2,1) !  2<Jx>_n<Jx>_p
                J2_PNP = J2_PNP + fac*2.d0*J_array(phi_n_index,1,2)*J_array(phi_p_index,2,2) !  2<Jy>_n<Jy>_p
                J2_PNP = J2_PNP + fac*2.d0*J_array(phi_n_index,1,3)*J_array(phi_p_index,2,3) !  2<Jz>_n<Jz>_p
                pJ2_PNP = pJ2_PNP + pfac*2.d0*pJ_array(phi_n_index,1,1)*pJ_array(phi_p_index,2,1) !  2<Jx>_n<Jx>_p
                pJ2_PNP = pJ2_PNP + pfac*2.d0*pJ_array(phi_n_index,1,2)*pJ_array(phi_p_index,2,2) !  2<Jy>_n<Jy>_p
                pJ2_PNP = pJ2_PNP + pfac*2.d0*pJ_array(phi_n_index,1,3)*pJ_array(phi_p_index,2,3) !  2<Jz>_n<Jz>_p
            end do
        end do
        deallocate(Jsquare_arry,pJsquare_arry,J_array,pJ_array)
    end subroutine

    subroutine calculate_Qlm_after_PNP(Q2m_PNP,pQ2m_PNP,cQ2m_PNP,pcQ2m_PNP)
        !---------------------------------------------------------------------------------
        !   Q2m_PNP: <q_1| Q_{lambda mu} R(alpha,beta,gamma)      P^{N} P^{Z}  |q_2 >      
        !  pQ2m_PNP: <q_1| Q_{lambda mu} R(alpha,beta,gamma)      P^{N} P^{Z} P|q_2 >
        !  cQ2m_PNP: <q_1| Q^{\dagger}_{2 mu} R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >
        ! pcQ2m_PNP: <q_1| Q^{\dagger}_{2 mu} R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !----------------------------------------------------------------------------------
        use Globals, only: projection_mesh,nucleus_attributes,mix
        use EM, only: calculate_Qlm
        integer :: L_n,L_p,phi_n_index, phi_p_index,it,m
        real(r64) :: phi_n,phi_p
        complex(r64) :: emiNphi,emiZphi,fac,pfac
        complex(r64), dimension(-2:2,2) :: Q2m,pQ2m,cQ2m,pcQ2m,Q2m_PNP,pQ2m_PNP,cQ2m_PNP,pcQ2m_PNP
        complex(r64), dimension(:,:,:), allocatable :: Q2m_array,pQ2m_array,cQ2m_array,pcQ2m_array

        Q2m_PNP = (0.d0, 0.d0)
        pQ2m_PNP = (0.d0, 0.d0)
        cQ2m_PNP = (0.d0, 0.d0)
        pcQ2m_PNP = (0.d0, 0.d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)

        allocate(Q2m_array(max(L_n,L_p),-2:2,2),pQ2m_array(max(L_n,L_p),-2:2,2),cQ2m_array(max(L_n,L_p),-2:2,2),pcQ2m_array(max(L_n,L_p),-2:2,2))
        do phi_n_index = 1, L_n
            it = 1
            call calculate_Qlm(2,phi_n_index,it,Q2m,pQ2m,cQ2m,pcQ2m)
            Q2m_array(phi_n_index,:,it) = Q2m(:,it)
            pQ2m_array(phi_n_index,:,it) = pQ2m(:,it)
            cQ2m_array(phi_n_index,:,it) = cQ2m(:,it)
            pcQ2m_array(phi_n_index,:,it) = pcQ2m(:,it)
        end do 
        do phi_p_index = 1, L_p
            it = 2
            call calculate_Qlm(2,phi_p_index,it,Q2m,pQ2m,cQ2m,pcQ2m)
            Q2m_array(phi_p_index,:,it) = Q2m(:,it)
            pQ2m_array(phi_p_index,:,it) = pQ2m(:,it)
            cQ2m_array(phi_p_index,:,it) = cQ2m(:,it)
            pcQ2m_array(phi_p_index,:,it) = pcQ2m(:,it)
        end do 

        do  phi_n_index = 1, L_n
            phi_n =  phi_n_index*projection_mesh%dphi(1)
            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
            do  phi_p_index = 1, L_p
                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                fac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2)
                pfac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)
                ! Q2: l=2
                it = 1 
                do m = -2, 2
                    Q2m_PNP(m,it) = Q2m_PNP(m,it) + fac*Q2m_array(phi_n_index,m,it)
                    pQ2m_PNP(m,it) = pQ2m_PNP(m,it) + pfac*pQ2m_array(phi_n_index,m,it)
                    cQ2m_PNP(m,it) = cQ2m_PNP(m,it) + fac*cQ2m_array(phi_n_index,m,it)
                    pcQ2m_PNP(m,it) = pcQ2m_PNP(m,it) + pfac*pcQ2m_array(phi_n_index,m,it)
                end do 
                it = 2
                do m = -2, 2
                    Q2m_PNP(m,it) = Q2m_PNP(m,it) + fac*Q2m_array(phi_p_index,m,it)
                    pQ2m_PNP(m,it) = pQ2m_PNP(m,it) + pfac*pQ2m_array(phi_p_index,m,it)
                    cQ2m_PNP(m,it) = cQ2m_PNP(m,it) + fac*cQ2m_array(phi_p_index,m,it)
                    pcQ2m_PNP(m,it) = pcQ2m_PNP(m,it) + pfac*pcQ2m_array(phi_p_index,m,it)
                end do
            end do
        end do
        deallocate(Q2m_array,pQ2m_array,cQ2m_array,pcQ2m_array)
    end subroutine

    subroutine calculate_r2_after_PNP(r2_PNP, pr2_PNP)
        !----------------------------------------------------------------
        !  r2_PNP: <q_1| r^2 R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >
        ! pr2_PNP: <q_1| r^2 R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !---------------------------------------------------------------

        use Globals, only: projection_mesh,nucleus_attributes,mix
        use EM, only: calculate_r2
        integer :: L_n,L_p,phi_n_index,phi_p_index,it
        real(r64) :: phi_n,phi_p
        complex(r64) :: emiNphi,emiZphi,fac,pfac
        complex(r64) :: r2,pr2,r2_PNP(2),pr2_PNP(2)
        complex(r64), dimension(:,:), allocatable :: r2_arry, pr2_arry

        r2_PNP = (0.d0, 0.d0)
        pr2_PNP = (0.d0, 0.d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)

        allocate(r2_arry(max(L_n,L_p),2),pr2_arry(max(L_n,L_p),2))
        do phi_n_index = 1, L_n
            it = 1
            call calculate_r2(phi_n_index,it,r2,pr2)
            r2_arry(phi_n_index,it) = r2
            pr2_arry(phi_n_index,it) = pr2
        end do 
        do phi_p_index = 1, L_p
            it = 2
            call calculate_r2(phi_p_index,it,r2,pr2)
            r2_arry(phi_p_index,it) = r2
            pr2_arry(phi_p_index,it) = pr2
        end do 

        do  phi_n_index = 1, L_n
            phi_n =  phi_n_index*projection_mesh%dphi(1)
            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
            do  phi_p_index = 1, L_p
                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                fac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2)
                pfac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)
                !
                it = 1
                r2_PNP(it) = r2_PNP(it) + fac*r2_arry(phi_n_index,it)
                pr2_PNP(it) = pr2_PNP(it) + pfac*pr2_arry(phi_n_index,it)
                it= 2
                r2_PNP(it) = r2_PNP(it) + fac*r2_arry(phi_p_index,it)
                pr2_PNP(it) = pr2_PNP(it) + pfac*pr2_arry(phi_p_index,it)
            end do
        end do
        deallocate(r2_arry,pr2_arry)
    end subroutine


    subroutine calculate_r2_2body_after_PNP(r2_2b_PNP, pr2_2b_PNP)
        !----------------------------------------------------------------
        !  r2_2b_PNP: <q_1| E_n R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >
        ! pr2_2b_PNP: <q_1| E_n R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !---------------------------------------------------------------

        use Globals, only: projection_mesh,nucleus_attributes,mix
        use Eccentricity, only: calculate_Eccentri_n_MM
        complex(r64),intent(out) :: r2_2b_PNP(5,3),pr2_2b_PNP(5,3)
        complex(r64) :: r2_2b_Each_Term_PNP(3,3),pr2_2b_Each_Term_PNP(3,3)
        integer :: L_n,L_p,phi_n_index,phi_p_index,it,mu
        real(r64) :: phi_n,phi_p
        complex(r64) :: emiNphi,emiZphi,fac,pfac
        integer,parameter :: n = 1
        complex(r64) :: Eccentri(2),pEccentri(2),Qn_mu(-n:n),pQn_mu(-n:n),Each_2B_Term(3),pEach_2B_Term(3)
        complex(r64), dimension(:,:,:),allocatable :: r2_2b_arry, pr2_2b_arry,Qn_mu_arry,pQn_mu_arry,Each_2B_Term_arry,pEach_2B_Term_arry
        r2_2b_PNP = (0.d0, 0.d0)
        pr2_2b_PNP = (0.d0, 0.d0)
        r2_2b_Each_Term_PNP = (0.d0, 0.d0)
        pr2_2b_Each_Term_PNP = (0.d0, 0.d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)

        allocate(r2_2b_arry(2,max(L_n,L_p),2),pr2_2b_arry(2,max(L_n,L_p),2), &
                 Qn_mu_arry(-n:n,max(L_n,L_p),2),pQn_mu_arry(-n:n,max(L_n,L_p),2), &
                 Each_2B_Term_arry(3,max(L_n,L_p),2), pEach_2B_Term_arry(3,max(L_n,L_p),2))
        do phi_n_index = 1, L_n
            it = 1
            call calculate_Eccentri_n_MM(n,phi_n_index,it,Eccentri,pEccentri,Qn_mu,pQn_mu,Each_2B_Term,pEach_2B_Term)
            r2_2b_arry(:,phi_n_index,it) = Eccentri(:)  ! For every phi_n, return 1b :1, same particle 2b:2
            pr2_2b_arry(:,phi_n_index,it) = pEccentri(:) 
            Qn_mu_arry(:,phi_n_index,it) = Qn_mu(:) 
            pQn_mu_arry(:,phi_n_index,it) = pQn_mu(:)
            Each_2B_Term_arry(:,phi_n_index,it) = Each_2B_Term(:)
            pEach_2B_Term_arry(:,phi_n_index,it) = pEach_2B_Term(:)
        end do
        do phi_p_index = 1, L_p
            it = 2
            call calculate_Eccentri_n_MM(n,phi_p_index,it,Eccentri,pEccentri,Qn_mu,pQn_mu,Each_2B_Term,pEach_2B_Term)
            r2_2b_arry(:,phi_p_index,it) = Eccentri(:)
            pr2_2b_arry(:,phi_p_index,it) = pEccentri(:)
            Qn_mu_arry(:,phi_p_index,it) = Qn_mu(:)
            pQn_mu_arry(:,phi_p_index,it) = pQn_mu(:)
            Each_2B_Term_arry(:,phi_p_index,it) = Each_2B_Term(:)
            pEach_2B_Term_arry(:,phi_p_index,it) = pEach_2B_Term(:)
        end do 

        do  phi_n_index = 1, L_n
            phi_n =  phi_n_index*projection_mesh%dphi(1)
            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
            do  phi_p_index = 1, L_p
                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                fac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2) ! 1/(L_n*L_p) * e^{-iN*phi_n} * e^{-iZ*phi_p} * norm_n * norm_p
                pfac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)
                ! 1B
                ! do it = 1, 2 ! 中子，质子
                it = 1 ! 中子

                r2_2b_PNP(1,it) = r2_2b_PNP(1,it) + fac*r2_2b_arry(1,phi_n_index,it) 
                pr2_2b_PNP(1,it) = pr2_2b_PNP(1,it) + pfac*pr2_2b_arry(1,phi_n_index,it)
                it = 2  ! 质子
                r2_2b_PNP(1,it) = r2_2b_PNP(1,it) + fac*r2_2b_arry(1,phi_p_index,it)
                pr2_2b_PNP(1,it) = pr2_2b_PNP(1,it) + pfac*pr2_2b_arry(1,phi_p_index,it)

                ! 2B
                it = 1
                ! 
                r2_2b_PNP(2,it) = r2_2b_PNP(2,it) + fac*r2_2b_arry(2,phi_n_index,it)
                pr2_2b_PNP(2,it) = pr2_2b_PNP(2,it) + pfac*pr2_2b_arry(2,phi_n_index,it)
                ! each teram of nn
                r2_2b_Each_Term_PNP(1,it) = r2_2b_Each_Term_PNP(1,it) + fac*Each_2B_Term_arry(1,phi_n_index,it)
                r2_2b_Each_Term_PNP(2,it) = r2_2b_Each_Term_PNP(2,it) + fac*Each_2B_Term_arry(2,phi_n_index,it)
                r2_2b_Each_Term_PNP(3,it) = r2_2b_Each_Term_PNP(3,it) + fac*Each_2B_Term_arry(3,phi_n_index,it)
                pr2_2b_Each_Term_PNP(1,it) = pr2_2b_Each_Term_PNP(1,it) + pfac*pEach_2B_Term_arry(1,phi_n_index,it)
                pr2_2b_Each_Term_PNP(2,it) = pr2_2b_Each_Term_PNP(2,it) + pfac*pEach_2B_Term_arry(2,phi_n_index,it)
                pr2_2b_Each_Term_PNP(3,it) = pr2_2b_Each_Term_PNP(3,it) + pfac*pEach_2B_Term_arry(3,phi_n_index,it)
                it = 2
                r2_2b_PNP(2,it) = r2_2b_PNP(2,it) + fac*r2_2b_arry(2,phi_p_index,it)
                pr2_2b_PNP(2,it) = pr2_2b_PNP(2,it) + pfac*pr2_2b_arry(2,phi_p_index,it)
                ! each teram of pp
                r2_2b_Each_Term_PNP(1,it) = r2_2b_Each_Term_PNP(1,it) + fac*Each_2B_Term_arry(1,phi_p_index,it)
                r2_2b_Each_Term_PNP(2,it) = r2_2b_Each_Term_PNP(2,it) + fac*Each_2B_Term_arry(2,phi_p_index,it)
                r2_2b_Each_Term_PNP(3,it) = r2_2b_Each_Term_PNP(3,it) + fac*Each_2B_Term_arry(3,phi_p_index,it)
                pr2_2b_Each_Term_PNP(1,it) = pr2_2b_Each_Term_PNP(1,it) +pfac*pEach_2B_Term_arry(1,phi_p_index,it)
                pr2_2b_Each_Term_PNP(2,it) = pr2_2b_Each_Term_PNP(2,it) + pfac*pEach_2B_Term_arry(2,phi_p_index,it)
                pr2_2b_Each_Term_PNP(3,it) = pr2_2b_Each_Term_PNP(3,it) + pfac*pEach_2B_Term_arry(3,phi_p_index,it)

                do mu = -n,n
                    r2_2b_PNP(2,3) = r2_2b_PNP(2,3) + fac*(-1)**mu* &
                            (Qn_mu_arry(mu,phi_n_index,1)*Qn_mu_arry(-mu,phi_p_index,2) + Qn_mu_arry(mu,phi_p_index,2)*Qn_mu_arry(-mu,phi_n_index,1))
                    pr2_2b_PNP(2,3) = pr2_2b_PNP(2,3) + pfac*(-1)**mu* &
                            (pQn_mu_arry(mu,phi_n_index,1)*pQn_mu_arry(-mu,phi_p_index,2) + pQn_mu_arry(mu,phi_p_index,2)*pQn_mu_arry(-mu,phi_n_index,1))
                    ! direct term should include `np and pn`
                    r2_2b_Each_Term_PNP(1,3) = r2_2b_Each_Term_PNP(1,3) + fac*(-1)**mu* &
                            (Qn_mu_arry(mu,phi_n_index,1)*Qn_mu_arry(-mu,phi_p_index,2) + Qn_mu_arry(mu,phi_p_index,2)*Qn_mu_arry(-mu,phi_n_index,1))
                    pr2_2b_Each_Term_PNP(1,3) = pr2_2b_Each_Term_PNP(1,3) + pfac*(-1)**mu* &
                            (pQn_mu_arry(mu,phi_n_index,1)*pQn_mu_arry(-mu,phi_p_index,2) + pQn_mu_arry(mu,phi_p_index,2)*pQn_mu_arry(-mu,phi_n_index,1))
                end do
            end do
        end do
        do it = 1, 3    ! ← 改为 3
            r2_2b_PNP(3,it) = r2_2b_Each_Term_PNP(1,it)
            r2_2b_PNP(4,it) = r2_2b_Each_Term_PNP(2,it)
            r2_2b_PNP(5,it) = r2_2b_Each_Term_PNP(3,it)
            pr2_2b_PNP(3,it) = pr2_2b_Each_Term_PNP(1,it)
            pr2_2b_PNP(4,it) = pr2_2b_Each_Term_PNP(2,it)
            pr2_2b_PNP(5,it) = pr2_2b_Each_Term_PNP(3,it)
        end do

        deallocate(r2_2b_arry,pr2_2b_arry,Qn_mu_arry,pQn_mu_arry,Each_2B_Term_arry,pEach_2B_Term_arry)
    end subroutine

    subroutine calculate_Eccentricity_after_PNP(Eccentri_PNP,pEccentri_PNP)
        !----------------------------------------------------------------
        !  Eccentri_PNP: <q_1| E_n R(alpha,beta,gamma) P^{N} P^{Z}  |q_2 >
        ! pEccentri_PNP: <q_1| E_n R(alpha,beta,gamma) P^{N} P^{Z} P|q_2 >
        !---------------------------------------------------------------

        use Globals, only: projection_mesh,nucleus_attributes,mix
        use Eccentricity, only: calculate_Eccentri_n_MM
        complex(r64),intent(out) :: Eccentri_PNP(5),pEccentri_PNP(5)
        complex(r64) :: Eccentri_Each_Term_PNP(3),pEccentri_Each_Term_PNP(3)
        integer :: L_n,L_p,phi_n_index,phi_p_index,it,mu ! Neutron/proton gauge angle Grid points
        real(r64) :: phi_n,phi_p ! Neutron/proton gauge angle
        complex(r64) :: emiNphi,emiZphi,fac,pfac
        integer,parameter :: n = 1
        complex(r64) :: Eccentri(2),pEccentri(2),Qn_mu(-n:n),pQn_mu(-n:n),Each_2B_Term(3),pEach_2B_Term(3)
        complex(r64), dimension(:,:,:),allocatable :: Eccentri_arry, pEccentri_arry,Qn_mu_arry,pQn_mu_arry,Each_2B_Term_arry,pEach_2B_Term_arry
        Eccentri_PNP = (0.d0, 0.d0)
        pEccentri_PNP = (0.d0, 0.d0)
        Eccentri_Each_Term_PNP = (0.d0, 0.d0)
        pEccentri_Each_Term_PNP = (0.d0, 0.d0)
        L_n = projection_mesh%nphi(1)
        L_p = projection_mesh%nphi(2)

        allocate(Eccentri_arry(2,max(L_n,L_p),2),pEccentri_arry(2,max(L_n,L_p),2),Qn_mu_arry(-n:n,max(L_n,L_p),2),pQn_mu_arry(-n:n,max(L_n,L_p),2), &
                 Each_2B_Term_arry(3,max(L_n,L_p),2), pEach_2B_Term_arry(3,max(L_n,L_p),2))
        do phi_n_index = 1, L_n
            it = 1
            call calculate_Eccentri_n_MM(n,phi_n_index,it,Eccentri,pEccentri,Qn_mu,pQn_mu,Each_2B_Term,pEach_2B_Term)
            Eccentri_arry(:,phi_n_index,it) = Eccentri(:)  ! For every phi_n, return 1b :1, same particle 2b:2
            pEccentri_arry(:,phi_n_index,it) = pEccentri(:) 
            Qn_mu_arry(:,phi_n_index,it) = Qn_mu(:) ! For every phi_n, return 1b :1, same particle 2b:2 sum_{pq} (Q_mu) rho_qp
            pQn_mu_arry(:,phi_n_index,it) = pQn_mu(:)
            Each_2B_Term_arry(:,phi_n_index,it) = Each_2B_Term(:)
            pEach_2B_Term_arry(:,phi_n_index,it) = pEach_2B_Term(:)
        end do 
        do phi_p_index = 1, L_p
            it = 2
            call calculate_Eccentri_n_MM(n,phi_p_index,it,Eccentri,pEccentri,Qn_mu,pQn_mu,Each_2B_Term,pEach_2B_Term)
            Eccentri_arry(:,phi_p_index,it) = Eccentri(:)
            pEccentri_arry(:,phi_p_index,it) = pEccentri(:)
            Qn_mu_arry(:,phi_p_index,it) = Qn_mu(:)
            pQn_mu_arry(:,phi_p_index,it) = pQn_mu(:)
            Each_2B_Term_arry(:,phi_p_index,it) = Each_2B_Term(:)
            pEach_2B_Term_arry(:,phi_p_index,it) = pEach_2B_Term(:)
        end do 

        do  phi_n_index = 1, L_n
            phi_n =  phi_n_index*projection_mesh%dphi(1)
            emiNphi = cdexp(-nucleus_attributes%neutron_number*cmplx(0,phi_n)) ! e^{-iN\phi_n}
            do  phi_p_index = 1, L_p
                phi_p =  phi_p_index*projection_mesh%dphi(2) 
                emiZphi = cdexp(-nucleus_attributes%proton_number*cmplx(0,phi_p)) ! e^{-iZ\phi_p}
                fac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%norm(phi_n_index,1)*mix%norm(phi_p_index,2) ! 1/(L_n*L_p) * e^{-iN*phi_n} * e^{-iZ*phi_p} * norm_n * norm_p
                pfac = 1.d0/(L_n*L_p)*emiNphi*emiZphi*mix%pnorm(phi_n_index,1)*mix%pnorm(phi_p_index,2)
                ! 1B
                it = 1
                Eccentri_PNP(1) = Eccentri_PNP(1) + fac*Eccentri_arry(1,phi_n_index,it)
                pEccentri_PNP(1) = pEccentri_PNP(1) + pfac*pEccentri_arry(1,phi_n_index,it)
                it = 2
                Eccentri_PNP(1) = Eccentri_PNP(1) + fac*Eccentri_arry(1,phi_p_index,it)
                pEccentri_PNP(1) = pEccentri_PNP(1) + pfac*pEccentri_arry(1,phi_p_index,it)

                ! 2B
                it = 1
                Eccentri_PNP(2) = Eccentri_PNP(2) - fac*Eccentri_arry(2,phi_n_index,it)
                pEccentri_PNP(2) = pEccentri_PNP(2) - pfac*pEccentri_arry(2,phi_n_index,it)
                ! each teram of nn
                Eccentri_Each_Term_PNP(1) = Eccentri_Each_Term_PNP(1) - fac*Each_2B_Term_arry(1,phi_n_index,it) 
                Eccentri_Each_Term_PNP(2) = Eccentri_Each_Term_PNP(2) - fac*Each_2B_Term_arry(2,phi_n_index,it)
                Eccentri_Each_Term_PNP(3) = Eccentri_Each_Term_PNP(3) - fac*Each_2B_Term_arry(3,phi_n_index,it)
                pEccentri_Each_Term_PNP(1) = pEccentri_Each_Term_PNP(1) - pfac*pEach_2B_Term_arry(1,phi_n_index,it)
                pEccentri_Each_Term_PNP(2) = pEccentri_Each_Term_PNP(2) - pfac*pEach_2B_Term_arry(2,phi_n_index,it)
                pEccentri_Each_Term_PNP(3) = pEccentri_Each_Term_PNP(3) - pfac*pEach_2B_Term_arry(3,phi_n_index,it)
                it = 2
                Eccentri_PNP(2) = Eccentri_PNP(2) - fac*Eccentri_arry(2,phi_p_index,it)
                pEccentri_PNP(2) = pEccentri_PNP(2) - pfac*pEccentri_arry(2,phi_p_index,it)
                ! each teram of nn
                Eccentri_Each_Term_PNP(1) = Eccentri_Each_Term_PNP(1) - fac*Each_2B_Term_arry(1,phi_p_index,it)
                Eccentri_Each_Term_PNP(2) = Eccentri_Each_Term_PNP(2) - fac*Each_2B_Term_arry(2,phi_p_index,it)
                Eccentri_Each_Term_PNP(3) = Eccentri_Each_Term_PNP(3) - fac*Each_2B_Term_arry(3,phi_p_index,it)
                pEccentri_Each_Term_PNP(1) = pEccentri_Each_Term_PNP(1) - pfac*pEach_2B_Term_arry(1,phi_p_index,it)
                pEccentri_Each_Term_PNP(2) = pEccentri_Each_Term_PNP(2) - pfac*pEach_2B_Term_arry(2,phi_p_index,it)
                pEccentri_Each_Term_PNP(3) = pEccentri_Each_Term_PNP(3) - pfac*pEach_2B_Term_arry(3,phi_p_index,it)

                do mu = -n,n
                    Eccentri_PNP(2) = Eccentri_PNP(2) + fac*(-1)**mu* &
                            (Qn_mu_arry(mu,phi_n_index,1)*Qn_mu_arry(-mu,phi_p_index,2) + Qn_mu_arry(mu,phi_p_index,2)*Qn_mu_arry(-mu,phi_n_index,1))
                    pEccentri_PNP(2) = pEccentri_PNP(2) + pfac*(-1)**mu* &
                            (pQn_mu_arry(mu,phi_n_index,1)*pQn_mu_arry(-mu,phi_p_index,2) + pQn_mu_arry(mu,phi_p_index,2)*pQn_mu_arry(-mu,phi_n_index,1))
                    ! direct term should include `np and pn`
                    Eccentri_Each_Term_PNP(1) = Eccentri_Each_Term_PNP(1) + fac*(-1)**mu* &
                            (Qn_mu_arry(mu,phi_n_index,1)*Qn_mu_arry(-mu,phi_p_index,2) + Qn_mu_arry(mu,phi_p_index,2)*Qn_mu_arry(-mu,phi_n_index,1))
                    pEccentri_Each_Term_PNP(1) = pEccentri_Each_Term_PNP(1) + pfac*(-1)**mu* &
                            (pQn_mu_arry(mu,phi_n_index,1)*pQn_mu_arry(-mu,phi_p_index,2) + pQn_mu_arry(mu,phi_p_index,2)*pQn_mu_arry(-mu,phi_n_index,1))
                end do
            end do
        end do
        Eccentri_PNP(3) = Eccentri_Each_Term_PNP(1)
        Eccentri_PNP(4) = Eccentri_Each_Term_PNP(2)
        Eccentri_PNP(5) = Eccentri_Each_Term_PNP(3)
        pEccentri_PNP(3) = pEccentri_Each_Term_PNP(1)
        pEccentri_PNP(4) = pEccentri_Each_Term_PNP(2)
        pEccentri_PNP(5) = pEccentri_Each_Term_PNP(3)

        deallocate(Eccentri_arry,pEccentri_arry,Qn_mu_arry,pQn_mu_arry,Each_2B_Term_arry,pEach_2B_Term_arry)
    end subroutine
END MODULE Kernel