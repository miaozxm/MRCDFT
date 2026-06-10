!==============================================================================!
! MODULE Observables                                                                 !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE GCM_Observables

    implicit none
    contains
    subroutine calculate_observables
        use GCM_Inout, only: write_GCM_observables
        call calculate_spectrum
        call calculate_reduced_transition_rate
        call write_GCM_observables
    end subroutine
    
    subroutine calculate_spectrum
        !-----------------------------------------------------------------------------
        ! <O> = \sum_{K1 q1 K2, q2 } f^J(K1,q1)*f^J(K2,q2)*N_KK(J,Pi(+/-),q1,K1,q2,K2)
        !----------------------------------------------------------------------------
        use Constants, only: r64
        use Globals, only: GCM_HWG,gcm_space,GCM_kernels,constraint,GCM_obser,Proj_option,GCM_basis
        integer :: J, parity,iM,qK1,q1,K1,qK2,q2,K2
        real(r64) :: beta2_aver,beta3_aver,be0,n_neu,n_pro,coef
        
        allocate(GCM_obser%E_ex(GCM_basis%N_max,0:gcm_space%Jmax,2),GCM_obser%beta2_aver(GCM_basis%N_max,0:gcm_space%Jmax,2),&
                 GCM_obser%beta3_aver(GCM_basis%N_max,0:gcm_space%Jmax,2),GCM_obser%N(GCM_basis%N_max,0:gcm_space%Jmax,2),&
                 GCM_obser%Z(GCM_basis%N_max,0:gcm_space%Jmax,2),   GCM_obser%rrms_p(GCM_basis%N_max,0:gcm_space%Jmax,2),source=0.d0)

        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            ! parity
            if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                parity = 1 ! +
            else
                parity = 2 ! -
            end if

            do iM = 1, GCM_HWG%M(J,parity)
                beta2_aver = 0.d0
                beta3_aver = 0.d0
                be0   = 0.d0
                n_neu = 0.d0
                n_pro = 0.d0
                do qK1 = 1, GCM_basis%N(J,parity)
                    q1 = GCM_basis%basis(1,qK1,J,parity)
                    K1 = GCM_basis%basis(2,qK1,J,parity)
                    beta2_aver = beta2_aver + GCM_HWG%gJKq(qK1,iM,J,parity)**2*constraint%betac(q1)
                    beta3_aver = beta3_aver + GCM_HWG%gJKq(qK1,iM,J,parity)**2*constraint%bet3c(q1)
                    do qK2 = 1, GCM_basis%N(J,parity)
                        q2 = GCM_basis%basis(1,qK2,J,parity)
                        K2 = GCM_basis%basis(2,qK2,J,parity)
                        coef = GCM_HWG%fJKq(qK1,iM,J,parity)*GCM_HWG%fJKq(qK2,iM,J,parity)
                        ! BE2 or <r^2>
                        be0   = be0   + coef*GCM_kernels%E0_KK(J,parity,q1,K1,q2,K2,2)
                        ! <N>
                        n_neu = n_neu + coef*GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,1)
                        ! <Z>
                        n_pro = n_pro + coef*GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,2)
                    end do
                end do
                GCM_obser%E_ex(iM,J,parity) = GCM_HWG%E(iM,J,parity) - GCM_HWG%E(1,0,1)
                GCM_obser%beta2_aver(iM,J,parity) = beta2_aver
                GCM_obser%beta3_aver(iM,J,parity) = beta3_aver
                GCM_obser%N(iM,J,parity) = n_neu
                GCM_obser%Z(iM,J,parity) = n_pro
                GCM_obser%rrms_p(iM,J,parity) = dsqrt(be0/(n_pro+1.E-10))
            end do
        end do 
    end subroutine

    subroutine calculate_reduced_transition_rate

    end subroutine
END MODULE