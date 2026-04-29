!==============================================================================!
! MODULE GCM_Inout                                                                  !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE GCM_Inout
    use Constants, only: u_start
    implicit none
    integer, private :: u_GCM = u_start + 21
    contains

    subroutine read_GCM_configuration(ifPrint)
        use Globals, only: InputFile,input_par,MPI_Infor
        integer :: i,is
        logical,intent(in),optional :: ifPrint
        character(len=*), parameter ::  format1 = "(10x,i5)", &
                                        format2 = "(10x,5e9.2)"
        open(u_GCM, file=InputFile%file_path_para, status='old')
        ! skip GCM parameters
        do i = 1, 49
            read(u_GCM,'(A)', iostat=is) 
            if (is /= 0) then
                print *, "Error reading GCM file"
                stop
            end if
        end do
        ! read
        read(u_GCM, format1) input_par%GCMType
        read(u_GCM, format1) input_par%kmax
        read(u_GCM, format2) input_par%zeta
        close(u_GCM)

        call set_GCM_parameters
        if(ifPrint.and. MPI_Infor%rank == 0) call printParameters

        contains

        subroutine set_GCM_parameters
            use Constants, only: Jmax_max
            use Globals, only: input_par,GCM_option,HWG
            integer :: J
            ! set GCM type
            GCM_option%GCMType = input_par%GCMType
            if(GCM_option%GCMType > 1 .or. GCM_option%GCMType < 0) stop 'GCMType wrong!'

            ! set kmax
            HWG%kmax = input_par%kmax
            if(input_par%kmax < 1) stop 'kmax cannot be less than 1.' 

            ! set cutoff
            do J = 0, Jmax_max
                if(J<=4) then
                    HWG%cutoff(J) = input_par%zeta(J+1)
                else 
                    HWG%cutoff(J) = input_par%zeta(5)
                end if 
            end do

        end subroutine
        subroutine printParameters
            use Globals, only: input_par,GCM_option
            use Tools, only: adjust_left
            integer :: Strlength = 40

            if(GCM_option%GCMType==0) then 
                write(*,"(5x,a,':   ',a)") adjust_left('GCM',Strlength),'only print the kernel'
            else if(GCM_option%GCMType==1) then
                write(*,"(5x,a,':   ',a)") adjust_left('GCM',Strlength),'perform GCM'
            end if 

            write(*,"(5x,a,':   ', 5(a,e9.2))") adjust_left('Zeta',Strlength),input_par%zeta
            write(*,"(a)") '=========================================================================================='
        end subroutine
    end subroutine

    subroutine read_kernels
        use Constants, only: r64
        use Proj_Inout, only: set_Proj_output_filename
        use Globals, only: gcm_space,Proj_outputfile,GCM_kernels,Proj_option
        integer :: q1, q2,J,K1_start,K1_end,K2_start,K2_end,K1,K2,q2_start,q2_end,parity
        integer :: read_J, read_K1, read_K2
        character(1) :: read_parity_char
        real(r64) :: read_N2, read_Z2, read_J2
        complex(r64) :: Q2_KK_21(2), E0_KK_21(2)
        character(1), dimension(2) :: ParityChar = ['+', '-']
        character(len=*), parameter ::  format1 = "(3i5,4x,a,4x,3f9.3)", &
                                        format2 = "(4e15.8)"
        ! N_KK(J,Pi(+/-),q1,K1,q2,K2), Norm kernel
        allocate(GCM_kernels%N_KK(0:gcm_space%Jmax,2,gcm_space%q1_start:gcm_space%q1_end,-gcm_space%Jmax:gcm_space%Jmax, &
                                    gcm_space%q2_start:gcm_space%q2_end,-gcm_space%Jmax:gcm_space%Jmax),source=(0.d0,0.d0))
        ! H_KK(J,Pi(+/-),q1,K1,q2,K2), Hamiltonian  kernel
        allocate(GCM_kernels%H_KK(0:gcm_space%Jmax,2,gcm_space%q1_start:gcm_space%q1_end,-gcm_space%Jmax:gcm_space%Jmax, &
                                    gcm_space%q2_start:gcm_space%q2_end,-gcm_space%Jmax:gcm_space%Jmax),source=(0.d0,0.d0))
        ! X_KK(J,Pi(+/-),q1,K1,q2,K2,it), Particle number kernel
        allocate(GCM_kernels%X_KK(0:gcm_space%Jmax,2,gcm_space%q1_start:gcm_space%q1_end,-gcm_space%Jmax:gcm_space%Jmax, &
                                    gcm_space%q2_start:gcm_space%q2_end,-gcm_space%Jmax:gcm_space%Jmax,2),source=(0.d0,0.d0))
        ! Q2_KK(J,Pi_i(+/-),q1,Ki,q2,Kf,it), <Ji+2 Kf q1 Pi_f ||Q2||Ji Ki q2 Pi_i>
        allocate(GCM_kernels%Q2_KK(0:gcm_space%Jmax,2,gcm_space%q1_start:gcm_space%q1_end,-gcm_space%Jmax:gcm_space%Jmax, &
                                    gcm_space%q2_start:gcm_space%q2_end,-gcm_space%Jmax:gcm_space%Jmax,2),source=(0.d0,0.d0))
        ! E0_KK(J,Pi_i(+/-),q1,Ki,q2,Kf,it), <J_f K_f q_1 Pi| r2 |J_i K_i q_2 Pi>
        allocate(GCM_kernels%E0_KK(0:gcm_space%Jmax,2,gcm_space%q1_start:gcm_space%q1_end,-gcm_space%Jmax:gcm_space%Jmax, &
                                    gcm_space%q2_start:gcm_space%q2_end,-gcm_space%Jmax:gcm_space%Jmax,2),source=(0.d0,0.d0))


        do q1 = gcm_space%q1_start, gcm_space%q1_end
            if(Proj_option%Kernel_Symmetry==0) then      ! All Kernels 
                q2_start = gcm_space%q2_start
                q2_end = gcm_space%q2_end
            else if(Proj_option%Kernel_Symmetry==1) then ! Triangular Kernels
                q2_start = q1
                q2_end = gcm_space%q2_end
            else if(Proj_option%Kernel_Symmetry==2) then ! Diagonal elements
                stop 'Only diagonal kernels cannot be used for GCM calculations.'
            else
                stop 'Kernel_Symmetry should be 0, 1, 2'
            end if
            do q2 = q2_start, q2_end
                call set_Proj_output_filename(q1,q2)
                open(Proj_outputfile%u_outputelem,file=Proj_outputfile%outputelem,form='unformatted',status='unknown')
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
                    do K1 = K1_start,K1_end
                        do K2 = K2_start,K2_end
                            ! In the axially symmetric case, the kernel is non-zero only when
                            ! the parity satisfies  Pi  = (-1)^J for  N_KK, H_KK, X_KK and E0_KK
                            ! the parity satisfies Pi_i = (-1)^J_i for  Q2_KK_12
                            if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                                parity = 1 ! +
                            else
                                parity = 2 ! -
                            end if
                            read(Proj_outputfile%u_outputelem,format1) read_J,read_K1,read_K2,read_parity_char,read_N2,read_Z2,read_J2
                            if(read_J/=J .or. read_K1/=K1 .or. read_K2/=K2 .or. read_parity_char/=ParityChar(parity)) then
                                write(*,*) Proj_outputfile%outputelem//': data is incorrect!' 
                                stop 'at read_kernels'
                            end if
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%N_KK(J,parity,q1,K1,q2,K2), GCM_kernels%H_KK(J,parity,q1,K1,q2,K2)
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,1), GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,2)
                            ! proton part
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%Q2_KK(J,parity,q1,K1,q2,K2,2), Q2_KK_21(2)
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%E0_KK(J,parity,q1,K1,q2,K2,2), E0_KK_21(2)
                            ! neutron part
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%Q2_KK(J,parity,q1,K1,q2,K2,1), Q2_KK_21(1)
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%E0_KK(J,parity,q1,K1,q2,K2,1), E0_KK_21(1)

                            GCM_kernels%H_KK(J,parity,q1,K1,q2,K2) = GCM_kernels%H_KK(J,parity,q1,K1,q2,K2)*GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                            GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,1) = GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,1)*GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                            GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,2) = GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,2)*GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                            if(Proj_option%Kernel_Symmetry==1 .and. q1/=q2) then
                                GCM_kernels%N_KK(J,parity,q2,K2,q1,K1) = GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                                GCM_kernels%H_KK(J,parity,q2,K2,q1,K1) = GCM_kernels%H_KK(J,parity,q1,K1,q2,K2)
                                GCM_kernels%X_KK(J,parity,q2,K2,q1,K1,:) = GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,:)
                                GCM_kernels%Q2_KK(J,parity,q2,K2,q1,K1,:) = Q2_KK_21(:)
                                GCM_kernels%E0_KK(J,parity,q2,K2,q1,K1,:) = E0_KK_21(:)
                            end if
                        end do 
                    end do 
                end do
                close(Proj_outputfile%u_outputelem)
            end do
        end do 
        call print_kernels
        contains 
        subroutine print_kernels
            use Globals, only: Proj_option,constraint
            integer :: J,parity,K1_start,K1_end,K2_start,K2_end,K1,q1,length_K,qK1,K2,q2,qK2,log_abs_norm
            real(r64) :: abs_norm
            character(1), dimension(2) :: ParityChar = ['+', '-']
            character(len=*), parameter ::  format11 = "(i3,a4,2(i3,2x),2f10.3,4f12.4)"
            write(*,*) 'Diagonal kernels:'
            do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
                if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                    parity = 1 ! +
                else
                    parity = 2 ! -
                end if
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
                write(*,*) 'J  Parity  K   q  beta2   beta3       E        n^J(q,q)   <N>     <Z>'
                do K1 = K1_start,K1_end
                    do q1 = gcm_space%q1_start, gcm_space%q1_end
                        write(*,format11) J, ParityChar(parity),K1,q1,constraint%betac(q1), constraint%bet3c(q2),&
                                    dreal(GCM_kernels%H_KK(J,parity,q1,K1,q1,K1)/GCM_kernels%N_KK(J,parity,q1,K1,q1,K1)),&
                                    dreal(GCM_kernels%N_KK(J,parity,q1,K1,q1,K1)), &
                                    dreal(GCM_kernels%X_KK(J,parity,q1,K1,q1,K1,1)), &
                                    dreal(GCM_kernels%X_KK(J,parity,q1,K1,q1,K1,2))
                    end do
                end do
                write(*,*) '--------------------------------------------------------------------'
                write(*,*) 'existing of norm matrix elements -- value 99 means missing, otherwise ln(overlap)'
                length_K = K1_end - K1_start + 1
                do K1 = K1_start,K1_end
                    do q1 = gcm_space%q1_start, gcm_space%q1_end
                        qK1 = q1+(K1-K1_start)*length_K 
                        write(*,'(2i3)',ADVANCE='NO') K1,q1
                        do K2 = K2_start,K2_end
                            do q2 = gcm_space%q2_start, gcm_space%q2_end
                                qK2 = q2+(K2-K2_start)*length_K 
                                if(qK2.le.qK1) then
                                    abs_norm = abs(dreal(GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)))
                                    if(abs_norm .gt. 0.d0) then 
                                        log_abs_norm = int(log(abs_norm))
                                        write(*,'(i3)',ADVANCE='NO') log_abs_norm
                                    else
                                        write(*,'(i3)',ADVANCE='NO') 99
                                    end if 
                                end if 
                            end do 
                        end do 
                        write(*,*) ''
                    end do
                end do
                write(*,*) '--------------------------------------------------------------------'
            end do
        end subroutine
    end subroutine
END MODULE