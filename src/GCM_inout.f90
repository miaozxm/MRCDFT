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
                                        format2 = "(10x,5e10.2)"
        open(u_GCM, file=InputFile%file_path_para, status='old')
        ! skip GCM parameters
        do i = 1, 50
            read(u_GCM,'(A)', iostat=is) 
            if (is /= 0) then
                print *, "Error reading GCM file"
                stop
            end if
        end do
        ! read
        read(u_GCM, format1) input_par%GCMType
        read(u_GCM, format1) input_par%Mmax
        read(u_GCM, format2) input_par%zeta
        close(u_GCM)

        call set_GCM_parameters
        if(ifPrint.and. MPI_Infor%rank == 0) call printParameters

        contains

        subroutine set_GCM_parameters
            use Constants, only: Jmax_max
            use Globals, only: input_par,GCM_option,GCM_HWG
            integer :: J
            ! set GCM type
            GCM_option%GCMType = input_par%GCMType
            if(GCM_option%GCMType > 1 .or. GCM_option%GCMType < 0) stop 'GCMType wrong!'

            ! set Mmax
            GCM_HWG%Mmax = input_par%Mmax
            if(input_par%Mmax < 0) stop 'Mmax cannot be less than 1.' 

            ! set cutoff
            do J = 0, Jmax_max
                if(J<=4) then
                    GCM_HWG%cutoff(J) = input_par%zeta(J+1)
                else 
                    GCM_HWG%cutoff(J) = input_par%zeta(5)
                end if 
            end do

        end subroutine
        subroutine printParameters
            use Globals, only: input_par,GCM_option,GCM_HWG
            use Tools, only: adjust_left
            integer :: Strlength = 40

            if(GCM_option%GCMType==0) then 
                write(*,"(5x,a,':   ',a)") adjust_left('GCM',Strlength),'GCM skipped.'
            else if(GCM_option%GCMType==1) then
                write(*,"(5x,A)") 'GCM:'
            end if 

            write(*,"(5x,a,': ', i3)") adjust_left('Mmax',Strlength),GCM_HWG%Mmax
            write(*,"(5x,a,': ', 5(e9.2))") adjust_left('Zeta',Strlength),input_par%zeta
            write(*,"(a)") '=========================================================================================='
        end subroutine
    end subroutine

    subroutine read_kernels
        use Constants, only: r64
        use Proj_Inout, only: set_Proj_output_filename
        use Globals, only: gcm_space,Proj_outputfile,GCM_kernels,Proj_option,GCM_outputfile
        integer :: q1, q2,J,K1_start,K1_end,K2_start,K2_end,K1,K2,q2_start,q2_end,parity
        integer :: read_J, read_K1, read_K2
        character(1) :: read_parity_char
        real(r64) :: read_N2, read_Z2, read_J2
        complex(r64) :: Q2_KK_21(2), E0_KK_21(2), r2_2b_dummy
        character(1), dimension(2) :: ParityChar = ['+', '-']
        character(len=*), parameter ::  format1 = "(3i5,4x,a,4x,3f12.3)", &
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
        ! r2_2b_KK(J,parity,q1,K1,q2,K2, term, channel)
        allocate(GCM_kernels%r2_2b_KK(0:gcm_space%Jmax,2,gcm_space%q1_start:gcm_space%q1_end,-gcm_space%Jmax:gcm_space%Jmax, &
                                gcm_space%q2_start:gcm_space%q2_end,-gcm_space%Jmax:gcm_space%Jmax,5,3),source=(0.d0,0.d0))

        write(*,"(5x,A)") 'Reading kernels for GCM calculations...'
        do q1 = gcm_space%q1_start, gcm_space%q1_end
            if(Proj_option%Kernel_Symmetry==0) then      ! All Kernels 
                q2_start = gcm_space%q2_start
                q2_end = gcm_space%q2_end
            else if(Proj_option%Kernel_Symmetry==1) then ! Triangular Kernels
                q2_start = q1
                q2_end = gcm_space%q2_end
            else if(Proj_option%Kernel_Symmetry==2) then ! Diagonal elements
                q2_start = q1
                q2_end = q1
                if(q1==gcm_space%q1_start) write(*,*) '!!! Warning: Only diagonal kernels cannot be used for GCM calculations!'
            else
                stop 'Kernel_Symmetry should be 0, 1, 2'
            end if
            do q2 = q2_start, q2_end
                call set_Proj_output_filename(q1,q2)
                write(GCM_outputfile%u_outGCM_standard,"(A,i3,A,i3,A,A)") 'Reading kernels for q1=',q1,' and q2=',q2,' from file: ', Proj_outputfile%outputelem
                open(Proj_outputfile%u_outputelem,file=Proj_outputfile%outputelem,form='formatted',status='old')
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
                            ! r²_2b 五项 (期望值格式 = ÷N × 4π/3，读取后 ×N 恢复)
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,1,2), &
                                                                         GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,1,1)    ! O1: 1B,p ; O2: 1B,n
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,2), &
                                                                         GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,3)    ! O3: 2B,pp ; O4: 2B,pn
                            read(Proj_outputfile%u_outputelem,format2) GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,1), &
                                                                         r2_2b_dummy                                          ! O5: 2B,nn ; padding


                            GCM_kernels%H_KK(J,parity,q1,K1,q2,K2) = GCM_kernels%H_KK(J,parity,q1,K1,q2,K2)*GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                            GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,1) = GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,1)*GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                            GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,2) = GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,2)*GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                            ! r²_2b ×N 恢复原始 kernel（只乘实际用到的5项）
                            GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,1,2) = GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,1,2) * GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)  ! O1: 1B,p
                            GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,1,1) = GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,1,1) * GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)  ! O2: 1B,n
                            GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,2) = GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,2) * GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)  ! O3: 2B,pp
                            GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,3) = GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,3) * GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)  ! O4: 2B,pn
                            GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,1) = GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,2,1) * GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)  ! O5: 2B,nn
                            if(Proj_option%Kernel_Symmetry==1 .and. q1/=q2) then
                                GCM_kernels%N_KK(J,parity,q2,K2,q1,K1) = GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)
                                GCM_kernels%H_KK(J,parity,q2,K2,q1,K1) = GCM_kernels%H_KK(J,parity,q1,K1,q2,K2)
                                GCM_kernels%X_KK(J,parity,q2,K2,q1,K1,:) = GCM_kernels%X_KK(J,parity,q1,K1,q2,K2,:)
                                GCM_kernels%Q2_KK(J,parity,q2,K2,q1,K1,:) = Q2_KK_21(:)
                                GCM_kernels%E0_KK(J,parity,q2,K2,q1,K1,:) = E0_KK_21(:)
                                GCM_kernels%r2_2b_KK(J,parity,q2,K2,q1,K1,:,:) = GCM_kernels%r2_2b_KK(J,parity,q1,K1,q2,K2,:,:)
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
            character(len=*), parameter ::  format11 = "(i3,a6,2x,(i3,1x),(i3,2x),2(f5.3,2x),(f7.4,2x),3(f12.4,2x))"
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
                
                write(GCM_outputfile%u_outGCM_standard,*) 'Diagonal kernels:'
                write(GCM_outputfile%u_outGCM_standard,*) ' J  Parity  K   q  beta2  beta3   n^J(q,q)       E            <N>           <Z>'
                do K1 = K1_start,K1_end
                    do q1 = gcm_space%q1_start, gcm_space%q1_end
                        write(GCM_outputfile%u_outGCM_standard,format11) J, ParityChar(parity),K1,q1,constraint%betac(q1), constraint%bet3c(q1),dreal(GCM_kernels%N_KK(J,parity,q1,K1,q1,K1)), &
                                        dreal(GCM_kernels%H_KK(J,parity,q1,K1,q1,K1)/(GCM_kernels%N_KK(J,parity,q1,K1,q1,K1)+1.0d-30)),&
                                        dreal(GCM_kernels%X_KK(J,parity,q1,K1,q1,K1,1)/(GCM_kernels%N_KK(J,parity,q1,K1,q1,K1)+1.0d-30)), &
                                        dreal(GCM_kernels%X_KK(J,parity,q1,K1,q1,K1,2)/(GCM_kernels%N_KK(J,parity,q1,K1,q1,K1)+1.0d-30))
                    end do
                end do  
                write(GCM_outputfile%u_outGCM_standard,*) 'Triangular norm kernels:  ln(overlap), value 99 means missing:'
                write(GCM_outputfile%u_outGCM_standard,*) ' K','  q'
                length_K = K1_end - K1_start + 1
                do K1 = K1_start,K1_end
                    do q1 = gcm_space%q1_start, gcm_space%q1_end
                        qK1 = q1+(K1-K1_start)*length_K 
                        write(GCM_outputfile%u_outGCM_standard,'(2i3,2x)',ADVANCE='NO') K1,q1
                        do K2 = K2_start,K2_end
                            do q2 = gcm_space%q2_start, gcm_space%q2_end
                                qK2 = q2+(K2-K2_start)*length_K 
                                if(qK2.le.qK1) then
                                    abs_norm = abs(dreal(GCM_kernels%N_KK(J,parity,q1,K1,q2,K2)))
                                    if(abs_norm .gt. 0.d0) then 
                                        log_abs_norm = int(log(abs_norm))
                                        write(GCM_outputfile%u_outGCM_standard,'(i4)',ADVANCE='NO') log_abs_norm
                                    else
                                        write(GCM_outputfile%u_outGCM_standard,'(i4)',ADVANCE='NO') 99
                                    end if 
                                end if 
                            end do 
                        end do 
                        write(GCM_outputfile%u_outGCM_standard,*)
                    end do
                end do
                write(GCM_outputfile%u_outGCM_standard,*) '-------------------------------------------------------------------------------'
                write(GCM_outputfile%u_outGCM_standard,*)
            end do
        end subroutine
    end subroutine

    subroutine set_GCM_output_filename
        use Globals, only: OUTPUT_PATH,r64,i16
        use Globals, only: BS,projection_mesh,nucleus_attributes,Proj_option,GCM_outputfile
        use Kernel, only: set_projection_mesh_points
        use Basis, only: set_Spherical_HO_basis_parameters
        use Tools, only: int2str
        integer :: q1,q2,AMPType,A
        integer(i16) :: AMP,name_nf1,name_nf2,nphi_1,nphi_2,nbeta_1,nbeta_2
        if(Proj_option%ProjectionType == 0 ) then 
            call set_projection_mesh_points ! set projection_mesh%nphi, projection_mesh%nbeta
            call set_Spherical_HO_basis_parameters(.False.) ! set BS%HO_sph%n0f
        end if 
        A = nucleus_attributes%mass_number_int
        !------
        AMPType = Proj_option%AMPType
        if(AMPType==0) AMP = 0 + 48
        if(AMPType==1) AMP = 1 + 48
        if(AMPType==3) AMP = 3 + 48
        !-----
        name_nf1 = mod(BS%HO_sph%n0f/10,10) + 48
        name_nf2 = mod(BS%HO_sph%n0f,10) + 48
        !-----
        nphi_1 = mod(projection_mesh%nphi(1)/10,10) + 48
        nphi_2 = mod(projection_mesh%nphi(1),10) + 48
        nbeta_1 = mod(projection_mesh%nbeta/10,10) + 48
        nbeta_2 = mod(projection_mesh%nbeta,10) + 48


        GCM_outputfile%outGCM_standard = OUTPUT_PATH//'GCM_'//int2str(A)//nucleus_attributes%name &
                            //'.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                            //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2)//'.out'
        GCM_outputfile%outGCM_HWG  = OUTPUT_PATH//'GCM_'//int2str(A)//nucleus_attributes%name &
                            //'_hwg'//'.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                            //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2)//'.dat'
        GCM_outputfile%outGCM_observables = OUTPUT_PATH//'GCM_'//int2str(A)//nucleus_attributes%name &
                            //'_obs'//'.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                            //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2)//'.dat'                     
        
    end subroutine

    subroutine write_GCM_observables
        use Globals, only: GCM_outputfile
        open(GCM_outputfile%u_outGCM_observables ,form='formatted',file=GCM_outputfile%outGCM_observables)
        write(*,"(5x,A)") 'Writing GCM observables to file: '//GCM_outputfile%outGCM_observables
        call print_spectrum
        close(GCM_outputfile%u_outGCM_observables)

        contains
        subroutine print_spectrum
            use Globals, only: gcm_space,Proj_option,GCM_HWG,GCM_obser
            integer :: J,parity,iM
            character(1), dimension(2) :: ParityChar = ['+', '-']
            character(len=*), parameter ::  format1 = "(i2,a1,i2,f12.4,3f10.4,f12.4,3f9.4,2f9.4)"
            write(GCM_outputfile%u_outGCM_observables,'(A,f10.4,A)') 'Ground State Energy:', GCM_HWG%E(1,0,1), 'MeV'
            write(GCM_outputfile%u_outGCM_observables,*) '------------------------------------------------------'
            write(GCM_outputfile%u_outGCM_observables,"(A)") 'Excitation Spectrum:'
            write(GCM_outputfile%u_outGCM_observables,"(A)") 'J^pi_i        E       E_ex    <beta2>   <beta3>       <N>       <Z>    rrms_p'
            do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
                ! parity
                if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                    parity = 1 ! +
                else
                    parity = 2 ! -
                end if
                do iM = 1, GCM_HWG%M(J,parity)
                    write(GCM_outputfile%u_outGCM_observables,format1) J,ParityChar(parity),iM,GCM_HWG%E(iM,J,parity),GCM_obser%E_ex(iM,J,parity),&
                                     GCM_obser%beta2_aver(iM,J,parity),GCM_obser%beta3_aver(iM,J,parity),&
                                     GCM_obser%N(iM,J,parity),GCM_obser%Z(iM,J,parity),GCM_obser%rrms_p(iM,J,parity)
                end do 
            end do
            write(GCM_outputfile%u_outGCM_observables,*) '------------------------------------------------------'
        end subroutine
    end subroutine
END MODULE