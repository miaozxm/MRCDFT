!==============================================================================!
! MODULE HWG                                                                   !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE HWG
    implicit none
    contains
    subroutine generator_coordinate_basis
        use Globals, only: constraint,gcm_space,Proj_option,GCM_basis
        integer :: J, parity,K1_start,K1_end,K2_start,K2_end,length_K,length_q,idx,K,q
        allocate(GCM_basis%basis(2,constraint%length*(2*gcm_space%Jmax+1),0:gcm_space%Jmax,2),source=999)
        
        GCM_basis%N_max = 0
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            ! parity
            if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                parity = 1 ! +
            else
                parity = 2 ! -
            end if
            ! K 
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
            length_K = K1_end - K1_start + 1
            ! q
            if(gcm_space%q1_start/=gcm_space%q2_start ) stop "q1_start must be the same as q2_start"
            if(gcm_space%q1_end/=gcm_space%q2_end ) stop "q1_end must be the same as q2_end"
            length_q = gcm_space%q1_end - gcm_space%q1_start + 1

            GCM_basis%N(J,parity) = length_q*length_K
            GCM_basis%N_max = max(GCM_basis%N_max, GCM_basis%N(J,parity))
            idx = 0
            do K = K1_start, K1_end
                do q = gcm_space%q1_start, gcm_space%q1_end
                    idx = idx + 1
                    GCM_basis%basis(1,idx,J,parity) = q ! q
                    GCM_basis%basis(2,idx,J,parity) = K ! K 
                end do 
            end do
            if(idx/=GCM_basis%N(J,parity)) stop 'Incorrect basis length!'
        end do 
    end subroutine

    subroutine solve_HWG_equation
        !--------------------------------------------------------------------
        ! Solve the Hill-Wheeler-Griffin equation :
        !         H * f = E * N * f
        ! where N is the norm (overlap) kernel and H is the Hamiltonian kernel
        !--------------------------------------------------------------------
        use Constants, only: r64
        use Globals, only: gcm_space,GCM_basis,Proj_option,GCM_kernels,GCM_HWG,GCM_outputfile
        use MathMethods, only: EIGENSOLVER_GEP
        real(r64), dimension(:,:), allocatable :: NN,HH,DD,GG,FF,WW,RR
        real(r64), dimension(:), allocatable :: E, EN
        real(r64) :: EPS
        integer :: J,parity,N,qK1,qK2,q1,K1,q2,K2,M,IS,IFL
        
        allocate(GCM_HWG%M(0:gcm_space%Jmax,2),GCM_HWG%E(1:GCM_basis%N_max,0:gcm_space%Jmax,2),&
                 GCM_HWG%fJKq(1:GCM_basis%N_max,1:GCM_basis%N_max,0:gcm_space%Jmax,2),&
                 GCM_HWG%gJKq(1:GCM_basis%N_max,1:GCM_basis%N_max,0:gcm_space%Jmax,2))
        
        write(*,"(5x,A)") 'Solving the HWG equation for each J^pi...'
        open(GCM_outputfile%u_outGCM_HWG ,form='formatted',file=GCM_outputfile%outGCM_HWG)

        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            ! parity
            if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                parity = 1 ! +
            else
                parity = 2 ! -
            end if

            N = GCM_basis%N(J,parity)
            ! store Norm matrix and Hamiltonian kernel
            allocate(NN(N,N),HH(N,N),source=0.d0)
            
            do qK1 = 1, N
                do qK2 = 1, N
                    q1 = GCM_basis%basis(1,qK1,J,parity)
                    K1 = GCM_basis%basis(2,qK1,J,parity)
                    q2 = GCM_basis%basis(1,qK2,J,parity)
                    K2 = GCM_basis%basis(2,qK2,J,parity)
                    NN(qK1,qK2) = dreal(GCM_kernels%N_KK(J,parity,q1,K1,q2,K2))
                    HH(qK1,qK2) = dreal(GCM_kernels%H_KK(J,parity,q1,K1,q2,K2))
                end do 
            end do
            call print_NN_HH

            ! solution of a generalized eigenvalue problem: HH*FF = E*NN*FF
            ! 1. Diagonalize NN: 
            !    NN * DD = EN * DD                         
            ! 2. Discard eigenpairs with EN/EN(1) < EPS                    
            ! 3. Construct reduced Hamiltonian:                            
            !    H_reduced = (DD*EN^{-1/2})^+ * HH * (DD*EN^{-1/2})
            ! 4. Diagonalize H_reduced to obtain eigenvalues E and GG      
            !    H_reduced * GG = E GG
            ! 5. Recover FF = (DD * EN^{-1/2}) * GG (IS>=1)     
            ! 6. Compute WW = NN * FF (IS>=2)                       
            ! 7. Compute orthogonal wave functions RR (IS>=3)       
            !    where RR = NN^{1/2}*FF = DD*EN^{1/2}*DD^+ FF  = DD * GG
            !    and satisfies RR^T * RR = I (identity)
            allocate(EN(N),E(N),DD(N,N),GG(N,N),FF(N,N),WW(N,N),RR(N,N),source=0.d0)
            ! Call the generalized eigenvalue solver
            ! Input:  NN (norm matrix), HH (Hamiltonian matrix), N (matrix size), EPS (cutoff)
            ! Output: EN (norm eigenvalues), DD (norm eigenvectors), 
            !         M (number of retained states after cutoff), E (energy eigenvalues),
            !         GG (reduced eigenvectors), FF (physical eigenvectors in original basis),
            !         WW (NN * FF), RR (collective wave function = NN^(1/2) * FF),
            !         IFL (error flag: 0 = success, 1=no eigenvalue retained )
            EPS = GCM_HWG%cutoff(J)
            IS = 3
            call EIGENSOLVER_GEP(N,N,NN,HH,EN,DD,M,E,GG,FF,WW,RR,EPS,IS,IFL)
            if(IFL/=0) stop 'No eigenvalue retained.'
            ! After discarding, the useful matrices and their dimensions are:
            ! EN(M), E(M)  
            ! DD: NxM , GG: MxM, FF: NxM, WW: NxM, RR: NxM

            ! store
            GCM_HWG%M   (        J,parity) = M
            GCM_HWG%E   (    1:M,J,parity) = E(1:M)      ! E^{J parity}
            GCM_HWG%fJKq(1:N,1:M,J,parity) = FF(1:N,1:M) ! f^{J parity}(K,q)
            GCM_HWG%gJKq(1:N,1:M,J,parity) = RR(1:N,1:M) ! g^{J parity}(K,q), collective wave function
            call print_HWG_Results
            deallocate(NN,HH,EN,E,DD,GG,FF,WW,RR)
        end do 
        close(GCM_outputfile%u_outGCM_HWG)
        contains

        subroutine print_NN_HH
            use Tools, only: int2str
            integer :: max_N,N,qK2,q2,K2,K1,q1,qK1
            character(len=*), parameter ::  format0 = "(3x,'[',I3,',',I3,']')", &
                                            format1 = "(a,i3,a,a)", &
                                            format2 = "('[',i3,',',i3,']',50f12.5)"
            character(1), dimension(2) :: ParityChar = ['+', '-']
            character(len=10000) :: header
            max_N = 50
            N = GCM_basis%N(J,parity)

            ! Create table headers
            header = ''
            do qK2 = 1, min(N,max_N)
                q2 = GCM_basis%basis(1,qK2,J,parity)
                K2 = GCM_basis%basis(2,qK2,J,parity)
                write(header(len_trim(header)+1:),format0) K2,q2
            end do
            if(N>max_N) then
                write(header(len_trim(header)+1:),*) '...'
            end if
            write(GCM_outputfile%u_outGCM_HWG,format1) '=================  J^pi:',J,ParityChar(parity),' ================='
            ! NN 
            write(GCM_outputfile%u_outGCM_HWG,"(A,3x,A)") 'NN( [K1,q1], [K2,q2] )', int2str(N)//'x'//int2str(N)
            write(GCM_outputfile%u_outGCM_HWG,"(10x,A)") trim(header)
            do qK1 = 1, min(N,max_N)
                q1 = GCM_basis%basis(1,qK1,J,parity)
                K1 = GCM_basis%basis(2,qK1,J,parity)
                write(GCM_outputfile%u_outGCM_HWG,format2) K1,q1,(NN(qK1,qK2), qK2=1,min(N,max_N)) 
            end do
            if(N>max_N) then
                write(GCM_outputfile%u_outGCM_HWG,*) '...'
            end if
            write(GCM_outputfile%u_outGCM_HWG,*) 
            ! HH
            write(GCM_outputfile%u_outGCM_HWG,"(A,3x,A)") 'HH( [K1,q1], [K2,q2] )',int2str(N)//'x'//int2str(N)
            write(GCM_outputfile%u_outGCM_HWG,"(10x,A)") trim(header)
            do qK1 = 1, min(N,max_N)
                q1 = GCM_basis%basis(1,qK1,J,parity)
                K1 = GCM_basis%basis(2,qK1,J,parity)
                write(GCM_outputfile%u_outGCM_HWG,format2) K1,q1,(HH(qK1,qK2), qK2=1,min(N,max_N)) 
            end do
            if(N>max_N) then
                write(GCM_outputfile%u_outGCM_HWG,*) '...'
            end if
            write(GCM_outputfile%u_outGCM_HWG,*) 
        end subroutine

        subroutine print_HWG_Results
            use Globals, only: GCM_basis,constraint
            use Tools, only: int2str
            integer :: max_M,i,K,q,qK
            character(len=500) :: header
            character(1), dimension(2) :: ParityChar = ['+', '-']
            character(len=*), parameter ::  format1 = "(50e11.3)", &
                                            format2 = '(2i3,1x,2f8.3,50(1x,f10.4))'

            max_M = 50
 
            write(GCM_outputfile%u_outGCM_HWG, '(A,e10.2)') ' Cutoff value of norm eigenvalues: ', EPS
            write(GCM_outputfile%u_outGCM_HWG,*) 'Number of retained states: ', M
            write(GCM_outputfile%u_outGCM_HWG,*)   'J^pi_i   eigenvalues_norm    eigenvalues_Hamiltonian'
            ! write(*, '(A)') '--------------------------------------------------' 
            do i = 1, M
                write(GCM_outputfile%u_outGCM_HWG, '(2x,A, 5X, E15.6, 5X, E15.6)')int2str(J)//ParityChar(parity)//'_'//int2str(i), EN(i), E(i)
            end do
            write(GCM_outputfile%u_outGCM_HWG,*) 
            ! f
            write(GCM_outputfile%u_outGCM_HWG,*) 'Eigenvectors f^{J pi}(K,q)'
            write(header, '(A, *(A11))') '  K  q    beta2   beta3', (trim(int2str(J)//ParityChar(parity)//'_'//int2str(i)), i=1, min(M, max_M))
            if(M>max_M) then
                write(header(len_trim(header)+1:),*) '...'
            end if
            write(GCM_outputfile%u_outGCM_HWG, '(A)') trim(header)
            do qK = 1, GCM_basis%N(J,parity)
                q = GCM_basis%basis(1,qK,J,parity)
                K = GCM_basis%basis(2,qK,J,parity)
                write(GCM_outputfile%u_outGCM_HWG,format2) K,q,constraint%betac(q), constraint%bet3c(q),(FF(qK,i), i=1,min(M,max_M))
            end do
            if(M>max_M) then
                write(GCM_outputfile%u_outGCM_HWG,*) '...'
            end if
            write(GCM_outputfile%u_outGCM_HWG,*) 
            ! g
            write(GCM_outputfile%u_outGCM_HWG,*) 'collective wave function:'
            write(GCM_outputfile%u_outGCM_HWG,*) 'Eigenvectors  g^{J pi}(K,q)'
            write(header, '(A, *(A11))') '  K  q    beta2   beta3', (trim(int2str(J)//ParityChar(parity)//'_'//int2str(i)), i=1, min(M, max_M))
            if(M>max_M) then
                write(header(len_trim(header)+1:),*) '...'
            end if
            write(GCM_outputfile%u_outGCM_HWG, '(A)') trim(header)
            do qK = 1, GCM_basis%N(J,parity)
                q = GCM_basis%basis(1,qK,J,parity)
                K = GCM_basis%basis(2,qK,J,parity)
                write(GCM_outputfile%u_outGCM_HWG,format2) K,q,constraint%betac(q), constraint%bet3c(q),(RR(qK,i), i=1,min(M,max_M))
            end do
            if(M>max_M) then
                write(GCM_outputfile%u_outGCM_HWG,*) '...'
            end if
            write(GCM_outputfile%u_outGCM_HWG,*) '---------------------------------------------'
            write(GCM_outputfile%u_outGCM_HWG,*) 
            write(GCM_outputfile%u_outGCM_HWG,*) 
        end subroutine
    end subroutine

END MODULE