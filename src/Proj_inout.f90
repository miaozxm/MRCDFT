!==============================================================================!
! MODULE Proj_Inout                                                            !
!                                                                              !
! This module contains functions and routines for reading and writing files.   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Proj_Inout
use Globals, only: Proj_outputfile
use Constants, only: i16,r64,u_start,pi,ngl,OUTPUT_PATH, Jmax_max
use CDFT_Inout, only: set_CDFT_output_filename
use Tools, only: int2str, adjust_left
implicit none
integer, private :: u_Proj = u_start + 11
logical :: first_kernel = .True.

contains

subroutine read_Proj_configuration(ifPrint)
    use Globals, only: InputFile,input_par,Proj_option,MPI_Infor
    integer :: i,is
    logical,intent(in),optional :: ifPrint
    character(len=*), parameter ::  format1 = "(10x,2f10.4)", &
                                    format2 = "(10x,i5)", &
                                    format3 = "(10x,2i5)" 

    open(u_Proj, file=InputFile%file_path_para, status='old')
    ! skip CDFT parameters
    do i = 1, 30
        read(u_Proj,'(A)', iostat=is) 
        if (is /= 0) then
            print *, "Error reading file"
            stop
        end if
    end do

    read(u_Proj, format2) input_par%ProjectionType
    read(u_Proj, format2) input_par%AMPType
    read(u_Proj, format2) input_par%PNPType
    read(u_Proj, format2) input_par%PPType
    read(u_Proj, format2) input_par%Kernel_Symmetry
    read(u_Proj, format3) input_par%q1_start, input_par%q1_end
    read(u_Proj, format3) input_par%q2_start, input_par%q2_end
    read(u_Proj, format2) input_par%Jmax
    read(u_Proj, format2) input_par%icm
    read(u_Proj, format2) input_par%nphi
    read(u_Proj, format2) input_par%Euler_Symmetry
    read(u_Proj, format2) input_par%nalpha
    read(u_Proj, format2) input_par%nbeta
    read(u_Proj, format2) input_par%ngamma
    read(u_Proj, format2) input_par%DsType
    read(u_Proj, format2) input_par%TDType
    read(u_Proj, format2) input_par%lambda_max
    read(u_Proj, format2) input_par%checkN2J2
    read(u_Proj, format2) input_par%EccentriType
    close(u_Proj)
    call set_Proj_parameters
    if(ifPrint.and. MPI_Infor%rank == 0) call printParameters
    contains
    subroutine set_Proj_parameters
        use Globals, only: input_par,Proj_option,gcm_space,constraint,TDs
        ! set projection type ( 0 : skip Proj;  1: do Proj)
        Proj_option%ProjectionType = input_par%ProjectionType
        if(Proj_option%ProjectionType > 1 .or. Proj_option%ProjectionType < 0) stop 'ProjectionType wrong!'

        ! set AMP type ((0) no (1) 1DAMP (2) 3DAMP)
        Proj_option%AMPtype = input_par%AMPType
        if(input_par%AMPType.ne.0 .and. input_par%AMPType.ne.1 .and. input_par%AMPType.ne.2) stop 'AMPType wrong!'
        if( input_par%AMPType==2) stop '(3DAMP) Not yet verified !'
        if(input_par%AMPtype==1 .and. mod(input_par%nbeta, 2) /= 0) stop 'nbeta must be an even number!' 

        ! set PNP type (0: no PNP; 1: PNP)
        Proj_option%PNPtype = input_par%PNPType
        if(input_par%PNPType.ne.0 .and. input_par%PNPType.ne.1 ) stop 'PNPType wrong!'
        if(input_par%PNPType.ne.0 .and. mod(input_par%nphi, 2) == 0) stop 'nphi must be an odd number!' 

        ! set PP type (0: no PP; 1: PP)
        Proj_option%PPtype = input_par%PPType
        if(input_par%PPType.ne.0 .and. input_par%PPType.ne.1 ) stop 'PNPType wrong!'

        ! Kernel Symmetry  (0: All ; 1: Triangular Matrix ; 2: Diagonal elements only)
        Proj_option%Kernel_Symmetry = input_par%Kernel_Symmetry
        if(Proj_option%Kernel_Symmetry > 2 .or. Proj_option%Kernel_Symmetry < 0) stop 'Kernel_Symmetry wrong!'

        ! Euler angles Symmetry (0: no, 1: Axially, 2: D2)
        Proj_option%Euler_Symmetry = input_par%Euler_Symmetry
        if(Proj_option%Euler_Symmetry .ne. 0 .and. Proj_option%Euler_Symmetry .ne. 1 .and. Proj_option%Euler_Symmetry .ne. 2) stop 'Euler_Symmetry wrong!'
        if(input_par%AMPType == 1  .and. (Proj_option%Euler_Symmetry .ne. 0 .and. Proj_option%Euler_Symmetry .ne. 1)) stop '1DAMP: Euler_Symmetry should be 0 or 1!'
        if(input_par%AMPType == 2  .and. (Proj_option%Euler_Symmetry .ne. 2)) stop '3DAMP: Euler_Symmetry should 2 !'

        ! center of mass correction (1: average  ; 2: HO approximation)
        Proj_option%icm = input_par%icm
        ! Set Norm overlap calculation method
        Proj_option%ihf = 3

        ! set GCM space
        gcm_space%Jmin = 0
        gcm_space%Jmax = input_par%Jmax
        if(gcm_space%Jmax > Jmax_max) stop 'Jmax_max too small!'
        gcm_space%Jstep = 1 
        gcm_space%q1_start = input_par%q1_start
        gcm_space%q1_end = input_par%q1_end
        if(gcm_space%q1_end == -1) gcm_space%q1_end = constraint%length
        gcm_space%q2_start = input_par%q2_start
        gcm_space%q2_end = input_par%q2_end
        if(gcm_space%q2_end == -1) gcm_space%q2_end = constraint%length
        if(gcm_space%q1_start < 1) stop 'q1_start less than 1 !'
        if(gcm_space%q1_end < gcm_space%q1_start) stop 'q1_end less than q1_start !'
        if(gcm_space%q1_end > constraint%length) stop 'q1_end greater than max length.'
        if(gcm_space%q2_start < 1) stop 'q2_start less than 1 !'
        if(gcm_space%q2_end < gcm_space%q2_start) stop 'q2_end less than q2_start !'
        if(gcm_space%q2_end > constraint%length) stop 'q2_end greater than max length.'

        ! set calculate density matrix element option
        Proj_option%DsType = input_par%DsType
        if(Proj_option%DsType<0 .or. Proj_option%DsType>3) stop "DsType should be 0 or 1 or 2 or 3"
        ! set calculate transition density matrix element option
        Proj_option%TDType = input_par%TDType
        if(Proj_option%TDType<0 .or. Proj_option%TDType>1) stop "TDType should be 0 or 1"
        TDs%lambda_max = input_par%lambda_max ! max lambda of 1B transition density
        ! set calculate N^2 and J^2 option
        Proj_option%checkN2J2 = input_par%checkN2J2
        if(Proj_option%checkN2J2<0 .or. Proj_option%checkN2J2>1) stop "checkN2J2 should be 0 or 1"
        ! set calculate eccentricity kernel option
        Proj_option%EccentriType = input_par%EccentriType
        if(Proj_option%EccentriType<0 .or. Proj_option%EccentriType>3) stop "EccentriType should be 0 or 1 or 2 or 3"
    end subroutine
    subroutine printParameters
        use Globals, only: input_par,Proj_option,gcm_space,TDs
        integer :: Strlength = 40
        character(len=5) :: AMP_char, PNP_char,PP_char
        if(Proj_option%AMPtype==0) then 
            AMP_char = 'noAMP'
        else if (Proj_option%AMPtype==1) then
            AMP_char = '1DAMP'
        else if (Proj_option%AMPtype==2) then
            AMP_char = '3DAMP'
        end if 
        if(Proj_option%PNPtype==0) then 
            PNP_char = 'noPNP'
        else if (Proj_option%PNPtype==1) then
            PNP_char = 'PNP'
        end if 
        if(Proj_option%PPtype==0) then 
            PP_char = 'noPP'
        else if (Proj_option%PNPtype==1) then
            PP_char = 'PP'
        end if 
        if(Proj_option%ProjectionType == 0) then
            write(*,'(5x,A)') 'Projection skipped.'
        else 
            write(*,'(5x,A)') AMP_char//'  +  '//PNP_char//'+  '//PP_char//':'
        end if 
        if(Proj_option%AMPtype /= 0) then
            write(*,"(5x,a,':   ',3(i2,a))") adjust_left('Number of euler angles',Strlength),input_par%nalpha,' (nalpha),  ',input_par%nbeta,' (nbeta),  ',input_par%ngamma,' (ngamma)'
            if(input_par%Euler_Symmetry==0) then
                write(*,"(5x,a,':   ',a)") adjust_left('Symmetry of euler angles',Strlength),'no'
            else if(input_par%Euler_Symmetry==1) then
                write(*,"(5x,a,':   ',a)") adjust_left('Symmetry of euler angles',Strlength),'Axially'
            else if(input_par%Euler_Symmetry==2) then
                write(*,"(5x,a,':   ',a)") adjust_left('Symmetry of euler angles',Strlength),'D2'
            end if 
        end if 
        if(Proj_option%PNPtype /= 0) then
            write(*,"(5x,a,':   ',i2,a)") adjust_left('Number gauge angles',Strlength),input_par%nphi,' (nphi)'
        end if

        if(Proj_option%Kernel_Symmetry==0) then 
            write(*,"(5x,a,':   ',a)") adjust_left('Kernels',Strlength),'All kernels'
        else if(Proj_option%Kernel_Symmetry==1) then
            write(*,"(5x,a,':   ',a)") adjust_left('Kernels',Strlength),'Upper triangular kernels'
        else if(Proj_option%Kernel_Symmetry==2) then
            write(*,"(5x,a,':   ',a)") adjust_left('Kernels',Strlength),'Diagonal kernels'
        end if 
        write(*,"(5x,a,':   ',a,i3,a,i3,a)") adjust_left('Quadratic constraint q1 range',Strlength), '[',gcm_space%q1_start,',',gcm_space%q1_end,' ]'
        write(*,"(5x,a,':   ',a,i3,a,i3,a)") adjust_left('Quadratic constraint q2 range',Strlength), '[',gcm_space%q2_start,',',gcm_space%q2_end,' ]'

        write(*,"(5x,a,':   ',i3)") adjust_left('Maximal J value',Strlength), gcm_space%Jmax

        if (Proj_option%checkN2J2 == 1) then
            write(*,"(5x,a,':   ',a)") adjust_left('Check N2 and J2',Strlength),'Yes'
        else 
            write(*,"(5x,a,':   ',a)") adjust_left('Check N2 and J2',Strlength),'No'
        end if

        if (Proj_option%DsType /= 0 ) then
            if(Proj_option%DsType==1) write(*,"(5x,a,':   ',a)") adjust_left('Calculate density ME',Strlength),'1B'
            if(Proj_option%DsType==2) write(*,"(5x,a,':   ',a)") adjust_left('Calculate density ME',Strlength),'2B'
            if(Proj_option%DsType==3) write(*,"(5x,a,':   ',a)") adjust_left('Calculate density ME',Strlength),'1B + 2B'
        else 
            write(*,"(5x,a,':   ',a)") adjust_left('Calculate density ME',Strlength),'No'
        end if 

        if (Proj_option%AMPtype > 0 .and. Proj_option%TDType /=0 ) then
            if(Proj_option%TDType==1)write(*,"(5x,a,':   ',a)") adjust_left('Calculate reduced transition density ME',Strlength),'1B'
            write(*,"(5x,a,':   ',i3)") adjust_left('Maximal lambda value(1BTD)',Strlength),TDs%lambda_max
        else 
            write(*,"(5x,a,':   ',a)") adjust_left('Calculate reduced transition density ME',Strlength),'No'
        end if 

        if (Proj_option%EccentriType /= 0) then
            if(Proj_option%EccentriType==1) write(*,"(5x,a,':   ',a)") adjust_left('Calculate eccentricity kernel',Strlength),'(1)'
            if(Proj_option%EccentriType==2) write(*,"(5x,a,':   ',a)") adjust_left('Calculate eccentricity kernel',Strlength),'(2)'
            if(Proj_option%EccentriType==3) write(*,"(5x,a,':   ',a)") adjust_left('Calculate eccentricity kernel',Strlength),'(3)'
        else 
            write(*,"(5x,a,':   ',a)") adjust_left('Calculate eccentricity kernel',Strlength),'No'
        end if 

        ! if(Proj_option%ihf == 1) then 
        !     write(*,"(5x,a,':   ',a)") adjust_left('Norm overlap formula',Strlength),'sqrt(det(D) det(R))'
        ! else if(Proj_option%ihf == 2) then
        !     write(*,"(5x,a,':   ',a)") adjust_left('Norm overlap formula',Strlength),'Robledo (2009) formula'
        ! else if(Proj_option%ihf == 3) then
        !     write(*,"(5x,a,':   ',a)") adjust_left('Norm overlap formula',Strlength),'Bertsch & Robledo (2011) formula'
        ! end if 
        write(*,*)
    end subroutine
end subroutine

subroutine read_wavefuntion_files(q1,q2)
    use Constants, only: nb_max 
    use Globals, only: wf1,wf2,constraint,outputfile
    integer :: q1,q2,nb,it,temp,ib
    integer, dimension(nb_max,2) :: kd
    call set_CDFT_output_filename(constraint%betac(q1),constraint%bet3c(q1))
    open(outputfile%u_outputwf,file=outputfile%outputwf,form='unformatted',status='unknown')      
    read(outputfile%u_outputwf) !dirac%ka
    read(outputfile%u_outputwf) kd
    read(outputfile%u_outputwf) wf1%v2
    read(outputfile%u_outputwf) wf1%v2d
    read(outputfile%u_outputwf) wf1%skk
    read(outputfile%u_outputwf) wf1%fg  
    read(outputfile%u_outputwf) nb
    read(outputfile%u_outputwf) wf1%ecm 
    read(outputfile%u_outputwf) !pairing%ibk
    close(outputfile%u_outputwf)
    do it=1,2
        temp = 0
        do ib=1,nb
            temp = temp + kd(ib,it)
        enddo
        wf1%nk(it) = temp
    enddo 
    call set_CDFT_output_filename(constraint%betac(q2),constraint%bet3c(q2))
    open(outputfile%u_outputwf,file=outputfile%outputwf,form='unformatted',status='unknown')      
    read(outputfile%u_outputwf) !dirac%ka
    read(outputfile%u_outputwf) kd
    read(outputfile%u_outputwf) wf2%v2
    read(outputfile%u_outputwf) wf2%v2d
    read(outputfile%u_outputwf) wf2%skk
    read(outputfile%u_outputwf) wf2%fg  
    read(outputfile%u_outputwf) nb
    read(outputfile%u_outputwf) wf2%ecm 
    read(outputfile%u_outputwf) !pairing%ibk
    close(outputfile%u_outputwf)
    do it=1,2
        temp = 0
        do ib=1,nb
            temp = temp + kd(ib,it)
        enddo
        wf2%nk(it) = temp
    enddo
end subroutine

subroutine set_Proj_Expectation_filename
    use Globals, only:BS,projection_mesh,nucleus_attributes,Proj_option
    integer :: AMPType,A
    integer(i16) :: AMP,name_nf1,name_nf2,nphi_1,nphi_2,nbeta_1,nbeta_2
    A = nucleus_attributes%mass_number_int
    !------
    AMPType = Proj_option%AMPType
    if(AMPType==0) AMP = 0 + 48
    if(AMPType==1) AMP = 1 + 48
    if(AMPType==3) AMP = 3 + 48
    !-----Nf
    name_nf1 = mod(BS%HO_sph%n0f/10,10) + 48
    name_nf2 = mod(BS%HO_sph%n0f,10) + 48
    !-----nphi,nbeta
    nphi_1 = mod(projection_mesh%nphi(1)/10,10) + 48
    nphi_2 = mod(projection_mesh%nphi(1),10) + 48
    nbeta_1 = mod(projection_mesh%nbeta/10,10) + 48
    nbeta_2 = mod(projection_mesh%nbeta,10) + 48
    Proj_outputfile%outExpectation = OUTPUT_PATH//'Proj_'//int2str(A)//nucleus_attributes%name &
                        //'_'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2) &
                        //'.'//char(nbeta_1)//char(nbeta_2)//'_Expectation.out'
end subroutine

subroutine write_Proj_output(q1,q2)
    use Globals, only: gcm_space,Proj_option
    use Eccentricity, only: calculate_Eccentricity_kernel_by_density_matrix_element
    integer,intent(in) :: q1,q2

    call set_Proj_output_filename(q1,q2)

    ! write kernels
    call write_kernels
    if(Proj_option%EccentriType/=0)then 
        if(Proj_option%EccentriType==2 .or. Proj_option%EccentriType==3) then
            call calculate_Eccentricity_kernel_by_density_matrix_element
        end if 
        call write_eccentricity_operators_kernels(q1,q2)
    end if 
    ! write r2 2body operators kernels
    call write_r2_2body_operators_kernels(q1,q2)
    ! write matrix elemets of operator
    if(q1== gcm_space%q1_start .and. q2==gcm_space%q2_start) then 
        call write_reduced_1B_multipole_matrix_elements
        call write_1B_operators_matrix_elements
    end if 

    ! write density and transition_density
    if (Proj_option%DsType > 0 )then
        if (Proj_option%DsType==1 .or. Proj_option%DsType==3) call write_1B_density_matrix_elements
    end if 
    if (Proj_option%TDType==1) then
        call write_reduced_1B_transition_density_matrix_elements(q1,q2)
    end if 

    ! write expectation
    if(q1==q2) then 
        call write_Proj_expectation(q1,q2)
    end if 
end subroutine

subroutine set_Proj_output_filename(q1,q2)
    use Globals, only:constraint,BS,projection_mesh,nucleus_attributes,Proj_option
    integer :: q1,q2,AMPType,A
    real(r64) :: beta2_1, beta3_1,beta2_2,beta3_2,abs2c1,abs3c1,abs2c2,abs3c2
    character :: signb21,signb31,signb22,signb32
    integer(i16), dimension(6) :: name1,name2
    integer(i16) :: AMP,name_nf1,name_nf2,nphi_1,nphi_2,nbeta_1,nbeta_2

    A = nucleus_attributes%mass_number_int
    beta2_1 = constraint%betac(q1)
    beta3_1 = constraint%bet3c(q1)
    beta2_2 = constraint%betac(q2)
    beta3_2 = constraint%bet3c(q2)
    if(beta2_1 >= 0.d0) signb21 = '+'
    if(beta2_1 < 0.d0)  signb21 = '-'
    if(beta3_1 >= 0.d0) signb31 = '+'
    if(beta3_1 < 0.d0)  signb31 = '-'
    if(beta2_2 >= 0.d0) signb22 = '+'
    if(beta2_2 < 0.d0)  signb22 = '-'
    if(beta3_2 >= 0.d0) signb32 = '+'
    if(beta3_2 < 0.d0)  signb32 = '-'
    !------
    AMPType = Proj_option%AMPType
    if(AMPType==0) AMP = 0 + 48
    if(AMPType==1) AMP = 1 + 48
    if(AMPType==3) AMP = 3 + 48
    !-------
    abs2c1 = abs(beta2_1)
    abs3c1 = abs(beta3_1)
    name1(1) = abs2c1 + 48 !In ASCII, character '0' start from 48. 
    name1(2) = mod(abs2c1*10,10.d0)+48
    name1(3) = mod(abs2c1*100,10.d0)+48
    name1(4) = abs3c1+48
    name1(5) = mod(abs3c1*10,10.d0)+48
    name1(6) = mod(abs3c1*100,10.d0)+48
    !------
    abs2c2 = abs(beta2_2)
    abs3c2 = abs(beta3_2)
    name2(1) = abs2c2 + 48 !In ASCII, character '0' start from 48. 
    name2(2) = mod(abs2c2*10,10.d0)+48
    name2(3) = mod(abs2c2*100,10.d0)+48
    name2(4) = abs3c2+48
    name2(5) = mod(abs3c2*10,10.d0)+48
    name2(6) = mod(abs3c2*100,10.d0)+48
    !-----
    name_nf1 = mod(BS%HO_sph%n0f/10,10) + 48
    name_nf2 = mod(BS%HO_sph%n0f,10) + 48
    !-----
    nphi_1 = mod(projection_mesh%nphi(1)/10,10) + 48
    nphi_2 = mod(projection_mesh%nphi(1),10) + 48
    nbeta_1 = mod(projection_mesh%nbeta/10,10) + 48
    nbeta_2 = mod(projection_mesh%nbeta,10) + 48
    Proj_outputfile%outputelem = OUTPUT_PATH//'Proj_'//int2str(A)//nucleus_attributes%name &
                        //'_kern.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) &
                        //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.elem'
    Proj_outputfile%outputDsME1B = OUTPUT_PATH//'Proj_'//int2str(A)//nucleus_attributes%name &
                        //'_D1B.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) &
                        //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.me'
    Proj_outputfile%outputTDME1B = OUTPUT_PATH//'Proj_'//int2str(A)//nucleus_attributes%name &
                        //'_TD1B.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) &
                        //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.me'
    Proj_outputfile%outputTDME1B_c = OUTPUT_PATH//'Proj_'//int2str(A)//nucleus_attributes%name &
                        //'_TD1B.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6)) &
                        //'_'//signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) //'.me'
    Proj_outputfile%outputEMme = OUTPUT_PATH//'jScheme_EM'//'_A'//int2str(A) &
                        //'_eMax'//char(name_nf1)//char(name_nf2)//'.me'       
    Proj_outputfile%outputm1Bme = OUTPUT_PATH//'mScheme_1B'//'_A'//int2str(A) &
                        //'_eMax'//char(name_nf1)//char(name_nf2)//'.me'
    Proj_outputfile%outputEccentricityKernel = OUTPUT_PATH//'Proj_'//int2str(A)//nucleus_attributes%name &
                        //'_Eccen.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) &
                        //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.elem'
    Proj_outputfile%outputR2bodyKernel = OUTPUT_PATH//'Proj_'//int2str(A)//nucleus_attributes%name &
                        //'_R2_2b.'//char(AMP)//'D'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2)//'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) &
                        //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.elem'

end subroutine

subroutine write_kernels
    use Globals, only: gcm_space,kernels,Proj_option
    integer :: J,K1,K2,parity
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(len=*), parameter ::  format1 = "(3i5,4x,a,4x,3f12.3)", &
                                    format2 = "(4e15.8)"
    open(Proj_outputfile%u_outputelem ,form='formatted',file=Proj_outputfile%outputelem)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            do K1 = -0,0
                do K2 = -0,0
                    ! In the axially symmetric case, the kernel is non-zero only when
                    ! the parity satisfies  Pi  = (-1)^J for  N_KK, H_KK, X_KK and E0_KK
                    ! the parity satisfies Pi_i = (-1)^J_i for  Q2_KK_12
                    if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                        parity = 1 ! +
                    else
                        parity = 2 ! -
                    end if
                    write(Proj_outputfile%u_outputelem,format1) J,K1,K2,ParityChar(parity),&
                                                                Real(kernels%N2_KK(J,K1,K2,1,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-6)), & ! N^2
                                                                Real(kernels%N2_KK(J,K1,K2,2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-6)), & ! Z^2
                                                                Real(kernels%J2_KK(J,K1,K2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-6))      ! J^2
                    write(Proj_outputfile%u_outputelem,format2) kernels%N_KK(J,K1,K2,parity), &
                                                                kernels%H_KK(J,K1,K2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)
                    write(Proj_outputfile%u_outputelem,format2) kernels%X_KK(J,K1,K2,1,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30),&
                                                                kernels%X_KK(J,K1,K2,2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)
                    ! ! proton part
                    write(Proj_outputfile%u_outputelem,format2) kernels%Q2_KK_12(J,K1,K2,2,parity),&
                                                                kernels%Q2_KK_21(J,K1,K2,2,parity)
                    write(Proj_outputfile%u_outputelem,format2) kernels%E0_KK(J,K1,K2,2,parity), &
                                                                kernels%E0_KK(J,K1,K2,2,parity)
                    ! ! neutron part
                    write(Proj_outputfile%u_outputelem,format2) kernels%Q2_KK_12(J,K1,K2,1,parity),&
                                                                kernels%Q2_KK_21(J,K1,K2,1,parity)
                    write(Proj_outputfile%u_outputelem,format2) kernels%E0_KK(J,K1,K2,1,parity), &
                                                                kernels%E0_KK(J,K1,K2,1,parity)
                end do
            end do
        end do 
    close(Proj_outputfile%u_outputelem)
end subroutine

subroutine write_r2_2body_operators_kernels(q1,q2)
    use Constants, only: pi
    use Globals, only: kernels,constraint,Proj_option
    integer :: q1,q2
    integer :: J,K1,K2,parity,ch
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(len=*), parameter ::  format1 = "(3i5,4x,a,4(4x,f5.3))", &
                                    format2 = "(3e15.8)"
    open(Proj_outputfile%u_outputR2bodyKernel,form='formatted',file=Proj_outputfile%outputR2bodyKernel)
        do J = 0,0 
            do K1 = -0,0
                do K2 = -0,0
                    ! In the axially symmetric case, the kernel is non-zero only when
                    ! the parity satisfies  Pi  = (-1)^J for  N_KK, H_KK, X_KK and E0_KK
                    ! the parity satisfies Pi_i = (-1)^J_i for  Q2_KK_12
                    if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                        parity = 1 ! +
                    else
                        parity = 2 ! -
                    end if
                    
                    ! 
                    ! write(Proj_outputfile%u_outputR2bodyKernel,'(A)') "# beta2_1 beta3_1 beta2_2 beta3_2  J  K1 K2 Parity   1B_KK      2B_KK      1B+2B"
                    
                    ! write(Proj_outputfile%u_outputR2bodyKernel,format1)  J,K1,K2,ParityChar(parity),constraint%betac(q1),constraint%bet3c(q1),constraint%betac(q2),constraint%bet3c(q2)

                    write(Proj_outputfile%u_outputR2bodyKernel,format2)   Real((4.0d0*pi/3.0d0) * SUM(kernels%r2_2b_KK(J,K1,K2,parity,1,:))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! 1B
                                                                                Real((4.0d0*pi/3.0d0) * SUM(kernels%r2_2b_KK(J,K1,K2,parity,2,:))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! 2B
                                     Real((4.0d0*pi/3.0d0)*(SUM(kernels%r2_2b_KK(J,K1,K2,parity,1,:)+kernels%r2_2b_KK(J,K1,K2,parity,2,:))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)))   ! 1B +2B
                    ! Each Term
                    write(Proj_outputfile%u_outputR2bodyKernel,format2)   Real((4.0d0*pi/3.0d0) * SUM(kernels%r2_2b_KK(J,K1,K2,parity,3,:))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! direct
                                                                                Real((4.0d0*pi/3.0d0)* SUM(kernels%r2_2b_KK(J,K1,K2,parity,4,:))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! exchange
                                                                                Real((4.0d0*pi/3.0d0)* SUM(kernels%r2_2b_KK(J,K1,K2,parity,5,:))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))    ! kappa
                    do ch=1,2,3 ! n,p,pn
                        write(Proj_outputfile%u_outputR2bodyKernel,format2)  Real((4.0d0*pi/3.0d0)* r2_2b_KK(J,K1,K2,parity,1,ch) /(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! 1B
                        Real((4.0d0*pi/3.0d0)* kernels%r2_2b_KK(J,K1,K2,parity,2,ch) /(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! 2B
                        Real((4.0d0*pi/3.0d0)*(kernels%r2_2b_KK(J,K1,K2,par ity,1,ch)+ kernels%r2_2b_KK(J,K1,K2,parity,2,ch)) /(kernels%N_KK(J,K1,K2,parity)+1.0d-30))  ! 1B+2B

                        write(Proj_outputfile%u_outputR2bodyKernel,format2)  Real((4.0d0*pi/3.0d0)* kernels%r2_2b_KK(J,K1,K2,parity,3,ch) /(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! direct
                        Real((4.0d0*pi/3.0d0)* kernels%r2_2b_KK(J,K1,K2,parity,4,ch) /(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! exchange
                        Real((4.0d0*pi/3.0d0)* kernels%r2_2b_KK(J,K1,K2,parity,5,ch) /(kernels%N_KK(J,K1,K2,parity)+1.0d-30))    ! kappa

                    end do


                    ! ! Separated n/p/np
                    ! write(Proj_outputfile%u_outputR2bodyKernel,format2)   Real((4.0d0*pi/3.0d0)* kernels%Ecc_sep_KK(J,K1,K2,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! n_1B
                    !                                                             Real((4.0d0*pi/3.0d0)* kernels%Ecc_sep_KK(J,K1,K2,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! p_1B
                    !                                                             Real((4.0d0*pi/3.0d0)* kernels%Ecc_sep_KK(J,K1,K2,parity,3)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))    ! np_2B
                
                    ! write(Proj_outputfile%u_outputR2bodyKernel,*) "By Density:"
                    ! ! ! 1B part
                    ! write(Proj_outputfile%u_outputR2bodyKernel,format2)   Real(kernels%r2_2b_KK_byDensity(J,K1,K2,1,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                    !                                                             Real(kernels%r2_2b_KK_byDensity(J,K1,K2,2,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                    !      Real((kernels%r2_2b_KK_byDensity(J,K1,K2,1,parity,1)+kernels%r2_2b_KK_byDensity(J,K1,K2,2,parity,1))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))
                    ! ! ! 2B part
                    ! write(Proj_outputfile%u_outputR2bodyKernel,format2)   Real(kernels%r2_2b_KK_byDensity(J,K1,K2,1,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                    !                                                             Real(kernels%r2_2b_KK_byDensity(J,K1,K2,2,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                    !      Real((kernels%r2_2b_KK_byDensity(J,K1,K2,1,parity,2)+kernels%r2_2b_KK_byDensity(J,K1,K2,2,parity,2))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))
                    ! ! ! np and pn part
                    ! write(Proj_outputfile%u_outputR2bodyKernel,format2)   Real(kernels%r2_2b_KK_byDensity(J,K1,K2,3,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)),& ! np part
                    !                                                             Real(kernels%r2_2b_KK_byDensity(J,K1,K2,3,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! pn part
                    !     Real((kernels%r2_2b_KK_byDensity(J,K1,K2,3,parity,1)+kernels%r2_2b_KK_byDensity(J,K1,K2,3,parity,2))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))
                    ! ! ! total 
                    ! write(Proj_outputfile%u_outputR2bodyKernel,format2)   Real((kernels%r2_2b_KK_byDensity(J,K1,K2,2,parity,1)+kernels%r2_2b_KK_byDensity(J,K1,K2,2,parity,2)&
                    !                                                                 + kernels%r2_2b_KK_byDensity(J,K1,K2,1,parity,1)+kernels%r2_2b_KK_byDensity(J,K1,K2,1,parity,2)&
                    !                                                                 + kernels%r2_2b_KK_byDensity(J,K1,K2,3,parity,1)+kernels%r2_2b_KK_byDensity(J,K1,K2,3,parity,2))&
                    !                                                                 /(kernels%N_KK(J,K1,K2,parity)+1.0d-30))  ! 1B +2B + np part                                       
                end do
            end do
        end do 
    close(Proj_outputfile%u_outputR2bodyKernel)
end subroutine

subroutine write_eccentricity_operators_kernels(q1,q2)
    use Globals, only: kernels,constraint,Proj_option
    integer :: q1,q2
    integer :: J,K1,K2,parity
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(len=*), parameter ::  format1 = "(3i5,4x,a,4(4x,f5.3))", &
                                    format2 = "(3e15.8)"
    open(Proj_outputfile%u_outputEccentricityKernel,form='formatted',file=Proj_outputfile%outputEccentricityKernel)
        do J = 0,0 
            do K1 = -0,0
                do K2 = -0,0
                    ! In the axially symmetric case, the kernel is non-zero only when
                    ! the parity satisfies  Pi  = (-1)^J for  N_KK, H_KK, X_KK and E0_KK
                    ! the parity satisfies Pi_i = (-1)^J_i for  Q2_KK_12
                    if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                        parity = 1 ! +
                    else
                        parity = 2 ! -
                    end if
                    write(Proj_outputfile%u_outputEccentricityKernel,format1)  J,K1,K2,ParityChar(parity),constraint%betac(q1),constraint%bet3c(q1),constraint%betac(q2),constraint%bet3c(q2)

                    write(Proj_outputfile%u_outputEccentricityKernel,format2)   Real(kernels%Eccentricity_KK(J,K1,K2,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! 1B
                                                                                Real(kernels%Eccentricity_KK(J,K1,K2,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! 2B
                                     Real((kernels%Eccentricity_KK(J,K1,K2,parity,1)+kernels%Eccentricity_KK(J,K1,K2,parity,2))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))   ! 1B +2B
                    ! Each Term
                    write(Proj_outputfile%u_outputEccentricityKernel,format2)   Real(kernels%Eccentricity_KK(J,K1,K2,parity,3)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! direct
                                                                                Real(kernels%Eccentricity_KK(J,K1,K2,parity,4)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! exchange
                                                                                Real(kernels%Eccentricity_KK(J,K1,K2,parity,5)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))    ! kappa
                    ! Separated n/p/np
                    write(Proj_outputfile%u_outputEccentricityKernel,format2)   Real(kernels%Ecc_sep_KK(J,K1,K2,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! n_1B
                                                                                Real(kernels%Ecc_sep_KK(J,K1,K2,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! p_1B
                                                                                Real(kernels%Ecc_sep_KK(J,K1,K2,parity,3)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))    ! np_2B
                
                    write(Proj_outputfile%u_outputEccentricityKernel,*) "By Density:"
                    ! ! 1B part
                    write(Proj_outputfile%u_outputEccentricityKernel,format2)   Real(kernels%Eccentricity_KK_byDensity(J,K1,K2,1,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                                                                                Real(kernels%Eccentricity_KK_byDensity(J,K1,K2,2,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                         Real((kernels%Eccentricity_KK_byDensity(J,K1,K2,1,parity,1)+kernels%Eccentricity_KK_byDensity(J,K1,K2,2,parity,1))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))
                    ! ! 2B part
                    write(Proj_outputfile%u_outputEccentricityKernel,format2)   Real(kernels%Eccentricity_KK_byDensity(J,K1,K2,1,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                                                                                Real(kernels%Eccentricity_KK_byDensity(J,K1,K2,2,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                         Real((kernels%Eccentricity_KK_byDensity(J,K1,K2,1,parity,2)+kernels%Eccentricity_KK_byDensity(J,K1,K2,2,parity,2))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))
                    ! ! np and pn part
                    write(Proj_outputfile%u_outputEccentricityKernel,format2)   Real(kernels%Eccentricity_KK_byDensity(J,K1,K2,3,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)),& ! np part
                                                                                Real(kernels%Eccentricity_KK_byDensity(J,K1,K2,3,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! pn part
                        Real((kernels%Eccentricity_KK_byDensity(J,K1,K2,3,parity,1)+kernels%Eccentricity_KK_byDensity(J,K1,K2,3,parity,2))/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))
                    ! ! total 
                    write(Proj_outputfile%u_outputEccentricityKernel,format2)   Real((kernels%Eccentricity_KK_byDensity(J,K1,K2,2,parity,1)+kernels%Eccentricity_KK_byDensity(J,K1,K2,2,parity,2)&
                                                                                    + kernels%Eccentricity_KK_byDensity(J,K1,K2,1,parity,1)+kernels%Eccentricity_KK_byDensity(J,K1,K2,1,parity,2)&
                                                                                    + kernels%Eccentricity_KK_byDensity(J,K1,K2,3,parity,1)+kernels%Eccentricity_KK_byDensity(J,K1,K2,3,parity,2))&
                                                                                    /(kernels%N_KK(J,K1,K2,parity)+1.0d-30))  ! 1B +2B + np part                                       
                end do
            end do
        end do 
    close(Proj_outputfile%u_outputEccentricityKernel)
end subroutine

! write expectation of different (q1,q2)
subroutine write_Proj_expectation(q1,q2)
    use Globals,only: gcm_space,Proj_outputfile,kernels,constraint,nucleus_attributes,Proj_option
    integer :: q1, q2,J,K1,K2,parity
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(len=*), parameter ::  format1 = "(4a9,3(a2,2x),a6,(1x,a9,2x),(a10,2x),(6x,a,5x),2(6x,a,3x),2(a9,2x),2(a7,2x))", &
                                    format2 = "(4(2x,f6.3,1x),3(i2,2x),(3x,a,3x),(f9.6,2x),2(f10.5,2x),2(f8.3,2x),2(f9.3,2x),(f7.3,2x),(f7.4,2x))"
    if(first_kernel) then
        write(Proj_outputfile%u_outExpectation,format1) "beta2_1","beta3_1","beta2_2","beta3_2","J","K1","K2","Parity",&
                                                        "N_Kernel","H_kernel","E","N","Z","N^2","Z^2","J^2", &
                                                        "r^2_p","r^2_2Bp"
        first_kernel = .False.
    endif
    do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
        do K1 = -0,0
            do K2 = -0,0
                if ((-1)**J == 1 .or. Proj_option%PPtype==0) then
                    parity = 1 ! +
                else
                    parity = 2 ! -
                end if
                write(Proj_outputfile%u_outExpectation,format2) constraint%betac(q1),constraint%bet3c(q1), &
                                    constraint%betac(q2),constraint%bet3c(q2),J,K1,K2,ParityChar(parity),  &
                                    Real(kernels%N_KK(J,K1,K2,parity)),Real(kernels%H_KK(J,K1,K2,parity)), &
                                    Real(kernels%H_KK(J,K1,K2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                                    Real(kernels%X_KK(J,K1,K2,1,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), &
                                    Real(kernels%X_KK(J,K1,K2,2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)),  &
                                    Real(kernels%N2_KK(J,K1,K2,1,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! N^2
                                    Real(kernels%N2_KK(J,K1,K2,2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)), & ! Z^2
                                    Real(kernels%J2_KK(J,K1,K2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30)),   & ! J^2
                                    (Real(kernels%E0_KK(J,K1,K2,2,parity)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))/nucleus_attributes%proton_number)**0.5d0, &!R_p  &
                                    (Real((4.0d0*pi/3.0d0) * kernels%r2_2b_KK(J,K1,K2,parity,1)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))/nucleus_attributes%proton_number)**0.5d0, & ! 1B
                                    (Real((4.0d0*pi/3.0d0) * kernels%r2_2b_KK(J,K1,K2,parity,2)/(kernels%N_KK(J,K1,K2,parity)+1.0d-30))/nucleus_attributes%proton_number)**0.5d0    ! 2B
            end do
        end do
    end do
    write(Proj_outputfile%u_outExpectation, '(/)')
end subroutine

! density matrix elements
subroutine write_1B_density_matrix_elements
    use Globals, only: gcm_space,Proj_option,BS,Proj_outputfile,Proj_densities
    integer :: J,K1_start,K1_end,K2_start,K2_end,K1,K2,iParity,Parity,ifg,m1,m2
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(1) :: Parity_c
    open(Proj_outputfile%u_outputDsME1B,form='formatted',file=Proj_outputfile%outputDsME1B)
    write(Proj_outputfile%u_outputDsME1B,*) "Pi   J  K1  K2  ifg  m1   m2    neutron            proton"
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
                    Parity_c = ParityChar(iParity)
                    if( Parity /= (-1)**J) cycle
                    do ifg = 1, 2
                    do m1 = 1, BS%HO_sph%idsp(1,ifg)
                        do m2 = 1, BS%HO_sph%idsp(1,ifg)
                        ! 1 Body
                        write(Proj_outputfile%u_outputDsME1B,"(1x,a1,4i4,2i5,2x,2f18.14)")Parity_c,J,K1,K2,ifg,m1,m2,&
                                    Real(Proj_densities%ME1B(J,K1,K2,iParity,1,ifg,m1,m2)),&
                                    Real(Proj_densities%ME1B(J,K1,K2,iParity,2,ifg,m1,m2))
                        end do
                    end do
                    end do 
                end do 
            end do 
        end do 
    end do
    close(Proj_outputfile%u_outputDsME1B)
end subroutine

subroutine write_reduced_1B_transition_density_matrix_elements(q1,q2)
    use Globals, only: Proj_outputfile,gcm_space,Proj_option,TDs
    integer, intent(in) :: q1,q2
    integer :: J,Ji,Jf,lambda,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,ifg,a,b,Parity_i,Parity_f,iPi,iPf
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(1) :: Parity_f_c,Parity_i_c
    ! q1-q2
    open(Proj_outputfile%u_outputTDME1B ,form='formatted',file=Proj_outputfile%outputTDME1B)
    write(Proj_outputfile%u_outputTDME1B,*) "Pf Pi  Jf  Ji  l   Kf  Ki  ifg a   b    neutron            proton"
    ! q2-q1
    if (Proj_option%Kernel_Symmetry == 1 .and. q1/=q2) then
        open(Proj_outputfile%u_outputTDME1B_c ,form='formatted',file=Proj_outputfile%outputTDME1B_c)
        write(Proj_outputfile%u_outputTDME1B_c,*) "Pf Pi  Jf  Ji  l   Kf  Ki  ifg a   b    neutron            proton"
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
                            do a = 1, TDs%nlj_length(ifg)
                                do b = 1, TDs%nlj_length(ifg)
                                    Parity_i = (-1)**Ji ! In the axially symmetric case, the kernel is non-zero only when Parity_i * (-1)^J_i = 1
                                    iPi = (3-Parity_i)/2 ! +1: 1, -1: 2
                                    Parity_i_c = ParityChar(iPi)
                                    do iPf = 1,2
                                        Parity_f = (-1)**(iPf+1) ! 1: +1 , 2: -1
                                        Parity_f_c = ParityChar(iPf)
                                        write(Proj_outputfile%u_outputTDME1B,"(1x,a1,2x,a1,8i4,2x,2f18.14)") Parity_f_c,Parity_i_c,Jf,Ji,lambda,Kf,Ki,ifg,a,b, &
                                            Real(TDs%reduced_TDME1B(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,1)),&
                                            Real(TDs%reduced_TDME1B(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,2))
                                        if (Proj_option%Kernel_Symmetry == 1 .and. q1/=q2) then
                                            write(Proj_outputfile%u_outputTDME1B_c,"(1x,a1,2x,a1,8i4,2x,2f18.14)")Parity_f_c,Parity_i_c,Jf,Ji,lambda,Kf,Ki,ifg,a,b, &
                                                Real(TDs%reduced_TDME1B_c(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,1)),&
                                                Real(TDs%reduced_TDME1B_c(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,2))
                                        end if 
                                    end do
                                end do 
                            end do 
                        end do
                    end do 
                end do
            end do
        end do 
    end do
    close(Proj_outputfile%u_outputTDME1B)
    if (Proj_option%Kernel_Symmetry == 1 .and. q1/=q2) then 
        close(Proj_outputfile%u_outputTDME1B_c)
    end if 
end subroutine

! matrix elements of operator
subroutine write_1B_operators_matrix_elements
    use Globals, only: BS, Proj_outputfile
    use Eccentricity, only: f_n,eccentricity_matrix_element_one_body
    use EM, only: rl_nl,ylm_ljm
    integer :: ifg,ndsp,i0sp,m1,m2,nr1,nl1,nj1,nm1,nr2,nl2,nj2,nm2
    real(r64) :: r2, r4, r2Y20, r4Y20, r4Y40, fn, e_1B
    open(Proj_outputfile%u_outputm1Bme ,form='formatted',file=Proj_outputfile%outputm1Bme) 
    write(Proj_outputfile%u_outputm1Bme,'(9A5,2A7,7A16)')  'ifg','m1','m2','n1','n2','l1','l2','2j1','2j2','2j_m1','2j_m2', &
                                                         'r^2','r^4',"r^2Y20","r^4Y20","r^4Y40",'f2','Eps1B'
    do ifg = 1, 2
        ndsp = BS%HO_sph%idsp(1,ifg)
        i0sp = BS%HO_sph%iasp(1,ifg)
        do m1 = 1, ndsp
            do m2= 1, ndsp
                nr1 = BS%HO_sph%nljm(i0sp+m1,1) ! n1  
                nl1 = BS%HO_sph%nljm(i0sp+m1,2) ! l1
                nj1 = BS%HO_sph%nljm(i0sp+m1,3) ! j1 + 1/2 
                nm1 = BS%HO_sph%nljm(i0sp+m1,4) ! m1 + 1/2
        
                nr2 = BS%HO_sph%nljm(i0sp+m2,1) ! n2
                nl2 = BS%HO_sph%nljm(i0sp+m2,2) ! l2
                nj2 = BS%HO_sph%nljm(i0sp+m2,3) ! j2 + 1/2
                nm2 = BS%HO_sph%nljm(i0sp+m2,4) ! m2 + 1/2

                ! r^n  ! <n1 l1 j1 m1 | r^n | n2 l2 j2 m2> = <n1 l1 | r^n | n2 l2> * <l1 j1 m1 | l2 j2 m2>
                if(nl1==nl2 .and. nj1==nj2 .and. nm1==nm2) then 
                    r2 =  rl_nl(nr1,nl1,2,nr2,nl2) !  <m1|r^2|m2>
                    r4 =  rl_nl(nr1,nl1,4,nr2,nl2) !  <m1|r^4|m2>
                else
                    r2 = 0.d0
                    r4 = 0.d0
                end if 
                ! r^2 Y20
                r2Y20 = rl_nl(nr1,nl1,2,nr2,nl2)*ylm_ljm(nl1,nj1,nm1,2,0,nl2,nj2,nm2) 
                ! r^4 Y20
                r4Y20 = rl_nl(nr1,nl1,4,nr2,nl2)*ylm_ljm(nl1,nj1,nm1,2,0,nl2,nj2,nm2)
                ! r^4 Y40
                r4Y40 = rl_nl(nr1,nl1,4,nr2,nl2)*ylm_ljm(nl1,nj1,nm1,4,0,nl2,nj2,nm2)
                ! eccentricity operator matrix element
                fn = f_n(ifg,m1,ifg,m2,2)
                call eccentricity_matrix_element_one_body(ifg,m1,ifg,m2,2,e_1B)
                write(Proj_outputfile%u_outputm1Bme,"(9i5,2i7,1x,7(f15.5,1x))") ifg,m1,m2,nr1,nr2,nl1,nl2,2*nj1-1,2*nj2-1,2*nm1-1,2*nm2-1,&
                                                                            r2,r4,r2Y20,r4Y20,r4Y40,fn,e_1B
            end do 
        end do 
    end do
end subroutine

subroutine  write_reduced_1B_multipole_matrix_elements
    use Globals, only: Proj_outputfile,BS,TDs,gcm_space
    use EM, only: reduced_multipole_matrix_elements,reduced_monopole_matrix_elements,rl_nl
    use TD, only: set_nlj_mapping
    integer :: ifg,lambda_start,lambda_end,lambda,a,b,i0sp,a_index,nra,nla,nja,b_index,nrb,nlb,njb
    real(r64),allocatable,dimension(:) :: Ql_ab
    real(r64) :: monopole_ab
    character(len=200) :: header, temp
    call set_nlj_mapping
    open(Proj_outputfile%u_outputEMme ,form='formatted',file=Proj_outputfile%outputEMme)
    ! store the nlj of a/b
    write(Proj_outputfile%u_outputEMme,*) "--------------------------------------------------------"
    write(Proj_outputfile%u_outputEMme,*) " n l j of a/b "
    write(Proj_outputfile%u_outputEMme,*) "--------------------------------------------------------"
    write(Proj_outputfile%u_outputEMme,*) " ifg a/b  n   l   j"   
    do ifg =1, 2 
        do a = 1, TDs%nlj_length(ifg)
            i0sp = BS%HO_sph%iasp(1,ifg)
            a_index = TDs%nlj_index(a,ifg)
            nra= BS%HO_sph%nljm(i0sp+a_index,1) ! n_r
            nla= BS%HO_sph%nljm(i0sp+a_index,2) ! l
            nja= BS%HO_sph%nljm(i0sp+a_index,3) ! j +1/2
            write(Proj_outputfile%u_outputEMme,"(5i4,'/2')") ifg, a, nra, nla, (nja*2-1)
        end do 
    end do

    ! store <a||Q_lambda||b>
    lambda_start = 1
    lambda_end = min(gcm_space%Jmax,TDs%lambda_max)
    header = '   a   b  sqrt(4pi)<a||r^2Y_0||b>'
    do lambda = lambda_start, lambda_end
        write(temp,'(A,I0,A)') '        <a||Q', lambda, '||b>'
        header = trim(header)//temp
    end do

    write(Proj_outputfile%u_outputEMme,*) "------------------------------------------------------------------------"
    write(Proj_outputfile%u_outputEMme,*) " reduced single-particle matrix element of the electromagnetic operator "
    write(Proj_outputfile%u_outputEMme,*) "-------------------------------------------------------------------------"
    write(Proj_outputfile%u_outputEMme,'(A)') trim(header)

    allocate(Ql_ab(lambda_start:lambda_end))
    ifg =  2 
    do a = 1, TDs%nlj_length(ifg)
        do b = 1, TDs%nlj_length(ifg)
            i0sp = BS%HO_sph%iasp(1,ifg)
            a_index = TDs%nlj_index(a,ifg)
            nra= BS%HO_sph%nljm(i0sp+a_index,1) ! n_r
            nla= BS%HO_sph%nljm(i0sp+a_index,2) ! l
            nja= BS%HO_sph%nljm(i0sp+a_index,3) ! j +1/2
            b_index = TDs%nlj_index(b,ifg)
            nrb= BS%HO_sph%nljm(i0sp+b_index,1) ! n_r
            nlb= BS%HO_sph%nljm(i0sp+b_index,2) ! l
            njb= BS%HO_sph%nljm(i0sp+b_index,3) ! j +1/2
            do lambda = lambda_start, lambda_end
                call reduced_multipole_matrix_elements(nra,nla,nja,lambda,nrb,nlb,njb,Ql_ab(lambda))
            end do
            call reduced_monopole_matrix_elements(nra,nla,nja,nrb,nlb,njb,monopole_ab)
            write(Proj_outputfile%u_outputEMme,'(2i4,1x,f17.10,8x,15(1x,f17.10))') a,b, monopole_ab,(Ql_ab(lambda), lambda=lambda_start,lambda_end)
        end do 
    end do
end subroutine

END MODULE Proj_Inout