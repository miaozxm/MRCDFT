!==============================================================================!
! MODULE CDFT                                                                  !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE CDFT
contains
    subroutine CDFT_Main
        use Globals, only: constraint, iteration, outputfile,expectations,option,pairing,MPI_Infor
        use CDFT_Inout, only: set_CDFT_Expectation_filename,set_CDFT_output_filename,write_result_DIR,adjust_left,int2str
        use Field, only: set_woodssaxon_parameters,calculate_meson_propagators,initial_potential_fields,calculate_fields
        use Forces, only : calculate_density_dependence_of_coupling_constants
        use Basis, only : set_Cylindrical_HO_basis,set_Spherical_HO_basis,transform_coefficients_form_cylindrical_to_spherical
        use Constraints, only: set_constraint_parameters,calculate_constraint_potential_coefficients
        use BCS, only: initial_pairing_field,calculate_occupation,set_block_level_of_KPi
        use DiracEquation, only: calculate_sigma_nabla,solve_dirac_equation
        use Density, only: calculate_densities_DIR,calculate_densities_RHB
        use Expectation, only: calculate_expectation_DIR,calculate_expectation_RHB
        use Broyden, only: set_broyden_parameters,mix_potentials_DIR,mix_matrix_elements_RHB
        use DeltaField, only: calculate_pairing_matrix_element_W,set_separable_pairing_parameters,initial_delta_field
        use RHBEquation, only: solve_RHB_equation

        implicit none
        logical :: ifPrint=.False.
        integer :: constraint_index,iteration_index
        character(len=150) :: format1,format2,format3,format4,format5

        if (MPI_Infor%rank == 0) write(*,*)'CDFT_Main: Start CDFT calculations' 
        format1 = "(112(1h=),/, 48x,'Iteration: ',i4,48x,/,(112(1h-)))"
        format2 = "(i3,'. It. si =',f15.10,'  E/A =',f15.10,'  R =',f15.10,'  b2 =',f15.10,'  b3 =',f15.10)"
        ! format2 = "(i3,'. It. si =',f30.25,'  E =',f30.25,'  R =',f30.25,'  b2 =',f30.25,'  b3 =',f30.25)"
        format3 = "(/,112(1h#),/,'#',28x ,'Iteration converged after',i4,' steps, si =',f13.10,28x,'#',/,112(1h#))"
        format4 = "('     Constraint: ',i4,'/',a, 'beta2=',f6.2,'   beat3=',f6.2)"
        format5 = "('     Iteration converged after',i4,' steps, si =',f13.8, ',  beta2 =',f6.3,',  beta3 =',f6.3,',',i4,'/',a)"

        call set_woodssaxon_parameters(.True.)
        call set_Cylindrical_HO_basis(.True.)
        call set_Spherical_HO_basis(.True.)
        call set_broyden_parameters

        call calculate_sigma_nabla(ifPrint .and. .False.)
        call calculate_meson_propagators(ifPrint .and. .False.)
        if (option%eqType .ne. 0) then ! for RHB case
            call set_separable_pairing_parameters
            call calculate_pairing_matrix_element_W 
        endif

        call set_CDFT_Expectation_filename
        open(outputfile%u_rotationalE, file=trim(outputfile%rotationalE)//'_'//int2str(MPI_Infor%rank), status='unknown')
        open(outputfile%u_outExpectation, file=trim(outputfile%outExpectation)//'_'//int2str(MPI_Infor%rank), status='unknown')

        do constraint_index = 1, constraint%length ! loop for different deformation parameters
            if (mod(constraint_index-1, MPI_Infor%nprocs) /= MPI_Infor%rank) cycle ! MPI process distribution
            constraint%index = constraint_index  
            call set_constraint_parameters
            call set_CDFT_output_filename(constraint%betac(constraint%index),constraint%bet3c(constraint%index))
            open(outputfile%u_outputf, file=outputfile%outputf, status='unknown')

            call initial_pairing_field(ifPrint .and. .True.) ! for BCS case (set initial fermi energy for RHB )
            call initial_delta_field ! for RHB case
            call initial_potential_fields(ifPrint .and. .True.)

            write(*,format4) constraint%index, adjust_left(constraint%length,4),constraint%betac(constraint%index),constraint%bet3c(constraint%index)
            do iteration_index = 1, iteration%iteration_max ! iteration loop
                iteration%ii = iteration_index
                write(outputfile%u_outputf,format1) iteration_index
                if(option%eqType .eq. 0) then
                    !  Dirac (BCS) equation case
                    ! solve Dirac equation and calculate densities
                    call solve_dirac_equation(ifPrint .and. .True.)
                    call calculate_occupation(ifPrint .and. .True.) ! BCS occupation probabilities
                    call calculate_densities_DIR(ifPrint .and. .True.)
                    ! calculate density dependence of coupling constants and fields by densities
                    call calculate_density_dependence_of_coupling_constants(ifPrint .and. .True.)
                    call calculate_fields
                    ! calculate expectation values
                    call calculate_expectation_DIR(ifPrint .and. .True.)
                    ! calculate new constraint potential coefficients
                    call calculate_constraint_potential_coefficients(ifPrint .and. .True.)
                    ! potentials for next iteration
                    call mix_potentials_DIR(ifPrint .and. .True.)
                else
                    ! RHB equation case
                    ! solve RHB equation and calculate densities
                    call solve_RHB_equation(ifPrint .and. .True.)
                    call calculate_densities_RHB(ifPrint .and. .True.)
                    ! calculate density dependence of coupling constants and fields by densities
                    call calculate_density_dependence_of_coupling_constants(ifPrint .and. .True.)
                    call calculate_fields
                    ! calculate expectation values
                    call calculate_expectation_RHB(ifPrint .and. .True.)
                    ! calculate new constraint potential coefficients
                    call calculate_constraint_potential_coefficients(ifPrint .and. .True.)
                    ! matrix elements for next iteration
                    call mix_matrix_elements_RHB
                endif
                write(outputfile%u_outputf,format2)iteration%ii,iteration%si,expectations%ea,expectations%rms,expectations%betg,expectations%beto

                ! check convergence and exit iteration loop
                if(option%eqType .eq. 0) then
                    ! For Dirac equation
                    ! Different blocking methods for odd-N and odd-Z nuclei
                    if(option%block_method == 2 .and. pairing%allow_block) exit ! exit iteration
                    if(iteration%ii.ge.2 .and. abs(iteration%si).lt.iteration%epsi) then
                        if(option%block_type==0 .or. pairing%allow_block .or. option%block_method==1) exit ! exit iteration
                        if(option%block_type==2) call set_block_level_of_KPi(.True.)
                        pairing%allow_block = .True.
                    endif
                else
                    ! For RHB equation
                    ! The RHB case currently does not support blocking for odd-N or odd-Z nuclei.
                    if(iteration%ii.ge.2 .and. abs(iteration%si).lt.iteration%epsi) exit
                end if
            end do
            write(outputfile%u_outputf,format3) iteration%ii,iteration%si
            write(*,format5) iteration%ii,iteration%si,expectations%betg,expectations%beto,constraint%index,adjust_left(constraint%length,4)

            ! after convergence 
            if(option%eqType .eq. 0) then
                ! transform wave function coefficients(F,G) from cylindrical to spherical basis
                call transform_coefficients_form_cylindrical_to_spherical(ifPrint .and. .True.)
                call write_result_DIR
            endif
            
            close(outputfile%u_outputf)
        end do

        close(outputfile%u_rotationalE)
        close(outputfile%u_outExpectation)
        if (MPI_Infor%rank == 0) write(*,*)'CDFT_Main: CDFT calculations completed'
    end subroutine
END MODULE CDFT

!==============================================================================!
! End of file                                                                  !
!==============================================================================!
