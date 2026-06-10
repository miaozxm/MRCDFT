!==============================================================================!
! MODULE GCM                                                                  !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE GCM
implicit none
contains
    subroutine GCM_main
        use Globals, only: GCM_outputfile
        use GCM_Inout, only: set_GCM_output_filename, read_kernels
        use HWG, only: generator_coordinate_basis,solve_HWG_equation
        use GCM_Observables, only: calculate_observables
        call set_GCM_output_filename
        open(GCM_outputfile%u_outGCM_standard ,form='formatted',file=GCM_outputfile%outGCM_standard)
        call read_kernels
        call generator_coordinate_basis
        call solve_HWG_equation
        call calculate_observables
        close(GCM_outputfile%u_outGCM_standard)
    end subroutine
END MODULE