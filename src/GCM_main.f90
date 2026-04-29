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
        use GCM_Inout, only: read_kernels
        use HWG, only: generator_coordinate_basis,solve_HWG_equation
        use GCM_Observables, only: calculate_observables
        call read_kernels
        call generator_coordinate_basis
        call solve_HWG_equation
        call calculate_observables
    end subroutine
END MODULE