!==============================================================================!
! PROGRAM xxx                                                                  !
!                                                                              !
! This program ...                                                             !
!                                                                              !
! Licence:                                                                     !
! Article describing the code:                                                 !
! GitHub repository:                                                           !
!==============================================================================!
PROGRAM MR_CDFT
    use omp_lib
    use CDFT
    use Proj
    use Globals, only: Proj_option,outputfile
    use MathMethods, only: math_gfv
    use CDFT_Inout, only: handle_input_config,read_file_b23,read_CDFT_configuration
    use Proj_Inout, only: read_Proj_configuration
    use Forces, only : set_force_parameters
    use Nucleus, only: set_nucleus_attributes
    implicit none
    integer :: num_threads
    real :: start_time,end_time,start_CPU_time,end_CPU_time,CPU_time_seconds,CPU_time_minutes,CPU_time_hours

    call handle_input_config
    ! Get the starting CPU time
    call CPU_TIME(start_CPU_time)
    start_time =  omp_get_wtime()
    ! Print number of threads
    !$OMP PARALLEL
    !$OMP MASTER
    num_threads = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL
    write(*,*) "OpenMP will create ", num_threads, " threads."

    write(*,*)"... PROGRAM Start!..."
    call read_file_b23
    call read_CDFT_configuration(.True.)
    call read_Proj_configuration(.True.)
    
    open(outputfile%u_config, file=outputfile%config, status='unknown')
    call math_gfv
    call set_nucleus_attributes(.True.)
    call set_force_parameters(.True.)
    if(Proj_option%ProjectionType == 0 .or. Proj_option%ProjectionType==1) then
        call CDFT_Main
    endif
    if(Proj_option%ProjectionType == 1 .or. Proj_option%ProjectionType==2) then
        call Proj_Main
    endif
    close(outputfile%u_config)
    write(*,*) " ...PROGRAM END!..."

    ! Get the ending CPU time
    call CPU_TIME(end_CPU_time)
    end_time =  omp_get_wtime()

    ! Calculate elapsed time
    CPU_time_seconds = end_CPU_time - start_CPU_time
    CPU_time_minutes = CPU_time_seconds / 60.0
    CPU_time_hours = CPU_time_seconds / 3600.0
    ! Print CPU time
    write(*,*) "CPU Time used (minutes): ", CPU_time_minutes
    write(*,*) "CPU Time used (hours): ", CPU_time_hours
    write(*,*) 'Real time (minutes)',(end_time - start_time)/60.0
END PROGRAM MR_CDFT