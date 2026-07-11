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
    use mpi
    use omp_lib
    use CDFT
    use Proj
    use GCM
    use Tools, only: adjust_left,int2str
    use Globals, only: option,Proj_option,GCM_option,outputfile,MPI_Infor
    use MathMethods, only: math_gfv
    use CDFT_Inout, only: handle_input_config,read_file_b23,read_CDFT_configuration
    use Proj_Inout, only: read_Proj_configuration
    use GCM_Inout, only: read_GCM_configuration
    use Forces, only : set_force_parameters
    use Nucleus, only: set_nucleus_attributes
    implicit none
    integer :: num_threads
    real :: start_time,end_time,start_CPU_time,end_CPU_time,CPU_time_seconds,CPU_time_minutes,CPU_time_hours
    ! MPI variable
    integer ::  ierr
    
    ! MPI initial
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, MPI_Infor%rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_Infor%nprocs, ierr)
    
    if (MPI_Infor%rank == 0) write(*,*)"... PROGRAM Start!..."
    if (MPI_Infor%rank == 0) then 
        num_threads = omp_get_max_threads()
        write(*,*) "MPI will create ", MPI_Infor%nprocs, " processes."
        write(*,*) "OpenMP will create ", num_threads, " threads."
        ! Get the starting CPU time
        call CPU_TIME(start_CPU_time)
        start_time =  omp_get_wtime()
    end if 
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call handle_input_config
    call read_file_b23
    call read_CDFT_configuration(.True.)
    call read_Proj_configuration(.True.)
    call read_GCM_configuration(.True.)
    open(outputfile%u_config, file=trim(outputfile%config)//'_'//int2str(MPI_Infor%rank), status='unknown')
    call math_gfv
    call set_nucleus_attributes(.True.)
    call set_force_parameters(.True.)

    ! CDFT
    if(option%CDFTType > 0) then
        if (MPI_Infor%rank == 0) write(*,*)'CDFT_Main: Start CDFT calculations'
        call CDFT_Main
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (MPI_Infor%rank == 0) write(*,*)'CDFT_Main: CDFT calculations completed'
    endif

    ! Projection
    if(Proj_option%ProjectionType > 0) then
        if (MPI_Infor%rank == 0) write(*,*)'Proj_Main: Start Proj calculations'
        call Proj_Main
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (MPI_Infor%rank == 0) write(*,*)'Proj_Main: Proj calculations completed'
    endif

    ! GCM
    if(GCM_option%GCMType > 0 ) then
        if (MPI_Infor%rank == 0) write(*,*)'GCM_Main: Start Proj calculations'
        if (MPI_Infor%rank == 0) call GCM_Main
        if (MPI_Infor%rank == 0) write(*,*)'GCM_Main: Proj calculations completed'
    endif

    close(outputfile%u_config)

    if (MPI_Infor%rank == 0) write(*,*) " ...PROGRAM END!..."
    if (MPI_Infor%rank == 0) then 
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
    end if 
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_Finalize(ierr)
END PROGRAM MR_CDFT