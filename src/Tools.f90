!==============================================================================!
! MODULE Tools                                                                  !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Tools
    implicit none

    interface adjust_left
        module procedure adjust_left_char
        module procedure adjust_left_int
    end interface

    contains
    function find_file(environment_name,file_name)
        !--------------------------------------------------------------------------------------!
        !function find_file:                                                                   ! 
        !returns directory that find file_name from input environment variable(directory_name).!
        !--------------------------------------------------------------------------------------!
        character(500) :: find_file
        character(*) :: file_name,environment_name
        character(1000) :: path
        integer :: istart,iend,i
        logical :: isthere,last_chance
        call GETENV(trim(environment_name),path)
            path=adjustl(path)
            istart = 1
            do while (.true.)
                i = istart
                do while (path(i:i).ne.':')
                    iend = i
                    i = i+ 1
                    if (path(i:i) == ' ') then
                        i = 0
                    exit
                    end if
                end do
                inquire(file=path(istart:iend)//'/'//trim(file_name),exist=isthere)
                if (isthere) then
                    find_file = path(istart:iend)//'/'
                    return
                else if (i == 0) then
                    inquire(file='./'//file_name,exist=last_chance)
                    if (last_chance) then
                        find_file = './'
                        return
                    else
                        print*,'FILE NOT FOUND: ',trim(adjustl(file_name))
                        print*,'CHECK ENVIRONMENT VARIABLE: ',trim(adjustl(environment_name))
                        stop
                    end if
                end if
                istart = iend+2
            end do
    end function find_file

    function int2str(i) result(str)
        implicit none
        integer, intent(in) :: i
        character(len=:), allocatable :: str
        character(len=32) :: buf
        write(buf,'(I0)') i
        str = trim(buf)
    end function int2str

    function real2str(x) result(str)
        implicit none
        real, intent(in) :: x
        character(len=:), allocatable :: str
        character(len=64) :: buf
        write(buf,'(G0)') x
        str = trim(buf)
    end function real2str

    FUNCTION adjust_left_char(s, length) RESULT(str)
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: s
        INTEGER, INTENT(IN) :: length
        CHARACTER(LEN=length) :: str
        INTEGER :: n
        n = LEN_TRIM(s)
        IF (n >= length) THEN
            str = s(1:length)
        ELSE
            str = ADJUSTL(s) // REPEAT(' ', length - n)
        END IF
    END FUNCTION

    function adjust_left_int(i, length) result(str)
        implicit none
        integer, intent(in) :: i, length
        character(len=length) :: str
        character(len=64) :: tmp

        write(tmp,'(I0)') i
        str = adjust_left_char(tmp, length)
    end function

    logical function file_exists(path)
        character(len=*), intent(in) :: path
        integer :: ios

        inquire(file=trim(path), exist=file_exists, iostat=ios)
        if (ios /= 0) file_exists = .false.
    end function file_exists

    subroutine make_directory(dir_path)
        character(len=*), intent(in) :: dir_path
        integer :: status
        logical :: dir_exists

        ! Check if directory exists (Modern Fortran 2003+)
        ! inquire(directory=trim(dir_path), exist=dir_exists)
        inquire(file=trim(dir_path)//'/.', exist=dir_exists)
        if (.not. dir_exists) then
            write(*,*) "Info: Creating directory: ", trim(dir_path)
            ! #ifdef _WIN32
            !     ! Windows implementation: Use 'cmd /c mkdir'
            !     ! cmd /c ensures the command is executed within the shell environment
            !     call execute_command_line("cmd /c mkdir " // trim(dir_path), exitstat=status)
            ! #else
                ! Unix/Linux/macOS implementation: Use 'mkdir -p'
                ! The -p flag creates parent directories and ignores existing ones
                call execute_command_line("mkdir -p " // trim(dir_path), exitstat=status)
            ! #endif
            if (status /= 0) then
                write(*,*) "Error: Failed to create directory. Exit Status:", status
            end if
        else
            ! Directory already exists, no action needed
            ! write(*,*) "Info: Directory exists:", trim(dir_path)
        end if
    end subroutine make_directory

END MODULE