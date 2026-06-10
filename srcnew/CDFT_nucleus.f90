!==============================================================================!
! MODULE Nucleus                                                               !
!                                                                              !
! This module contains the variables and routines to set the target nucleus.   !                                                                      !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Nucleus
use Constants, only: r64
use Globals, only : nucleus_attributes, outputfile
implicit none

contains

subroutine set_nucleus_attributes(ifPrint)
!-------------------------------------------------------------------------------------!
!Purpose : Set nucleus name, proton number, neutron number of the target nucleus.     !
!                                                                                     !
!-------------------------------------------------------------------------------------!
    use Globals, only: input_par
    logical,intent(in),optional :: ifPrint
    integer, parameter :: MAXZ = 150
    character(len=2*MAXZ +2) :: T
    integer :: np
    integer :: proton_number
    character(len=3) :: nucleus_name
    ! set input parameters
    nucleus_attributes%name = input_par%nucleus_name
    nucleus_attributes%mass_number_int = input_par%nucleus_mass_number
    ! set other nucleus attributes
    T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
    T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
    T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
    T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
    T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
    T(201:240) = 'FmMdNoLrRfHaSgNsHsMr10111213141516171819'
    T(241:280) = '2021222324252627282930313233343536373839'
    T(281:282) = '40'
    proton_number = nucleus_attributes%proton_number
    nucleus_name = nucleus_attributes%name
    if (proton_number /= 0 ) then                                                 
        if (proton_number.lt.0.or.proton_number.gt.MAXZ) stop '[Nucleus] neutron_number wrong!'
        nucleus_name = T(2*proton_number+1 : 2*proton_number+2)                                                                                            
    else if(nucleus_name /= '') then                                                              
        do np = 0,maxz                                           
            if (nucleus_name.eq.T(2*np+1 : 2*np+2)) then                            
                proton_number = np                                                
                if (proton_number.eq.139) proton_number = 140                                                                                
            endif                                                 
        end do  
    else
        stop '[Nucleus] proton_number or nucleus_name should be given!'                                                                              
    endif
    nucleus_attributes%name = nucleus_name
    nucleus_attributes%proton_number = proton_number
    nucleus_attributes%mass_number = nucleus_attributes%mass_number_int  ! int -> float 
    nucleus_attributes%neutron_number = nucleus_attributes%mass_number_int - nucleus_attributes%proton_number
    
    if(ifPrint) call printNucleusAttribute
    contains
    subroutine printNucleusAttribute
        character(len=*),parameter :: format1='(a20,5x,a5)', &
                                      format2='(a20,5x,f6.2)'
        write(outputfile%u_config,outputfile%format_title) 'set_nucleus_attributes'
        write(outputfile%u_config,format1) 'Name:', nucleus_attributes%name
        write(outputfile%u_config,format2) 'A(mass number):', nucleus_attributes%mass_number
        write(outputfile%u_config,format2) 'Z(proton number):', nucleus_attributes%proton_number
        write(outputfile%u_config,format2) 'N(neutron number):', nucleus_attributes%neutron_number
    end subroutine printNucleusAttribute 
end subroutine set_nucleus_attributes



END MODULE Nucleus