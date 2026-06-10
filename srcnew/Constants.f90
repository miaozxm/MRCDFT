!==============================================================================!
! MODULE Constants                                                             !
!                                                                              !
! This module contains the variables related to the physical, mathematical and !
! numerical constants.                                                         !
!==============================================================================!
MODULE Constants

use Iso_fortran_env

implicit none
public
!!! Definition of input/output units (for portability)
integer, parameter :: u_input  = input_unit,  &
                      u_output = output_unit, &
                      u_start  = u_input + u_output, &!the beginning of the unit of files
                      u_config = u_start + 999 ! unit of config.out
! outputfile path
character(len=*), parameter :: OUTPUT_PATH = './output/'

!!! Definition of kind parameters (for portability)
integer, parameter :: i8  = int8,   & ! integer  8 bits
                      i16 = int16,  & ! integer 16 bits 
                      i32 = int32,  & ! integer 32 bits 
                      i64 = int64,  & ! integer 64 bits 
                      r32 = real32, & ! real 32 bits (single precision)
                      r64 = real64    ! real 64 bits (double precision)

!!! Definition of simple names for numerical values
real(r64) ::    zero  = 0.0d0, &
                one   = 1.0d0, &
                two   = 2.0d0, &
                half  = 0.5d0, &
                third = 0.333333333333333333d0
                        


!!! Definition of physical constants
real(r64) ::    pi  = 3.141592653589793d0, & ! 3.14159265...
                hbc = 197.328284d0,        & ! hbar*c in Mev.fm
                radius_r0 = 1.2,         & ! radius factor
                alphi = 137.03602        !& 

!!! Constants used in some subroutine or functions(to declare the length of arrays)
! maximal number for GFV                                                      
integer, parameter :: igfv   = 100, &
                      igfvbc = 100                                            
! number of gauss-meshpoints
integer, parameter :: ngh   = 32,&
                      ngl   = 16,&
                      nghl  = ngh * ngl 

! number for cylindrical harmonic oscillator basis
integer, parameter :: nt_max = 2023, &! ntx! maximal number of all levels for protons or neutrons(fermions)
                      ntb_max = 132, & !nox !maximal number of all levels for bosons
                      nb_max = 42, &! nbx! maximal number of nb(nb is the number of K-parity-blocks)
                      nz_max = 21, &!nzx ! maximal nz-quantum number for fermions
                      nr_max = 10, &!nrx ! maximal nr-quantum number for fermions
                      ml_max = 21, &!mlx ! maximal ml-quantum number for fermions 
                      nz_max_boson= 20, &!nzbx! maximal nz-quantum number for bosons  
                      nr_max_boson= 10, &!nrbx! maximal nr-quantum number for bosons
                      nfx = 242, & ! maximal dimension F of one k-block
                      ngx = 264, & ! maximal dimension G of one k-block
                      nhx = nfx  + ngx, &
                      nfgx = nfx * ngx

! number for spherical harmonic oscillator basis
integer, parameter ::   nt3x = nt_max*3,& ! dimension of all |n_r l j m_j > quantum number 
                        nh2x  = 3000
integer, parameter ::   nr_max_sph = 16, &
                        l_max_sph = 17

! maximal Dirac states number
integer, parameter :: nkx = 1000 !2771 

! fixed text
character(len=1) :: tp(2) = ['+','-'], &
                    tis(2) = ['n','p'], &
                    tl(0:20) = ['s','p','d','f','g','h','i','j','k','l','m','n','o','P','q','r','S','t','u','v','w']
character(len=8) :: tit(2) = ['Neutron:','Proton: ']

! ! broyden's mixing method
! integer, parameter :: nn = 6*nghl+3 ! 6 = vps_n + vps_p + vms_n + vms_p + delq_n + delq_p; 3 = c1x + c2x + c3x
! integer, parameter :: mm = 7 ! M: mix M previous iterations

!!!Option
integer :: icou = 1, & ! Coulomb-field: not at all (0), direct term (1), plus exchange (2), ususally set to 1.
           icm  = 0, & ! Center of mass:  not the usual sense of center of mass correction, make sure icm=0 is used
           itx  = 2    ! 1: only neutrons;  2: neutrons and protons                 

           
!!! RHB parameters           
! maximal oscillator quantum number for fermions
integer, parameter :: N0FX  = 14, &
                      nxx   = N0FX/2


integer, parameter :: NNNX  = (N0FX+1)*(N0FX+1)
integer, parameter :: NBX   = nb_max, &
                      NB2X  = NBX + NBX
integer, parameter :: NHHX  = nhx * nhx

! RHB Matrix dimension
integer,parameter :: NHBX   = nhx + nhx, &
                     NHBQX  = NHBX*NHBX 

! working space 
integer, parameter :: MVX1  = ( (nxx+1)**4+(nxx+1)**2)/2, &
                      MVX2a = nxx*(nxx+1)*(2*nxx+1), &
                      MVX2  = MVX2a*(3*nxx**2+3*nxx-1)/15, &
                      MVX3  = nxx*(nxx+1)*(2*nxx+1)/2, &
                      MVX4  = nxx**2*(nxx+1)**2/2, &
                      MVX5  = nxx*(nxx+1)/2, &
                      MVX   = MVX1+MVX2+MVX3+MVX4+MVX5

integer, parameter :: MSUM1 = nxx*(nxx+1)/2, &
                      MSUM2 = nxx*(nxx+1)*(2*nxx+1)/6, &
                      MSUM3 = nxx**2*(nxx+1)**2/4, &
                      MSUM4 = nxx*(nxx+1)*(2*nxx+1)*(3*nxx**2+3*nxx-1)/30, &
                      MVTX1 = 6+7*nxx+4*MSUM4+16*MSUM3+27*MSUM2+22*MSUM1, &
                      MVTX2 = nxx+4*MSUM4+8*MSUM3+9*MSUM2+5*MSUM1, &
                      MVTX  = MVTX1+MVTX2

!***************** Proj parameters *****************
! number of meshpoints
integer, parameter ::   ngr = 16, &
                        ntheta = 16, &
                        nphi = 24
integer, parameter :: Jmax_max = 10

real(r64) :: eps_occupation = 5.E-9
END MODULE Constants
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
