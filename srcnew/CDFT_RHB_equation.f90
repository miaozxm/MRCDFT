!==============================================================================!
! MODULE                                                                  !
!                                                                              !
! This module solves the RHB in cylindrical oscillator basis                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
Module RHBEquation
use Constants, only: r64
use Globals, only: RHB,RHB_pairing
implicit none
contains
subroutine solve_RHB_equation(ifPrint)
    !--------------------------------------------------------------------------
    !   solves the RHB-Equation in cylindrical oscillator basis 
    !   
    !  RHB Matrix:
    !   ( A - lambba        B       Delta_ff        0    )
    !   (   B^T       -C - lambda      0            0    )
    !   ( Delta_ff       0        -A + lambda      -B    )
    !   (    0           0       -B^T        -C + lambda )
    ! 
    !    ( h_D - lambda    Delta~   )
    !  = (                          )
    !    ( Delta~     -h_D + lambda )
    !  where 
    !         ( A     B  )                              (Delta_ff   Delta_fg)
    !   h_D = (          )  Delta~ = keep Delta_ff of   (                   )
    !         ( B     -C )                              (Delta_gf   Delta_gg)
    !--------------------------------------------------------------------------
    use MathMethods, only: sdiag
    use Constants, only: zero,half,one,NHBQX,NBX, NHBX
    use Globals, only: dirac, pairing,BS,iteration
    use Broyden, only: set_initial_matrix_elements_RHB

    logical,intent(in),optional :: ifPrint
    integer :: it,lit,klp,kla, ib,mul,nf,ng,nh,nhb,m,n2,n1,k,n
    real(r64) :: dl,al,xh,xl,snold,sn,v2,alold
    real(r64), dimension(NHBQX) :: hb
    real(r64), dimension(NHBX) :: e,ez
    real(r64) :: epsi = 1.d-8
    integer :: maxl = 200

    ! the initial RHB matrix elements are calculated in the first iteration
    ! the subsequent matrix elements are calculated by mix_matrix_elements_RHB.
    if(iteration%ii .eq. 1) call set_initial_matrix_elements_RHB

    do it=1,2
        dl = 100.d0
        al = pairing%ala(it)
        xh = pairing%ala(it) + dl
        xl = pairing%ala(it) - dl
        do lit = 1,maxl ! loop over lambda-iteration
            snold = sn
            sn  = zero
            klp = 0 ! record the number of particles states(all blocks, positive states )
            kla = 0 ! record the number of anti-particles states(all blocks, negative states )
            do ib = 1,BS%HO_cyl%nb ! loop over differnt blocks
                mul  = BS%HO_cyl%mb(ib)
                nf   = BS%HO_cyl%id(ib,1)
                ng   = BS%HO_cyl%id(ib,2)
                nh   = nf + ng
                nhb  = nh + nh
                m    = ib + (it-1)*NBX
                !! calculation of the RHB-Matrix:
                ! lower triangular of RHB Matrix
                do n2 = 1,nh
                    do n1 = n2,nh
                        hb(   n1+(   n2-1)*nhb) =  dirac%hh(n1+(n2-1)*nh,m) ! lower triangular of  h_D
                        hb(nh+n1+(nh+n2-1)*nhb) = -dirac%hh(n1+(n2-1)*nh,m) ! lower triangular of -h_D
                        hb(nh+n1+(   n2-1)*nhb) =  RHB_pairing%delta(n1+(n2-1)*nh,m) ! lower triangular of Delta_ff
                        hb(nh+n2+(   n1-1)*nhb) =  RHB_pairing%delta(n2+(n1-1)*nh,m) ! upper triangular of Delta_ff
                    enddo
                    hb(   n2+(   n2-1)*nhb) =  hb(n2+(n2-1)*nhb) - al ! h_D - lambda
                    hb(nh+n2+(nh+n2-1)*nhb) = -hb(n2+(n2-1)*nhb) ! h_D + lambda
                enddo
                !! Diagonalization of this symmetric matrix 
                ! call the subroutine with the parameter '+1,' it makes the energies in ascending order
                call sdiag(nhb,nhb,hb,e,hb,ez,+1)

                !! store eigenvalues and wave functions
                ! particles
                RHB%ka(ib,it) = klp ! begining of positive states in ib block
                do k = 1,nf  ! nf: number of positive energy state 
                    klp = klp + 1
                    RHB%equ(klp,it) = e(nh+k)
                    do n = 1,nhb
                        RHB%fguv(n,klp,it) = hb(n+(nh+k-1)*nhb)
                    enddo
                    ! 
                    v2 = zero
                    do n = 1,nh
                        v2 = v2 + RHB%fguv(nh+n,klp,it)**2
                    enddo
                    if (v2.lt.zero) v2 = zero
                    if (v2.gt.one)  v2 = one
                    sn = sn + v2*mul
                enddo
                RHB%kd(ib,it) = klp - RHB%ka(ib,it)
                ! anti-particles
                RHB%ka(ib,it+2) = kla
                do k = 1,ng ! ng: number of negative energy state   
                    kla = kla + 1
                    RHB%equ(kla,it+2) = e(ng-k+1) 
                    do n = 1,nhb
                        RHB%fguv(n,kla,it+2) = hb(n+(ng-k)*nhb)
                    enddo
                    v2 = zero
                    do n = 1,nh  ! no-sea approximation
                        v2 = v2 + RHB%fguv(nh+n,kla,it+2)**2
                    enddo ! no-sea approximation
                    sn = sn + v2*mul
                enddo
                RHB%kd(ib,it+2) = kla - RHB%ka(ib,it+2)
            enddo
            ! calculate new chemical potential( al ): change the value of al
            call calculate_new_chemical_potential
            if (abs(al-alold).lt.epsi) goto 30
        enddo
        stop '[solve_RHB_equation]: Lambda-Iteration interupted !'
        
        ! Lambda-Iteration success!
        30 pairing%ala(it) = al

        if(ifPrint) call printRHB
    enddo
    
    contains

    subroutine calculate_new_chemical_potential
        use Globals, only: nucleus_attributes
        real(r64) :: particle_number,dd,dn
        particle_number = nucleus_attributes%neutron_number
        if(it.eq.2) particle_number = nucleus_attributes%proton_number
        if (lit.gt.1) dd = (sn - snold)/(al - alold)
        ! calculation of a new lambda-value
        alold = al
        dn    = sn - particle_number
        if (dn.lt.zero) then
            xl = al
        else
            xh = al
        endif
        if (lit.eq.1) then
            if(dabs(dn).le.0.1d0) then
                al = al - dn
            else
                al = al - 0.1d0*sign(one,dn)
            endif
        else
            ! secant method
            if (dd.eq.zero) dd = 1.d-20
            al = al - dn/dd
            if (al.lt.xl.or.al.gt.xh) then
                ! bisection
                al = half*(xl+xh)
            endif
        endif
        return
    end subroutine calculate_new_chemical_potential

    subroutine printRHB
        use Globals,only: outputfile
        character(len=*), parameter :: format1= "(i4,a,i4,3f13.8)"
        write(outputfile%u_outputf,*) '*************************BEGIN solve_RHB_equation ********************'                                                          
        write(outputfile%u_outputf,format1) lit,'. Lambda-Iteration successful:',it,al,sn
        write(outputfile%u_outputf,*) '*************************END solve_RHB_equation**********************'  
    end subroutine printRHB
end subroutine solve_RHB_equation


end Module RHBEquation