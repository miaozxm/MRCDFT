!====================================================================================!
! MODULE Broyden                                                                      !
!                                                                                    !
! This module contains subroutines related to  
!                                                                                    !
!                                                                                    !
! List of subroutines and functions:                                                 !
! - subroutine                                                                       !
!====================================================================================!
Module Broyden

use Constants, only: r64,zero,half,ngh,ngl,nghl,NBX
use Globals, only: fields,constraint
implicit none
integer :: nn ! dimension of elements to be mixed
integer :: mm ! M: mix M previous iterations
real(r64),dimension(:),allocatable :: vin 
real(r64),dimension(:),allocatable :: vou
real(r64),dimension(:,:),allocatable :: df ! \Delta F, but df(i,inex) sotre latest F 
real(r64),dimension(:,:),allocatable :: dv ! \Delta V, but dv(i,inex) sotre latest V_{in}

contains
subroutine set_broyden_parameters
    use Constants, only:  MVTX, MVX
    use Globals, only: option

    if (option%eqType .eq. 0) then ! dirac 
        nn = 6*nghl+3 ! 6 = vps_n + vps_p + vms_n + vms_p + delq_n + delq_p; 3 = c1x + c2x + c3x
        mm = 7 
    else ! RHB
        nn = 2*MVTX+2*MVX+2+3
        mm = 7
    endif

    allocate(vin(nn))
    allocate(vou(nn))
    allocate(df(nn,mm))
    allocate(dv(nn,mm))
end subroutine set_broyden_parameters

subroutine set_initial_matrix_elements_RHB
    use Constants, only: NHHX,NB2X
    use Globals, only: BS,dirac,RHB_pairing,pairing
    use DiracEquation, only: calculate_all_dirac_matrix
    integer :: ipos,it,ib,nf,ng,nh,m,i2,i1
    ! calculate dirac matrix elements: dirac%hh
    allocate(dirac%hh(NHHX,NB2X))
    call calculate_all_dirac_matrix
    ! store initial matrix elements for broyden mixing
    ipos = 0
    do it = 1,2
        ! store h_D(ff and gg, fg and gf remain unchanged in the iteration) and delta
        do ib = 1,BS%HO_cyl%nb
            nf = BS%HO_cyl%id(ib,1)
            ng = BS%HO_cyl%id(ib,2)
            nh = nf + ng
            m = ib + (it-1)*NBX
            do i2 = 1,nf
                do i1 = i2,nf
                    ipos = ipos + 1
                    vin(ipos) = dirac%hh(i1+(i2-1)*nh,m) ! lower triangular of A (ff)
                enddo
            enddo
            do i2 = nf+1,nh
                do i1 = i2,nh
                    ipos = ipos+1                    
                    vin(ipos) = dirac%hh(i1+(i2-1)*nh,m) ! lower triangular of C (gg)
                enddo
            enddo                    
            do i2 = 1,nf
                do i1 = i2,nf
                    ipos = ipos+1
                    vin(ipos) = RHB_pairing%delta(i1+(i2-1)*nh,m) ! lower triangular of \Delta_ff
                enddo
            enddo           
        enddo
        ipos = ipos+1
        vin(ipos) = pairing%ala(it) ! lambda
    enddo
    ipos = ipos+1
    vin(ipos) = constraint%c1x ! constraint potential
    ipos = ipos + 1
    vin(ipos) = constraint%c2x ! constraint potential Quardrupole
    ipos = ipos + 1
    vin(ipos) = constraint%c3x ! constraint potential Octupole
end subroutine set_initial_matrix_elements_RHB

subroutine mix_matrix_elements_RHB
    use Globals, only: BS,dirac,RHB_pairing,pairing,constraint
    use DiracEquation, only: calculate_all_dirac_matrix
    integer :: ipos,it,ib,nf,ng,nh,m,i2,i1
    ! store new matrix elements(vou)
    call calculate_matrix_elements_RHB
    ! mix matrix elements by broyden's method 
    call broyden_mixing
    ! set the new mixed matrix elements
    ipos = 0
    do it = 1,2
        do ib = 1,BS%HO_cyl%nb
            nf = BS%HO_cyl%id(ib,1)
            ng = BS%HO_cyl%id(ib,2)
            nh = nf + ng
            m = ib + (it-1)*NBX
            do i2 = 1,nf
                do i1 = i2,nf
                    ipos = ipos + 1
                    dirac%hh(i1+(i2-1)*nh,m) = vin(ipos) ! lower triangular of A (ff)
                    dirac%hh(i2+(i1-1)*nh,m) = vin(ipos) ! upper triangular of A (ff)
                enddo
            enddo
            do i2 = nf+1,nh
                do i1 = i2,nh
                    ipos = ipos+1                    
                    dirac%hh(i1+(i2-1)*nh,m) = vin(ipos) ! lower triangular of C (gg)
                    dirac%hh(i2+(i1-1)*nh,m) = vin(ipos) ! upper triangular of C (gg)
                enddo
            enddo                    
            do i2 = 1,nf
                do i1 = i2,nf
                    ipos = ipos+1
                    RHB_pairing%delta(i1+(i2-1)*nh,m) = vin(ipos) ! lower triangular of \Delta_ff
                    RHB_pairing%delta(i2+(i1-1)*nh,m) = vin(ipos) ! upper triangular of \Delta_ff
                enddo
            enddo           
        enddo
        ipos = ipos+1
        pairing%ala(it) = vin(ipos)  ! lambda
    enddo
    ipos = ipos+1
    constraint%c1x = vin(ipos) ! constraint potential
    ipos = ipos + 1
    constraint%c2x = vin(ipos) ! constraint potential Quardrupole
    ipos = ipos + 1
    constraint%c3x = vin(ipos) ! constraint potential Octupole
end subroutine mix_matrix_elements_RHB

subroutine calculate_matrix_elements_RHB
    use Constants, only: half,hbc,zero
    use Globals, only: force,densities,BS,dirac,RHB_pairing,pairing
    use DiracEquation, only: calculate_all_dirac_matrix
    use DeltaField, only: calculate_delta_field
    real(r64),dimension(4) :: glt
    real(r64) :: re,s1,s2,v1,v2,vcn
    integer :: ih,il,i,ipos,it,ib,nf,ng,nh,m,i2,i1
    ! calculate potentials
    do il = 1,ngl
        do ih = 1,ngh
            i = ih + ngh*(il-1)
            ! meson-exchange models
            if (force%option%ipc.eq.0) then
                glt(1) = force%couplg%ggsig * force%ff(i,1,1) * fields%sigma(ih,il) ! sigma
                glt(2) = force%couplg%ggome * force%ff(i,2,1) * fields%omega(ih,il) ! omega
                glt(3) = force%couplg%ggdel * force%ff(i,3,1) * fields%delta(ih,il) ! delta
                glt(4) = force%couplg%ggrho * force%ff(i,4,1) * fields%rho(ih,il)   ! rho
                ! rearangement field
                if (force%option%idd.eq.1) then
                    re =  force%couplg%ggsig * force%ff(i,1,2) * fields%sigma(ih,il) * densities%ro(i,1) &
                        + force%couplg%ggome * force%ff(i,2,2) * fields%omega(ih,il) * densities%ro(i,2) &
                        + force%couplg%ggdel * force%ff(i,3,2) * fields%delta(ih,il) * densities%ro(i,3) &
                        + force%couplg%ggrho * force%ff(i,4,2) * fields%rho(ih,il)   * densities%ro(i,4)                   
                    glt(2) = glt(2) + re
                endif
            ! point-coupling models
            elseif (force%option%ipc.eq.1) then
                glt(1) = force%couplg%ggsig * force%ff(i,1,1) * densities%ro(i,1) ! S
                glt(2) = force%couplg%ggome * force%ff(i,2,1) * densities%ro(i,2) ! V
                glt(3) = force%couplg%ggdel * force%ff(i,3,1) * densities%ro(i,3) ! TS
                glt(4) = force%couplg%ggrho * force%ff(i,4,1) * densities%ro(i,4) ! TV
                ! derivative terms
                glt(1) = glt(1) + force%coupld%ddsig*densities%dro(i,1)
                glt(2) = glt(2) + force%coupld%ddome*densities%dro(i,2)
                glt(3) = glt(3) + force%coupld%dddel*densities%dro(i,3)
                glt(4) = glt(4) + force%coupld%ddrho*densities%dro(i,4)
                ! rearangement field
                if (force%option%inl.eq.1) then
                    glt(1) = glt(1) + half * force%couplg%ggsig * force%ff(i,1,2) * densities%ro(i,1)**2
                    glt(2) = glt(2) + half * force%couplg%ggome * force%ff(i,2,2) * densities%ro(i,2)**2
                    glt(3) = glt(3) + half * force%couplg%ggdel * force%ff(i,3,2) * densities%ro(i,3)**2
                    glt(4) = glt(4) + half * force%couplg%ggrho * force%ff(i,4,2) * densities%ro(i,4)**2
                else if (force%option%inl.eq.0) then
                    re =  half * force%couplg%ggsig * force%ff(i,1,2) * densities%ro(i,1)**2 &
                        + half * force%couplg%ggome * force%ff(i,2,2) * densities%ro(i,2)**2 &
                        + half * force%couplg%ggdel * force%ff(i,3,2) * densities%ro(i,3)**2 &
                        + half * force%couplg%ggrho * force%ff(i,4,2) * densities%ro(i,4)**2       
                    glt(2) = glt(2) + re
                endif
            else
                stop 'calculate_potentials: ipc not properly defined'
            endif

            s1 = glt(1) - glt(3) ! neutron scalar
            s2 = glt(1) + glt(3) ! proton  scalar
            v1 = glt(2) - glt(4) ! neutron vector
            v2 = glt(2) + glt(4) + fields%coulomb(ih,il) ! proton  vector

            fields%vps(ih,il,1) = (v1+s1)*hbc ! V+S for neutron
            fields%vps(ih,il,2) = (v2+s2)*hbc ! V+S for proton
            fields%vms(ih,il,1) = (v1-s1)*hbc ! V-S for neutron
            fields%vms(ih,il,2) = (v2-s2)*hbc ! V-S for proton
            vcn = constraint%c1x*constraint%vc(i,1) + & ! total constraint potentials
                  constraint%c2x*constraint%vc(i,2) + &
                  constraint%c3x*constraint%vc(i,3)
            fields%vpstot(ih,il,1) = fields%vps(ih,il,1) + vcn ! total V+S for neutron
            fields%vpstot(ih,il,2) = fields%vps(ih,il,2) + vcn ! total V+S for porton
            fields%vmstot(ih,il,1) = fields%vms(ih,il,1) + vcn ! total V-S for neutron
            fields%vmstot(ih,il,2) = fields%vms(ih,il,2) + vcn ! total V-S for porton
        enddo
    enddo
    
    ! calculate dirac matrix elements: dirac%hh
    call calculate_all_dirac_matrix

    ! calculate matrix elements \Delta_ff
    call calculate_delta_field
    ! store matix elements for broyden mixing
    ipos = 0
    do it = 1,2
        ! store h_D(ff and gg, fg and gf remain unchanged in the iteration) and delta
        do ib = 1,BS%HO_cyl%nb
            nf = BS%HO_cyl%id(ib,1)
            ng = BS%HO_cyl%id(ib,2)
            nh = nf + ng
            m = ib + (it-1)*NBX
            do i2 = 1,nf
                do i1 = i2,nf
                    ipos = ipos + 1
                    vou(ipos) = dirac%hh(i1+(i2-1)*nh,m) ! lower triangular of A (ff)
                enddo
            enddo
            do i2 = nf+1,nh
                do i1 = i2,nh
                    ipos = ipos+1                    
                    vou(ipos) = dirac%hh(i1+(i2-1)*nh,m) ! lower triangular of C (gg)
                enddo
            enddo                    
            do i2 = 1,nf
                do i1 = i2,nf
                    ipos = ipos+1
                    vou(ipos) = RHB_pairing%delta(i1+(i2-1)*nh,m) ! lower triangular of \Delta_ff
                    if(pairing%del(it) .lt. 1.d-5) then
                        vou(ipos) = zero
                        vin(ipos) = zero
                    endif
                enddo
            enddo           
        enddo
        ipos = ipos+1
        vou(ipos) = pairing%ala(it) ! lambda
    enddo
    ipos = ipos+1
    vou(ipos) = constraint%c1x ! constraint potential
    ipos = ipos + 1
    vou(ipos) = constraint%c2x ! constraint potential Quardrupole
    ipos = ipos + 1
    vou(ipos) = constraint%c3x ! constraint potential Octupole
end subroutine calculate_matrix_elements_RHB

subroutine mix_potentials_DIR(ifPrint)
    !--------------------------------------------------------------------------
    !   Calculation of the potentials at gauss-meshpoints
    ! 
    !-------------------------------------------------------------------------
    use Globals,only: fields,pairing
    logical,intent(in),optional :: ifPrint

    real(r64) :: vcn
    integer :: i,il,ih
    ! store new potentials(vou) and old potentials(vin) 
    call calculate_potentials_DIR
    ! mix potentials by broyden's method
    call broyden_mixing
    ! set the new mixed potentials
    constraint%c1x  = vin(1+nn-3)
    constraint%c2x  = vin(2+nn-3)
    constraint%c3x  = vin(3+nn-3)
    do il = 1, ngl
        do ih = 1,ngh
            i = ih + ngh*(il-1)
            fields%vps(ih,il,1) = vin( i )          ! V+S for neutron
            fields%vps(ih,il,2) = vin( i +   nghl ) ! V+S for proton
            fields%vms(ih,il,1) = vin( i + 2*nghl ) ! V+S for neutron
            fields%vms(ih,il,2) = vin( i + 3*nghl ) ! V+S for proton
            pairing%delq(i,1)= vin( i + 4*nghl )    ! pairing potential for neutron
            pairing%delq(i,2)= vin( i + 5*nghl )    ! pairing potential for proton
            vcn = constraint%c1x*constraint%vc(i,1) + & ! total constraint potentials
                  constraint%c2x*constraint%vc(i,2) + &
                  constraint%c3x*constraint%vc(i,3)
            fields%vpstot(ih,il,1) = fields%vps(ih,il,1) + vcn ! total V+S for neutron
            fields%vpstot(ih,il,2) = fields%vps(ih,il,2) + vcn ! total V+S for porton
            fields%vmstot(ih,il,1) = fields%vms(ih,il,1) + vcn ! total V-S for neutron
            fields%vmstot(ih,il,2) = fields%vms(ih,il,2) + vcn ! total V-S for porton
        enddo
    enddo
    constraint%c1xold = constraint%c1x ! constraint potential coefficient(Dipole)
    constraint%c2xold = constraint%c2x ! constraint potential coefficient(Quardrupole)
    constraint%c3xold = constraint%c3x ! constraint potential coefficient(Octupole)

    if(ifPrint) call print_potentials
    contains
    subroutine print_potentials
        use Field, only: prigh
        use Constants, only: one
        use Globals, only: outputfile
        write(outputfile%u_outputf,*) '*************************BEGIN print_potentials ********************'
        call prigh(outputfile%u_outputf,1,fields%vps(1,1,1),one,'V+S  n')
        call prigh(outputfile%u_outputf,1,fields%vms(1,1,1),one,'V-S  n')
        call prigh(outputfile%u_outputf,1,fields%vps(1,1,2),one,'V+S  p')
        call prigh(outputfile%u_outputf,1,fields%vms(1,1,2),one,'V-S  p')
        call prigh(outputfile%u_outputf,1,fields%vpstot(1,1,1),one,'V+S Total  n')
        call prigh(outputfile%u_outputf,1,fields%vmstot(1,1,1),one,'V-S Total  n')
        call prigh(outputfile%u_outputf,1,fields%vpstot(1,1,2),one,'V+S Total  p')
        call prigh(outputfile%u_outputf,1,fields%vmstot(1,1,2),one,'V-S Total  p')        
        call prigh(outputfile%u_outputf,1,pairing%delq(1,1),one,'Del  n')
        call prigh(outputfile%u_outputf,1,pairing%delq(1,2),one,'Del  p')
        write(outputfile%u_outputf,*) '*************************END print_potentials **********************'
    end subroutine print_potentials
end subroutine mix_potentials_DIR

subroutine calculate_potentials_DIR
    !----------------------------------------------------------------
    !   calculate new potentials 
    !   store new potentials(vou) and old potentials(vin)
    !----------------------------------------------------------------
    use Constants, only: half,hbc
    use Globals, only: force,densities,pairing
    real(r64),dimension(4) :: glt
    real(r64) :: re,s1,s2,v1,v2,vp1,vp2
    integer :: ih,il,i
    do il = 1,ngl
        do ih = 1,ngh
            i = ih + ngh*(il-1)
            ! meson-exchange models
            if (force%option%ipc.eq.0) then
                glt(1) = force%couplg%ggsig * force%ff(i,1,1) * fields%sigma(ih,il) ! sigma
                glt(2) = force%couplg%ggome * force%ff(i,2,1) * fields%omega(ih,il) ! omega
                glt(3) = force%couplg%ggdel * force%ff(i,3,1) * fields%delta(ih,il) ! delta
                glt(4) = force%couplg%ggrho * force%ff(i,4,1) * fields%rho(ih,il)   ! rho
                ! rearangement field
                if (force%option%idd.eq.1) then
                    re =  force%couplg%ggsig * force%ff(i,1,2) * fields%sigma(ih,il) * densities%ro(i,1) &
                        + force%couplg%ggome * force%ff(i,2,2) * fields%omega(ih,il) * densities%ro(i,2) &
                        + force%couplg%ggdel * force%ff(i,3,2) * fields%delta(ih,il) * densities%ro(i,3) &
                        + force%couplg%ggrho * force%ff(i,4,2) * fields%rho(ih,il)   * densities%ro(i,4)                   
                    glt(2) = glt(2) + re
                endif
            ! point-coupling models
            elseif (force%option%ipc.eq.1) then
                glt(1) = force%couplg%ggsig * force%ff(i,1,1) * densities%ro(i,1) ! S
                glt(2) = force%couplg%ggome * force%ff(i,2,1) * densities%ro(i,2) ! V
                glt(3) = force%couplg%ggdel * force%ff(i,3,1) * densities%ro(i,3) ! TS
                glt(4) = force%couplg%ggrho * force%ff(i,4,1) * densities%ro(i,4) ! TV
                ! derivative terms
                glt(1) = glt(1) + force%coupld%ddsig*densities%dro(i,1)
                glt(2) = glt(2) + force%coupld%ddome*densities%dro(i,2)
                glt(3) = glt(3) + force%coupld%dddel*densities%dro(i,3)
                glt(4) = glt(4) + force%coupld%ddrho*densities%dro(i,4)
                ! rearangement field
                if (force%option%inl.eq.1) then
                    glt(1) = glt(1) + half * force%couplg%ggsig * force%ff(i,1,2) * densities%ro(i,1)**2
                    glt(2) = glt(2) + half * force%couplg%ggome * force%ff(i,2,2) * densities%ro(i,2)**2
                    glt(3) = glt(3) + half * force%couplg%ggdel * force%ff(i,3,2) * densities%ro(i,3)**2
                    glt(4) = glt(4) + half * force%couplg%ggrho * force%ff(i,4,2) * densities%ro(i,4)**2
                else if (force%option%inl.eq.0) then
                    re =  half * force%couplg%ggsig * force%ff(i,1,2) * densities%ro(i,1)**2 &
                        + half * force%couplg%ggome * force%ff(i,2,2) * densities%ro(i,2)**2 &
                        + half * force%couplg%ggdel * force%ff(i,3,2) * densities%ro(i,3)**2 &
                        + half * force%couplg%ggrho * force%ff(i,4,2) * densities%ro(i,4)**2       
                    glt(2) = glt(2) + re
                endif
            else
                stop 'calculate_potentials: ipc not properly defined'
            endif

            s1 = glt(1) - glt(3) ! neutron scalar
            s2 = glt(1) + glt(3) ! proton  scalar
            v1 = glt(2) - glt(4) ! neutron vector
            v2 = glt(2) + glt(4) + fields%coulomb(ih,il) ! proton  vector

            ! pairing field
            if (pairing%ide(1).eq.4) then
                vp1 = pairing%vpair(1) * densities%rkapp(i,1) / 2.d0 ! neutron pairing field
                vp2 = pairing%vpair(2) * densities%rkapp(i,2) / 2.d0 ! proton pairing field
            endif

    !-----------------------------------------------------------------------------
    ! store the potentials in 1D arrays for easier handling of broyden's mixing
    ! vin = old pots
    ! vou = new pots
    !-----------------------------------------------------------------------------
            vin( i          ) = fields%vps(ih,il,1) ! V+S for neutron(old)
            vin( i +   nghl ) = fields%vps(ih,il,2) ! V+S for proton(old)
            vin( i + 2*nghl ) = fields%vms(ih,il,1) ! V-S for neutron(old)
            vin( i + 3*nghl ) = fields%vms(ih,il,2) ! V-S for proton(old)
            vou( i          ) = (v1+s1)*hbc ! V+S for neutron(calculated)
            vou( i +   nghl ) = (v2+s2)*hbc ! V+S for proton(calculated)
            vou( i + 2*nghl ) = (v1-s1)*hbc ! V-S for neutron(calculated)
            vou( i + 3*nghl ) = (v2-s2)*hbc ! V-S for proton(calculated)
            if (pairing%ide(1).eq.4) then
                vin( i + 4*nghl ) = pairing%delq(i,1) ! pairing potential for neutron(old)
                vin( i + 5*nghl ) = pairing%delq(i,2) ! pairing potential for proton(old)
                vou( i + 4*nghl ) = vp1 ! pairing potential for neutron(calculated)
                vou( i + 5*nghl ) = vp2 ! pairing potential for porton(calculated)
            else
                vin( i + 4*nghl ) = zero
                vin( i + 5*nghl ) = zero
                vou( i + 4*nghl ) = zero
                vou( i + 5*nghl ) = zero
            endif
        enddo
    enddo
    vin(1+nn-3) = constraint%c1xold ! constraint potential(old)
    vin(2+nn-3) = constraint%c2xold ! constraint potential Quardrupole(old)
    vin(3+nn-3) = constraint%c3xold ! constraint potential Octupole(old)
    vou(1+nn-3) = constraint%c1x ! constraint potential(calculated)
    vou(2+nn-3) = constraint%c2x ! constraint potential Quardrupole(calculated)
    vou(3+nn-3) = constraint%c3x ! constraint potential Octupole(calculated)
end subroutine calculate_potentials_DIR

subroutine broyden_mixing
    !----------------------------------------------------------------
    !                      Broyden's Method
    !
    !     more information:
    !     Broyden's Method in Nuclear Structure Calculations, by
    !     Andrzej Baran, Aurel Bulgac, Michael McNeil Forbes, 
    !     Gaute Hagen, Witold Nazarewicz, Nicolas Schunck, 
    !     and Mario V. Stoitsov
    !     https://doi.org/10.1103/PhysRevC.78.014318
    !     
    !     Note:
    !     Linear mixing is used in the first 3 iterations
    !     to ensure the stability of Broytden's method.
    !----------------------------------------------------------------
    use Globals, only: iteration
    use MathMethods, only: dnrm2,dscal,ddot,dsytrf,dsytri,dcopy
    integer :: i,ii,iuse,ipos,inex,j,info,k
    real(r64) :: dnorm,xmi,gamma

    real(r64) :: bw0 = 0.01d0  ! \omega_0
    real(r64),dimension(mm,mm) :: bbeta ! \beta_{kn}
    real(r64),dimension(mm)    :: bwork ! c_k^m
    real(r64),dimension(nn)    :: curv ! alpha * F^m
    integer,dimension(mm)    :: ibwork

    ii = iteration%ii
    xmi = iteration%xmix
    ! calculate F=V_{out}-V_{in} and sotre in vou
    do i = 1,nn
        vou(i) = vou(i) - vin(i) ! F = V_{out} -V_{in}
    enddo
    ! Maximum interval value of F
    iteration%si  = zero
    do i = 1,nn
        iteration%si = max(iteration%si,abs(vou(i)))
    enddo

    ! linear mixing for 3 steps
    if (mm.eq.0 .or. ii.le.3) then  
         do i = 1, nn
            vin(i) = vin(i) + iteration%xmix*vou(i)
         enddo
    ! broyden mixing
    else
        iuse = min( ii-1 - 1, mm )
        ipos = ii - 2 - ( (ii-3)/mm )*mm
        inex = ii - 1 - ( (ii-2)/mm )*mm

        ! calculate \Delta F and \Delta V
        if( ii .eq. 4 ) then
            do j = 1, mm
                do i = 1, nn
                    df(i,j) = zero !\Delta F
                    dv(i,j) = zero ! \Delta V
                enddo
            enddo
        else
            do i = 1, nn
                df(i,ipos) = vou(i) - df(i,ipos) ! F^{n+1} - F^{n}
                dv(i,ipos) = vin(i) - dv(i,ipos) ! V_{in}^{n+1} - V_{in}^{n}  
            enddo
            dnorm = sqrt( dnrm2(nn,df(1,ipos),1)**2.0d0 )
            call dscal( nn, 1.0d0 / dnorm, df(1,ipos), 1 ) ! \Delta F
            call dscal( nn, 1.0d0 / dnorm, dv(1,ipos), 1 ) ! \Delta V
        endif 

        ! calculate matrix \beta_{kn} and c_k^m
        do i = 1, iuse
            do j = i+1, iuse
               bbeta(i,j) = ddot( nn, df(1,j), 1, df(1,i), 1 )
            enddo
            bbeta(i,i) = 1.d0 + bw0*bw0
        enddo
        call dsytrf( 'U', iuse, bbeta, mm, ibwork, bwork, mm, info )
        if( info .ne. 0 ) stop 'broyden_mixing: info at DSYTRF V+S '        
        call dsytri( 'U', iuse, bbeta, mm, ibwork, bwork, info )
        if( info .ne. 0 ) stop 'broyden_mixing: info at DSYTRI V+S '
        do i = 1, iuse
            do j = i+1, iuse
               bbeta(j,i) = bbeta(i,j) ! Matrix \beta is symmetric
            enddo
            bwork(i) = ddot( nn, df(1,i), 1, vou, 1 ) !  c_k^m =\omega_k (\Delata F_{k})'*F^(m), where \omega_k = 1, m is this iteration.
        enddo
        ! calculate \alpha * F^m, linear part
        do i = 1, nn
            curv(i) = xmi * vou(i) ! \alpha * F^m = \alpha ( V_{out}^m - V_{in}^m )
        enddo
        ! calculate \alpha * F^m - \sum^{m-1}_{n=max(1,m-M)}\omega_n \gamma_{mn}u^n, where \omega_n = 1 
        do i = 1, iuse ! sum n
            gamma = 0.0d0
            do j = 1, iuse ! sum k
               gamma = gamma + bbeta(j,i) * bwork(j) ! \gamma_{mn} =c_k^m \beta_{kn}
            enddo
            do k = 1, nn
               curv(k) = curv(k) - gamma * ( dv(k,i) + xmi * df(k,i) ) ! \alpha *F^m - \gamma_{mn}*u^n, u^n = \alpha \Delta F^n + \Delta V^n 
            enddo
        enddo
        ! store F and V_{in} for next iteration
        call dcopy( nn, vou, 1, df(1,inex), 1 ) ! store F of this iteration to df(1,interation_index)
        call dcopy( nn, vin, 1, dv(1,inex), 1 ) ! sotre V_{in} of this iteration to dv(1,interation_index
        ! V_{in}^{m+1}
        do i = 1, nn
            vin(i) = vin(i) + curv(i)
        enddo
    endif
end subroutine broyden_mixing

end Module Broyden