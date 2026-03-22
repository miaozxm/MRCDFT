!==============================================================================!
! MODULE Occupation                                                            !
!                                                                              !
! This module calculate the occupation.                                        !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE BCS 
use Constants, only: r64,nkx
use Globals, only: pairing,dirac
implicit none

integer :: nx1 ! the number of sigle particle states
real(r64), dimension(nkx) :: ek ! energy of kth single particle state
real(r64), dimension(nkx) :: dk ! pairing gaps
real(r64), dimension(nkx) :: sk ! cutoff weights
real(r64) :: alx ! fermi energy
real(r64) :: tzx ! neutron or proton numbers N
real(r64) :: ecut ! cut-off energy in pairing window
integer :: block_level ! block_level for neutron or proton

contains

subroutine initial_pairing_field(ifPrint)
    !--------------------------------------------------------------------
    !   set the initial `pairing potential`, `fermi energy`, `cutoff energy`
    !   they will be updated after each iteration
    !--------------------------------------------------------------------
    use Globals, only: input_par,nucleus_attributes,nghl,constraint,option
    logical,intent(in),optional :: ifPrint
    integer :: it,ii
    ! set input parameters
    pairing%ide = input_par%pairing_ide
    pairing%dec = input_par%pairing_dec
    pairing%ga  = input_par%pairing_ga
    pairing%del = input_par%pairing_del
    pairing%vpair = input_par%pairing_vpair
    ! set other parameters
    do it = 1,2
        pairing%gg(it) = pairing%ga(it)/nucleus_attributes%mass_number + 1.d-10
        ! initial pairing potential
        if (pairing%ide(it) == 4) then
            do ii =1,nghl
                pairing%delq(ii,it) = pairing%del(it)
            enddo
        endif
        ! only initialize the first loop, and use the calculated value from the previous loop for subsequent loops.
        if(constraint%index==1) then
            ! initial fermi energy
            pairing%ala(it) = -7.0
            ! inital cutoff energy
            pairing%ecut(it) = 5.d0
        endif
    enddo
    ! set block
    pairing%block_level = input_par%block_level
    pairing%block_K = input_par%K
    pairing%block_Pi = input_par%Pi
    if(option%block_type == 1 .and. option%block_method == 1) then
        pairing%allow_block = .True.
    else
        pairing%allow_block = .False.
    endif

    if(ifPrint) call printPairingProperties
    contains
    subroutine printPairingProperties
        use Globals, only: outputfile,BS
        character(len=*), parameter :: format1 = "(a20,5x,2f10.6)", &
                                       format2 = "(a,2(f8.4,2h/A),2f10.6,/)"
        write(outputfile%u_outputf,*) '*************************BEGIN set_pairing_parameters ********************'
        write(outputfile%u_outputf,format1) 'Gap parameter(dec) :', pairing%dec
        write(outputfile%u_outputf,format1) 'Gap parameter(del) :', pairing%del
        write(outputfile%u_outputf,format1) 'Pairing Window     :', (pairing%pwi+7.d0)/BS%HO_cyl%hom,pairing%pwi
        write(outputfile%u_outputf,format2) 'Pairing const.     :', pairing%ga,pairing%gg
        write(outputfile%u_outputf,"(a,/)") '*************************END set_pairing_parameters ********************'
    end subroutine printPairingProperties
end subroutine initial_pairing_field

subroutine calculate_occupation(ifPrint)
    use Constants, only: itx
    logical,intent(in),optional :: ifPrint
    integer :: it

    do it =1,itx ! loop over neutron-proton 
        if (pairing%ide(it).eq.1) then                                         
            stop 'ide=1 not implemented' 
        endif

        !  2: Frozen gap BCS 
        if (pairing%ide(it).eq.2) then
            stop 'ide=2 not implemented' 
        endif

        ! 3: Constant G BCS 
        if (pairing%ide(it).eq.3) then
             stop 'ide=3 not implemented'                                       
        endif
        
        ! 4: Delta force BCS 
        if (pairing%ide(it).eq.4) then
            call delta_force_bcs(it)
        endif
    enddo

    if(ifPrint) call printPairing
    contains
    subroutine printPairing
        use Constants,only: zero
        use Globals,only: outputfile,nucleus_attributes
        integer :: it
        character(len=*), parameter :: format1= "(a,3f30.25)"
        character(len=*), parameter :: format2= "(a,15x,2f30.25)"
        write(outputfile%u_outputf,*) '*************************BEGIN print_pairing ********************'
        write(outputfile%u_outputf,'(/,28x,a,8x,a,9x,a)') 'neutron','proton','total'
        ! Lambda                                                             
        write(outputfile%u_outputf,format1) ' lambda(fermi energy) ',pairing%ala
        ! Delta                                                             
        write(outputfile%u_outputf,format1) ' Delta ...............',pairing%del
        ! trace of kappa                                                    
        write(outputfile%u_outputf,format1) ' spk .................',pairing%spk
        do it = 1,2
            if (pairing%spk(it).ne.0) then
                pairing%ga(it) = nucleus_attributes%mass_number * pairing%del(it) / pairing%spk(it)
            else
                pairing%ga(it) = zero   
            endif
        enddo
        write(outputfile%u_outputf,format1) ' effective G*A .......',pairing%ga(1),pairing%ga(2)
        write(outputfile%u_outputf,format1) ' Cut off param .......',pairing%ecut(1),pairing%ecut(2)
        write(outputfile%u_outputf,*) '*************************END print_pairing **********************'
    end subroutine printPairing
end subroutine calculate_occupation

subroutine delta_force_bcs(it)
    !-------------------------------------------------------------------------------------------------------------------------------------------------------------!
    !   Input:                                                                    
    !        it:   1 for neutron; 2 for proton                                
    !-------------------------------------------------------------------------------------------------------------------------------------------------------------
    !   Formulas:
    !        1. \Delta_k = \int \psi^{\dagger}_{k}(\boldsymbol{r}) \Delta_{\tau}(\boldsymbol{r}) \psi_{k}(\boldsymbol{r}) d\boldsymbol{r}
    !        2. f_{k}=\frac{1}{1+\exp \left[\left(\epsilon_k-\epsilon_F-\Delta E\right) / (\Delta E/10)\right]}
    !        3. v_k^2=\frac{1}{2}\left(1-\frac{\epsilon_k-\epsilon_F}{\sqrt{\left(\epsilon_k-\epsilon_F\right)^{2}+(f_k\Delta_k)^{2}}}\right), \quad v_k^2+u_k^2=1
    !        4. \sum_{k}v_k^2 = N
    !        5. 2\sum_{k}f_k = N + 1.65N^{2/3}
    !------------------------------------------------------------------------------------------------------------------------------------------------------------     
    !   Variables:
    !        el = \epsilon_k -\epsilon_F 
    !        dd = f_k \Delta_k
    !        ei = \frac{1}{2} \frac{1}{\sqrt{(\epsilon_k -\epsilon_F)^2+(f_k\Delta_k)^2}} 
    !        v2 = v_k^2 = \frac{1}{2}(1 - \frac{\epsilon_k - \epsilon _F}{\sqrt{(\epsilon_k -\epsilon_F)^2+(f_k\Delta_k)^2}}) 
    !        dev2 = \sum_{k} v_k^2f_k\Delta_k
    !        fv2  = \sum_{k} f_k v_k^2
    !        uv = 2 u_k v_k = \frac{f_k\Delta_k}{\sqrt{(\epsilon_k -\epsilon_F)^2+(f_k\Delta_k)^2}}
    !        se = -2\sum _{k} f_k u_k v_k \Delta_k
    !        sp = 2\sum_{k}f_k u_k v_k
    !------------------------------------------------------------------------------------------------------------------------------------------------------------ 
    !   Outputs:
    !        pairing gaps:                de(k,it)  = \Delta_k
    !        occupation probability:      vv(k,it)  = 2v_k^2
    !        cutoff weights:              skk(k,it) = f_k
    !        fermi energy:                ala(it)   = \epsilon_F
    !        cutoff energy:               ecut(it)  = \Delta E
    !        particle number fluctuation: disp(it)  = <(\Delta N) ^2> =  <N^2>-<N>^2 = 4\sum_k u^2_k v^2_k
    !        trace of tensor kappa:       spk(it)   = \sum_{k} f_k u_k v_k
    !        average of gaps:             dev2(it)  = <\Delta>\equiv \frac{\sum_{k} f_{k} v_{k}^{2} \Delta_{k}}{\sum_{k} f_{k} v_{k}^{2}}
    !        average of gaps2:            del(it)   = \Delta^{uv}\equiv \frac{\sum_{k} f_{k} u_kv_{k} \Delta_{k}}{\sum_{k} f_{k} u_kv_{k}}
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------! 
    use Constants, only: half,one,zero
    integer, intent(in) :: it
    integer :: k
    real(r64) :: se,sp,dev2,fv2,disp,el,dd,ei,v2,uv

    ! caculate pairing gaps(store in pairing%de)
    call calculate_pairing_gap(it) !calculate \Delta_k

    ! set Module Occupation's global variables
    call set_occupation_global_variable(it)

    ! calculate cutoff weights (caculate sk())
    call calculate_cutoff_weights

    ! calculate new fermi surface (update alx)
    call calculate_new_fermi_energy
    ! calculate new cutoff energy (update ecut)
    call calculate_new_cutoff_energy

    !----------------------------------------------------------------
    ! caculate and store pairing variable
    !----------------------------------------------------------------
    se      = 0.0d0
    sp      = 0.0d0
    dev2    = 0.0d0
    fv2     = 0.0d0
    disp    = 0.0d0
    do k = 1, nx1
        if( isBlocked(k)) then
            ! pairing%vv(k,it) = 1.0d0
            ! cycle
            v2 = 0.5d0
            goto 100
        endif
        el = ek(k) - alx
        dd = sk(k)*dk(k)
        ei = half/sqrt(el**2 + dd**2)
        v2 = half - el*ei ! v_k^2
        dev2 = dev2 + dd*v2
        fv2  = fv2 + sk(k)*v2
        uv = 2*dd*ei ! 2*u_k*v_k
        se = se - dd*uv
        sp = sp + sk(k)*uv
        disp = disp + 4.d0*v2*(one-v2)
        if(v2.gt.1.0d0) v2 = one
        if(v2.lt.0.0d0) v2 = zero

    100 pairing%de(k,it) = dk(k)
        pairing%vv(k,it) = 2*v2
        pairing%skk(k,it)= sk(k)
    enddo

    ! pairing%epair(it)= 0.5d0*se
    pairing%ala(it)  = alx  ! fermi energy
    pairing%ecut(it) = ecut ! cut-off energy
    pairing%disp(it) = disp ! fluctuations of particle number: 4*sum_(k>0) u^2_k * v^2_k
    pairing%spk(it)  = 0.5d0 * sp ! trace of kappa: \sum_{k} f_k u_k v_k
    pairing%dev2(it) = dev2/fv2 !average of single-particle pairing gaps
    pairing%del(it)  = -se/sp ! average of single-particle pairing gaps : \Delta^{uv}\equiv \frac{\sum_{k} f_{k} u_kv_{k} \Delta_{k}}{\sum_{k} f_{k} u_kv_{k}}

end subroutine delta_force_bcs

subroutine calculate_pairing_gap(it)
    !----------------------------------------------------------------------
    ! For delta force: calculate the pairing gaps from pairing tensor.
    ! 
    ! Pairing gaps is calculated using the following formula:
    !   \Delta_{k} = \int \psi^{\dagger}_{k}(\boldsymbol{r}) \Delta_{\tau}(\boldsymbol{r}) \psi_{k}(\boldsymbol{r}) d\boldsymbol{r}
    ! where \psi_{k}(\boldsymbol{r}) is a single particle state, \Delta_{\tau}(\boldsymbol{r}) is the pairing potential.
    ! Beacuse, the single particle state: \psi_{k}(\boldsymbol{r}) is expanded in the oscillator basis,
    ! actually,the pairing gaps is calculated using the following formula:
    !       \Delta_{k} = \sum_{n n'} \Delta_{n n'}^{\tau}f^{(k)}_{n}f^{(k)}_{n'} + \sum_{n n'} \Delta_{n n'}^{\tau}g^{(k)}_{n}g^{(k)}_{n'}
    ! where \Delta_{n n'}=  \int  \Psi_{n} \Psi_{n'} \Delta_{\tau}(\boldsymbol{r})  d\boldsymbol{r}
    ! and \Psi_{\alpha} is one of the oscillator basis, 
    ! and f^{(k)}_{n} is the coefficients for expanding the single-particle state's lager component in the oscillator basis
    ! and g^{(k)}_{n} is the coefficients for expanding the single-particle state's small component in the oscillator basis.
    !
    !--------------------------------------------------------------------------------------------------------------------
    ! output: 
    !       pairing%de(k,it)    gairing gaps: \Delta_{k} 
    !----------------------------------------------------------------------
    use Globals, only: BS,iteration
    use Constants, only:nhx,nb_max
    use DiracEquation, only: pot
    integer, intent(in) :: it

    integer :: ib,nf,ng,nh,i0f,i0g,i2,i1,k1,k2,k,n2,n1
    real(r64),dimension(:,:,:), allocatable :: del12
    real(r64) :: gapmax,ss
    allocate(del12(nhx*nhx,nb_max,2))

    do ib = 1,BS%HO_cyl%nb ! number of K-parity-blocks 
        nf  = BS%HO_cyl%id(ib,1)   ! dimension large components of block b                                                                                   
        ng  = BS%HO_cyl%id(ib,2)   ! dimension small components of block b                                           
        nh  = nf + ng                                               
        i0f = BS%HO_cyl%ia(ib,1)   ! begin of the large components of block b is ia(b,1)+1                                            
        i0g = BS%HO_cyl%ia(ib,2)   ! begin of the small components of block b is ia(b,2)+1
        
        ! initialize to zero
        do i2=1,nh
            do i1=1,nh
                del12(i1+(i2-1)*nh,ib,it)=0d0
            enddo
        enddo

        ! calculate \Delta_{\alpha \alpha'}^{\tau} = \int  \Psi_{\alpha} \Psi_{\alpha'} \Delta_{\tau}(\boldsymbol{r})  d\boldsymbol{r} 
        ! where \Psi_{\alpha} is one of the oscillator basis.
        ! Note that only the upper triangular matrix elements are computed and stored in del12 by subroutine pot.
        call pot(i0f,nf,nh,pairing%delq(1,it),del12(1,ib,it)) ! delq = \Delta_{tau}(\boldsymbol(r))
        call pot(i0g,ng,nh,pairing%delq(1,it),del12(nf+1+nf*nh,ib,it))

        ! symmetrize del12 because only the upper triangular matrix elements are computed by subroutine pot
        do i2 = 1,nh
            do i1 = i2+1,nh
                del12(i2+(i1-1)*nh,ib,it) = del12(i1+(i2-1)*nh,ib,it)
            enddo
        enddo
    enddo

    gapmax = 0.0d0
    do ib=1,BS%HO_cyl%nb
        k1 = dirac%ka(ib,it) + 1 ! begining position of energy state of k block
        k2 = dirac%ka(ib,it) + dirac%kd(ib,it) ! dimension of energy state of k block
        do k=k1,k2
            nf  = BS%HO_cyl%id(ib,1)   ! dimension large components of block b                                                                                   
            ng  = BS%HO_cyl%id(ib,2)   ! dimension small components of block b
            nh = nf + ng
            ss = 0.0d0
            ! large components' contribution 
            do n2=1,nf
                do n1=1,nf
                    ss = ss + dirac%fg(n1,k,it)*del12((n2-1)*nh+n1,ib,it)*dirac%fg(n2,k,it)
                enddo
            enddo
            ! small components' contribution
            ! do n2=1,ng
            !     do n1=1,ng
            !         ss = ss + dirac%fg(n1+nf,k,it)*del12((n2+nf-1)*nh+nf+n1,ib,it)*dirac%fg(n2+nf,k,it)
            !     enddo
            ! enddo
            
            pairing%de(k,it) = ss+0.01d0 !to avoid pairing collapse (plus a small number)
            if(iteration%ii.eq.1) pairing%de(k,it)=pairing%del(it)
            gapmax = max(gapmax, pairing%de(k,it))
        enddo

    !     if(gapmax.lt.0.1d0) then
    !         write(*,*) '...... pairing collapse .....',it
    !         do k=k1,k2
    !             pairing%de(k,it) = pairing%de(k,it)+0.1d0
    !         enddo
    !    endif 
    enddo
    return
end subroutine calculate_pairing_gap

subroutine set_occupation_global_variable(it)
    !----------------------------------------------------------------
    ! set Module Occupation's global variables 
    !----------------------------------------------------------------
    use Globals, only: constraint,nucleus_attributes
    integer, intent(in) :: it
    integer :: k
    real(r64) :: vk1,uk1
    ! block level for Odd-N/Z nuclei
    if(it==1) then
        block_level = pairing%block_level(1)
    else
        block_level = pairing%block_level(2)
    endif

    ! neutron or proton numbers(store tzx)
    if(it==1) then
        tzx=nucleus_attributes%neutron_number
    else if(it==2) then
        tzx=nucleus_attributes%proton_number
    endif

    nx1 = dirac%nk(it) ! the number of sigle particle states
    do k = 1,nx1 ! the number of sigle particle states
        ! if(isBlocked(k)) cycle
        ek(k) = dirac%ee(k,it) ! energy of kth single particle state
        vk1   = sqrt(pairing%vv(k,it)/2.d0) ! v_k
        uk1   = sqrt(1.d0 -pairing%vv(k,it)/2.d0) ! u_k
        dk(k) = pairing%de(k,it) + 4.d0*constraint%clam2(constraint%index)*vk1*uk1
    enddo
    alx = pairing%ala(it) ! fermi energy
    ecut = pairing%ecut(it) ! cut-off energy

end subroutine set_occupation_global_variable

subroutine calculate_cutoff_weights
    !----------------------------------------------------------------------- 
    !   Output:
    !        sk:  cutoff weights 
    !----------------------------------------------------------------------
    use Constants, only: one
    real(r64) :: efermi,ecutoff,swx
    integer :: k

    efermi = alx ! fermi energy
    ecutoff = ecut ! cut-off energy
    swx = ecutoff/10.d0
    do k = 1, nx1
        sk(k) = 1.0d0/(1.0d0 + dexp((ek(k) - efermi -ecutoff)*(one/swx)))
    enddo
end subroutine calculate_cutoff_weights

subroutine calculate_new_fermi_energy
    !---------------------------------------------------------------
    ! calculate new fermi energy
    !
    !---------------------------------------------------------------
    use MathMethods, only:brak0,rtbrent
    real(r64) :: xl,xh,sl,sh

    ! roughly calculate the range of fermi energy
    call brak0(pnumb3,alx,xl,xh,sl,sh,1.0d0)
    ! more accurately calculate the fermi energy
    alx = rtbrent(pnumb3,xl,xh,sl,sh,1.d-10) 

end subroutine calculate_new_fermi_energy

function pnumb3(alam)
    !----------------------------------------------------------------------------
    !   caculate
    !        (2\sum_{k}v_k^2 - N)
    !
    !   where v_k^2 = 1/2 * (1 - \frac{\epsilon_k - \epsilon_F}{\sqrt{(\epsilon_k - \epsilon_F)^2 + f_k^2\Delta_k^2} })
    !   \epsilon_F is chemical potential(or fermi energy)
    !   \epsilon_k is single particle energy
    !      N       is particle(neutron or proton) number
    !----------------------------------------------------------------------------
    !   INPUT:
    !        alam:  \epsion_F   
    !----------------------------------------------------------------------------
    real(r64),intent(in) :: alam
    real(r64) :: pnumb3
    
    real(r64) :: s,el
    integer :: k

    s   = - tzx
    do k = 1,nx1
        if(isBlocked(k)) then
            s = s + 1.0d0
            cycle
        endif
        el = ek(k) - alam
        s  = s + (1.0d0-el/sqrt(el**2 + (sk(k)*dk(k))**2))
    enddo
    pnumb3 = s
    return
end function pnumb3

subroutine calculate_new_cutoff_energy
    !---------------------------------------------------------------
    ! calculate new cutoff energy
    !
    !---------------------------------------------------------------
    use MathMethods, only:rtbrent
    real(r64) :: ecl,ech,sl,sh

    ! roughly calculate the range of cutoff energy
    ecl = ecut
  1 continue
    if(ecl.gt.0.d0 .and. sfk(ecl).gt.0.d0) then
        ecl = ecl-0.5d0
        goto 1
    endif
    ech = ecut
  2 continue            
    if(sfk(ech).lt.0.d0) then
        ech = ech+0.5d0
        goto 2
    endif
    sl = sfk(ecl)
    sh = sfk(ech)
    ! more accurately calculate the fermi energy
    ecut = rtbrent(sfk,ecl,ech,sl,sh,1.d-10)
end subroutine calculate_new_cutoff_energy

function sfk(pw)
    !----------------------------------------------------------------------------
    !  calcuate 
    !    2*sum_sk - N - 1.65*N^{2/3} 
    !  for smooth pairing cutoff
    !----------------------------------------------------------------------------
    !  INPUT:
    !         pw: cutoff energy
    !---------------------------------------------------------------------------
    use Constants, only: two,third
    real(r64), intent(in) :: pw
    real(r64) :: sfk

    real(r64) :: swx,wd,pw1,s,tmp
    integer :: k

    swx = 0.1d0 * abs(pw)
    wd = 1.0d0 / swx
    pw1 = pw + alx
    s = - tzx - 1.65d0*tzx**(two * third)
    do k =1,nx1
        ! if(isBlocked(k)) then
        !     s = s + tzx + 1.65d0*tzx**(two * third)
        !     s = s - (tzx-1) - 1.65d0*(tzx-1)**(two * third)
        !     ! cycle 
        ! endif
        if( ((ek(k)-pw1)*wd) .le. 500d0 ) then
            s = s + 2.d0/(1.0d0+dexp((ek(k)-pw1)*wd))
        endif
    enddo
    sfk = s
end function sfk

function isBlocked(k)
    !-------------------------------------------------------
    !   If the neutron(it=1) or proton (it=2) is odd number
    !and k energy level is the blocked level, return True
    !------------------------------------------------------
    integer :: k
    logical :: isBlocked
    isBlocked = .False.
    if(MOD(int(tzx),2).eq.1 .and. k.eq.block_level .and. pairing%allow_block) then
        isBlocked = .True.
    endif
end function isBlocked

subroutine set_block_level_of_KPi(ifPrint)
    !--------------------------------------------------------------------
    !  Assign the energy level of the minimum quasiparticle energy
    !  under the given  K^\Pi  to `block_level`.
    !------------------------------------------------------
    use Globals, only: BS
    logical,intent(in),optional :: ifPrint
    real(r64) :: min_qusiparticle_Energy,smax,s,qusiparticle_Energy
    integer :: it,ib,parity,min_qusiparticle_level,nf,i0f,kk,k1,k2,n,imax,kk_parity
    ! find the level
    do it = 1,2
        ib = pairing%block_K(it) ! 1: 1/2, 2: 3/2
        parity = pairing%block_Pi(it) ! 1 for '+'; -1 for '-'
        min_qusiparticle_Energy = 50.0
        min_qusiparticle_level = 0
        if(ib > 0) then 
            nf = BS%HO_cyl%id(ib,1)
            i0f = BS%HO_cyl%ia(ib,1)
            k1 = dirac%ka(ib,it) + 1
            k2 = dirac%ka(ib,it) + dirac%kd(ib,it)
            do kk = k1,k2
                if (dirac%ee(kk,it)-pairing%ala(it).gt.30.0) cycle
                ! search for main oscillator component
                smax = 0.d0
                do n = 1,nf
                    s = abs(dirac%fg(n,kk,it))
                    if (s.gt.smax) then
                        smax = s
                        imax = n
                    endif
                enddo

                ! quasiparticle energy
                qusiparticle_Energy = sqrt((dirac%ee(kk,it)-pairing%ala(it))**2+(pairing%skk(kk,it)*pairing%de(kk,it))**2)
                
                ! parity
                kk_parity = (-1) ** (BS%HO_cyl%nz(i0f+imax) + BS%HO_cyl%ml(i0f+imax))
                ! find min quasiparticle energy with given parity(Pi)
                if(qusiparticle_Energy < min_qusiparticle_Energy  .and. kk_parity == parity ) then
                    min_qusiparticle_Energy = qusiparticle_Energy
                    min_qusiparticle_level = kk
                endif
            enddo   
        endif 
        ! assign block_level
        pairing%block_level(it) = min_qusiparticle_level
    enddo

    if(ifPrint) call printKPi
    contains
    subroutine printKPi
        use Globals,only: outputfile,option
        integer :: K  
        character :: parity_str
        character(len=*), parameter :: format1 = "('   K^\pi:',13x,i2,'/2',a)"
        write(outputfile%u_outputf,*) '*************************BEGIN  printKPi********************'
        if(option%block_method==2) then 
            write(outputfile%u_outputf,*) "Block Method: convergence -> block"
        else if(option%block_method==3) then
            write(outputfile%u_outputf,*) "Block Method: convergence -> block -> convergence"
        endif
        ! Neutron
        K = pairing%block_K(1)
        if(pairing%block_Pi(1)==1) then 
            parity_str = '+'
        else if(pairing%block_Pi(1)==-1) then
            parity_str = '-'
        endif
        if(K>0) then
            write(outputfile%u_outputf,*) "Block Neutron:"
            write(outputfile%u_outputf,format1) 2*K-1, parity_str
            write(outputfile%u_outputf,*) "  block level:", pairing%block_level(1)
        endif
        ! Proton
        K = pairing%block_K(2)
        if(pairing%block_Pi(2)==1) then 
            parity_str = '+'
        else if(pairing%block_Pi(2)==-1) then
            parity_str = '-'
        endif
        if(K>0) then
            write(outputfile%u_outputf,*) "Block Proton:"
            write(outputfile%u_outputf,format1) 2*K-1, parity_str
            write(outputfile%u_outputf,*) "  block level:", pairing%block_level(2)
        endif
        write(outputfile%u_outputf,*) '*************************END printKPi **********************'
    end subroutine printKPi
end subroutine set_block_level_of_KPi

END MODULE BCS