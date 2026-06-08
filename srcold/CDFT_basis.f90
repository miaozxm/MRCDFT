!==============================================================================!
! MODULE Basis                                                                 !
!                                                                              !
! This module contains the variables and routines related to the Harmonic Osc- !
! illator model space.                                                         !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_basis                                                       !
!==============================================================================!
MODULE Basis
use Globals, only: BS
implicit none

contains
! Cylindrical
subroutine set_Cylindrical_HO_basis(ifPrint)
    logical, intent(in), optional :: ifPrint
    call set_Cylindrical_HO_basis_parameters(ifPrint .and. .True.) 
    ! calculate quantum number of basis
    call calculate_Cylindrical_HO_basis_quantum_number_Fermions
    call calculate_Cylindrical_HO_basis_quantum_number_Bosons
    ! call calculate_Spherical_HO_basis_quantum_number(ifPrint .and. .False.)
    ! set gauss points in Cylindrical coordinates
    call set_Cylindrical_gauss_points(ifPrint .and. .True.) ! set_gauss_points should after set_basis_parameters
    ! calculate the basis in Gauss points
    call calculate_Cylindrical_HO_basis_In_Gauss_points(ifPrint .and. .False.)
end subroutine set_Cylindrical_HO_basis

subroutine set_Cylindrical_HO_basis_parameters(ifPrint)
    !------------------------------------------------------------------------------!
    ! subroutine set_basis_parameters                                              !
    !                                                                              !
    ! Set parameters, ...            !
    !------------------------------------------------------------------------------!
    use Constants, only: hbc,two,third,zero
    use Globals, only: input_par,force,nucleus_attributes
    logical,intent(in),optional :: ifPrint
    !----------------set for cylindrical HO basis---------------------
    ! set input parameters
    BS%HO_cyl%n0f = input_par%basis_n0f
    BS%HO_cyl%n0b = input_par%basis_n0b
    BS%HO_cyl%b0  = input_par%basis_b0
    BS%HO_cyl%q   = input_par%basis_q
    BS%HO_cyl%beta0 = input_par%basis_beta0
    ! set other basis parameters
    BS%HO_cyl%hb0 = hbc/(two * force%masses%amu)
    BS%HO_cyl%hom = 41.0* nucleus_attributes%mass_number**(-third)
    if(BS%HO_cyl%b0 <= zero) then
        BS%HO_cyl%b0 = sqrt(two * BS%HO_cyl%hb0/BS%HO_cyl%hom)
    endif
    BS%HO_cyl%bp = BS%HO_cyl%q ** ((-1.d0/6.d0))
    BS%HO_cyl%bz = BS%HO_cyl%q ** ((1.d0/3.d0))
    
    if(ifPrint) call printBasisParameters
    contains
    subroutine printBasisParameters
        use Globals, only: outputfile
        character(len=*),parameter :: format1="(2(a16,i5))", &
                                      format2="(a16,f12.5)"
        write(outputfile%u_config,outputfile%format_title)'set_Cylindrical_HO_basis_parameters'
        write(outputfile%u_config,format1) 'n0f    =',BS%HO_cyl%n0f, 'n0b = ', BS%HO_cyl%n0b
        write(outputfile%u_config,format2) 'hb0    =',BS%HO_cyl%hb0
        write(outputfile%u_config,format2) 'hom    =',BS%HO_cyl%hom
        write(outputfile%u_config,format2) 'b0     =',BS%HO_cyl%b0
        write(outputfile%u_config,format2) 'beta0  =',BS%HO_cyl%beta0
        write(outputfile%u_config,format2) 'b_z/b0(bz) =',BS%HO_cyl%bz
        write(outputfile%u_config,format2) 'b_z    =',BS%HO_cyl%b0*BS%HO_cyl%bz
        write(outputfile%u_config,format2) 'b_prep/b0(bp) =',BS%HO_cyl%bz
        write(outputfile%u_config,format2) 'b_prep =',BS%HO_cyl%b0*BS%HO_cyl%bp
    end subroutine printBasisParameters

end subroutine set_Cylindrical_HO_basis_parameters

subroutine calculate_Cylindrical_HO_basis_quantum_number_Fermions
    !------------------------------------------------------------------------------!
    ! Cylindrical Oscillator-Base for Fermions                                                 !
    !------------------------------------------------------------------------------!
    use Constants, only: nb_max,nhx,nt_max,nz_max,nr_max,ml_max,tp
    integer :: il,k,ib,ip,ilk,im,ir,iz,nn,ipf,ipg,i1,i2,nz1,nr1,ml1,nz2,nr2,ml2,&
               nn1,nza2,nze2,ipa,iba,ie,i3
    integer,dimension(nhx,2) :: nzz,nrr,mll
    integer,dimension(2) :: idd,idds
    logical :: loz,lor,lom 
    BS%HO_cyl%nzm = 0
    BS%HO_cyl%nrm = 0
    BS%HO_cyl%mlm = 0
    il = 0 ! all level numbers
    do k = 1,nb_max !loop over K  
        !-----------------------------------------------------------------------
        ! Basis for the large components f  
        ! 2*n_r + n_z +  m_l  should not be larger than N_F( HO%n0f )                        
        ! k=Omega+1/2, ip=1 for parity '+', ip=2 for parity '-'
        ! Omega = m_l + m_s
        ! ib is an index to order different compinations of Omega and parity                                              
        !-----------------------------------------------------------------------
        ib = k 
        ip = 1
        BS%HO_cyl%ikb(ib) = k 
        BS%HO_cyl%ipb(ib) = ip
        BS%HO_cyl%mb(ib)  = 2
        write(BS%HO_cyl%txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  K = ',k+k-1,'/2',tp(ip)
        ilk = 0
        ! the following three loops make sure that the quantum numbers satisfy nn=2*ir+im+iz 
        ! and pi=(-1)^{n} and k=Omega+1/2=m or m+1
        do im = k-1,k ! loop over quantum number ml
            do ir=0, (BS%HO_cyl%n0f -im)/2 ! loop over quantum number nr
                do iz=0, BS%HO_cyl%n0f  ! loop over quantum number nz
                    nn = 2*ir + im +iz
                    if(nn > BS%HO_cyl%n0f) exit 
                    ilk = ilk + 1
                    if(ilk >nhx) stop '[Basis]: nhx too small'
                    nzz(ilk,ip) = iz
                    nrr(ilk,ip) = ir
                    mll(ilk,ip) = im
                enddo
            enddo
        enddo
        idd(ip) = ilk
        !----------------------------------------------------------------------
        ! Basis for the small components g                            
        ! N_F(for small components) should be larger than N_F(for larger components) by 1, 
        ! the extra states are called spurious states                                            
        !------------------------------------------------
        ipf = ip
        ipg = ip
        ilk = idd(ipg)
        do i1 = 1,idd(ipf)
            nz1 = nzz(i1,ipf)
            nr1 = nrr(i1,ipf)
            ml1 = mll(i1,ipf)

            loz = .true.
            lor = .true.
            lom = .true.

            do i2 = 1, ilk
                nz2 = nzz(i2,ipg)
                nr2 = nrr(i2,ipg)
                ml2 = mll(i2,ipg)
                !nz1+1
                if(loz .and. nz2.eq.nz1+1 .and. nr2.eq.nr1 .and. ml2.eq.ml1) loz = .false.
                ! nr1+1, ml1-1
                if(lor .and. nr2.eq.nr1+1 .and. nz2.eq.nz1 .and. ml2.eq.ml1-1) lor = .false.
                !ml1+1
                if(lom .and. ml2.eq.ml1+1 .and. nr2.eq.nr1 .and. nz2.eq.nz1) lom = .false. 
            enddo
            !add them to the basis of small components
            !the following three 'if' statements check whether the above mentioned relations among the quantum numbers are still satisfied
            if(loz) then
                nn1 = nz1 + 2*nr1 + ml1 + 1
                nza2= mod(nn1-ml1,2)/2                                
                nze2= (nn1-ml1)/2 
                if (nza2.le.(nz1+1)/2.and.(nz1+1)/2.le.nze2) then ! always true
                    ilk = ilk + 1                                      
                    nzz(ilk,ipg) = nz1 + 1                             
                    nrr(ilk,ipg) = nr1                                 
                    mll(ilk,ipg) = ml1
                endif
            endif
            if (lor) then                                            
                 nn1= nz1 + 2*nr1 + ml1 +1                             
                 if (k-1.le.ml1-1.and.ml1-1.le.min(k,nn1)) then  ! ml1==k   
                     ilk = ilk + 1                                      
                     nzz(ilk,ipg) = nz1                                 
                     nrr(ilk,ipg) = nr1 + 1                             
                     mll(ilk,ipg) = ml1 - 1                             
                 endif                                                 
            endif
            if (lom) then                                            
                nn1= nz1 + 2*nr1 + ml1 +1                             
                if (k-1.le.ml1+1.and.ml1+1.le.min(k,nn1)) then  !ml1==k−1      
                    ilk = ilk + 1                                      
                    nzz(ilk,ipg) = nz1                                 
                    nrr(ilk,ipg) = nr1                                 
                    mll(ilk,ipg) = ml1 + 1    
                endif                                                 
            endif
            
            if(ilk > nhx) stop '[Basis]: nhx too small'
        enddo
        idds(ipg)= ilk - idd(ipg)
        !--------------------------------------------------------------------------------------------
        ! Reordering and construction of the fields ia, id, iag0
        ! Attention:  note that through ipa=3-ip; ib=2*(k-1)+ip, iba=2*(k-1)+ipa ???
        !             and id(ib,2)=idd(ipa)+idds(ipa)
        !             and idg0(ib)=idds(ipa)
        !             and ia(iba,2)=il
        !             and iag0(iba)=il+idd(ip)
        !             the index has been redirected, so the ordering described at the head is ahieved.
        !          id(ib,1) ia(ib,1) for large component f 
        !          id(ib,2) ia(ib,2) for small component g
        !--------------------------------------------------------------------------------------------
        ipa = ip
        iba = ib
        BS%HO_cyl%id(ib,1) = idd(ip)
        BS%HO_cyl%id(ib,2) = idd(ipa) + idds(ipa)
        BS%HO_cyl%idg0(ib) = idds(ipa)

        BS%HO_cyl%ia(ib,1) = il
        BS%HO_cyl%ia(iba,2) = il
        BS%HO_cyl%iag0(iba) = il + idd(ip)

        ie = idd(ip) + idds(ip)
        if(il+ie > nt_max) stop '[Basis]: nt_max too small'
        do i3=1,ie
            il = il + 1
            BS%HO_cyl%nz(il) = nzz(i3,ip)
            BS%HO_cyl%nr(il) = nrr(i3,ip)
            BS%HO_cyl%ml(il) = mll(i3,ip)                                       
            BS%HO_cyl%ms(il) = (-1)**(k-BS%HO_cyl%ml(il)+1) !ms=1: spin up Omega=m+1/2, ms=-1: spin down, Omega=m-1/2
            nn = BS%HO_cyl%nz(il) +2*BS%HO_cyl%nr(il) + BS%HO_cyl%ml(il)
            if(nn < 10) then
                write(BS%HO_cyl%tb(il),"(i2,a1,'[',3i1,']')") 2*k-1,tp(ip),nn,BS%HO_cyl%nz(il),BS%HO_cyl%ml(il)
            else
                write(BS%HO_cyl%tb(il),"(4i2)") 2*k-1,nn,BS%HO_cyl%nz(il),BS%HO_cyl%ml(il)
            endif
            BS%HO_cyl%nzm = max0(BS%HO_cyl%nzm,BS%HO_cyl%nz(il))
            BS%HO_cyl%nrm = max0(BS%HO_cyl%nrm,BS%HO_cyl%nr(il))
            BS%HO_cyl%mlm = max0(BS%HO_cyl%mlm,BS%HO_cyl%ml(il))
        enddo
    enddo 
    BS%HO_cyl%nt = il
    BS%HO_cyl%nb = ib
    if(BS%HO_cyl%mlm > ml_max) stop '[Basis]: ml_max too small'
    if(BS%HO_cyl%nrm > nr_max) stop '[Basis]: nr_max too small'
    if(BS%HO_cyl%nzm > nz_max) stop '[Basis]: nz_max too small'

end subroutine calculate_Cylindrical_HO_basis_quantum_number_Fermions

subroutine calculate_Cylindrical_HO_basis_quantum_number_Bosons
    !------------------------------------------------------------------------------!
    ! Cylindrical Oscillator-Base for Bosons                                                   !
    !------------------------------------------------------------------------------!
    use Globals, only: ntb_max,nr_max_boson,nz_max_boson
    integer :: il,iz,ir
    ! if (mod(BS%HO_cyl%n0b,2).ne.0) stop ' [Basis]: n0b must be even'
    BS%HO_cyl%nzbm = 0
    BS%HO_cyl%nrbm = 0
    il = 0
    do iz =0,BS%HO_cyl%n0b !loop over nz-quantum number
        do ir=0, (BS%HO_cyl%n0b-iz)/2 !loop over nr-quantum number
            if (iz + 2*ir > BS%HO_cyl%n0b) exit
            il = il + 1
            if(il > ntb_max) stop '[Basis]: ntb_max too small'
            BS%HO_cyl%nrb = ir
            BS%HO_cyl%nzb = iz
            BS%HO_cyl%nrbm = max0(BS%HO_cyl%nrbm,ir)
            BS%HO_cyl%nzbm = max0(BS%HO_cyl%nzbm,iz)
        enddo
    enddo
    BS%HO_cyl%ntb = il
    if(BS%HO_cyl%nrbm > nr_max_boson) stop '[Basis]: nr_max_boson too small'
    if(BS%HO_cyl%nzbm > nz_max_boson) stop '[Basis]: nz_max_boson too small'
end subroutine calculate_Cylindrical_HO_basis_quantum_number_Bosons

subroutine set_Cylindrical_gauss_points(ifPrint)
    use Globals, only: gauss
    use Constants, only: ngh,ngl,pi
    use MathMethods, only: GaussLaguerre, GaussHermite
    logical, intent(in), optional :: ifPrint
    integer :: ih,il
    gauss%nh = ngh
    gauss%nl = ngl
    ! set gauss integral parameters
    call GaussHermite(gauss%nh,gauss%xh,gauss%w_hermite)
    call GaussLaguerre(gauss%nl,gauss%xl,gauss%w_laguerre)
    !set z-axis for both - and +
    do ih=1,gauss%nh
        gauss%zb(ih) = gauss%xh(ih)*BS%HO_cyl%b0*BS%HO_cyl%bz
        gauss%wh(ih) = gauss%w_hermite(ih)*exp(gauss%xh(ih)**2)
    enddo
    !set r-axis
    do il=1,gauss%nl
        gauss%rb(il) = sqrt(gauss%xl(il))*BS%HO_cyl%b0*BS%HO_cyl%bp
        gauss%wl(il) = gauss%w_laguerre(il)*exp(gauss%xl(il))
        do ih=1,gauss%nh
            gauss%wdcor(ih,il) = BS%HO_cyl%b0**3 * BS%HO_cyl%bz*BS%HO_cyl%bp**2 * pi * gauss%wh(ih) * gauss%wl(il)
        enddo
    enddo

    if(ifPrint) call printGaussPoints
    contains
    subroutine printGaussPoints
        use Globals, only: outputfile
        character(len=*), parameter :: format1= "(a30,15f9.5,/,30x,15f9.5,/,30x,5f9.5)"
        write(outputfile%u_config,outputfile%format_title) 'set_Cylindrical_gauss_points'
        write(outputfile%u_config,*) 'Gauss-Hermite:  b_z = ',BS%HO_cyl%b0*BS%HO_cyl%bz
        write(outputfile%u_config,format1) 'points(zeta)               :',gauss%xh
        write(outputfile%u_config,format1) 'weight(w_h)                :',gauss%w_hermite
        write(outputfile%u_config,format1) 'zb(z=zeta*b_z)             :',gauss%zb
        write(outputfile%u_config,format1) 'wh(w_h*e^{zeta^2})         :',gauss%wh
        write(outputfile%u_config,*) 'Gauss-Lagueere: b_prep = ',BS%HO_cyl%b0*BS%HO_cyl%bp
        write(outputfile%u_config,format1) 'points(eta)                :',gauss%xl
        write(outputfile%u_config,format1) 'weight(w_l)                :',gauss%w_laguerre
        write(outputfile%u_config,format1) 'rb(r_prep=sqrt(eta)*b_prep):',gauss%rb
        write(outputfile%u_config,format1) 'wl(w_l*e^{eta})            :',gauss%wh    
    end subroutine printGaussPoints
end subroutine set_Cylindrical_gauss_points

subroutine calculate_Cylindrical_HO_basis_In_Gauss_points(ifPrint)
    !--------------------------------------------------------------------------------!
    ! Calculates the wavefunctions for the cylindrical oscillator.      
    !     they are given as:                                                 
    !                                                                       
    !     Phi(zz,rr,phi) = 1/((b0**3 bz*bp*bp)^(1/2)) * psi_nz(zz) *        
    !                      psi_(L,nr)(rr) * exp(i*L*phi)/sqrt(2*pi)         
    !                                                                       
    !     zz is z-coordinate in units of fm                                 
    !     rr is perpendicular coordinate in units of fm                    
    !                                                                       
    !     psi_nz(zz)     = N_nz * H_nz(z) * exp(-z*z/2)                                                                                            
    !     psi_(L,nr)(rr) = N_(nr,L) * sqrt(2) * eta^(L/2) * L_nr^L(eta) * exp(-eta/2)   
    !                                                                       
    !     z = zz/(bz*b0),    r = rr/(bp*b0),       eta = r*r                                                                                      
    !     N_nz     = 1/sqrt(sqrt(pi) * 2^nz * nz!)                                                                                                 
    !     N_(nr,L) = sqrt( nr! / (nr+L)! )                                  
    !                                                                       
    !                                                                       
    !     the contribution to the density from the level i is:                                                                                     
    !     rho_k(zz,rr)= 1/(2*pi b0*3 bz*bp*bp) * ( psi_nz(zz) * psi_(L,nr)(rr) )^2                 
    !                                                                       
    !---- The z-function at meshpoint xh(ih) is stored in QH(nz,ih)         
    !     such that QH is normalized in such way that the                   
    !     norm integral reads                                               
    !                                                                       
    !     \int dz_(-inf)^(inf) (psi_nz)**2 = 1 = \sum_i QH(nz,i)**2         
    !                                                                       
    !     this means, that QH contains the following factors:               
    !                                                                       
    !     a)  the z-part of the wavefunction psi_nz(zz)                     
    !     b)  the Gaussian weight sqrt( WH(i) ):                            
    !         \inf_(-inf)^inf f(z) dz = \sum_i f(x_i) * WH(i)               
    !                                                                       
    !     having QH(nz,i) we get the z-wavefunction:                                                                                               
    !     psi_nz(zz) =  QH(nz,i) / sqrt( WH(i) )                            
    !                                                                       
    !---- The r-function at meshpoint XL(il) is stored in QL(nr,L,il)       
    !     such that QL is normalized in such way that the                   
    !     2-dimensional norm integral reads                                 
    !                                                                       
    !     \int_0^inf r dr (phi_(L,nr)**2 = 1 = \sum_i QL(nr,L,i)**2         
    !                                                                       
    !     this means, that QL contains the following factors:              
    !                                                                       
    !     a)  the part of the wavefunction psi_(nr,L)(rr)                   
    !     b)  a factor sqrt(1/2) from the transformation from r to eta      
    !     c)  the Gaussian weight sqrt( WL(i) ):                            
    !         \inf_0^inf f(eta) d(eta) = \sum_i f(XL(i)) * WL(i)            
    !                                                                       
    !     having QL(nr,L,i) we get the r-wavefunction:                                                                                             
    !     psi_(nr,L)(rr) =  QL(nr,L,i) * sqrt( 2 / WL(i) )                  
    !                                                                       
    !                                                                       
    !---- the density contribution from the level k                         
    !                                                                       
    !     rho_k(zz,rr)= (QH(nz,ih) * QL(nr,L,il))^2 /                       
    !                         (  pi * WH(ih)*WL(il)* b0**3 * bz*bp*pb)      
    !                                                                       
    !----------------------------------------------------------------------
    !                                                                       
    !     QH1 contains the z-derivatives in the following form:             
    !                                                                       
    !     d/dz psi_nz(zz) = QH1(nz,i) / sqrt( WH(i) )                       
    !                                                                      
    !     QL1 contains the r-derivatives in the following form:             
    !                                                                       
    !     r d/dr psi_(nr,L)(rr) = QL1(nr,L,i) * sqrt( 2 / WL(i) )           
    !     note this difference with the text                                
    !----------------------------------------------------------------------
    !                                                                       
    !     QHB(nz,i) is the z-function for the expansion of the mesonfields  
    !     QLB(nr,i) is the r-function for the expansion of the mesonfields  
    !                                                                       
    !     QHB(nz,i) = psi_nz(zz) / (b0b*bz)^(1/2)                          
    !     QLB(nr,i) = psi_(nr,L=0)(rr) / ( sqrt(2*pi) * b0b*bp)             
    !                                                                       
    !----------------------------------------------------------------------
    ! 
    !Output:
    ! QH (0:nz_max,  1:ngh)  : z-function at meshpoint xh(ih) with sqrt(wh)
    ! QH1(0:nz_max,  1:ngh)  : the first-order derivative of z-funtion with sqrt(wh)
    ! QL (0:nr_max, 0:ml_max, 1:ngl)  : r-function at meshpoint xl(il) with sqrt(wl/2) 
    ! QL1(0:nr_max, 0:ml_max, 1:ngl)  : the first-order derivative of r-funtion with r*sqrt(wl/2);NOTE: r d/dr psi_(nr,L)(rr) = QL1(nr,L,i) * sqrt( 2 / WL(i) )
    ! QHB(0:nz_max_boson, 1:ngh)  : the z-function for the expansion of the mesonfields
    ! QLB(0:nr_max_boson, 1:ngl)  : the r-function for the expansion of the mesonfields
    !----------------------------------------------------------------------
    use Constants, only:  r64,nz_max,nr_max,ml_max,nr_max_boson, nz_max_boson,pi,one,half
    use Globals, only: gauss,gfv
    logical, intent(in), optional :: ifPrint
    real(r64) :: w4pii,z,psiH_0,qh_0,cb,zb,psiHB_0, &
                 x,psiL_00,ql_l0,tmp0,tmp1,tmp2,tmp3,tmp4, xb,psiLB_00
    integer :: ih,n,il,l
    allocate(BS%HO_cyl%qh (0:nz_max, 1:gauss%nh))
    allocate(BS%HO_cyl%qh1(0:nz_max, 1:gauss%nh))
    allocate(BS%HO_cyl%ql (0:nr_max, 0:ml_max, 1:gauss%nl))
    allocate(BS%HO_cyl%ql1(0:nr_max, 0:ml_max, 1:gauss%nl))
    allocate(BS%HO_cyl%qhb (0:nz_max_boson, 1:gauss%nh))
    allocate(BS%HO_cyl%qlb (0:nr_max_boson, 1:gauss%nl))

    !----------------------------------------------------------------------
    ! z-dependence
    !----------------------------------------------------------------------
    w4pii = pi**(-0.25d0)
    cb = one
    do ih = 1, gauss%nh
        !------------basis for the fermions----------------
        z = gauss%xh(ih)
        psiH_0 = w4pii*exp(-half*z**2)
        qh_0 = psiH_0*sqrt(gauss%wh(ih))
        BS%HO_cyl%qh(0,ih) = qh_0
        BS%HO_cyl%qh(1,ih) = gfv%sq(2) * qh_0 * z
        BS%HO_cyl%qh1(0,ih) = - qh_0 * z
        BS%HO_cyl%qh1(1,ih) = gfv%sq(2) * qh_0 * (one -z**2)
        do n = 2, BS%HO_cyl%nzm
            BS%HO_cyl%qh(n,ih) = gfv%sqi(n)*(gfv%sq(2)*z*BS%HO_cyl%qh(n-1,ih)-gfv%sq(n-1)*BS%HO_cyl%qh(n-2,ih))
            BS%HO_cyl%qh1(n,ih) = gfv%sq(n+n) * BS%HO_cyl%qh(n-1,ih) - z * BS%HO_cyl%qh(n,ih)
        enddo
        !------------basis for the bosons---------------- 
        zb = z*cb
        psiHB_0 = w4pii*exp(-half*zb**2)/sqrt(BS%HO_cyl%b0*BS%HO_cyl%bz)
        BS%HO_cyl%qhb(0,ih) = psiHB_0
        BS%HO_cyl%qhb(1,ih) = gfv%sq(2) * psiHB_0 * zb
        do n = 2, nz_max_boson
            BS%HO_cyl%qhb(n,ih) = gfv%sqi(n)*(gfv%sq(2)*zb*BS%HO_cyl%qhb(n-1,ih)-gfv%sq(n-1)*BS%HO_cyl%qhb(n-2,ih))
        enddo
    enddo
    !----------------------------------------------------------------------
    ! r-dependence
    ! note here x is actually eta=r*r
    !----------------------------------------------------------------------
    do il = 1, gauss%nl
        !------------basis for the fermions----------------
        x = gauss%xl(il)
        psiL_00 = gfv%sq(2) * exp(-half*x)
        do l = 0,BS%HO_cyl%mlm
              ql_l0 = psiL_00*sqrt(half*gauss%wl(il)*x**l)
              BS%HO_cyl%ql(0,l,il) = ql_l0 * gfv%wfi(l)
              BS%HO_cyl%ql(1,l,il) = ql_l0 * (l+1-x) * gfv%wfi(l+1)
              BS%HO_cyl%ql1(0,l,il) = ql_l0 * (l-x) * gfv%wfi(l) ! why not ql_10 * (l-x)/(sqrt(x)) ?
              BS%HO_cyl%ql1(1,l,il) = ql_l0 * (l*l+l-x*(l+l+3)+x*x) * gfv%wfi(l+1)
            do n = 2,BS%HO_cyl%nrm
                tmp0 = gfv%sqi(n) * gfv%sqi(n+l)
                tmp1 = n+n+l-1-x
                tmp2 = gfv%sq(n-1) * gfv%sq(n-1+l)
                tmp3 = n+n+l-x
                tmp4 = 2 * gfv%sq(n) * gfv%sq(n+l)
                BS%HO_cyl%ql(n,l,il) = (tmp1*BS%HO_cyl%ql(n-1,l,il) - tmp2*BS%HO_cyl%ql(n-2,l,il))*tmp0
                BS%HO_cyl%ql1(n,l,il) =  tmp3*BS%HO_cyl%ql(n,l,il)  -tmp4*BS%HO_cyl%ql(n-1,l,il)
            enddo
        enddo
        !------------basis for the bosons----------------
        xb = gauss%xl(il) * cb **2 ! bosons
        psiLB_00 = gfv%sq(2)*exp(-half*xb) / (sqrt(2*pi)*BS%HO_cyl%b0*BS%HO_cyl%bp) 
        BS%HO_cyl%qlb(0,il) = psiLB_00
        BS%HO_cyl%qlb(1,il) = psiLB_00 * (1-xb)
        do n = 2,nr_max_boson
            BS%HO_cyl%qlb(n,il) = ((2*n-1-xb) * BS%HO_cyl%qlb(n-1,il) - (n-1) * BS%HO_cyl%qlb(n-2,il)) / n
        enddo
    enddo

    if(ifPrint .and. .false.) call printCylindricalHOBasis
    if(ifPrint) call test_integration_of_Cylindrical_HO_basis
    contains
    subroutine printCylindricalHOBasis
        use Globals, only: outputfile
        character(len=*),parameter :: format1 = "(a,i3,30f15.8)", &
                                      format2 = "(a,i4,i3,30f15.8)"
        integer :: n,l
        integer :: ix = 10
        write(outputfile%u_config,outputfile%format_title) 'calculate_Cylindrical_HO_basis_In_Gauss_points'
        do n = 0,min(BS%HO_cyl%nzm,BS%HO_cyl%nzbm)
            write(outputfile%u_config,'(a,i2,a)') ' QH(nz=',n,',ih=1...)'
            write(outputfile%u_config,format1) ' H     xh',n,(BS%HO_cyl%qh(n,ih),ih=1,ix)              
            write(outputfile%u_config,format1) ' dH/dx xh',n,(BS%HO_cyl%qh1(n,ih),ih=1,ix)             
            write(outputfile%u_config,format1) ' H  cb*xh',n,(BS%HO_cyl%qhb(n,ih)*sqrt(gauss%wh(ih)),ih=1,ix)
        enddo
        do l=0, BS%HO_cyl%mlm
            write(outputfile%u_config,*) '      nr ml    QL(nr,l,il=1,...)'
            do n=0,BS%HO_cyl%nrm
                write(outputfile%u_config,format2)'ql  ', n,l,(BS%HO_cyl%ql(n,l,il),il=1,ix) 
                write(outputfile%u_config,format2)'ql1 ', n,l,(BS%HO_cyl%ql1(n,l,il),il=1,ix)
                if(l==0) then
                write(outputfile%u_config,format2)'qlb ', n,l,(BS%HO_cyl%qlb(n,il),il=1,ix) 
                endif
            enddo
        enddo
    end subroutine printCylindricalHOBasis
end subroutine calculate_Cylindrical_HO_basis_In_Gauss_points

subroutine test_integration_of_Cylindrical_HO_basis
    use Constants, only: r64,nz_max,nz_max_boson,ml_max,nr_max,pi
    use Globals, only: gauss,outputfile
    real(r64):: s,sb,s2
    integer :: n1,n2,ih,l,il
    character(len=*),parameter :: format1 = "(' G.Hermit: n1 =',i3,' n2 =',i3, 3f15.8)",&
                                  format2 = "(' G.Lague.: l =' ,i2,' n1 =',i3,' n2 =',i3, 3f15.8)"
    write(outputfile%u_config,outputfile%format_title) 'test_integration_of_Cylindrical_HO_basis'
    ! test for hermite integration
    do n1 = 0,min(BS%HO_cyl%nzm,BS%HO_cyl%nzbm)
         do n2 = 0,n1
               s  = 0.0d0
               sb = 0.0d0
               s2 = 0.0d0
               do ih = 1,gauss%nh
                  s = s  + BS%HO_cyl%qh(n1,ih)*BS%HO_cyl%qh(n2,ih)
                  if (n1.le.nz_max_boson) sb = sb + BS%HO_cyl%qhb(n1,ih)*BS%HO_cyl%qhb(n2,ih)*gauss%wh(ih)
               enddo
               write(outputfile%u_config, format1) n1,n2,s,sb*BS%HO_cyl%b0*BS%HO_cyl%bz
         enddo
    enddo
    ! test for Laguerre integration
    do l = 0,BS%HO_cyl%mlm
        do n1 = 0,min(BS%HO_cyl%nrm,BS%HO_cyl%nrbm)
            do n2 = 0,n1
                s  = 0.0d0
                sb = 0.0d0
                s2 = 0.0d0
                do il = 1,gauss%nl
                    s = s + BS%HO_cyl%ql(n1,l,il)*BS%HO_cyl%ql(n2,l,il)
                    if (l.eq.0) sb = sb + BS%HO_cyl%qlb(n1,il)*BS%HO_cyl%qlb(n2,il)*gauss%wl(il)
                enddo
                write(outputfile%u_config,format2) l,n1,n2,s,sb*pi*(BS%HO_cyl%b0*BS%HO_cyl%bp)**2
            enddo                                                       
        enddo                                                       
    enddo 
end subroutine test_integration_of_Cylindrical_HO_basis

! Cylindrical -> Spherical
subroutine transform_coefficients_form_cylindrical_to_spherical(ifPrint)
    !----------------------------------------------------------
    !   calculate the expansion coefficients  of s.p. states on 
    ! spherical oscillator basis
    !----------------------------------------------------------
    use Constants, only: r64,nh2x
    use Globals, only: dirac
    logical,intent(in) :: ifPrint
    integer :: it,ib,k1,k2,ibsp,k,ifg,ndsp,i0sp,nfgsp,nd,i0,nfg,&
               il,nr1,nl1,nj1,nm1,nnsp,n,nz2,nr2,ml2,ms2,nncyl
    real(r64) :: sm,coef
    do it = 1,2 ! loop over neutrons and protons
        do ib = 1,BS%HO_cyl%nb ! loop over K and parity-blocks
            k1 = dirac%ka(ib,it) + 1
            k2 = dirac%ka(ib,it) + dirac%kd(ib,it)
            ibsp=1
            do k = k1,k2 ! loop over S.P.States
                do ifg = 1,2 ! loop over contributions from large and small components
                    ! sph.
                    ndsp  = BS%HO_sph%idsp(ibsp,ifg)
                    i0sp  = BS%HO_sph%iasp(ibsp,ifg)
                    nfgsp = (ifg-1)*BS%HO_sph%idsp(ibsp,1)
                    ! cyl.
                    nd  = BS%HO_cyl%id(ib,ifg)
                    i0  = BS%HO_cyl%ia(ib,ifg)
                    nfg = (ifg-1)*BS%HO_cyl%id(ib,1)
                    do il= 1,ndsp ! loop of spher. oscillator basis states
                        nr1= BS%HO_sph%nljm(i0sp+il,1)
                        nl1= BS%HO_sph%nljm(i0sp+il,2)
                        nj1= BS%HO_sph%nljm(i0sp+il,3)
                        nm1= BS%HO_sph%nljm(i0sp+il,4) 
                        nnsp = 2*(nr1-1) + nl1 !N = 2(N_r - 1) + l
                        sm = 0.d0
                        ! if (abs(nm1).ne.abs(ikb(ib))) goto 112
                        do n = 1,nd ! loop of cylin. oscillator basis states n
                            nz2 = BS%HO_cyl%nz(n+i0)
                            nr2 = BS%HO_cyl%nr(n+i0)
                            ml2 = BS%HO_cyl%ml(n+i0)
                            ms2= BS%HO_cyl%ms(n+i0)  ! ms=-1(dowm) and ms=+1 (up)
                            if(ms2.eq.-1) ms2=0 ! set ms= 0(dowm) and ms=+1 (up)
                            nncyl = nz2 + 2*nr2 + ml2 ! N = 2n_p + m_l + n_z
                            ! nss2 + ml2 = abs(kb(ib)) = Omega > 0
                            if (nncyl.ne.nnsp) cycle
                            coef = transsph(nz2,nr2,ml2,ms2,nr1,nl1,nj1,nm1)
                            sm   = sm + dirac%fg(nfg+n, k, it)*coef
                        enddo !n
                        ! 112 continue
                        BS%HO_sph%fg(nfgsp+il,k,it) = sm
                        if(nfgsp+il>nh2x) stop 'nh2x too small !'
                    enddo ! il 
                enddo ! ifg
            enddo !k
        enddo ! ib 
    enddo ! it

    if(ifPrint) call test_integration_of_Spherical_and_Cylindrical_HO_basis
end subroutine transform_coefficients_form_cylindrical_to_spherical

subroutine test_integration_of_Spherical_and_Cylindrical_HO_basis
    use Constants, only: r64,nkx,nbx
    use Globals, only: outputfile,dirac,BS
    integer :: it,ib,k1,k2,k,ifg,ndsp,i0sp,nfgsp,nd,i0,nfg,il,n
    real(r64),dimension(nkx,nbx,2) :: ormfg,ormsphfg
    character(len=*),parameter :: format1 = "(a,i3,a,i2,a,i1,a,2f12.8)"
    write(outputfile%u_outputf,*) '************* BEGIN test_integration_of_Spherical_and_Cylindrical_HO_basis **************'
    write(outputfile%u_outputf,*) '(  k,ib,it) fg^2_{Cart}  fg^2_{sph.}'
    ! loop over neutrons and protons
    do it = 1,2
        ! loop over simplex-parity-blocks
        do ib = 1,BS%HO_cyl%nb
            k1 = dirac%ka(ib,it) + 1
            k2 = dirac%ka(ib,it) + dirac%kd(ib,it)
            ! loop over S.P.States
            do k = k1,k2
                ormfg(k,ib,it)=0
                ormsphfg(k,ib,it)=0
                !loop over contributions from large and small components
                do ifg = 1,2
                    ! the block of spher. HO basis states
                    ndsp  = BS%HO_sph%idsp(1,ifg)
                    i0sp  = BS%HO_sph%iasp(1,ifg)
                    nfgsp = (ifg-1)*BS%HO_sph%idsp(1,1)
                    ! the block of Cartesian HO basis states
                    nd  = BS%HO_cyl%id(ib,ifg)
                    i0  = BS%HO_cyl%ia(ib,ifg)
                    nfg = (ifg-1)*BS%HO_cyl%id(ib,1)

                    !loop of spher. oscillator basis states------------
                    do il=1,ndsp
                        ! nr= BS%HO_sph%nljm(i0sp+il,1)
                        ! nl= BS%HO_sph%nljm(i0sp+il,2)
                        ! nj= BS%HO_sph%nljm(i0sp+il,3)
                        ! nm= BS%HO_sph%nljm(i0sp+il,4)
                        ormsphfg(k,ib,it)=ormsphfg(k,ib,it)+BS%HO_sph%fg(nfgsp+il,k,it)*BS%HO_sph%fg(nfgsp+il,k,it)
                    enddo !il

                    ! loop of Cartesian oscillator basis states n
                    do n =  1,nd
                        ! nz2 = BS%HO_cyl%nz(n+i0)
                        ! nr2 = BS%HO_cyl%nr(n+i0)
                        ! ml2 = BS%HO_cyl%ml(n+i0)
                        ormfg(k,ib,it)=ormfg(k,ib,it)+dirac%fg(nfg+n,k,it)*dirac%fg(nfg+n,k,it)
                    enddo !n
                enddo !ifg
                write(outputfile%u_outputf,format1) '(', k,', ', ib,', ',it,')=',ormfg(k,ib,it), ormsphfg(k,ib,it)
            enddo !k
        enddo !ib
    enddo !it

    write(outputfile%u_outputf,*) 'fg(k,ib,it)*fg(k+1,ib,it)    fgsph(k,ib,it)*fgsph(k+1,ib,it) '
    ! loop over neutrons and protons
    do it = 1,2
        ! loop over simplex-parity-blocks
        do ib = 1,BS%HO_cyl%nb
            k1 = dirac%ka(ib,it) + 1
            k2 = dirac%ka(ib,it) + dirac%kd(ib,it)
            ! loop over S.P.States
            do k = k1,k2-1
                ormfg(k,ib,it)=0
                ormsphfg(k,ib,it)=0
                ! loop over contributions from large and small components
                do ifg = 1,2
                    ! the block of spher. HO basis states
                    ndsp  = BS%HO_sph%idsp(1,ifg)
                    i0sp  = BS%HO_sph%iasp(1,ifg)
                    nfgsp = (ifg-1)*BS%HO_sph%idsp(1,1)
                    ! the block of Cartesian HO basis states
                    nd  = BS%HO_cyl%id(ib,ifg)
                    i0  = BS%HO_cyl%ia(ib,ifg)
                    nfg = (ifg-1)*BS%HO_cyl%id(ib,1)
                    ! loop of spher. oscillator basis states------------
                    do il=1,ndsp
                        ! nr= BS%HO_sph%nljm(i0sp+il,1)
                        ! nl= BS%HO_sph%nljm(i0sp+il,2)
                        ! nj= BS%HO_sph%nljm(i0sp+il,3)
                        ! nm= BS%HO_sph%nljm(i0sp+il,4)
                        ormsphfg(k,ib,it)=ormsphfg(k,ib,it)+BS%HO_sph%fg(nfgsp+il,k,it)*BS%HO_sph%fg(nfgsp+il,k+1,it)
                    enddo !il

                    ! loop of Cartesian oscillator basis states n
                    do n =  1,nd 
                        ! nz2 = BS%HO_cyl%nz(n+i0)
                        ! nr2 = BS%HO_cyl%nr(n+i0)
                        ! ml2 = BS%HO_cyl%ml(n+i0)
                        ormfg(k,ib,it)=ormfg(k,ib,it)+dirac%fg(nfg+n,k,it)*dirac%fg(nfg+n,k+1,it)
                    enddo !n
                enddo !ifg
                write(outputfile%u_outputf,format1) '(', k,', ', ib,', ',it,')=',ormfg(k,ib,it), ormsphfg(k,ib,it)
            enddo !k
        enddo !ib
    enddo !it

    write(outputfile%u_outputf,*) '************END test_integration_of_Spherical_and_Cylindrical_HO_basis ******************'
end subroutine test_integration_of_Spherical_and_Cylindrical_HO_basis

function transsph(nz,nr,ml,ms,n,l,j,m) 
    !----------------------------------------------------------------
    !  m=ml+ms, ms=0 (down);1(up)
    !  Transformation from cylindrical oscillator with spin | ms >
    !  to spherical oscillator |n s l j m> with spin | ms >
    !
    !----------------------------------------------------------------
    implicit real*8(a-h,o-z)
    implicit integer(i-n)
    double precision transsph
    transsph = 0.0d0
    if(m.ne.ml+ms) return
    transsph = clebslj(l,j,ms,m)*cylsph(nz,nr,ml,n,l,ml)
end function transsph

function clebslj(l,j,s,m)
    !---------------------------------------------------------------------
    !     Calculates  Clebsch-Gordan coefficients <l 1/2; m_l m_s| j m_j>     
    !     l = L, m_l = M-S, m_s = S-1/2, j = J-1/2, m_j = M -1/2
    !     Table 1: https://scipp.ucsc.edu/~haber/ph215/clebsch18.pdf
    !---------------------------------------------------------------------
    use Globals, only: gfv
    double precision clebslj
    integer :: L,J,S,M
    clebslj = 0.0d0
    ! spin parallel orbital angular momentum
    if (J.eq.L+1) then
        if (S.eq.1) clebslj = gfv%sq(L+M)*gfv%sqi(2*L+1)
        if (S.eq.0) clebslj = gfv%sq(L-M+1)*gfv%sqi(2*L+1)
    endif
    ! spin anti-parallel orbital angular momentum
    if (J.eq.L) then
        if (S.eq.1) clebslj = -gfv%sq(L-M+1)*gfv%sqi(2*L+1)
        if (S.eq.0) clebslj =  gfv%sq(L+M)*gfv%sqi(2*L+1)
    endif
end function clebslj

function cylsph(nz,nr,ml,n,l,m) 
    !------------------------------------------------------------------
    !     Transformation from cylindrical oscillator to spherical oscillator
    !     yields <nz, nr, ml | n, l, ml >  
    !------------------------------------------------------
    !     Ref.: R.R.Chasman and S. Wahlborn, NPA90(1967)401
    !
    !     | nz,nr,ml > = |nz> | nr,ml > exp( i ml phi )/sqrt(2 pi)
    !       where | nr,ml > = sqrt( 2 nr! / (nr+|ml|)!) r_p ** |ml|
    !                         *  L^|ml|_nr ( r_p**2 ) exp( - r_p**2/2)
    !     np(olar) = 2*nr + |ml|
    !
    !     | n l m >  = sqrt( 2 (n-1)! / ( n+l-1/2)! ) r**l L^(l+1/2)_(n-1) (r**2)
    !                  exp( -r**2/2 )  Y_lm (theta,phi)
    !-----------------------------------------------------------------------
    implicit real*8(a-h,o-z)
    implicit integer(i-n)
    double precision cylsph
    cylsph = 0.d0
    if(ml.ne.m) return
    nn = 2*(n-1) + l
    do nx=0,nn-nz
        ny=nn-nz-nx
        c1 = carsph(nx,ny,nz,n,l,m)
        c2 = carzyl(nx,ny,nz,nz,nr,ml)
        if(c2.ne.0.d0) cylsph = cylsph + c1*c2
    enddo
end function cylsph

function carsph(nx,ny,nz,n,l,m) 

    !-----------------------------------------------------
    !     Transformation from cartesian oscillator to spherical oscillator
    !     yields <nx, ny, nzl | n, l, m >  
    !
    !     | nx,ny,nz > = | nx > i**ny | ny > | nz >
    !       where | n > = pi**(-1/4) (2**n n!)**(-1/2) H_N(x) exp(-x**2/2)
    !
    !     | n l m >  = sqrt( 2 (n-1)! / ( n+l-1/2)! ) r**l L^(l+1/2)_(n-1) (r**2)
    !                  exp( -r**2/2 )  Y_lm (theta,phi)
    !-----------------------------------------------------------------------
    use Constants, only: igfvbc
    use Globals, only: gfv
    implicit real*8(a-h,o-z)
    implicit integer(i-n)
    double precision carsph
    carsph = 0.0d0
    ma = iabs(m)
    np = nx + ny
    nr = np - ma
    if (nr.lt.0.or.mod(nr,2).ne.0) return
    nr = nr/2
    nn = 2*(n-1) + l
    if (nn.ne.nz+np) return
    if (nn.ge.igfvbc) stop '[Basis CARSPH]: igfvbc too small'
    s1 = 0.0d0
    do k = 0,l
        if (2*k.ge.l+ma.and.k.lt.n+l) then   
            s1 = s1 + gfv%iv(k)*gfv%ibc(l,k)*gfv%ibc(2*k,l+ma)*gfv%ibc(n-1+l-k,nr) 
        endif
    enddo
    if(s1.eq.0.d0) return
    s2 = 0.0d0
    do k = 0,nr
        k2 = 2*k
        if (ny-ma.le.k2.and.k2.le.ny) then 
            s2 = s2 + gfv%iv(k)*gfv%ibc(nr,k)*gfv%ibc(ma,ny-k2)
        endif
    enddo
    if(s2.eq.0.d0) return
    aa = gfv%iv(n-1)*sqrt(2**(n+l+1)/(gfv%fak(n-1)*gfv%fad(n-1+l)))
    cc = gfv%iv(l-ma)*gfv%sqh(l)*gfv%wf(l+m)*gfv%wf(l-m)*gfv%fi(l)/2**l
    zz = gfv%wf(nx)*gfv%wf(ny)*gfv%wf(nz)*gfv%sqh(0)**(nn+1)
    carsph = aa*cc*zz*s1*s2
    if (m.lt.0) carsph = gfv%iv(ma+ny)*carsph
end function carsph

function carzyl(nx,ny,nz,nzz,nr,ml)    
    !-----------------------------------------------------------------------
    !     Transformation from cartesian oscillator to zylindrical oscillator
    !     yields <nx, ny, nz | nzz, nr, ml >

    !     | nx,ny,nz > = | nx > i**ny | ny > | nz >
    !       where | n > = pi**(-1/4) (2**n n!)**(-1/2) H_n(x) exp(-x**2/2)
    !
    !     | nz,nr,ml > = |nz> | nr,ml > exp( i ml phi )/sqrt(2 pi)
    !       where | nr,ml > = sqrt( 2 nr! / (nr+|ml|)!) r_p ** |ml|
    !                         *  L^|ml|_nr ( r_p**2 ) exp( - r_p**2/2)
    !----------------------------------------------------------------------
    use Constants, only: igfvbc
    use Globals, only: gfv
    implicit real*8(a-h,o-z)
    implicit integer(i-n)
    double precision carzyl
    carzyl = 0.0d0
    if (nz.ne.nzz) return
    ma = iabs(ml)
    np = 2*nr + ma
    if (np.ne.nx+ny) return
    if (np.ge.igfvbc) stop '[Basis CARZYL]: igfvbc too small'
    s = 0.0d0
    do k = 0,nr
        k2 = 2*k
        if (ny-ma.le.k2.and.k2.le.ny) s = s + gfv%iv(k)*gfv%ibc(nr,k)*gfv%ibc(ma,ny-k2) 
    enddo
    carzyl = gfv%iv(nr)*gfv%sqh(0)**np*gfv%wf(nx)*gfv%wf(ny)*gfv%wfi(nr)*gfv%wfi(nr+ma)*s
    if(ml.lt.0) carzyl = gfv%iv(ny)*carzyl
end function carzyl

! Spherical
subroutine set_Spherical_HO_basis(ifPrint)
    logical, intent(in), optional :: ifPrint
    call set_Spherical_HO_basis_parameters(ifPrint .and. .True.) 
    ! calculate quantum number of basis
    call calculate_Spherical_HO_basis_quantum_number(ifPrint .and. .False.)
    ! set gauss points in spherical coordinates
    call set_Spherical_gauss_points(ifPrint .and. .True.)
    ! radial and angular part wave function of spherical oscillator in Gauss points
    ! R_{nl}(x),  Y_{lm}(\theta,\phi) and their first-order derivation
    call calculate_Spherical_HO_basis_In_Gauss_points(ifPrint .and. .False.)

end subroutine set_Spherical_HO_basis

subroutine set_Spherical_HO_basis_parameters(ifPrint)
    !------------------------------------------------------------------------------!
    ! subroutine set_Spherical_HO_basis_parameters                                 !
    !                                                                              !
    ! Set parameters, ...                                                          !
    !------------------------------------------------------------------------------!
    use Constants, only: hbc, two, third, zero
    use Globals, only: input_par,force,nucleus_attributes
    logical,intent(in),optional :: ifPrint
    BS%HO_sph%n0f = input_par%basis_n0f
    BS%HO_sph%b0  = input_par%basis_b0
    BS%HO_sph%hb0 = hbc/(two * force%masses%amu)
    BS%HO_sph%hom = 41.0* nucleus_attributes%mass_number**(-third)   
    if (BS%HO_sph%b0 .le. zero) then
        BS%HO_sph%b0 = sqrt(two * BS%HO_sph%hb0/BS%HO_sph%hom)
    endif
    if(ifPrint) call printSphBasisParameters
    contains
    subroutine printSphBasisParameters
        use Globals, only: outputfile
        character(len=*),parameter :: format1="(1(a16,i5))", &
                                      format2="(a16,f12.5)"
        write(outputfile%u_config,outputfile%format_title)'set_Spherical_HO_basis_parameters'
        write(outputfile%u_config,format1) 'n0f    =',BS%HO_sph%n0f
        write(outputfile%u_config,format2) 'hb0    =',BS%HO_sph%hb0
        write(outputfile%u_config,format2) 'hom    =',BS%HO_sph%hom
        write(outputfile%u_config,format2) 'b0     =',BS%HO_sph%b0
    end subroutine printSphBasisParameters
end subroutine set_Spherical_HO_basis_parameters

subroutine calculate_Spherical_HO_basis_quantum_number(ifPrint)
    !----------------------------------------------------------------------------! 
    ! set Spherical Harmonic Oscillator Basis quantum number                     !
    ! |n_r l j m_j >
    ! Major quantum number:  N = 2*(n_r -1) + l                                  !
    ! nljm(:,1): n_r
    ! nljm(:,2): l
    ! nljm(:,3): j + 1/2 ; while |l - 1/2 |<= j <= |l + 1/2|
    ! nljm(:,4): m_j + 1/2; while -j <= m_j <= j
    !----------------------------------------------------------------------------!
    use Constants, only: nt3x
    logical, intent(in), optional :: ifPrint
    integer :: nrms,nlms,njms,nmms,nfx0,ngx0
    integer :: il,ik,ib,ifg,nn1,nr1,nl1,nj1,nm1
    character(len=*), parameter ::  format1 = "(1h[,4i2,1h],a)", &
                                    format2 = "(4i2)"
    nrms = 0
    nlms = 0
    njms = 0
    nmms = 0
    nfx0 = 0
    ngx0 = 0

    il = 0    
    ib = 1   
    do ifg = 1,2 ! loop over large and small components
        BS%HO_sph%iasp(ib,ifg) = il ! start position of this ib
        do nn1 = 0, BS%HO_sph%n0f+ifg-1, 1 !loop over major quantum number N.
            do nr1 = 1, nn1/2+1     !loop over radial quantum number nr.
                do nl1 = 0, nn1     !loop over orbital angular momentum l.
                    if(2*(nr1-1)+nl1 .ne. nn1) cycle
                    do nj1 = nl1, nl1+1 ! j=l+1/2 or l-1/2 ! nj=j+1/2=(l +or-1/2 + 1/2)=nl+1 or nl
                        do nm1 = -nj1+1, nj1 ! magnetic quantum number m_j
                            il = il+1
                            if(il.gt.nt3x) stop 'BASE : ntx too small'
                            BS%HO_sph%nljm(il,1)=nr1
                            BS%HO_sph%nljm(il,2)=nl1
                            BS%HO_sph%nljm(il,3)=nj1
                            BS%HO_sph%nljm(il,4)=nm1
                            BS%HO_sph%nps(il)   =ifg
                            if (nn1.lt.10) then
                                write(BS%HO_sph%tts(il),format1) nr1,nl1,nj1,nm1,' '
                            else
                                write(BS%HO_sph%tts(il),format2) nr1,nl1,nj1,nm1
                            endif
                            nrms = max0(nrms,nr1)
                            nlms = max0(nlms,nl1)
                            njms = max0(njms,nj1)
                            nmms = max0(nmms,nm1)
                        enddo !nm
                    enddo !nj
                enddo !nl
            enddo !nr
        enddo !nn
        BS%HO_sph%idsp(ib,ifg) = il - BS%HO_sph%iasp(ib,ifg) ! dimension length of this ib
    enddo !ifg 
    nfx0 = max(nfx0,BS%HO_sph%idsp(ib,1)) ! max dimension of these blocks in large component
    ngx0 = max(ngx0,BS%HO_sph%idsp(ib,2)) ! max dimension of these blocks in small component
    
    if(ifPrint) call printSphBasis
    contains
    subroutine printSphBasis
        use Globals, only: outputfile
        use Constants, only: nt_max
        integer :: nf,ng,nh,mf,i
        character(len=*), parameter ::  format3 = "(a,2i8)", &
                                        format4 = "(a,i8,a,i8)",&
                                        format5 = "(i4,a,i2,a,i2,a,i2,a,i3,a,i3,a)"
        write(outputfile%u_config,outputfile%format_title) 'calculate_Spherical_HO_basis_quantum_number'
        write(outputfile%u_config,format3) ' Maximal nr:       nrms = ',nrms
        write(outputfile%u_config,format3) ' Maximal nl:       nlms = ',nlms
        write(outputfile%u_config,format3) ' Maximal nj:       njms = ',njms
        write(outputfile%u_config,format3) ' Maximal nm:       nmms = ',nmms
        write(outputfile%u_config,format4) ' Dimension basis(large) = ',nfx0,"  <",nt_max*2
        write(outputfile%u_config,format4) ' Dimension basis(small) = ',ngx0,"  <",nt_max*2
        nf  = BS%HO_sph%idsp(ib,1)
        ng  = BS%HO_sph%idsp(ib,2)
        nh  = nf + ng
        mf  = BS%HO_sph%iasp(ib,1)
        do i = mf+1,mf+nh
            write(outputfile%u_config,format5) i,'   nr = ',BS%HO_sph%nljm(i,1),'   nl = ',BS%HO_sph%nljm(i,2), &
                                                '   nj = ',BS%HO_sph%nljm(i,3),'   nm = ',BS%HO_sph%nljm(i,4)
            if (i.eq.mf+nf) write(outputfile%u_config,'(3x,61(1h-))')
        enddo
    end subroutine printSphBasis                                                          
end  subroutine calculate_Spherical_HO_basis_quantum_number

subroutine set_Spherical_gauss_points(ifPrint)
    use Constants, only:  r64, pi, ngr, ntheta, nphi
    use MathMethods, only: GaussHermite,GaussLegendre_X1toX2
    use Globals, only : spatial_grid
    implicit none
    logical, intent(in), optional :: ifPrint
    integer :: ix,ithe,iphi,i
    real(r64) :: xx, theta, measure
    real(r64), dimension(ngr*2) :: x, wx

    ! gauss points 
    call GaussHermite(ngr*2,x,wx)
    do ix = 1, ngr
        spatial_grid%x(ix) = x(ix)
        spatial_grid%wx(ix) = wx(ix) * dexp(x(ix)**2)  
    end do
    call GaussLegendre_X1toX2(0.d0,   pi, spatial_grid%theta, spatial_grid%wtheta, ntheta)
    call GaussLegendre_X1toX2(0.d0, 2*pi, spatial_grid%phi,   spatial_grid%wphi,   nphi)

    ! weighted volume elements in spherical coordinates
    ! wwsp = b0^3 *x^2 * sin(theta) * wx * wtheta * wphi
    !      = r^2 sin(theta) * wr * wtheta * wphi
    do ix=1,ngr
        xx = spatial_grid%x(ix)**2
        do ithe=1,ntheta 
            theta=spatial_grid%theta(ithe)
            measure = xx*dsin(theta)
            do iphi=1,nphi
                i = ix + (ithe-1)*ngr + (iphi-1)*ngr*ntheta
                spatial_grid%wwsp(i) = BS%HO_sph%b0**3*spatial_grid%wx(ix)*spatial_grid%wtheta(ithe) &
                                        *spatial_grid%wphi(iphi)*measure
            enddo
        enddo
    enddo
    if(ifPrint) call printGaussPoints
    contains
    subroutine printGaussPoints
        use Globals, only: outputfile
        character(len=*), parameter :: format1= "(a30,15f9.5,/,30x,15f9.5,/,30x,15f9.5)"
        write(outputfile%u_config,outputfile%format_title) 'set_Spherical_gauss_points'
        write(outputfile%u_config,*) 'Gauss-Hermite:'
        write(outputfile%u_config,format1) 'x               :',spatial_grid%x
        write(outputfile%u_config,format1) 'w_x             :',spatial_grid%wx
        write(outputfile%u_config,*) 'Gauss-Legendree:'
        write(outputfile%u_config,format1) 'theta           :',spatial_grid%theta
        write(outputfile%u_config,format1) 'w_theta         :',spatial_grid%wtheta
        write(outputfile%u_config,format1) 'phi             :',spatial_grid%phi
        write(outputfile%u_config,format1) 'w_phi           :',spatial_grid%wphi    
    end subroutine printGaussPoints
end subroutine set_Spherical_gauss_points

subroutine calculate_Spherical_HO_basis_In_Gauss_points(ifPrint)
    !--------------------------------------------------------------------------------------
    ! Calculates the wavefunctions for the spherical oscillator.
    !     the are given as:  
    !     phi(r,theta,phi) = R_nl(r) * Y_lm(theta,phi) 
    !
    !     R_nl(r) = b^(-3/2) N_nl * x**l * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
    !     N_nl    = sqrt(2 * (n-1)!/Gamma(n+l+1/2)!)      
    !     where x=r/b, n=1,2,3,...
    !     R_nl is normalized in such a way that the norm integral reads
    !     \int dr r**2 R_nl(r)^2 = 1 
    !
    !    Y_lm (theta,phi) = sqrt((2l+1)/4*pi)*sqrt((l-m)!/(l+m)!)P_l^m(cos(theta)) e^(im*phi)
    !                      = sqrt((2l+1)/4*pi)d^l_{-m0}(theta) e^(im*phi)
    !    where d^l_{mm'}(theta) = <lm|e^{i*theta*J_y}|lm'>
    !-----------------------------------------------------------------------------------------
    !    Note: 
    !    1) In the program we store R_nl(x) not R_nl(r)
    !       R_nl(x) = N_nl * x**l * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
    !       (without b^{-3/2})
    !    2) R_nl and Y_lm are multiplied by sqrt(wx) and sqrt(wthet*wphi)
    !
    !    Output:
    !    BS%HO_sph%Rnl(x,n,l) = R_nl(x) * sqrt(wx)
    !    BS%HO_sph%Rnl1(x,n,l) = R'_nl(x) * sqrt(wx)
    !    BS%HO_sph%Ylm(phi,theta,l,m) = Y^_lm (theta,phi) * sqrt(wthet*wphi)
    !    BS%HO_sph%Ylm1(phi,theta,l,m) = d(Y^_lm (theta,phi))/d(theta) * sqrt(wthet*wphi)
    !-----------------------------------------------------------------------------------------
    use Constants, only: r64,half,ngr,ntheta,nphi,pi,zero,nr_max_sph,l_max_sph
    use Globals, only: spatial_grid, gfv, BS
    implicit none
    logical,intent(in) :: ifPrint
    integer :: nrx,lx,ix,l,n,itheta,iphi,m
    real(r64) :: x,xx,wx,theta,phi,cos_theta,cot_theta,wtp
    complex(r64) :: ctheta,cphi,ctemp1,ctemp2
    real(r64), dimension(:), allocatable :: rnl
    nrx = nr_max_sph
    lx = l_max_sph
    allocate(BS%HO_sph%Rnl(1:ngr,1:nrx,0:lx))
    allocate(BS%HO_sph%Rnl1(1:ngr,1:nrx,0:lx))
    allocate(BS%HO_sph%Ylm(1:nphi,1:ntheta,0:lx,-lx:lx))
    allocate(BS%HO_sph%Ylm1(1:nphi,1:ntheta,0:lx,-lx:lx))
    allocate(rnl(nrx))
    !----------------------------------------------------------------------------------------------
    !   x (r/b) dependence, contain the weight (but no r^2)
    !   BS%HO_sph%Rnl(x,n,l) = R_nl(x) * sqrt(w)
    !   BS%HO_sph%Rnl1(x,n,l) = R'_nl(x) * sqrt(w)
    !----------------------------------------------------------------------------------------------
    do ix = 1, ngr
        x  = spatial_grid%x(ix)
        wx = spatial_grid%wx(ix)
        xx = x**2 
        BS%HO_sph%Rnl1(ix,1,0)= gfv%sq(2)*gfv%wgi(1)*(-x)*exp(-half*xx)*dsqrt(wx)
        do l = 0, lx
            if(l.ge.1) then 
                BS%HO_sph%Rnl1(ix,1,l) = gfv%sq(2)*gfv%wgi(l+1)*(l*x**(l-1)-x**(l+1))*exp(-half*xx)*dsqrt(wx)   
            endif
            call R_nl(nrx,l,x,rnl)
            do n = 1, nrx
                BS%HO_sph%Rnl(ix,n,l) = rnl(n)*dsqrt(wx)
            enddo !n 
        enddo !l
        do l = 0, lx-1 
            do n = 2, nrx
                BS%HO_sph%Rnl1(ix,n,l)= (l/x-x)*BS%HO_sph%Rnl(ix,n,l) - 2.d0*gfv%sq(n-1)*BS%HO_sph%Rnl(ix,n-1,l+1)
            enddo !n 
        enddo !l
        ! l=lx case 
        call R_nl(nrx,lx+1,x,rnl)
        do n = 2, nrx  
            BS%HO_sph%Rnl1(ix,n,lx)= (lx/x-x)*BS%HO_sph%Rnl(ix,n,lx) - 2.d0*gfv%sq(n-1)*rnl(n-1)*dsqrt(wx) 
        enddo
    enddo
    !----------------------------------------------------------------------------------------------
    !   theta, phi dependence, contain the weight (but no sin(theta))
    !   BS%HO_sph%Ylm(phi,theta,l,m) = Y^_lm (theta,phi) * sqrt(wthet*wphi)
    !   BS%HO_sph%Ylm1(phi,theta,l,m) = d(Y^_lm (theta,phi))/d(theta) * sqrt(wthet*wphi)
    !----------------------------------------------------------------------------------------------
    do itheta  = 1,ntheta
        theta     = spatial_grid%theta(itheta)
        cos_theta = dcos(theta)
        cot_theta = cos_theta/dsin(theta)
        ctheta    = cmplx(0,theta) ! complex i*theta
        do iphi = 1,nphi
            phi  = spatial_grid%phi(iphi)
            cphi = cmplx(0,phi) ! complex i*phi
            wtp  = spatial_grid%wtheta(itheta)*spatial_grid%wphi(iphi) 
            do l= 0, lx
                do m= -l, l
                    BS%HO_sph%Ylm(iphi,itheta,l,m)=dsqrt((2*l+1)/(4*pi))*djmk(l,-m,0,cos_theta,0)*cdexp(m*cphi)*dsqrt(wtp)
                enddo !m
            enddo !l
            do l = 0, lx
                do m = -l, l
                    ctemp1 = m*cot_theta*BS%HO_sph%Ylm(iphi,itheta,l,m)
                    if(m.lt.l) ctemp2 = dsqrt(l*(l+1.d0)-m*(m+1.d0))*BS%HO_sph%Ylm(iphi,itheta,l,m+1)*cdexp(-cphi)
                    if(m.eq.l) ctemp2 = zero
                    BS%HO_sph%Ylm1(iphi,itheta,l,m)= ctemp1 + ctemp2
                enddo !m
            enddo !l
        enddo !phi
    enddo !theta
    if (ifPrint) call test_integration_of_Spherical_HO_basis
end subroutine calculate_Spherical_HO_basis_In_Gauss_points

subroutine test_integration_of_Spherical_HO_basis
    use Constants, only: r64, nr_max_sph,l_max_sph, ngr, ntheta, nphi, zero,pi
    use Globals, only: outputfile,BS,spatial_grid
    implicit none
    character(len=*), parameter ::  format1 = "(' G.Lag.: l =',i2,' n1 =',i3,'  n2 =',i3,f15.8)",&
                                    format2 = "(9x,'l=',i2, ' n2 =',i3,' t1=',f15.8,' t2=',f15.8)", &
                                    format3 = "('Y_(',i4,i4,')*Y_(',i4,i4,')=',2f15.8)"
    integer :: lx,nrx,l,n1,n2,ix,l1,l2,m1,m2,ithe,iphi
    real(r64) :: s,t1,t2,x,wx,theta,cos_theta,phi,wtp
    complex(r64) :: sy,cphi,test1,test2
    lx = l_max_sph 
    nrx = nr_max_sph 
    write(outputfile%u_config,outputfile%format_title) 'test_integration_of_Spherical_HO_basis'
    ! test for radial integration
    do l = 0,lx
        do n1 = 1,nrx
            do n2 = 1,n1
                s  = zero  
                t1 = zero
                t2 = zero
                do ix = 1,ngr ! try to use ngr = 39 to check
                    x = spatial_grid%x(ix)
                    wx = spatial_grid%wx(ix)
                    s  = s + x**2*BS%HO_sph%Rnl(ix,n1,l)*BS%HO_sph%Rnl(ix,n2,l) 
                    t1 = t1 + dcos(x)*BS%HO_sph%Rnl(ix,n2,l)*dsqrt(wx)
                    t2 = t2 - dsin(x)*BS%HO_sph%Rnl1(ix,n2,l)*dsqrt(wx)
                enddo !ix
                write(outputfile%u_config,format1) l,n1,n2,s
                ! write(outputfile%u_config,format2) l,n2,t1,t2   
            enddo !n2
        enddo !n1
    enddo 
    ! test angle part integration Ylm
    do l1=0,lx
    do l2=0,lx
    do m1=-l1,l1
    do m2=-l2,l2
            sy = zero 
            do ithe  = 1,ntheta ! try to use ntheta = 60 to check
                theta = spatial_grid%theta(ithe)
                do iphi = 1,nphi ! try to use nphi = 60 to check
                    sy = sy+dsin(theta)*conjg(BS%HO_sph%Ylm(iphi,ithe,l1,m1))*BS%HO_sph%Ylm(iphi,ithe,l2,m2)
                enddo !phi
            enddo ! theta
            if(abs(sy).gt.0.1) write(outputfile%u_config,format3) l1,m1,l2,m2,sy
    enddo !m2
    enddo !m1
    enddo !l2
    enddo !l1
    ! check for the derivative of ylm1 
    do ithe  = 1,ntheta
        theta = spatial_grid%theta(ithe)
        cos_theta= dcos(theta) 
        iphi = nphi
        phi  = spatial_grid%phi(iphi)
        cphi = cmplx(0,phi)
        wtp  = spatial_grid%wtheta(ithe)*spatial_grid%wphi(iphi) 
        test1 = 3*dsqrt(35.d0/(2*pi))*4*dsin(theta)**3*dcos(theta)*cdexp(4*cphi)/16.d0
        test2 = BS%HO_sph%Ylm1(iphi,ithe,4,4)/dsqrt(wtp)  
        write(outputfile%u_config,*) 'costhe=',cos_theta, 'test1=',test1, 'test2=',test2  
    enddo ! 
end subroutine test_integration_of_Spherical_HO_basis

subroutine R_nl(n,l,x,rnl)
    !--------------------------------------------------------------------------------------
    ! Calculates the radial part of the spherical harmonic oscillator wavefunction.
    ! The wavefunction is given by:
    !     R_nl(x) = N_nl * x**l * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
    ! where x = r/b, n=1,2,3,..., l=0,1,2,...
    ! N_nl    = sqrt(2 * (n-1)!/Gamma(n+l+1/2)!)
    !--------------------------------------------------------------------------------------
    use Constants, only: r64,one,half
    use Globals, only: gfv
    implicit none
    integer, intent(in) :: n, l
    real(r64), intent(in) :: x
    real(r64), intent(out) :: rnl(n)
    real(r64) :: xx,xl
    integer :: i
    xx = x*x 
    if (l.eq.0) then
        xl = one
    else
        xl = x**l
    endif
    rnl(1) = gfv%sq(2)*gfv%wgi(l+1)*exp(-half*xx)*xl
    rnl(2) = rnl(1)*(l+1.5d0-xx)*gfv%shi(l+1)
    do i = 3,n
       rnl(i) = ((2*i+l-2.5d0-xx)*rnl(i-1)-gfv%sq(i-2)*gfv%sqh(i-2+l)*rnl(i-2))*gfv%sqi(i-1)*gfv%shi(i-1+l)
    enddo
end subroutine R_nl

double precision function djmk(J,M,K,COSBET,IS) 
    !--------------------------------------------------------------------------------------
    ! Calculates the Wigner-d matrix element d^j_{mm'}(\theta) = <jm|e^{i*\theta*J_y}|jm'>
    ! where J_y is the y-component of the angular momentum operator.
    !
    ! IS = 0: integer values for     j = J,      m = M,     m' = K
    ! IS = 1: half integer valus for j = J-1/2 , m = M-1/2, m' = K-1/2
    ! COSBET = cos(beta)
    !--------------------------------------------------------------------------------------
    use Constants, only: r64
    use Globals, only: gfv
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    integer :: J, M, K, IS
    real(r64) :: COSBET
    djmk = 0.d0
    if (iabs(m).gt.j.or.iabs(k).gt.j) return
    if (abs(COSBET).gt.1.d0) stop 'DJMK: cos(beta) is larger than one'
    if (COSBET.eq.1.d0) then
        if (m.eq.k) djmk = 1.d0
        return
    endif
    if (COSBET.eq.-1.d0) then
        if (m.eq.is-k) djmk = dfloat(gfv%iv(j-m))
        return
    endif
    JMM = J-M
    JMK = J-K
    JPM = J+M
    JPK = J+K
    MPK = M+K
    IF (IS.EQ.1) THEN
        JPM = JPM-1
        JPK = JPK-1
        MPK = MPK-1
    ENDIF
    C2 = (1.d0+COSBET)/2
    S2 = (1.d0-COSBET)/2
    C  = SQRT(C2)
    S  = SQRT(S2)
    CS = C2/S2
    IA = MAX(0,-MPK)
    IE = MIN(JMM,JMK)
    X  = gfv%IV(JMM+IA)*C**(2*IA+MPK)*S**(JMM+JMK-2*IA)*gfv%FI(JMM-IA)*gfv%FI(JMK-IA)*gfv%FI(MPK+IA)*gfv%FI(IA)
    do i = ia,ie
        DJMK = DJMK + X
        X = -X*CS*(JMM-I)*(JMK-I)/((I+1)*(MPK+I+1))
    enddo
    DJMK = DJMK*gfv%WF(JPM)*gfv%WF(JMM)*gfv%WF(JPK)*gfv%WF(JMK)
    return
end function djmk

subroutine set_spherical_oscillator_wave_function
    !-----------------------------------------------------------------------------------------------------
    !
    !  calculates and store the wave functions for the spherical oscillator.
    !
    !  -----
    !  The wave functions of the spherical oscillator basis is 
    !    \Phi_{m}(r, \theta, \phi, \sigma) = |n_r l j m_j> 
    !        = \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s} R_{n_r l}(r) Y_{l m_l}(\theta,\phi) \chi_{m_{s}}^{\sigma}
    !   
    !   Am(m, i, m_s, ifg) = C^{j m_j}_{l m_l 1/2 m_s} R_{n_r l}(r) Y_{l m_l}(\theta,\phi)
    !      m=(n_r,l,j,m_j); 
    !      i=(r,theta,phi);
    !      m_s = 0 (m_s = -1/2), 1 (m_s = +1/2) 
    !      ifg= 1 (large component), 0 (small component);
    !  Note: m_j = m_l + m_s => m_l is fixed given m_j and m_s, removing the m_l and m_s summation.
    ! -----
    !  The second derivative of the wave function is 
    !   \nabla^2 \Phi_{m}(r, \theta, \phi, \sigma) 
    !               = \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s} 
    !                 \nabla^2(R_{n_r l}(r) Y_{l m_l}(\theta,\phi)) \chi_{m_{s}}^{\sigma}
    !               = \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s} 
    !                 2/b^2*[x^2/2 - (2n_r+l-1/2)]R_{n_r l}(r) Y_{l m_l}(\theta,\phi)) \chi_{m_{s}
    !   where x = r/b
    !
    !   Dm(m, i, m_s, ifg) = C^{j m_j}_{l m_l 1/2 m_s}*2/b^2*[x^2/2-(2n_r+l-1/2)]R_{n_r l}(r) Y_{l m_l}(\theta,\phi))
    !
    ! ----
    ! The first-order derivative of the wave function is
    !   \nabla \Phi_{m}(r, \theta, \phi, \sigma)
    !             = \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s}\nabla(R_{n_r l}(r) Y_{l m_l}(\theta,\phi)) \chi_{m_{s}}^{\sigma}
    !             = \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s} dR_{n_r l}_dr Y_{lm} \chi_{m_{s}}^{\sigma} e_r
    !             + \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s} 1/r*R_{n_r l} dY_{lm})_dtheta \chi_{m_{s}}^{\sigma} e_theta
    !             + \sum_{m_l m_s} C^{j m_j}_{l m_l 1/2 m_s} (im_l)/(r*sin\theta)R_{n_r l} Y_{lm})\chi_{m_{s}}^{\sigma} e_phi
    !
    !  Y1(m, i, m_s, ifg) = C^{j m_j}_{l m_l 1/2 m_s} dR_{n_r l}_dr Y_{lm}
    !  Y2(m, i, m_s, ifg) = C^{j m_j}_{l m_l 1/2 m_s} 1/r*R_{n_r l} dY_{lm})_dtheta
    !  Y3(m, i, m_s, ifg) = C^{j m_j}_{l m_l 1/2 m_s} (m_l)/(r*sin\theta)R_{n_r l} Y_{lm}) (divided by i)
    ! -----
    ! Note: 
    ! 1) In the program, R_nl(r) corresponds to R_{nl}(x) without the b^{-3/2} normalization factor. 
    !    This allows the substitution of r^2dr with x^2dx when computing integrals involving R^2_{nl}.
    ! 2) The measure  dsqrt(x**2*dsin(theta)) is mutiplied to the wave function.
    !-----------------------------------------------------------------------------------------------------
    use Constants, only: one,two,r64,nt3x,ngr,ntheta,nphi
    use Globals, only: BS,spatial_grid,option
    integer :: ifg,ib,begin_pos,basis_dim,m,nr,nl,nj,nm,ix,itheta,iphi,i
    real(r64) :: b0,cg0,cg1,x,theta,measure,temp1,temp2
    
    ! set R_{nl}(x),  Y_{lm}(\theta,\phi) and their first-order derivation
    if(option%CDFTType==0) then
        call set_Spherical_HO_basis(.False.)
    endif
    allocate(BS%HO_sph%Am(1:nt3x,1:ngr*ntheta*nphi,0:1,1:2))
    allocate(BS%HO_sph%Dm(1:nt3x,1:ngr*ntheta*nphi,0:1,1:2))
    allocate(BS%HO_sph%Y1(1:nt3x,1:ngr*ntheta*nphi,0:1,1:2))
    allocate(BS%HO_sph%Y2(1:nt3x,1:ngr*ntheta*nphi,0:1,1:2))
    allocate(BS%HO_sph%Y3(1:nt3x,1:ngr*ntheta*nphi,0:1,1:2))
    ! set spherical oscillator wave function
    b0 = BS%HO_sph%b0
    do ifg=1,2
        ib = 1
        begin_pos =  BS%HO_sph%iasp(1,ifg) 
        basis_dim =  BS%HO_sph%idsp(1,ifg)
        if(basis_dim > nt3x) stop 'nt3x is too small!'
        do m = 1, basis_dim ! loop of oscillator basis states m
            nr = BS%HO_sph%nljm(begin_pos+m,1) ! n_r
            nl = BS%HO_sph%nljm(begin_pos+m,2) ! l
            nj = BS%HO_sph%nljm(begin_pos+m,3) ! j + 1/2
            nm = BS%HO_sph%nljm(begin_pos+m,4) ! m_j + 1/2

            ! if(nm.gt.nl) then ! m_j + 1/2 > l
            !     ! m_l = l, m_s = +1/2, j = l+1/2, m_j = l+1/2
            !     cg0 = zero
            !     cg1 = clebslj(nl,nj,1,nm)
            ! elseif(nm.lt.-nl) then ! m_j + 1/2 < -l
            !     ! m_l =-l, m_s = -1/2, j = l-1/2, m_j =-l-1/2
            !     cg0 = clebslj(nl,nj,0,nm)
            !     cg1 = zero
            ! else
            !     cg0 = clebslj(nl,nj,0,nm) ! ms = -1/2
            !     cg1 = clebslj(nl,nj,1,nm) ! ms = +1/2
            ! endif

            cg0 = clebslj(nl,nj,0,nm) ! ms = -1/2
            cg1 = clebslj(nl,nj,1,nm) ! ms = +1/2
            do ix = 1, ngr
                x = spatial_grid%x(ix)
                temp1 = two/(b0**2)*(x**2/2.d0-(2*nr+nl-1.d0/2.d0)) ! 2/b^2[x^2/2-(2n_r +l-1/2)]
                do itheta = 1, ntheta
                    theta = spatial_grid%theta(itheta)
                    measure = dsqrt(x**2*dsin(theta))  
                    temp2 = one/(x*b0*dsin(theta)) ! 1/(r*sin(theta))
                    do iphi = 1, nphi 
                        i = ix + (itheta-1)*ngr + (iphi-1)*ngr*ntheta
                        BS%HO_sph%Am(m,i,0,ifg) = measure*cg0*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm) !ms=-1/2 
                        BS%HO_sph%Am(m,i,1,ifg) = measure*cg1*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm-1) !ms=+1/2
                        BS%HO_sph%Dm(m,i,0,ifg) = measure*cg0*temp1*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm) !ms=-1/2 
                        BS%HO_sph%Dm(m,i,1,ifg) = measure*cg1*temp1*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm-1) !ms=+1/2
                        BS%HO_sph%Y1(m,i,0,ifg) = measure*cg0/b0*BS%HO_sph%Rnl1(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm) !ms=-1/2
                        BS%HO_sph%Y1(m,i,1,ifg) = measure*cg1/b0*BS%HO_sph%Rnl1(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm-1) !ms=+1/2
                        BS%HO_sph%Y2(m,i,0,ifg) = measure*cg0/(b0*x)*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm1(iphi,itheta,nl,nm) !ms=-1/2
                        BS%HO_sph%Y2(m,i,1,ifg) = measure*cg1/(b0*x)*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm1(iphi,itheta,nl,nm-1) !ms=+1/2
                        BS%HO_sph%Y3(m,i,0,ifg) = measure*cg0*nm*temp2*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm) !ms=-1/2 ! divided by i
                        BS%HO_sph%Y3(m,i,1,ifg) = measure*cg1*(nm-1)*temp2*BS%HO_sph%Rnl(ix,nr,nl) &
                                                    *BS%HO_sph%Ylm(iphi,itheta,nl,nm-1) !ms=-1/2 ! divided by i
                    enddo !iphi
                enddo !ithe
            enddo !ix
        enddo !m
    enddo !ifg
end subroutine

END MODULE Basis