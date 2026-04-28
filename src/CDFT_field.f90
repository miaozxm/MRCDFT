!====================================================================================!
! MODULE Field                                                                       !
!                                                                                    !
! This module contains subroutines related to                                        !       
!                                                                                    !
!                                                                                    !
! List of subroutines and functions:                                                 !
! - subroutine                                                                       !
!====================================================================================!
MODULE Field
use Constants, only: r64,zero,one,third,half,pi,icou,hbc,alphi,nghl,ngh,ngl
use Globals, only: woodssaxon,gauss,constraint,nucleus_attributes,fields,&
                   outputfile
use CDFT_Inout, only: read_fields
implicit none    
 

contains
subroutine set_woodssaxon_parameters(ifPrint)
    use Globals,only: input_par
    logical, intent(in), optional :: ifPrint
    ! set input parameters
    woodssaxon%qs    = input_par%woodssaxon_qs
    woodssaxon%beta2 = input_par%woodssaxon_beta2
    woodssaxon%beta3 = input_par%woodssaxon_beta3
    ! Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
    ! W. Koepf and P. Ring Z. Phys. A 339, 81 (1990)
    ! 1: neutron and 2: proton
    woodssaxon%v0  = -71.28d0
    woodssaxon%akv = 0.4616d0
    woodssaxon%r0v = [1.2334d0,1.2496d0]
    woodssaxon%av  = [0.615d0,0.6124d0]
    woodssaxon%vso = [11.1175d0,8.9698d0]
    woodssaxon%rso = [1.1443d0,1.1401d0]
    woodssaxon%aso = [0.6476d0,0.6469d0]

    if(ifPrint) call printSaxonWoodParameters
    contains
    subroutine printSaxonWoodParameters
        character(len=*),parameter :: format1 = "(a,2f30.25)"
        write(outputfile%u_config,outputfile%format_title)  'set_woodssaxon_parameters'
        ! print parameters
        write(outputfile%u_config,*) 'Saxon-Woods Parameters : '
        write(outputfile%u_config,format1) ' v0     = ',woodssaxon%v0                                  
        write(outputfile%u_config,format1) ' kappa  = ',woodssaxon%akv                                 
        write(outputfile%u_config,format1) ' lambda = ',woodssaxon%vso                                 
        write(outputfile%u_config,format1) ' r0     = ',woodssaxon%r0v                                 
        write(outputfile%u_config,format1) ' a      = ',woodssaxon%av                                  
        write(outputfile%u_config,format1) ' r0-so  = ',woodssaxon%rso                                 
        write(outputfile%u_config,format1) ' a-so   = ',woodssaxon%aso                                 
        write(outputfile%u_config,format1) ' beta2  = ',woodssaxon%beta2                              
        write(outputfile%u_config,format1) ' beta3  = ',woodssaxon%beta3
    end subroutine printSaxonWoodParameters
end subroutine set_woodssaxon_parameters

subroutine initial_potential_fields(ifPrint)
    !-------------------------------------------------------------------------------------!
    ! given_initial_field                                                                       !
    !Purpose : initializes the potential.                                                 !
    !          init=0: reads fields from a file that already exists.                      !
    !          init=1: calculates fields in saxon-woods form.                             !
    !                                                                                     !
    !-------------------------------------------------------------------------------------!
    use Globals, only: input_par
    logical, intent(in), optional :: ifPrint
    ! set input parameters
    fields%inin = input_par%inin
    
    ! set initial fields
    if (fields%inin==0) then 
        ! read from input files
        call read_fields
    else if(fields%inin==1) then
        ! saxon-woods fields
        call calculate_saxonwoods_fields(ifPrint .and. .True.)
    else
        stop '[Potential]: inin must be 0 or 1'
    endif 
end subroutine initial_potential_fields

subroutine calculate_saxonwoods_fields(ifPrint)
    use Globals, only: BS
    logical, intent(in), optional :: ifPrint
    integer :: it,ih,il
    real(r64) :: z,zz,r,rr,cos_theta,Y20,Y30,facb, u,w,rc,tmp_c
    real(r64),dimension(2) :: rrv, rls, vp, vls 

    do it =1,2
        rrv(it) = woodssaxon%r0v(it)*nucleus_attributes%mass_number**third 
        rls(it) = woodssaxon%rso(it)*nucleus_attributes%mass_number**third 
        vp(it)  = woodssaxon%v0 * (one - woodssaxon%akv &
                *(-1)**it * ( nucleus_attributes%proton_number - nucleus_attributes%neutron_number) &
                /nucleus_attributes%mass_number)
        vls(it) = vp(it) * woodssaxon%vso(it)
        do ih =1,gauss%nh
            z  = gauss%zb(ih)
            zz = z**2
            do il = 1,gauss%nl
                rr = (gauss%xl(il)*(BS%HO_cyl%bp * BS%HO_cyl%b0)**2 + zz)
                r  = sqrt(rr)
                cos_theta = z/r ! cos(\theta)
                Y20 = sqrt(5/(16*pi))*(3*cos_theta**2-one)
                Y30 = sqrt(7/(16*pi))*(5*cos_theta**3-3*cos_theta)
                if(constraint%icstr == 0) then
                    facb = one + 0.75*(woodssaxon%beta2*Y20 + woodssaxon%beta3*Y30)
                elseif(constraint%icstr == 1) then
                    facb = one + 0.75*(constraint%betac(constraint%index)*Y20 + woodssaxon%beta3*Y30)
                elseif(constraint%icstr == 2) then
                    facb = one + 0.75*(constraint%betac(constraint%index)*Y20 + constraint%bet3c(constraint%index)*Y30)
                elseif(constraint%icstr == 3) then
                    facb = one + 0.75*(woodssaxon%beta2*Y20 + constraint%bet3c(constraint%index)*Y30)
                else 
                    stop "[Potential]: constraint%icstr must be 0,1,2,3"
                endif
                u= vp(it)/(one+exp((r-rrv(it)*facb)/woodssaxon%av(it))) 
                w= -vls(it)/(one+exp((r-rls(it)*facb)/woodssaxon%aso(it)))
            
                ! set fields
                fields%vps(ih,il,it) = u
                fields%vms(ih,il,it) = w 
                fields%vpstot(ih,il,it) = u
                fields%vmstot(ih,il,it) = w
                !  pairing%delq(ih,il,it) = pairing%del(it) ! samething have been done in moudule_preparation's set_pairing_parameters

                ! Coulomb, approximate the coulomb potential
                fields%coulomb(ih,il) = zero
                if(icou /= 0 .and. it==2) then ! only for proton
                    rc = rrv(2)
                    if(r < rc ) then
                        tmp_c = half * (3/rc-r*r/(rc**3))
                    else 
                        tmp_c = one/r
                    endif
                    fields%coulomb(ih,il) = tmp_c * nucleus_attributes%proton_number/alphi
                    fields%vps(ih,il,2) = fields%vps(ih,il,2) + hbc*fields%coulomb(ih,il)
                    fields%vms(ih,il,2) = fields%vms(ih,il,2) + hbc*fields%coulomb(ih,il)
                    fields%vpstot(ih,il,2) = fields%vpstot(ih,il,2) + hbc*fields%coulomb(ih,il)
                    fields%vmstot(ih,il,2) = fields%vmstot(ih,il,2) + hbc*fields%coulomb(ih,il)
                endif
            enddo
        enddo
    enddo

    if(ifPrint) call printSaxonWoodsField
    contains
    subroutine printSaxonWoodsField
        write(outputfile%u_outputf,*)   '***********BEGIN generate_initial_fields/calculate_saxonwoods_fields ******************'
        write(outputfile%u_outputf,*) 'Calculates fields in saxon-woods form !'
        ! coulomb force
        if(icou == 0) then
            write(outputfile%u_outputf,*) 'Without Coulomb force'
        else if (icou == 1) then
            write(outputfile%u_outputf,*) 'With Coulomb force'
        else if (icou == 2) then
            write(outputfile%u_outputf,*) 'With Coulomb force with exchange'
        else 
            write(outputfile%u_outputf,*) '[Preparation]: Wrong icou !'
        endif
        ! print fields
        write(outputfile%u_outputf,"(/,a)") 'Caculated Fields:'
        call prigh(outputfile%u_outputf,1,fields%coulomb,one,'Coulomb')
        call prigh(outputfile%u_outputf,1,fields%vps,one,'VPS-n')
        call prigh(outputfile%u_outputf,1,fields%vps(1,1,2),one,'VPS-p')
        call prigh(outputfile%u_outputf,1,fields%vms,one,'VMS-n')
        call prigh(outputfile%u_outputf,1,fields%vms(1,1,2),one,'VMS-p')
        call prigh(outputfile%u_outputf,1,fields%vpstot,one,'VPStotal-n')
        write(outputfile%u_outputf,"(a,/)")  '***********END generate_initial_fields/calculate_saxonwoods_fields *****************'
    end subroutine printSaxonWoodsField
end subroutine calculate_saxonwoods_fields
subroutine prigh(u_output,is,ff,f,text) 
    !---------------------------------------------------------------!
    ! IS = 1:  prints f*ff(z,r) at gauss-meshpoints                 !
    ! IS = 2:  prints f*ff(z,r)/wdcor at gauss-meshpoints           !
    ! (is=2) is used to output density                              !
    !  wdcor is the overall normalization factor equal to,          !
    !           wdcor=pi*WH(ih)*WL(il)* b0**3 * bz*bp*bp            !
    !---------------------------------------------------------------!
    integer :: u_output,is
    real(r64) :: f
    character(len=*) :: text
    real(r64),dimension(gauss%nh,gauss%nl) :: ff

    integer :: ixh =7, ixl = 3 ! change it to print more data; ixh should less than gauss%nh, ixl should less than gauss%nl
    integer :: ih,il
    real(r64) :: r 
    write(u_output,'(/,1x,a6,12f30.25)') text, (gauss%zb(ih),ih=1,ixh)
    do il = 1,ixl
        r = gauss%rb(il)
        if(is==1) then
            write(u_output,'(1x,f6.3,12f30.25)') r, (f*ff(ih,il),ih=1,ixh)
        elseif(is==2) then
            write(u_output,'(1x,f6.3,12f30.25)') r, (f*ff(ih,il)/gauss%wdcor(ih,il),ih=1,ixh)
        else
            write(u_output,*)'Wrong is in prigh()'
        endif
    enddo
end subroutine prigh

subroutine calculate_fields
    !----------------------------------------------------------------
    !       calculation of the coulomb fiels and meson-fields in the 
    ! oscillator basis. The fields are given in (fm^-1).
    !    meson fields:  sig(i) = phi(i,1)*ggsig/gsig
    !                   ome(i) = phi(i,2)*ggome/gome
    !                   del(i) = phi(i,3)*ggdel/gdel
    !                   rho(i) = phi(i,4)*ggrho/grho
    !----------------------------------------------------------------
    use Globals, only: densities, force 
    integer :: imes, i, is, ih, il
    real(r64), dimension(nghl) :: so, ph
    real(r64) :: f, sb, sc, si1, si2, si3, ss, sv, so1, so2, so3

    integer :: maxs = 1000
    real(r64) :: xmixsig = 1.d0, epss = 1.d-9

    ! calculate coulomb field
    if(icou.eq.1) call coulom(densities%drvp, fields%coulomb)
    ! point-coupling models does not require solving K-G equations
    if(force%option%ipc.eq.1) return

    !meson-exchange models need to solve K-G equationss
    ! sigma-meson
    imes = 1
    if (force%option%inl.eq.0) then
        do i = 1,nghl
            so(i) = force%ff(i,imes,1)*densities%ro(i,imes)
        enddo
        call gordon(imes,so,fields%sigma)
    else
        f  = force%couplg%ggsig / force%couplm%gsig
        sb = force%nonlin%g2 / force%couplm%gsig
        sc = force%nonlin%g3 / force%couplm%gsig
        ! loop of sigma-iteration
        do is = 0,maxs
            do il = 1,ngl
                do ih =1,ngh
                    i = ih+(il-1)*ngh
                    si1 = fields%sigma(ih,il) *f
                    si2 = si1*si1
                    si3 = si2*si1
                    so(i) = densities%ro(i,imes) + sb*si2 + sc*si3
                enddo
            enddo
            call gordon(imes,so,ph)
            ss = zero
            do il = 1, ngl
                do ih = 1, ngh
                    i = ih+(il-1)*ngh
                    sv = ph(i) - fields%sigma(ih,il)
                    fields%sigma(ih,il) = fields%sigma(ih,il) + xmixsig*sv
                    ss = dmax1(ss,abs(sv))
                enddo
            enddo
            if (ss.lt.epss) goto 10
        enddo
        stop 'FIELD: sigma-iteration has not converged'
    endif
    10 continue

    ! omega-meson
    imes = 2
    if (force%option%inl.eq.0 .or. abs(force%nonlin%c3).lt.1.d-8) then 
        do i = 1,nghl
        so(i) = force%ff(i,imes,1)*densities%ro(i,imes)
        enddo
        call gordon(imes,so,fields%omega)
    else
        f  = force%couplg%ggome / force%couplm%gome
        sc = -force%nonlin%c3 / force%couplm%gome
        ! loop of omega-iteration
        do is = 0,maxs
            do il = 1,ngl
                do ih = 1,ngh
                    i = ih+(il-1)*ngh
                    so1 = fields%omega(ih,il)*f
                    so2 = so1*so1
                    so3 = so2*so1
                    so(i) = densities%ro(i,imes) + sc*so3
                enddo
            enddo
            call gordon(imes,so,ph)
            ss = zero
            do il = 1,ngl
                do ih =1,ngh
                    i = ih+(il-1)*ngh
                    sv       = ph(i) - fields%omega(ih,il)
                    fields%omega(ih,il) = fields%omega(ih,il) + xmixsig*sv
                    ss = dmax1(ss,abs(sv))
                enddo
            enddo
            if (ss.lt.epss) goto 11
        enddo
        stop 'FIELD: omega-iteration has not converged'
    endif
    11 continue 

    ! delta-meson
    imes = 3
    do i = 1,nghl
        so(i) = force%ff(i,imes,1)*densities%ro(i,imes)
    enddo
    call gordon(imes,so,fields%delta)

    ! rho-meson
    imes = 4
    do i = 1,nghl
        so(i) = force%ff(i,imes,1)*densities%ro(i,imes)
    enddo
    call gordon(imes,so,fields%rho)
    return
end subroutine calculate_fields

subroutine calculate_meson_propagators(ifPrint) 
    !-------------------------------------------------------------------
    ! Calculates the meson-propagators GG.                              
    ! DD =/nabla^2 is the Laplace operator in oscillator space                   
    ! GG = m**2/(-DD+m**2)
    ! Actually, what is calculated is (-DD+m**2)/m**2 in Oscillator Basis
    !--------------------------------------------------------------------
    use Constants, only: ntb_max,zero,one,half
    use Globals, only: force,gfv,matrices,BS
    logical,intent(in),optional :: ifPrint
    integer :: i2,nz2,nr2,i1,nz1,nr1,no,imes,k,i,ifl
    real(r64),dimension(4) :: mass
    real(r64),dimension(ntb_max*ntb_max) ::dd,aa
    real(r64),dimension(ntb_max*ntb_max,4) :: gg
    real(r64) :: az2,ap2,t,m2i,d

    ! point-coupling models does not require solving K-G equations, thus no meson propagators
    if(force%option%ipc.eq.1) return
    
    mass(1) = force%masses%amsig
    mass(2) = force%masses%amome
    mass(3) = force%masses%amdel
    mass(4) = force%masses%amrho

    az2 = one/(BS%HO_cyl%b0*BS%HO_cyl%bz)**2
    ap2 = one/(BS%HO_cyl%b0*BS%HO_cyl%bp)**2

    no = BS%HO_cyl%ntb
    do i2 = 1,BS%HO_cyl%ntb
        nz2 = BS%HO_cyl%nzb(i2)
        nr2 = BS%HO_cyl%nrb(i2)
        dd(i2+(i2-1)*no) = -(az2*(nz2+half)+ap2*(nr2+nr2+1))
        do i1 = 1, i2-1
            nz1 = BS%HO_cyl%nzb(i1)
            nr1 = BS%HO_cyl%nrb(i1)
            t   = zero
            if (nr1.eq.nr2) then
               if (nz2.eq.nz1+2) t = az2*half*gfv%sq(nz1+1)*gfv%sq(nz2)
               if (nz2.eq.nz1-2) t = az2*half*gfv%sq(nz1)*gfv%sq(nz2+1)
            endif
            if (nz1.eq.nz2) then
               if (nr2.eq.nr1+1) t = -ap2*nr2
               if (nr2.eq.nr1-1) t = -ap2*nr1
            endif
            dd(i1+(i2-1)*no) = t
            dd(i2+(i1-1)*no) = t
        enddo
    enddo
                                                    
    do imes = 1,4                                                     
        m2i = (one/(mass(imes)+1.d-10))**2                                     
        do k = 1,no                                                    
            do i = 1,no                                                 
                aa(i+(k-1)*no) = -dd(i+(k-1)*no)*m2i              
                gg(i+(k-1)*no,imes) = zero                               
            enddo                                                       
            aa(k+(k-1)*no) = aa(k+(k-1)*no) + one                  
            gg(k+(k-1)*no,imes) = one                                
        enddo                                                          
        call lingd(no,no,no,no,aa,gg(1,imes),d,ifl)
    enddo
    if(ifl == -1) stop '[Field]: error in subroutine lingd!'
    matrices%meson_propagators(:,:) = gg(:,:)

    if(ifPrint) call printMesonPropagator
    contains

    subroutine printMesonPropagator
        use Globals,only: outputfile
        use DiracEquation, only: aprint
        integer :: i
        write(outputfile%u_config,outputfile%format_title) 'calculate_meson_propagators'
        do imes =1,4
            call aprint(outputfile%u_config,1,1,1,no,no,no,gg(1,imes),(/('',i=1,no)/),(/('',i=1,no)/),'Meson-Propagator')
        enddo
    end subroutine printMesonPropagator

end subroutine calculate_meson_propagators
subroutine lingd(ma,mx,n,m,a,x,d,ifl)
    !----------------------------------------------------------------------
    ! Solves the system of linear equations A*X = B .
    ! At the beginning, the matrix B is stored in X. During the calculation 
    ! it will be overwritten. D is the determinant of A .
    !----------------------------------------------------------------------
    integer :: ma, mx, n, m, ifl
    real(r64),dimension(ma,m) :: a
    real(r64),dimension(mx,m) :: x
    real(r64) :: d

    integer :: i,j,k,l,k1,n1
    real(r64) :: p,q,tol,cp,cq

    ! constant
    real(r64) :: tollim = 1.d-10, one = 1.d0, zero = 0.d0
    ifl = 1
    p = zero
    do i = 1,n
        q = zero
        do j = 1,n
        q = q + abs(a(i,j))
        enddo
        if (q.gt.p)   p=q
    enddo
    tol = tollim*p
    d   = one
    do k = 1,n
        p = zero
        do j = k,n
            q = abs(a(j,k))
            if (q.ge.p) then
                p = q
                i = j
            endif
        enddo
        if (p.le.tol) then
            ! write (*,100) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
            ! 100  format(1x,80a1,' *****  ERROR IN LINGD , TOLERANZ =', e10.4,' VALUE OF A(',i3,',',i3,') IS ',e10.4/1x,80a1) ! check it !
            ifl = -1
            return
        endif
        cp = one/a(i,k)
        if (i.ne.k) then
            d = -d
            do l = 1,m
                cq     = x(i,l)
                x(i,l) = x(k,l)
                x(k,l) = cq
            enddo
            do l = k,n
                cq     = a(i,l)
                a(i,l) = a(k,l)
                a(k,l) = cq
            enddo
        endif
        d = d*a(k,k)
        if (k.eq.n) goto 1
        k1 = k+1
        do i = k1,n
            cq = a(i,k)*cp
            do l = 1,m
                x(i,l) = x(i,l) - cq*x(k,l)
            enddo
            do l = k1,n
                a(i,l) = a(i,l) - cq*a(k,l)
            enddo
        enddo
    enddo
  1 do l = 1,m
        x(n,l) = x(n,l)*cp
    enddo
    if (n.gt.1) then
        n1 = n-1
        do k=1,n1
            cp = one/a(n-k,n-k)
            do l=1,m
                cq = x(n-k,l)
                do i = 1,k
                    cq = cq - a(n-k,n+1-i)*x(n+1-i,l)
                enddo
                x(n-k,l) = cq*cp
            enddo
        enddo
    endif
    return 
end subroutine lingd

subroutine coulom(dro,cou)
    !---------------------------------------------------------------------
    !                                                                       
    !     Coulom-Field (direct part) in zylindric coordinates               
    !     in units of fm^(-1)                                              
    !                                                                       
    !     cou(r) = int d^3r' rho_p(r') / |r-r'|
    !            = 1/2 int d^3r |r-r'| Delta rho(r')
    !   
    !    The expansion coefficient of the elliptic integral get from:
    !     https://personal.math.ubc.ca/~cbm/aands/frameindex.htm                         
    !----------------------------------------------------------------------
    use Constants, only: nghl,two,ngh,ngl
    real(r64), dimension(nghl) :: dro, cou
    real(r64) :: f, r, r4r, z, s, r1, r2r, r4s, z1, drop , y, x, xx, e
    integer :: il, ih, kl, kh


    f = one/(alphi*two*pi)
    do il = 1,ngl
        r  = gauss%rb(il)
        r4r = -4*r
        do ih = 1,ngh
            z = gauss%zb(ih)
            s = zero
            do kl = 1,ngl
                r1  = gauss%rb(kl)
                r2r = (r+r1)**2
                r4s = r4r*r1
                do kh = 1,ngh
                    z1 = gauss%zb(kh)
                    drop = dro(kh+(kl-1)*ngh)
                    y = (z-z1)**2+r2r
                    x  = one + r4s/y
                    xx = x*x
                    ! expansion of the eliptic integral
                    e = one + 0.44325141463*x + 0.06260601220*xx + 0.04757383546*xx*x + 0.01736506451*xx*xx
                    if (x.gt.1d-9) then 
                        e = e - dlog(x) *( 0.24998368310*x + 0.09200180037*xx + 0.04069697526*xx*x + 0.00526449639*xx*xx)
                    endif
                    s = s + drop*e*sqrt(y)
                enddo
            enddo
            cou(ih+(il-1)*ngh) = 2.d0*f*s
      enddo
    enddo
    return
end subroutine coulom

subroutine gordon(imes, so, phi)
    !----------------------------------------------------------------------
    !                                                                       
    !     solution of the Klein-Gordon-equation.                                 
    !     by expansion in zylindrical oscillator
    !     
    !     Input:                              
    !       imes:  number of meson                                            
    !       so:    source
    !     Output:                                                     
    !       phi:   meson field                                                
    !                                                                       
    !----------------------------------------------------------------------
    use Constants, only: ntb_max
    use Globals, only: BS, gauss, matrices
    integer, intent(in) :: imes
    real(r64), dimension(nghl), intent(in) :: so
    real(r64), dimension(nghl), intent(out) :: phi

    real(r64), dimension(ntb_max) :: rn,pn
    integer :: n, nr, nz, il, ih, k
    real(r64) :: s

    ! transformation of the inhomogeneous part to the oscillator base
    do n = 1,ntb_max
        nr = BS%HO_cyl%nrb(n)
        nz = BS%HO_cyl%nzb(n)
        s = zero
        do il = 1,ngl
            do ih = 1,ngh
                s = s + BS%HO_cyl%qhb(nz,ih) * BS%HO_cyl%qlb(nr,il) * so(ih+(il-1)*ngh) * gauss%wdcor(ih,il)
            enddo
        enddo
        rn(n) = s
    enddo

    ! multiplication with the Greens function
    do n = 1,ntb_max
        s = zero
        do k = 1,ntb_max
            s = s + matrices%meson_propagators(n+(k-1)*ntb_max,imes)*rn(k)
        enddo
        pn(n) = s
    enddo

    ! transformation back to the coordinate space
    do il = 1,ngl
        do ih = 1,ngh
            s = zero
            do n = 1,ntb_max
                nz = BS%HO_cyl%nzb(n)
                nr = BS%HO_cyl%nrb(n)
                s = s + pn(n)*BS%HO_cyl%qhb(nz,ih)*BS%HO_cyl%qlb(nr,il)
            enddo
            phi(ih+(il-1)*ngh) = s
        enddo
    enddo
    return
end subroutine

END MODULE Field