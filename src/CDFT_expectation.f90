!==============================================================================!
! MODULE Expectation                                                           !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Expectation
use Constants, only: r64,zero,itx,nghl,ngh,ngl,hbc
use Globals, only: expectations,gauss,densities,constraint,nucleus_attributes,BS,&
                   force,dirac,pairing,RHB_pairing
implicit none
real(r64),dimension(3) :: xn ! integral of scalar density: \int \rho_v(r) dr
real(r64),dimension(3) :: xs ! integral of vector density: \int \rho_s(r) dr
contains
subroutine calculate_expectation_DIR(ifPrint)
    !----------------------------------------------------------------
    !   calculates expectation values of operators
    !----------------------------------------------------------------
    logical,intent(in),optional :: ifPrint
    call multipole(ifPrint .and. .True.)
    call energy_DIR(ifPrint .and. .True.)
end subroutine

subroutine calculate_expectation_RHB(ifPrint)
    logical,intent(in),optional :: ifPrint
    call multipole(ifPrint .and. .True.)
    call energy_RHB
end subroutine 

subroutine multipole(ifPrint)
    !----------------------------------------------------------------
    !    calculate multipole moments
    !----------------------------------------------------------------
    use Constants, only: pi,radius_r0,third
    logical,intent(in),optional :: ifPrint
    integer :: it,ih,il
    real(r64) :: A_,N_,Z_,z,zz,rrp,rr,rq,ro,rh,x,fac0,fac0n,fac0p,r01,fac1,r02,fac2,r03,fac3,r04,fac4
    real(r64), dimension(3) :: r2,rd,qq,oo,hh,bet2,oot,het
    A_ = nucleus_attributes%mass_number
    N_ = nucleus_attributes%neutron_number
    Z_ = nucleus_attributes%proton_number
    do it = 1,3
        xn(it) = zero  ! \int \rho_v(r) dr
        xs(it) = zero  ! \int \rho_s(r) dr
        r2(it) = zero  ! root mean square radius
        rd(it) = zero  ! dipole
        qq(it) = zero  ! quadrupole
        oo(it) = zero  ! octupole
        hh(it) = zero  ! hexadecupole
    enddo 
    do ih = 1,ngh
        z  = gauss%zb(ih)
        zz = z**2
        do il = 1,ngl
            rrp = gauss%rb(il)**2                                                                        
            ! r^2                                        
            rr = zz + rrp                                                           
            ! for quadrupole moment                                              
            rq = 3*zz - rr                                                                                                            
            ! for octupole moment
            ro = 2*z*zz-3*z*rrp
            ! for hexadecupole moment                                            
            rh = 8*zz**2 - 24*zz*rrp + 3*rrp**2
                                    
            do it = 1,itx                                                
                x      = densities%rv(ih+(il-1)*ngh,it)                              
                xn(it) = xn(it) + x                                       
                r2(it) = r2(it) + x*rr                                      
                qq(it) = qq(it) + x*rq ! Q
                oo(it) = oo(it) + x*ro ! O
                hh(it) = hh(it) + x*rh ! H
                xs(it) = xs(it) + densities%rs(ih+(il-1)*ngh,it)                      
                rd(it) = rd(it) + x*z
            enddo                                                          
        enddo                                                             
    enddo                                                                                                                              
    do it = 1,itx                                                   
        r2(it) = sqrt(r2(it)/xn(it))                                   
    enddo                                                             
                                                                       
    xn(3) = xn(1) + xn(2)                                             
    xs(3) = xs(1) + xs(2)                                             
    r2(3) = sqrt((N_*r2(1)**2 + Z_*r2(2)**2)/A_)
    rd(3) = rd(1) + rd(2)                                             
    qq(3) = qq(1) + qq(2)     
    oo(3) = oo(1) + oo(2)
    hh(3) = hh(1) + hh(2)

    fac0  = 4*pi/(3*A_) 
    fac0n = 4*pi/(3*N_)
    fac0p = 4*pi/(3*Z_)

    r01   = radius_r0*A_**third !R0
    fac1  = sqrt(3/(4*pi))
    ! Q
    r02   = r01**2
    fac2  = sqrt(5/(16*pi))
    bet2(1) = fac0n*fac2*qq(1)/r02
    bet2(2) = fac0p*fac2*qq(2)/r02
    bet2(3) = fac0 *fac2*qq(3)/r02
    ! O
    r03    = r01**3
    fac3   = sqrt(7/(16*pi))
    oot(1) = fac0n*fac3*oo(1)/r03
    oot(2) = fac0p*fac3*oo(2)/r03
    oot(3) = fac0 *fac3*oo(3)/r03
    !H
    r04    = r01**4
    fac4   = sqrt(9/(256*pi))
    het(1) = fac0n*fac4*hh(1)/r04
    het(2) = fac0p*fac4*hh(2)/r04
    het(3) = fac0 *fac4*hh(3)/r04

    expectations%rc    = sqrt(r2(2)**2 + 0.64)  ! charge radius                                   
    expectations%rms   = r2(3) ! root mean square
    constraint%calq1 = fac0*fac1*rd(3)/r01 !
    expectations%beta2   = bet2(3) ! deformation parameter \beta
    constraint%calq2   = expectations%beta2  !
    expectations%beta3  = oot(3) !
    constraint%calq3  = expectations%beta3 !
    expectations%betg  = bet2(3) !
    expectations%beto  = oot(3) !
    expectations%qq2p  = fac2*qq(2) ! Q20_p
    expectations%qq3p  = fac3*oo(2) !                      
    expectations%dd0 = N_/A_*rd(2)-Z_/A_*rd(1) ! dipole moment: D0=e*(N/A*zp-Z/A*zn)
    expectations%ddz = rd(2)

    if(ifPrint) call printMultipole
    contains
    subroutine printMultipole
        use Constants,only: radius_r0,third
        use Globals,only: outputfile,nucleus_attributes
        integer :: i
        character(len=*), parameter :: format1 = "(a,3f15.6)"
        character(len=*), parameter :: format2 = "(a,15x,2f15.6)"
        write(outputfile%u_outputf,"(//,1x,A,/,1x,33(1h-))") 'Multipole:'
        write(outputfile%u_outputf,'(/,28x,a,8x,a,9x,a)') 'neutron','proton','total'
        ! particle number                                                   
        write(outputfile%u_outputf,format1) ' particle number .....',xn
        write(outputfile%u_outputf,format1) ' trace scalar density ',xs
        ! rms-Radius                                                        
        write(outputfile%u_outputf,format1) ' rms-Radius ..........',r2
        ! charge-Radius
        write(outputfile%u_outputf,format2) ' charge-Radius, R0....',expectations%rc, &
              radius_r0*nucleus_attributes%mass_number**third
        write(outputfile%u_outputf,*) ''

        ! quadrupole-moment                                                 
        write(outputfile%u_outputf,format1) ' quadrupole moment ...',qq
        write(outputfile%u_outputf,format1) '    <Q20> ............',(qq(i)*fac2,i=1,3)
        ! quadrupole deformation                                            
        write(outputfile%u_outputf,format1) ' beta  ...............',bet2
        write(outputfile%u_outputf,*) ''

        ! octupole moment                                             
        write(outputfile%u_outputf,format1) ' octupole moment .....',oo
        write(outputfile%u_outputf,format1) '    <Q30> ............',(oo(i)*fac3,i=1,3)
        ! octupole deformation
        write(outputfile%u_outputf,format1) ' oota ................',oot
        write(outputfile%u_outputf,*) ''

        ! hexadecupole moment
        write(outputfile%u_outputf,format1) ' hexadecupole moment .',hh
        write(outputfile%u_outputf,format1) '    <Q40> ............',(hh(i)*fac4,i=1,3)
        ! hexadecupole-deformation
        write(outputfile%u_outputf,format1) ' heta ................',het

        write(outputfile%u_outputf,*) '*************************END print_multipole **********************'
    end subroutine printMultipole
end subroutine multipole

subroutine energy_DIR(ifPrint)
    !----------------------------------------------------------------
    !   calculate energy
    !----------------------------------------------------------------
    use Globals, only: matrices
    logical,intent(in),optional :: ifPrint
    integer :: it,ih,il,k,ib,nf,ng,k1,k2,i2,i1
    real(r64) :: s,ekin
    real(r64) :: ept(3),epart(3),ekt(3),emes(4),enl,ecou,esig,eome,edel,erho,ecm
    ! single particle energy                                         
    do it = 1,itx                                                         
        epart(it) = zero
        do k = 1, dirac%nk(it) ! nk: the number of single particle eigenenergy                                                 
            epart(it) = epart(it) + dirac%ee(k,it)*pairing%vv(k,it)
        enddo
        epart(it) = epart(it) + constraint%c1x*constraint%calq1
        epart(it) = epart(it) + constraint%c2x*constraint%calq2
        epart(it) = epart(it) + constraint%c3x*constraint%calq3

    enddo                                                          
    epart(3) = epart(1) + epart(2)
                                                    
    ! E_kin: alpha*p
    do it=1,itx                                               
        ekin = zero                                                    
        do ib = 1,BS%HO_cyl%nb ! number of K blocks                                                   
            nf = BS%HO_cyl%id(ib,1) ! dimension of large components in ib block                                               
            ng = BS%HO_cyl%id(ib,2) ! dimension of small components in ib block                                               
            k1 = dirac%ka(ib,it) + 1 ! begining position of energies(or eigenstate) of ib block in dirac%ee                                         
            k2 = dirac%ka(ib,it) + dirac%kd(ib,it) !ending postion; kd: dimension of energies(or eigenstate) of ib block in dirac%ee                              
            if (k1.le.k2) then ! The beginning must be less than or equal to the end                                    
                do i2 = 1,nf                                             
                    do i1 = 1,ng                                          
                        s = zero                                           
                        do k = k1,k2                                       
                        s = s + dirac%fg(nf+i1,k,it)*dirac%fg(i2,k,it)*pairing%vv(k,it)     
                        enddo                          
                        ekin = ekin + s*matrices%sp(i1+(i2-1)*ng,ib)           
                    enddo                                                 
                enddo                                                    
            endif                                                       
        enddo                         
        ekt(it) = 2*ekin  + force%masses%amu*(xs(it)-xn(it))*hbc                   
    enddo                                                                                              
    ekt(3)   = ekt(1) + ekt(2)     

    ! field energies
    call efield(emes,enl,ecou)
    esig  = emes(1)
    eome  = emes(2)
    edel  = emes(3)
    erho  = emes(4)

    ! pairing energy                                                    
    do it = 1,itx
        ept(it) = zero 
        do il=1,ngl
            do ih=1,ngh
                k = ih+ngh*(il-1)
                ept(it) = ept(it) + densities%rkapp(k,it)*pairing%delq(k,it)*gauss%wdcor(ih,il)
            enddo
        enddo
        ept(it) = -ept(it)/2                                 
    enddo
    ept(3) = ept(1) + ept(2)                                          
       
    ! center–of–mass correction
    ecm = -0.75d0*BS%HO_cyl%hom   ! the estimate formulation from the simple harmonic oscillator shell model
    if(ifPrint .and. .True.) write (*,*) '1E-cm=',ecm
    if(ifPrint .and. .True.) call centmas(ecm) ! microscopic c.m. correction (calculated only at the final step)
    
    if(ifPrint .and. .True.) write (*,*) '2E-cm=',ecm
                                                                                                                              
    expectations%etot = ekt(3) + esig + eome + erho + ecou + enl + ept(3) + ecm ! E                 
    expectations%ea   = expectations%etot/nucleus_attributes%mass_number  ! E/A 
    expectations%ecm   = ecm !E_{CM}

    if(ifPrint) call printEnergy
    contains
    subroutine printEnergy
        use Globals,only: outputfile
        character(len=*), parameter :: format1 = "(a,3f15.6)"
        character(len=*), parameter :: format2 = "(a,30x,2f15.6)"
        write(outputfile%u_outputf,"(//,1x,A,/,1x,33(1h-))") 'Enegy:'
        write(outputfile%u_outputf,'(/,28x,a,8x,a,9x,a)') 'neutron','proton','total'
        ! single-particle energy
        write(outputfile%u_outputf,format1) ' Particle Energy .....',epart
        ! kinetic energy
        write(outputfile%u_outputf,format1) ' Kinetic Energy ......',ekt
        ! sigma energy
        write(outputfile%u_outputf,format2) ' E-sigma .............',esig
        ! nonlinear part sigma energy
        write(outputfile%u_outputf,format2) ' E-sigma non linear ..',enl
        ! omega energy
        write(outputfile%u_outputf,format2) ' E-omega .............',eome
        ! rho-energy
        write(outputfile%u_outputf,format2) ' E-rho ...............',erho                                                                                            
        ! Coulomb energy
        write(outputfile%u_outputf,format2) ' Coulomb direct ......',ecou
        ! pairing energy
        write(outputfile%u_outputf,format1) ' Pairing Energy ......',ept  
        ! center of mass correction
        write(outputfile%u_outputf,format2) ' E-cm ................',expectations%ecm
        
        ! total energy                                                      
        write(outputfile%u_outputf,format2) ' Total Energy ........',expectations%etot                                                              
        ! energy per particle                                               
        write(outputfile%u_outputf,format2) ' E/A .................',expectations%ea

    end subroutine printEnergy
end subroutine energy_DIR

subroutine efield(emes,enl,ecou)
    !----------------------------------------------------------------
    !   calculate field energies
    !   Ooutputs:
    !       emes(4) : energy expectation of free terms of mesons energies
    !       enl     : energy expectation of non-linear terms of (sigma)mesons
    !       ecou    : energy expectation of coulomb potential
    !----------------------------------------------------------------
    use Globals, only: force,fields
    use Constants, only: half,icou
    real(r64),intent(out) :: emes(4), enl, ecou
    integer :: m,i,il,ih
    real(r64) :: sigma,omega,delta,rho,f,b2,b3,s,rho_c

    ! meson-exchange module
    if (force%option%ipc.eq.0) then
       ! free terms
        sigma = zero
        omega = zero
        delta = zero
        rho   = zero
        do il = 1,ngl
            do ih = 1,ngh
                i = ih + ngh*(il-1)
                sigma = sigma + force%couplg%ggsig*force%ff(i,1,1)* fields%sigma(ih,il) *densities%ro(i,1) *gauss%wdcor(ih,il)
                omega = omega + force%couplg%ggome*force%ff(i,2,1)* fields%omega(ih,il) *densities%ro(i,2) *gauss%wdcor(ih,il)
                delta = delta + force%couplg%ggdel*force%ff(i,3,1)* fields%delta(ih,il) *densities%ro(i,3) *gauss%wdcor(ih,il)
                rho   = rho   + force%couplg%ggrho*force%ff(i,4,1)* fields%rho(ih,il)   *densities%ro(i,4) *gauss%wdcor(ih,il)
            enddo
        enddo
        emes(1) = half*hbc*sigma
        emes(2) = half*hbc*omega
        emes(3) = half*hbc*delta
        emes(4) = half*hbc*rho
        
        ! non-linear terms in sigma
        if (force%option%inl.gt.0) then
            f  = force%couplg%ggsig/force%couplm%gsig
            b2 = f**3*force%nonlin%g2/3
            b3 = f**4*force%nonlin%g3/2
            s = zero
            do il = 1,ngl
                do ih =1,ngh
                    s = s - (b2*fields%sigma(ih,il)**3 + b3*fields%sigma(ih,il)**4)*gauss%wdcor(ih,il)
                enddo
            enddo
            enl = half*hbc*s
        endif

    ! point-coupling module
    elseif (force%option%ipc.eq.1) then
        sigma = zero
        omega = zero
        delta = zero
        rho   = zero
        do il = 1,ngl
            do ih =1,ngh
                i = ih + ngh*(il-1)
                sigma = sigma + ( force%couplg%ggsig*force%ff(i,1,1) *densities%ro(i,1)**2 ) * gauss%wdcor(ih,il)
                omega = omega + ( force%couplg%ggome*force%ff(i,2,1) *densities%ro(i,2)**2 ) * gauss%wdcor(ih,il)
                delta = delta + ( force%couplg%ggdel*force%ff(i,3,1) *densities%ro(i,3)**2 ) * gauss%wdcor(ih,il)
                rho   = rho   + ( force%couplg%ggrho*force%ff(i,4,1) *densities%ro(i,4)**2 ) * gauss%wdcor(ih,il)
            enddo
        enddo
        ! derivative terms
        do il = 1,ngl
            do ih =1,ngh
                i = ih + ngh*(il-1)
                sigma = sigma + (force%coupld%ddsig*densities%ro(i,1)*densities%dro(i,1) ) * gauss%wdcor(ih,il)
                omega = omega + (force%coupld%ddome*densities%ro(i,2)*densities%dro(i,2) ) * gauss%wdcor(ih,il)
                delta = delta + (force%coupld%dddel*densities%ro(i,3)*densities%dro(i,3) ) * gauss%wdcor(ih,il)
                rho   = rho   + (force%coupld%ddrho*densities%ro(i,4)*densities%dro(i,4) ) * gauss%wdcor(ih,il)
            enddo
        enddo
        emes(1) = half*hbc*sigma
        emes(2) = half*hbc*omega
        emes(3) = half*hbc*delta
        emes(4) = half*hbc*rho
        enl = zero ! without non-linear terms
    else
        stop 'in EFIELD: ipc not properly defined'
    endif

    ! Coulomb energy
    ecou  = zero
    if (icou.ne.0) then
        do il =1,ngl
            do ih =1,ngh
                i = ih + ngh*(il-1)
                rho_c = half*(densities%ro(i,2) + densities%ro(i,4)) ! densities%rv(i,2)
                ecou  = ecou + fields%coulomb(ih,il)* rho_c *gauss%wdcor(ih,il)
            enddo
        enddo 
    endif
    ecou  = hbc*ecou/2
    return
end subroutine efield

subroutine centmas(ecm)
    !----------------------------------------------------------------
    !   calculates the center-of-mass correction 
    !----------------------------------------------------------------
    use Constants, only: one
    use Globals, only: gfv
    real(r64),intent(out) :: ecm
    integer :: it,ib,ifg,nd,i0,nfg,n2,nz2,nr2,ml2,ll,n1,nz1,nr1,ml1,i12,&
               k1,k2,k,il,jb,l1,l2,l,md,j0,mfg,ki,kj,ms2,ms1,ivu
               
    real(r64) :: ap1,ap2,az1,az2,htem,hcent,dterm,sa,s,tdx,tdy,tdz,&
                 r,rrp,a,b,td1a,td1b,td1c,dtemp,eta,etb,vv1,vv2,uu1,uu2,uv1,uv2,w,&
                 eaxp,eayp,eazp,ebxp,etax,etay,etaz,etbx,etaxa,etaxb,eterm

    ecm = zero

    ap1  = one/(BS%HO_cyl%b0*BS%HO_cyl%bp)
    ap2  = one/(BS%HO_cyl%b0*BS%HO_cyl%bp)**2
    az1  = one/(BS%HO_cyl%b0*BS%HO_cyl%bz)
    az2  = one/(BS%HO_cyl%b0*BS%HO_cyl%bz)**2
    ! loop over neutrons and protons
    do it =1,2 
        htem  = -one/(2. * force%masses%amu * nucleus_attributes%mass_number*hbc)
        hcent = htem*hbc**2

        ! canculate the dirac term   
        dterm = zero
        do ib = 1,BS%HO_cyl%nb ! loop over K-parity-blocks
            do ifg = 1,2 ! loop over contributions from large and small components
                nd  = BS%HO_cyl%id(ib,ifg)
                i0  = BS%HO_cyl%ia(ib,ifg)
                nfg = (ifg-1)* BS%HO_cyl%id(ib,1)
                do n2 = 1,nd ! loop of oscillator basis states n2
                    nz2 = BS%HO_cyl%nz(n2+i0)
                    nr2 = BS%HO_cyl%nr(n2+i0)
                    ml2 = BS%HO_cyl%ml(n2+i0)
                    ll  = ml2**2
                    do n1 = n2,nd ! loop of oscillator basis states n1
                        nz1 = BS%HO_cyl%nz(n1+i0)
                        nr1 = BS%HO_cyl%nr(n1+i0)
                        ml1 = BS%HO_cyl%ml(n1+i0)
                        if (ml1.ne.ml2) goto 20
                        i12 = 2-n2/n1

                        ! loop over different eigenstates of Dirac equation k
                        sa = zero
                        k1 = dirac%ka(ib,it) + 1 ! begining in k block
                        k2 = dirac%ka(ib,it) + dirac%kd(ib,it) ! ending in k block
                        do k = k1,k2
                            sa = sa + dirac%fg(nfg+n1,k,it) * dirac%fg(nfg+n2,k,it) *pairing%vv(k,it)/2.
                        enddo
                        s = sa*i12

                        ! loop over the mesh-points
                        tdx  = zero
                        tdy  = zero
                        tdz  = zero
                        if (nz1.eq.nz2) then
                            do il = 1,ngl
                                r = gauss%rb(il)     
                                rrp = r**2
                                tdy = tdy + BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql(nr2,ml2,il)*ll/rrp

                                a = gfv%sq(nr2+1)*gfv%sq(nr2+ml2+1)
                                b = gfv%sq(nr2)*gfv%sq(nr2+ml2)
                                td1a = a*BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql1(nr2+1,ml2,il)/rrp
                                td1b = BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql1(nr2,ml2,il)/rrp
                                if(nr2>0) then
                                    td1c = b*BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql1(nr2-1,ml2,il)/rrp
                                else
                                    td1c = 0
                                endif
                                tdx = tdx + (td1a - td1b - td1c)
                            enddo
                            if (nr1.eq.nr2) tdz = -az2*(2*nz2+1)/2.

                        else if (nz1.eq.nz2+2) then 
                            if (nr1.eq.nr2) tdz = az2*gfv%sq(nz2+1)*gfv%sq(nz2+2)/2.
                            
                        else if (nz1.eq.nz2-2) then
                            if (nr1.eq.nr2)  tdz = az2*gfv%sq(nz2)*gfv%sq(nz2-1)/2.
                        endif 
                        dtemp = (tdx - tdy + tdz)*s
                        dterm = dterm + dtemp
                20  enddo
                enddo
            enddo
        enddo
        ! end the calculation of the dirac term

        ! calculate the exchange terms     
        eta  = zero
        etb  = zero
        ! loop over K,K'-parity-blocks
        do ib = 1,BS%HO_cyl%nb
            do jb =1,BS%HO_cyl%nb
                ! loop over different eigenstates of Dirac eq k and k' in ib and jb block
                k1 = dirac%ka(ib,it) + 1
                k2 = dirac%ka(ib,it) + dirac%kd(ib,it)
                l1 = dirac%ka(jb,it) + 1
                l2 = dirac%ka(jb,it) + dirac%kd(jb,it)
                do k = k1,k2
                    do l = l1,l2
                        vv1  = pairing%vv(k,it)/2.
                        vv2  = pairing%vv(l,it)/2.
                        uu1  = one - vv1
                        uu2  = one - vv2
                        uv1  = sqrt(vv1*uu1)
                        uv2  = sqrt(vv2*uu2)
                        w    = vv1*vv2 + uv1*uv2

                        eaxp = zero
                        eayp = zero
                        eazp = zero
                        ebxp = zero

                        ! loop over contributions from large and small componets
                        do ifg = 1,2
                            nd  = BS%HO_cyl%id(ib,ifg)
                            md  = BS%HO_cyl%id(jb,ifg)
                            i0  = BS%HO_cyl%ia(ib,ifg)
                            j0  = BS%HO_cyl%ia(jb,ifg)
                            nfg = (ifg-1)*BS%HO_cyl%id(ib,1)
                            mfg = (ifg-1)*BS%HO_cyl%id(jb,1)
                            ki  = BS%HO_cyl%ikb(ib)
                            kj  = BS%HO_cyl%ikb(jb)

                            ! loop of oscillator basis states n2 and n1
                            do n2 = 1,md
                                nz2 = BS%HO_cyl%nz(n2+j0)
                                nr2 = BS%HO_cyl%nr(n2+j0)
                                ml2 = BS%HO_cyl%ml(n2+j0)
                                ms2 = BS%HO_cyl%ms(n2+j0)
                                do n1 = 1,nd
                                    nz1 = BS%HO_cyl%nz(n1+i0)
                                    nr1 = BS%HO_cyl%nr(n1+i0)
                                    ml1 = BS%HO_cyl%ml(n1+i0)
                                    ms1 = BS%HO_cyl%ms(n1+i0)

                                    etax = zero
                                    etay = zero
                                    etaz = zero
                                    etbx = zero          

                                    sa  =  dirac%fg(nfg+n1,k,it) * dirac%fg(mfg+n2,l,it)
                                    ! calculate the first term of exchange eta: same parity (++) or (--)
                                    if (ms1.eq.ms2) then
                                        ! etax
                                        if (ki.eq.kj) then
                                            etaxa = zero
                                            etaxb = zero      
                                            do il = 1,ngl
                                                etaxa = etaxa + BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql(nr2,ml2,il)
                                            enddo
                                            if (nz1.eq.nz2-1) etaxb = az1*gfv%sq(nz2)*gfv%sqi(2)
                                            if (nz1.eq.nz2+1) etaxb = -az1*gfv%sq(nz2+1)*gfv%sqi(2)
                                            etax = etaxa*etaxb
                                        ! etay
                                        else if (ki.eq.kj+1) then
                                            if (nz1.eq.nz2) then
                                                do il = 1,ngl
                                                    r = gauss%rb(il)
                                                    etay = etay + (BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql1(nr2,ml2,il)*gfv%sqi(2)/r&
                                                           -BS%HO_cyl%ql(nr1,ml1,il) * BS%HO_cyl%ql(nr2,ml2,il)*ml2*gfv%sqi(2)/r)
                                                enddo
                                            endif
                                        ! etaz
                                        else if (ki.eq.kj-1) then
                                            if (nz1.eq.nz2) then
                                                do il = 1,ngl
                                                    r = gauss%rb(il)
                                                    etaz = etaz + (BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql1(nr2,ml2,il)*gfv%sqi(2)/r&
                                                           +BS%HO_cyl%ql(nr1,ml1,il) * BS%HO_cyl%ql(nr2,ml2,il)*ml2*gfv%sqi(2)/r)
                                                enddo
                                            endif
                                        endif
                                        eaxp = eaxp + etax*sa
                                        eayp = eayp + etay*sa
                                        eazp = eazp + etaz*sa
                                        ! end the calculation of first term
                                    ! calculate the first term of exchange eta: opposite parity (+-) or (+-)
                                    else
                                        ivu =  -gfv%iv(ifg)*ms1
                                        ! etbx
                                        if (ki.eq.1.and.kj.eq.1) then
                                            if (nz1.eq.nz2) then
                                                do il = 1,ngl
                                                    r = gauss%rb(il)
                                                    etbx=etbx+ivu*(BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql1(nr2,ml2,il)*gfv%sqi(2)/r&
                                                         + BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql(nr2,ml2,il)*ml2*gfv%sqi(2)/r)
                                                enddo
                                            endif
                                        endif
                                        ebxp = ebxp + etbx*sa
                                    endif
                                enddo !n1
                            enddo !n2
                        enddo !ifg

                        eta = eta + w*(eaxp**2 + eayp**2 + eazp**2)
                        etb = etb + w*(ebxp**2)    
                    enddo !l
                enddo !k
            enddo !jb 
        enddo!ib
        eterm = eta + etb
        ! end the calculation of the exchange terms

        ecm = ecm - 2*hcent*(dterm + eterm)
    enddo !it
    return
end subroutine centmas

subroutine energy_RHB
    use Globals, only: matrices
    use Constants, only: NBX,NHHX
    integer :: it,ib,nf,ng,nh,m,n2,n1,il
    real(r64) :: emcc2,ekin
    real(r64),dimension(NHHX) :: h0
    real(r64) :: ept(3),epart(3),ekt(3),emes(4),enl,ecou,esig,eome,edel,erho,ecm
 
    ! single particle energy
    do it = 1,itx                                                         
        epart(it) = zero
        do ib = 1,BS%HO_cyl%nb
            nf  = BS%HO_cyl%id(ib,1)
            ng  = BS%HO_cyl%id(ib,2)
            nh  = nf + ng
            m   = ib + (it-1)*NBX
            epart(it) = epart(it) + trabt(nh,nh,nh,nh,dirac%hh(1,m),densities%rosh(1,m))
        enddo
        ! ! constrained energy
        epart(it) = epart(it) + constraint%c1x*constraint%calq1
        epart(it) = epart(it) + constraint%c2x*constraint%calq2
        epart(it) = epart(it) + constraint%c3x*constraint%calq3
    enddo                                                  
    epart(3) = epart(1) + epart(2)

    ! kinetic energy E_kin: alpha*p
    do it = 1,itx 
        ekin = zero
        do ib = 1,BS%HO_cyl%nb
            nf  = BS%HO_cyl%id(ib,1)
            ng  = BS%HO_cyl%id(ib,2)
            nh  = nf + ng
            m   = ib + (it-1)*NBX
            ! construction of the free Dirac-operator H0
            emcc2 = 2*force%masses%amu*hbc
            do n2 = 1,nf
                do n1 = 1,nf
                    h0(n1+(n2-1)*nh) = zero
                enddo
                do n1 = 1,ng
                    h0(nf+n1+(n2-1)*nh) = matrices%sp(n1+(n2-1)*ng,ib)
                    h0(n2+(nf+n1-1)*nh) = h0(nf+n1+(n2-1)*nh)
                enddo
            enddo
            do n2 = nf+1,nh
                do n1 = nf+1,nh
                    h0(n1+(n2-1)*nh) = zero
                enddo
                h0(n2+(n2-1)*nh) = - emcc2
            enddo
            ekin = ekin + trabt(nh,nh,nh,nh,h0,densities%rosh(1,m))
        enddo
        ekt(it) = ekin
    enddo
    ekt(3) = ekt(1) + ekt(2)

    ! field energies
    call efield(emes,enl,ecou)
    esig  = emes(1)
    eome  = emes(2)
    edel  = emes(3)
    erho  = emes(4)

    ! pairing energy
    do it = 1,itx
        ept(it) = zero
        il = 0
        do ib = 1,BS%HO_cyl%nb
            nf = BS%HO_cyl%id(ib,1)
            ng = BS%HO_cyl%id(ib,2)
            nh = nf + ng
            m  = ib + (it-1)*NBX
            do n2 = 1,nf
                do n1 = n2,nf
                    il = il + 1
                    ept(it)  = ept(it) + RHB_pairing%delta(n1+(n2-1)*nh,m)*RHB_pairing%kappa(il,it)
                enddo
            enddo
        enddo 
        ept(it) = -ept(it)/2
    enddo
    ept(3) = ept(1) + ept(2)
    
    ! center–of–mass correction
    ecm = -0.75d0*BS%HO_cyl%hom   ! the estimate formulation from the simple harmonic oscillator shell model

    expectations%etot = ekt(3) + esig + eome + erho + ecou + enl + ept(3) + ecm ! E                 
    expectations%ea   = expectations%etot/nucleus_attributes%mass_number  ! E/A 
    expectations%ecm   = ecm !E_{CM}

end subroutine energy_RHB

function trabt(ma,mb,n1,n2,aa,bb)
    real(r64) :: trabt
    integer :: ma,mb,n1,n2
    real(r64),dimension(ma,n2) :: aa
    real(r64),dimension(mb,n2) :: bb
    integer :: i,k
    real(r64) :: s
    s = 0.d0
    do i = 1,n1
        do k = 1,n2
            s = s + aa(i,k)*bb(i,k)
        enddo
    enddo
    trabt = s
end function

END MODULE Expectation