!==============================================================================!
! MODULE Density                                                               !
!                                                                              !
! This module calculates the densities in r-space at Gauss-meshpoints          !                                                 !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Density
use Constants, only: r64
use Globals, only: densities
implicit none
contains
subroutine calculate_densities_DIR(ifPrint)
    use Constants, only: zero,one,nghl,itx,ngl,ngh
    use Globals, only: BS,gfv,dirac,pairing,gauss
    logical,intent(in),optional :: ifPrint
    integer :: i,it,ib,nh,ifg,nd,i0,nfg,ivv,&
               n2,nz2,nr2,ml2,n1,nz1,nr1,ml1,&
               il2,k1,k2,k,il,ih,ihl
    real(r64) :: ap2,az2,s,sc,ben,qlab,qltab,bfl,qhab,sro,fgr,qh1ab,stau,sdro,sdt
    real(r64),dimension(2) :: rshell, cshell
    real(r64), dimension(nghl,2) :: rs  ! \rho_s
    real(r64), dimension(nghl,2) :: drs ! \Delta \rho_s
    real(r64), dimension(nghl,2) :: rv  ! \rho_v
    real(r64), dimension(nghl,2) :: drv ! \Delta \rho_v
    real(r64), dimension(nghl,2) :: rkapp ! bcs \kappa
    ap2 = one/(BS%HO_cyl%b0 * BS%HO_cyl%bp)**2
    az2 = one/(BS%HO_cyl%b0 * BS%HO_cyl%bz)**2
    ! initialize to zero
    do i=1,nghl
        do it =1,2
            rs(i,it)  = zero
            drs(i,it) = zero
            rv(i,it)  = zero
            drv(i,it) = zero
            rkapp(i,it) = zero ! bcs kappa
        enddo
    enddo

    do ib=1,BS%HO_cyl%nb ! loop over K-parity-blocks
        nh = BS%HO_cyl%id(ib,1) + BS%HO_cyl%id(ib,2)
        do ifg = 1,2 ! loop over contributions form large and small components
            nd = BS%HO_cyl%id(ib,ifg)
            i0 = BS%HO_cyl%ia(ib,ifg)
            nfg = (ifg-1)*BS%HO_cyl%id(ib,1)
            ivv = gfv%iv(ifg) ! iv(n) = (-1)**n
            do n2 = 1,nd
                nz2 = BS%HO_cyl%nz(n2+i0)
                nr2 = BS%HO_cyl%nr(n2+i0)
                ml2 = BS%HO_cyl%ml(n2+i0)
                do n1 = n2,nd ! due to the contribution is symmetric and this loop only cover half 
                    nz1 = BS%HO_cyl%nz(n1+i0)
                    nr1 = BS%HO_cyl%nr(n1+i0)
                    ml1 = BS%HO_cyl%ml(n1+i0)
                    if (ml1 .ne. ml2) cycle
                    il2 = 2 - n2/n1 ! if n2==n1, il2 = 1, otherwise il2=2 (symetric contribution,so twice for n1!=n2)
                    do it=1,itx !loop over neutrons and protons
                        s = zero
                        sc = zero
                        k1 = dirac%ka(ib,it) + 1 ! ka(k,itx): begining index of eneries in k block
                        k2 = dirac%ka(ib,it) + dirac%kd(ib,it) ! ending index of energies in k block; kd(k,itx): dimension of eneries in k block                       
                        do k =k1,k2 ! loop over different eigenstates of Dirac equation k
                            s  = s  + dirac%fg(nfg+n1,k,it)*dirac%fg(nfg+n2,k,it)*pairing%vv(k,it)
                            sc = sc + dirac%fg(nfg+n1,k,it)*dirac%fg(nfg+n2,k,it)*pairing%skk(k,it) &
                                      *sqrt(pairing%vv(k,it)*(2d0-pairing%vv(k,it)))
                        enddo
                        rshell(it) = s*il2
                        if(ifg.eq.1) then
                            cshell(it) = sc*il2
                        else
                            cshell(it) = 0.0d0
                        endif
                    enddo
                    !loop over the mesh-points
                    ben = az2*(nz1+nz2+1) + ap2*2*(nr1+nr2+ml1+1)
                    do il =1,ngl
                        qlab = BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql(nr2,ml1,il)
                        qltab = ap2*(BS%HO_cyl%ql1(nr1,ml1,il)*BS%HO_cyl%ql1(nr2,ml1,il) + ml1**2*qlab)/gauss%xl(il)
                        bfl = ap2*gauss%xl(il) - ben
                        do ih =1,ngh
                            qhab = BS%HO_cyl%qh(nz1,ih)*BS%HO_cyl%qh(nz2,ih)                            
                            sro = qlab*qhab
                            ihl = ih + (il-1)*ngh
                            do it=1,itx
                                ! scalar and vector density
                                fgr = rshell(it)*sro
                                rs(ihl,it) = rs(ihl,it) - ivv*fgr ! scalar
                                rv(ihl,it) = rv(ihl,it) + fgr ! vector
                                rkapp(ihl,it) = rkapp(ihl,it) + cshell(it)*sro ! kappa
                                ! delta scalar and vector density
                                qh1ab = az2*BS%HO_cyl%qh1(nz1,ih)*BS%HO_cyl%qh1(nz2,ih)
                                stau  = qh1ab*qlab + qhab*qltab                      
                                sdro  = sro*(az2*gauss%xh(ih)**2 + bfl)
                                sdt   = 2*rshell(it)*(sdro+stau)     
                                drs(ihl,it) = drs(ihl,it) - ivv*sdt ! Delta scalar 
                                drv(ihl,it) = drv(ihl,it) + sdt ! Delta vector
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    if(ifPrint) call printDensityResult
    do i =1,nghl
        il = (i-1)/ngh + 1
        ih = mod(i-1,ngh)+1 
        s = 1.d0/gauss%wdcor(ih,il) ! remove the weight factors(wh and wl, which are included in qh and ql) of Gaussian points ! multiply 1/(b_z*(b_\perp)^2\pi)
        densities%rs  = rs
        densities%rv  = rv
        densities%drs = drs
        densities%drv = drv
        densities%ro(i,1)  = s*(rs(i,1) + rs(i,2))  ! rho_s
        densities%ro(i,2)  = s*(rv(i,1) + rv(i,2))  ! rho_v
        densities%ro(i,3)  = s*(-rs(i,1) + rs(i,2)) ! rho_TS ?
        densities%ro(i,4)  = s*(-rv(i,1) + rv(i,2)) ! rho_3 /rho_TV
        densities%dro(i,1) = s*(drs(i,1) + drs(i,2)) 
        densities%dro(i,2) = s*(drv(i,1) + drv(i,2)) 
        densities%dro(i,3) = s*(-drs(i,1) + drs(i,2))  
        densities%dro(i,4) = s*(-drv(i,1) + drv(i,2))
        densities%rkapp(i,1)=s*rkapp(i,1) ! kappa
        densities%rkapp(i,2)=s*rkapp(i,2)
        densities%drvp(i)   =drv(i,2)
    enddo
    contains
    subroutine printDensityResult
        use Globals,only: outputfile
        integer :: it,i
        real(r64) :: s
        ! check, whether integral over drv vanishes
        write(outputfile%u_outputf,*) "******************BEGIN DENSITY********************************"
        do it=1,itx
            s = zero
            do i =1,nghl
                s = s+ drv(i,it)
            enddo
            if(abs(s).gt.1.0d-5) then
                write(outputfile%u_outputf, *) 'drv is not zero in densit',it,s
            endif
            write(outputfile%u_outputf,'(a,i3,2f15.8)') 'Integral over dro:',it,s
        enddo

        ! check the particle number
        do it=1,itx
            s  = zero                                                      
            do i = 1,nghl                                                  
                s  =  s + rv(i,it)                                          
            enddo                                                          
            write(outputfile%u_outputf,'(a,i3,2f15.8)') 'Integral over rv :',it,s
        enddo
        write(outputfile%u_outputf,'(a,3i6)') 'The gauss mesh point: x, y, tot', ngh, ngl, nghl
        write(outputfile%u_outputf,*) "******************END DENSITY************************************"
    end subroutine printDensityResult
end subroutine calculate_densities_DIR

subroutine calculate_densities_RHB(ifPrint)
    use Constants, only: zero,one,nghl,itx,ngl,ngh,NBX
    use Globals, only: BS,gfv,gauss
    logical,intent(in),optional :: ifPrint
    integer :: i,it,ib,nh,ifg,nd,i0,nfg,ivv,&
               n2,nz2,nr2,ml2,n1,nz1,nr1,ml1,&
               il2,il,ih,ihl
    real(r64) :: ap2,az2,s,ben,qlab,qltab,bfl,qhab,sro,fgr,qh1ab,stau,sdro,sdt
    real(r64),dimension(2) :: rshell
    real(r64), dimension(nghl,2) :: rs  ! \rho_s
    real(r64), dimension(nghl,2) :: drs ! \Delta \rho_s
    real(r64), dimension(nghl,2) :: rv  ! \rho_v
    real(r64), dimension(nghl,2) :: drv ! \Delta \rho_v
    
    !! ------calculate single-particle density and the pairing tensor-----------
    call calculate_density_and_kappa_RHB

    !! ------calculate densities: rho_s, rho_V, rho_TS, rho_TV ----------------
    ap2 = one/(BS%HO_cyl%b0 * BS%HO_cyl%bp)**2
    az2 = one/(BS%HO_cyl%b0 * BS%HO_cyl%bz)**2
    ! initialize to zero
    do i=1,nghl
        do it =1,2
            rs(i,it)  = zero
            drs(i,it) = zero
            rv(i,it)  = zero
            drv(i,it) = zero
            ! 

        enddo
    enddo

    do ib=1,BS%HO_cyl%nb ! loop over K-parity-blocks
        nh = BS%HO_cyl%id(ib,1) + BS%HO_cyl%id(ib,2)
        do ifg = 1,2 ! loop over contributions form large and small components
            nd = BS%HO_cyl%id(ib,ifg)
            i0 = BS%HO_cyl%ia(ib,ifg)
            nfg = (ifg-1)*BS%HO_cyl%id(ib,1)
            ivv = gfv%iv(ifg) ! iv(n) = (-1)**n
            do n2 = 1,nd
                nz2 = BS%HO_cyl%nz(n2+i0)
                nr2 = BS%HO_cyl%nr(n2+i0)
                ml2 = BS%HO_cyl%ml(n2+i0)
                do n1 = n2,nd ! due to the contribution is symmetric and this loop only cover half 
                    nz1 = BS%HO_cyl%nz(n1+i0)
                    nr1 = BS%HO_cyl%nr(n1+i0)
                    ml1 = BS%HO_cyl%ml(n1+i0)
                    if (ml1 .ne. ml2) cycle
                    il2 = 2 - n2/n1 ! if n2==n1, il2 = 1, otherwise il2=2 (symetric contribution,so twice for n1!=n2)
                    !
                    rshell(1) = densities%rosh(nfg+n1 + (nfg+n2-1)*nh, ib)*il2
                    rshell(2) = densities%rosh(nfg+n1 + (nfg+n2-1)*nh, ib+NBX)*il2
                    !loop over the mesh-points
                    ben = az2*(nz1+nz2+1) + ap2*2*(nr1+nr2+ml1+1)
                    do il =1,ngl
                        qlab = BS%HO_cyl%ql(nr1,ml1,il)*BS%HO_cyl%ql(nr2,ml1,il)
                        qltab = ap2*(BS%HO_cyl%ql1(nr1,ml1,il)*BS%HO_cyl%ql1(nr2,ml1,il) + ml1**2*qlab)/gauss%xl(il)
                        bfl = ap2*gauss%xl(il) - ben
                        do ih =1,ngh
                            qhab = BS%HO_cyl%qh(nz1,ih)*BS%HO_cyl%qh(nz2,ih)                            
                            sro = qlab*qhab
                            ihl = ih + (il-1)*ngh
                            do it=1,itx
                                ! scalar and vector density
                                fgr = rshell(it)*sro
                                rs(ihl,it) = rs(ihl,it) - ivv*fgr ! scalar
                                rv(ihl,it) = rv(ihl,it) + fgr ! vector
                                !

                                ! delta scalar and vector density
                                qh1ab = az2*BS%HO_cyl%qh1(nz1,ih)*BS%HO_cyl%qh1(nz2,ih)
                                stau  = qh1ab*qlab + qhab*qltab                      
                                sdro  = sro*(az2*gauss%xh(ih)**2 + bfl)
                                sdt   = 2*rshell(it)*(sdro+stau)     
                                drs(ihl,it) = drs(ihl,it) - ivv*sdt ! Delta scalar 
                                drv(ihl,it) = drv(ihl,it) + sdt ! Delta vector
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    if(ifPrint) call printDensityResult
    do i =1,nghl
        il = (i-1)/ngh + 1
        ih = mod(i-1,ngh)+1 
        s = 1.d0/gauss%wdcor(ih,il) ! remove the weight factors(wh and wl, which are included in qh and ql) of Gaussian points ! multiply 1/(b_z*(b_\perp)^2\pi)
        densities%rs  = rs
        densities%rv  = rv
        densities%drs = drs
        densities%drv = drv
        densities%ro(i,1)  = s*(rs(i,1) + rs(i,2))  ! rho_s
        densities%ro(i,2)  = s*(rv(i,1) + rv(i,2))  ! rho_v
        densities%ro(i,3)  = s*(-rs(i,1) + rs(i,2)) ! rho_TS ?
        densities%ro(i,4)  = s*(-rv(i,1) + rv(i,2)) ! rho_3 /rho_TV
        densities%dro(i,1) = s*(drs(i,1) + drs(i,2)) 
        densities%dro(i,2) = s*(drv(i,1) + drv(i,2)) 
        densities%dro(i,3) = s*(-drs(i,1) + drs(i,2))  
        densities%dro(i,4) = s*(-drv(i,1) + drv(i,2))
        !
        densities%drvp(i)   =drv(i,2)
    enddo
    contains
    subroutine printDensityResult
        use Globals,only: outputfile
        integer :: it,i
        real(r64) :: s
        ! check, whether integral over drv vanishes
        write(outputfile%u_outputf,*) "******************BEGIN DENSITY********************************"
        do it=1,itx
            s = zero
            do i =1,nghl
                s = s+ drv(i,it)
            enddo
            if(abs(s).gt.1.0d-5) then
                write(outputfile%u_outputf, *) 'drv is not zero in densit',it,s
            endif
            write(outputfile%u_outputf,'(a,i3,2f15.8)') 'Integral over dro:',it,s
        enddo

        ! check the particle number
        do it=1,itx
            s  = zero                                                      
            do i = 1,nghl                                                  
                s  =  s + rv(i,it)                                          
            enddo                                                          
            write(outputfile%u_outputf,'(a,i3,2f15.8)') 'Integral over rv :',it,s
        enddo
        write(outputfile%u_outputf,*) "******************END DENSITY************************************"
    end subroutine printDensityResult
end subroutine calculate_densities_RHB

subroutine calculate_density_and_kappa_RHB
    use Constants, only: zero,NBX,NHHX,NB2X
    use Globals, only: pairing,BS,RHB,RHB_pairing
    real(r64) :: sp,sr,sk
    integer :: it,il,ib,nf,ng,nh,k1,ke,k1a,kea,mul,m,n2,n1,k,i0,ml1,ml2,i12
    allocate(densities%rosh(NHHX,NB2X) )
    do it = 1,2
        sp  = zero
        il = 0
        do ib = 1,BS%HO_cyl%nb ! loop over the blocks
            nf  = BS%HO_cyl%id(ib,1)
            ng  = BS%HO_cyl%id(ib,2)
            nh  = nf + ng
            k1  = RHB%ka(ib,it) + 1
            ke  = RHB%ka(ib,it) + RHB%kd(ib,it)
            k1a = RHB%ka(ib,it+2) + 1
            kea = RHB%ka(ib,it+2) + RHB%kd(ib,it+2)
            mul = BS%HO_cyl%mb(ib)
            m   = ib + (it-1)*NBX
            ! calculation of rho_{n,n'}
            do n2 = 1,nh
                do n1 =  n2,nh
                    sr = zero
                    do k = k1,ke
                        sr = sr + RHB%fguv(nh+n1,k,it)*RHB%fguv(nh+n2,k,it)
                    enddo
                    do k = k1a,kea ! no-sea approximation
                        sr = sr + RHB%fguv(nh+n1,k,it+2)*RHB%fguv(nh+n2,k,it+2)
                    enddo
                    sr = mul*sr
                    densities%rosh(n1+(n2-1)*nh,m) = sr
                    densities%rosh(n2+(n1-1)*nh,m) = sr
                enddo
            enddo
            ! contributions of large components f*f to kappa
            i0  = BS%HO_cyl%ia(ib,1)
            do n2 = 1,nf
                do n1 =  n2,nf
                    ml1 = BS%HO_cyl%ml(i0+n1)
                    ml2 = BS%HO_cyl%ml(i0+n2)
                    i12 = 2 - n2/n1
                    il  = il + 1
                    sk = zero
                    do k = k1,ke
                        sk = sk + RHB%fguv(nh+n1,k,it)*RHB%fguv(n2,k,it)
                    enddo
                    do k = k1a,kea ! no-sea approximation
                        sk = sk + RHB%fguv(nh+n1,k,it+2)*RHB%fguv(n2,k,it+2)
                    enddo
                    sk = mul*sk
                    if (ml1.ne.ml2) sk = zero  ! remove question ???
                    RHB_pairing%kappa(il,it) = i12*sk ! pairing tensor: kappa
                    if (n1.eq.n2) sp = sp + RHB_pairing%kappa(i12,it) ! trace of kappa
                enddo
            enddo
            pairing%spk(it) = sp
        enddo
    enddo
    return
end subroutine calculate_density_and_kappa_RHB

END MODULE Density