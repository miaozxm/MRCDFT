!==============================================================================!
! MODULE Dirac                                                                 !
!                                                                              !
! This module solves the Dirac-Equation in cylindrical oscillator basis        !                                                 !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE DiracEquation
use Constants, only: r64
use Globals, only: dirac,BS
implicit none

contains

subroutine solve_dirac_equation(ifPrint)
    !------------------------------------------------------------------------------!
    !                                                                              !
    ! solves the Dirac-Equation in cylindrical oscillator basis                    !
    ! units: Mev                                                                   !
    !                                                                              !
    !------------------------------------------------------------------------------!
    use MathMethods, only: sdiag
    use Constants, only: itx,nhx,nkx
    use Globals, only: nucleus_attributes
    logical,intent(in),optional :: ifPrint
    integer :: it,ib,kl,nf,ng,nh,k1,k,i
    real(r64), dimension(nhx*nhx) :: hh
    real(r64), dimension(nhx) :: e, ez
    character(len=8), dimension(nhx) :: tfg
    logical :: endContinue
    do it = 1,itx !loop over protons and neutrons
        kl = 0
        do ib = 1,BS%HO_cyl%nb
            dirac%ka(ib,it) = kl
            ! calculate dirac matrix of ib block
            call calculate_block_dirac_matrix(it,ib,hh,nf,ng,endContinue)
            nh = nf + ng
            if(endContinue) cycle
            if(.False.) call printDiracMatrix
            !! diagonalization
            ! call the subroutine with the parameter '+1,' it makes the energies in ascending order
            call sdiag(nh,nh,hh,e,hh,ez,+1)
            !! store energies and wavefunctions
            ! no sea approximation: only nf states (positive energy states) are considered
            k1 = kl
            do k = 1,nf
               kl = kl + 1
               dirac%nblo(kl,it) = ib
               dirac%ee(kl,it)   = e(ng+k)  !consider only positive energy
               do i = 1,nf+ng
                  dirac%fg(i,kl,it) = hh(i+(ng+k-1)*(nf+ng))
               enddo
            enddo
            dirac%kd(ib,it) = nf
            if(ifPrint .and. .True.) call printEnergyAndCoefficientOfEachBlock                                                
        enddo
        dirac%nk(it) = kl
        if(kl>nkx) stop '[DiracEquation]: nkx too small'
        if (it==1 .and. 2*kl .lt. nucleus_attributes%neutron_number) stop '[DiracEquation]: pwi too small'
        if (it==2 .and. 2*kl .lt. nucleus_attributes%proton_number) stop '[DiracEquation]: pwi too small'
    enddo

    contains
    subroutine printDiracMatrix
        use Globals, only: outputfile
        integer :: i0f,i0g
        i0f = BS%HO_cyl%ia(ib,1)   ! begin of the large components of block b is ia(b,1)+1 
        i0g = BS%HO_cyl%ia(ib,2)   ! begin of the small components of block b is ia(b,2)+1
        do i = 1,nf
            tfg(i) = BS%HO_cyl%tb(i+i0f)
        enddo
        do i = 1,ng
            tfg(i+nf) = BS%HO_cyl%tb(i+i0g)
        enddo
        write(outputfile%u_outputf,'(/,a)') BS%HO_cyl%txb(ib)
        call aprint(outputfile%u_outputf,2,2,6,nh,nh,nh,hh,tfg,(/('',i=1,nh)/),'HH')
    end subroutine printDiracMatrix
    subroutine printEnergyAndCoefficientOfEachBlock
        use Globals, only: outputfile
        call aprint(outputfile%u_outputf,1,1,6,1,1,nf,dirac%ee(k1+1,it),(/('',i=1,1)/),(/('',i=1,nf)/),'E')
        call aprint(outputfile%u_outputf,1,2,6,nh,nh,nf,dirac%fg(1,k1+1,it),tfg,(/('',i=1,nf)/),'FG')
    end subroutine printEnergyAndCoefficientOfEachBlock
end subroutine solve_dirac_equation

subroutine calculate_all_dirac_matrix
    use Constants, only: NBX
    integer :: it, ib, m ,nf, ng
    logical :: endContinue
    ! calculate initial Hamiltonian matrix element for RHB equation
    do it = 1,2
        do ib = 1,BS%HO_cyl%nb
            m = ib + (it-1)*NBX
            call calculate_block_dirac_matrix(it,ib,dirac%hh(1,m),nf,ng,endContinue)
            ! if(endContinue) cycle
        enddo
    enddo
end subroutine calculate_all_dirac_matrix

subroutine calculate_block_dirac_matrix(it,ib,hh,nf,ng,endContinue)
    !-----------------------------------------------------------------
    !   calculate the dirac matrix of block ib
    !----------------------------------------------------------------
    use Constants, only: nhx,hbc
    use Globals, only: force,matrices,fields
    integer :: it,ib,nf,ng
    logical :: endContinue
    real(r64), dimension(nhx*nhx) :: hh
    integer :: nh,i0f,i0g,i2,i1,i
    real(r64) :: emcc2
    endContinue = .False.
    emcc2 = 2*force%masses%amu*hbc                             
    nf  = BS%HO_cyl%id(ib,1)   ! dimension large components of block b                                            
    if (nf.eq.0) then
        endContinue = .True. 
        return                                     
    endif
    ng  = BS%HO_cyl%id(ib,2)   ! dimension small components of block b                                           
    nh  = nf + ng                                               
    i0f = BS%HO_cyl%ia(ib,1)   ! begin of the large components of block b is ia(b,1)+1                                            
    i0g = BS%HO_cyl%ia(ib,2)   ! begin of the small components of block b is ia(b,2)+1
    
    !! calculation of the Dirac-matrix
    ! the hamitonian is a 2*2 matrix as defined in the text
    !         (A, B )  i.e. (ff, fg)
    !         (B, -C)       (gf, gg)

    ! calculate the section (gf)
    do i2 = 1,nf                                                
        do i1 = 1,ng                                             
            hh(nf+i1+(i2-1)*nh) = matrices%sp(i1+(i2-1)*ng,ib) ! B already been computed by subroutine calculate_sigma_nabla in Matrix module 
        enddo                                                    
    enddo
    ! calculate the section v+s (ff)
    call pot(i0f, nf, nh, fields%vpstot(1,1,it), hh)                 
    ! calculate the section v-s (gg)
    call pot(i0g, ng, nh, fields%vmstot(1,1,it), hh(nf+1+nf*nh)) 
    ! the matrix has been shifted by amu therefore the original v+s+amu
    ! is now v+s, while the orignal v-s-amu is now (v-s)-2*amu
    do i = nf+1,nh                                              
        hh(i+(i-1)*nh) = hh(i+(i-1)*nh) - emcc2                  
    enddo
    ! the Dirac Matrix is a symmetric matrix, 
    ! set (fg) by making the matrix hh symmetric
    do i1=1,nh
        do i2=i1,nh
            hh(i1+(i2-1)*nh)=hh(i2+(i1-1)*nh)
        enddo
    enddo
end subroutine calculate_block_dirac_matrix

subroutine pot(i0,n,nh,v,tt)
    !----------------------------------------------------------------------
    !                                                                       
    !     integrate 
    !           \int  \Psi_{\alpha} \Psi_{\alpha'} f(\boldsymbol{r})  d\boldsymbol{r} 
    !     in the oscillator basis:
    !           \Psi_{\alpha} = \Phi_{n_z}(z)\Phi_{n_r}^{m_l}(r_{\perp})\frac{1}{\sqrt{2\pi}}e^{im_l\varphi} \chi_{ms}(s)
    ! 
    !     So, the integral expression is: 
    !       \delta_{m_s m_s'} \delta_{m_lm_l'}  \int \Phi_{n_z'}(z) \Phi_{n_r'}^{m_l'}(r_{\perp})  f(\boldsymbol{r}) \Phi_{n_z}(z) \Phi_{n_r}^{m_l}(r_{\perp}) d\boldsymbol{r}
    !       =  \delta_{m_s m_s'} \delta_{m_lm_l'} \sum_{z_i} \sum_{r_j}   \Phi_{n_z'}(z_i)\Phi_{n_z}(z_i)wh(z_i) *\Phi_{n_r'}^{m_l'}(r_j)\Phi_{n_r}^{m_l}(r_j)wl(r_j)*f(z_i,r_j)
    !       =  \delta_{m_s m_s'} \delta_{m_lm_l'} \sum_{z_i} \sum_{r_j}   qh(n_z',z_i)*qh(n_z,z_i)*ql(n_r',r_j)*ql(n_r,r_j)*f(z_j,r_j)
    !    where, z_i, r_j is the gauss points.
    !    If \Psi_{\alpha} and \Psi_{\alpha'} in the same block, m_s==m_s' is automatically satisfied while m_l==m_l', because  K = m_l + m_s 
    !------------------------------------------------------------------------------   
    !                                                                       
    !     i0   initial point for the quantum number array                   
    !     n    dimension of \alpha
    !     nh   dimension of array tt
    !     v    the function f() in the above integral                                                   
    !     tt   big array to store the integral results                                                                
    !----------------------------------------------------------------------
    !NOTE:
    !  1. "For each loop, i1 starts from i2,  so only the upper triangular matrix elements are computed."
    !  2. In calculate_Cylindrical_HO_basis_In_Gauss_points we allocate qh (and qh1,ql...) as qh (0:nz_max, 1:nh)), 
    !     the first index start from 0, however in this subroutine we start i2 from 1 because i0 which from ia(b,1)
    !     (or ia(b,2)) and ia(b,1)+1 is the begin of the large(small) components of block b. 
    !----------------------------------------------------------------
    use Constants, only: ngh,ngl
    integer, intent(in) :: i0,n,nh
    real(r64), intent(in),dimension(ngh,ngl) :: v
    real(r64), dimension(nh,nh) :: tt

    integer :: nz2,nr2,ml2,nz1,nr1,ml1,i2,i1,il,ih
    real(r64) :: t, s
    do i2 = 1,n                                                      
        nz2 = BS%HO_cyl%nz(i0+i2)                                                
        nr2 = BS%HO_cyl%nr(i0+i2)                                                
        ml2 = BS%HO_cyl%ml(i0+i2)                                                
        do i1 = i2,n                                                   
        nz1 = BS%HO_cyl%nz(i0+i1)                                             
        nr1 = BS%HO_cyl%nr(i0+i1)                                             
        ml1 = BS%HO_cyl%ml(i0+i1)                                                                                                                    
        if (ml1.eq.ml2) then                                        
            t = 0.0d0                                                
            do il = 1,ngl                                                                                                                
                s = 0.0d0                                              
                do ih = 1,ngh                                         
                    s = s + v(ih,il) * BS%HO_cyl%qh(nz1,ih) * BS%HO_cyl%qh(nz2,ih)             
                enddo
                t = t + s*BS%HO_cyl%ql(nr1,ml1,il) * BS%HO_cyl%ql(nr2,ml2,il)               
            enddo                                                    
            tt(i1,i2) = t        
        else                                                        
            tt(i1,i2) = 0.0d0     
        endif
        enddo                                                          
    enddo
end subroutine pot

subroutine calculate_sigma_nabla(ifPrint)
    !-----------------------------------------------------------------------------------------------------
    ! calculates single particle matrix elements:
    ! \int d\boldsymbol{r} \Psi_{\alpha} \boldsymbol{\sigma}\cdot\boldsymbol{\nabla} \Psi_{\tilde{\alpha}}
    ! for Fermions in the cylindrical oscillator basis
    !-----------------------------------------------------------------------------------------------------
    use Constants, only: hbc,zero
    use Globals, only: gfv,gauss,matrices
    logical,intent(in),optional :: ifPrint
    integer :: ib,nf,ng,i0f,i0g,i2,i1, nz2,nr2,ml2,ms2,nz1,nr1,ml1,ms1, m,il
    real(r64) :: fz,fp,s
    fz = hbc/(BS%HO_cyl%b0*BS%HO_cyl%bz)
    fp = hbc/(BS%HO_cyl%b0*BS%HO_cyl%bp)
    do ib =1,BS%HO_cyl%nb
        nf = BS%HO_cyl%id(ib,1) ! dimension of large components in block ib
        ng = BS%HO_cyl%id(ib,2) ! dimension of small components in block ib
        i0f = BS%HO_cyl%ia(ib,1) ! begin of large components in block ib
        i0g = BS%HO_cyl%ia(ib,2) ! begin of small components in block ib
        do i2 = 1,nf
            nz2 = BS%HO_cyl%nz(i0f+i2)
            nr2 = BS%HO_cyl%nr(i0f+i2)
            ml2 = BS%HO_cyl%ml(i0f+i2)
            ms2 = BS%HO_cyl%ms(i0f+i2)
            do i1 = 1,ng
                nz1 = BS%HO_cyl%nz(i0g+i1)
                nr1 = BS%HO_cyl%nr(i0g+i1)
                ml1 = BS%HO_cyl%ml(i0g+i1)
                ms1 = BS%HO_cyl%ms(i0g+i1)
                
                s = zero
                if(ms1.eq.ms2) then ! also ml1 equal ml2 with the same k
                  if (nr1.eq.nr2) then                                     
                    if (nz1.eq.nz2+1) s = -ms1*fz*gfv%sq(nz1)*gfv%sqi(2)          
                    if (nz1.eq.nz2-1) s = +ms1*fz*gfv%sq(nz2)*gfv%sqi(2)
                  endif
                else
                    if (nz1.eq.nz2) then
                        if (ml1.eq.ml2+1) then                                
                            m = -ml2                                           
                        else                                                  
                            m = +ml2                                           
                        endif                                                 
                        do il = 1,gauss%nl                                         
                            s=s+BS%HO_cyl%ql(nr1,ml1,il)*(BS%HO_cyl%ql1(nr2,ml2,il)+m*BS%HO_cyl%ql(nr2,ml2,il))/sqrt(gauss%xl(il))
                        enddo                                                
                        s = s*fp
                    endif
                endif
                matrices%sp(i1+(i2-1)*ng,ib) = -s 
            enddo
        enddo
    enddo

    if(ifPrint) call printSIgmaNablaMatrices
    contains

    subroutine printSIgmaNablaMatrices
        use Constants, only: nfx,ngx
        use Globals, only: outputfile
        integer :: ib,i,nf,ng,i0f,i0g
        character(len=8) :: tp(nfx),tm(ngx)
        write(outputfile%u_config,outputfile%format_title) 'calculate_sigma_nabla'
        do ib = 1,BS%HO_cyl%nb
            nf = BS%HO_cyl%id(ib,1) ! dimension of large components in block ib
            ng = BS%HO_cyl%id(ib,2) ! dimension of small components in block ib
            i0f = BS%HO_cyl%ia(ib,1) ! begin of large components in block ib
            i0g = BS%HO_cyl%ia(ib,2) ! begin of small components in block ib
            write(outputfile%u_config,"(/,a)") BS%HO_cyl%txb(ib)
            do i= 1,nf
                tp(i) = BS%HO_cyl%tb(i+i0f)
            enddo
            do i = 1,ng
                tm(i) = BS%HO_cyl%tb(i+i0g)
            enddo
            call aprint(outputfile%u_config,1,3,6,ng,ng,nf,matrices%sp(1,ib),tm,tp,'Sigma * P')
        enddo
    end subroutine printSIgmaNablaMatrices
end subroutine calculate_sigma_nabla

subroutine aprint(u_write,is,it,ns,ma,n1,n2,a,t1,t2,text)
    !----------------------------------------------------------------!
    !                                                                !
    !     IS = 1    Full matrix                                      !
    !          2    Lower diagonal matrix                            !
    !          3    specially stored symmetric matrix                !
    !                                                                !
    !     IT = 1    numbers for rows and columns                     !
    !          2    text for rows and numbers for columns            !
    !          3    text for rows and columns                        !
    !                                                                !
    !     NS = 1     FORMAT   8F8.4     80 Coulums                   !
    !     NS = 2     FORMAT   8f8.2     80 Coulums                   !
    !     NS = 3     FORMAT  17F4.1     80 Coulums                   !
    !     NS = 4     FORMAT  30F4..1    120 Coulums                  !
    !     NS = 5     FORMAT  5F12.8     80 Coulums                   !
    !     NS = 6     FORMAT  5F12.4     80 Coulums                   !
    !     NS = 7     FORMAT  4E13.6     80 Coulums                   !
    !     NS = 8     FORMAT  8E15.8    130 Coulums                   !
    !----------------------------------------------------------------!
    implicit double precision (a-h,o-z)
    implicit integer (i-n)
    integer :: u_write,is,it,ns,ma,n1,n2
    real(r64),dimension(ma*n2),intent(in) :: a
    character(len=8) :: t1(n1), t2(n2)
    character(len=*) :: text

    character(len=30) :: fmt1, fmt2
    character*20 :: fti,ftt,fmt(8),fmti(8),fmtt(8)
    integer:: nsp(8),nspalt

    data nsp/8,8,17,30,5,5,4,8/
    data fmt /'8f8.4)',            '8F8.2)',   &         
              '17f4.1)',           '30f4.1)',  &         
              '5f12.8)',           '5f12.4)',  &        
              '4e13.6)',           '8e15.8)'/           
    data fmti/'(11x,8(i4,4x))',    '(11x,8(i4,4x))',    &
              '(11x,17(1x,i2,1x))','(11x,30(1x,i2,1x))',&
              '(11x,6(i4,8x))',    '(11x,10(i4,8x))',   &
              '(11x,5(i4,9x))',    '(11x,8(i4,11x))'/   
    data fmtt/'(11x,8a8)',         '(11x,8a8)',       &  
              '(11x,17a4)',        '(11x,30a4)',      &  
              '(11x,6(a8,2x))',    '(11x,10(4x,a8))', &  
              '(11x,5(a8,5x))',    '(11x,8(a8,7x))'/

    fmt1   = '(4x,i3,4x,' // fmt(ns)
    fmt2   = '(1x,a8,2x' // fmt(ns) 
    fti    = fmti(ns)               
    ftt    = fmtt(ns)               
    nspalt = nsp(ns)
    write(u_write,'(//,3x,a)') text
    ka = 1                           
    ke = nspalt  ! end coulum number to display in first Part                      
    nteil = n2/nspalt 
    if (nteil*nspalt.ne.n2) nteil = nteil + 1
    do 10  nt = 1,nteil                                            
        if (n2.gt.nspalt)  write(u_write,100)  nt                       
 100    format(//, 10x,'Part',i5,' of the Matrix',/)                
        if (nt.eq.nteil) ke = n2    ! in last Part, coulum end to n2                                  
        if (it.lt.3) then                                           
            write(u_write,fti) (k,k=ka,ke) ! number for coulums                               
        else                                                        
            write(u_write,ftt) (t2(k),k=ka,ke)  ! text for coulums                          
        endif  
                                                    
        do 20  i=1,n1                                                  
            kee=ke                                                      
            if (is.eq.2.and.ke.gt.i) kee=i                              
            if (ka.gt.kee) goto 20                                      
            if (is.eq.3) then                                           
                if (it.eq.1) then                                        
                    write(u_write,fmt1) i,(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
                else                                                     
                    write(u_write,fmt2) t1(i),(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
                endif                                                    
            else                                                        
                if (it.eq.1) then                                        
                    write(u_write,fmt1) i,(a(i+(k-1)*ma),k=ka,kee)             
                else                                                     
                    write(u_write,fmt2) t1(i),(a(i+(k-1)*ma),k=ka,kee)         
                endif                                                    
            endif                                                       
    20  continue                                                       

    ka=ka+nspalt                                                   
    ke=ke+nspalt                                                   
 10 continue
end subroutine aprint

END MODULE DiracEquation