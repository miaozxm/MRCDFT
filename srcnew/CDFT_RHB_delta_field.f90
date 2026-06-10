!==============================================================================!
! MODULE DeltaField                                                               !
!                                                                              !
! This module calculate the pairing delta field                                !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!

Module DeltaField
use Constants, only: r64,NB2X,NHHX,NNNX,MVX,NBX,N0FX
use Globals, only : RHB_pairing,BS
implicit none
integer :: nnmax ! total number of (N_z , N_p)
integer :: mv ! upper triangular elements number of pairing matrix(all blocks and only large components)
contains

subroutine set_separable_pairing_parameters
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------

    ! set mv
    call set_variable_mv

    ! set separabel force's two parameters: a and G
    RHB_pairing%gl(1) = -728.d0 ! -G of neutron, Mev fm^3
    RHB_pairing%gl(2) = -728.d0 ! -G of proton, Mev fm^3
    RHB_pairing%ga    = 0.415d0 ! a^2 ( a = 0.644204936 fm)

end subroutine set_separable_pairing_parameters

subroutine set_variable_mv
    integer :: ib,nf,n2,n1
    ! count mv
    mv = 0
    do ib = 1,BS%HO_cyl%nb ! all blocks
         nf = BS%HO_cyl%id(ib,1) ! only large components
         do n2 =  1,nf
            do n1 = n2,nf ! only triangular
                mv = mv + 1
            enddo
         enddo
    enddo
    return
end subroutine set_variable_mv

subroutine initial_delta_field
    !--------------------------------------------------------------------------
    !   intialze pairing field \Delta as pairing%del 
    !   Note: only \Delta_{ff} set as pairing%del
    !--------------------------------------------------------------------------
    use Constants, only: zero,one
    use Globals, only: pairing
    integer :: i,j,ib,nf,ng,nh,m1,m2,it,ne,n
    real(r64) :: delt
    ! set zeros
    do i =1,NB2X
        do j=1,NHHX
            RHB_pairing%delta(j,i) = 0.d0
        enddo
    enddo
    ! initialize delta field
    do ib = 1,BS%HO_cyl%nb
        nf  = BS%HO_cyl%id(ib,1)
        ng  = BS%HO_cyl%id(ib,2)
        nh = nf + ng
        m1 = ib
        m2 = ib + NBX
        do it = 1,2
            if (abs(pairing%del(it)).lt.1.d-5) then
                pairing%del(it)=1.d-5
                pairing%spk(it) = zero
            else
                pairing%spk(it) = one
            endif         
        enddo

        ne = nf ! only \Delta_{ff} 
        delt = pairing%del(1)
        do n = 1,ne
            RHB_pairing%delta(n+(n-1)*nh,m1) = delt
        enddo

        ne = nf ! ! only \Delta_{ff} 
        delt = pairing%del(2)
        do n = 1,ne
            RHB_pairing%delta(n+(n-1)*nh,m2) = delt
        enddo
    enddo
    return
end subroutine initial_delta_field

subroutine calculate_pairing_matrix_element_W
    !-------------------------------------------------------------------------
    !     calculate paring matrix element W_{12}^{N_zN_p}
    ! where W_{12}^{N_zN_p} = V_{12}^{N_z} V_{12}^{N_p}
    !-------------------------------------------------------------------------
    use Constants, only: zero
    integer :: nnzm,nnrm,nn,nnz,nnr,i
    real(r64),dimension(MVX) :: wn
    real(r64),dimension(1:MVX, 0:N0FX) :: vnz, vnr
    real(r64) :: smax, eps = 1.d-20
    ! calculate V_{12}^{N_z} and V_{12}^{N_p}
    call calculate_pairing_matrix_element_VNz(nnzm,vnz) ! calculate V_{12}^{N_z}
    call calculate_pairing_matrix_element_VNp(nnrm,vnr) ! calculate V_{12}^{N_p}
    ! calculate the single particle matrix elements W_{12}^{N_zN_p}
    nn = 0
    do nnz = 0,nnzm ! N_z
        do nnr = 0,nnrm ! N_p
            nn = nn + 1
            if (nn.gt.NNNX) stop '[delta_field]: NNNX too small'
            smax = zero
            do i = 1,mv
                wn(i) = vnz(i,nnz)*vnr(i,nnr) ! V_{12}^{N_z} * V_{12}^{N_p}
                smax = max(smax,abs(wn(i)))
            enddo
            if (smax.lt.eps) then ! too small
                nn = nn - 1
            else
                do i = 1,mv
                    RHB_pairing%wnn(i,nn) = wn(i) ! W_{12}^{N_zN_p}
                enddo
            endif 
        enddo
    enddo
    nnmax = nn
end subroutine calculate_pairing_matrix_element_W

subroutine calculate_pairing_matrix_element_VNz(nnzm,vnz)
    !-------------------------------------------------------------------------
    !   calculate V_{12}^{N_z}
    !-------------------------------------------------------------------------
    real(r64),dimension(1:MVX, 0:N0FX),intent(out) :: vnz
    integer, intent(out) :: nnzm
    integer :: k,i,il,ib,nf,i0f,n2,nz2,n1,nz1,nn,n,nh,nnh,nx,n0f
    real(r64) :: b0z
    real(r64),dimension(0:N0FX) :: pnosc
    n0f = BS%HO_cyl%n0f
    ! initialize to zero
    do k = 0,N0FX
        do i = 1,mv
            vnz(i,k) = 0.d0
        enddo
    enddo
    b0z  = BS%HO_cyl%b0 * BS%HO_cyl%bz
    call pnoscz(n0f,sqrt(RHB_pairing%ga),b0z,pnosc)
    il = 0
    do ib = 1,BS%HO_cyl%nb
        nf = BS%HO_cyl%id(ib,1)
        i0f = BS%HO_cyl%ia(ib,1)
        do n2 = 1,nf
            nz2 = BS%HO_cyl%nz(i0f+n2)
            do n1 = n2,nf ! only triangular
                nz1 = BS%HO_cyl%nz(i0f+n1)
                il  = il + 1                   
                if (mod(nz1+nz2,2).eq.0) then 
                    do nn = 0,nz1+nz2,2 ! N_z
                        n  = nz1 + nz2 - nn ! n_z = n_{z_1} + n_{z_2} - N_z
                        nh  = n/2
                        nnh = nn/2
                        if (nh .gt.n0f) stop '[calculate_pairing_matrix_element_VNz]: n too large'
                        if (nnh.gt.n0f) stop '[calculate_pairing_matrix_element_VNz]: nn too large'
                        vnz(il,nnh) = pnosc(nh)*Talmi_Moshinsky_1D(nz1,nz2,nn,n) ! V_{12}^{N_z}
                    enddo
                endif
            enddo
        enddo
    enddo
    ! calculate the maximal nnz:  nnzm 
    nnzm =0
    do nn = 0,n0f
        nx = 0
        do i = 1,mv
            if (abs(vnz(i,nn)).gt.1.d-20) then
               nx = nn
            endif
        enddo
        nnzm = max(nnzm,nx)
    enddo
end subroutine calculate_pairing_matrix_element_VNz

subroutine pnoscz(nm,a,b,pnosc)
    use Constants, only: pi
    use Globals, only: gfv
    integer,intent(in) :: nm
    real(r64),intent(in) :: a,b
    real(r64),dimension(0:N0FX),intent(out) :: pnosc
    real(r64) :: s0,s1
    integer :: nh,n
    s0 = (2*pi)**0.25d0
    s0 = s0*sqrt(b/(a*a+b*b))
    s1 = (a*a-b*b)/(a*a+b*b)
    do nh = 0,nm
        n = 2*nh
        pnosc(nh) = s0 * s1**nh * gfv%wf(n) / 2**nh * gfv%fi(nh)
    enddo
end subroutine pnoscz

function Talmi_Moshinsky_1D(n1,n2,n3,n4)
    !-----------------------------------------------------------------------------------
    !   1-dimensional Talmi-Moshinsky bracket: <n1 , n2 | n3 , n4 >
    ! quantum number n_i start from zero: n = 0,1,2,....
    ! 
    ! M_{N_z,n_z}^{n_{z_1}n_{z_2}}
    !     =\frac{1}{\sqrt{2^{N_{z}+n_{z}}}}\sqrt{\frac{n_{z_{1}}!n_{z_{2}}!}{N_{z}!n_{z}!}}
    !     \delta_{n_{z_{1}}+n_{z_{2}},N_{z}+n_{z}}\sum_{s=0}^{n_{z}}(-)^{s}
    !    \left(\begin{array}{c}N_{z}\\n_{z_{1}}-n_{z}+s\end{array}\right.\right)
    !    \left(\begin{array}{c}n_{z}\\s\end{array}\right)
    ! where n1=n_{z_1}, n2= n_{z_2}, n3=N_z, n4=n_z
    ! ----------------------------------------------------------------------------------
    !   wf(n)   = sqrt(n!) 
    !   wfi(n)  = 1/sqrt(n!)
    !   iv(n)   = (-1)**n
    !   ibc(m,n)  = m!/(n!(m-n)!) 
    !-----------------------------------------------------------------------------------
    use Globals, only: gfv
    use Constants, only: zero, two
    double precision ::  Talmi_Moshinsky_1D
    integer, intent(in) :: n1,n2,n3,n4
    real(r64) :: f,s
    integer :: m4,m3
    Talmi_Moshinsky_1D = zero
    if (n1+n2.ne.n3+n4) return
    f = gfv%wf(n1) * gfv%wfi(n3) * gfv%wf(n2) * gfv%wfi(n4) / sqrt(two**(n3+n4))
    s = zero
    do m4 = 0,n4
        m3 = n1 - n4 + m4
        if (m3.ge.0.and.n3.ge.m3) then 
            s = s + gfv%iv(m4) * gfv%ibc(n3,m3) * gfv%ibc(n4,m4)
        endif
    enddo
    Talmi_Moshinsky_1D = f*s
    return
end function Talmi_Moshinsky_1D

subroutine calculate_pairing_matrix_element_VNp(nnrm,vnr)
    !-------------------------------------------------------------------------
    !   calculate V_{12}^{N_p}
    !-------------------------------------------------------------------------
    real(r64),dimension(1:MVX,0:N0FX),intent(out) :: vnr
    integer, intent(out) :: nnrm
    integer :: k,i,il,ib,nf,i0f,n2,nr2,ml2,n1,nr1,ml1,nn,n,nx,n0f
    real(r64) :: b0p
    real(r64),dimension(0:N0FX) :: pnosc
    n0f = BS%HO_cyl%n0f
    ! initialize to zero
    do k = 0, N0FX
        do i= 1, mv
            vnr(i,k) = 0.d0
        enddo
    enddo
    b0p  = BS%HO_cyl%b0 * BS%HO_cyl%bp
    call pnoscp(n0f,sqrt(RHB_pairing%ga),b0p,pnosc)

    il = 0
    do ib = 1,BS%HO_cyl%nb
        nf = BS%HO_cyl%id(ib,1)
        i0f= BS%HO_cyl%ia(ib,1)
        do n2 = 1,nf
            nr2 = BS%HO_cyl%nr(i0f+n2)
            ml2 = BS%HO_cyl%ml(i0f+n2)
            do n1 = n2,nf
                nr1 = BS%HO_cyl%nr(i0f+n1)
                ml1 = BS%HO_cyl%ml(i0f+n1)
                il  = il + 1        
                if (ml1.eq.ml2) then
                    do nn = 0,nr1+nr2+ml1 ! N_p
                        n = nr1+nr2+ml1-nn ! n_r = n_{r_1} + n_{r_2} + |m_{l_1}| - N_p
                        if (nn.gt.n0f) stop '[calculate_pairing_matrix_element_VNp]: nn too large'
                        if (n .gt.n0f) stop '[calculate_pairing_matrix_element_VNp]: n  too large'
                        vnr(il,nn) = pnosc(n)*Talmi_Moshinsky_2D(ml1,nr1,nr2,0,nn,n) ! V_{12}^{N_p}
                    enddo   ! nn
                endif
            enddo
        enddo
    enddo
    ! calculate the maximal nnr:  nnrm 
    nnrm = 0
    do nn = 0,n0f
        nx = 0
        do i = 1,mv
            if (abs(vnr(i,nn)).gt.1.d-20) then
                nx = nn
            endif
        enddo
        nnrm = max(nnrm,nx)
    enddo
end subroutine calculate_pairing_matrix_element_VNp

subroutine pnoscp(nm,a,b,pnosc)
    use Constants, only: pi
    integer, intent(in) :: nm
    real(r64),intent(in) :: a,b
    real(r64),dimension(0:nm) ::  pnosc
    integer :: n
    real(r64) :: f1,f2
    if(nm.gt.BS%HO_cyl%n0f) stop '[PNOSCP]: nm too large'
    do n = 0,nm
        f1=b/(b**2+a**2)
        f2=((b**2-a**2)/(b**2+a**2))**n
        pnosc(n) = f1*f2/2/pi
    enddo
    return
end subroutine pnoscp

function Talmi_Moshinsky_2D(im1,in1,in2,im3,in3,in4)
    !-------------------------------------------------------------------------------------------------------------------------------------
    !   2-dimensional-Talmi-Moshinsky bracket: <n1 m1, n2 -m1 | n3 m3, n4 -m3>
    ! radial quantum number start from zero: n=0,1,2,....	
    ! orbital angular momentum m1=-m2>0,m3=-m4>0
    !
    !  \begin{aligned} 
    !    & M_{N_p M_p n_p m_p}^{n_{r_1} m_{l_1} n_{r_2} m_{l_2}}\\
    !    &=\frac{(-)^{N_p+n_p-n_{r_1}-n_{r_2}}}{\sqrt{2^{2 N_p+2 n_p+\left|M_p\right|+\left|m_p\right|}}} 
    !       \sqrt{\frac{\left(n_{r_1}\right) !\left(n_{r_1}+\left|m_{l_1}\right|\right)!
    !                   \left(n_{r_2}\right) !\left(n_{r_2}+\left|m_{l_2}\right|\right)!}
    !           {\left(N_p\right) !\left(N_p+\left|M_p\right|\right) !\left(n_p\right)!
    !            \left(n_p+\left|m_p\right|\right) !}} \\ 
    !    & \times \delta_{2 n_{r_1}+\left|m_{l_1}\right|+2 n_{r_2}+\left|m_{l_2}\right|, 2 N_p+\left|M_p\right|+2 n_p+\left|m_p\right|} 
    !             \delta_{m_{l_1}+m_{l_2}, M_p+m_p} \\ 
    !    & \times \sum_{Q, R, S=0}^{N_p} \sum_{T=0}^{M_p} \sum_{q, r, s=0}^{n_p} \sum_{t=0}^{m_p}(-)^{r+s+t}
    !             \left(\begin{array}{cccc}N_p & & \\ N_p-Q-R-S & Q & R & S\end{array}\right)
    !             \left(\begin{array}{c}M_p \\ T\end{array}\right)
    !             \left(\begin{array}{cccc}n_p & & \\ n_p-q-r-s & q & r & s\end{array}\right)
    !             \left(\begin{array}{c}m_p \\ t\end{array}\right) \\ 
    !  \end{aligned}
    ! where n1=n_{r_1}, m1=m_{l_1}, n2=n_{r_2}, m2=m_{l_2}, n3=N_p, m3=M_p, n4=n_p, m4=m_p
    !-----------------------------------------------------------------------------------------------------------------------------------
    !     iv(n)  =  (-1)**n
    !     fak(n) =  n!
    !     fi(n)  =  1/n!
    !     wf(n)  =  sqrt(n!)
    !     wfi(n) =  1/sqrt(n!)
    !-------------------------------------------------------------------------
    use Constants, only: zero,two
    use Globals, only: gfv
    double precision ::  Talmi_Moshinsky_2D
    integer, intent(in) :: im1,in1,in2,im3,in3,in4
    integer :: im2, im4, nn12,nn34, n1,n2,n3,n4,m1,m2,m3,m4, &
                nn1,nn2,nn3,nn4,i3,j3,k3,l3,i4,j4,k4,l4,it3,it4
    real(r64) :: s34, prout,sn3,sn4,sm3,sm4

    im2 = -im1
    im4 = -im3
    if (im1.lt.0)     stop '[Talmi_Moshinsky_2D]: m1 < 0'
    if (im3.lt.0)     stop '[Talmi_Moshinsky_2D]: m3 < 0'

    Talmi_Moshinsky_2D = zero
    nn12 = 2*in1+iabs(im1)+2*in2+iabs(im2)
    nn34 = 2*in3+iabs(im3)+2*in4+iabs(im4)
    if (nn12.ne.nn34) return
    s34 = zero
    if (im1.lt.im2) then
        n1 = in2
        n2 = in1
        m1 = im2
        m2 = im1
    else
        n1 = in1
        n2 = in2
        m1 = im1
        m2 = im2
    endif
    n3 = in3
    n4 = in4
    m3 = im3
    m4 = im4

    nn1 = 2*n1 + abs(m1)
    nn2 = 2*n2 + abs(m2)
    nn3 = 2*n3 + abs(m3)
    nn4 = 2*n4 + abs(m4)
    if (m3.gt.m4) then
        if (n1+n2.ne.n3+n4+m2-m4) return
    else
        if (n1+n2.ne.n3+n4+m2-m3) return
    endif

    prout = gfv%iv(n3+n4-n1-n2) / sqrt(two**(nn3+nn4))* &
            gfv%wf(n1)*gfv%wf(n1+abs(m1))* &
            gfv%wf(n2)*gfv%wf(n2+abs(m2))* &
            gfv%wfi(n3)*gfv%wfi(n3+abs(m3))* &
            gfv%wfi(n4)*gfv%wfi(n4+abs(m4))

    sn3 = zero     
    sn4 = zero  
    sm3 = zero
    sm4 = zero
    do i3 = 0,n3
    do j3 = 0,n3
    do k3 = 0,n3
        l3 = n3 - i3 - j3 - k3
        do i4 = 0,n4
        do j4 = 0,n4
        do 20 k4 = 0,n4
            l4 = n4 - i4 - j4 - k4
            if (l3.lt.0.or.l4.lt.0) goto 20 
            do it3 = 0,abs(m3)
                do 10 it4 = 0,abs(m4)
                    if (m3.gt.m4) then
                        if (i3+i4+j3+j4+it3.ne.n2) goto 10
                        if (j3+j4.ne.m2+k3+k4-it3+it4) goto 10
                        if (l3+l4.ne.n3+n4-n2-k3-k4+it3) goto 10
                    else
                        if (i3+i4+j3+j4+it4.ne.n2) goto 10
                        if (j3+j4.ne.m2+k3+k4+it3-it4) goto 10
                        if (l3+l4.ne.n3+n4-n2-k3-k4+it4) goto 10
                    endif
                    sn3 = gfv%fak(n3) * gfv%fi(l3) * gfv%fi(i3) * gfv%fi(j3) * gfv%fi(k3)
                    sn4 = gfv%iv(j4+k4) * gfv%fak(n4) * gfv%fi(l4) * gfv%fi(i4) * gfv%fi(j4) * gfv%fi(k4)
                    sm3 = gfv%fak(abs(m3)) * gfv%fi(it3) * gfv%fi(abs(m3)-it3)
                    sm4 = gfv%iv(it4) * gfv%fak(abs(m4)) * gfv%fi(it4) * gfv%fi(abs(m4)-it4)
                    s34 = s34 + sn3*sn4*sm3*sm4
                10 continue
            enddo
        20 continue
        enddo
        enddo
    enddo
    enddo
    enddo
    Talmi_Moshinsky_2D = s34*prout
    return
end function Talmi_Moshinsky_2D

subroutine calculate_delta_field
    !-------------------------------------------------------------------------
    !    calculate pairing field \Delta  in a separable form of the Gogny force.
    ! where \Delta_{12} = -G\sum_{N_z}sum_{N_p} W_{12}^{N_zN_p*} P_{N_zN_p} 
    ! and P_{N_zN_p} = 1/2\sum_{12}W_{12}^{N_zN_p}\kappa_{12}
    !-------------------------------------------------------------------------
    use Constants, only: half,zero
    integer :: it,i,i12,ib,i0,nf,ng,nh,m,n2,n1,nn
    real(r64) :: s,g
    real(r64), dimension(NNNX) :: pnn ! P_{N_zN_p}
    do it = 1 ,2
        ! calculate P_{N_zN_p}
        do nn = 1,nnmax
            s = zero
            do i = 1,mv
                s = s + RHB_pairing%wnn(i,nn) * RHB_pairing%kappa(i,it)
            enddo
            pnn(nn) = s
        enddo
        ! calculate \Delta_{12}
        g = half * RHB_pairing%gl(it) ! 1/2*G
        i12 = 0
        do ib = 1,BS%HO_cyl%nb
            i0 = BS%HO_cyl%ia(ib,1)
            nf = BS%HO_cyl%id(ib,1)
            ng = BS%HO_cyl%id(ib,2)
            nh = nf + ng
            m  = ib + (it-1)*NBX
            do n2 =  1,nf
                do n1 = n2,nf
                    i12 = i12 + 1
                    s = zero
                    do nn = 1,nnmax
                        s = s + RHB_pairing%wnn(i12,nn)*pnn(nn)
                    enddo 
                    RHB_pairing%delta(n1+(n2-1)*nh,m) = -g*s
                    RHB_pairing%delta(n2+(n1-1)*nh,m) = -g*s
                enddo
            enddo
        enddo
    enddo
end subroutine calculate_delta_field

END Module DeltaField