!==============================================================================!
! MODULE MathMethods                                                           !
!                                                                              !
! This module contains the variables and routines related to mathematical fun- !
! ctions and numerical methods.                                                !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE MathMethods
    ! BLAS external
    implicit none
    external :: zgemv,zdotu,zdotc,zGEMM

contains

subroutine blas_set_num_threads_local(n)
    integer :: n 
    call mkl_set_num_threads_local(n)
end subroutine

subroutine print_blas_num_threads
    integer :: mkl_get_max_threads
    print *, "Global target threads:", mkl_get_max_threads()
end subroutine

function zTrace(A, N) result(tr)
    implicit none
    complex*16, intent(in) :: A(:, :)
    integer, intent(in) :: N
    complex :: tr
    integer :: i, max_dim
    max_dim = min(N, size(A, 1), size(A, 2))
    tr = 0.0
    do i = 1, max_dim
        tr = tr + A(i, i)
    end do
end function


subroutine zGEMM_Trace(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    implicit none
    character(len=1), intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex*16, intent(in) :: alpha, beta
    complex*16, intent(in) :: A(lda, *)
    complex*16, intent(in) :: B(ldb, *)
    complex*16, intent(inout) :: C(ldc, *)

    integer :: i, l
    complex*16 :: temp
    
    !$omp parallel do private(i,temp,l)
    do i = 1, min(m,n)

        if ((transa=='T' .or. transa=='t') .and. &
            (transb=='N' .or. transb=='n')) then
            ! C_ii = sum_l A(l,i) * B(l,i)
            call zdotu(temp, k, A(1,i), 1, B(1,i), 1)

        else if ((transa=='C' .or. transa=='c') .and. &
                 (transb=='N' .or. transb=='n')) then
            ! C_ii = sum_l conj(A(l,i)) * B(l,i)
            call zdotc(temp, k, A(1,i), 1, B(1,i), 1)

        else
            temp = (0.0d0,0.0d0)
            do l = 1, k
                select case (transa)
                case ('N','n')
                    temp = temp + A(i,l) * select_b(transb, B, i, l)
                case ('T','t')
                    temp = temp + A(l,i) * select_b(transb, B, i, l)
                case ('C','c')
                    temp = temp + conjg(A(l,i)) * select_b(transb, B, i, l)
                end select
            end do
        end if
        C(i,i) = alpha * temp + beta * C(i,i)
    end do
    !$omp end parallel do

    contains

    pure function select_b(transb, B, i, l) result(val)
        character(len=1), intent(in) :: transb
        complex*16, intent(in) :: B(ldb,*)
        integer, intent(in) :: i, l
        complex*16 :: val

        select case (transb)
        case ('N','n')
            val = B(l,i)
        case ('T','t')
            val = B(i,l)
        case ('C','c')
            val = conjg(B(i,l))
        end select
    end function

end subroutine

subroutine sort(n, e, descending)
    !-----------------------------------------------------------------------
    !
    !     Sorts a set of numbers according to their size
    !
    !     Input:
    !       n          - number of elements to sort
    !       e          - array to be sorted
    !       descending - logical flag:
    !                    .true.  = sort in descending order (large to small)
    !                    .false. = sort in ascending order (small to large)
    !-----------------------------------------------------------------------
    integer, intent(in) :: n
    double precision, intent(inout) :: e(n)
    logical, intent(in) :: descending

    integer :: i, j, k
    double precision :: p

    do i = 1, n-1
        k = i
        p = e(i)
        do j = i+1, n
            if (descending) then
                ! Look for larger elements for descending sort
                if (e(j) > p) then
                    k = j
                    p = e(j)
                end if
            else
                ! Look for smaller elements for ascending sort
                if (e(j) < p) then
                    k = j
                    p = e(j)
                end if
            end if
        end do
        if (k /= i) then
            e(k) = e(i)
            e(i) = p
        end if
    end do
end subroutine sort

subroutine GaussLaguerre(N,X,W)
!-------------------------------------------------------------------------------------!
!Purpose : Compute the zeros of Laguerre polynomial Ln(x)in the interval [0,inf],     !
!          and the corresponding weighting coefficients for Gauss-Laguerr integration.!
!                                                                                     !                                                                                                  
!Input   : n    --- Order of the Laguerre polynomial                                  !
!          X(n) --- Zeros of the Laguerre polynomial                                  !
!          W(n) --- Corresponding weighting coefficients                              !
!                                                                                     ! 
!Integral: WL  =  W * exp(X)                                                          !
!          \int_0^\infty  f(z) exp(-z) dz  =   \sum_i f(X(i)) W(i)                    !
!          \int_0^\infty  f(z) dz          =   \sum_i f(X(i)) WL(i)                    !
!-------------------------------------------------------------------------------------!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    DIMENSION X(N),W(N)
    HN=1.0D0/N
    DO 35 NR=1,N
        IF (NR.EQ.1) Z=HN
        IF (NR.GT.1) Z=X(NR-1)+HN*NR**1.27
        IT=0
10      IT=IT+1
        Z0=Z
        P=1.0D0
        DO 15 I=1,NR-1
            P=P*(Z-X(I))
15      CONTINUE
        F0=1.0D0
        F1=1.0D0-Z
        DO 20 K=2,N
            PF=((2.0D0*K-1.0D0-Z)*F1-(K-1.0D0)*F0)/K
            PD=K/Z*(PF-F1)
            F0=F1
            F1=PF
20        CONTINUE
        FD=PF/P
        Q=0.0D0
        DO 30 I=1,NR-1
            WP=1.0D0
            DO 25 J=1,NR-1
                IF (J.EQ.I) GO TO 25
                WP=WP*(Z-X(J))
25          CONTINUE
            Q=Q+WP
30      CONTINUE
        GD=(PD-Q*FD)/P
        Z=Z-FD/GD
        IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
        X(NR)=Z
        W(NR)=1.0D0/(Z*PD*PD)           
35  CONTINUE
    RETURN
end subroutine GaussLaguerre

subroutine GaussHermite(N,X,W)
!-------------------------------------------------------------------------------------!
!Purpose : Compute the zeros of Hermite polynomial Ln(x) in the interval [-inf,inf],  ! 
!          and the corresponding weighting coefficients for Gauss-Hermite integration.!
!                                                                                     !
!Input   : n    --- Order of the Hermite polynomial                                   !  
!          X(n) --- Zeros of the Hermite polynomial                                   !
!          W(n) --- Corresponding weighting coefficients                              !
!                                                                                     !
!Integral: WH  =  W * exp(X**2)                                                       !
!          \int_-infty^\infty  f(z) exp(-z**2) dz  =   \sum_i f(X(i)) W(i)            !
!          \int_-infty^\infty  f(z) dz             =   \sum_i f(X(i)) WH(i)           !
!-------------------------------------------------------------------------------------!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    DIMENSION X(N),W(N)
    if(N.ge.80) stop "[GaussHermite] N is too large, and the generated grid points are incorrect!"
    HN=1.0D0/N
    ZL=-1.1611D0+1.46D0*N**0.5
    DO 40 NR=1,N/2
        IF (NR.EQ.1) Z=ZL
        IF (NR.NE.1) Z=Z-HN*(N/2+1-NR)
        IT=0
10      IT=IT+1
        Z0=Z
        F0=1.0D0
        F1=2.0D0*Z
        DO 15 K=2,N
            HF=2.0D0*Z*F1-2.0D0*(K-1.0D0)*F0
            HD=2.0D0*K*F1
            F0=F1
            F1=HF
15      CONTINUE
        P=1.0D0
        DO 20 I=1,NR-1
            P=P*(Z-X(I))
20      CONTINUE
        FD=HF/P
        Q=0.0D0
        DO 30 I=1,NR-1
            WP=1.0D0
            DO 25 J=1,NR-1
                IF (J.EQ.I) GO TO 25
                WP=WP*(Z-X(J))
25          CONTINUE
            Q=Q+WP
30      CONTINUE
        GD=(HD-Q*FD)/P
        Z=Z-FD/GD
        IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
        X(NR)=Z
        X(N+1-NR)=-Z
        R=1.0D0
        DO 35 K=1,N
            R=2.0D0*R*K
35      CONTINUE
        W(NR)=3.544907701811D0*R/(HD*HD)
        W(N+1-NR)=W(NR)
40  CONTINUE
    IF (N.NE.2*INT(N/2)) THEN
        R1=1.0D0
        R2=1.0D0
        DO 45 J=1,N
            R1=2.0D0*R1*J
            IF (J.GE.(N+1)/2) R2=R2*J
45      CONTINUE
        W(N/2+1)=0.88622692545276D0*R1/(R2*R2)
        X(N/2+1)=0.0D0
    ENDIF
    RETURN
end  subroutine GaussHermite

subroutine GaussLegendre(N,X,W)
!-------------------------------------------------------------------------------------!
!Purpose : Compute the zeros of Legendre polynomial Pn(x) in the interval [-1,1], and ! 
!          the corresponding weighting coefficients for Gauss-Legendre integration.   !
!                                                                                     !
!Input   : n    --- Order of the Legendre polynomial                                  !  
!          X(n) --- Zeros of the Legendre polynomial                                  !
!          W(n) --- Corresponding weighting coefficients                              !
!                                                                                     !
!Integral:                                                                            !
!          \int_{-1}^{+1}     f(z) dz  =   \sum_i  f(X(i)) * W(i)                     !
!-------------------------------------------------------------------------------------!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    DIMENSION X(N),W(N)
    N0=(N+1)/2
    DO 45 NR=1,N0
        Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10      Z0=Z
        P=1.0D0
        DO 15 I=1,NR-1
            P=P*(Z-X(I))
15      CONTINUE
        F0=1.0D0
        IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
        F1=Z
        DO 20 K=2,N
            PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
            PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
            F0=F1
            F1=PF
20      CONTINUE
        IF (Z.EQ.0.0) GO TO 40
        FD=PF/P
        Q=0.0D0
        DO 35 I=1,NR
            WP=1.0D0
            DO 30 J=1,NR
                IF (J.NE.I) WP=WP*(Z-X(J))
30          CONTINUE
            Q=Q+WP
35      CONTINUE
        GD=(PD-Q*FD)/P
        Z=Z-FD/GD
        IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40      X(NR)=Z
        X(N+1-NR)=-Z
        W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
        W(N+1-NR)=W(NR)
45  CONTINUE
    RETURN
end subroutine GaussLegendre

subroutine GaussLegendre_X1toX2(X1,X2,X,W,N)
!-------------------------------------------------------------------------------------!
!Purpose : Compute the zeros of Legendre polynomial Pn(x) in the interval [X1,X2], and! 
!          the corresponding weighting coefficients for Gauss-Legendre integration.   !
!                                                                                     !
!Input   : n    --- Order of the Legendre polynomial                                  !  
!          X(n) --- Zeros of the Legendre polynomial                                  !
!          W(n) --- Corresponding weighting coefficients                              !
!                                                                                     !
!Integral:                                                                            !
!          \int_{-1}^{+1}     f(z) dz  =   \sum_i  f(X(i)) * W(i)                     !
!-------------------------------------------------------------------------------------!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    dimension w(n),x(n) 
    PARAMETER (EPS=3.d-14)
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do 12 i=1,m
    z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
        p1=1.d0
        p2=0.d0
        do 11 j=1,n
        p3=p2
        p2=p1
        p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
    if(abs(z-z1).gt.EPS)goto 1
    x(i)=xm-xl*z
    x(n+1-i)=xm+xl*z
    w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
    w(n+1-i)=w(i)
12    continue
    return
end subroutine GaussLegendre_X1toX2

subroutine math_gfv
!-------------------------------------------------------------------------------------!
!Purpose : Calculates sign, sqrt, factorials, etc. of integers and half int.          !
!                                                                                     !
!Input   :                                  !  
!                                                                                     !
!Variable Meaning:                                                                    !
!          iv(n)  =  (-1)**n                                                          !
!          sq(n)  =  sqrt(n)                                                          ! 
!          sqi(n) =  1/sqrt(n)                                                        !
!          sqh(n) =  sqrt(n+1/2)                                                      !
!          shi(n) =  1/sqrt(n+1/2)                                                    !
!          fak(n) =  n!                                                               !
!          ibc(m,n) = m!/(n!(m-n)!)                                                   !
!          fad(n) =  (2*n+1)!!                                                        !
!          fdi(n) =  1/(2*n+1)!!                                                      !
!          fi(n)  =  1/n!                                                             !
!          wf(n)  =  sqrt(n!)                                                         !
!          wfi(n) =  1/sqrt(n!)                                                       !
!          wfd(n) =  sqrt((2*n+1)!!)                                                  !
!          gm2(n) =  gamma(n+1/2)                                                     !
!          gmi(n) =  1/gamma(n+1/2)                                                   !
!          wg(n)  =  sqrt(gamma(n+1/2))                                               !
!          wgi(n) =  1/sqrt(gamma(n+1/2))                                             !
!-------------------------------------------------------------------------------------!
    use Constants, only: igfv,igfvbc,zero,one,half,pi
    use Globals, only: gfv
    integer :: i,m,n
    gfv%iv(0)  = +1
    gfv%sq(0)  = zero
    gfv%sqi(0) = 1.d30
    gfv%sqh(0) = sqrt(half)
    gfv%shi(0) = 1/gfv%sqh(0)
    gfv%fak(0) = one
    gfv%fad(0) = one
    gfv%fi(0)  = one
    gfv%fdi(0) = one
    gfv%wf(0)  = one
    gfv%wfi(0) = one
    gfv%wfd(0) =  one
    !gm2(0) = Gamma(1/2) = sqrt(pi)
    gfv%gm2(0) =  sqrt(pi)
    gfv%gmi(0) =  1/gfv%gm2(0)
    gfv%wg(0)  =  sqrt(gfv%gm2(0))
    gfv%wgi(0) =  1/gfv%wg(0)
    do i = 1,igfv
            gfv%iv(i)  = -gfv%iv(i-1)
            gfv%iv(-i) = gfv%iv(i)
            gfv%sq(i)  = dsqrt(dfloat(i))
            gfv%sqi(i) = one/gfv%sq(i)
            gfv%sqh(i) = sqrt(i+half)
            gfv%shi(i) = one/gfv%sqh(i)
            gfv%fak(i) = i*gfv%fak(i-1)
            gfv%fad(i) = (2*i+1)*gfv%fad(i-1)
            gfv%fi(i)  = one/gfv%fak(i)
            gfv%fdi(i) = one/gfv%fad(i)
            gfv%wf(i)  = gfv%sq(i)*gfv%wf(i-1)
            gfv%wfi(i) = one/gfv%wf(i)
            gfv%wfd(i) = sqrt(gfv%fad(i))
            gfv%gm2(i) = (i-half)*gfv%gm2(i-1)
            gfv%gmi(i) = one/gfv%gm2(i)
            gfv%wg(i)  = gfv%sqh(i-1)*gfv%wg(i-1)
            gfv%wgi(i) = one/gfv%wg(i)
    enddo
    !THE ARRAY OF BINOMIAL COEFFICIENTS
    gfv%ibc(0,0)= one
    do m=1,igfvbc
        do n=0,m
            gfv%ibc(m,n)=gfv%fak(m)/(gfv%fak(n)*gfv%fak(m-n)) 
        enddo   
    enddo  
    return
end subroutine math_gfv

! recursive function factorial(n) result(facto)
! !------------------------------------------------------------------------------!
! ! function factorial                                                           !
! !                                                                              !
! ! Computes the factorial: n! = n * (n-1) * ... * 1                             !
! !------------------------------------------------------------------------------!
! integer, intent(in) :: n
! real(r64) :: facto 
! if ( n <= 0 ) then 
!   facto = one  
! else
!   facto = n * factorial(n-1)
! endif
! end function


subroutine sdiag(nmax,n,a,d,x,e,is)                                                                      
!=======================================================================
!       diagonalization of a symmetric matrix   
!                                                                  
!   A   matrix to be diagonalized                                       
!   D   eigenvalues                                                     
!   X   eigenvectors                                                    
!   E   auxiliary field                                                 
!   IS = 1  eigenvalues are ordered and major component of X is positive
!        0  eigenvalues are not ordered                                 
!-----------------------------------------------------------------------
    implicit double precision (a-h,o-z)
    implicit integer(i-n)                                                                                                      
    dimension a(nmax,nmax),x(nmax,nmax),e(n),d(n)                                                                                            
    data tol,eps/1.e-32,1.e-10/,zero/0.d0/,one/1.d0/,half/0.5d0/

    if (n.eq.1) then                                                  
        d(1)   = a(1,1)                                                
        x(1,1) = one                                                   
        return                                                         
    endif                                                                                                                                    
    do i = 1,n                                                        
        do j = 1,i                                                     
        x(i,j)=a(i,j)                                               
        enddo                                                          
    enddo                                                                                                                                   
    ! Householder-reduction                                             
    do i = n,2,-1                                                     
        l = i - 2                                                      
        f = x(i,i-1)                                                   
        g = f                                                          
        h = zero                                                       
        do k = 1,l                                                     
        h = h + x(i,k)*x(i,k)                                       
        enddo                                                          
        s = h + f*f                                                    
        if (s.lt.tol) h = zero                                         
        if (h.gt.zero) then                                            
        l = l+1                                                     
        g = dsqrt(s)                                                
        if (f.ge.zero) g = -g                                       
        h = s - f*g                                                 
        hi = one/h                                                  
        x(i,i-1) = f - g                                            
        f = zero                                                    
        do j = 1,l                                                  
            x(j,i) = x(i,j)*hi                                       
            s = zero                                                 
            do k = 1,j                                               
                s = s + x(j,k)*x(i,k)                                 
            enddo                                                    
            j1 = j+1                                                 
            do k = j1,l                                              
                s = s + x(k,j)*x(i,k)                                 
            enddo                                                    
            e(j) = s*hi                                              
            f = f + s*x(j,i)                                         
        enddo                                                       
        f = f*hi*half                                               
        do j = 1,l                                                  
            s    = x(i,j)                                            
            e(j) = e(j) - f*s                                        
            p    = e(j)                                              
            do k = 1,j                                               
                x(j,k) = x(j,k) - s*e(k) - x(i,k)*p                   
            enddo                                                    
        enddo                                                       
        endif                                                          
        d(i)   = h                                                     
        e(i-1) = g                                                     
    enddo                                                             
                                                                       
    ! transformation matrix                                             
    d(1) = zero                                                       
    e(n) = zero                                                       
    b    = zero                                                       
    f    = zero                                                       
    do i = 1,n                                                        
        l = i-1                                                        
        if (d(i).ne.zero) then                                         
        do j = 1,l                                                  
            s = zero                                                 
            do k = 1,l                                               
                s = s + x(i,k)*x(k,j)                                 
            enddo                                                    
            do k = 1,l                                               
                x(k,j) = x(k,j) - s*x(k,i)                            
            enddo                                                    
        enddo                                                       
        endif                                                          
        d(i)   = x(i,i)                                                
        x(i,i) = one                                                   
        do j = 1,l                                                     
        x(i,j) = zero                                               
        x(j,i) = zero                                               
        enddo                                                          
    enddo                                                             
                                                                       
    ! diagonalizition of tri-diagonal-matrix                            
    do l = 1,n                                                        
        h = eps*( abs(d(l))+ abs(e(l)))                                
        if (h.gt.b) b = h

        ! test for splitting                                             
        do j = l,n                                                     
            if ( abs(e(j)).le.b) goto 10                                
        enddo                                                          
                                                                    
    ! test for convergence                                           
    10  if (j.gt.l) then                                               
    20      p  = (d(l+1)-d(l))/(2*e(l))                                 
            r  = dsqrt(p*p+1.d0)                                        
            pr = p + r                                                  
            if (p.lt.zero) pr = p - r                                   
            h = d(l) - e(l)/pr                                          
            do i = l,n                                                  
                d(i) = d(i) - h                                          
            enddo                                                       
            f = f + h                                                   
                                                                    
            ! QR-transformation                                           
            p = d(j)                                                    
            c = one                                                     
            s = zero                                                    
            do i = j-1,l,-1                                             
                g = c*e(i)                                               
                h = c*p                                                  
                if ( abs(p).lt.abs(e(i))) then                           
                    c = p/e(i)                                            
                    r = dsqrt(c*c+one)                                    
                    e(i+1) = s*e(i)*r                                     
                    s = one/r                                             
                    c = c/r                                               
                else                                                     
                    c = e(i)/p                                            
                    r = dsqrt(c*c+one)                                    
                    e(i+1) = s*p*r                                        
                    s = c/r                                               
                    c = one/r                                             
                endif                                                    
                p = c*d(i) - s*g                                         
                d(i+1) = h + s*(c*g+s*d(i))                              
                do k = 1,n                                               
                    h        = x(k,i+1)                                   
                    x(k,i+1) = x(k,i)*s + h*c                             
                    x(k,i)   = x(k,i)*c - h*s                             
                enddo                                                    
            enddo                                                       
            e(l) = s*p                                                  
            d(l) = c*p                                                  
            if ( abs(e(l)).gt.b) goto 20                                
        endif                                                          
                                                                        
        ! convergence                                                    
        d(l) = d(l) + f                                                
    enddo                                                             
                                                                    
    if (is.eq.0) return                                               
    ! ordering of eigenvalues                                           
    do i = 1,n                                                        
        k  = i                                                         
        p  = d(i)                                                      
        do j = i+1,n                                                   
            if (d(j).lt.p) then                                         
                k = j                                                   
                p = d(j)                                                
            endif                                                       
        enddo                                                          
        if (k.ne.i) then                                               
            d(k) = d(i)                                                 
            d(i) = p                                                    
            do j = 1,n                                                  
                p      = x(j,i)                                          
                x(j,i) = x(j,k)                                          
                x(j,k) = p                                               
            enddo                                                       
        endif                                                          
    enddo                                                             
                                                                    
    ! signum                                                            
    do  k = 1,n                                                       
        s = zero                                                       
        do i = 1,n                                                     
            h = abs(x(i,k))                                             
            if (h.gt.s) then                                            
                s  = h                                                   
                im = i                                                   
            endif                                                       
        enddo                                                          
        if (x(im,k).lt.zero) then                                      
            do i = 1,n                                                  
                x(i,k) = -x(i,k)                                         
            enddo                                                       
        endif                                                          
    enddo                                                             
                                                                    
    return                                                                                                               
end subroutine sdiag                                                               

subroutine brak0(func,x0,x1,x2,f1,f2,step)
!----------------------------------------------------------------------
!     subroutine for bracketing a root of a function
!
!     given a function monotonous groving f(x) = FUNC and 
!     given an initial point X0 this routine searches for an interval
!     such that ther is root of f(x) between x1 and x2
!     x1 < x2
!
!     INPUT:  X0   starting point
!             FUNC monotonously growing function of x
!             STEP inital step-size
!     OUTPUT: X1,X2 an intervall, which brakets a zero of f(x)
!             f1,f2 values of the function f(x) at x1 and x2
!
!----------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
    implicit integer(i-n)
    external func
    data maxit/100/
    x  = x0
    dx = step 
    s  = func(x)
    if (s.lt.0.0) then
        x1 = x
        f1 = s
        do i = 1,maxit
            x = x + dx 
            s = func(x)
            if (s.gt.0.0) then
                x2 = x
                f2 = s
                goto 20
            else
                x1 = x
                f1 = s
            endif
            dx = 2*dx
        enddo
        stop 'BRAK0 no success in bracketing'
    else
        x2 = x
        f2 = s
        do i = 1,maxit
            x = x - dx 
            s = func(x)
            if (s.lt.0.0) then
                x1 = x
                f1 = s
                goto 20
            else 
                x2 = x
                f2 = s
            endif
            dx = 2*dx
        enddo
        stop 'BRAK0 no success in bracketing'
    endif
 20 continue      
    return
end subroutine brak0

function rtbrent(func,x1,x2,y1,y2,tol)      
!-----------------------------------------------------------------------
!     Using Brent's method find the root of the function FUNC(X) 
!     known to lie between X1 and X2.
!     Y1 = FUNC(X1), Y2 = FUNC(X2) 
!     The root will be returned as RTBRENT with accuracy TOL
!
!     from: NUMERICAL RECIPIES, 9.3
!
!----------------------------------------------------------------------
    implicit real*8(a-h,o-z)
    implicit integer(i-n)
    double precision rtbrent

    external func
    parameter (itmax = 100, eps = 1.d-12) !maximum number of iterations, machine floating point precision
    data zero/0.0/,one/1.d0/,half/0.5d0/
    a  = x1
    b  = x2
    fa = y1
    fb = y2
    if (fa*fb.gt.zero) stop 'RTBRENT: root must be bracketed'
    fc = fb
    do 10 iter = 1,itmax
        if (fb*fc.gt.zero) then
            ! rename a,b,c and adjust bounding interval
            c  = a              
            fc = fa
            d  = b - a
            e  = d
        endif
        if (abs(fc).lt.abs(fb)) then
            a  = b
            b  = c
            c  = a
            fa = fb
            fb = fc
            fc = fa
        endif
        ! convergence check
        tol1 = 2*eps*abs(b)+half*tol
        xm = half*(c-b)
        if (abs(xm).le.tol1 .or. fb.eq.zero) then
            rtbrent = b
            return
        endif
        if (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
            ! attempt inverse quadratic interpolation
            s = fb/fa
            if (a.eq.c) then
               p = 2*xm*s
               q = one - s
            else
               q = fa/fc
               r = fb/fc 
               p = s*(2*xm*q*(q-r) - (b-a)*(r-one))
               q = (q-one)*(r-one)*(s-one)
            endif
            if (p.gt.zero) q = -q
            ! check whether in bounds
            p = abs(p)
            if (2*p.lt.dmin1(3*xm*q-abs(tol1*q),abs(e*q))) then
                ! accept interpolation
                e = d
                d = p/q
            else
                ! interpolation failed, use besection
                d = xm
                e = d
            endif
        else
            ! bounds decreasing too slowly, use bisection
            d = xm
            e = d
        endif
        ! move last best guess to a
        a  = b
        fa = fb
        if (abs(d).gt.tol1) then
            ! evaluate new trial root
            b = b + d
        else
            b = b + sign(tol1,xm)
        endif
        fb = func(b)
 10 continue
    stop 'RTBRENT: exceeding maximum number of iterations'
end function rtbrent

function DNRM2(N,X,INCX)
!----------------------------------------------------------------
!  DNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!           DNRM2 := sqrt( x'*x )
!
!  This version written on 25-October-1982.
!  Modified on 14-October-1993 to inline the call to DLASSQ.
!  Sven Hammarling, Nag Ltd.
!----------------------------------------------------------------
    DOUBLE PRECISION dnrm2
    !.. Scalar Arguments ..
    INTEGER INCX,N
    !.. Array Arguments ..
    DOUBLE PRECISION X(*)
    !.. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    !.. Local Scalars ..
    DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
    INTEGER IX
    !.. Intrinsic Functions ..
    INTRINSIC ABS,SQRT

    IF (N.LT.1 .OR. INCX.LT.1) THEN
        NORM = ZERO
    ELSE IF (N.EQ.1) THEN
        NORM = ABS(X(1))
    ELSE
        SCALE = ZERO
        SSQ = ONE
        ! The following loop is equivalent to this call to the LAPACK
        ! auxiliary routine:
        ! CALL DLASSQ( N, X, INCX, SCALE, SSQ )
        DO 10 IX = 1,1 + (N-1)*INCX,INCX
            IF (X(IX).NE.ZERO) THEN
                ABSXI = ABS(X(IX))
                IF (SCALE.LT.ABSXI) THEN
                    SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                    SCALE = ABSXI
                ELSE
                    SSQ = SSQ + (ABSXI/SCALE)**2
                END IF
            END IF
        10 CONTINUE
        NORM = SCALE*SQRT(SSQ)
    END IF
    DNRM2 = NORM
    RETURN
end function DNRM2

subroutine DSCAL(N,DA,DX,INCX)
!----------------------------------------------------------------
!   DSCAL scales a vector by a constant.
!
!   uses unrolled loops for increment equal to one.
!   jack dongarra, linpack, 3/11/78.
!   modified 3/93 to return if incx .le. 0.
!   modified 12/3/93, array(1) declarations changed to array(*)   
!
!----------------------------------------------------------------
    !.. Scalar Arguments ..
    DOUBLE PRECISION DA
    INTEGER INCX,N
    !.. Array Arguments ..
    DOUBLE PRECISION DX(*)
    !.. Local Scalars ..
    INTEGER I,M,MP1,NINCX
    !.. Intrinsic Functions ..
    INTRINSIC MOD
    IF (N.LE.0 .OR. INCX.LE.0) RETURN
    IF (INCX.EQ.1) GO TO 20
    ! code for increment not equal to 1
    NINCX = N*INCX
    DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
    10 CONTINUE
    RETURN
    ! code for increment equal to 1
    20 M = MOD(N,5)
        IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
            DX(I) = DA*DX(I)
    30 CONTINUE
        IF (N.LT.5) RETURN
    40 MP1 = M + 1
        DO 50 I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
    50 CONTINUE
    RETURN
end subroutine DSCAL

function DDOT(N,DX,INCX,DY,INCY)
!----------------------------------------------------------------
!  DDOT forms the dot product of two vectors.
!
!  uses unrolled loops for increments equal to one.
!  jack dongarra, linpack, 3/11/78.
!  modified 12/3/93, array(1) declarations changed to array(*)
!----------------------------------------------------------------
    DOUBLE PRECISION DDOT
    !.. Scalar Arguments ..
    INTEGER INCX,INCY,N
    !.. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
    !.. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY,M,MP1
    !.. Intrinsic Functions ..
    INTRINSIC MOD
    DDOT = 0.0d0
    DTEMP = 0.0d0
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    ! code for unequal increments or equal increments not equal to 1
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
    10 CONTINUE
    DDOT = DTEMP
    RETURN
    ! ccode for both increments equal to 1
    20 M = MOD(N,5)
        IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
            DTEMP = DTEMP + DX(I)*DY(I)
    30 CONTINUE
        IF (N.LT.5) GO TO 60
    40 MP1 = M + 1
        DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) &
                + DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
    50 CONTINUE
    60 DDOT = DTEMP
    RETURN
end function DDOT

subroutine DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!----------------------------------------------------------------
!                -- LAPACK routine (version 3.1) --
!    Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!    November 2006
!
!Purpose:
!    DSYTRF computes the factorization of a real symmetric matrix A using
!    the Bunch-Kaufman diagonal pivoting method.  The form of the
!    factorization is
!   
!       A = U*D*U**T  or  A = L*D*L**T
!   
!    where U (or L) is a product of permutation and unit upper (lower)
!    triangular matrices, and D is symmetric and block diagonal with
!    1-by-1 and 2-by-2 diagonal blocks.
!   
!    This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!Arguments:
!    UPLO    (input) CHARACTER*1
!            = 'U':  Upper triangle of A is stored;
!            = 'L':  Lower triangle of A is stored.
!  
!    N       (input) INTEGER
!            The order of the matrix A.  N >= 0.
!  
!    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!            N-by-N upper triangular part of A contains the upper
!            triangular part of the matrix A, and the strictly lower
!            triangular part of A is not referenced.  If UPLO = 'L', the
!            leading N-by-N lower triangular part of A contains the lower
!            triangular part of the matrix A, and the strictly upper
!            triangular part of A is not referenced.
!  
!            On exit, the block diagonal matrix D and the multipliers used
!            to obtain the factor U or L (see below for further details).
!  
!    LDA     (input) INTEGER
!            The leading dimension of the array A.  LDA >= max(1,N).
!  
!    IPIV    (output) INTEGER array, dimension (N)
!            Details of the interchanges and the block structure of D.
!            If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!            interchanged and D(k,k) is a 1-by-1 diagonal block.
!            If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!            columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!            is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!            IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!  
!    WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!  
!    LWORK   (input) INTEGER
!            The length of WORK.  LWORK >=1.  For best performance
!            LWORK >= N*NB, where NB is the block size returned by ILAENV.
!  
!            If LWORK = -1, then a workspace query is assumed; the routine
!            only calculates the optimal size of the WORK array, returns
!            this value as the first entry of the WORK array, and no error
!            message related to LWORK is issued by XERBLA.
!  
!    INFO    (output) INTEGER
!            = 0:  successful exit
!            < 0:  if INFO = -i, the i-th argument had an illegal value
!            > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
!                  has been completed, but the block diagonal matrix D is
!                  exactly singular, and division by zero will occur if it
!                  is used to solve a system of equations.
!
!Further Details:
!    If UPLO = 'U', then A = U*D*U', where
!       U = P(n)*U(n)* ... *P(k)U(k)* ...,
!    i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!    1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!    and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!    defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!    that if the diagonal block D(k) is of order s (s = 1 or 2), then
!  
!               (   I    v    0   )   k-s
!       U(k) =  (   0    I    0   )   s
!               (   0    0    I   )   n-k
!                  k-s   s   n-k
!  
!    If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!    If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!    and A(k,k), and v overwrites A(1:k-2,k-1:k).
!  
!    If UPLO = 'L', then A = L*D*L', where
!       L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!    i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!    n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!    and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!    defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!    that if the diagonal block D(k) is of order s (s = 1 or 2), then
!  
!               (   I    0     0   )  k-1
!       L(k) =  (   0    I     0   )  s
!               (   0    v     I   )  n-k-s+1
!                  k-1   s  n-k-s+1
!  
!    If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!    If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!    and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!----------------------------------------------------------------
    !.. Scalar Arguments ..
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LWORK, N
    !.. Array Arguments ..
    INTEGER            IPIV( * )
    DOUBLE PRECISION   A( LDA, * ), WORK( * )
    !.. Local Scalars ..
    LOGICAL            LQUERY, UPPER
    INTEGER            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
    !.. External Functions ..
    !LOGICAL            LSAME
    !INTEGER            ILAENV
    !EXTERNAL           LSAME, ILAENV
    !.. External Subroutines ..
    !EXTERNAL           DLASYF, DSYTF2, XERBLA
    !.. Intrinsic Functions ..
    INTRINSIC          MAX

    !.. Executable Statements ..
    ! Test the input parameters.
    INFO = 0
    UPPER = LSAME( UPLO, 'U' )
    LQUERY = ( LWORK.EQ.-1 )
    IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
        INFO = -1
    ELSE IF( N.LT.0 ) THEN
        INFO = -2
    ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
        INFO = -4
    ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
        INFO = -7
    END IF

    IF( INFO.EQ.0 ) THEN
        ! Determine the block size
        NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
        LWKOPT = N*NB
        WORK( 1 ) = LWKOPT
    END IF
    IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DSYTRF', -INFO )
        RETURN
    ELSE IF( LQUERY ) THEN
        RETURN
    END IF
    NBMIN = 2
    LDWORK = N
    IF( NB.GT.1 .AND. NB.LT.N ) THEN
        IWS = LDWORK*NB
        IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DSYTRF', UPLO, N, -1, -1, -1 ) )
        END IF
    ELSE
        IWS = 1
    END IF

    IF( NB.LT.NBMIN ) NB = N

    IF( UPPER ) THEN
        ! Factorize A as U*D*U' using the upper triangle of A
        !
        ! K is the main loop index, decreasing from N to 1 in steps of
        ! KB, where KB is the number of columns factorized by DLASYF;
        ! KB is either NB or NB-1, or K for the last block
        K = N
        10  CONTINUE
        ! If K < 1, exit from loop
        IF( K.LT.1 ) GO TO 40
        IF( K.GT.NB ) THEN
            ! Factorize columns k-kb+1:k of A and use blocked code to
            ! update columns 1:k-kb
            CALL DLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK, IINFO )
        ELSE
            ! Use unblocked code to factorize columns 1:k of A
            CALL DSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
        END IF
        ! Set INFO on the first occurrence of a zero pivot
        IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO
        ! Decrease K and return to the start of the main loop
        K = K - KB
        GO TO 10
    ELSE
        ! Factorize A as L*D*L' using the lower triangle of A
        ! K is the main loop index, increasing from 1 to N in steps of
        ! KB, where KB is the number of columns factorized by DLASYF;
        ! KB is either NB or NB-1, or N-K+1 for the last block
         K = 1
        20  CONTINUE
        ! If K > N, exit from loop
        IF( K.GT.N ) GO TO 40
        IF( K.LE.N-NB ) THEN
            ! Factorize columns k:k+kb-1 of A and use blocked code to
            ! update columns k+kb:n
            CALL DLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ), WORK, LDWORK, IINFO )
        ELSE
            ! Use unblocked code to factorize columns k:n of A
            CALL DSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
        END IF
        ! Set INFO on the first occurrence of a zero pivot
        IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + K - 1
        ! Adjust IPIV
        DO 30 J = K, K + KB - 1
            IF( IPIV( J ).GT.0 ) THEN
                IPIV( J ) = IPIV( J ) + K - 1
            ELSE
                IPIV( J ) = IPIV( J ) - K + 1
            END IF
        30  CONTINUE
        ! Increase K and return to the start of the main loop
        K = K + KB
        GO TO 20
    END IF
    40 CONTINUE
    WORK( 1 ) = LWKOPT
    RETURN
end subroutine DSYTRF

FUNCTION LSAME(CA,CB)
!----------------------------------------------------------------
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!Purpose:
!   LSAME returns .TRUE. if CA is the same letter as CB regardless of case.
!
!Arguments:
!   CA  (input) CHARACTER*1
!   CB  (input) CHARACTER*1
!   CA and CB specify the single characters to be compared.
!----------------------------------------------------------------
    LOGICAL LSAME
    !.. Scalar Arguments ..
    CHARACTER CA,CB
    !.. Intrinsic Functions ..
    INTRINSIC ICHAR
    !.. Local Scalars ..
    INTEGER INTA,INTB,ZCODE
    ! Test if the characters are equal
    LSAME = CA .EQ. CB
    IF (LSAME) RETURN
    ! Now test for equivalence if both characters are alphabetic.
    ZCODE = ICHAR('Z')
    ! Use 'Z' rather than 'A' so that ASCII can be detected on Prime
    ! machines, on which ICHAR returns a value with bit 8 set.
    ! ICHAR('A') on Prime machines returns 193 which is the same as
    ! ICHAR('A') on an EBCDIC machine.
    INTA = ICHAR(CA)
    INTB = ICHAR(CB)
    IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
        ! ASCII is assumed - ZCODE is the ASCII code of either lower or
        ! upper case 'Z'.
        IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
        IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
    ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
        ! EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
        ! upper case 'Z'.
        IF (INTA.GE.129 .AND. INTA.LE.137 .OR. &
            INTA.GE.145 .AND. INTA.LE.153 .OR. &
            INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
        IF (INTB.GE.129 .AND. INTB.LE.137 .OR. &
            INTB.GE.145 .AND. INTB.LE.153 .OR. &
            INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
    ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
        ! ASCII is assumed, on Prime machines - ZCODE is the ASCII code
        ! plus 128 of either lower or upper case 'Z'.
        IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
        IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
    END IF
    LSAME = INTA .EQ. INTB
    RETURN
END FUNCTION LSAME

subroutine XERBLA(SRNAME,INFO)
!----------------------------------------------------------------
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!Purpose:
!   XERBLA  is an error handler for the LAPACK routines.
!   It is called by an LAPACK routine if an input parameter has an
!   invalid value.  A message is printed and execution stops.
!
!   Installers may consider modifying the STOP statement in order to
!   call system-specific exception-handling facilities.
!
!Arguments:
!   SRNAME  (input) CHARACTER*6
!           The name of the routine which called XERBLA.
!
!   INFO   (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!----------------------------------------------------------------
    !.. Scalar Arguments ..
    CHARACTER*6        SRNAME
    INTEGER            INFO
    WRITE( *, FMT = 9999 )SRNAME, INFO
    STOP
    9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',&
                'an illegal value' )
end subroutine XERBLA

function ILAENV(ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!----------------------------------------------------------------
!  -- LAPACK auxiliary routine (version 3.1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     January 2007
!
!Purpose:
!     ILAENV is called from the LAPACK routines to choose problem-dependent
!     parameters for the local environment.  See ISPEC for a description of
!     the parameters.
!   
!     ILAENV returns an INTEGER
!     if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!     if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!   
!     This version provides a set of parameters which should give good,
!     but not optimal, performance on many of the currently available
!     computers.  Users are encouraged to modify this subroutine to set
!     the tuning parameters for their particular machine using the option
!     and problem size information in the arguments.
!   
!     This routine will not function correctly if it is converted to all
!     lower case.  Converting it to all upper case is allowed.
!   
!Arguments:
!     ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines (DEPRECATED)
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR method
!               for nonsymmetric eigenvalue problems (DEPRECATED)
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!          12 <= ISPEC <= 16:
!               xHSEQR or one of its subroutines,
!               see IPARMQ for detailed explanation
!
!     NAME     (input) CHARACTER*(*)
!              The name of the calling subroutine, in either upper case or
!              lower case.
!
!     OPTS    (input) CHARACTER*(*)
!             The character options to the subroutine NAME, concatenated
!             into a single character string.  For example, UPLO = 'U',
!             TRANS = 'T', and DIAG = 'N' for a triangular routine would
!             be specified as OPTS = 'UTN'.
!
!    N1       (input) INTEGER
!    N2       (input) INTEGER
!    N3       (input) INTEGER
!    N4       (input) INTEGER
!             Problem dimensions for the subroutine NAME; these may not all
!             be required.
!
!Further Details:
!     The following conventions have been used when calling ILAENV from the
!     LAPACK routines:
!     1)  OPTS is a concatenation of all of the character options to
!         subroutine NAME, in the same order that they appear in the
!         argument list for NAME, even if they are not used in determining
!         the value of the parameter specified by ISPEC.
!     2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!         that they appear in the argument list for NAME.  N1 is used
!         first, N2 second, and so on, and unused problem dimensions are
!         passed a value of -1.
!     3)  The parameter value returned by ILAENV is checked for validity in
!         the calling subroutine.  For example, ILAENV is used to retrieve
!         the optimal blocksize for STRTRI as follows:
!   
!         NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!         IF( NB.LE.1 ) NB = MAX( 1, N )
!----------------------------------------------------------------
    INTEGER ILAENV
    !.. Scalar Arguments ..
    CHARACTER*( * )    NAME, OPTS
    INTEGER            ISPEC, N1, N2, N3, N4
    !.. Local Scalars ..
    INTEGER            I, IC, IZ, NB, NBMIN, NX
    LOGICAL            CNAME, SNAME
    CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
    !.. Intrinsic Functions ..
    INTRINSIC          CHAR, ICHAR, INT, MIN, REAL

    !.. External Functions ..
    !INTEGER            IEEECK, IPARMQ
    !EXTERNAL           IEEECK, IPARMQ
    !.. Executable Statements ..
    GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
            130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
    !Invalid value for ISPEC
    ILAENV = -1
    RETURN

    10 CONTINUE
    ! Convert NAME to upper case if the first character is lower case.
    ILAENV = 1
    SUBNAM = NAME
    IC = ICHAR( SUBNAM( 1: 1 ) )
    IZ = ICHAR( 'Z' )
    IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
        ! ASCII character set
        IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) SUBNAM( I: I ) = CHAR( IC-32 )
            20  CONTINUE
        END IF
    ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
        ! EBCDIC character set
        IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
            ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
            ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
                SUBNAM( 1: 1 ) = CHAR( IC+64 )
                DO 30 I = 2, 6
                    IC = ICHAR( SUBNAM( I: I ) )
                    IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                        ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                        ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
                30  CONTINUE
        END IF
    ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
        ! Prime machines:  ASCII+128
        IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )SUBNAM( I: I ) = CHAR( IC-32 )
            40  CONTINUE
        END IF
    END IF
    C1 = SUBNAM( 1: 1 )
    SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
    CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
    IF( .NOT.( CNAME .OR. SNAME ) ) RETURN
    C2 = SUBNAM( 2: 3 )
    C3 = SUBNAM( 4: 6 )
    C4 = C3( 2: 3 )
    GO TO ( 50, 60, 70 )ISPEC

    50 CONTINUE
    ! ISPEC = 1:  block size
    ! In these examples, separate code is provided for setting NB for
    ! real and complex.  We assume that NB will take the same value in
    ! single or double precision.
    NB = 1
    IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
        ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
        ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
        END IF
    ELSE IF( C2.EQ.'PO' ) THEN
        IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
        END IF
    ELSE IF( C2.EQ.'SY' ) THEN
        IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
        ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
        END IF
    ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRF' ) THEN
            NB = 64
        ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
        ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
        END IF
    ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NB = 32
            END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
        END IF
    ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &   
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NB = 32
            END IF
        END IF
    ELSE IF( C2.EQ.'GB' ) THEN
        IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
                IF( N4.LE.64 ) THEN
                    NB = 1
                ELSE
                    NB = 32
               END IF
            ELSE
                IF( N4.LE.64 ) THEN
                    NB = 1
                ELSE
                    NB = 32
               END IF
            END IF
        END IF
    ELSE IF( C2.EQ.'PB' ) THEN
        IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
                IF( N2.LE.64 ) THEN
                    NB = 1
                ELSE
                    NB = 32
                END IF
            ELSE
                IF( N2.LE.64 ) THEN
                    NB = 1
                ELSE
                    NB = 32
                END IF
            END IF
        END IF
    ELSE IF( C2.EQ.'TR' ) THEN
        IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
                NB = 64
            ELSE
                NB = 64
            END IF
        END IF
    ELSE IF( C2.EQ.'LA' ) THEN
        IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
        END IF
    ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
        IF( C3.EQ.'EBZ' ) THEN
            NB = 1
        END IF
    END IF
    ILAENV = NB
    RETURN

    60 CONTINUE
    ! ISPEC = 2:  minimum block size
    NBMIN = 2
    IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
                NBMIN = 2
            ELSE
                NBMIN = 2
            END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
                NBMIN = 2
            ELSE
                NBMIN = 2
            END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
                NBMIN = 2
            ELSE
                NBMIN = 2
            END IF
        ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
                NBMIN = 2
            ELSE
                NBMIN = 2
            END IF
        END IF
    ELSE IF( C2.EQ.'SY' ) THEN
        IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
                NBMIN = 8
            ELSE
                NBMIN = 8
            END IF
        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
                NBMIN = 2
        END IF
    ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
        END IF
    ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NBMIN = 2
            END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NBMIN = 2
            END IF
        END IF
    ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NBMIN = 2
            END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NBMIN = 2
            END IF
        END IF
    END IF
    ILAENV = NBMIN
    RETURN

    70 CONTINUE
    ! ISPEC = 3:  crossover point
    NX = 0
    IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
                NX = 128
            ELSE
                NX = 128
            END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
                NX = 128
            ELSE
                NX = 128
            END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
                NX = 128
            ELSE
                NX = 128
            END IF
        END IF
    ELSE IF( C2.EQ.'SY' ) THEN
        IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
        END IF
    ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRD' ) THEN
            NX = 32
        END IF
    ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NX = 128
            END IF
        END IF
    ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                NX = 128
            END IF
        END IF
    END IF
    ILAENV = NX
    RETURN

    80 CONTINUE
    ! ISPEC = 4:  number of shifts (used by xHSEQR)
    ILAENV = 6
    RETURN

    90 CONTINUE
    ! ISPEC = 5:  minimum column dimension (not used)
    ILAENV = 2
    RETURN

    100 CONTINUE
    ! ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
    ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
    RETURN

    110 CONTINUE
    ! ISPEC = 7:  number of processors (not used)
    ILAENV = 1
    RETURN

    120 CONTINUE
    ! ISPEC = 8:  crossover point for multishift (used by xHSEQR)
    ILAENV = 50
    RETURN

    130 CONTINUE
    ! ISPEC = 9:  maximum size of the subproblems at the bottom of the
    !             computation tree in the divide-and-conquer algorithm
    !             (used by xGELSD and xGESDD)
    !
    ILAENV = 25
    RETURN

    140 CONTINUE
    ! ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
    ! ILAENV = 0
    ILAENV = 1
    IF( ILAENV.EQ.1 ) THEN
        ILAENV = IEEECK( 0, 0.0, 1.0 )
    END IF
    RETURN

    150 CONTINUE
    ! ISPEC = 11: infinity arithmetic can be trusted not to trap
    ! ILAENV = 0
    ILAENV = 1
    IF( ILAENV.EQ.1 ) THEN
        ILAENV = IEEECK( 1, 0.0, 1.0 )
    END IF
    RETURN

    160 CONTINUE
    ! 12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
    ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
    RETURN
end function ILAENV

subroutine DLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!----------------------------------------------------------------
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!Purpose:
!     DLASYF computes a partial factorization of a real symmetric matrix A
!     using the Bunch-Kaufman diagonal pivoting method. The partial
!     factorization has the form:
!   
!     A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
!           ( 0  U22 ) (  0   D  ) ( U12' U22' )
!   
!     A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
!           ( L21  I ) (  0  A22 ) (  0    I   )
!   
!     where the order of D is at most NB. The actual order is returned in
!     the argument KB, and is either NB or NB-1, or N if N <= NB.
!   
!     DLASYF is an auxiliary routine called by DSYTRF. It uses blocked code
!     (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
!     A22 (if UPLO = 'L').
!
!Arguments:
!     UPLO    (input) CHARACTER*1
!             Specifies whether the upper or lower triangular part of the
!             symmetric matrix A is stored:
!             = 'U':  Upper triangular
!             = 'L':  Lower triangular
!   
!     N       (input) INTEGER
!             The order of the matrix A.  N >= 0.
!   
!     NB      (input) INTEGER
!             The maximum number of columns of the matrix A that should be
!             factored.  NB should be at least 2 to allow for 2-by-2 pivot
!             blocks.
!   
!     KB      (output) INTEGER
!             The number of columns of A that were actually factored.
!             KB is either NB-1 or NB, or N if N <= NB.
!   
!     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!             On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!             n-by-n upper triangular part of A contains the upper
!             triangular part of the matrix A, and the strictly lower
!             triangular part of A is not referenced.  If UPLO = 'L', the
!             leading n-by-n lower triangular part of A contains the lower
!             triangular part of the matrix A, and the strictly upper
!             triangular part of A is not referenced.
!             On exit, A contains details of the partial factorization.
!   
!     LDA     (input) INTEGER
!             The leading dimension of the array A.  LDA >= max(1,N).
!   
!     IPIV    (output) INTEGER array, dimension (N)
!             Details of the interchanges and the block structure of D.
!             If UPLO = 'U', only the last KB elements of IPIV are set;
!             if UPLO = 'L', only the first KB elements are set.
!   
!             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!             interchanged and D(k,k) is a 1-by-1 diagonal block.
!             If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!             columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!             is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!             IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!             interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!   
!     W       (workspace) DOUBLE PRECISION array, dimension (LDW,NB)
!   
!     LDW     (input) INTEGER
!             The leading dimension of the array W.  LDW >= max(1,N).
!   
!     INFO    (output) INTEGER
!             = 0: successful exit
!             > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!                  has been completed, but the block diagonal matrix D is
!                  exactly singular.
!----------------------------------------------------------------   

    !.. Scalar Arguments ..
    CHARACTER          UPLO
    INTEGER            INFO, KB, LDA, LDW, N, NB
    !.. Array Arguments ..
    INTEGER            IPIV( * )
    DOUBLE PRECISION   A( LDA, * ), W( LDW, * )

    !.. Parameters ..
    DOUBLE PRECISION   ZERO, ONE
    PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
    DOUBLE PRECISION   EIGHT, SEVTEN
    PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
    !.. Local Scalars ..
    INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP, KSTEP, KW
    DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1, ROWMAX, T
    !.. External Functions ..
    !LOGICAL            LSAME
    !INTEGER            IDAMAX
    !EXTERNAL           LSAME, IDAMAX
    !.. External Subroutines ..
    !EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
    !.. Intrinsic Functions ..
    INTRINSIC          ABS, MAX, MIN, SQRT
    !.. Executable Statements ..
    INFO = 0
    !Initialize ALPHA for use in choosing pivot block size.
    ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
    IF( LSAME( UPLO, 'U' ) ) THEN
        ! Factorize the trailing columns of A using the upper triangle
        ! of A and working backwards, and compute the matrix W = U12*D
        ! for use in updating A11
        !
        ! K is the main loop index, decreasing from N in steps of 1 or 2
        !
        ! KW is the column of W which corresponds to column K of A
        K = N
        10  CONTINUE
        KW = NB + K - N
        ! Exit from loop
        IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 ) GO TO 30
        ! Copy column K of A to column KW of W and update it
        CALL DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
        IF( K.LT.N ) CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, &
            W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
        KSTEP = 1
        ! Determine rows and columns to be interchanged and whether
        ! a 1-by-1 or 2-by-2 pivot block will be used
        ABSAKK = ABS( W( K, KW ) )
        ! IMAX is the row-index of the largest off-diagonal element in
        ! column K, and COLMAX is its absolute value
        IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
        ELSE
            COLMAX = ZERO
        END IF
        IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            ! Column K is zero: set INFO and continue
            IF( INFO.EQ.0 ) INFO = K
            KP = K
        ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
                ! no interchange, use 1-by-1 pivot block
               KP = K
            ELSE
                ! Copy column IMAX to column KW-1 of W and update it
                CALL DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
                CALL DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA, W( IMAX+1, KW-1 ), 1 )
                IF( K.LT.N ) CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), &
                    LDA, W( IMAX, KW+1 ), LDW, ONE, W( 1, KW-1 ), 1 )
                ! JMAX is the column-index of the largest off-diagonal
                ! element in row IMAX, and ROWMAX is its absolute value
                JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
                ROWMAX = ABS( W( JMAX, KW-1 ) )
                IF( IMAX.GT.1 ) THEN
                    JMAX = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                    ROWMAX = MAX( ROWMAX, ABS( W( JMAX, KW-1 ) ) )
                END IF
                IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                    ! no interchange, use 1-by-1 pivot block
                    KP = K
                ELSE IF( ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX ) THEN
                    ! interchange rows and columns K and IMAX, use 1-by-1
                    ! pivot block
                    KP = IMAX
                    ! copy column KW-1 of W to column KW
                    CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
                ELSE
                    ! interchange rows and columns K-1 and IMAX, use 2-by-2
                    ! pivot block
                    KP = IMAX
                    KSTEP = 2
                END IF
            END IF
            KK = K - KSTEP + 1
            KKW = NB + KK - N
            ! Updated column KP is already stored in column KKW of W
            IF( KP.NE.KK ) THEN
                ! Copy non-updated column KK to column KP
                A( KP, K ) = A( KK, K )
                CALL DCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),LDA )
                CALL DCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
                ! Interchange rows KK and KP in last KK columns of A and W
                CALL DSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
                CALL DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW )
            END IF
            IF( KSTEP.EQ.1 ) THEN
                ! 1-by-1 pivot block D(k): column KW of W now holds
                ! W(k) = U(k)*D(k)
                ! where U(k) is the k-th column of U
                ! Store U(k) in column k of A
                CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
                R1 = ONE / A( K, K )
                CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
                ! 2-by-2 pivot block D(k): columns KW and KW-1 of W now hold
                ! ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
                ! where U(k) and U(k-1) are the k-th and (k-1)-th columns of U
                IF( K.GT.2 ) THEN
                    ! Store U(k) and U(k-1) in columns k and k-1 of A
                    D21 = W( K-1, KW )
                    D11 = W( K, KW ) / D21
                    D22 = W( K-1, KW-1 ) / D21
                    T = ONE / ( D11*D22-ONE )
                    D21 = T / D21
                    DO 20 J = 1, K - 2
                        A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                        A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
                    20  CONTINUE
                END IF
                ! Copy D(k) to A
                A( K-1, K-1 ) = W( K-1, KW-1 )
                A( K-1, K ) = W( K-1, KW )
                A( K, K ) = W( K, KW )
            END IF
        END IF
        ! Store details of the interchanges in IPIV
        IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
        ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
        END IF
        ! Decrease K and return to the start of the main loop
        K = K - KSTEP
        GO TO 10

        30 CONTINUE
        ! Update the upper triangle of A11 (= A(1:k,1:k)) as
        ! A11 := A11 - U12*D*U12' = A11 - U12*W'
        ! computing blocks of NB columns at a time
        DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
            ! Update the upper triangle of the diagonal block
            DO 40 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', JJ-J+1, N-K, -ONE, &
                          A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE, &
                          A( J, JJ ), 1 )
            40  CONTINUE
            ! Update the rectangular superdiagonal block
            CALL DGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE,&
                       A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE, &
                       A( 1, J ), LDA )
        50  CONTINUE
        ! Put U12 in standard form by partially undoing the interchanges
        ! in columns k+1:n
        J = K + 1
        60  CONTINUE
        JJ = J
        JP = IPIV( J )
        IF( JP.LT.0 ) THEN
            JP = -JP
            J = J + 1
        END IF
        J = J + 1
        IF( JP.NE.JJ .AND. J.LE.N )CALL DSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
        IF( J.LE.N ) GO TO 60
        ! Set KB to the number of columns factorized
        KB = N - K
    ELSE
        ! Factorize the leading columns of A using the lower triangle
        ! of A and working forwards, and compute the matrix W = L21*D
        ! for use in updating A22
        ! K is the main loop index, increasing from 1 in steps of 1 or 2
        K = 1
        70 CONTINUE
        ! Exit from loop
        IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N ) GO TO 90
        ! Copy column K of A to column K of W and update it
        CALL DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
        CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA,&
                    W( K, 1 ), LDW, ONE, W( K, K ), 1 )
        KSTEP = 1
        ! Determine rows and columns to be interchanged and whether
        ! a 1-by-1 or 2-by-2 pivot block will be used

        ABSAKK = ABS( W( K, K ) )

        ! IMAX is the row-index of the largest off-diagonal element in
        ! column K, and COLMAX is its absolute value
        IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
        ELSE
            COLMAX = ZERO
        END IF
        IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            ! Column K is zero: set INFO and continue
            IF( INFO.EQ.0 ) INFO = K
                KP = K
        ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
            ! no interchange, use 1-by-1 pivot block
               KP = K
            ELSE
                ! Copy column IMAX to column K+1 of W and update it
                CALL DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1)
                CALL DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ), 1)
                CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), &
                          LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 )
                ! JMAX is the column-index of the largest off-diagonal
                ! element in row IMAX, and ROWMAX is its absolute value
                JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
                ROWMAX = ABS( W( JMAX, K+1 ) )
                IF( IMAX.LT.N ) THEN
                    JMAX = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                    ROWMAX = MAX( ROWMAX, ABS( W( JMAX, K+1 ) ) )
                END IF
                IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                ! no interchange, use 1-by-1 pivot block
                    KP = K
                ELSE IF( ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX ) THEN
                    ! interchange rows and columns K and IMAX, use 1-by-1
                    ! pivot block
                    KP = IMAX
                    ! copy column K+1 of W to column K
                    CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
                ELSE
                    ! interchange rows and columns K+1 and IMAX, use 2-by-2
                    ! pivot block
                    KP = IMAX
                    KSTEP = 2
                END IF
            END IF
            KK = K + KSTEP - 1
            ! Updated column KP is already stored in column KK of W
            IF( KP.NE.KK ) THEN
                ! Copy non-updated column KK to column KP
                A( KP, K ) = A( KK, K )
                CALL DCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
                CALL DCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
                ! Interchange rows KK and KP in first KK columns of A and W
                CALL DSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
                CALL DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
            IF( KSTEP.EQ.1 ) THEN
                ! 1-by-1 pivot block D(k): column k of W now holds
                ! W(k) = L(k)*D(k)
                ! where L(k) is the k-th column of L
                ! Store L(k) in column k of A
                CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
                IF( K.LT.N ) THEN
                    R1 = ONE / A( K, K )
                    CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
                END IF
            ELSE
                ! 2-by-2 pivot block D(k): columns k and k+1 of W now hold
                ! ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
                ! where L(k) and L(k+1) are the k-th and (k+1)-th columns
                ! of L
                IF( K.LT.N-1 ) THEN
                    ! Store L(k) and L(k+1) in columns k and k+1 of A
                    D21 = W( K+1, K )
                    D11 = W( K+1, K+1 ) / D21
                    D22 = W( K, K ) / D21
                    T = ONE / ( D11*D22-ONE )
                    D21 = T / D21
                    DO 80 J = K + 2, N
                        A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                        A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
                    80 CONTINUE
                END IF
                ! Copy D(k) to A
                A( K, K ) = W( K, K )
                A( K+1, K ) = W( K+1, K )
                A( K+1, K+1 ) = W( K+1, K+1 )
            END IF
        END IF
        ! Store details of the interchanges in IPIV
        IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
        ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
        END IF
        ! Increase K and return to the start of the main loop
        K = K + KSTEP
        GO TO 70

        90  CONTINUE
        ! Update the lower triangle of A22 (= A(k:n,k:n)) as
        ! A22 := A22 - L21*D*L21' = A22 - L21*W'
        ! computing blocks of NB columns at a time
        DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
            ! Update the lower triangle of the diagonal block
            DO 100 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', J+JB-JJ, K-1, -ONE, &
                          A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE, &
                          A( JJ, JJ ), 1 )
            100  CONTINUE
            ! Update the rectangular subdiagonal block
            IF( J+JB.LE.N ) CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB, &
                                        K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW, &
                                        ONE, A( J+JB, J ), LDA )
        110  CONTINUE
        ! Put L21 in standard form by partially undoing the interchanges
        ! in columns 1:k-1
         J = K - 1
        120  CONTINUE
        JJ = J
        JP = IPIV( J )
        IF( JP.LT.0 ) THEN
            JP = -JP
            J = J - 1
        END IF
        J = J - 1
        IF( JP.NE.JJ .AND. J.GE.1 ) CALL DSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
        IF( J.GE.1 ) GO TO 120
        ! Set KB to the number of columns factorized
         KB = K - 1
    END IF
    RETURN
end subroutine DLASYF

subroutine DCOPY(N,DX,INCX,DY,INCY)
!---------------------------------------------------------------
!Purpose:
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!---------------------------------------------------------------
    !.. Scalar Arguments ..
    INTEGER INCX,INCY,N
    !.. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)

    !.. Local Scalars ..
    INTEGER I,IX,IY,M,MP1
    !.. Intrinsic Functions ..
    INTRINSIC MOD
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20

    ! code for unequal increments or equal increments
    !not equal to 1
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
    10 CONTINUE
    RETURN
    ! code for both increments equal to 1
    ! clean-up loop
    20 M = MOD(N,7)
    IF (M.EQ.0) GO TO 40
    DO 30 I = 1,M
        DY(I) = DX(I)
    30 CONTINUE
    IF (N.LT.7) RETURN
    40 MP1 = M + 1
    DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
    50 CONTINUE
    RETURN
end subroutine DCOPY

subroutine DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!----------------------------------------------------------------
!Purpose:
!     DGEMV  performs one of the matrix-vector operations
!        y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!     where alpha and beta are scalars, x and y are vectors and A is an
!     m by n matrix.
!   
!Arguments:   
!     TRANS  - CHARACTER*1.
!              On entry, TRANS specifies the operation to be performed as
!              follows:
!   
!                 TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!   
!                 TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!   
!                 TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!   
!              Unchanged on exit.
!   
!     M      - INTEGER.
!              On entry, M specifies the number of rows of the matrix A.
!              M must be at least zero.
!              Unchanged on exit.
!   
!     N      - INTEGER.
!              On entry, N specifies the number of columns of the matrix A.
!              N must be at least zero.
!              Unchanged on exit.
!   
!     ALPHA  - DOUBLE PRECISION.
!              On entry, ALPHA specifies the scalar alpha.
!              Unchanged on exit.
!   
!     A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!              Before entry, the leading m by n part of the array A must
!              contain the matrix of coefficients.
!              Unchanged on exit.
!   
!     LDA    - INTEGER.
!              On entry, LDA specifies the first dimension of A as declared
!              in the calling (sub) program. LDA must be at least
!              max( 1, m ).
!              Unchanged on exit.
!   
!     X      - DOUBLE PRECISION array of DIMENSION at least
!              ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!              and at least
!              ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!              Before entry, the incremented array X must contain the
!              vector x.
!              Unchanged on exit.
!   
!     INCX   - INTEGER.
!              On entry, INCX specifies the increment for the elements of
!              X. INCX must not be zero.
!              Unchanged on exit.
!   
!     BETA   - DOUBLE PRECISION.
!              On entry, BETA specifies the scalar beta. When BETA is
!              supplied as zero then Y need not be set on input.
!              Unchanged on exit.
!   
!     Y      - DOUBLE PRECISION array of DIMENSION at least
!              ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!              and at least
!              ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!              Before entry with BETA non-zero, the incremented array Y
!              must contain the vector y. On exit, Y is overwritten by the
!              updated vector y.
!   
!     INCY   - INTEGER.
!              On entry, INCY specifies the increment for the elements of
!              Y. INCY must not be zero.
!              Unchanged on exit.
!   
!   
!     Level 2 Blas routine.
!   
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs
!----------------------------------------------------------------
    !.. Scalar Arguments ..
    DOUBLE PRECISION ALPHA,BETA
    INTEGER INCX,INCY,LDA,M,N
    CHARACTER TRANS
    !.. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*),Y(*)

    !.. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    !.. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
    !.. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    !.. External Subroutines ..
    !EXTERNAL XERBLA
    !.. Intrinsic Functions ..
    INTRINSIC MAX

    ! Test the input parameters.
    INFO = 0
    IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.LSAME(TRANS,'C')) THEN
        INFO = 1
    ELSE IF (M.LT.0) THEN
        INFO = 2
    ELSE IF (N.LT.0) THEN
        INFO = 3
    ELSE IF (LDA.LT.MAX(1,M)) THEN
        INFO = 6
    ELSE IF (INCX.EQ.0) THEN
        INFO = 8
    ELSE IF (INCY.EQ.0) THEN
        INFO = 11
    END IF
        IF (INFO.NE.0) THEN
        CALL XERBLA('DGEMV ',INFO)
        RETURN
    END IF

    ! Quick return if possible.
    IF ((M.EQ.0) .OR. (N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN

    ! Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    ! up the start points in  X  and  Y.
    IF (LSAME(TRANS,'N')) THEN
        LENX = N
        LENY = M
    ELSE
        LENX = M
        LENY = N
    END IF
    IF (INCX.GT.0) THEN
        KX = 1
    ELSE
        KX = 1 - (LENX-1)*INCX
    END IF
    IF (INCY.GT.0) THEN
        KY = 1
    ELSE
        KY = 1 - (LENY-1)*INCY
    END IF

    ! Start the operations. In this version the elements of A are
    ! accessed sequentially with one pass through A.

    ! First form  y := beta*y.
    IF (BETA.NE.ONE) THEN
        IF (INCY.EQ.1) THEN
            IF (BETA.EQ.ZERO) THEN
                DO 10 I = 1,LENY
                    Y(I) = ZERO
                10 CONTINUE
            ELSE
                DO 20 I = 1,LENY
                    Y(I) = BETA*Y(I)
                20 CONTINUE
            END IF
        ELSE
            IY = KY
            IF (BETA.EQ.ZERO) THEN
                DO 30 I = 1,LENY
                    Y(IY) = ZERO
                    IY = IY + INCY
                30 CONTINUE
            ELSE
                DO 40 I = 1,LENY
                    Y(IY) = BETA*Y(IY)
                    IY = IY + INCY
                40 CONTINUE
            END IF
        END IF
    END IF
    IF (ALPHA.EQ.ZERO) RETURN
    IF (LSAME(TRANS,'N')) THEN
        ! Form  y := alpha*A*x + y.
        JX = KX
        IF (INCY.EQ.1) THEN
            DO 60 J = 1,N
                IF (X(JX).NE.ZERO) THEN
                    TEMP = ALPHA*X(JX)
                    DO 50 I = 1,M
                        Y(I) = Y(I) + TEMP*A(I,J)
                    50 CONTINUE
                END IF
                    JX = JX + INCX
            60 CONTINUE
        ELSE
            DO 80 J = 1,N
                IF (X(JX).NE.ZERO) THEN
                    TEMP = ALPHA*X(JX)
                    IY = KY
                    DO 70 I = 1,M
                        Y(IY) = Y(IY) + TEMP*A(I,J)
                        IY = IY + INCY
                    70 CONTINUE
                END IF
                    JX = JX + INCX
            80 CONTINUE
        END IF
    ELSE
        ! Form  y := alpha*A'*x + y.
        JY = KY
        IF (INCX.EQ.1) THEN
            DO 100 J = 1,N
                TEMP = ZERO
                DO 90 I = 1,M
                    TEMP = TEMP + A(I,J)*X(I)
                90 CONTINUE
                    Y(JY) = Y(JY) + ALPHA*TEMP
                    JY = JY + INCY
            100 CONTINUE
        ELSE
            DO 120 J = 1,N
                TEMP = ZERO
                IX = KX
                DO 110 I = 1,M
                    TEMP = TEMP + A(I,J)*X(IX)
                    IX = IX + INCX
                110 CONTINUE
                Y(JY) = Y(JY) + ALPHA*TEMP
                JY = JY + INCY
            120 CONTINUE
        END IF
    END IF
    RETURN
end subroutine DGEMV

function IDAMAX(N,DX,INCX)
!----------------------------------------------------------------
! Purpose:
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!----------------------------------------------------------------
    INTEGER IDAMAX
    ! .. Scalar Arguments ..
    INTEGER INCX,N
    ! .. Array Arguments ..
    DOUBLE PRECISION DX(*)

    ! .. Local Scalars ..
    DOUBLE PRECISION DMAX
    INTEGER I,IX
    ! .. Intrinsic Functions ..
    INTRINSIC DABS

    IDAMAX = 0
    IF (N.LT.1 .OR. INCX.LE.0) RETURN
    IDAMAX = 1
    IF (N.EQ.1) RETURN
    IF (INCX.EQ.1) GO TO 20
    ! code for increment not equal to 1
    IX = 1
    DMAX = DABS(DX(1))
    IX = IX + INCX
    DO 10 I = 2,N
        IF (DABS(DX(IX)).LE.DMAX) GO TO 5
        IDAMAX = I
        DMAX = DABS(DX(IX))
    5   IX = IX + INCX
    10 CONTINUE
    RETURN
    ! code for increment equal to 1
    20 DMAX = DABS(DX(1))
        DO 30 I = 2,N
            IF (DABS(DX(I)).LE.DMAX) GO TO 30
            IDAMAX = I
            DMAX = DABS(DX(I))
    30 CONTINUE
    RETURN
end function IDAMAX

subroutine DSWAP(N,DX,INCX,DY,INCY)
!----------------------------------------------------------------
! Purpose:
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!----------------------------------------------------------------
    ! .. Scalar Arguments ..
    INTEGER INCX,INCY,N
    ! .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)

    ! .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY,M,MP1
    ! .. Intrinsic Functions ..
    INTRINSIC MOD

    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20

    ! code for unequal increments or equal increments not equal
    ! to 1
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
    10 CONTINUE
    RETURN
    ! code for both increments equal to 1
    ! clean-up loop
    20 M = MOD(N,3)
    IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
        30 CONTINUE
        IF (N.LT.3) RETURN
    40 MP1 = M + 1
    DO 50 I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
    50 CONTINUE
    RETURN
end subroutine DSWAP

subroutine DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!----------------------------------------------------------------
!Purpose:
!   DGEMM  performs one of the matrix-matrix operations
!
!       C := alpha*op( A )*op( B ) + beta*C,
!
!   where  op( X ) is one of
!       op( X ) = X   or   op( X ) = X',
!
!   alpha and beta are scalars, and A, B and C are matrices, with op( A )
!   an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!Arguments:
!   TRANSA - CHARACTER*1.
!       On entry, TRANSA specifies the form of op( A ) to be used in
!       the matrix multiplication as follows:
!           TRANSA = 'N' or 'n',  op( A ) = A.
!
!           TRANSA = 'T' or 't',  op( A ) = A'.
!
!           TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!   TRANSB - CHARACTER*1.
!       On entry, TRANSB specifies the form of op( B ) to be used in
!       the matrix multiplication as follows:
!
!           TRANSB = 'N' or 'n',  op( B ) = B.
!
!           TRANSB = 'T' or 't',  op( B ) = B'.
!
!           TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!   M      - INTEGER.
!            On entry,  M  specifies  the number  of rows  of the  matrix
!            op( A )  and of the  matrix  C.  M  must  be at least  zero.
!            Unchanged on exit.
!
!   N      - INTEGER.
!             On entry,  N  specifies the number  of columns of the matrix
!             op( B ) and the number of columns of the matrix C. N must be
!             at least zero.
!             Unchanged on exit.
!   
!   K      - INTEGER.
!            On entry,  K  specifies  the number of columns of the matrix
!            op( A ) and the number of rows of the matrix op( B ). K must
!            be at least  zero.
!            Unchanged on exit.
!   
!   ALPHA  - DOUBLE PRECISION.
!            On entry, ALPHA specifies the scalar alpha.
!            Unchanged on exit.
!   
!   A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!            k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!            Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!            part of the array  A  must contain the matrix  A,  otherwise
!            the leading  k by m  part of the array  A  must contain  the
!            matrix A.
!            Unchanged on exit.
!
!   LDA    - INTEGER.
!            On entry, LDA specifies the first dimension of A as declared
!            in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!            LDA must be at least  max( 1, m ), otherwise  LDA must be at
!            least  max( 1, k ).
!            Unchanged on exit.
!
!   B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!            n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!            Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!            part of the array  B  must contain the matrix  B,  otherwise
!            the leading  n by k  part of the array  B  must contain  the
!            matrix B.
!            Unchanged on exit.
!
!   LDB    - INTEGER.
!            On entry, LDB specifies the first dimension of B as declared
!            in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!            LDB must be at least  max( 1, k ), otherwise  LDB must be at
!            least  max( 1, n ).
!            Unchanged on exit.
!
!   BETA   - DOUBLE PRECISION.
!            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!            supplied as zero then C need not be set on input.
!            Unchanged on exit.
!
!   C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!            Before entry, the leading  m by n  part of the array  C must
!            contain the matrix  C,  except when  beta  is zero, in which
!            case C need not be set on entry.
!            On exit, the array  C  is overwritten by the  m by n  matrix
!            ( alpha*op( A )*op( B ) + beta*C ).
!
!   LDC    - INTEGER.
!            On entry, LDC specifies the first dimension of C as declared
!            in  the  calling  (sub)  program.   LDC  must  be  at  least
!            max( 1, m ).
!            Unchanged on exit.
!
!
!   Level 3 Blas routine.
!
!   -- Written on 8-February-1989.
!      Jack Dongarra, Argonne National Laboratory.
!      Iain Duff, AERE Harwell.
!      Jeremy Du Croz, Numerical Algorithms Group Ltd.
!      Sven Hammarling, Numerical Algorithms Group Ltd.
!----------------------------------------------------------------
    ! .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA,BETA
    INTEGER K,LDA,LDB,LDC,M,N
    CHARACTER TRANSA,TRANSB
    ! .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)

    ! .. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    ! .. External Subroutines ..
    !EXTERNAL XERBLA
    ! .. Intrinsic Functions ..
    INTRINSIC MAX
    ! .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
    LOGICAL NOTA,NOTB
    ! .. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

    ! Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    ! transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    ! and  columns of  A  and the  number of  rows  of  B  respectively.
    NOTA = LSAME(TRANSA,'N')
    NOTB = LSAME(TRANSB,'N')
    IF (NOTA) THEN
        NROWA = M
        NCOLA = K
    ELSE
        NROWA = K
        NCOLA = M
    END IF
    IF (NOTB) THEN
        NROWB = K
    ELSE
        NROWB = N
    END IF

    ! Test the input parameters.
    INFO = 0
    IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. (.NOT.LSAME(TRANSA,'T'))) THEN
        INFO = 1
    ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. (.NOT.LSAME(TRANSB,'T'))) THEN
        INFO = 2
    ELSE IF (M.LT.0) THEN
        INFO = 3
    ELSE IF (N.LT.0) THEN
        INFO = 4
    ELSE IF (K.LT.0) THEN
        INFO = 5
    ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
        INFO = 8
    ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
        INFO = 10
    ELSE IF (LDC.LT.MAX(1,M)) THEN
        INFO = 13
    END IF
    IF (INFO.NE.0) THEN
        CALL XERBLA('DGEMM ',INFO)
        RETURN
    END IF
    ! Quick return if possible.
    IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
    ! And if  alpha.eq.zero.
    IF (ALPHA.EQ.ZERO) THEN
        IF (BETA.EQ.ZERO) THEN
            DO 20 J = 1,N
                DO 10 I = 1,M
                    C(I,J) = ZERO
                10 CONTINUE
            20 CONTINUE
        ELSE
            DO 40 J = 1,N
                DO 30 I = 1,M
                    C(I,J) = BETA*C(I,J)
                30 CONTINUE
            40 CONTINUE
        END IF
        RETURN
    END IF
    ! Start the operations.
    IF (NOTB) THEN
        IF (NOTA) THEN
        ! Form  C := alpha*A*B + beta*C.
            DO 90 J = 1,N
                IF (BETA.EQ.ZERO) THEN
                    DO 50 I = 1,M
                        C(I,J) = ZERO
                    50 CONTINUE
                ELSE IF (BETA.NE.ONE) THEN
                    DO 60 I = 1,M
                        C(I,J) = BETA*C(I,J)
                    60 CONTINUE
                END IF
                DO 80 L = 1,K
                    IF (B(L,J).NE.ZERO) THEN
                        TEMP = ALPHA*B(L,J)
                        DO 70 I = 1,M
                            C(I,J) = C(I,J) + TEMP*A(I,L)
                        70 CONTINUE
                    END IF
                80 CONTINUE
            90 CONTINUE
        ELSE
        ! Form  C := alpha*A'*B + beta*C
            DO 120 J = 1,N
                DO 110 I = 1,M
                    TEMP = ZERO
                    DO 100 L = 1,K
                        TEMP = TEMP + A(L,I)*B(L,J)
                    100 CONTINUE
                    IF (BETA.EQ.ZERO) THEN
                        C(I,J) = ALPHA*TEMP
                    ELSE
                        C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                    END IF
                110 CONTINUE
            120 CONTINUE
        END IF
    ELSE
        IF (NOTA) THEN
        ! Form  C := alpha*A*B' + beta*C
            DO 170 J = 1,N
                IF (BETA.EQ.ZERO) THEN
                    DO 130 I = 1,M
                        C(I,J) = ZERO
                    130 CONTINUE
                ELSE IF (BETA.NE.ONE) THEN
                    DO 140 I = 1,M
                        C(I,J) = BETA*C(I,J)
                    140 CONTINUE
                END IF
                DO 160 L = 1,K
                    IF (B(J,L).NE.ZERO) THEN
                        TEMP = ALPHA*B(J,L)
                        DO 150 I = 1,M
                            C(I,J) = C(I,J) + TEMP*A(I,L)
                        150 CONTINUE
                    END IF
                160 CONTINUE
            170 CONTINUE
        ELSE
        ! Form  C := alpha*A'*B' + beta*C
            DO 200 J = 1,N
                DO 190 I = 1,M
                    TEMP = ZERO
                    DO 180 L = 1,K
                        TEMP = TEMP + A(L,I)*B(J,L)
                    180 CONTINUE
                    IF (BETA.EQ.ZERO) THEN
                        C(I,J) = ALPHA*TEMP
                    ELSE
                        C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                    END IF
                190 CONTINUE
            200 CONTINUE
        END IF
    END IF
    RETURN
end subroutine DGEMM

subroutine DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
!-------------------------------------------------------------------
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!Purpose:
!     DSYTF2 computes the factorization of a real symmetric matrix A using
!     the Bunch-Kaufman diagonal pivoting method:
!   
!        A = U*D*U'  or  A = L*D*L'
!   
!     where U (or L) is a product of permutation and unit upper (lower)
!     triangular matrices, U' is the transpose of U, and D is symmetric and
!     block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!   
!     This is the unblocked version of the algorithm, calling Level 2 BLAS.
!   
!Arguments:   
!     UPLO    (input) CHARACTER*1
!             Specifies whether the upper or lower triangular part of the
!             symmetric matrix A is stored:
!             = 'U':  Upper triangular
!             = 'L':  Lower triangular
!   
!     N       (input) INTEGER
!             The order of the matrix A.  N >= 0.
!   
!     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!             On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!             n-by-n upper triangular part of A contains the upper
!             triangular part of the matrix A, and the strictly lower
!             triangular part of A is not referenced.  If UPLO = 'L', the
!             leading n-by-n lower triangular part of A contains the lower
!             triangular part of the matrix A, and the strictly upper
!             triangular part of A is not referenced.
!   
!             On exit, the block diagonal matrix D and the multipliers used
!             to obtain the factor U or L (see below for further details).
!   
!     LDA     (input) INTEGER
!             The leading dimension of the array A.  LDA >= max(1,N).
!   
!     IPIV    (output) INTEGER array, dimension (N)
!             Details of the interchanges and the block structure of D.
!             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!             interchanged and D(k,k) is a 1-by-1 diagonal block.
!             If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!             columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!             is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!             IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!             interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!   
!     INFO    (output) INTEGER
!             = 0: successful exit
!             < 0: if INFO = -k, the k-th argument had an illegal value
!             > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!                  has been completed, but the block diagonal matrix D is
!                  exactly singular, and division by zero will occur if it
!                  is used to solve a system of equations.
!   
!Further Details:
!     09-29-06 - patch from
!       Bobby Cheng, MathWorks
!   
!       Replace l.204 and l.372
!            IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!       by
!            IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
!   
!     01-01-96 - Based on modifications by
!       J. Lewis, Boeing Computer Services Company
!       A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!     1-96 - Based on modifications by J. Lewis, Boeing Computer Services
!            Company
!   
!     If UPLO = 'U', then A = U*D*U', where
!        U = P(n)*U(n)* ... *P(k)U(k)* ...,
!     i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!     1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!     and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!     defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!     that if the diagonal block D(k) is of order s (s = 1 or 2), then
!   
!                (   I    v    0   )   k-s
!        U(k) =  (   0    I    0   )   s
!                (   0    0    I   )   n-k
!                   k-s   s   n-k
!   
!     If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!     If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!     and A(k,k), and v overwrites A(1:k-2,k-1:k).
!   
!     If UPLO = 'L', then A = L*D*L', where
!        L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!     i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!     n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!     and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!     defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!     that if the diagonal block D(k) is of order s (s = 1 or 2), then
!   
!                (   I    0     0   )  k-1
!        L(k) =  (   0    I     0   )  s
!                (   0    v     I   )  n-k-s+1
!                   k-1   s  n-k-s+1
!   
!     If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!     If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!     and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!------------------------------------------------------------------
    !.. Scalar Arguments ..
    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    !.. Array Arguments ..
    INTEGER            IPIV( * )
    DOUBLE PRECISION   A( LDA, * )
    !.. Parameters ..
    DOUBLE PRECISION   ZERO, ONE
    PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
    DOUBLE PRECISION   EIGHT, SEVTEN
    PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )

    !.. Local Scalars ..
    LOGICAL            UPPER
    INTEGER            I, IMAX, J, JMAX, K, KK, KP, KSTEP
    DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1,&
                       ROWMAX, T, WK, WKM1, WKP1
    !.. External Functions ..
    !LOGICAL            LSAME, DISNAN
    !INTEGER            IDAMAX
    !EXTERNAL           LSAME, IDAMAX, DISNAN
    !.. External Subroutines ..
    !EXTERNAL           DSCAL, DSWAP, DSYR, XERBLA
    !.. Intrinsic Functions ..
    INTRINSIC          ABS, MAX, SQRT

    !.. Executable Statements ..
    ! Test the input parameters.
    INFO = 0
    UPPER = LSAME( UPLO, 'U' )
    IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
        INFO = -1
    ELSE IF( N.LT.0 ) THEN
        INFO = -2
    ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
        INFO = -4
    END IF
    IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DSYTF2', -INFO )
        RETURN
    END IF
    ! Initialize ALPHA for use in choosing pivot block size.
    ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
    IF( UPPER ) THEN
        ! Factorize A as U*D*U' using the upper triangle of A
        ! K is the main loop index, decreasing from N to 1 in steps of
        ! 1 or 2
        K = N
        10  CONTINUE
        !  If K < 1, exit from loop
        IF( K.LT.1 ) GO TO 70
        KSTEP = 1
        ! Determine rows and columns to be interchanged and whether
        ! a 1-by-1 or 2-by-2 pivot block will be used
        ABSAKK = ABS( A( K, K ) )

        ! IMAX is the row-index of the largest off-diagonal element in
        ! column K, and COLMAX is its absolute value
        IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
        ELSE
            COLMAX = ZERO
        END IF
        IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
            ! Column K is zero or contains a NaN: set INFO and continue
            IF( INFO.EQ.0 ) INFO = K
            KP = K
        ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
                ! no interchange, use 1-by-1 pivot block
                KP = K
            ELSE
                ! JMAX is the column-index of the largest off-diagonal
                ! element in row IMAX, and ROWMAX is its absolute value
                JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
                ROWMAX = ABS( A( IMAX, JMAX ) )
                IF( IMAX.GT.1 ) THEN
                    JMAX = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                    ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
                END IF
                IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                    ! no interchange, use 1-by-1 pivot block
                    KP = K
                ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
                    ! interchange rows and columns K and IMAX, use 1-by-1
                    ! pivot block
                    KP = IMAX
                ELSE
                    ! interchange rows and columns K-1 and IMAX, use 2-by-2
                    ! pivot block
                    KP = IMAX
                    KSTEP = 2
                END IF
            END IF
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
                ! Interchange rows and columns KK and KP in the leading
                ! submatrix A(1:k,1:k)
                CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
                CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ),LDA )
                T = A( KK, KK )
                A( KK, KK ) = A( KP, KP )
                A( KP, KP ) = T
                IF( KSTEP.EQ.2 ) THEN
                    T = A( K-1, K )
                    A( K-1, K ) = A( KP, K )
                    A( KP, K ) = T
                END IF
            END IF
            ! Update the leading submatrix
            IF( KSTEP.EQ.1 ) THEN
                ! 1-by-1 pivot block D(k): column k now holds
                ! W(k) = U(k)*D(k)
                ! where U(k) is the k-th column of U
                ! Perform a rank-1 update of A(1:k-1,1:k-1) as
                ! A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
                R1 = ONE / A( K, K )
                CALL DSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
                ! Store U(k) in column k
                CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
                ! 2-by-2 pivot block D(k): columns k and k-1 now hold
                ! ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
                ! where U(k) and U(k-1) are the k-th and (k-1)-th columns
                ! of U
                !
                ! Perform a rank-2 update of A(1:k-2,1:k-2) as
                ! A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
                !    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
                IF( K.GT.2 ) THEN
                    D12 = A( K-1, K )
                    D22 = A( K-1, K-1 ) / D12
                    D11 = A( K, K ) / D12
                    T = ONE / ( D11*D22-ONE )
                    D12 = T / D12
                    DO 30 J = K - 2, 1, -1
                        WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                        WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                        DO 20 I = J, 1, -1
                            A( I, J ) = A( I, J ) - A( I, K )*WK - &
                                        A( I, K-1 )*WKM1
                        20  CONTINUE
                        A( J, K ) = WK
                        A( J, K-1 ) = WKM1
                    30 CONTINUE
                END IF
            END IF
        END IF
        ! Store details of the interchanges in IPIV
            IF( KSTEP.EQ.1 ) THEN
                IPIV( K ) = KP
            ELSE
                IPIV( K ) = -KP
                IPIV( K-1 ) = -KP
            END IF
        ! Decrease K and return to the start of the main loop
        K = K - KSTEP
        GO TO 10
    ELSE
        ! Factorize A as L*D*L' using the lower triangle of A
        ! K is the main loop index, increasing from 1 to N in steps of
        ! 1 or 2
        K = 1
        40  CONTINUE
        ! If K > N, exit from loop
        IF( K.GT.N ) GO TO 70
        KSTEP = 1
        ! Determine rows and columns to be interchanged and whether
        ! a 1-by-1 or 2-by-2 pivot block will be used
        ABSAKK = ABS( A( K, K ) )
        ! IMAX is the row-index of the largest off-diagonal element in
        ! column K, and COLMAX is its absolute value
        IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
        ELSE
            COLMAX = ZERO
        END IF
        IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
            ! Column K is zero or contains a NaN: set INFO and continue
            IF( INFO.EQ.0 ) INFO = K
            KP = K
        ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
                ! no interchange, use 1-by-1 pivot block
                KP = K
            ELSE
                ! JMAX is the column-index of the largest off-diagonal
                ! element in row IMAX, and ROWMAX is its absolute value
                JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
                ROWMAX = ABS( A( IMAX, JMAX ) )
                IF( IMAX.LT.N ) THEN
                    JMAX = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                    ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
                END IF
                IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                    ! no interchange, use 1-by-1 pivot block
                    KP = K
                ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
                    ! interchange rows and columns K and IMAX, use 1-by-1
                    ! pivot block
                    KP = IMAX
                ELSE
                    ! interchange rows and columns K+1 and IMAX, use 2-by-2
                    ! pivot block
                    KP = IMAX
                    KSTEP = 2
                END IF
            END IF
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
                ! Interchange rows and columns KK and KP in the trailing
                ! submatrix A(k:n,k:n)
                IF( KP.LT.N ) CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
                CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),LDA )
                T = A( KK, KK )
                A( KK, KK ) = A( KP, KP )
                A( KP, KP ) = T
                IF( KSTEP.EQ.2 ) THEN
                    T = A( K+1, K )
                    A( K+1, K ) = A( KP, K )
                    A( KP, K ) = T
                END IF
            END IF
            ! Update the trailing submatrix
            IF( KSTEP.EQ.1 ) THEN
                ! 1-by-1 pivot block D(k): column k now holds
                ! W(k) = L(k)*D(k)
                ! where L(k) is the k-th column of L
                IF( K.LT.N ) THEN
                    ! Perform a rank-1 update of A(k+1:n,k+1:n) as
                    ! A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
                    D11 = ONE / A( K, K )
                    CALL DSYR( UPLO, N-K, -D11, A( K+1, K ), 1, A( K+1, K+1 ), LDA )
                    ! Store L(k) in column K
                    CALL DSCAL( N-K, D11, A( K+1, K ), 1 )
               END IF
            ELSE
                ! 2-by-2 pivot block D(k)
                IF( K.LT.N-1 ) THEN
                    ! Perform a rank-2 update of A(k+2:n,k+2:n) as
                    ! A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
                    ! where L(k) and L(k+1) are the k-th and (k+1)-th
                    ! columns of L
                    D21 = A( K+1, K )
                    D11 = A( K+1, K+1 ) / D21
                    D22 = A( K, K ) / D21
                    T = ONE / ( D11*D22-ONE )
                    D21 = T / D21
                    DO 60 J = K + 2, N
                        WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                        WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
                        DO 50 I = J, N
                            A( I, J ) = A( I, J ) - A( I, K )*WK - A( I, K+1 )*WKP1
                        50  CONTINUE
                        A( J, K ) = WK
                        A( J, K+1 ) = WKP1
                    60  CONTINUE
                END IF
            END IF
        END IF
        ! Store details of the interchanges in IPIV
        IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
        ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
        END IF
        ! Increase K and return to the start of the main loop
        K = K + KSTEP
        GO TO 40
    END IF
    70 CONTINUE
    RETURN
end subroutine DSYTF2

subroutine DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
!----------------------------------------------------------------
!Purpose:
!   DSYR   performs the symmetric rank 1 operation
!      A := alpha*x*x' + A,
!   where alpha is a real scalar, x is an n element vector and A is an
!   n by n symmetric matrix.
! 
!   Arguments: 
!   UPLO   - CHARACTER*1.
!            On entry, UPLO specifies whether the upper or lower
!            triangular part of the array A is to be referenced as
!            follows:
! 
!               UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                   is to be referenced.
! 
!               UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                   is to be referenced.
! 
!            Unchanged on exit.
! 
!   N      - INTEGER.
!            On entry, N specifies the order of the matrix A.
!            N must be at least zero.
!            Unchanged on exit.
! 
!   ALPHA  - DOUBLE PRECISION.
!            On entry, ALPHA specifies the scalar alpha.
!            Unchanged on exit.
! 
!   X      - DOUBLE PRECISION array of dimension at least
!            ( 1 + ( n - 1 )*abs( INCX ) ).
!            Before entry, the incremented array X must contain the n
!            element vector x.
!            Unchanged on exit.
! 
!   INCX   - INTEGER.
!            On entry, INCX specifies the increment for the elements of
!            X. INCX must not be zero.
!            Unchanged on exit.
! 
!   A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!            Before entry with  UPLO = 'U' or 'u', the leading n by n
!            upper triangular part of the array A must contain the upper
!            triangular part of the symmetric matrix and the strictly
!            lower triangular part of A is not referenced. On exit, the
!            upper triangular part of the array A is overwritten by the
!            upper triangular part of the updated matrix.
!            Before entry with UPLO = 'L' or 'l', the leading n by n
!            lower triangular part of the array A must contain the lower
!            triangular part of the symmetric matrix and the strictly
!            upper triangular part of A is not referenced. On exit, the
!            lower triangular part of the array A is overwritten by the
!            lower triangular part of the updated matrix.
! 
!   LDA    - INTEGER.
!            On entry, LDA specifies the first dimension of A as declared
!            in the calling (sub) program. LDA must be at least
!            max( 1, n ).
!            Unchanged on exit.
! 
! 
!   Level 2 Blas routine.
! 
!   -- Written on 22-October-1986.
!      Jack Dongarra, Argonne National Lab.
!      Jeremy Du Croz, Nag Central Office.
!      Sven Hammarling, Nag Central Office.
!      Richard Hanson, Sandia National Labs.
!----------------------------------------------------------------
    ! .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA
    INTEGER INCX,LDA,N
    CHARACTER UPLO
    ! .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*)

    ! .. Parameters ..
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D+0)
    ! .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,IX,J,JX,KX
    ! .. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    ! .. External Subroutines ..
    !EXTERNAL XERBLA
    ! .. Intrinsic Functions ..
    INTRINSIC MAX    

    ! Test the input parameters.
    INFO = 0
    IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
        INFO = 1
    ELSE IF (N.LT.0) THEN
        INFO = 2
    ELSE IF (INCX.EQ.0) THEN
        INFO = 5
    ELSE IF (LDA.LT.MAX(1,N)) THEN
        INFO = 7
    END IF
    IF (INFO.NE.0) THEN
        CALL XERBLA('DSYR  ',INFO)
        RETURN
    END IF
    ! Quick return if possible.
    IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
    ! Set the start point in X if the increment is not unity.
    IF (INCX.LE.0) THEN
        KX = 1 - (N-1)*INCX
    ELSE IF (INCX.NE.1) THEN
        KX = 1
    END IF
    ! Start the operations. In this version the elements of A are
    ! accessed sequentially with one pass through the triangular part
    ! of A.
    IF (LSAME(UPLO,'U')) THEN
    ! Form  A  when A is stored in upper triangle.
        IF (INCX.EQ.1) THEN
            DO 20 J = 1,N
                IF (X(J).NE.ZERO) THEN
                    TEMP = ALPHA*X(J)
                    DO 10 I = 1,J
                        A(I,J) = A(I,J) + X(I)*TEMP
                    10 CONTINUE
                END IF
            20 CONTINUE
        ELSE
            JX = KX
            DO 40 J = 1,N
                IF (X(JX).NE.ZERO) THEN
                    TEMP = ALPHA*X(JX)
                    IX = KX
                    DO 30 I = 1,J
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
                    30 CONTINUE
                END IF
                JX = JX + INCX
            40 CONTINUE
        END IF
    ELSE
    ! Form  A  when A is stored in lower triangle.
        IF (INCX.EQ.1) THEN
            DO 60 J = 1,N
                IF (X(J).NE.ZERO) THEN
                    TEMP = ALPHA*X(J)
                    DO 50 I = J,N
                        A(I,J) = A(I,J) + X(I)*TEMP
                    50 CONTINUE
                END IF
            60 CONTINUE
        ELSE
            JX = KX
            DO 80 J = 1,N
                IF (X(JX).NE.ZERO) THEN
                    TEMP = ALPHA*X(JX)
                    IX = JX
                    DO 70 I = J,N
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
                    70 CONTINUE
                END IF
                JX = JX + INCX
            80 CONTINUE
        END IF
    END IF
    RETURN
end subroutine DSYR

function DISNAN(DIN)
!----------------------------------------------------------------
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!Purpose:
!   DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!   otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!   future.
!Arguments:
!   DIN     (input) DOUBLE PRECISION
!           Input to test for NaN.
!----------------------------------------------------------------
    LOGICAL DISNAN
    ! .. Scalar Arguments ..
    DOUBLE PRECISION DIN
    ! .. External Functions ..
    !LOGICAL DLAISNAN
    !EXTERNAL DLAISNAN
    ! .. Executable Statements ..
    DISNAN = DLAISNAN(DIN,DIN)
    RETURN
end function DISNAN

function DLAISNAN(DIN1,DIN2)
!----------------------------------------------------------------
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!Purpose:
!   This routine is not for general use.  It exists solely to avoid
!   over-optimization in DISNAN.
!   
!   DLAISNAN checks for NaNs by comparing its two arguments for
!   inequality.  NaN is the only floating-point value where NaN != NaN
!   returns .TRUE.  To check for NaNs, pass the same variable as both
!   arguments.
!   
!   Strictly speaking, Fortran does not allow aliasing of function
!   arguments. So a compiler must assume that the two arguments are
!   not the same variable, and the test will not be optimized away.
!   Interprocedural or whole-program optimization may delete this
!   test.  The ISNAN functions will be replaced by the correct
!   Fortran 03 intrinsic once the intrinsic is widely available.
!   
!Arguments:  
!   DIN1     (input) DOUBLE PRECISION
!   DIN2     (input) DOUBLE PRECISION
!            Two numbers to compare for inequality.
!----------------------------------------------------------------
    LOGICAL DLAISNAN
    ! .. Scalar Arguments ..
    DOUBLE PRECISION DIN1,DIN2

    ! .. Executable Statements ..
    DLAISNAN = (DIN1.NE.DIN2)
    RETURN
end function DLAISNAN

subroutine DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
!----------------------------------------------------------------
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!Purpose:   
!     DSYTRI computes the inverse of a real symmetric indefinite matrix
!     A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!     DSYTRF.
!   
!Arguments:   
!     UPLO    (input) CHARACTER*1
!             Specifies whether the details of the factorization are stored
!             as an upper or lower triangular matrix.
!             = 'U':  Upper triangular, form is A = U*D*U**T;
!             = 'L':  Lower triangular, form is A = L*D*L**T.
!   
!     N       (input) INTEGER
!             The order of the matrix A.  N >= 0.
!   
!     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!             On entry, the block diagonal matrix D and the multipliers
!             used to obtain the factor U or L as computed by DSYTRF.
!   
!             On exit, if INFO = 0, the (symmetric) inverse of the original
!             matrix.  If UPLO = 'U', the upper triangular part of the
!             inverse is formed and the part of A below the diagonal is not
!             referenced; if UPLO = 'L' the lower triangular part of the
!             inverse is formed and the part of A above the diagonal is
!             not referenced.
!   
!     LDA     (input) INTEGER
!             The leading dimension of the array A.  LDA >= max(1,N).
!   
!     IPIV    (input) INTEGER array, dimension (N)
!             Details of the interchanges and the block structure of D
!             as determined by DSYTRF.
!   
!     WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!   
!     INFO    (output) INTEGER
!             = 0: successful exit
!             < 0: if INFO = -i, the i-th argument had an illegal value
!             > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!                  inverse could not be computed.
!----------------------------------------------------------------
    !.. Scalar Arguments ..
    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    !.. Array Arguments ..
    INTEGER            IPIV( * )
    DOUBLE PRECISION   A( LDA, * ), WORK( * )

    !.. Parameters ..
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
    !.. Local Scalars ..
    LOGICAL            UPPER
    INTEGER            K, KP, KSTEP
    DOUBLE PRECISION   AK, AKKP1, AKP1, D, T, TEMP
    !.. External Functions ..
    !LOGICAL            LSAME
    !DOUBLE PRECISION   DDOT
    !EXTERNAL           LSAME, DDOT
    !.. External Subroutines ..
    !EXTERNAL           DCOPY, DSWAP, DSYMV, XERBLA
    !.. Intrinsic Functions ..
    INTRINSIC          ABS, MAX
    !.. Executable Statements ..

    ! Test the input parameters.
    INFO = 0
    UPPER = LSAME( UPLO, 'U' )
    IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
        INFO = -1
    ELSE IF( N.LT.0 ) THEN
        INFO = -2
    ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
        INFO = -4
    END IF
    IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DSYTRI', -INFO )
        RETURN
    END IF
    ! Quick return if possible
    IF( N.EQ.0 ) RETURN
    !  Check that the diagonal matrix D is nonsingular.
    IF( UPPER ) THEN
    ! Upper triangular storage: examine D from bottom to top
        DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
        10 CONTINUE
    ELSE
    ! Lower triangular storage: examine D from top to bottom.
        DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
        20 CONTINUE
    END IF
    INFO = 0
    IF( UPPER ) THEN
    ! Compute inv(A) from the factorization A = U*D*U'.
    ! K is the main loop index, increasing from 1 to N in steps of
    ! 1 or 2, depending on the size of the diagonal blocks.
    K = 1
    30 CONTINUE
    ! If K > N, exit from loop.
    IF( K.GT.N ) GO TO 40
    IF( IPIV( K ).GT.0 ) THEN
        ! 1 x 1 diagonal block
        ! Invert the diagonal block.
        A( K, K ) = ONE / A( K, K )
        ! Compute column K of the inverse.
        IF( K.GT.1 ) THEN
            CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
            CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )
            A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ), 1 )
        END IF
        KSTEP = 1
    ELSE
        ! 2 x 2 diagonal block
        ! Invert the diagonal block.
        T = ABS( A( K, K+1 ) )
        AK = A( K, K ) / T
        AKP1 = A( K+1, K+1 ) / T
        AKKP1 = A( K, K+1 ) / T
        D = T*( AK*AKP1-ONE )
        A( K, K ) = AKP1 / D
        A( K+1, K+1 ) = AK / D
        A( K, K+1 ) = -AKKP1 / D
        ! Compute columns K and K+1 of the inverse.
        IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ), 1 )
               A( K, K+1 ) = A( K, K+1 ) - DDOT( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL DCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K+1 ), 1 )
               A( K+1, K+1 ) = A( K+1, K+1 ) - DDOT( K-1, WORK, 1, A( 1, K+1 ), 1 )
        END IF
            KSTEP = 2
        END IF

        KP = ABS( IPIV( K ) )
        IF( KP.NE.K ) THEN
            ! Interchange rows and columns K and KP in the leading
            ! submatrix A(1:k+1,1:k+1)
            CALL DSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL DSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF
        END IF
        K = K + KSTEP
        GO TO 30
    40 CONTINUE
    ELSE
        ! Compute inv(A) from the factorization A = L*D*L'.
        ! K is the main loop index, increasing from 1 to N in steps of
        ! 1 or 2, depending on the size of the diagonal blocks.
        K = N
        50 CONTINUE
        ! If K < 1, exit from loop.
        IF( K.LT.1 ) GO TO 60
        IF( IPIV( K ).GT.0 ) THEN
            ! 1 x 1 diagonal block
            ! Invert the diagonal block.
            A( K, K ) = ONE / A( K, K )
            ! Compute column K of the inverse.
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ), 1 )
            END IF
            KSTEP = 1
        ELSE
            ! 2 x 2 diagonal block
            ! Invert the diagonal block.
            T = ABS( A( K, K-1 ) )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
            ! Compute columns K-1 and K of the inverse.
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ), 1 )
               A( K, K-1 ) = A( K, K-1 ) - DDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 )
               CALL DCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K-1 ), 1 )
               A( K-1, K-1 ) = A( K-1, K-1 ) - DDOT( N-K, WORK, 1, A( K+1, K-1 ), 1 )
            END IF
            KSTEP = 2
        END IF
        KP = ABS( IPIV( K ) )
        IF( KP.NE.K ) THEN
            !  Interchange rows and columns K and KP in the trailing
            ! submatrix A(k-1:n,k-1:n)
            IF( KP.LT.N ) CALL DSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL DSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
                TEMP = A( K, K-1 )
                A( K, K-1 ) = A( KP, K-1 )
                A( KP, K-1 ) = TEMP
            END IF
        END IF
        K = K - KSTEP
        GO TO 50
        60 CONTINUE
    END IF
    RETURN
end subroutine DSYTRI

subroutine DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!----------------------------------------------------------------
!Purpose:
!   DSYMV  performs the matrix-vector  operation
!      y := alpha*A*x + beta*y,
!   where alpha and beta are scalars, x and y are n element vectors and
!   A is an n by n symmetric matrix.
! 
!Arguments:
!   UPLO   - CHARACTER*1.
!            On entry, UPLO specifies whether the upper or lower
!            triangular part of the array A is to be referenced as
!            follows:
! 
!               UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                   is to be referenced.
! 
!               UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                   is to be referenced.
! 
!            Unchanged on exit.
! 
!   N      - INTEGER.
!            On entry, N specifies the order of the matrix A.
!            N must be at least zero.
!            Unchanged on exit.
! 
!   ALPHA  - DOUBLE PRECISION.
!            On entry, ALPHA specifies the scalar alpha.
!            Unchanged on exit.
! 
!   A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!            Before entry with  UPLO = 'U' or 'u', the leading n by n
!            upper triangular part of the array A must contain the upper
!            triangular part of the symmetric matrix and the strictly
!            lower triangular part of A is not referenced.
!            Before entry with UPLO = 'L' or 'l', the leading n by n
!            lower triangular part of the array A must contain the lower
!            triangular part of the symmetric matrix and the strictly
!            upper triangular part of A is not referenced.
!            Unchanged on exit.
! 
!   LDA    - INTEGER.
!            On entry, LDA specifies the first dimension of A as declared
!            in the calling (sub) program. LDA must be at least
!            max( 1, n ).
!            Unchanged on exit.
! 
!   X      - DOUBLE PRECISION array of dimension at least
!            ( 1 + ( n - 1 )*abs( INCX ) ).
!            Before entry, the incremented array X must contain the n
!            element vector x.
!            Unchanged on exit.
! 
!   INCX   - INTEGER.
!            On entry, INCX specifies the increment for the elements of
!            X. INCX must not be zero.
!            Unchanged on exit.
! 
!   BETA   - DOUBLE PRECISION.
!            On entry, BETA specifies the scalar beta. When BETA is
!            supplied as zero then Y need not be set on input.
!            Unchanged on exit.
! 
!   Y      - DOUBLE PRECISION array of dimension at least
!            ( 1 + ( n - 1 )*abs( INCY ) ).
!            Before entry, the incremented array Y must contain the n
!            element vector y. On exit, Y is overwritten by the updated
!            vector y.
! 
!   INCY   - INTEGER.
!            On entry, INCY specifies the increment for the elements of
!            Y. INCY must not be zero.
!            Unchanged on exit.
! 
! 
!   Level 2 Blas routine.
! 
!   -- Written on 22-October-1986.
!      Jack Dongarra, Argonne National Lab.
!      Jeremy Du Croz, Nag Central Office.
!      Sven Hammarling, Nag Central Office.
!      Richard Hanson, Sandia National Labs.
!----------------------------------------------------------------
    ! .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA,BETA
    INTEGER INCX,INCY,LDA,N
    CHARACTER UPLO
    ! .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*),Y(*)

    ! .. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    ! .. Local Scalars ..
    DOUBLE PRECISION TEMP1,TEMP2
    INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
    ! .. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    ! .. External Subroutines ..
    !EXTERNAL XERBLA
    ! .. Intrinsic Functions ..
    INTRINSIC MAX

    ! Test the input parameters.
    INFO = 0
    IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
        INFO = 1
    ELSE IF (N.LT.0) THEN
        INFO = 2
    ELSE IF (LDA.LT.MAX(1,N)) THEN
        INFO = 5
    ELSE IF (INCX.EQ.0) THEN
        INFO = 7
    ELSE IF (INCY.EQ.0) THEN
        INFO = 10
    END IF
    IF (INFO.NE.0) THEN
        CALL XERBLA('DSYMV ',INFO)
        RETURN
    END IF
    ! Quick return if possible.
    IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
    ! Set up the start points in  X  and  Y.
    IF (INCX.GT.0) THEN
        KX = 1
    ELSE
        KX = 1 - (N-1)*INCX
    END IF
    IF (INCY.GT.0) THEN
        KY = 1
    ELSE
        KY = 1 - (N-1)*INCY
    END IF

    ! Start the operations. In this version the elements of A are
    ! accessed sequentially with one pass through the triangular part
    ! of A.

    ! First form  y := beta*y.
    IF (BETA.NE.ONE) THEN
        IF (INCY.EQ.1) THEN
            IF (BETA.EQ.ZERO) THEN
                DO 10 I = 1,N
                    Y(I) = ZERO
                10 CONTINUE
            ELSE
                DO 20 I = 1,N
                    Y(I) = BETA*Y(I)
                20 CONTINUE
            END IF
        ELSE
            IY = KY
            IF (BETA.EQ.ZERO) THEN
                DO 30 I = 1,N
                    Y(IY) = ZERO
                    IY = IY + INCY
                30 CONTINUE
            ELSE
                DO 40 I = 1,N
                    Y(IY) = BETA*Y(IY)
                    IY = IY + INCY
                40 CONTINUE
            END IF
        END IF
    END IF
    IF (ALPHA.EQ.ZERO) RETURN
    IF (LSAME(UPLO,'U')) THEN
    ! Form  y  when A is stored in upper triangle.
        IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
            DO 60 J = 1,N
                TEMP1 = ALPHA*X(J)
                TEMP2 = ZERO
                DO 50 I = 1,J - 1
                    Y(I) = Y(I) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(I)
                50 CONTINUE
                Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
            60 CONTINUE
        ELSE
            JX = KX
            JY = KY
            DO 80 J = 1,N
                TEMP1 = ALPHA*X(JX)
                TEMP2 = ZERO
                IX = KX
                IY = KY
                DO 70 I = 1,J - 1
                    Y(IY) = Y(IY) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(IX)
                    IX = IX + INCX
                    IY = IY + INCY
                70 CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
            80 CONTINUE
        END IF
    ELSE
        ! Form  y  when A is stored in lower triangle.
        IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
            DO 100 J = 1,N
                TEMP1 = ALPHA*X(J)
                TEMP2 = ZERO
                Y(J) = Y(J) + TEMP1*A(J,J)
                DO 90 I = J + 1,N
                    Y(I) = Y(I) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(I)
                90 CONTINUE
                Y(J) = Y(J) + ALPHA*TEMP2
            100 CONTINUE
        ELSE
            JX = KX
            JY = KY
            DO 120 J = 1,N
                TEMP1 = ALPHA*X(JX)
                TEMP2 = ZERO
                Y(JY) = Y(JY) + TEMP1*A(J,J)
                IX = JX
                IY = JY
                DO 110 I = J + 1,N
                    IX = IX + INCX
                    IY = IY + INCY
                    Y(IY) = Y(IY) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(IX)
                110 CONTINUE
                Y(JY) = Y(JY) + ALPHA*TEMP2
                JX = JX + INCX
                JY = JY + INCY
            120 CONTINUE
        END IF
    END IF
    RETURN
end subroutine DSYMV

function IEEECK( ISPEC, ZERO, ONE )
!----------------------------------------------------------------
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!Purpose:
!    IEEECK is called from the ILAENV to verify that Infinity and
!    possibly NaN arithmetic is safe (i.e. will not trap).
!  
!Arguments:
!    ISPEC   (input) INTEGER
!            Specifies whether to test just for inifinity arithmetic
!            or whether to test for infinity and NaN arithmetic.
!            = 0: Verify infinity arithmetic only.
!            = 1: Verify infinity and NaN arithmetic.
!  
!    ZERO    (input) REAL
!            Must contain the value 0.0
!            This is passed to prevent the compiler from optimizing
!            away this code.
!  
!    ONE     (input) REAL
!            Must contain the value 1.0
!            This is passed to prevent the compiler from optimizing
!            away this code.
!  
!    RETURN VALUE:  INTEGER
!            = 0:  Arithmetic failed to produce the correct answers
!            = 1:  Arithmetic produced the correct answers
!----------------------------------------------------------------
    INTEGER IEEECK
    ! .. Scalar Arguments ..
    INTEGER            ISPEC
    REAL               ONE, ZERO 
    !.. Local Scalars ..
    REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, NEGZRO, NEWZRO, POSINF

    !.. Executable Statements ..
    IEEECK = 1

    POSINF = ONE / ZERO
    IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
    END IF

    NEGINF = -ONE / ZERO
    IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
    END IF

    NEGZRO = ONE / ( NEGINF+ONE )
    IF( NEGZRO.NE.ZERO ) THEN
        IEEECK = 0
        RETURN
    END IF

    NEGINF = ONE / NEGZRO
    IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
    END IF

    NEWZRO = NEGZRO + ZERO
    IF( NEWZRO.NE.ZERO ) THEN
        IEEECK = 0
        RETURN
    END IF

    POSINF = ONE / NEWZRO
    IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
    END IF

    NEGINF = NEGINF*POSINF
    IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
    END IF

    POSINF = POSINF*POSINF
    IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
    END IF


    ! Return if we were only asked to check infinity arithmetic
    IF( ISPEC.EQ.0 ) RETURN

    NAN1 = POSINF + NEGINF
    NAN2 = POSINF / NEGINF
    NAN3 = POSINF / POSINF
    NAN4 = POSINF*ZERO
    NAN5 = NEGINF*NEGZRO
    NAN6 = NAN5*0.0

    IF( NAN1.EQ.NAN1 ) THEN
        IEEECK = 0
        RETURN
    END IF

    IF( NAN2.EQ.NAN2 ) THEN
        IEEECK = 0
        RETURN
    END IF

    IF( NAN3.EQ.NAN3 ) THEN
        IEEECK = 0
        RETURN
    END IF

    IF( NAN4.EQ.NAN4 ) THEN
        IEEECK = 0
        RETURN
    END IF

    IF( NAN5.EQ.NAN5 ) THEN
        IEEECK = 0
        RETURN
    END IF

    IF( NAN6.EQ.NAN6 ) THEN
        IEEECK = 0
        RETURN
    END IF

    RETURN
end function IEEECK

function IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!----------------------------------------------------------------
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!Purpose:
!       This program sets problem and machine dependent parameters
!       useful for xHSEQR and its subroutines. It is called whenever 
!       ILAENV is called with 12 <= ISPEC <= 16
!
!Arguments:
!   ISPEC  (input) integer scalar
!           ISPEC specifies which tunable parameter IPARMQ should
!           return.
!
!           ISPEC=12: (INMIN)  Matrices of order nmin or less
!                   are sent directly to xLAHQR, the implicit
!                   double shift QR algorithm.  NMIN must be
!                   at least 11.
!
!           ISPEC=13: (INWIN)  Size of the deflation window.
!                   This is best set greater than or equal to
!                   the number of simultaneous shifts NS.
!                   Larger matrices benefit from larger deflation
!                   windows.
!
!           ISPEC=14: (INIBL) Determines when to stop nibbling and
!                   invest in an (expensive) multi-shift QR sweep.
!                   If the aggressive early deflation subroutine
!                   finds LD converged eigenvalues from an order
!                   NW deflation window and LD.GT.(NW*NIBBLE)/100,
!                   then the next QR sweep is skipped and early
!                   deflation is applied immediately to the
!                   remaining active diagonal block.  Setting
!                   IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!                   multi-shift QR sweep whenever early deflation
!                   finds a converged eigenvalue.  Setting
!                   IPARMQ(ISPEC=14) greater than or equal to 100
!                   prevents TTQRE from skipping a multi-shift
!                   QR sweep.
!
!           ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!                   a multi-shift QR iteration.
!
!           ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!                   following meanings.
!                   0:  During the multi-shift QR sweep,
!                       xLAQR5 does not accumulate reflections and
!                       does not use matrix-matrix multiply to
!                       update the far-from-diagonal matrix
!                       entries.
!                   1:  During the multi-shift QR sweep,
!                       xLAQR5 and/or xLAQRaccumulates reflections and uses
!                       matrix-matrix multiply to update the
!                       far-from-diagonal matrix entries.
!                   2:  During the multi-shift QR sweep.
!                       xLAQR5 accumulates reflections and takes
!                       advantage of 2-by-2 block structure during
!                       matrix-matrix multiplies.
!                   (If xTRMM is slower than xGEMM, then
!                   IPARMQ(ISPEC=16)=1 may be more efficient than
!                   IPARMQ(ISPEC=16)=2 despite the greater level of
!                   arithmetic work implied by the latter choice.)
!
!   NAME    (input) character string
!           Name of the calling subroutine
!
!   OPTS    (input) character string
!           This is a concatenation of the string arguments to
!           TTQRE.
!
!   N       (input) integer scalar
!           N is the order of the Hessenberg matrix H.
!
!   ILO     (input) INTEGER
!   IHI     (input) INTEGER
!           It is assumed that H is already upper triangular
!           in rows and columns 1:ILO-1 and IHI+1:N.
!
!   LWORK   (input) integer scalar
!           The amount of workspace available.
!
! Further Details:
!      Little is known about how best to choose these parameters.
!      It is possible to use different values of the parameters
!      for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!
!      It is probably best to choose different parameters for
!      different matrices and different parameters at different
!      times during the iteration, but this has not been
!      implemented --- yet.
!
!
!      The best choices of most of the parameters depend
!      in an ill-understood way on the relative execution
!      rate of xLAQR3 and xLAQR5 and on the nature of each
!      particular eigenvalue problem.  Experiment may be the
!      only practical way to determine which choices are most
!      effective.
!
!      Following is a list of default values supplied by IPARMQ.
!      These defaults may be adjusted in order to attain better
!      performance in any particular computational environment.
!
!      IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!                       Default: 75. (Must be at least 11.)
!
!      IPARMQ(ISPEC=13) Recommended deflation window size.
!                       This depends on ILO, IHI and NS, the
!                       number of simultaneous shifts returned
!                       by IPARMQ(ISPEC=15).  The default for
!                       (IHI-ILO+1).LE.500 is NS.  The default
!                       for (IHI-ILO+1).GT.500 is 3*NS/2.
!
!      IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!
!      IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!                       a multi-shift QR iteration.
!
!                       If IHI-ILO+1 is ...
!
!                       greater than      ...but less    ... the
!                       or equal to ...      than        default is
!
!                               0               30       NS =   2+
!                              30               60       NS =   4+
!                              60              150       NS =  10
!                             150              590       NS =  **
!                             590             3000       NS =  64
!                            3000             6000       NS = 128
!                            6000             infinity   NS = 256
!
!                   (+)  By default matrices of this order are
!                        passed to the implicit double shift routine
!                        xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!                        values of NS are used only in case of a rare
!                        xLAHQR failure.
!
!                   (**) The asterisks (**) indicate an ad-hoc
!                        function increasing from 10 to 64.
!
!      IPARMQ(ISPEC=16) Select structured matrix multiply.
!                       (See ISPEC=16 above for details.)
!                       Default: 3.
!----------------------------------------------------------------
    INTEGER IPARMQ
    ! .. Scalar Arguments ..
    INTEGER            IHI, ILO, ISPEC, LWORK, N
    CHARACTER          NAME*( * ), OPTS*( * )
    ! .. Parameters ..
    INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
    PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, ISHFTS = 15, IACC22 = 16 )
    INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
    PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500 )
    REAL               TWO
    PARAMETER          ( TWO = 2.0 )
    ! .. Local Scalars ..
    INTEGER            NH, NS
    !.. Intrinsic Functions ..
    INTRINSIC          LOG, MAX, MOD, NINT, REAL

    ! .. Executable Statements ..
    IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. ( ISPEC.EQ.IACC22 ) ) THEN
        ! Set the number simultaneous shifts
        NH = IHI - ILO + 1
        NS = 2
        IF( NH.GE.30 ) NS = 4
        IF( NH.GE.60 ) NS = 10
        IF( NH.GE.150 ) NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
        IF( NH.GE.590 ) NS = 64
        IF( NH.GE.3000 ) NS = 128
        IF( NH.GE.6000 ) NS = 256
        NS = MAX( 2, NS-MOD( NS, 2 ) )
    END IF

    IF( ISPEC.EQ.INMIN ) THEN
        ! Matrices of order smaller than NMIN get sent to xLAHQR, the classic double shift algorithm.
        ! This must be at least 11.
        IPARMQ = NMIN
    ELSE IF( ISPEC.EQ.INIBL ) THEN
        ! INIBL: skip a multi-shift qr iteration and whenever aggressive early deflation finds
        ! at least (NIBBLE*(window size)/100) deflations. 
        IPARMQ = NIBBLE
    ELSE IF( ISPEC.EQ.ISHFTS ) THEN
        ! NSHFTS: The number of simultaneous shifts =====
        IPARMQ = NS
    ELSE IF( ISPEC.EQ.INWIN ) THEN
        ! NW: deflation window size. 
        IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
        ELSE
            IPARMQ = 3*NS / 2
        END IF
    ELSE IF( ISPEC.EQ.IACC22 ) THEN
        ! IACC22: Whether to accumulate reflections before updating the far-from-diagonal elements
        ! and whether to use 2-by-2 block structure while doing it.  A small amount of work could be saved
        ! by making this choice dependent also upon the NH=IHI-ILO+1.
        IPARMQ = 0
        IF( NS.GE.KACMIN ) IPARMQ = 1
        IF( NS.GE.K22MIN ) IPARMQ = 2
    ELSE
        ! invalid value of ispec
        IPARMQ = -1
    END IF
end function IPARMQ

subroutine clingd(ma,mx,n,m,a,x,d,ifl)      
    !----------------------------------------------------------------------       
    !     solves the system A*X=B, where B is at the beginning on X
    !     it will be overwritten later on, d is the determinant of A    
    !----------------------------------------------------------------------
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    COMPLEX*16 A,X,D,CP,CQ
    DIMENSION  A(MA,N),X(MX,M)
    character(len=4) :: CC
    CC = '----'
    TOLLIM = 1.E-16
    IFL=1
    P=0.D0
    DO 10 I=1,N
        Q=0.D0
        DO 20 J=1,N
            Q=Q+CDABS(A(I,J))
        20 CONTINUE
        IF (Q.GT.P) P=Q
    10   CONTINUE
    TOL=TOLLIM*P
    D=CMPLX(1.D0,0.D0)
    DO 30 K=1,N
        P=0.D0
        DO 40 J=K,N
            Q=CDABS(A(J,K))
            IF (Q.LT.P) GOTO 40
            P=Q
            I=J
        40   CONTINUE
        IF (P.GT.TOL) GOTO 70
        WRITE (*,200) (CC,J=1,22),TOL,I,K,A(I,K),(CC,J=1,22)
        200 FORMAT (1X,22A4/' *****  ERROR IN LINGM, TOLERANZ =',E11.4,' WERT VON A(',I3,',',I3,') IST ',E11.4/1X,22A4)
        IFL=-1
        RETURN
        70 CP=1./A(I,K)
        IF (I.EQ.K) GOTO 90
        D=-D
        DO 81 L=1,M
            CQ=X(I,L)
            X(I,L)=X(K,L)
            X(K,L)=CQ
        81 CONTINUE
        DO 80 L=K,N
            CQ=A(I,L)
            A(I,L)=A(K,L)
            A(K,L)=CQ
        80 CONTINUE
        90 D=D*A(K,K)
        IF (K.EQ.N) GOTO 1
        K1=K+1

        DO I=K1,N
            CQ=A(I,K)*CP
            DO L=1,M
                X(I,L)=X(I,L)-CQ*X(K,L)
            END DO
            DO L=K1,N
                A(I,L)=A(I,L)-CQ*A(K,L)
            END DO
        END DO
    30 CONTINUE  
    1 DO 126 L=1,M
        X(N,L)=X(N,L)*CP
    126 CONTINUE
    IF (N.EQ.1) RETURN
    N1=N-1

    DO K=1,N1
        CP=1./A(N-K,N-K)
        DO L=1,M
            CQ=X(N-K,L)
            DO I=1,K
                CQ=CQ-A(N-K,N+1-I)*X(N+1-I,L)
            END DO
            X(N-K,L)=CQ*CP
        END DO
    END DO
    return
end 

subroutine cmatmul_ABC(A, B, C, D, m, k, p, n)
    !-------------------------------------------------------------
    !   Compute the product of three complex matrices:
    !       D = A * B * C
    !
    ! Inputs:
    !   A(m, k) - complex*16, first matrix
    !   B(k, p) - complex*16, second matrix
    !   C(p, n) - complex*16, third matrix
    !
    ! Output:
    !   D(m, n) - complex*16, the result matrix
    !-------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, k, p, n
    complex*16, intent(in) :: A(m, k)
    complex*16, intent(in) :: B(k, p)
    complex*16, intent(in) :: C(p, n)
    complex*16, intent(out) :: D(m, n)

    complex*16, allocatable :: Temp(:,:)
    integer :: i, j, l

    allocate(Temp(m, p))
    Temp = (0.0d0, 0.0d0)

    ! Compute Temp = A * B
    do i = 1, m
        do j = 1, p
            do l = 1, k
                Temp(i, j) = Temp(i, j) + A(i, l) * B(l, j)
            end do
        end do
    end do

    D = (0.0d0, 0.0d0)

    ! Compute D = Temp * C
    do i = 1, m
        do j = 1, n
            do l = 1, p
                D(i, j) = D(i, j) + Temp(i, l) * C(l, j)
            end do
        end do
    end do

    deallocate(Temp)

end subroutine cmatmul_ABC
    
Subroutine ZPfaffianF(SK,LDS,N,Ipiv,Pf)
    !-------------------------------------------------------------------------
    !     subroutine to calculate the Pfaffian of a skew-symmetric matrix from 
    !     C. Gonzalez-Ballestero, L. M. Robledo and G. F. Bertsch,
    !     Computer Physics Communications (2011)         
    !==========================================================================
    !+---------------------------------------------------------------------+
    !   Computes the Pfaffian of a skew-symmetric matrix SK using         |
    !   the ideas of Aitken's block diagonalization formula.              |
    !   Full pivoting is implemented to make the whole procedure stable   |
    !+---------------------------------------------------------------------+
    !   Only the upper triangle of SK is used                             |
    !   In output the lower triangle contains the transformation matrix   | 
    !+---------------------------------------------------------------------+
    !   SK .....Skew symmetric input matrix                               |
    !   LDS ....Leading dimension of matrix SK                            |
    !   N  .....Dimension of SK                                           |
    !   Ipiv ...Integer vector to hold the pivoted columns                |
    !+---------------------------------------------------------------------+
    ! This program may be freely used providing such use cites the source,|
    !      Numeric and symbolic evaluation of the pfaffian of general     |
    !      skew-symmetric matrices                                        |
    !                                                                     |
    !      C. Gonzalez-Ballestero, L.M.Robledo and G.F. Bertsch           |
    !      Computer Physics Communications (2011)                         |
    !+---------------------------------------------------------------------+
    Implicit none
    Double Complex SK(LDS,N)
    Double Complex Pf,SS,SW,one,zero
    Double Precision epsln,phas,big
    Integer Ipiv(N,2),N,LDS,NB,IB,i,j,k,ip,NR,NC,I1,I2
    one = dcmplx(1.0d+00,0.0d+00)
    zero= dcmplx(0.0d+00,0.0d+00)
    epsln = 1.0d-13   ! smallest number such that 1+epsln=1
    if(mod(N,2).eq.1) then ! N odd
        Pf = zero
    else
        NB= N / 2   ! Number of 2x2 blocks
        if (NB.eq.1) then 
        ! The pfaffian of a 2x2 matrix is Pf=Sk(1,2)
        Pf = SK(1,2)
        else   ! NB.gt.1
            Pf = one
            do IB=1,NB-1
                NR = IB*2 - 1 ! row numb of the 1,2 element of the 2x2 block
                NC = NR + 1
                big = 0.0d+00
                I1 = NR
                I2 = NC
                do i=NR,N-1 ! all rows
                    do j=i+1,N
                        if(abs(SK( i , j )).gt. big) then
                            big = abs(SK( i , j ))
                            I1 = i
                            I2 = j 
                        end if
                    end do
                end do
                Ipiv( IB , 1) = I1 ! to initialize
                Ipiv( IB , 2) = I2 ! to initialize
                phas = 1.0d+00
                if(I1.eq.NR) phas = -phas
                if(I2.eq.NC) phas = -phas 
                ! Pivoting of element NR,NC with ip1,ip2
                ! 
                ! Pivoting for a skew-symmetric matrix (Upper)
                !
                if(I1.ne.NR) Call Zexch(SK,LDS,N,NR,I1)
                if(I2.ne.NC) Call Zexch(SK,LDS,N,NC,I2)
                ss = Sk( NR , NC )
                if(abs(ss).gt.epsln) then
                    !  Updating the Schur complement  matrix
                    do i = 2*IB+1 , N-1
                        do j = i+1 , N
                            Sk(i,j) = Sk(i,j) + (Sk(NC,i)*Sk(NR,j)-Sk(NR,i)*Sk(NC,j))/SS
                        end do
                    end do
                    ! Storing X and Y vectors in the lower part of the matrix    
                    do i = 2*IB+1 , N
                        Sk( i , NR ) =  -Sk( NC , i )/SS
                        Sk( i , NC ) =   Sk( NR , i )/SS
                    end do
                    if (IB.ge.2) then  ! swap
                        ! Swapping      
                        do j= 1 , 2*IB
                            SW           = Sk( NR , j )
                            Sk( NR , j ) = Sk( I1 , j ) 
                            Sk( I1 , j ) = SW
                            SW           = Sk( NC , j )
                            Sk( NC , j ) = Sk( I2 , j ) 
                            Sk( I2 , j ) = SW
                        end do       
                    end if
                else
                    Pf = zero
                    return
                end if ! dabs(ss).gt.epsln
                Pf = Pf * SS * phas
            end do ! IB
            Pf = Pf * Sk(N-1,N)
        end if !    NB.eq.1
    end if !    mod(N,2).eq.1
    return
end

Subroutine Zexch(SK,LDS,N,I1,I2)
    !  Exchange of rows i1 i2 and columns i1 i2
    !  SK is assumed upper triangular
    !  i1 < i2
    Implicit none
    Double Complex SK(LDS,N)
    Double Complex SS
    Integer LDS,N,I1,I2,I
    SK(i1,i2) = -SK(i1,i2)
    if(I1.ne.1) then
        do i=1,I1-1
            SS           = SK( i , I1 )
            SK( i , I1 ) = SK( i , I2 )
            SK( i , I2 ) = SS
        end do
    end if
    if(I2.ne.N) then
        do i=I2+1,N
            SS           = SK( I1 , i )
            SK( I1 , i ) = SK( I2 , i )
            SK( I2 , i ) = SS
        end do
    end if
    if(I2.ge.I1+2) then
        do i=I1+1,I2-1
            SS           =  SK( I1 , i )
            SK( I1 , i ) = -SK( i , I2 )
            SK( i , I2 ) = -SS
        end do
    end if
    return
end
END MODULE MathMethods