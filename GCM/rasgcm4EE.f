c_______________________________________________________________________
       program rasgcm 

c     Version v08
c     - Kernels can be calculated with analytic expressions
c     Version v06
c     - neutron part of matrix element <Q_L> is also calculate.
c     - gcmr2(iq) is set to zero all the time.
c     Version v05
c     - name of the element file is renamed.
c     Version v04
c     jmy(2015/2/28)
c     Version v03
c     jmy(2015/2/9)
c     - subroutine readkernel() is changed by including E3 matrix elements
c    
c     Version v01
c     jmy(2015/1/26)    
c     -configuration mixing of parity, particle number and angular
c      momentum projected axial-octupole deformed states 
!-----------------------------------------------------------------------  
      implicit real*8 (a-h,o-z)
      include 'dic.par'

      parameter (mgcm = 5)
      parameter (jsize   = jmax/2+1)
      parameter (jdimax2 = 2*jmax+1)
!
      logical lpr,iexist      
      character tp*1,tis*1,tit*8,tl*1                           ! textex 
      character*2 nucnam                                        ! nucnuc  
      complex*16 cnum,cden 
      character*58 olpf,elem,dens
      character*28 outputwf,outputwf1,outputwf2                 ! ttttt1
      dimension iexst(maxq*jmax*maxq*jmax,maxq*jmax*maxq*jmax)
      character(500) find_file
c
      common /ttttt4/ outputwf,outputwf1,outputwf2
      common /textex/ tp(2),tis(2),tit(2),tl(0:30) 
      common /mathco/ zero,one,two,half,third,pi  
      common /gfviv / iv(-igfv:igfv)
      common /gfvsq / sq(0:igfv)  
      common /physco/ hbc,alphi,r0
      common /nucnuc/ amas,nama,nneu,npro,nucnam 
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gcmcut/ zeta(5),jmx,icase,kmx 
      common /gcmsta/ nmaxdi(0:jmax),jkmax(0:jmax)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp,lrpa 
      common /eulera/ nbet,nphi 
      common /pkcout/ nubd,iden,nubw
      common /trdens/ transrho(ngr,mgcm,mgcm,jsize,jsize,jdimax2),
     &                scalerho(ngr,mgcm,mgcm,jsize,jsize,jdimax2)
c..............................................................................
   99 format (/,'  __________________________________________________ ',
     1        /,' |                                                  |',
     2        /,' |  program rasgcm4EE for ee nuclei                 |',
     3        /,' |                                                  |',
     4        /,' |  23 Jan  2025                                    |',
     5        /,' |  Copyright NTG@SYSU                              |',
     6        /,' |                                                  |',
     6        /,' |  If you have found any problem in the code       |',
     6        /,' |                                                  |',
     9        /,' |  please contact us by sending email.             |',
     9        /,' |__________________________________________________|',
     9        /)
  100 format ('      No Projections (No AMP & No PNP)  ')
  101 format(a,i2)
  102 format(/,1x,74('-'),/)
  103 format ('      Pure Particle Number   Projection (PNP)  ')
  104 format ('      Pure Angular  Momentum Projection (AMP)  ')
  200 format ('      Both PN and AM projection (PNAMP)  ')
  400 format(/,' existing matrix elements',
     1         ' -- value 99 means missing, otherwise ln(overlap)',/)
  401 format(2x,2i3,1x,55i3)
  402 format ('    K  iq')
  403 format ('                                  ')
      write(*,99)
c     ........................................................ default value
      call default() 
c     ........................................................ read parameters
      call reader()
      if(nbet.eq.1.and.nphi.eq.1) write(*,100)               
      if(nbet.eq.1.and.nphi.gt.1) write(*,103)               
      if(nbet.gt.1.and.nphi.eq.1) write(*,104)               
      if(nbet.gt.1.and.nphi.gt.1) write(*,200)               
      write(*,102)
c     ........................................................ initialization
      jmn  = 0              
      jdf  = 1
      if(icase.eq.0) then
        name = 1+48 
      else if(icase.eq.1) then
        name = 3+48 
      else
        stop ' icase should be either 0 or 1'
      endif
!     determination of the location of the elements
        WFS_DIR = find_file("GCM_FILES_DIR",'HFB.wf')
        !print *, 'WFS_DIR=',WFS_DIR
        inquire (file=trim(WFS_DIR)//'HFB.wfs',exist=iexist)
        if(.not.iexist) stop ' pls create the HFB.wf in the dir...'
 
!     ........................................................ initialization
      do i=1,maxq*jmax*maxq*jmax 
         iexst(i,1) = 99
      enddo
      monopole = 0   ! 1: with constraint on monopole; 0: no
!     ......................................................... loop over q1
      do iq1=1, maxmp                           
          betac1  = gcmbet(iq1) 
          gammac1 = gcmgam(iq1)
          r2c1    = gcmr2(iq1)
          iq2m    = maxmp  
          if(nbet.eq.1.and.monopole.eq.0) iq2m = iq1               
!     ......................................................... loop over q2
          do iq2 =iq1, iq2m        
             betac2  = gcmbet(iq2) 
             gammac2 = gcmgam(iq2)  
             r2c2    = 0.d0 ! gcmr2(iq2)
            call re_name(betac1,gammac1,r2c1,betac2,gammac2,r2c2,name,elem,dens) 
c              write(*,*) elem
              ! print *, 'file=',trim(WFS_DIR)//elem
              open(lou,file=trim(WFS_DIR)//elem,err=900)
  900         call readkernel(icase,iq1,iq2,jmn,jmx,jdf,elem) 
              call kernels(jmn,jmx,jdf,iq1,iq2,icase,nmaxdi,.false.) 
              close(lou)
      enddo !iq2    
      enddo !iq1
      write(*,*) 'Finished reading and calculate kernels.'
c    ......................................................... loop over angular momentum J
      do jproj = jmn,jmx,jdf 
          iki0 = 0
c         if(jproj.eq.1) goto 89
c         if(iv(jproj).gt.zero)  iki0    = 0  
c         if(iv(jproj).lt.zero)  iki0    = 2 
          kmax = icase*jproj   
!------------------------------------------------------------------------    
              do k1=iki0,kmax,2   ! K-value
              do k2=iki0,kmax,2 
                 k1m = k1-iki0   ! mesh point
                 k2m = k2-iki0
c    ......................................................... print the diagonal element 
                 call table(k1,k2,k1m,k2m,jproj,iexst)
              enddo
              enddo  
      
c    ......................................................... print the norm overlap 
      print 400
      print 402
c     print 403
      do k1=iki0,kmax,2   ! K-value
         k1m = k1-iki0    ! mesh point
      do iq1 = 1,maxmp
         iqk1=iq1+k1m/2*maxmp
         write(*,401) k1,iq1,(iexst(iqk1,iqk2),iqk2=1,iqk1)
      enddo
      enddo  
      nmaxdj = nmaxdi(jproj)  
c     ........................................................... solution of HWG equation
      call hweq(jproj,nmaxdj,maxdj,.false.) 
c      call hweq(jproj,nmaxdj,maxdj,.true.) 
      jkmax(jproj) = maxdj 
      write(*,102)
!------------------------------------------------------------------------  
   89 enddo !jproj  
c     ............................................................ calculate observables
      write(*,*) ' Startprint  calculate observables '
      call obsers(jmn,jmx,jdf) 
c     ............................................................ print spectrum  
      write(*,*) ' Startprint spectrum '
      call spect(jmn,jdf)
c     .................................................... print table with trans. strengths
      call trans(jmn,jdf) 
!     ................................................. projected density
      if(iden.gt.0.and.icase.eq.0) then
!     ..................................................initialization
      do i=1,ngr*mgcm*mgcm*jsize*jsize*jdimax2
         transrho(i,1,1,1,1,1) = zero
         scalerho(i,1,1,1,1,1) = zero
      enddo ! i
      do iq1=1, maxmp
          betac1  = gcmbet(iq1)
          gammac1 = gcmgam(iq1)
          r2c1    = gcmr2(iq1)
          do iq2=1, maxmp
             betac2  = gcmbet(iq2)
             gammac2 = gcmgam(iq2)
             r2c2    = gcmr2(iq2)
             call re_name(betac1,gammac1,r2c1,betac2,gammac2,r2c2, 
     &                   name,elem,dens) 
              write(22,*) dens 
             open(nubd,file=dens,err=800)
c             call readen(jmn,jdf)
c             call caden(iq1,iq2,jmn,jdf)
          enddo ! iq2
      enddo ! iq1
c      call hwdens(jmn,jdf) 
 800  endif  ! iden
c     .................................................... end
      stop 'alright' 
      END  

c______________________________________________________________________________
      subroutine table(k1,k2,k1m,k2m,jproj,iexst)

c..............................................................................
c     This subroutine prints out a table for the diagonal elements. 
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'dic.par'
      character*2 nucnam                                        ! nucnuc  
      complex*16 njqqkk,hjqqkk,xnkk 
      complex*16 qp,qcp,qpred,qcpred,q0p,q0pc
      complex*16 q3p,q3cp,q3pred,q3cpred,q1pred,q1cpred
      complex*16 q3n,q3cn,q3nred,q3cnred
      complex*16 qn,qcn,qnred,qcnred,q0n,q0nc
      dimension iexst(maxq*jmax,maxq*jmax)


      common /nucnuc/ amas,nama,nneu,npro,nucnam  
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax)
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp

      common /big   / qp    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcp   (maxq,maxq,0:jmax,0:jmax,0:jmax), 
     2                qpred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0p   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0pc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                xnkk  (maxq,maxq,0:jmax,0:jmax,0:jmax,2)  
      common /nme_oct/ q3p    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                 q3cp   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                 q3pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                 q3cpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                 q1pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                 q1cpred(maxq,maxq,0:jmax,0:jmax,0:jmax) 
      common /neutron/qn    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                qnred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcnred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0nc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                q3n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     7                q3cn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     8                q3nred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     9                q3cnred(maxq,maxq,0:jmax,0:jmax,0:jmax)       

      common /qmoment/ qtot(maxq)
      common /mathco/ zero,one,two,half,third,pi  
      data r0 /1.2d0/ 


  201 format ('  J  K1 K2 iq     beta2   beta3    qtot       E       ',
     1        '   n^J(q,q)   <N>     <Z>')
  202 format (4i3,3f10.3,2f12.4,2f9.4)
c   
      write(*,201) 
c     ................................. initialization
      r00 = r0*amas**third
      pi  = 4.0d0 * atan2(1.0d0,1.0d0)
      fac = dsqrt(16*pi/5)
c     ................................ print diagonal (in q space) matrix elements
      do iq1 = 1,maxmp
	   ee  = dreal(hjqqkk(jproj,iq1+k1m/2*maxmp,iq1+k2m/2*maxmp)) 
	   gg  = dreal(njqqkk(jproj,iq1+k1m/2*maxmp,iq1+k2m/2*maxmp)) 
	   xnn = dreal(xnkk(iq1,iq1,jproj,k1,k2,1)) 
	   xpp = dreal(xnkk(iq1,iq1,jproj,k1,k2,2)) 
	   qtot(iq1) = (3*amas*r00**2)/(4*pi)*gcmbet(iq1)*fac
           write(*,202),jproj,k1,k2,iq1,gcmbet(iq1),
     &               gcmgam(iq1),qtot(iq1),ee,gg,xnn,xpp
      enddo ! iq1
c     ................................ print matrix of existing matrix elements
      do iq1 = 1,maxmp
      do iq2 = 1,maxmp
	   iqk1=iq1+k1m/2*maxmp
	   iqk2=iq2+k2m/2*maxmp
	   gg  = dreal(njqqkk(jproj,iqk1,iqk2))  
          if (abs(gg).gt.0.0d0) then
            iexst(iqk1,iqk2) = log(abs(gg)) 
          endif
      enddo ! iq1
      enddo ! iq2

      END  
c______________________________________________________________________________
      subroutine readkernel(icase,iq1,iq2,jprmi,jprma,jdf,elem)

c..............................................................................
c     write matrix elements for Hill-Wheeler-Griffin solver to file lou      .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'dic.par'
! 
      complex*16 nkk,hkk,xnkk 
      complex*16 qp,qcp,qpred,qcpred,q0p,q0pc
      complex*16 q3p,q3cp,q3pred,q3cpred,q1pred,q1cpred
      complex*16 q3n,q3cn,q3nred,q3cnred
      complex*16 qn,qcn,qnred,qcnred,q0n,q0nc
      character*58 elem

      common /big   / qp    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcp   (maxq,maxq,0:jmax,0:jmax,0:jmax), 
     2                qpred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0p   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0pc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                xnkk  (maxq,maxq,0:jmax,0:jmax,0:jmax,2)  
      common /nme_oct/ q3p    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                 q3cp   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                 q3pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                 q3cpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                 q1pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                 q1cpred(maxq,maxq,0:jmax,0:jmax,0:jmax) 
      common /neutron/qn    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                qnred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcnred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0nc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                q3n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     7                q3cn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     8                q3nred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     9                q3cnred(maxq,maxq,0:jmax,0:jmax,0:jmax)       
          
           
      common /kerne/  nkk(0:jmax,0:jmax,0:jmax),
     &                hkk(0:jmax,0:jmax,0:jmax) 
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax)
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gfviv / iv(-igfv:igfv)
      common /mathco/ zero,one,two,half,third,pi 
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp,lrpa  
 
  401 format (10i5)
  402 format (8e15.8)
      character(len=*), parameter ::  format1 = "(3i5,4x)", 
     &                                format2 = "(4e15.8)"

!      iho=1
!      if(iho.eq.1) then
!      call kernel_hoa(iq1,iq2,jprma)      
!      return
!      endif
      do jproj = jprmi,jprma,jdf 
      do k1  =  0, jproj, 2 
      do k2  =  0, jproj, 2 
         nkk(jproj,k1,k2)           = 0.d0
         hkk(jproj,k1,k2)           = 0.d0
         qp(iq2,iq1,jproj,k1,k2)    = 0.d0
         qp(iq1,iq2,jproj,k1,k2)    = 0.d0
         qpred(iq2,iq1,jproj,k1,k2) = 0.d0
         qpred(iq1,iq2,jproj,k1,k2) = 0.d0
         q3p(iq1,iq2,jproj,k1,k2)    = 0.d0
         q3p(iq2,iq1,jproj,k1,k2)    = 0.d0
         q3pred(iq2,iq1,jproj,k1,k2) = 0.d0
         q3pred(iq1,iq2,jproj,k1,k2) = 0.d0
         q1pred(iq2,iq1,jproj,k1,k2) = 0.d0
         q1pred(iq1,iq2,jproj,k1,k2) = 0.d0

         q0p(iq2,iq1,jproj,k1,k2)   = 0.d0
         q0p(iq1,iq2,jproj,k1,k2)   = 0.d0
         xnkk(iq1,iq2,jproj,k1,k2,1)= 0.d0
         xnkk(iq2,iq1,jproj,k1,k2,1)= 0.d0
         xnkk(iq2,iq1,jproj,k1,k2,2)= 0.d0
         xnkk(iq1,iq2,jproj,k1,k2,2)= 0.d0
       enddo ! k2 
       if(icase.ne.0) then
         qpred(iq1,iq2,jproj,k1,k1max+2)  = 0.d0
         qpred(iq2,iq1,jproj,k1,k1max+2)  = 0.d0
         qpred(iq1,iq2,jproj,k1,-k1max-2) = 0.d0
         qpred(iq2,iq1,jproj,k1,-k1max-2) = 0.d0
      endif
       enddo ! k1 
      enddo ! jproj

c    ......................................................... 
      icheck = 0
      do jpr = jprmi,jprma,jdf 
c       .................................... write matrix elements to tape 
      !   read (lou,401,end=403) jproj   
        read (lou,format1) jproj, k1r, k2r 
      !   read (lou,402) bet1,gam1,bet2,gam2
c        if(abs(gcmbet(iq1)-bet1).gt.0.01) stop ' beta1 is not correct !'
c        if(abs(gcmgam(iq1)-gam1).gt.1.0)  stop ' gam1  is not correct !'
c        if(abs(gcmbet(iq2)-bet2).gt.0.01) stop ' beta2 is not correct !'
c        if(abs(gcmgam(iq2)-gam2).gt.1.0)  stop ' gam2  is not correct !'
c    .........................................................
c        if(iv(jproj).gt.zero) then         ! for even spin 
c           ik0     = 0 
c        endif
c        if(iv(jproj).lt.zero) then       ! for odd spin 
c           ik0     = 2
c        endif 
c    .........................................................
         ik0     = 0 
         k1max   = 0 !  icase*jproj
         do k1  =  ik0, k1max, 1 
         do k2  =  ik0, k1max, 1 
            !    read (lou,401)  k1r,k2r  
               read (lou,format2)  nkk(jproj,k1,k2), hkk(jproj,k1,k2)   
               read (lou,format2)  xnkk(iq1,iq2,jproj,k1,k2,1),
     &                         xnkk(iq1,iq2,jproj,k1,k2,2) 
!                read (lou,402)  qp(iq1,iq2,jproj,k1,k2),
!      &                         qcp(iq1,iq2,jproj,k1,k2) 
               read (lou,format2)  qpred(iq1,iq2,jproj,k1,k2),
     &                         qcpred(iq1,iq2,jproj,k1,k2) 
!        if(jproj.eq.0) write(*,*) 'qpred',qpred(iq1,iq2,jproj,k1,k2)   
               read (lou,format2)  q0p(iq1,iq2,jproj,k1,k2),
     &                         q0pc(iq1,iq2,jproj,k1,k2) 
!                read (lou,402)  q3p(iq1,iq2,jproj,k1,k2),
!      &                         q3cp(iq1,iq2,jproj,k1,k2)
!                read (lou,402)  q3pred(iq1,iq2,jproj,k1,k2),
!      &                         q3cpred(iq1,iq2,jproj,k1,k2)
!                read (lou,402)  q1pred(iq1,iq2,jproj,k1,k2),
!      &                         q1cpred(iq1,iq2,jproj,k1,k2)
c     ....................................... qcp, qcpred were set to zero for iq1=iq2 cases in PKC code
c          q0p    = <J,K1,q1 ||e*r^2|| J,K2,q2>
c          qp     = <J,K1,q1 ||e*Q_2|| J,K2,q2>
c          qcp    = <J,K1,q2 ||e*Q_2|| J,K2,q1>
c          qpred  = <J,K1,q1 ||e*Q_2|| J+2,K2,q2>
c          qcpred = <J,K1,q2 ||e*Q_2|| J+2,K2,q1>
c          q3p    = <J,K1,q1  ||e*Q_3|| J+1,K2,q2>
c          q3cp   = <J,K1,q2  ||e*Q_3|| J+1,K2,q1>
c          q3pred  = <J,K1,q1 ||e*Q_3|| J+3,K2,q2>
c          q3cpred = <J,K1,q2 ||e*Q_3|| J+3,K2,q1>
c          q1pred  = <J,K1,q1 ||e*D_1|| J+1,K2,q2>
c          q1cpred = <J,K1,q2 ||e*D_1|| J+1,K2,q1>
c     ........................................ 
          if(iq1.ne.iq2) then 
            ! qp(iq2,iq1,jproj,k1,k2)    = qcp(iq1,iq2,jproj,k1,k2) 
            qpred(iq2,iq1,jproj,k1,k2) = qcpred(iq1,iq2,jproj,k1,k2) 
            ! q3p(iq2,iq1,jproj,k1,k2)   = q3cp(iq1,iq2,jproj,k1,k2) 
            ! q3pred(iq2,iq1,jproj,k1,k2)= q3cpred(iq1,iq2,jproj,k1,k2)
            ! q1pred(iq2,iq1,jproj,k1,k2)= q1cpred(iq1,iq2,jproj,k1,k2) 
            q0p(iq2,iq1,jproj,k1,k2)   = q0pc(iq1,iq2,jproj,k1,k2) 
            xnkk(iq2,iq1,jproj,k1,k2,1)= xnkk(iq1,iq2,jproj,k1,k2,1)
            xnkk(iq2,iq1,jproj,k1,k2,2)= xnkk(iq1,iq2,jproj,k1,k2,2) 
          endif 
c             goto 98
c    ............................................... neutron part     
!           read (lou,402) qn(iq1,iq2,jproj,k1,k2),
!      &                   qcn(iq1,iq2,jproj,k1,k2) 
          read (lou,format2) qnred(iq1,iq2,jproj,k1,k2),
     &                   qcnred(iq1,iq2,jproj,k1,k2) 
          read (lou,format2) q0n(iq1,iq2,jproj,k1,k2),
     &                   q0nc(iq1,iq2,jproj,k1,k2)  
!           read (lou,402) q3nred(iq1,iq2,jproj,k1,k2),
!      &                   q3cnred(iq1,iq2,jproj,k1,k2)   
!           read (lou,402) q3n(iq1,iq2,jproj,k1,k2),
!      &                   q3cn(iq1,iq2,jproj,k1,k2)           


          if(iq1.ne.iq2) then 
            ! qn(iq2,iq1,jproj,k1,k2)    = qcn(iq1,iq2,jproj,k1,k2) 
            qnred(iq2,iq1,jproj,k1,k2) = qcnred(iq1,iq2,jproj,k1,k2) 
            ! q3n(iq2,iq1,jproj,k1,k2)   = q3cn(iq1,iq2,jproj,k1,k2) 
            ! q3nred(iq2,iq1,jproj,k1,k2)= q3cnred(iq1,iq2,jproj,k1,k2) 
            q0n(iq2,iq1,jproj,k1,k2)   = q0nc(iq1,iq2,jproj,k1,k2)  
          endif 
                           
 98       enddo ! k2 
c     ....................... additional E2 matrix elements <J,K1 | Q_2 | J+2,K2>
      !     if(icase.ne.0)
!      &    read (lou,402)  qpred(iq1,iq2,jproj,k1,k1max+2),
!      &                    qpred(iq2,iq1,jproj,k1,k1max+2)
          enddo ! k1   		                
      icheck = icheck + 1
  403 continue    
  99  enddo ! jproj
      write(22,*) elem
c     ...................................... if the elem is not calculated
      if(icheck.eq.0) then
      write(12,*) elem, 'does not exist ...'
      do jproj = jprmi,jprma,jdf 
         do k1  =  0, jproj, 2 
         do k2  =  0, jproj, 2 
            nkk(jproj,k1,k2)           = 0.d0
            hkk(jproj,k1,k2)           = 0.d0
            qp(iq2,iq1,jproj,k1,k2)    = 0.d0
            qp(iq1,iq2,jproj,k1,k2)    = 0.d0
            qpred(iq2,iq1,jproj,k1,k2) = 0.d0
            qpred(iq1,iq2,jproj,k1,k2) = 0.d0
            q3p(iq1,iq2,jproj,k1,k2)    = 0.d0
            q3p(iq2,iq1,jproj,k1,k2)    = 0.d0
            q3pred(iq2,iq1,jproj,k1,k2) = 0.d0
            q3pred(iq1,iq2,jproj,k1,k2) = 0.d0
            q1pred(iq2,iq1,jproj,k1,k2) = 0.d0
            q1pred(iq1,iq2,jproj,k1,k2) = 0.d0

            q0p(iq2,iq1,jproj,k1,k2)   = 0.d0
            q0p(iq1,iq2,jproj,k1,k2)   = 0.d0
            xnkk(iq1,iq2,jproj,k1,k2,1)= 0.d0
            xnkk(iq2,iq1,jproj,k1,k2,1)= 0.d0
            xnkk(iq2,iq1,jproj,k1,k2,2)= 0.d0
            xnkk(iq1,iq2,jproj,k1,k2,2)= 0.d0
          enddo ! k2 
          if(icase.ne.0) then
            qpred(iq1,iq2,jproj,k1,k1max+2)  = 0.d0
            qpred(iq2,iq1,jproj,k1,k1max+2)  = 0.d0
            qpred(iq1,iq2,jproj,k1,-k1max-2) = 0.d0
            qpred(iq2,iq1,jproj,k1,-k1max-2) = 0.d0
         endif
          enddo ! k1 
      enddo ! jproj
      endif
      return
      end subroutine readkernel
c______________________________________________________________________________
      subroutine reader()

c     .........................................................................
c     reader sets from input file     
c     .........................................................................
      implicit real*8 (a-h,o-z)
!
      include 'dic.par'
!   
      character nucnam*2                                        ! nucnuc 
!  

      common /mathco/ zero,one,two,half,third,pi 
      common /nucnuc/ amas,nama,nneu,npro,nucnam  
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp,lrpa 
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gcmcut/ zeta(5),jmx,icase,kmx 
      common /guagea/ mphi(2),dphi(2),phi(2),iphi
      common /eulera/ nbet,nphi 
      common /basnnn/ n0f,n0b
!............................
  80  format(10x,i5)
  81  format(10x,i6)
  90  format(10x,5e9.2)
!----------------------------------
!    reading from data
!............................         
      open(lin,file='data',status='old')  
!---- Nucleus under consideration
      read(lin,'(a2,i4)') nucnam,nama
      amas = nama       
      call nucleus(2,npro,nucnam)
      nneu = nama - npro
       write(*,*) nucnam,nama,npro
      read(lin,80) ioe
      read(lin,80) n0f 
      read(lin,80) icase 
      read(lin,80) nphi
      mphi(1) = nphi 
      mphi(2) = nphi 
      read(lin,80) nbet
      read(lin,80) jmx 
      if(nbet.eq.1) jmx = 0
      read(lin,80) kmx 
      read(lin,90) zeta    
!............................
!    reading from betagam
!............................ 
      open(10,file='betgam.dat',status='old')  
      read(10,81) maxmp
      do iq1 = 1, maxmp 
         read(10,*) gcmbet(iq1),gcmgam(iq1) !,gcmr2(iq1)  
      enddo ! iq1
      return
      end
c________________________________________________________________________
      subroutine re_name(betac1,gammac1,r2c1,betac2,gammac2,r2c2,name,elem,dens)
c.......................................................................
      implicit double precision (a-h,o-z)
      include 'dic.par'

      character*28 outputwf,outputwf1,outputwf2                 ! ttttt1
      character*58 olpf,elem,dens
      character*1 sign1,sign2
      character*2 nucnam                                        ! nucnuc
      character*3 name0
      character :: signb21,signb31,signb22,signb32
      integer, dimension(6) :: name1,name2
!---------------------------
      common /nucnuc/ amas,nama,nneu,npro,nucnam
      common /guagea/ mphi(2),dphi(2),phi(2),iphi
      common /eulera/ nbet,nphi 
      common /basnnn/ n0f,n0b
!---------------------------
              name0 = './'
              if(betac1.ge.0.d0) sign1='+'
              if(betac1.lt.0.d0) sign1='-'
              if(betac2.ge.0.d0) sign2='+'
              if(betac2.lt.0.d0) sign2='-'
              ab2c1  = abs(betac1)
              ab3c1  = abs(gammac1)
              name11 = ab2c1+48
              name21 = mod(ab2c1*10,10.d0)  +48
              name31 = mod(ab2c1*100,10.d0) +48
              name01 = ab3c1+48
              name41 = mod(ab3c1*10,10.d0)+48
              name51 = mod(ab3c1*100,10.d0)   +48
              name61 = r2c1+48
              name71 = mod(r2c1*10,10.d0) +48
              namer1 = mod(r2c1*100.0001,10.d0) +48
!---------------------------
              ab2c2  = abs(betac2)
              ab3c2  = abs(gammac2)
              name12 = ab2c2+48
              name22 = mod(ab2c2*10,10.d0)  +48
              name32 = mod(ab2c2*100,10.d0) +48
              name02 = ab3c2+48
              name42 = mod(ab3c2*10,10.d0)+48
              name52 = mod(ab3c2*100,10.d0)   +48
              name62 = r2c2+48
              name72 = mod(r2c2*10,10.d0) +48
              namer2 = mod(r2c2*100.0001,10.d0) +48
!---------------------------  
              jphi   = mphi(1) 
              name81  = mod(jphi/10,10) + 48 
              name82  = mod(jphi,10) + 48 
              name91  = mod(nbet/10,10) + 48 
              name92  = mod(nbet,10) + 48 

              name_nf1 = mod(n0f/10,10) + 48
              name_nf2 = mod(n0f,10) + 48


!---------------------------
              beta2_1 = betac1 
              beta3_1 = gammac1
              beta2_2 = betac2
              beta3_2 = gammac2
              if(beta2_1 >= 0.d0) signb21 = '+'
              if(beta2_1 < 0.d0)  signb21 = '-'
              if(beta3_1 >= 0.d0) signb31 = '+'
              if(beta3_1 < 0.d0)  signb31 = '-'
              if(beta2_2 >= 0.d0) signb22 = '+'
              if(beta2_2 < 0.d0)  signb22 = '-'
              if(beta3_2 >= 0.d0) signb32 = '+'
              if(beta3_2 < 0.d0)  signb32 = '-'
              !-------
              abs2c1 = abs(beta2_1)
              abs3c1 = abs(beta3_1)
              name1(1) = abs2c1 + 48 !In ASCII, character '0' start from 48. 
              name1(2) = mod(abs2c1*10,10.d0)+48
              name1(3) = mod(abs2c1*100,10.d0)+48
              name1(4) = abs3c1+48
              name1(5) = mod(abs3c1*10,10.d0)+48
              name1(6) = mod(abs3c1*100,10.d0)+48
              !------
              abs2c2 = abs(beta2_2)
              abs3c2 = abs(beta3_2)
              name2(1) = abs2c2 + 48 !In ASCII, character '0' start from 48. 
              name2(2) = mod(abs2c2*10,10.d0)+48
              name2(3) = mod(abs2c2*100,10.d0)+48
              name2(4) = abs3c2+48
              name2(5) = mod(abs3c2*10,10.d0)+48
              name2(6) = mod(abs3c2*100,10.d0)+48
              !-----
            !   name_nf1 = mod(BS%HO_sph%n0f/10,10) + 48
            !   name_nf2 = mod(BS%HO_sph%n0f,10) + 48
              !-----
            !   nphi_1 = mod(projection_mesh%nphi(1)/10,10) + 48
            !   nphi_2 = mod(projection_mesh%nphi(1),10) + 48
            !   nbeta_1 = mod(projection_mesh%nbeta/10,10) + 48
            !   nbeta_2 = mod(projection_mesh%nbeta,10) + 48
c              if(mphi(1).lt.10)
c     &         name6  = '0'//char(jphi)
!---------------------------  
       if(gammac1.lt.0.d0) then
        outputwf1 =name0//'/dio'//sign1//char(name11)
     &           //char(name21)//char(name31)
     &       //'.-'//char(name01)//char(name41)//char(name51)//'.wf'
       else
        outputwf1 =name0//'/dio'//sign1//char(name11)
     &           //char(name21)//char(name31)
     &      //'.+'//char(name01)//char(name41)//char(name51)//'.wf'
       endif
!
       if(gammac2.lt.0.d0)then
        outputwf2 =name0//'/dio'//sign2//char(name12)
     &           //char(name22)//char(name32)
     &      //'.-'//char(name02)//char(name42)//char(name52)//'.wf'
       else
        outputwf2 =name0//'/dio'//sign2//char(name12)
     &           //char(name22)//char(name32)
     &          //'.+'//char(name02)//char(name42)//char(name52)//'.wf'
       endif

!--------------------- 
      olpf =name0//'/ovlp.'//char(name)//'D.'
     &          //char(name11)//char(name21)//char(name31)
     &          //char(name41)//char(name51)//'.'//char(name61)
     &          //char(name71)//char(namer1)//
     &          '_'
     &          //char(name12)//char(name22)//char(name32)
     &          //char(name42)//char(name52)//'.'//char(name62)
     &          //char(name72)//char(namer2)//'.dat'   
!--------------------- 
!       elem ='kern.'//char(name)//'D'
!      &          //'_eMax'//char(name_nf1)//char(name_nf2)
!      &          //'.'//char(name81)//char(name82)
!      &          //'.'//char(name91)//char(name92)//'.'
!      &          //sign1//char(name11)//char(name21)
!      &          //char(name31)//char(name41)//char(name51)//'.'
!      &          //char(name61)//char(name71)//char(namer1)//
!      &          '_'
!      &          //sign2//char(name12)//char(name22)//char(name32)
!      &          //char(name42)//char(name52)//'.'
!      &          //char(name62)//char(name72)//char(namer2)//'.elem'
     
       elem = 'Proj_.'//int2str(nama)//nucnam
     &    //'_kern.'//char(name)//'D' 
     &     //'_eMax'//char(name_nf1)//char(name_nf2) 
     &     //'.'//char(name81)//char(name82) 
     &     //'.'//char(name91)//char(name92) 
     &     //signb21//char(name1(1))//char(name1(2))//char(name1(3)) 
     &     //signb31//char(name1(4))//char(name1(5))//char(name1(6)) 
     &     //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) 
     &     //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.elem' 

!--------------------- 
      dens ='/Proj_kern.'//char(name)//'D'
     &          //'.'//char(name81)//char(name82)
     &          //'.'//char(name91)//char(name92)//'.'
     &          //sign1//char(name11)//char(name21)
     &          //char(name31)//char(name41)//char(name51)//'.'
     &          //char(name61)//char(name71)//char(namer1)//
     &          '_'
     &          //sign2//char(name12)//char(name22)//char(name32)
     &          //char(name42)//char(name52)//'.'
     &          //char(name62)//char(name72)//char(namer2)//'.dens'
      return
      end
c______________________________________________________________________________
      function int2str(i) result(str)
            implicit none
            integer, intent(in) :: i
            character(len=:), allocatable :: str
            character(len=32) :: buf
            write(buf,'(I0)') i
            str = trim(buf)
      end function int2str

      subroutine kernels(jmn,jmx,jdf,iq1,iq2,icase,nmaxdi,lpr) 

c..............................................................................  
c     copy kernels (nkk,hkk) to (njqqkk,hjqqkk) for HW equation
c.............................................................................. 
      implicit real*8 (a-h,o-z)
      include 'dic.par'
      logical lpr
c 
      COMPLEX*16 on,etot,calp,cgam,cpi
      COMPLEX*16 ctemp1,njqqkk,hjqqkk   
      complex*16 nikk,hikk,nkk,hkk,novlpi,hovlpi 
      complex*16 onang,etotang,qlmang,xnang   
      complex*16 q00p,bet0ang,be0 
      dimension nmaxdi(0:jmax),xn(2)
      character*8 name
c  
      common /kerne/  nkk(0:jmax,0:jmax,0:jmax),
     &                hkk(0:jmax,0:jmax,0:jmax) 
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax)
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gfviv / iv(-igfv:igfv)
      common /gfvsq / sq(0:igfv)  
      common /mathco/ zero,one,two,half,third,pi  
      common /averpn/ ala0(2)
      common /nucnuc/ amas,nama,npr(2),nucnam
      data   eps /1.d-5/ 
! 
   30 format('nkk(',2i2,2x,2i2,')=',2f15.8,2x,'hkk=',2f15.8) 
!--------------------------------------------- 
 111  if(icase.eq.0) then  
!--- loop over total angular momentum 
      k1  = 0 
      k2  = 0 
      k1m = -k1
      k2m = -k2
      do iis = jmn, jmx, jdf
         nmaxdi(iis)  = maxmp      
         hjqqkk(iis,iq1,iq2) = hkk(iis,0,0)   
         njqqkk(iis,iq1,iq2) = nkk(iis,0,0)  
!------ symmetry in Kernels 
         hjqqkk(iis,iq2,iq1) = hkk(iis,0,0)   
         njqqkk(iis,iq2,iq1) = nkk(iis,0,0) 
         if(lpr) then
            write(*,30) iq1,iq2,k1m,k2m,nkk(iis,0,0),hkk(iis,0,0)
         endif
  119 enddo ! iis	  
      return 
      endif
!--------------------------------------------- 
 222   if(icase.eq.1) then 
       do iis= jmn, jmx, jdf
          if(iis.eq.1) goto 120
c----- initialization  
          fac05 = one 
          fac3 = fac05*(2*iis+1)/(8*pi1**2)  
c----------------------------------   
	        if(iv(iis).gt.zero) then       ! for even spin
                 nmaxdi(iis)  = (iis/2+1)*maxmp
                     iki0     = 0 
              endif
	        if(iv(iis).lt.zero) then       ! for odd spin
                 nmaxdi(iis)  = (iis-1)/2*maxmp
                      iki0    = 2
              endif   
c     ..................................................................................
!      K1, K2 start from 0 or 2 and only the matrix elements with even K are non-zero
!      For even spin J, K1 and K2 start from iki0 =0
!      For odd  spin J, K1 and K2 start from iki0 =2
!      the mesh point i is given by: i = iq + (K-iki0)/2*maxmp  
c     ..................................................................................
!      For example:
!**************************************************************************************************
!      1)  J=4, iki0=0
!          K=   iki0,                     2,                               4
!          i= 1,2,...,maxmp,  maxmp+1,maxmp+2,...,2*maxmp,  2*maxmp+1,2*maxmp+2,...,2*maxmp+maxmp
!
!      2)  J=7, iki0=2
!          K=   iki0,                     4,                               6
!          i= 1,2,...,maxmp,  maxmp+1,maxmp+2,...,2*maxmp,  2*maxmp+1,2*maxmp+2,...,2*maxmp+maxmp
!**************************************************************************************************
              do k1=iki0,iis,2   ! K-value
              do k2=iki0,iis,2 
	           k1m = k1-iki0   ! mesh point
	           k2m = k2-iki0
                 hjqqkk(iis,iq1+k1m/2*maxmp,iq2+k2m/2*maxmp) 
     &           = hkk(iis,k1,k2)
                 njqqkk(iis,iq1+k1m/2*maxmp,iq2+k2m/2*maxmp) 
     &           = nkk(iis,k1,k2)  
              enddo ! k2
              enddo ! k1 
c---- symmetry in kernels 
             do k1=iki0,iis,2
             do k2=iki0,iis,2
	          k1m = k1-iki0   ! mesh point
	          k2m = k2-iki0
                hjqqkk(iis,iq2+k2m/2*maxmp,iq1+k1m/2*maxmp)
     &           = hkk(iis,k1,k2)
                njqqkk(iis,iq2+k2m/2*maxmp,iq1+k1m/2*maxmp)
     &           = nkk(iis,k1,k2) 
              enddo !
              enddo ! 
!----------output
       if(lpr) then
           do k1=iki0,iis,2
           do k2=iki0,iis,2
                k1m = k1-iki0   ! mesh point
                k2m = k2-iki0
                write(*,30) iq1,iq2,k1m,k2m, 
     &                      nkk(iis,k1,k2),hkk(iis,k1,k2) 
           enddo !
           enddo !
       endif 
!------------------------		   
  120 enddo ! iis
!------------------------		   
      endif
      return
!---- End
      end  
c______________________________________________________________________________ 
      subroutine hweq(jtot,nmaxd,kbon,lpr)  

c.............................................................................. 
!           Solution of Hill-Wheel equation 
!       The standard method is adopted to solve
!       the HW equation. More details can be found
!       in <The nuclear many-body problem> by P.Ring & P.Schuck
c.............................................................................. 
      implicit real*8 (a-h,o-z)
      include 'dic.par' 
      real*8 NN,HH
      complex*16 njqqkk,hjqqkk,xnkk 
      complex*16 qp,qcp,qpred,qcpred,q0p,q0pc
      complex*16 q3p,q3cp,q3pred,q3cpred,q1pred,q1cpred
      complex*16 q3n,q3cn,q3nred,q3cnred
      complex*16 qn,qcn,qnred,qcnred,q0n,q0nc
      logical lpr
!---------------------------------------------------------------------
      dimension HH(nmaxd,nmaxd),NN(nmaxd,nmaxd),GG(nmaxd,nmaxd),
     &          DD(nmaxd,nmaxd),FF(nmaxd,nmaxd)     
      DIMENSION WW(nmaxd,nmaxd),RR(nmaxd,nmaxd)
      dimension E(nmaxd),EN(nmaxd),Z(nmaxd) 
!--------------------------------------------------------------
      common /mathco/ zero,one,two,half,third,pi  
      common /gfviv / iv(-igfv:igfv)
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax) 
      common /heelwh/ FJKq(maxq*jmax,maxq*jmax,0:jmax),
     &                EJKq(maxq*jmax,0:jmax)  
      common /probab/ FN(maxq*jmax,maxq*jmax,0:jmax)
      common /gcmexp/ bet_aver(maxq*jmax,0:jmax), 
     &                gam_aver(maxq*jmax,0:jmax),
     &                qspec(maxq*jmax,0:jmax),
     &                xrp(maxq*jmax,0:jmax),
     &                q1(maxq*jmax,0:jmax),
     &                q2(maxq*jmax,0:jmax)   
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gcmcut/ zeta(5),jmx,icase,kmx 
      common /qmoment/ qtot(maxq)

      common /big   / qp    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcp   (maxq,maxq,0:jmax,0:jmax,0:jmax), 
     2                qpred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0p   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0pc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                xnkk  (maxq,maxq,0:jmax,0:jmax,0:jmax,2)  
      common /gcmobs/ xnn(maxq*jmax,jmax) 
      data IS/3/ 
c     ...........................................
   10 format(3x,5f18.12)
  100 format ('cden(',2i4,')=',2f12.8)
  101 format(3i4,f12.5,5f9.5)     
  200 format (//,' eigenvalues of the norm ',/)
  201 format (1h ,1p7e11.3)
  202 format (/,' cutoff value of the norm eigenvalues ',1pe10.2)
  203 format (//,' eigenvalues of the hamiltonian ',/)
  401 format(2x,2i3,1x,3f8.3,10f10.5) 
  402 format ('    K  iq   beta2  beta3    rms R  ' 
     1   '      g1       g2       g3 ')
  404 format ('    K  iq   beta2  beta3    rms R  '
     1   '     fn1      fn2      fn3 ')

  403 format (//,'  eigenvectors ',/) 
c     ............................................. set cutoff in eigenvalue of norm kernel
      if(jtot.le.4) EPS=zeta(jtot+1)
      if(jtot.gt.4) EPS=zeta(5)
c     .............................................. initialization  
      N = nmaxd 
        write(66,*) '----- matrix HH(iq2,iq1) ---- for J=', jtot 
      do iq1=1,N
         do iq2=1,N
            NN(iq2,iq1) 
     &      = dreal(njqqkk(jtot,iq2,iq1)) 
            HH(iq2,iq1) 
     &      = dreal(hjqqkk(jtot,iq2,iq1)*njqqkk(jtot,iq2,iq1))  
          if(abs(NN(iq2,iq1)).lt.1.d-10) then
             NN(iq2,iq1) = 0.d0
             HH(iq2,iq1) = 0.d0
           endif
         enddo !q2

        write(66,'(50f12.5)') (HH(iq2,iq1),iq2=1,N) 
      enddo !q1  
!    ------------------------------------------------------------ print out NN matrix
        write(66,*) '----- matrix NN(iq2,iq1) ---- for J=', jtot 
      do iq1=1,N
        write(66,'(50f12.5)') (NN(iq2,iq1),iq2=1,N) 
      enddo !q1  
!-----------------------------------------------------------------
!
!     solution of a generalized eigenvalue problem: HH*FF = E*NN*FF
!
!-----------------------------------------------------------------
      write(*,*) 'begin solution of HWG eq. ....for ..',jtot
      call HNDIAG(NMAXD,N,M,HH,NN,E,EN,FF,DD,GG,RR,WW,
     &                       Z,EPS,IS,IFL) 
!-----------------------------------------------------------------
       if(M.LE.0) return
       kbon = M
       do k=1,kbon 
          do iqk=1,N     
             FJKq(iqk,k,jtot) = FF(iqk,k)  
c             if(jtot.eq.0.and.k.eq.1)
c     &       write(21,'(i3,2f12.6)') iqk, FF(iqk,k)  
          enddo ! iqk 
            EJKq(k,jtot) = E(k)   
       enddo !k
c    ............................................... norm eigenvalue
      print 200
      print 201, (EN(i),i=1,kbon) 
c     .............................................. cutoff for norm eigenvalues
      print 202,EPS
c     .............................................  hamiltonian eigenvalues
      print 203
      print 201, (E(i),i=1,kbon)  

c     .............................................  hamiltonian eigenvectors
      print 403
      print 402
!     ............................................  RR: g(K,q)
      do iqk1=1,N
         k1 =2*Int((iqk1-1)/maxmp)          ! k value
         if(iv(jtot).lt.zero) k1 = k1 + 2   ! if jtot = 3, 5, 7, ...; k value
         iq1 = mod(iqk1,maxmp)              ! mesh point only in q-space
         if(iq1.eq.0) iq1 = maxmp  
         write(*,401) k1,iq1,gcmbet(iq1),gcmgam(iq1),gcmr2(iq1),
     &                (RR(iqk1,ik),ik=1,min(kbon,kmx))
      enddo ! iqk1     
c     ......................................... average value of beta and gamma
       do ik=1,kbon 
	    betav = zero 
	    gamav = zero  
          do iqk1=1,N   
             iq1 = mod(iqk1,maxmp)             ! mesh point only in q-space
             if(iq1.eq.0) iq1 = maxmp 
	       ivv   =  1
	       if(abs(gcmgam(iq1)-180.d0).lt.1) ivv = -1
	       betav = betav + RR(iqk1,ik)**2 * gcmbet(iq1)*ivv
	       gamav = gamav + RR(iqk1,ik)**2 * gcmgam(iq1)
          enddo ! iqk1  
	       bet_aver(ik,jtot) = betav 
         gam_aver(ik,jtot) = gamav 
	  enddo ! k 
c     ............................................ if lpr, test
       if(lpr) then
       do k=1,kbon
	    cden  = zero 
	    xnn(jtot,k) = zero
          do iqk1=1,N   
             k1 =2*Int((iqk1-1)/maxmp)          ! k value
             if(iv(jtot).lt.zero) k1 = k1 + 2  ! if jtot = 3, 5, 7, ...; k value
             iq1 = mod(iqk1,maxmp)             ! mesh point only in q-space
             if(iq1.eq.0) iq1 = maxmp  
             fns = zero
          do iqk2=1,N 
             k2 =2*Int((iqk2-1)/maxmp)          ! k value
             if(iv(jtot).lt.zero) k2 = k2 + 2  ! if jtot = 3, 5, 7, ...; k value
             iq2 = mod(iqk2,maxmp)             ! mesh point only in q-space
             if(iq2.eq.0) iq2 = maxmp 

               fns  = fns +  dreal(njqqkk(jtot,iqk2,iqk1))
     &                     *  FJKq(iqk2,k,jtot)  
               coef = FJKq(iqk1,k,jtot)* FJKq(iqk2,k,jtot)  
	       cden = cden + coef * dreal(njqqkk(jtot,iqk2,iqk1)) 
	       xnn(jtot,k) = xnn(jtot,k) 
     &                     + coef * dreal(xnkk(iq1,iq2,ji,k1,k2,1))
     &                            * dreal(njqqkk(jtot,iqk2,iqk1))  
          enddo ! iqk2 
               FN(iqk1,k,jtot) = fns            ! <q|k>       
          enddo ! iqk1 
c	    write(*,100) jtot,k,cden,xnn(jtot,k)
	  enddo ! k 
      print 404
        do iqk1=1,N   
         k1 =2*Int((iqk1-1)/maxmp)          ! k value
         if(iv(jtot).lt.zero) k1 = k1 + 2   ! if jtot = 3, 5, 7, ...; k value
         iq1 = mod(iqk1,maxmp)              ! mesh point only in q-space
         if(iq1.eq.0) iq1 = maxmp
         write(*,401) k1,iq1,gcmbet(iq1),gcmgam(iq1),gcmr2(iq1),
     &                (FN(iqk1,ik,jtot),ik=1,min(kbon,kmx))
         enddo ! iqk1 

	endif
!---- END  
      return 
      end 
c______________________________________________________________________________
      subroutine obsers(jmn,jmx,jdf)

c..............................................................................
c     calculate spectroscopic quadupole moment                                . 
c..............................................................................
!
      implicit real*8(a-h,o-z)
      include 'dic.par'
!
      complex*16 njqqkk,hjqqkk,xnkk 
      complex*16 qp,qcp,qpred,qcpred,q0p,q0pc
      complex*16 q3p,q3cp,q3pred,q3cpred,q1pred,q1cpred
      complex*16 q3n,q3cn,q3nred,q3cnred
      complex*16 qn,qcn,qnred,qcnred,q0n,q0nc
      complex*16 aa(maxq*jmax,maxq*jmax),bb(maxq*jmax,maxq*jmax),
     &           s(maxq*jmax,maxq*jmax,0:jmax) 
      complex*16 y(maxq*jmax),zz(maxq*jmax),w,x 

      common /mathco/ zero,one,two,half,third,pi
      common /gfviv / iv(-igfv:igfv)
      common /gfvsq / sq(0:igfv)  
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax)
      common /heelwh/ FJKq(maxq*jmax,maxq*jmax,0:jmax),
     &                EJKq(maxq*jmax,0:jmax) 
      common /gcmexp/ bet_aver(maxq*jmax,0:jmax), 
     &                gam_aver(maxq*jmax,0:jmax),
     &                qspec(maxq*jmax,0:jmax),
     &                xrp(maxq*jmax,0:jmax),
     &                q1(maxq*jmax,0:jmax),
     &                q2(maxq*jmax,0:jmax)   
      common /gcmbe3/ q3spec(maxq*jmax,0:jmax)   
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gcmsta/ nmaxdi(0:jmax),jkmax(0:jmax) 
      common /big   / qp    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcp   (maxq,maxq,0:jmax,0:jmax,0:jmax), 
     2                qpred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0p   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0pc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                xnkk  (maxq,maxq,0:jmax,0:jmax,0:jmax,2)  
      common /nme_oct/ q3p    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                 q3cp   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                 q3pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                 q3cpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                 q1pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                 q1cpred(maxq,maxq,0:jmax,0:jmax,0:jmax) 
      common /neutron/qn    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                qnred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcnred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0nc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                q3n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     7                q3cn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     8                q3nred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     9                q3cnred(maxq,maxq,0:jmax,0:jmax,0:jmax)
     
      common /groust/ egs
      common /eltran/ be2fi(maxq*jmax,maxq*jmax,0:jmax),
     &                be0fi(maxq*jmax,maxq*jmax,0:jmax)
      common /nucnuc/ amas,nama,nneu,npro,nucnam
!     
      write(*,*) ' Obser subroutine finished '
  203 format (i3,3f12.6)
      do ji = jmn,jmx,jdf 
         ji2   = 2*ji + 1
         kboni = jkmax(ji) 
c         if(ji.eq.1) goto 89
c     ............................................ print out (n^J(q,q))^(1/2)
         if(ji.eq.0) then
            do iq1=1,maxmp
	             gg  = dreal(njqqkk(ji,iq1,iq1)) 
               write(23,203) iq1,sqrt(gg),FJKq(iq1,1,ji),FJKq(iq1,2,ji)
            enddo ! iq1
         endif
c     ............................................ initialization
        im   = nmaxdi(ji)
        do iqk1 = 1, im
        do iqk2 = 1, im
           aa(iqk1,iqk2)  = hjqqkk(ji,iqk1,iqk2)*njqqkk(ji,iqk1,iqk2)
           bb(iqk1,iqk2)  = njqqkk(ji,iqk1,iqk2)
        enddo !  iqk1
        enddo !  iqk2
c
         jf   = ji 
         do li= 1,kboni
            lf= li  
c     ............................................ initialization
            do iqk1 = 1, im
               s(iqk1,li,ji)  = FJKq(iqk1,li,ji)  
            enddo !  iqk1
            be2 = zero 
            be3 = zero 
            be0 = zero 
            do iqki = 1,nmaxdi(ji)                !summation over q_i and K_i 
               iqi  = iqki  ! mod(iqki,maxmp)             ! mesh point only in q-space
               iki  = 0 ! 2*Int((iqki-1)/maxmp)       ! k value
!                                                   if(iqi.eq.0)       iqi = maxmp 
            do iqkf = 1,nmaxdi(jf)                !summation over q_i and K_i 
               iqf  = iqkf !mod(iqkf,maxmp)             ! mesh point only in q-space
               ikf  = 0    ! 2*Int((iqkf-1)/maxmp)       ! k value
!                                                        if(iqf.eq.0)       iqf = maxmp 
!    ..............................................
               coef  =  FJKq(iqki,li,ji)*FJKq(iqkf,lf,jf) 
               be2 = be2 + coef  * qp  (iqf,iqi,ji,iki,ikf)   
               be3 = be3 + coef  * q3p  (iqf,iqi,ji,iki,ikf)   
               be0 = be0 + coef  * q0p (iqi,iqf,ji,iki,ikf)     
c 
            enddo ! iqf 
          enddo ! iqi    
          qspec(li,ji) 
     &          = be2*dsqrt(16*pi/5)*wignei(ji,2,ji,ji,0,-ji)   
          q3spec(li,ji) 
     &          = be3*dsqrt(16*pi/7)*wignei(ji,3,ji,ji,0,-ji)   
          xrp(li,ji) 
     &          = dsqrt(be0)  ! /dsqrt(sq(ji2)) Aug. 17, 2013 (jmy) 
c       ........................................... accuracy of the calculation  
        ek = EJKq(li,ji) 
        do i=1,im
           y (i) = 0.0d0
           zz(i) = 0.0d0
           do j=1,im
              x = 0.0
            do l=1,im
               x = x + (aa(j,l)-ek*bb(j,l)) * s(l,li,ji)
            enddo  ! l
               w    = aa(i,j)-ek*bb(i,j)
               y(i) = y(i) + w*x  ! (H -EN)^2 f
               zz(i)= zz(i)+ w*s(j,li,ji) ! (H -EN)f
          enddo ! j
        enddo  ! i
        q1(li,ji) = 0.0d0
        q2(li,ji) = 0.0d0
        do i=1,im
          q1(li,ji) = q1(li,ji) + dreal(y (i)*y (i))
          q2(li,ji) = q2(li,ji) + dreal(zz(i)*zz(i))
        enddo ! i  
!        write(*,*) 'q1=',q1(li,ji),'  q2=',q2(li,ji)
c................................................. 
        enddo !li     
   89 enddo ! ji
!     
      return
      end
c______________________________________________________________________________
      subroutine spect(jmn,jdf)
c..............................................................................
c     calculate spectrum                                                      .
c                                                                             .
c     iqk = iq+k1m/2*maxmp                                                    .
c     k1m = k1-iki0, where iki0=0 (even J) or 2 (odd J)                       .
c     iq = number of points in q space                                        .
c     ik = number of points in k space                                        .
c     k  = number of eigenvalues of the norm kernel                           . 
c                                                                             .
c.............................................................................. 
c     EJKq(k,J)    : energy for state J_k                                     .
c     FJKq(iqk,k,J): wave function for state J_k                              .
c..............................................................................
!
      implicit real*8(a-h,o-z)
      include 'dic.par'
!
      complex*16 njqqkk,hjqqkk,xnkk 
      character*1 py(0:1)
      complex*16 qp,qcp,qpred,qcpred,q0p,q0pc
      complex*16 q3p,q3cp,q3pred,q3cpred,q1pred,q1cpred
      complex*16 q3n,q3cn,q3nred,q3cnred
      complex*16 qn,qcn,qnred,qcnred,q0n,q0nc

      common /mathco/ zero,one,two,half,third,pi
      common /gfviv / iv(-igfv:igfv)
      common /gfvsq / sq(0:igfv)  
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax)
      common /heelwh/ FJKq(maxq*jmax,maxq*jmax,0:jmax),
     &                EJKq(maxq*jmax,0:jmax) 
      common /gcmexp/ bet_aver(maxq*jmax,0:jmax), 
     &                gam_aver(maxq*jmax,0:jmax),
     &                qspec(maxq*jmax,0:jmax),
     &                xrp(maxq*jmax,0:jmax),
     &                q1(maxq*jmax,0:jmax),
     &                q2(maxq*jmax,0:jmax)  
      common /gcmbe3/ q3spec(maxq*jmax,0:jmax) 
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gcmsta/ nmaxdi(0:jmax),jkmax(0:jmax) 
      common /big   / qp    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcp   (maxq,maxq,0:jmax,0:jmax,0:jmax), 
     2                qpred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0p   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0pc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                xnkk  (maxq,maxq,0:jmax,0:jmax,0:jmax,2)  
      common /nme_oct/ q3p    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                 q3cp   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                 q3pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                 q3cpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                 q1pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                 q1cpred(maxq,maxq,0:jmax,0:jmax,0:jmax) 
      common /neutron/qn    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                qnred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcnred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0nc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                q3n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     7                q3cn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     8                q3nred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     9                q3cnred(maxq,maxq,0:jmax,0:jmax,0:jmax)
     
      common /groust/ egs
      common /gcmcut/ zeta(5),jmx,icase,kmx
      common /eltran/ be2fi(maxq*jmax,maxq*jmax,0:jmax),
     &                be0fi(maxq*jmax,maxq*jmax,0:jmax)
      common /nucnuc/ amas,nama,nneu,npro,nucnam
      data py(0:1)/'+','-'/
!     
  102 format(1x,92('-'),/)
  103 format('   Ground State energy  ',/)
  104 format('   E_g.s.= ',f10.4,' MeV')
  200 format('   excitation spectrum  ',/)
  201 format (' J^pi_i      E       E_ex      <beta2>   <beta3>    ' 
     &  ' Q^s_2     <N>       <Z>    rrms_p      accuracy')
  203 format (86(' '), '<H-<H>>^2')
  204 format (92(' '), '<H^2-<H>^2>^2 ')
  202 format (i2,a1,i2,f12.4,3f10.4,f12.4,3f9.4,2f9.4) 
c..................................................initialization
      r00 = 1.2d0*amas**third
      egs = EJKq(1,0)  
      print 103
      write(*,104) egs
      print 102
      print 200
      print 201
      print 203
      print 204
      print 102
c    ...........................................
      do ji   = jmn,jmx,jdf   
            kbon  = min(jkmax(ji),kmx)  
         do li    = 1,kbon 
            xei   = EJKq(li,ji) - egs 
c    ............................................. initialization 
            qtot  = zero 
            xn    = zero
            xp    = zero 
            cden  = 0.d0
            do iqki = 1,nmaxdi(ji)                !summation over q_i and K_i 
               iqi  = mod(iqki,maxmp)             ! mesh point only in q-space
               iki  = 2*Int((iqki-1)/maxmp)       ! k value
               if(iqi.eq.0)       iqi = maxmp
!               if(iv(ji).lt.zero) iki = iki + 2   ! if jtot = 3, 5, 7, ...; k value 
            do iqkf = 1,nmaxdi(ji)                !summation over q_i and K_i 
               iqf  = mod(iqkf,maxmp)             ! mesh point only in q-space
               ikf  = 2*Int((iqkf-1)/maxmp)       ! k value
               if(iqf.eq.0)       iqf = maxmp
!               if(iv(ji).lt.zero) ikf = ikf + 2   ! if jtot = 3, 5, 7, ...; k value
!    ..............................................
               coef = FJKq(iqki,li,ji)*FJKq(iqkf,li,ji)   
               cden= cden + coef* dreal(njqqkk(ji,iqkf,iqki))   
!    .............................................. calculations of <N> & <Z>  
              xn = xn + coef * dreal(xnkk(iqi,iqf,ji,iki,ikf,1))
     &                       * dreal(njqqkk(ji,iqkf,iqki))   
              xp = xp + coef * dreal(xnkk(iqi,iqf,ji,iki,ikf,2))
     &                       * dreal(njqqkk(ji,iqkf,iqki))    
            enddo ! iqkf 
          enddo ! iqki 
          xrp(li,ji) = xrp(li,ji)/dsqrt(xp+0.0000001d0)          ! normalized to proton number
          write(*,202) ji, py(mod(ji,2)), li, EJKq(li,ji), xei, 
     &              bet_aver(li,ji), 
     &              gam_aver(li,ji),qspec(li,ji),
     &              xn,xp,xrp(li,ji), 
     &              q1(li,ji),q2(li,ji)    
        enddo !li    
   89 enddo ! ji
!
      print 102
      return
      end
c______________________________________________________________________________
      subroutine trans(jmn,jdf)

c..............................................................................
c     calculate reduced transition rates                                      .
c                                                                             .
c     iqk = iq+k1m/2*maxmp
c     k1m = k1-iki0, where iki0=0 (even J) or 2 (odd J) 
c     iq = number of points in q space                                        .
c     ik = number of points in k space                                        .
c     k  = number of eigenvalues of the norm kernel                           . 
c                                                                             .
c.............................................................................. 
c     EJKq(k,J)    : energy for state J_k                                     .
c     FJKq(iqk,k,J): wave function for state J_k                              .
c..............................................................................
!
      implicit real*8(a-h,o-z)
      include 'dic.par'
!
      complex*16 njqqkk,hjqqkk,xnkk 
      complex*16 qp,qcp,qpred,qcpred,q0p,q0pc
      complex*16 q3p,q3cp,q3pred,q3cpred,q1pred,q1cpred
      complex*16 q3n,q3cn,q3nred,q3cnred
      complex*16 qn,qcn,qnred,qcnred,q0n,q0nc

      common /mathco/ zero,one,two,half,third,pi
      common /gfviv / iv(-igfv:igfv)
      common /gfvsq / sq(0:igfv)  
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax)
      common /heelwh/ FJKq(maxq*jmax,maxq*jmax,0:jmax),
     &                EJKq(maxq*jmax,0:jmax) 
      common /gcmexp/ bet_aver(maxq*jmax,0:jmax), 
     &                gam_aver(maxq*jmax,0:jmax),
     &                qspec(maxq*jmax,0:jmax),
     &                xrp(maxq*jmax,0:jmax),
     &                q1(maxq*jmax,0:jmax),
     &                q2(maxq*jmax,0:jmax)   
      common /gcmbe3/ q3spec(maxq*jmax,0:jmax)  
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gcmsta/ nmaxdi(0:jmax),jkmax(0:jmax) 
      common /big   / qp    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcp   (maxq,maxq,0:jmax,0:jmax,0:jmax), 
     2                qpred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0p   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0pc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                xnkk  (maxq,maxq,0:jmax,0:jmax,0:jmax,2)  
      common /nme_oct/ q3p    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                 q3cp   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                 q3pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                 q3cpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                 q1pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                 q1cpred(maxq,maxq,0:jmax,0:jmax,0:jmax) 
      common /neutron/qn    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                qnred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcnred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0nc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                q3n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     7                q3cn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     8                q3nred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     9                q3cnred(maxq,maxq,0:jmax,0:jmax,0:jmax)
      common /groust/ egs
      common /gcmcut/ zeta(5),jmx,icase,kmx 
      common /eltran/ be2fi(maxq*jmax,maxq*jmax,0:jmax),
     &                be0fi(maxq*jmax,maxq*jmax,0:jmax)
      common /nucnuc/ amas,nama,nneu,npro,nucnam
!   
c..............................................................................
   20 format(2i2,3f12.6) 
  203 format (4i5,1f10.2,1f12.3,1f8.3,2f10.3)
 

  202 format (' i1 Ji i2 Jf   Ei    Ef      DeltaE',
     1        '  E2:(J_i->J_f) (J_f -> J_i)   |Q(t)|   ',
     1        ' beta(t)   <Ji||E2||Jf>   m(E0)   rho^2(E0)',
     1        '     M_n/M_p    eta ',/,
     2        '               MeV   MeV     MeV',
     3        '       e^2 fm^4    e^2 fm^4    e fm^2 ',
     4         '                  e fm^2      e fm^2       ') 
    
   
  204 format (4i3,2f12.3,1f12.3,9f12.3)
  205 format (' i1 Ji i2 Jf   Ei    Ef    DeltaE   ',
     1        'E3:(J_i->J_f)(J_f -> J_i) <Ji||E3||Jf>  ',
     1        '|Q_3(t)|    beta_3(t)     M_n/M_p   eta',/,
     2        '               MeV   MeV     MeV   ',
     3        '     e^2 fm^6    e^2 fm^6   ',
     4        '  e fm^3      e fm^3 ')   
 
 
  207 format (/,' ........   E2 transitions .......',/)
  208 format (/,' ........   E3 transitions .......',/)
  209 format (/,' ........   E1 transitions .......',/)
  210 format (' i1 Ji i2 Jf   Ei    Ef    DeltaE',
     1        '  E1:(J_i->J_f)(J_f -> J_i) <Ji||E1||Jf> |Q_1(t)|',
     1        '   |beta_1(t)|  ',/,
     2        '               MeV   MeV     MeV',
     3        '     e^2 fm^2    e^2 fm^2   ',
     4        '     e fm       e fm ')  
  
  211 format (4i3,2f15.3,1f15.3,9f15.3)
  
c..................................................initialization
      r00 = 1.2d0*amas**third
      fac = 4.d0*pi/(3.d0*npro*r00**2)
      fac3= 4.d0*pi/(3.d0*npro*r00**3)
      fac1= 4.d0*pi/(3.d0*npro*r00**1)
      egs = EJKq(1,0)  
      pniv= (amas-npro)/npro   ! N/Z 
c    ............................................. E2
      Lif = 2  
      print 207
      write(*,*) ' W.u./e^2fm^4=',fm2wu(amas,2)
      do ji = jmn,jmx,jdf 
         ji2   = 2*ji + 1
         kboni = min(jkmax(ji),kmx)   
            jfmax = min(ji+Lif,jmx)
            if(ji.eq.jmn) print 202
         do jf = ji, jfmax, Lif 
            jf2   = 2*jf + 1
            kbonf = min(jkmax(jf),kmx)  
            fac2  = wignei(ji,2,jf,0,0,0)*sq(jf2)*iv(ji-2)   
            fac22 = fac2**2   
         do li=1,kboni 
            xei  = EJKq(li,ji) - egs
         do lf=1,kbonf
            xef1 = EJKq(lf,jf) - egs
c    ............................................. initialization
            be2  = zero
            be0  = zero 
            b2n  = zero                            ! <Q_2> for neutrons
            ratio_np =zero
            do iqki = 1,nmaxdi(ji)                !summation over q_i and K_i 
               iqi  = iqki                         ! mod(iqki,maxmp)             ! mesh point only in q-space
               iki  = 0                            ! 2*Int((iqki-1)/maxmp)       ! k value
!                                                  if(iqi.eq.0)  iqi = maxmp 
            do iqkf = 1,nmaxdi(jf)                !summation over q_i and K_i 
               iqf  = iqkf                         ! mod(iqkf,maxmp)             ! mesh point only in q-space
               ikf  = 0                            ! 2*Int((iqkf-1)/maxmp)       ! k value
!                                                  if(iqf.eq.0)   iqf = maxmp 
!    .............................................. 
               coef  =  FJKq(iqki,li,ji)*FJKq(iqkf,lf,jf) 
!    .............................................. calculations of BE0 and spectroscopic Q
               if(jf.eq.ji) then  
               be2 = be2 + coef * qp  (iqf,iqi,ji,iki,ikf) 
               be0 = be0 + coef * q0p (iqf,iqi,ji,iki,ikf)
               b2n = b2n + coef * qn  (iqf,iqi,ji,iki,ikf) 
               endif
!    .............................................. calculations of BE2 
               if(jf.eq.ji+2) then
                  be2 = be2 + coef*qpred(iqf,iqi,ji,iki,ikf)
                  b2n = b2n + coef*qnred(iqf,iqi,ji,iki,ikf) 
               endif  
            enddo ! iqf 
           enddo ! iqi    
            if(abs(be2).gt.1.d-10) ratio_np=b2n/be2 
            eta   = ratio_np/pniv 
            be2s  = be2**2/ji2
            be22  = be2s*ji2/jf2
            beta2 = zero
            if(fac2.ne.0.d0)   
     &      beta2 = fac * dsqrt(be2s)/fac2   
            qt    = abs(beta2)*dsqrt(16*pi/5)/fac 
            if(jf.eq.ji.and.li.eq.lf)
     &      qspec(li,ji) 
     &            = be2s*dsqrt(16*pi/5)*wignei(ji,2,ji,ji,0,-ji)
            xme0  = be0     ! /sq(ji2)  is removed 2013/08/17
            rho0  = be0/r00**2
             rho02 = rho0**2 ! /ji2  is removed 2013/08/17
c     .............................................. print out 
          deif = xef1 - xei      
          if(be2.lt.0.d0) isgn=-1
          if(be2.ge.0.d0) isgn=+1
          print 204,li,ji,lf,jf,xei,xef1,deif,
     &        be2s,be22,qt,isgn*abs(beta2),be2,xme0,rho02,ratio_np,eta   
c     ..............................................
          be2fi(lf,li,ji) = be2 
          be0fi(lf,li,ji) = rho02    ! rho^2(E0)
        enddo !lf  	   
        enddo !li   
        enddo ! jf   
   89 enddo ! ji

c    ............................................. E1: J -> J+1
      Lif = 1
       print 209
      write(*,*) ' W.u./e^2fm^2=',fm2wu(amas,1)
      do ji = jmn,jmx,jdf
         ji2   = 2*ji + 1
         kboni = min(jkmax(ji),kmx)
            jfmax = min(ji+Lif,jmx)
            jf = ji + Lif
            jf2   = 2*jf + 1
            kbonf = min(jkmax(jf),kmx)
            fac11 = wignei(ji,1,jf,0,0,0)*sq(jf2)*iv(ji-1)
            fac12 = fac1**2 
            if(ji.eq.jmn) print 210
         do li=1,kboni
            xei  = EJKq(li,ji) - egs
         do lf=1,kbonf
            xef1 = EJKq(lf,jf) - egs
c    ............................................. initialization
            be1 = zero
            do iqki = 1,nmaxdi(ji)                !summation over q_i and K_i
               iqi  = mod(iqki,maxmp)             ! mesh point only in q-space
               iki  = 2*Int((iqki-1)/maxmp)       ! k value
               if(iqi.eq.0)       iqi = maxmp
            do iqkf = 1,nmaxdi(jf)                !summation over q_i and K_i
               iqf  = mod(iqkf,maxmp)             ! mesh point only in q-space
               ikf  = 2*Int((iqkf-1)/maxmp)       ! k value
               if(iqf.eq.0)       iqf = maxmp
!    ..............................................
               coef  =  FJKq(iqki,li,ji)*FJKq(iqkf,lf,jf)
               be1 = be1 + coef*q1pred(iqf,iqi,ji,iki,ikf)
c
            enddo ! iqf
          enddo ! iqi
            be1s   = be1**2/ji2
            be12  = be1s*ji2/jf2
            beta1 = zero
            if(fac11.ne.0.d0)
     &      beta1 = fac1* dsqrt(be1s)/fac11
            q1t    = abs(beta1)*dsqrt(16*pi/3)/fac1
c     .............................................. print out
          deif = xef1 - xei
          print 211,li,ji,lf,jf,xei,xef1,deif,
     &              be1s,be12,be1,q1t,abs(beta1)
c     ..............................................
        enddo !lf
        enddo !li
   69 enddo ! ji
      
c    ............................................. E3: J -> J+1
      Lif = 1
       print 208
      write(*,*) ' W.u./e^2fm^6=',fm2wu(amas,3)
      do ji = jmn,jmx,jdf
         ji2   = 2*ji + 1
         kboni = min(jkmax(ji),kmx)
            jfmax = min(ji+Lif,jmx)
            jf = ji + Lif
            jf2   = 2*jf + 1
            kbonf = min(jkmax(jf),kmx)
            fac31 = wignei(ji,3,jf,0,0,0)*sq(jf2)*iv(ji-3)
            fac32 = fac3**2 
            if(ji.eq.jmn) print 205
         do li=1,kboni
            xei  = EJKq(li,ji) - egs
         do lf=1,kbonf
            xef1 = EJKq(lf,jf) - egs
c    ............................................. initialization
            be3 = zero
            b3n = zero
            ratio_np =zero
            do iqki = 1,nmaxdi(ji)                !summation over q_i and K_i
               iqi  = mod(iqki,maxmp)             ! mesh point only in q-space
               iki  = 2*Int((iqki-1)/maxmp)       ! k value
               if(iqi.eq.0)       iqi = maxmp
            do iqkf = 1,nmaxdi(jf)                !summation over q_i and K_i
               iqf  = mod(iqkf,maxmp)             ! mesh point only in q-space
               ikf  = 2*Int((iqkf-1)/maxmp)       ! k value
               if(iqf.eq.0)       iqf = maxmp
!    ..............................................
               coef  =  FJKq(iqki,li,ji)*FJKq(iqkf,lf,jf)
               be3 = be3 + coef*q3p(iqf,iqi,ji,iki,ikf)
               b3n = b3n + coef*q3n(iqf,iqi,ji,iki,ikf)
c
            enddo ! iqf
          enddo ! iqi
            if(abs(be3).gt.1.d-10) ratio_np = b3n/be3
            eta   = ratio_np/pniv
            be3s   = be3**2/ji2
            be32  = be3s*ji2/jf2
            beta3 = zero
            if(fac31.ne.0.d0)
     &      beta3 = fac3* dsqrt(be3s)/fac31
            q3t    = abs(beta3)*dsqrt(16*pi/7)/fac3
c     .............................................. print out
          deif = xef1 - xei
          if(be3.lt.0.d0) isgn=-1
          if(be3.ge.0.d0) isgn=+1
          print 211,li,ji,lf,jf,xei,xef1,deif,
     &              be3s,be32,be3,q3t,isgn*abs(beta3),ratio_np,eta
c     ..............................................
        enddo !lf
        enddo !li
   79 enddo ! ji
   
   
c    ............................................. E3: J -> J+3
      Lif = 3 
       print 208
      write(*,*) ' W.u./e^2fm^6=',fm2wu(amas,3)
      do ji = jmn,jmx,jdf 
         ji2   = 2*ji + 1
         kboni = min(jkmax(ji),kmx)   
            jfmax = min(ji+Lif,jmx)
            jf = ji + Lif
            jf2   = 2*jf + 1
            kbonf = min(jkmax(jf),kmx)  
            fac31 = wignei(ji,3,jf,0,0,0)*sq(jf2)*iv(ji-3)   
            fac32 = fac3**2    
            if(ji.eq.jmn) print 205
         do li=1,kboni 
            xei  = EJKq(li,ji) - egs
         do lf=1,kbonf
            xef1 = EJKq(lf,jf) - egs
c    ............................................. initialization
            be3 = zero 
            b3n = zero 
            ratio_np =zero
            do iqki = 1,nmaxdi(ji)                !summation over q_i and K_i 
               iqi  = mod(iqki,maxmp)             ! mesh point only in q-space
               iki  = 2*Int((iqki-1)/maxmp)       ! k value
               if(iqi.eq.0)       iqi = maxmp 
            do iqkf = 1,nmaxdi(jf)                !summation over q_i and K_i 
               iqf  = mod(iqkf,maxmp)             ! mesh point only in q-space
               ikf  = 2*Int((iqkf-1)/maxmp)       ! k value
               if(iqf.eq.0)       iqf = maxmp 
!    .............................................. 
               coef  =  FJKq(iqki,li,ji)*FJKq(iqkf,lf,jf) 
               be3 = be3 + coef*q3pred (iqf,iqi,ji,iki,ikf)
               b3n = b3n + coef*q3nred (iqf,iqi,ji,iki,ikf)
c 
            enddo ! iqf 
          enddo ! iqi    
            if(abs(be3).gt.1.d-10) ratio_np = b3n/be3
            eta    = ratio_np/pniv
            be3s   = be3**2/ji2
            be32  = be3s*ji2/jf2
            beta3 = zero
            if(fac31.ne.0.d0)   
     &      beta3 = fac3* dsqrt(be3s)/fac31   
            q3t    = abs(beta3)*dsqrt(16*pi/7)/fac3   
c     .............................................. print out 
          deif = xef1 - xei
          print 211,li,ji,lf,jf,xei,xef1,deif,
     &              be3s,be32,be3,q3t,abs(beta3),ratio_np,eta  
c     ..............................................
        enddo !lf  	   
        enddo !li    
   99 enddo ! ji
!
      return
      end
c______________________________________________________________________________
      subroutine default()

c..............................................................................
      implicit real*8 (a-h,o-z)
!
      include 'dic.par' 
      common /baspar/ hom,hb0,b0  
      common /nucnuc/ amas,nama,nneu,npro,nucnam  
      common /initia/ vin,rin,ain,inin,inink 
      common /mathco/ zero,one,two,half,third,pi  
      common /physco/ hbc,alphi,r0  
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp,lrpa  
      common /pkcout/ nubd,iden,nubw
      common /numass/ amu
c     ........................................................... signs and factorials 
      call gfv()
      call binom()  
!     ............... density
      iden = 0 ! 1: calculate transition density; 0: no 
!     ............... constant
      amu    =  939.d0                 ! MeV
!---- tapes
      l6   = 10
      lin  = 3
      lou  = 26
      lwin = 1
      lwou = 2
      lplo = 11
      laka = 12
      lvpp = 13
      lrpa = 14 
      nubd = 21
      nubw = 22
    
      return
!-end-DEFAULT
      end
c______________________________________________________________________________
      blockdata block1

c..............................................................................
      implicit real*8 (a-h,o-z)
!
      character tp*1,tis*1,tit*8,tl*1                    ! common textex
!
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
!
!
!---- fixed texts
      data tp/'+','-'/,tis/'n','p'/,tit/'Neutron:','Proton: '/
      data tl/'s','p','d','f','g','h','i','j','k','l','m',
     &            'n','o','P','q','r','S','t','u','v','w',
     &            'x','y','z','0','0','0','0','0','0','0'/
!
!
!---- physical constants
      data hbc/197.328284d0/,r0/1.2d0/,alphi/137.03602/
!     data hbc/197.3289d0/,r0/1.2/,alphi/137.03602d0/
!
!---- mathemathical constants
!     are determined in GFV
      data zero/0.0d0/,one/1.d0/,two/2.d0/
      data half/0.5d0/,third/0.333333333333333333d0/
      data pi/3.141592653589793d0/ 
!-end-BLOCK1 
      end

c______________________________________________________________________________
      subroutine nucleus(is,npro,te)

c..............................................................................
!
!     is = 1 determines the symbol for a given proton number npro
!          2 determines the proton number for a given symbol te
!
c..............................................................................
!
      PARAMETER (MAXZ=140)
!
      CHARACTER TE*2,T*(2*MAXZ+2)
!
      T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
      T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
      T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
      T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
      T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
      T(201:240) = 'FmMdNoLrRfHaSgNsHsMr10111213141516171819'
      T(241:280) = '2021222324252627282930313233343536373839'
      T(281:282) = '40'
!
! ... Rf is called also as Ku (kurchatovium)
! ... Ha: IUPAC calls it as dubnium (Db). J.Chem.Educ. 1997, 74, 1258
! ... Ha is called also as Db (Dubnium)
!
      if (is.eq.1) then
         if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong' 
         te = t(2*npro+1:2*npro+2)
         return
      else
!
         do np = 0,maxz
            if (te.eq.t(2*np+1:2*np+2)) then
               npro = np
	       if (npro.gt.maxz) write(6,100) TE
               return
            endif
         enddo
!
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
      endif
!
      stop
!-END-NUCLEUS
      END
!
c______________________________________________________________________________
      subroutine gfv

c..............................................................................
!
!     Calculates sign, sqrt, factorials, etc. of integers and half int.
!
!     iv(n)  =  (-1)**n
!     sq(n)  =  sqrt(n)
!     sqi(n) =  1/sqrt(n)
!     sqh(n) =  sqrt(n+1/2)
!     shi(n) =  1/sqrt(n+1/2)
!     fak(n) =  n!
!     ibc(m,n) = m!/(n!(m-n)!)  
!     fad(n) =  (2*n+1)!!
!     fdi(n) =  1/(2*n+1)!!
!     fi(n)  =  1/n!
!     wf(n)  =  sqrt(n!)
!     wfi(n) =  1/sqrt(n!)
!     wfd(n) =  sqrt((2*n+1)!!)
!     gm2(n) =  gamma(n+1/2)
!     gmi(n) =  1/gamma(n+1/2)
!     wg(n)  =  sqrt(gamma(n+1/2))
!     wgi(n) =  1/sqrt(gamma(n+1/2))
!
c..............................................................................
      implicit double precision (a-h,o-z)
!
      include 'dic.par'
!      parameter (igfv = 100)
!
      common /gfviv / iv(-igfv:igfv)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvibc/ ibc(0:igfvbc,0:igfvbc)
      common /gfvfak/ fak(0:igfv)
      common /gfvfad/ fad(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvfdi/ fdi(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfvwfd/ wfd(0:igfv)
      common /gfvgm2/ gm2(0:igfv)
      common /gfvgmi/ gmi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
!
!---- mathemathical constants
!     data zero/0.0d0/,one/1.d0/,two/2.d0/
!     data half/0.5d0/,third/0.333333333333333333d0/
!     data pi/3.141592653589793d0/
!
      zero  = 0.d0
      one   = 1.d0
      two   = 2.d0
      half  = one/two
      third = one/3.d0
      pi    = 4*atan(one)
!
      iv(0)  = +1
      sq(0)  =  zero
      sqi(0) =  1.d30
      sqh(0) =  sqrt(half)
      shi(0) =  1/sqh(0)
      fak(0) =  one
      fad(0) =  one
      fi(0)  =  one
      fdi(0) =  one
      wf(0)  =  one
      wfi(0) =  one
      wfd(0)=  one
!     gm2(0) = Gamma(1/2) = sqrt(pi)
      gm2(0) =  sqrt(pi)
      gmi(0) =  1/gm2(0)
      wg(0)  =  sqrt(gm2(0))
      wgi(0) =  1/wg(0)
      do i = 1,igfv
         iv(i)         = -iv(i-1)
         iv(-igfv+i-1) = -iv(i)
         sq(i)  = dsqrt(dfloat(i))
         sqi(i) = one/sq(i)
         sqh(i) = sqrt(i+half)
         shi(i) = one/sqh(i)
         fak(i) = i*fak(i-1)
         fad(i) = (2*i+1)*fad(i-1)
         fi(i)  = one/fak(i)
         fdi(i) = one/fad(i)
         wf(i)  = sq(i)*wf(i-1)
         wfi(i) = one/wf(i)
         wfd(i) = sqrt(fad(i))
         gm2(i) = (i-half)*gm2(i-1)
         gmi(i) = one/gm2(i)
         wg(i)  = sqh(i-1)*wg(i-1)
         wgi(i) = one/wg(i)
      enddo
      ibc(0,0)= one
      do m=1,igfvbc
            do n=0,m
                  ibc(m,n)=fak(m)/(fak(n)*fak(m-n))
            enddo
      enddo   
      return
!-end-GFV
      end 
c______________________________________________________________________________
      subroutine binom()
c
c..............................................................................
C     THE ARRAY OF BINOMIAL COEFFICIENTS
C     BIN(I,J)= = I!/J!/(I-J)! 
c..............................................................................
      implicit double precision (a-h,o-z)
c
      parameter (IGFV = 100)
c
      common /bin0/ bin(0:IGFV,0:IGFV)
c
      do i = 0,IGFV
         do k = 0,IGFV
            bin(i,k) = 0.d0
         enddo   ! k
         bin(i,0)=1.d0
         bin(i,i)=1.d0
         do k = 1,i/2
            bin(i,k)   = dnint(bin(i,k-1)/dfloat(k)*dfloat(i-k+1))
            bin(i,i-k) = bin(i,k)
         enddo
      enddo
!      ngfv = IGFV
      return
c-end-BINOM
      end
c_______________________________________________________________________
      SUBROUTINE HNDIAG(NMAX,N,M,HH,NN,E,EN,FF,DD,GG,RR,WW,Z,EPS,IS,IFL)
!
c......................................................................C 
!     diagonalizes    HH*FF = E * NN*FF                                C
!     diagonalizes    NN*DD = EN * DD                                  C
!     eigenvalus of the norm EN<EPS are neglected                      C
!     IS = 0   only eigenvalues E,EN  and eigenvectors GG,DD           C
!     IS = 1   also eigenvectors  FF (possibly redundant)              C
!     IS = 2   also covariant components  WW = NN*FF                   C
!     IS = 3   also orthogonal wave functions RR = DD * GG             C
!              NN**(-1/2) * HH * NN**(-1/2) * RR  =  E  RR             C
!                                                                      C
!     HH(N,N),GG(N,N),NN(N,N),DD(N,N),FF(N,N),WW(N,N),RR(N,N)          C      
!     E(N)  eigenvalues of the generalized eigenvalue problem          C      
!     EN(N) eigenvalues of the norm matrix NN                          C
!     Z(N)  auxiliary field                                            C       
!                                                                      C
c......................................................................C
      implicit real*8 (a-h,o-z)
      real*8 NN
!                                                                             
      dimension HH(NMAX,N),NN(NMAX,N),GG(NMAX,NMAX),
     &          DD(NMAX,N),FF(NMAX,N) 
      DIMENSION WW(NMAX,N),RR(NMAX,N)
      dimension E(N),E1(N),EN(N),Z(N)
!------------------------------------                                                      
      IFL = 0   
      do K=1,N
         do I=1,N
            DD(I,K)=-NN(I,K)    
        enddo
      enddo
!-----diagonalization of NN
      CALL SDIAG(NMAX,N,DD,EN,DD,Z,1) 
!      
      do I=1,N
         EN(I)=-EN(I)+1.d-8
      enddo
!--- find out all the eigenvalues which are large enough compared with the largest one,
!    and drop the small ones
      do 1 I=1,N
         IF (EN(I)/EN(1).LT.EPS) GOTO 2
    1    CONTINUE
    2 M=I-1 
      IF (M.LE.0) THEN                                              
         WRITE(6,*) ' SUBR. HNDIAG: NO EIGENVALUE OF N LARGER THAN EPS'
         IFL=1
         RETURN
      ENDIF
!----- up here, all left eigenvalues are large enough, the number is M
!      parepare the new Hamiltonian in the collective subspace
!-----
!----- DD = u
      do K=1,M
         S=1./SQRT(EN(K))
         do I=1,N
           RR(I,K)=DD(I,K)*S
         enddo !I
      enddo !K
!-------------
      do K=1,M
         do I=1,N
            S=0.
            do L=1,N
               S=S+HH(I,L)*RR(L,K)
            enddo !L
            E1(I)=S
         enddo !I
         do J=K,M
            S=0.
            do L=1,N
               S=S+RR(L,J)*E1(L)
            enddo ! L
               GG(J,K)=S
         enddo !J 
       enddo !K
!================
  100 format(5x,12f10.5)
     
!---- diagonalization of the new Hamiltonian
! 
      CALL sdiag(NMAX,M,GG,E,GG,Z,1) ! dim(E) = M 
!
!     calculation of the contra-variant components FF
      IF (IS.EQ.0) RETURN
      CALL MAB(NMAX,NMAX,NMAX,N,M,M,RR,GG,FF,1,0)
!
!     calculation of the co-variant components WW (overlaps)
      IF (IS.LT.2) RETURN 
      CALL MAB(NMAX,NMAX,NMAX,N,N,M,NN,FF,WW,1,0)
!
!     BESTIMMUNG DER ORTHOGONALEN WELLENFUNKTIONEN RR
      IF (IS.LT.3) RETURN
!     calculation of the orthogonal wave functions RR
      CALL MAB(NMAX,NMAX,NMAX,N,M,M,DD,GG,RR,1,0)
!-----
      RETURN
!-end-HNDIAG
      END
c_______________________________________________________________________
      subroutine mab(ma,mb,mc,n1,n2,n3,a,b,c,iph,is)
!
c......................................................................
!
      implicit real*8 (a-h,o-z)
!
      dimension a(ma,ma),b(mb,mb),c(mc,mc)
!      dimension a(n1,n2),b(n2,n3),c(n1,n3)
!
      if (is.eq.0) then
         if (iph.lt.0) then
            do i = 1,n1
            do k = 1,n3
               s = 0.0d0
               do l = 1,n2
                  s = s+a(i,l)*b(l,k)
               enddo
               c(i,k) = - s
            enddo
            enddo
         else
            do i = 1,n1
            do k = 1,n3
               s = 0.0d0
               do l = 1,n2
                  s = s+a(i,l)*b(l,k)
               enddo
               c(i,k) = s
            enddo
            enddo
         endif
      else
         if (iph.lt.0) then
            do i = 1,n1
            do k = 1,n3
               s = c(i,k)
               do l = 1,n2
                  s = s-a(i,l)*b(l,k)
               enddo
               c(i,k) = s
            enddo
            enddo
         else
            do i = 1,n1
            do k = 1,n3
               s = c(i,k)
               do l = 1,n2
                  s = s+a(i,l)*b(l,k)
               enddo
               c(i,k)=s
            enddo
            enddo
         endif
      endif
      return
!-end-MAB
      end 
c_______________________________________________________________________
      subroutine sdiag(nmax,n,a,d,x,e,is)

c......................................................................
!
!     A   matrix to be diagonalized
!     D   eigenvalues    
!     X   eigenvectors
!     E   auxiliary field
!     IS = 1  eigenvalues are ordered and major component of X is positiv
!          0  eigenvalues are not ordered            
c......................................................................
      implicit double precision (a-h,o-z)
!
      dimension a(nmax,nmax),x(nmax,nmax),e(n),d(n)
!
      data tol,eps/1.e-32,1.e-10/                           
!
      if (n.eq.1) then
         d(1)=a(1,1)  
         x(1,1)=1.
         return
      endif
!
      do 10 i=1,n 
      do 10 j=1,i 
   10    x(i,j)=a(i,j)
!
!cc   householder-reduktion
      i=n
   15 if (i-2) 200,20,20
   20 l=i-2
      f=x(i,i-1)
      g=f            
      h=0  
      if (l) 31,31,32
   32 do 30 k=1,l
   30 h=h+x(i,k)*x(i,k)
   31 s=h+f*f         
      if (s-tol) 33,34,34              
   33 h=0                             
      goto 100                       
   34 if (h) 100,100,40             
   40 l=l+1                        
      g= dsqrt(s)
      if (f.ge.0.) g=-g        
      h=s-f*g                 
      hi=1.d0/h                
      x(i,i-1)=f-g          
      f=0.0                 
      if (l) 51,51,52     
   52 do 50 j=1,l        
      x(j,i)=x(i,j)*hi  
      s=0.0             
      do 55 k=1,j     
   55 s=s+x(j,k)*x(i,k)                      
      j1=j+1                                
      if (l-j1) 57,58,58                   
   58 do 59 k=j1,l                        
   59 s=s+x(k,j)*x(i,k)                  
   57 e(j)=s*hi                         
   50 f=f+s*x(j,i)                     
   51 f=f*hi*.5d0                      
!                                    
      if (l) 100,100,62             
   62 do 60 j=1,l                  
      s=x(i,j)                    
      e(j)=e(j)-f*s              
      p=e(j)                    
      do 65 k=1,j              
   65 x(j,k)=x(j,k)-s*e(k)-x(i,k)*p        
   60 continue                            
  100 continue                           
      d(i)=h                            
      e(i-1)=g                         
      i=i-1                           
      goto 15            
!            
!cc   Bereitstellen der Transformationmatrix 
  200 d(1)=0.0                               
      e(n)=0.0                              
      b=0.0                                
      f=0.0                               
      do 210 i=1,n                      
      l=i-1                            
      if (d(i).eq.0.) goto 221        
      if (l) 221,221,222             
  222 do 220 j=1,l                  
      s=0.0                         
      do 225 k=1,l                
  225 s=s+x(i,k)*x(k,j)          
      do 226 k=1,l              
  226 x(k,j)=x(k,j)-s*x(k,i)   
  220 continue                
  221 d(i)=x(i,i)            
      x(i,i)=1              
      if (l) 210,210,232   
  232 do 230 j=1,l        
      x(i,j)=0.0          
  230 x(j,i)=0.0         
  210 continue         
!
!cc   Diagonalisieren der Tri-Diagonal-Matrix
      DO 300 L=1,N                     
      h=eps*( abs(d(l))+ abs(e(l)))
      if (h.gt.b) b=h             
!
!cc   Test fuer Splitting        
      do 310 j=l,n              
      if ( abs(e(j)).le.b) goto 320
  310 continue                 
!
!cc   test fuer konvergenz    
  320 if (j.eq.l) goto 300   
  340 p=(d(l+1)-d(l))/(2*e(l))          
      r= dsqrt(p*p+1.d0)
      pr=p+r                           
      if (p.lt.0.) pr=p-r             
      h=d(l)-e(l)/pr                 
      do 350 i=l,n                  
  350 d(i)=d(i)-h                  
      f=f+h                       
!
!cc   QR-transformation          
      p=d(j)                    
      c=1.d0                     
      s=0.0                    
      i=j                    
  360 i=i-1                 
      if (i.lt.l) goto 362 
      g=c*e(i)            
      h=c*p              
      if ( abs(p)- abs(e(i))) 363,364,364
  364 c=e(i)/p                          
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*p*r                     
      s=c/r                           
      c=1.d0/r                         
      goto 365                      
  363 c=p/e(i)                     
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*e(i)*r             
      s=1.d0/r                      
      c=c/r                     
  365 p=c*d(i)-s*g             
      d(i+1)=h+s*(c*g+s*d(i)) 
      do 368 k=1,n           
         h=x(k,i+1)            
         x(k,i+1)=x(k,i)*s+h*c
  368    x(k,i)=x(k,i)*c-h*s 
      goto 360           
  362 e(l)=s*p          
      d(l)=c*p         
      if ( abs(e(l)).gt.b) goto 340
!
!cc   konvergenz      
  300 d(l)=d(l)+f    
!
      if (is.eq.0) return
!cc   ordnen der eigenwerte    
      do 400 i=1,n            
      k=i                    
      p=d(i)                
      j1=i+1               
      if (j1-n) 401,401,400   
  401 do 410 j=j1,n          
      if (d(j).ge.p) goto 410 
      k=j                    
      p=d(j)                
  410 continue             
  420 if (k.eq.i) goto 400
      d(k)=d(i)          
      d(i)=p            
      do 425 j=1,n     
      p=x(j,i)        
      x(j,i)=x(j,k)  
  425 x(j,k)=p      
  400 continue     
!                 
!     signum
      do k = 1,n
         s = 0.0d0
         do i = 1,n
            h = abs(x(i,k))
            if (h.gt.s) then
               s  = h
               im = i
            endif
         enddo   ! i
         if (x(im,k).lt.0.0d0) then
            do i = 1,n
               x(i,k) = - x(i,k)
       	    enddo
         endif
      enddo   ! k
! 
      return
!-end-SDIAG
      end 
c______________________________________________________________________________
      double precision function wignei(j1,j2,j3,m1,m2,m3)

c......................................................................
C
C     Calculates Wigner-coefficients for integer j- and m-values
C
c......................................................................
      implicit double precision (a-h,o-z)
      include 'dic.par'
C
      common /gfvfi / fi(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfviv / iv(-igfv:igfv)
c
      wignei=0.d0           
      if (m1+m2+m3.ne.0) return
      i0 = j1+j2+j3+1
      if (igfv.lt.i0) stop 'in wignei: igfv too small'
      i1 = j1+j2-j3
      i2 = j1-m1
      i3 = j2+m2
      i4 = j3-j2+m1
      i5 = j3-j1-m2
      n2 = min0(i1,i2,i3) 
      n1 = max0(0,-i4,-i5)
      if (n1.gt.n2) return
      do 11 n = n1,n2
   11 wignei = wignei + iv(n)*
     &         fi(n)*fi(i1-n)*fi(i2-n)*fi(i3-n)*fi(i4+n)*fi(i5+n)
      wignei = iv(i2+i3)*wignei*wfi(i0)*wf(i1)*wf(i2+i4)*wf(i3+i5)*
     &         wf(i2)*wf(j1+m1)*wf(j2-m2)*wf(i3)*wf(j3-m3)*wf(j3+m3)
C
      return
c-end-WIGNEI
      end  
!==================================================================
      double precision function djmk(j,m,k,cosbet,is)
!================================================================== 
!     Calculates the Wigner-Functions  d_mm'^j (beta) 
!     IS = 0: integer values for j = J,  m = M, m' = K
!
!     IS = 1: half integer valus for j = J-1/2 , m = M-1/2, m' = K-1/2
!
!     COSBET = cos(beta)
!
!     The Wigner-Functions and the phases are defined as in Edmonds.
!
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'dic.par'
      common /gfviv / iv(-igfv:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvfi / fi(0:igfv)
C
      djmk = 0.d0
      if (iabs(m).gt.j.or.iabs(k).gt.j) return
      if (abs(cosbet).gt.1.d0)
     &   stop 'in DJMK: cos(beta) is larger than one'
      if (cosbet.eq.1.d0) then
         if (m.eq.k) djmk = 1.d0
         return
      endif
      if (cosbet.eq.-1.d0) then
         if (m.eq.is-k) djmk = dfloat(iv(j-m))
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
      X  = IV(JMM+IA) * C**(2*IA+MPK) * S**(JMM+JMK-2*IA)
     &     *FI(JMM-IA)*FI(JMK-IA)*FI(MPK+IA)*FI(IA)

      do i = ia,ie
         DJMK = DJMK + X
         X = -X*CS*(JMM-I)*(JMK-I)/((I+1)*(MPK+I+1))
      enddo
      DJMK = DJMK*WF(JPM)*WF(JMM)*WF(JPK)*WF(JMK)
C
      return
C-end-djmk
      end
c______________________________________________________________________________
      function ssum (n,a)

      double precision zero,a(n),ssum

      parameter (zero=0.0d0)

      ssum = zero
      do i=1,n
        ssum = ssum + a(i)
      enddo

      return
      end function ssum
c..............................................................................
c______________________________________________________________________________
      function fm2wu(amass,lambda)
c
c     transformation from e^2fm^(2*lambda) to W.u.

      double precision fm2wu,el,fac,amass
      integer lambda 
      
      fm2wu =(3.d0/(3.d0+lambda))**2
     &     *amass**(2*lambda/3.d0)
     &     *(1.2d0)**(2*lambda)/(4.d0*3.1415926) 

      return
      end function fm2wu 
c..............................................................................
c______________________________________________________________________________
      subroutine kernel_hoa(iq1,iq2,jprma)

c..............................................................................
c     Kernels in HO approximation (Ring & Schuck, The nuclear many-body problem).
c     page 409
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'dic.par'
! 
      complex*16 nkk,hkk,xnkk 
      complex*16 qp,qcp,qpred,qcpred,q0p,q0pc
      complex*16 q3p,q3cp,q3pred,q3cpred,q1pred,q1cpred
      complex*16 q3n,q3cn,q3nred,q3cnred
      complex*16 qn,qcn,qnred,qcnred,q0n,q0nc
      character*58 elem

      common /big   / qp    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcp   (maxq,maxq,0:jmax,0:jmax,0:jmax), 
     2                qpred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0p   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0pc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                xnkk  (maxq,maxq,0:jmax,0:jmax,0:jmax,2)  
      common /nme_oct/ q3p    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                 q3cp   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                 q3pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                 q3cpred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                 q1pred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                 q1cpred(maxq,maxq,0:jmax,0:jmax,0:jmax) 
      common /neutron/qn    (maxq,maxq,0:jmax,0:jmax,0:jmax),
     1                qcn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     2                qnred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     3                qcnred(maxq,maxq,0:jmax,0:jmax,0:jmax),
     4                q0n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     5                q0nc  (maxq,maxq,0:jmax,0:jmax,0:jmax),
     6                q3n   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     7                q3cn   (maxq,maxq,0:jmax,0:jmax,0:jmax),
     8                q3nred (maxq,maxq,0:jmax,0:jmax,0:jmax),
     9                q3cnred(maxq,maxq,0:jmax,0:jmax,0:jmax)       
          
           
      common /kerne/  nkk(0:jmax,0:jmax,0:jmax),
     &                hkk(0:jmax,0:jmax,0:jmax) 
      common /gcmhov/ njqqkk(0:jmax,maxq*jmax,maxq*jmax),
     &                hjqqkk(0:jmax,maxq*jmax,maxq*jmax)
      common /gcmmes/ gcmbet(100),gcmgam(100),gcmr2(100),step,maxmp
      common /gfviv / iv(-igfv:igfv)
      common /mathco/ zero,one,two,half,third,pi 
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp,lrpa  
      
100   format('(',2f6.3,')',4f12.6)   
       write(*,*) '.....Kernels are determined with HO appro......'
      if(jprma.ne.0) stop ' In HO approximation, jprma should be zero'
c     ........................... initial
        jproj = 0
        k1=0
        k2=0
        s2 = 0.002342     ! s^2
        h2d2m = 4.496     ! h^2/2M
        cmomeg2d2 = 311.d0    ! M*h*ome^2/2
        
        betac1  = gcmbet(iq1) 
        gammac1 = gcmgam(iq1) 
        betac2  = gcmbet(iq2)  !  beta3
        gammac2 = gcmgam(iq2)  !  beta3'  
       if(abs(betac1-beta2).gt.1.d-2)  
     & stop ' In HO approximation, beta2 is not equal'
        dgam  = gammac2-gammac1
        gam12 = gammac1+gammac2
        nkk(jproj,k1,k2) = dexp(-dgam**2/(4.d0*s2))
        hkk(jproj,k1,k2) =  h2d2m*(1.d0/(2*s2)-dgam**2/(4*s2))
     &                    + cmomeg2d2*(s2/2.d0 + gam12**2/4.d0)
        hkk(jproj,k1,k2)= (hkk(jproj,k1,k2)-2596.45)*nkk(jproj,k1,k2)
        
       write(*,100)gammac1,gammac2,nkk(jproj,k1,k2),hkk(jproj,k1,k2)
      return
      end 


            character(500) function find_file(directory_name,file_name)
            character(*) :: file_name,directory_name
            character(1000) :: path
            integer :: istart,iend,i
            logical :: isthere,last_chance

! file_name = file_name)
! directory_name = adjustl(directory_name)

             call GETENV(trim(directory_name),path)
             path=adjustl(path)

             istart = 1
              do while (.true.)
              i = istart
       do while (path(i:i).ne.':')
                iend = i
                 i = i+ 1
                if (path(i:i) == ' ') then
                 i = 0
                 exit
                end if
                end do
        inquire(file=path(istart:iend)//'/'//trim(file_name),
     &   exist=isthere)
        if (isthere) then
         find_file = path(istart:iend)//'/'
         return
        else if (i == 0) then
        inquire(file='./'//file_name,exist=last_chance)
        if (last_chance) then
           find_file = './'
           return
           else
           PRINT*, 'FILE NOT FOUND: ',trim(adjustl(file_name))
           print*, 'CHECK ENVIRONMENT VARIABLE: ',
     &      trim(adjustl(directory_name))
           STOP
           end if
        end if
         istart = iend+2
        end do
        end function find_file

