!=======================================================================
      Module hd_core
!=======================================================================
! ... contain B-spline representation for one-electron Dirac hamiltonian
! ... including the interaction with closed core shells
!-----------------------------------------------------------------------
      Use DBS_grid,     only: ns,ks,ms
      Use df_orbitals,  only: nbf,kbs,p

      Implicit none
   
      Real(8), allocatable :: hd (:,:,:)    
      Real(8), allocatable :: hdc(:,:,:)    
      Integer :: nkap = 0
      Integer, allocatable :: kap_list(:)

      Integer :: nbk = 0       ! current number of coefficients
      Integer :: mbk = 0       ! maximum dimension
      Integer :: ibk = 2**10   ! initial dimension
      Integer :: kbk = 0       ! max. multipole index
      Integer :: kcore = 0     ! core size

! ... coefficients and their attributes:

      Real(8), allocatable :: cbk(:)   
      Integer, allocatable :: kr1(:),kr2(:),kr3(:),kr4(:)

      End Module hd_core


!======================================================================
      Subroutine Alloc_hd_core(ncore)
!======================================================================
!     allocation arrays in hd_core Module
!----------------------------------------------------------------------
      Use hd_core

      Implicit none
      Integer, intent(in) :: ncore
      Integer :: i,k,k1,k2
      Integer, external :: ipointer
      
      if(allocated(hd)      ) Deallocate(hd)
      if(allocated(hdc)     ) Deallocate(hdc)
      if(allocated(kap_list)) Deallocate(kap_list)
      nkap = 0
      kcore = 0
      if(ms.le.0) Return

      k1=minval(kbs(1:nbf))
      k2=maxval(kbs(1:nbf))

      nkap=0
      Do k=k1,k2
       if(ipointer(nbf,kbs,k).ne.0) nkap=nkap+1 
      End do

      Allocate(kap_list(nkap))
      Allocate(hd(ms,ms,nkap))
      if(ncore.gt.0) Allocate(hdc(ms,ms,nkap))

      i=0
      Do k=k1,k2
       if(ipointer(nbf,kbs,k).eq.0) Cycle
       i=i+1; kap_list(i)=k 
      End do

      kcore = ncore

      End Subroutine Alloc_hd_core


!======================================================================
      Subroutine Gen_hd_core (ncore,mbreit,kappa)
!======================================================================
!     Generates the HD-matrixes and integrals 
!----------------------------------------------------------------------
      Use hd_core 

      Implicit none
      Integer, intent(in) :: ncore,mbreit,kappa 
      Integer :: k,ik

      if(nkap.eq.0) Call Alloc_hd_core(ncore)

!... set up hd matrix:

      Do ik = 1,nkap;  k=kap_list(ik)
       if(kappa.ne.0.and.kappa.ne.k) Cycle 
       Call Gen_hd(k,hd(1,1,ik))
      End do   
   
! ... core contribution:

      if(ncore.gt.0) then

       if(nbk.eq.0) Call Add_int_core(ncore,mbreit)

       Do ik = 1,nkap;  k=kap_list(ik)
        if(kappa.ne.0.and.kappa.ne.k) Cycle 
        Call Load_int_core(ik)
       End do

      end if

      End Subroutine Gen_hd_core


!=====================================================================
      Subroutine Add_int_core(ncore,mbreit)
!=====================================================================
!     modify the hd-matrix due to interaction with the closed core
!     Breit corrections only from exchange interaction?
!---------------------------------------------------------------------
      Use hd_core
      
      Implicit none
      Integer, intent(in) :: ncore,mbreit 
      Integer :: i, kc,lc,jc, k,k1,k2, v, kk,ll,jj, ik, int
      Integer, External ::  l_kappa, j_kappa, itra

      Real(8) :: C,CC,CR, S(8)
      Real(8), external :: CJKJ, SMU

      Do ik=1,nkap; kk=kap_list(ik); ll=l_kappa(kk); jj=j_kappa(kk)

      Do i=1,ncore; kc=kbs(i); lc=l_kappa(kc); jc=j_kappa(kc)

! ... direct interaction:

       C = dble(jc+1);  Call Add_hd_coef(0,1,i,ik,C)  

! ... exchange interaction:

       k1=iabs(jc-jj)/2;  k2=iabs(jc+jj)/2
 
       Do k=k1,k2

        C = -CJKJ(jc,k,jj)**2 / (jj+1); 
        if(C.eq.0.d0) Cycle
 
        CR = C; if(mod(lc+k+ll,2).eq.1.or.itra(lc,k,ll).eq.0) CR = 0.d0

        if(CR.ne.0.d0)  Call Add_hd_coef(k,4,i,ik,C)  

        if(mbreit.eq.0) Cycle

        Do v=k-1,k+1; if(v.lt.0) Cycle

         if(mod(lc+v+ll,2).ne.1) Cycle
         if(SMU(kk,kc,kc,kk,k,v,S).eq.0.d0) Cycle

         int=53;  CC = C*(S(1)+S(4))
                  Call Add_hd_coef(v,int,i,ik,CC)  ! Sk(Pi Pc; Qc Qi)
         int=54;  CC = C*(S(2)+S(3))
                  Call Add_hd_coef(v,int,i,ik,CC)  ! Sk(Pc Pi; Qi Qc)
         int=63;  CC = C*(S(5)+S(6))
                  Call Add_hd_coef(v,int,i,ik,CC)  ! Sk(Pi Qc; Qc Pi)
         int=64;  CC = C*(S(7)+S(8))
                  Call Add_hd_coef(v,int,i,ik,CC)  ! Sk(Pc Qi; Qj Pc)
        
        End do   ! over v
       End do    ! over k
      End do     ! over i
      End do     ! over ik

      kbk = maxval(kr1(1:nbk))


      End Subroutine add_int_core


!=====================================================================
      Subroutine Load_int_core(ik)
!=====================================================================
      Use hd_core

      Implicit none
      Integer :: j, ii,jj, k, met, ik,io
      Real(8) :: x(ns,ns),xx(ns,ns),d(ns,ks),dd(ns,ks)

      hdc(:,:,ik) = hd(:,:,ik)

      Do k = 0,kbk

      Do ii=1,2; Do jj=1,2

      Select Case (10*ii+jj)      
       Case(11);  Call mrk_pppp(k)
       Case(12);  Call mrk_pqpq(k)
       Case(21);  Call mrk_qpqp(k)
       Case(22);  Call mrk_qqqq(k)
      End Select

! ... add direct integrals:

      if(k.eq.0) then
       dd = 0.d0
       Do j = 1,nbk   
        if(kr4(j).ne.ik) Cycle
        if(cbk(j).eq.0.d0) Cycle
        if(kr1(j).ne.0) Cycle
        if(kr2(j).ne.1) Cycle
        io = kr3(j)
        Call density(ns,ks,d,p(1,jj,io),p(1,jj,io),'s')
        dd = dd + cbk(j)*d
       End do
    
       Call Convol(ns,ks,d,dd,1,'s','s')
       Call UPDATE_HS(ms,hdc(1,1,ik),ii,ii,ns,ks,d,'s')          
      end if

! ... exchange contribution:

      xx = 0.d0; met=0
      Do j = 1,nbk   
       if(kr4(j).ne.ik) Cycle
       if(cbk(j).eq.0.d0) Cycle
       if(kr1(j).ne.k) Cycle
       if(kr2(j).ne.4) Cycle
       io = kr3(j)
       Call density(ns,ks,x,p(1,ii,io),p(1,jj,io),'x')
       xx = xx + cbk(j)*x
       met=1
      End do

      if(met.eq.1) then
       Call Convol(ns,ks,x,xx,4,'s','s')
       Call UPDATE_HS(ms,hdc(1,1,ik),ii,jj,ns,ks,x,'x')          
      end if

! ... Breit interaction: not added yet
       
      End do; End do
      End do  ! over k

      End Subroutine Load_int_core        


!======================================================================
      Subroutine Alloc_hd_coefs(m)
!======================================================================
! ... allocate, deallocate, or reallocate the data
!----------------------------------------------------------------------
      Use hd_core

      Implicit none
      Integer :: m
      Integer, allocatable :: iar(:)
      Real(8), allocatable :: rar(:)

      if(m.le.0) then
       if(allocated(cbk)) Deallocate (cbk,kr1,kr2,kr3,kr4)
       nbk = 0; mbk = 0 
      elseif(.not.allocated(cbk)) then
       mbk = m; nbk = 0
       Allocate(cbk(mbk),kr1(mbk),kr2(mbk),kr3(mbk),kr4(mbk))
      elseif(m.le.mbk) then
       Return
      elseif(nbk.eq.0) then
       Deallocate (cbk,kr1,kr2,kr3,kr4)
       mbk = m
       Allocate(cbk(mbk),kr1(mbk),kr2(mbk),kr3(mbk),kr4(mbk))
      else
       Allocate(rar(nbk))
       rar=cbk; Deallocate(cbk); Allocate(cbk(m)); cbk(1:nbk)=rar
       Deallocate(rar)
       Allocate(iar(nbk))
       iar=kr1; Deallocate(kr1); Allocate(kr1(m)); kr1(1:nbk)=iar
       iar=kr2; Deallocate(kr2); Allocate(kr2(m)); kr2(1:nbk)=iar
       iar=kr3; Deallocate(kr3); Allocate(kr3(m)); kr3(1:nbk)=iar
       iar=kr4; Deallocate(kr4); Allocate(kr4(m)); kr4(1:nbk)=iar
       Deallocate(iar)
       mbk = m
      end if

      End Subroutine Alloc_hd_coefs


!======================================================================
      Subroutine Add_hd_coef(k1,k2,k3,k4,C)
!======================================================================
!     add new data to the list 
!----------------------------------------------------------------------
      Use hd_core

      Implicit none
      Integer, intent(in) :: k1,k2,k3,k4
      Real(8), intent(in) :: C
      Integer :: i,k,l,m

      if(mbk.eq.0) Call Alloc_hd_coefs(ibk)

! ... search position (k) for new integral

      k=1; l=nbk
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (k1.lt.kr1(m)) then;       l = m - 1
      elseif(k1.gt.kr1(m)) then;       k = m + 1
      else
       if    (k2.lt.kr2(m)) then;      l = m - 1
       elseif(k2.gt.kr2(m)) then;      k = m + 1
       else
        if    (k3.lt.kr3(m)) then;     l = m - 1
        elseif(k3.gt.kr3(m)) then;     k = m + 1
        else
         if    (k4.lt.kr4(m)) then;    l = m - 1
         elseif(k4.gt.kr4(m)) then;    k = m + 1
         else
          cbk(m)=cbk(m)+C
          Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue 

! ... shift the rest data up:

      Do i=nbk,k,-1
       m = i + 1
       cbk(m)=cbk(i)
       kr1(m)=kr1(i); kr2(m)=kr2(i); kr3(m)=kr3(i); kr4(m)=kr4(i)
      End do

! ... add new integral:

      cbk(k)=C; kr1(k)=k1; kr2(k)=k2; kr3(k)=k3; kr4(k)=k4; nbk=nbk+1
      if(nbk.eq.mbk) Call Alloc_hd_coefs(mbk+ibk) 

      End Subroutine Add_hd_coef


!======================================================================
      Subroutine Get_hd_core (kappa,dhl)
!======================================================================
!     Generates the HD-matrixes and integrals 
!----------------------------------------------------------------------
      Use hd_core 

      Implicit none
      Integer :: kappa,ik, Ipointer
      Real(8) :: dhl(ms,ms)

      ik = Ipointer(nkap,kap_list,kappa)
      if(ik.eq.0) Stop 'Get_hd_core: hd_core is not initialized'

      if(kcore.eq.0) then
       dhl(:,:)=hd(:,:,ik)
      else
       dhl(:,:)=hdc(:,:,ik)
      end if
      
      End Subroutine Get_hd_core


!======================================================================
      Subroutine Get_hd (kappa,dhl)
!======================================================================
!     Generates the HD-matrixes and integrals 
!----------------------------------------------------------------------
      Use hd_core 

      Implicit none
      Integer :: kappa,ik, Ipointer
      Real(8) :: dhl(ms,ms)

      ik = Ipointer(nkap,kap_list,kappa)
      if(ik.eq.0) Stop 'Get_hd: hd_core is not initialized'

      dhl(:,:)=hd(:,:,ik)
      
      End Subroutine Get_hd


!======================================================================
      Real(8) Function dh_value(i,j)
!======================================================================
!     dhl(i,j) integral
!----------------------------------------------------------------------
      Use hd_core

      Implicit none
      Integer, intent(in) :: i,j
      Real(8) :: v(ms),vi(ms),vj(ms),dhl(ms,ms)

      dh_value = 0.d0
      if(kbs(i).ne.kbs(j)) Return

      Call Get_hd(kbs(i),dhl)  
      Call Get_pv_df(i,vi)
      Call Get_pv_df(j,vj)
      v = MATMUL(dhl,vj)
      dh_value = SUM(vi*v)
    
      End Function dh_value


!======================================================================
      Real(8) Function dhc_value(i,j)
!======================================================================
!     dhl(i,j) integral
!----------------------------------------------------------------------
      Use hd_core

      Implicit none
      Integer, intent(in) :: i,j
      Real(8) :: v(ms),vi(ms),vj(ms),dhl(ms,ms)

      dhc_value = 0.d0
      if(kbs(i).ne.kbs(j)) Return

      Call Get_hd_core(kbs(i),dhl)      !  associate ???
      Call Get_pv_df(i,vi)
      Call Get_pv_df(j,vj)
      v = MATMUL(dhl,vj)
      dhc_value = SUM(vi*v)
    
      End Function dhc_value


