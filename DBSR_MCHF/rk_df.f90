
!======================================================================
      Real(8) Function rk_df (i1,j1,i2,j2,k) 
!======================================================================
!     Returns  rk_df (i1, j1; i2, j2), base on the assembling two-electron
!     B-spline integrals (see module DBS_integral)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use df_orbitals 
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8), external :: SUM_AmB

      icallrk = icallrk + 1

      if(i1.eq.ii1.and.i2.eq.ii2.and.j1.eq.jj1.and.j2.eq.jj2.and.k.eq.kk) then
       rk_df = rk; Return
      end if
      if(kk.ne.k) kk=-100

      if(i1.ne.ii1.or.i2.ne.ii2) then
       Call density (ns,ks,dens1,p(1,1,i1),p(1,1,i2),'s')
       Call density (ns,ks,dens3,p(1,2,i1),p(1,2,i2),'s')
       ii1=i1; ii2=i2; kk=-100
       iden1 = iden1 + 1
      end if

      if(j1.ne.jj1.or.j2.ne.jj2) then
       Call density (ns,ks,dens2,p(1,1,j1),p(1,1,j2),'s')
       Call density (ns,ks,dens4,p(1,2,j1),p(1,2,j2),'s')
       jj1=j1; jj2=j2
       iden2 = iden2+ 1
      end if


      if(kk.eq.-100) then
       Call mrk_pppp(k);  Call convol(ns,ks,conv1,dens1,2,'s','s')
       Call mrk_qqqq(k);  Call convol(ns,ks,conv2,dens3,2,'s','s')
       Call mrk_qpqp(k);  Call convol(ns,ks,conv3,dens3,2,'s','s')
       Call mrk_pqpq(k);  Call convol(ns,ks,conv4,dens1,2,'s','s')
       kk = k
       iconv=iconv+1
      end if
 
      rk =  SUM_AmB(ns,ks,conv1,dens2,'s') + &
            SUM_AmB(ns,ks,conv2,dens4,'s') + &
            SUM_AmB(ns,ks,conv3,dens2,'s') + &
            SUM_AmB(ns,ks,conv4,dens4,'s')

      rk_df = rk

      End Function rk_df

!======================================================================
      Real(8) Function rk_pppp_df (i1,j1,i2,j2,k) 
!======================================================================
!     Returns  Rk (P_i1, P_j1; P_i2, P_j2) 
!----------------------------------------------------------------------
      Use DBS_grid,      only: ns,ks
      Use df_orbitals,   only: p  
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_pppp(k)
      Call density (ns,ks,dens,p(1,1,i1),p(1,1,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,1,j1),p(1,1,j2),'s')
      rk_pppp_df = SUM_AmB(ns,ks,conv,dens,'s')
  
      End Function rk_pppp_df


!======================================================================
      Real(8) Function rk_qqqq_df (i1,j1,i2,j2,k) 
!======================================================================
!     Returns  Rk (Q_i1, Q_j1; Q_i2, Q_j2) 
!----------------------------------------------------------------------
      Use DBS_grid,      only: ns,ks
      Use df_orbitals,   only: p  
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_qqqq(k)
      Call density (ns,ks,dens,p(1,2,i1),p(1,2,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,2,j1),p(1,2,j2),'s')
      rk_qqqq_df = SUM_AmB(ns,ks,conv,dens,'s')

      End Function rk_qqqq_df


!======================================================================
      Real(8) Function rk_qpqp_df (i1,j1,i2,j2,k) 
!======================================================================
!     Returns  Rk (Q_i1, P_j1; Q_i2, P_j2) 
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use df_orbitals,      only: p
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_qpqp(k)
      Call density (ns,ks,dens,p(1,2,i1),p(1,2,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,1,j1),p(1,1,j2),'s')
      rk_qpqp_df = SUM_AmB(ns,ks,conv,dens,'s')

      End Function rk_qpqp_df


!======================================================================
      Real(8) Function rk_pqpq_df (i1,j1,i2,j2,k) 
!======================================================================
!     Returns  Rk (P_i1, Q_j1; P_i2, Q_j2) 
!----------------------------------------------------------------------
      Use DBS_grid,      only: ns,ks
      Use df_orbitals,   only: p 
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_pqpq(k)
      Call density (ns,ks,dens,p(1,1,i1),p(1,1,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,2,j1),p(1,2,j2),'s')
      rk_pqpq_df = SUM_AmB(ns,ks,conv,dens,'s')

      End Function rk_pqpq_df


!======================================================================
      Real(8) Function sk_ppqq_df (i1,j1,i2,j2,k)
!======================================================================
!     Returns  S^k (i1, j1; i2, j2), base on the assembling two-electron
!     B-spline integrals (see module DBS_integrals).
!----------------------------------------------------------------------
      Use DBS_grid,          only: ns,ks
      Use df_orbitals,       only: p
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: di(ns,ns),dj(ns,ns),dc(ns,ns)
      Real(8), external :: SUM_AmB

      if(k.lt.0) Stop 'Stop in Sk_ppqq: k < 0'                 ! ???

      Call msk_ppqq(k)
      Call density (ns,ks,di,p(1,1,i1),p(1,2,i2),'n')
      Call convol  (ns,ks,dc,di,2,'n','n')
      Call density (ns,ks,dj,p(1,1,j1),p(1,2,j2),'n')
      sk_ppqq_df = SUM_AmB(ns,ks,dc,dj,'n')

      End Function sk_ppqq_df


!======================================================================
      Real(8) Function sk_pqqp_df (i1,j1,i2,j2,k)
!======================================================================
!     Returns  S^k (i1, j1; i2, j2), base on the assembling two-electron
!     B-spline integrals (see module DBS_integrals).
!----------------------------------------------------------------------
      Use DBS_grid,          only: ns,ks
      Use df_orbitals,       only: p
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: di(ns,ns),dj(ns,ns),dc(ns,ns)
      Real(8), external :: SUM_AmB

      if(k.lt.0) Stop 'Stop in sk_pqqp: k < 0'                   

      Call msk_pqqp(k)
      Call density (ns,ks,di,p(1,1,i1),p(1,2,i2),'n')
      Call convol  (ns,ks,dc,di,2,'n','n')
      Call density (ns,ks,dj,p(1,2,j1),p(1,1,j2),'n')
      sk_pqqp_df = SUM_AmB(ns,ks,dc,dj,'n')

      End Function sk_pqqp_df


!=======================================================================
      Subroutine pri_int (nu,icase,k,i1,i2,i3,i4,S,SS)
!=======================================================================
!     print the integral
!-----------------------------------------------------------------------
      Use df_orbitals

      Implicit none
      Integer, intent(in) :: nu,icase,k,i1,i2,i3,i4
      Real(8), intent(in) :: S,SS
      Character AINT(0:3)/'O','L','R','S'/

      if(icase.eq.1) then
       write(nu,'(a,a,a,a,a,a,f20.10,f10.5)') &        
        AINT(icase),'(',EBS(i1),',',EBS(i2),')=',S,SS  
      else
       write(nu,'(a,i2,a,a,a,a,a,a,a,a,a,f20.10,f10.5)') &
        AINT(icase),k,'(',EBS(i1),',',EBS(i2),';', &    
        EBS(i3),',',EBS(i4),')=',S,SS  
      end if

      End Subroutine pri_int



