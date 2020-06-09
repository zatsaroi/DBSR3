!======================================================================
      Subroutine SUB_nlLS(C)
!======================================================================
!     get LS configurations 
!----------------------------------------------------------------------
      Use jj2ls; Use conf_jj; Use conf_LS, only: LS

      Implicit none

      Integer :: mt(msh)       ! the number of term in shell i
      Integer :: nt(msh)       ! the term inder consideration
      Integer :: i,kt,ii
      Real(8) :: C,CT,T(msh)
      Integer, External :: Jterm, Itrans, Iterm_LS

      Do i=1,n
       mt(i)  = Itrans(ln1(i),iq1(i),jq1(i),jq2(i),kt1(i),kt2(i),jot(i),-1,CT)
      if(mt(i).eq.0) Return
       moments(i+2*n) = jot(i)+1
       moments(i+5*n) = JT(i)+1
      End do

      nt(1)=1; i=1

    1 kt = Itrans(ln1(i),iq1(i),jq1(i),jq2(i),kt1(i),kt2(i),jot(i),nt(i),CT)
      ii = Iterm_LS(ln1(i),iq1(i),kt,LS(i,1),LS(i,2),LS(i,3))
      if(i.eq.1) then; T(1)=C*CT; else; T(i)=T(i-1)*CT; end if
      if(i.lt.n) then
       i=i+1; nt(i)=1; go to 1
      else
       CALL Sub_Iterm(T(N),Jtotal,n,ln1,iq1)
      end if

    2 nt(i)=nt(i)+1
      if(nt(i).gt.mt(i)) then
       if(i.eq.1) go to 3
       i=i-1; go to 2
      end if
      go to 1

    3 Return

      End  Subroutine SUB_nlLS


!----------------------------------------------------------------------
      Subroutine Sub_Iterm(C,Jtotal,n,ll,qq)
!----------------------------------------------------------------------
!     exhaustion of intermediate terms
!----------------------------------------------------------------------

      USE conf_LS, only: msh,LS

      Implicit none 

      Integer :: Jtotal,n,ll(n),qq(n), i,i1,i2, j1,j2
      Integer :: IL_min(msh),IL_max(msh),IS_min(msh),IS_max(msh)
      Real(8) :: C

      LS(1,4)=LS(1,2)
      LS(1,5)=LS(1,3)

      if(n.eq.1) then;  CALL Add_recup(C,Jtotal,n,ll,qq); Return; end if

      i1=2                         ! i1 - low  limit
      i2=n                         ! i2 - high limit in array LS(...)

      i=i1
    1 j1=i-1; j2=i

      IL_min(i)=IABS(LS(j1,4)-LS(j2,2))+1
      IL_max(i)=     LS(j1,4)+LS(j2,2) -1
      IS_min(i)=IABS(LS(j1,5)-LS(j2,3))+1
      IS_max(i)=     LS(j1,5)+LS(j2,3) -1

      LS(i,4)=IL_min(i)
      LS(i,5)=IS_min(i)

    2 if(i.lt.i2) then
       i=i+1; go to 1
      else
       Call Add_recup(C,Jtotal,n,ll,qq)
      end if

    3 if(LS(i,5).lt.IS_max(i)) then
       LS(i,5)=LS(i,5)+2
       go to 2
      elseif(LS(i,4).lt.IL_max(i)) then
       LS(i,4)=LS(i,4)+2
       LS(i,5)=IS_min(i)
       go to 2
      else
       if(i.eq.i1) go to 4
       i=i-1; go to 3
      end if

    4 Return

      End  ! Subroutine Sub_Iterm


!================================================================
      Subroutine Add_recup(C,Jtot,n,ll,qq)
!================================================================

      USE conf_LS
      USE jj2ls, m => n

      Implicit none

      Integer :: Jtot,n,ll(n),qq(n), i,Ltot,Stot
      Integer, external :: ITRI,Iadd_symc_LS,Iadd_symt_LS 
      Real(8) :: C,Z
      Real(8), external :: ZRECUP

      Ltot = LS(n,4)
      Stot = LS(n,5) 

      if(ITRI(Ltot,Stot,Jtot+1).eq.0) Return

      Do i=1,n
       moments(i    ) = LS(i,2) 
       moments(i+n  ) = LS(i,3)
       moments(i+3*n) = LS(i,4)
       moments(i+4*n) = LS(i,5)
      End do

      if(n.eq.1) then
       Z = one
      else
       Z = ZRECUP(nmom,moments)
      end if

      Z=Z*C;  if(Z.eq.0.d0) Return
 
      iconf = Iadd_symc_LS(Ltot,Stot,n,qq,ll)
      iterm = Iadd_symt_LS(iconf,n,LS)

      if(iterm.gt.LS_nterms) Stop 'Add_recup: iterm > LS_nterms'

      LS_term = iterm

      C_term(JJ_term,LS_term) = C_term(JJ_term,LS_term) + Z

      End Subroutine Add_recup


