!=====================================================================
      Subroutine get_LS_terms
!======================================================================
!     generate all possible LS states from list of configurations 
!--------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: ic,i

      Do ic=1,ncfg

       Call Get_cfg_jj(ic)

       no1 = 1; ln1(1)=ln(1); iq1(1)=iq(1)
       Do i=2,no
        if(nn(i).eq.nn(i-1).and.ln(i).eq.ln(i-1)) then    
         iq1(no1)=iq1(no1)+iq(i)
        else
         no1=no1+1; ln1(no1)=ln(i); iq1(no1)=iq(i)
        end if
       End do

       Call  gen_LS_term (Jtotal,no1,ln1,iq1)

      End do

      End Subroutine get_LS_terms


!======================================================================
      Subroutine gen_LS_term(Jtot,ino,iln,iiq)
!======================================================================
!     exhaustion of shell-terms
!----------------------------------------------------------------------

      Use conf_LS

      Implicit none
      Integer, intent(in) :: Jtot,ino,iln(*),iiq(*)
      Integer :: mt(msh)       ! the number of term in shell i
      Integer :: nt(msh)       ! the term inder consideration

      Integer :: i,i1,i2, ii,IA,IL,IS
      Integer, external :: Iterm_LS

      no = ino; ln(1:no)=iln(1:no); iq(1:no)=iiq(1:no)

      i1=1                     ! i1 - low  limit of shells
      i2=no                    ! i2 - high limit of shells
      Do i=i1,i2
       mt(i) = Iterm_LS(ln(i),iq(i),-1,IA,IL,IS)
      End do


      i=i1                     ! first shell under consideration
      nt(i)=1

    1 ii = Iterm_LS(ln(i),iq(i),nt(i),LS(i,1),LS(i,2),LS(i,3))
      if(i.lt.i2) then
       i=i+1; nt(i)=1; go to 1
      else
       CALL Sum_Iterm(Jtot)
      end if

    2 nt(i)=nt(i)+1
      if(nt(i).gt.mt(i)) then
        if(i.eq.i1) go to 3
        i=i-1; go to 2
        end if
      go to 1

    3 Return

      End  Subroutine Gen_LS_term


!----------------------------------------------------------------------
      Subroutine Sum_Iterm(Jtot)
!----------------------------------------------------------------------
!     exhaustion of intermediate terms
!----------------------------------------------------------------------
      USE conf_LS
 
      Implicit none

      Integer :: IL_min(msh),IL_max(msh),IS_min(msh),IS_max(msh)
      Integer :: Jtot, i,i1,i2,j1,j2


      LS(1,4)=LS(1,2)
      LS(1,5)=LS(1,3)

      if(no.eq.1) then;  CALL Add_term(Jtot); Return; end if

      i1=2                         ! i1 - low  limit
      i2=no                        ! i2 - high limit in array LS(...)

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
         Call Add_term(Jtot)
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

      End  ! Subroutine Sum_Iterm


!======================================================================
      Subroutine Add_term(Jtot)
!======================================================================
      USE conf_LS

      Implicit none
      Integer :: Jtot, Ltot,Stot
      Integer, external :: Itri, Iadd_symc_LS, Iadd_symt_LS

! ... check the total term: 

      Ltot=LS(no,4)
      Stot=LS(no,5)

      if(Itri(Jtot+1,Ltot,Stot).eq.0) Return

! ... add term: 

      iconf = Iadd_symc_LS(Ltot,Stot,no,iq,ln)
      iterm = Iadd_symt_LS(iconf,no,LS)

      End Subroutine Add_term



