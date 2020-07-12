!======================================================================
      Subroutine Pre_det_exp
!======================================================================
!     defines the determinant expansions and write the information
!     in scratch files 'nud'
!-----------------------------------------------------------------------
      Use dbsr_breit; Use conf_jj;  Use symc_list; Use symt_list
      Use term_exp;   Use nljm_orbitals

      Implicit none
      Integer, allocatable :: IP_kt(:),IP_det(:,:),JTs(:,:),JTi(:,:),IPs(:,:)
      Real(8), allocatable :: C_det(:), CC_det(:,:)
      Integer :: i,j, k,kt,kdt,ktm, it,it1,it2, JW,JQ
      Real(8) :: S
      Integer, external :: Ndets_jq, Jterm, mj_value

      Call Alloc_nljm(ne,msh)
      rewind(nud)

! ... loop over conf. symmeteries:

      ic_case = 0
      Do ic=1,nsymc;  if(IC_need(ic).eq.0) Cycle

! ... define configuration ic:

       Call Get_symc(ic,Jtotal,no,nn,kn,ln,jn,iq,in)

       k = 1
       Do i=1,no
        ipn(i)=k; k=k+iq(i); md(i)=Ndets_jq(jn(i),iq(i))
        if(mj_max.lt.jn(i)) mj_max=jn(i)
       End do

! ... define relevant angular symmetries:

       it1=IC_term1(ic); it2=IC_term2(ic); ktm=it2-it1+1

       if(allocated(IP_kt)) Deallocate(IP_kt); Allocate(IP_kt(ktm ))
       if(allocated(JTs  )) Deallocate(JTs  ); Allocate(JTs(no,ktm))
       if(allocated(JTi  )) Deallocate(JTi  ); Allocate(JTi(no,ktm))
       if(allocated(IPs  )) Deallocate(IPs  ); Allocate(IPs(no,ktm))
       kt = 0
       Do k =it1,it2
        it=JP_term(k); if(IT_need(it).eq.0) Cycle
        kt = kt + 1; IP_kt(kt) = it
        Call Get_symt(it,ic,no,Jshell,Vshell,Jintra)
        Do i=1,no
         JTs(i,kt) = Jshell(i)
         JTi(i,kt) = Jintra(i)
         IPs(i,kt) = Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ)
        End do
       End do

      if(kt.eq.0) then; IC_need(ic)=0; Cycle; end if

      if(allocated(C_det)) Deallocate(C_det);  Allocate(C_det(kt))

! ... define the det.expansions:

       rewind(nua)
       Call Det_expn_jj;  if(kdt.eq.0) Stop 'Pre_detexp: kdt = 0'

! ... record results (re-write from 'nua' to 'nud'):

       write(nud) ic,kt,kdt
       write(nud) IP_kt(1:kt)
       rewind(nua)
       Allocate(CC_det(kt,kdt),IP_det(ne,kdt))
       Do i = 1,kdt;  read(nua) IP_det(:,i),CC_det(:,i); End do
       write(nud) IP_det
       write(nud) CC_det
       Deallocate(CC_det,IP_det)
       ic_case = ic_case + 1

      End do    ! over ic

      if(allocated(IP_kt)) Deallocate(IP_kt)
      if(allocated(JTs  )) Deallocate(JTs  )
      if(allocated(JTi  )) Deallocate(JTi  )
      if(allocated(IPs  )) Deallocate(IPs  )
      if(allocated(C_det)) Deallocate(C_det)

      if(allocated(mj_orb)) Deallocate(mj_orb)
      Allocate(mj_orb(mj_max+1))
      Do i=1,mj_max+1;  mj_orb(i) = mj_value(i);  End do

Contains

!======================================================================
      Subroutine Det_expn_jj
!======================================================================
!     find all possible determinants for given configurations.
!     The determinants and their coefficients are recoded to unit 'nua'
!
!     Calls: Det_sh_jq, DETC_jq, Clebsh2
!----------------------------------------------------------------------
      Implicit none
      Real(8) :: C
      Real(8), external :: DETC_jq, Clebsh2

      kdt=0; i=1; nd(i)=1
    1 Call DET_sh_jq(jn(i),iq(i),nd(i),MJs(i),Idet(ipn(i)))

      if(i.eq.1) then
       MJi(1)=MJs(1)
      else
       MJi(i) = MJi(i-1)+MJs(i)
      end if

      if(i.lt.no) then;  i=i+1;  nd(i)=1;  go to 1; end if

!--------------------------------------------------------------------
! ... expansion coefficient:

      C_det = 0.d0; k = 0
      Do it=1,kt
       if(MJi(no).ne.JTi(no,it)) Cycle
       C = DETC_jq(jn(1),iq(1),IPs(1,it),nd(1))
       if(C.eq.0.d0) Cycle

       Do j=2,no
        C = C * DETC_jq(jn(j),iq(j),IPs(j,it),nd(j))
        if(C.eq.0.d0) Exit
        C = C * Clebsh2(JTi(j-1,it),MJi(j-1), &
                        JTs(j  ,it),MJs(j  ), &
                        JTi(j  ,it),MJi(j  ))
        if(C.eq.0.d0) Exit
       End do

       if(C.ne.0.d0) then; C_det(it)=C; k=1; end if
      End do

      if(k.ne.0) then
       kdt=kdt+1; write(nua) Idet(1:ne),C_det(1:kt)
      end if

!--------------------------------------------------------------------
    2 nd(i)=nd(i)+1                ! selecting the next case

      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          ! to end
       i=i-1; go to 2
      end if
      go to 1

    3 Continue

      End Subroutine Det_expn_jj

      End Subroutine Pre_det_exp
