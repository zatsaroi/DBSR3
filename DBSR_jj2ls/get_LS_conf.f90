!=====================================================================
      Subroutine get_LS_conf
!======================================================================
!     generate all possible LS states from list of JJ configurations 
!--------------------------------------------------------------------
      Use jj2ls,    JTOT => JT
      Use conf_jj

      Implicit none
      Integer :: jc,jt,ic,it, i
      Integer, external :: Jadd_cfg_LS

      Do jc=1,ncfg

       Call Get_cfg_jj(jc); jt=iterm
 
       no1 = 1; nn1(1)=nn(1); in1(1)=in(1)
       Do i=2,no
        if(nn(i).eq.nn(i-1).and.ln(i).eq.ln(i-1)) Cycle
        no1=no1+1; nn1(no1)=nn(i); in1(no1)=in(i)
       End do                
       Do it=1,LS_nterms; if(C_term(jt,it).eq.zero) Cycle

        ic = Jadd_cfg_LS(it,no1,nn1,in1)       
        write(nua) jc,ic,C_term(jt,it)

       End do
      End do

      End Subroutine get_LS_conf


!======================================================================
      Integer Function Jadd_cfg_LS(it,no_i,nn_i,kn_i) 
!======================================================================
!     add new CAS to cfg_list with given term 'it'
!----------------------------------------------------------------------

      USE conf_LS

      Implicit none 
      Integer, intent(in) :: no_i,nn_i(*),kn_i(*)

      Integer :: i,ic,it,ip
      Integer, External :: Ifind_nlk

      Jadd_cfg_LS = 0;  if(no_i.le.0) Return
      
      if(mcfg.eq.0) Call Alloc_cfg_LS(icfg)

      Call Get_symt_LS(it,ic,no,LS)
      Call Get_symc_LS(ic,LS(no,4),LS(no,5),no,nn,ln,iq,kn)

      if(no.ne.no_i) then
       write(*,*) 'no,no_i', no,no_i
       Call prj_conf_LS(0,0.d0)
       Stop 'Jadd_cfg_LS: no <> no_i'
      end if

      nn(1:no) = nn_i(1:no)
      kn(1:no) = kn_i(1:no)
      Do i = 1,no; np(i)=Ifind_nlk(nn(i),ln(i),kn(i),2); End do

! ... check if we already have such state:

      Do ic = 1,ncfg
       if(IC_term(ic).ne.it) Cycle
       ip = ip_state(ic); Jadd_cfg_LS = ic
       Do i = 1,no; ip=ip+1
        if(np(i).ne.IP_orb(ip)) then; Jadd_cfg_LS=0; Exit; end if
       End do
       if(Jadd_cfg_LS.ne.0) Return
      End do

      ncfg=ncfg+1
      if(ncfg.eq.mcfg.or.lcfg+no.gt.kcfg) Call Alloc_cfg_LS(mcfg+icfg)

      IC_term(ncfg)=it
      ip_state(ncfg)=lcfg
      Do i=1,no; lcfg=lcfg+1; IP_orb(lcfg)=np(i); End do
      Jadd_cfg_LS = ncfg

      End Function Jadd_cfg_LS


!=====================================================================
      Subroutine write_LS_conf(nuc,core1)
!======================================================================
!     record LS configurations in unit 'nuc'
!--------------------------------------------------------------------
      Use conf_LS

      Character(*) :: core1

      write(nuc,*)
      write(nuc,'(a)') trim(core1)
      Do ic=1,ncfg
       Call Get_cfg_LS(ic)
       Call Incode_c
       write(nuc,'(a)') trim(config)
       write(nuc,'(a)') trim(couple)
      End do
      write(nuc,'(a)') '*'

      End Subroutine write_LS_conf


!=====================================================================
      Subroutine write_LS_conf_c(nuc,core1,C,ipt)
!======================================================================
!     record  LS configurations
!--------------------------------------------------------------------

      Use conf_LS

      Character(*) :: core1
      Real(8) :: C(*)
      Integer :: ipt(*)

      write(nuc,*)
      write(nuc,'(a)') trim(core1)
      Do i=1,ncfg; ic=ipt(i)
       Call Get_cfg_LS(ic)
       Call Incode_c
       write(nuc,'(a,T65,f11.8)') trim(config),C(ic)
       write(nuc,'(a)') trim(couple)
      End do
      write(nuc,'(a)') '*'

      End Subroutine write_LS_conf_c



!=====================================================================
      Subroutine get_LS_core(core_LS)
!======================================================================
!     define LS core
!--------------------------------------------------------------------
      Use conf_jj

      Character(*) :: core_LS
      Character(4), External :: ELF4
 
      core_LS = ' '; if(ncore.eq.0) Return
      core_LS(1:4) = ELF4(nn_core(1),l_core(1),0)
      ii = 1
      Do i=2,ncore
       if(nn_core(i-1).eq.nn_core(i).and.l_core(i-1).eq.l_core(i)) Cycle
       ii=ii+1; j = (ii-1)*4 + 1
       core_LS(j:j+3) = ELF4(nn_core(i),l_core(i),0)
      End do      

      End Subroutine get_LS_core

