!======================================================================      
      Subroutine SUB_JJ
!======================================================================
!     define channels orbitals in case of JJ-coupling
!----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none

      Real(8) :: CT
      Integer :: j,l,n,kappa,ip,iset,it,ic,jc,ic1,ic2,JT 
      Integer :: lch_min,lch_max
      Character(5), external :: ELi
      Integer, external :: New_index, Ifind_jjorb, Iadd_cfg_jj, &
                           kappa_lj, ich_del

      nch = 0
      Do it=1,ntarg

       JT = jtarg(it)
       j_min=iabs(Jpar-JT);  if(mod(j_min,2).eq.0) j_min=j_min+1
       j_max=iabs(Jpar+JT);  if(mod(j_max,2).eq.0) j_max=j_max+1
       if(j_min.gt.j_max) Cycle

       Do j = j_min,j_max,2

        lch_min = (j-1)/2; lch_max = (j+1)/2

        Do l=lch_min,lch_max

         if(max_ll.ge.0.and.l.gt.max_ll) Cycle
         if(min_ll.ge.0.and.l.lt.min_ll) Cycle

         if(ptarg(it)*(-1)**l.ne.ipar) Cycle

         kappa = kappa_lj(l,j)

         if(max_ka.ne.0.and.kappa.gt.max_ka) Cycle
         if(min_ka.ne.0.and.kappa.lt.min_ka) Cycle

         iset = New_index(kappa,ksmax,nwf,KEF,IEF) 
         n = ichar('k')

         ip=Ifind_jjorb(n,kappa,iset,2)

         if(ich_del(ilsp,kappa,it).ne.0) Cycle

         if(nch.eq.mch) Call Alloc_channel_jj(mch+imch)
         nch=nch+1; kch(nch)=kappa; iptar(nch)=it
                    lch(nch)=l;  jjch(nch) = j
         ipch(nch) = ip 
         ELC(nch)=ELi(n,kch(nch),iset)

! ... define scattering configurations:      

         ic1=ic_targ(it-1)+1; ic2=ic_targ(it)
         CT = 0.d0

         Do ic=ic1,ic2

          Call Get_cfg_jj(ic);  CT = CT + WC(ic)**2

          no=no+1
          nn(no)=n; kn(no)=kappa; ln(no)=l; iq(no)=1; jn(no)=j; in(no)=iset
          Jshell(no)=j; Vshell(no)=1; Jintra(no)=Jpar; Jtotal=Jintra(no)

          jc = Iadd_cfg_jj('add');   WC(jc) = WC(ic)

         End do    ! over target configuration (ic)

         ipconf(nch)=ncfg-ncfg_targ

!         if(abs(CT-1.d0).gt.c_norm) &
!         write(pri,'(i5,a20,a20,a,F10.3,a)') &
!                  nch, BFT(it),AFT(it),ELC(nch),CT, '  - check normalization'
  
        End do     ! over small l
       End do      ! over small j

      End do       ! over targets
  
      End Subroutine SUB_JJ
