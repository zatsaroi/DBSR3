!======================================================================
      Subroutine Read_conf_jj
!======================================================================
!     Read the configuration list from c-file (unit 'nuc').
!     Prepare the angular-coefficienta arrays.
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use orb_jj;  Use symc_list;  Use symt_list

      Implicit none
      Integer :: i,k,it,jt,ij,nsymc0,nsymt0
      Integer, external :: Jdef_ncfg, Iadd_cfg_jj

! ... Check c-file:  name.c or c=...

      AF_cfg = TRIM(name)//'.c'
      Call Read_aarg('c',AF_cfg)
      if(Icheck_file(AF_cfg).eq.0) then
       write(scr,'(/a,a,a/)') 'c-file ',trim(AF_cfg),' is absent'
       write(log,'(/a,a,a/)') 'c-file ',trim(AF_cfg),' is absent'
       Stop 'DBSR_MCHF: stop in Read_conf_jj routine'
      else
       Open(nuc,file=AF_cfg)
       write(log,'(/a,T35,a)') 'Configuration file: ',trim(AF_cfg)
      end if

! ... Check angular-bank file: name.bnk, int_bnk  or bnk=...

      AF_bnk = TRIM(name)//'.bnk'
      Call Read_aarg('bnk',AF_bnk)
      if(Icheck_file(AF_bnk).ne.0) then
       Open(nub,file=AF_bnk,form='UNFORMATTED')
      elseif(Icheck_file(BF_bnk).ne.0) then
       AF_bnk = BF_bnk
       Open(nub,file=AF_bnk,form='UNFORMATTED')
      else 
       write(scr,'(/a/)') 'angular-coef-bank file is absent'
       write(log,'(/a/)') 'angular-coef-bank file is absent'
       Stop 'DBSR_MCHF: stop in Read_conf_jj routine'
      end if

      write(log,'(/a,T35,a)') 'Angular-coefficients file: ',trim(AF_bnk)

! ... define ncfg:

      ncfg = Jdef_ncfg(nuc)
      if(ncfg.eq.0) then
       write(scr,'(/a/)') 'Number of configuration ncfg = 0'
       write(log,'(/a/)') 'Number of configuration ncfg = 0'
       Stop 'DBSR_MCHF: stop in Read_conf_jj routine'
      end if

      Call Read_core_jj(nuc)
      Call Read_peel_orb(nuc)

! ... read bank information:

      Call Read_symc(nub); nsymc0 = nsymc
      Call Read_symt(nub); nsymt0 = nsymt

! ... check symmetries in c-file:

      rewind(nuc); ne=0; parity=0
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'***') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ

      Call Decode_cj;  Call Test_cj

      i = Iadd_cfg_jj('detect')
      if(i.lt.0) Stop 'R_confj: repeated states in c-file?'
      go to 1
    2 Continue     

      if(nsymt.gt.nsymt0.or.nsymc.gt.nsymc0) &
       Stop 'bnk-file is not complete - run dbsr_breit first!'
!----------------------------------------------------------------------
! ... define connection between states and symmetries:

      if(allocated(IS_order )) Deallocate(IS_order)
      if(allocated(IT_state1)) Deallocate(IT_state1)
      if(allocated(IT_state2)) Deallocate(IT_state2)
      Allocate(IS_order(ncfg),IT_state1(nsymt),IT_state2(nsymt))
      IT_state1=0; IT_state2=0

      Call SORTI (ncfg,IS_term,IS_order)

      Do i=1,ncfg
       it=IS_term(IS_order(i))
       if(IT_state1(it).eq.0) IT_state1(it)=i; IT_state2(it)=i 
      End do

!----------------------------------------------------------------------
! ... define if we need additional angular coefficients:

      if(allocated(IT_done)) Deallocate(IT_done)
                             Allocate(IT_done(nsymt*(nsymt+1)/2))
      Call Read_done(nub)

      k = 0
      Do it = 1,nsymt;  if(IT_state1(it).eq.0) Cycle
       Do jt = 1,it;    if(IT_state1(jt).eq.0) Cycle
        ij = (it-1)*it/2+jt
        if(IT_done(ij).gt.0) Cycle
        k=1; Exit
       End do
       if(k.eq.1) Exit
      End do
      if(k.eq.1) Stop 'bnk-file is not complete - run dbsr_breit first!'

      write(log,'(/80(''-''))')
      write(log,*)
      write(log,'(a)')         'c-file data: '
      write(log,*)
      write(log,'(a,T28,i8)')  'number of electrons:',ne
      write(log,'(a,T28,i8)')  'number of configurations:',ncfg
      write(log,'(a,T28,i8)')  'core orbitals:',ncore
      if(LEN_TRIM(core).gt.0)  &
      write(log,'(a,a)')       'core: ',trim(core)
      write(log,'(a,T28,i8)')  'number of peel orbitals:',nwf-ncore
      write(log,*)
      write(log,'(10a8)')  (ELF(i),i=ncore+1,nwf)

      kmin = 0
      kmax = maxval(JEF(1:nwf))

      Call AV_energy_coef(ncore,nwf,lef,jef)

      End Subroutine read_conf_jj


