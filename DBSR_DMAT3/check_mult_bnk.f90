!======================================================================
      Subroutine Check_mult_bnk
!======================================================================
!     Reads the configuration list from c-files (unit 'nuc'),
!     defines involved angular symmetries and compares them with the
!     data in the mult-bnk.
!     Prepares the angular-coefficient arrays.
!----------------------------------------------------------------------
      Use dbsr_dmat;  Use conf_jj;  Use symc_list  
                      Use orb_jj;   Use symt_list

      Implicit none
      Character(1) :: kt
      Integer :: i,j,k,it,jt, ne1,ne2, nsymc0,nsymt0
      Integer, external :: DEF_ij

! ... read bank information:

      read(nub) kt,k

      if(kt.ne.ktype) Stop 'DBSR_DMAT: different ktype in the bnk-file'
      if(k .ne.kpol ) Stop 'DBSR_DMAT: different kpol  in the bnk-file'
      Call Read_symc(nub);  nsymc0=nsymc
      Call Read_symt(nub);  nsymt0=nsymt

!      write(pri,'(/a)') 'mult_bnk:'
!      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
!      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt

! ... define angular symmetries from c-file:

      kset1 = 1
      Call Read_conf_jj(nuc1,kset1,'add','check')
      ne1=ne; parity1=parity; ncfg1=ncfg; nwf1=nwf

!      write(pri,'(/a)') 'first set:'
!      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
!      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt
!      write(pri,'(a,i5)') 'number of atomic states   = ',ncfg1

      kset2 = kset1 + maxval(ief(1:nwf))
      Call Read_conf_jj(nuc2,kset2,'add','check')
      ne2=ne; parity2=parity; ncfg2=ncfg-ncfg1; nwf2=nwf-nwf1

!      write(pri,'(/a)') 'second set:'
!      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
!      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt
!      write(pri,'(a,i5)') 'number of atomic states   = ',ncfg2

      if(ne1.ne.ne2) Stop 'DBSR_DMAT: ne1 <> ne2'

      if(ktype(1:1).eq.'E'.and.mod(kpol,2).eq.1.and. &
         parity1.eq.parity2) Stop 'DBSR_mult: parity ?'
      if(ktype(1:1).eq.'M'.and.mod(kpol,2).eq.0.and. &
         parity1.eq.parity2) Stop 'DBSR_mult: parity ?'

      if(nsymc.gt.nsymc0) then
       write(*,*) 'nsymc0, nsymc', nsymc0, nsymc 
       Stop 'DBSR_DMAT: run DBSR_MULT first'
      end if
 
      if(nsymt.gt.nsymt0) then
       write(*,*) 'nsymt0, nsymt', nsymt0, nsymt 
       Stop 'DBSR_DMAT: run DBSR_MULT first'
      end if

!----------------------------------------------------------------------
!                                          define IP_stat and IT_state:

      if(allocated(IS_order )) Deallocate(IS_order)
                               Allocate  (IS_order(ncfg))
      if(allocated(IT_state1)) Deallocate(IT_state1)
                               Allocate  (IT_state1(nsymt))
      if(allocated(IT_state2)) Deallocate(IT_state2)
                               Allocate  (IT_state2(nsymt))

      Call SORTI(ncfg,IS_term,IS_order)

      IT_state1=0; IT_state2=0
      Do i=1,ncfg
       it=IS_term(IS_order(i))
       if(IT_state1(it).eq.0) IT_state1(it)=i; IT_state2(it)=i 
      End do

!----------------------------------------------------------------------
! ... define if we have all needed coefficients:

      if(allocated(IT_done)) Deallocate(IT_done)
                             Allocate(IT_done(nsymt*(nsymt+1)/2))
! ... skip some imformation:

      Call Read_done(nub)

      k = 1
      Do i =       1,ncfg1; it = IS_term(i)
      Do j = ncfg1+1,ncfg ; jt = IS_term(j) 
       if(IT_done(DEF_ij(it,jt)).eq.0) k=0
       if(k.eq.0) Exit
      End do; if(k.eq.0) Exit; End do 
      if(k.eq.0) then
       write(*,*) ' mult_bnk is not full '
       Stop 'DBSR_mult: run DBSR_MULT first'
      end if

      Deallocate(IT_done)

      End Subroutine Check_mult_bnk



!======================================================================
      Subroutine read_conf_jj(muc,kshift,job,check)
!======================================================================
!     read and add configurations to the list "conf_jj"
!     job  =  'add'     -  just add
!          =  'detect'  -  return -ic if exist and 
!          =   others   -  add if not exist 
!     check = 'check'   -  check the configurations for number of 
!                          electrons and parity 
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Character(*), intent(in) :: job, check 
      Integer, intent(in) :: muc,kshift
      Integer, external :: Iadd_cfg_jj
      Integer :: nuc,i,ic

      nuc=iabs(muc); if(muc.gt.0) rewind(nuc)
      if(check.eq.'check') then; ne=0; parity=0; end if
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'*') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ
      Call Decode_cj
      in = in + kshift
      if(check.eq.'check') Call Test_cj
      ic = Iadd_cfg_jj(job)
      if(ic.lt.0) Stop 'Read_conf_jj: repeated states?'
      WC(ic)=0.d0
      i=INDEX(CONFIG,')',BACK=.TRUE.)+1
      if(LEN_TRIM(CONFIG(i:)).ne.0) Read(CONFIG(i:),*) WC(ic)
      go to 1
    2 Continue

      End Subroutine read_conf_jj
