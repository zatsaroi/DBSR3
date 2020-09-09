!======================================================================
      Subroutine Check_mult_bnk
!======================================================================
!     Read the configuration list from c-files (units 'nu1' and 'nu2'),
!     define there angular symmetries and compare with existing ones
!     in dataset mult_bnk (unit nub).
!     Prepare the angular arrays.
!----------------------------------------------------------------------
      Use dbsr_mult
      Use conf_jj;     Use symc_list;    Use symt_list

      Implicit none
      Integer :: i,k, ic,jc, ik,jk, it,jt, ij,ijc, iort_c, &
                 ne1,ne2, parity1,parity2
      Integer, external :: Iort_conf_jj, DEF_ij, ITRA, Jparity

! ... read old information, if any:

      if(new) then
       Call Alloc_symc(0);    Call Alloc_symt(0)
      else
       Call Read_symc(nub);   Call Read_symt(nub)
      end if

      write(pri,'(/a)') 'mult_bnk:'
      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt

!---------------------------------------------------------------------
! ... define new symmetries from c-file:

      Call Read_conf_jj(nu1,0,'detect','check'); ne1=ne; parity1=parity
      
      write(pri,'(/a)') 'first set:'
      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt

      Call Read_conf_jj(nu2,0,'detect','check'); ne2=ne; parity2=parity

      write(pri,'(/a)') 'second set:'
      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt

      if(ne1.ne.ne2) Stop 'DBSR_mult: ne1 <> ne2'
      if(ktype.eq.'E'.and.mod(kpol,2).eq.1.and. &
         parity1.eq.parity2) write(*,*) 'DBSR_mult: parity ?'
      if(ktype.eq.'M'.and.mod(kpol,2).eq.0.and. &
         parity1.eq.parity2) write(*,*) 'DBSR_mult: parity ?'

!----------------------------------------------------------------------
! ... define IP_term and IC_term pointers:

      if(Allocated(IT_stat )) Deallocate(IT_stat)
                              Allocate  (IT_stat (nsymt))
      if(Allocated(JP_term )) Deallocate(JP_term)
                              Allocate  (JP_term (nsymt))
      if(Allocated(IC_term1)) Deallocate(IC_term1)
                              Allocate  (IC_term1(nsymc))
      if(Allocated(IC_term2)) Deallocate(IC_term2)
                              Allocate  (IC_term2(nsymc))

      IT_stat = 1
      Do i=1,nsymt
       if(IT_conf(i).gt.0) Cycle
       IT_conf(i) = iabs(IT_conf(i));  IT_stat(i) = 0
      End do

      Call SortI(nsymt,IT_conf,JP_term)

      IC_term1=0
      Do i=1,nsymt
       ic=IT_conf(JP_term(i))
       if(IC_term1(ic).eq.0) IC_term1(ic)=i
                             IC_term2(ic)=i
      End do

!----------------------------------------------------------------------
! ... define if we need additional calculations

      if(allocated(IC_need)) Deallocate(IC_need)
                             Allocate  (IC_need(nsymc))
      if(allocated(JC_need)) Deallocate(JC_need)
                             Allocate  (JC_need(nsymc*(nsymc+1)/2))
      if(allocated(IT_need)) Deallocate(IT_need)
                             Allocate  (IT_need(nsymt))
      if(allocated(IT_done)) Deallocate(IT_done)
                             Allocate  (IT_done(nsymt*(nsymt+1)/2))

      IT_done=0    
      if(.not.new) Call Read_done(nub)
      icalc=.FALSE.; IC_need=0; JC_need=0; IT_need=0
      
       Do ic = 1,nsymc
        Call Get_symc(ic,Jtotal1,no1,nn1,kn1,ln1,jn1,iq1,in1) 
        parity1 = Jparity(no1,ln1,iq1)
       Do jc = ic,nsymc
        Call Get_symc(jc,Jtotal2,no2,nn2,kn2,ln2,jn2,iq2,in2)      
        parity2 = Jparity(no2,ln2,iq2)

         iort_c = Iort_conf_jj(1)
         if(ITRA(Jtotal1,Jtotal2,kpol+kpol).eq.0) iort_c=1
         if(ktype.eq.'E'.and.mod(kpol,2).eq.1.and. &
          parity1.eq.parity2) iort_c=1
         if(ktype.eq.'M'.and.mod(kpol,2).eq.0.and. &
          parity1.eq.parity2) iort_c=1

         k = 0 
         Do ik=IC_term1(ic),IC_term2(ic);  it=JP_term(ik)
         Do jk=IC_term1(jc),IC_term2(jc);  jt=JP_term(jk)
           ij = DEF_ij(it,jt)
           if(IT_done(ij).eq.0) IT_done(ij) = iort_c
           if(IT_stat(it)*IT_stat(jt).eq.0.and.IT_done(ij).eq.0) &    
              IT_done(ij)=-1
           if(IT_done(ij).eq.0) k = 1
         End do
         End do
         if(k.eq.0) Cycle

         ijc=DEF_ij(ic,jc); JC_need(ijc)=1; IC_need(ic)=1; IC_need(jc)=1
         icalc=.TRUE.

       End do
       End do

       ! ... define IT_need:

       Do it = 1,nsymt; k = 0
        Do jt = 1,nsymt
         ij=DEF_ij(it,jt); if(IT_done(ij).ne.0) Cycle; k=1; Exit
        End do         
        IT_need(it) = k
       End do

       if(icalc) then
        write(pri,'(/a)') 'Need of additional calculations --> yes '
       else
        write(pri,'(/a)') 'Need of additional calculations --> no  '
       end if

      Deallocate(IT_stat)

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
      go to 1
    2 Continue

      End Subroutine read_conf_jj
