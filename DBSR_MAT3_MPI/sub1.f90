!======================================================================
      Subroutine SUB1
!======================================================================
!     run calculations for one partial wave
!----------------------------------------------------------------------
      Use dbsr_mat                      
      Use c_data           
      Use DBS_integrals, only: memory_DBS_integrals

      Implicit none
      Integer :: i, ich,jch, k, it,jt, is,js, nelc_core
      Real(8) :: C,t0,t1,t2
      Integer, external :: Ifind_channel_jj, no_ic_jj

      Call CPU_time(t0)

! ... read configuration expansion and orbitals information: 

      Call Read_data

! ... L-integrals for bound orbitals:

      Call Gen_dhl_core(ncore,mbreit,0)

! ... initialize arrays:

      if(allocated(CBUF)) Deallocate(CBUF,itb,jtb,intb,idfb)
      Allocate(CBUF(mcbuf), itb(mcbuf), jtb(mcbuf), intb(mcbuf), idfb(mcbuf))
      mem_buffer = 24.d0 * mcbuf /(1024*1024)

      if(allocated(ich_state)) Deallocate(ich_state)
      Allocate(ich_state(ncfg))
      Do is=1,ncfg;   ich_state(is) = Ifind_channel_jj(is);  End do
      memory_conf_jj = memory_conf_jj + 4.d0*ncfg/(1024*1024)

      if(allocated(no_state)) Deallocate(no_state)
      Allocate(no_state(ncfg))
      Do is=1,ncfg;   no_state(is) = no_ic_jj(is);  End do
      memory_conf_jj = memory_conf_jj + 4.d0*ncfg/(1024*1024)

      Call CPU_time(t1)
      write(pri,'(/a,T40,f10.2,a)') 'Preparations:',(t1-t0)/60,' min'

!----------------------------------------------------------------------
    1 Continue

      Call CPU_time(t1)

      idiag=0;  Call Alloc_dbsr_matrix

      write(pri,'(/a/  )') 'Main dimensions in dbsr_matrix module:'
      write(pri,'(a,i8,a)') 'nch    = ',nch,  '  - number of channels'
      write(pri,'(a,i8,a)') 'npert  = ',npert,'  - number of perturbers'
      write(pri,'(a,i8,a)') 'ncp    = ',ncp,  '  - number of perturber config.s'
      write(pri,'(a,i8,a)') 'ms     = ',ms,   '  - number of splines' 
      mhm = nch*ms+npert
      write(pri,'(a,i8,a)') 'mhm    = ',mhm,  '  - matrix dimension' 

      write(pri,'(/a/)') 'Main memory consumings:'
      write(pri,'(a,T40,f10.2,a)') 'memory of conf_jj module:', memory_conf_jj,' Mb'
      write(pri,'(a,T40,f10.2,a)') 'memory of DBS_gauss module:', memory_DBS_gauss,' Mb'
      write(pri,'(a,T40,f10.2,a)') 'memory of DBS_orbitals:', memory_DBS_orbitals,' Mb'
      write(pri,'(a,T40,f10.2,a)') 'memory of DBS_integrals:', memory_DBS_integrals,' Mb'
      write(pri,'(a,T40,f10.2,a)') 'memory of cdata module:', mem_cdata,' Mb'
      write(pri,'(a,T40,f10.2,a)') 'memory of bufer:', mem_buffer,' Mb'
      write(pri,'(a,T40,f10.2,a)') 'matrix memory:', mem_mat,' Mb'
      C = memory_DBS_gauss + memory_DBS_orbitals + memory_DBS_integrals + &
          mem_cdata + mem_buffer + mem_mat + memory_conf_jj
      write(pri,'(a,T40,f10.2,a)') 'total_estimations:', C,' Mb'

! ... update <.|p> vectors:

      Call Get_v_ch(1)

!-----------------------------------------------------------------------
! ... full overlap matrix:

      icase=0; Call Alloc_c_data(ntype_O,0,mk,mblock,nblock,kblock,eps_c) 
      Do ich=1,nch;  Call UPDATE_HX(ich,ich,fppqq,1.d0);  End do

                     Call State_res

! ... save overlaps diagonal blocks:

      Do ich=1,nch; i=icc(ich,ich); diag(:,:,ich)=hch(:,:,i); End do

      Call CPU_time(t2)
      write(pri,'(/a,T40,f10.2,a)') 'O-integrals:',(t2-t1)/60,' min'

! ... diagonal blocks for the Hamiltonian matrix:

      Call CPU_time(t1)

      idiag = 1

! ... core-energy shift for diagonal blocks:

      Do ich=1,nch; i=icc(ich,ich); hch(:,:,i)=hch(:,:,i)*Ecore; End do

! ... L-integrals:

      icase=1; Call Alloc_c_data(ntype_L,0,0,mblock,nblock,kblock,eps_c) 
      Call State_res   

      Call CPU_time(t2)
      write(pri,'(/a,T40,f10.2,a)') 'Diag. L-integrals:',(t2-t1)/60,' min'

! ... R-integrals:

      Call CPU_time(t1)

      icase=2; Call Alloc_c_data(ntype_R,0,mk,mblock,nblock,kblock,eps_c) 
      Call State_res  

      Call CPU_time(t2)
      write(pri,'(/a,T40,f10.2,a)') 'Diag. R-integrals:',(t2-t1)/60,' min'

! ... S-integrals:

      if(mbreit.eq.1) then
       Call CPU_time(t1)
       icase=3; Call Alloc_c_data(ntype_S,0,mk+1,mblock,nblock,kblock,eps_c) 
       Call State_res  
       Call CPU_time(t2)
       write(pri,'(/a,T40,f10.2,a)') 'Diag. S-integrals:',(t2-t1)/60,' min'
      end if

! ... target energy part:

      Do ich = 1,nch;  it=iptar(ich); C=Etarg(it)-Ecore
       Call UPDATE_HX(ich,ich,fppqq,C)
      End do

! ... orthogonal conditions:

      Call DBS_ORTH 

! ... diagonalize the diagonal blocks:

      Call CPU_time(t1);  Call Diag_channels;  Call CPU_time(t2)

      write(pri,'(/a,T40,f10.2,a)') 'diagonalization:',(t2-t1)/60,' min '

      nsol = SUM(ipsol)
      write(pri,'(/72("-"))')
      write(pri,'(a,i6,a)') 'nsol =',nsol,'  -  number of channel solutions (new basis)'
      write(pri,'(72("-")/)')

! ... record overlap matrix:      

      if(allocated(eval)) Deallocate(eval); Allocate(eval(nsol))
      k = 0
      Do ich = 1,nch; Do is = 1,ipsol(ich)
       k=k+1; eval(k) = diag(is,ms,ich) 
      End do; End do

      if(allocated(jpsol)) Deallocate(jpsol); Allocate(jpsol(0:nch))
      jpsol = ipsol
      Do i=1,nch;  jpsol(i)=jpsol(i-1)+jpsol(i); End do

      rewind(nui)
      write(nui) ns,nch,npert,nsp,nsq    
      write(nui) nsol
      write(nui) jpsol
      write(nui) eval
      write(nui) (((diag(i,js,ich),i=1,ms),js=1,ipsol(ich)),ich=1,nch)

! ... tranform overlap matrix:      

      Call CPU_time(t1);  icase=0;  Call transform_matrix

! ... check big overlaps:

      k = 0; Call Check_mat(k) 
      if(k.gt.0) &
      write(*,'(a,a,i3)') 'Checking overlap matrix:','  k =',k
      write(pri,'(/a,f6.3,a,i3)') &
         'Checking overlap matrix: s_ovl =', s_ovl,'   k =',k
      if(k.gt.0) go to 1

      if(allocated(overlaps)) Deallocate(overlaps)

      Call Record_matrix

      Call CPU_time(t2)
      write(pri,'(/a,T40,f10.2,a)') 'diagonal blocks:',(t2-t0)/60,' min '
      write(*  ,'( a,T40,f10.2,a)') 'diagonal blocks:',(t2-t0)/60,' min '

!----------------------------------------------------------------------
! ... non-doagonal blocks:

      Call CPU_time(t1)

      idiag = -1

! ... recalculate the overlap non-diagonal blocks:    ???

      hch = 0.d0 
      if(npert.gt.0) then; hcp = 0.d0; hp = 0.d0; end if
      Do it=1,ntarg; Do jt=1,it
        if(it.eq.jt) Cycle; k=it*(it-1)/2+jt; otarg(k)=0.d0
      End do; End do

      icase=0; Call Alloc_c_data(ntype_O,0,mk,mblock,nblock,kblock,eps_c) 
      Call State_res

! ... core-energy shift:

      hch = hch  * Ecore      
      if(npert.gt.0) then; hcp = hcp * Ecore; hp = hp * Ecore; end if   

      Call CPU_time(t2)
      write(pri,'(/a,T40,f10.2,a)') 'None-diag. O-integrals:',(t2-t1)/60,' min'

! ... L-integrals:

      Call CPU_time(t1)

      icase=1; Call Alloc_c_data(ntype_L,0,0,mblock,nblock,kblock,eps_c) 
      Call State_res   

      Call CPU_time(t2)
      write(pri,'(/a,T40,f10.2,a)') 'None-diag. L-integrals:',(t2-t1)/60,' min'

! ... R-integrals:

      Call CPU_time(t1)

      icase=2; Call Alloc_c_data(ntype_R,0,mk,mblock,nblock,kblock,eps_c) 
      Call State_res  

      Call CPU_time(t2)
      write(pri,'(/a,T40,f10.2,a)') 'None-diag. R-integrals:',(t2-t1)/60,' min'

! ... S-integrals:

      if(mbreit.eq.1) then
       Call CPU_time(t1)
       icase=3; Call Alloc_c_data(ntype_S,0,mk+1,mblock,nblock,kblock,eps_c) 
       Call State_res  
       Call CPU_time(t2)
       write(pri,'(/a,T40,f10.2,a)') 'None-diag. S-integrals:',(t2-t1)/60,' min'
      end if

! ... orthogonal conditions:

      Call DBS_ORTH 
  
! ... check the diagonal asymptotic coef.s:

      nelc_core=0; Do i=1,nclosed; nelc_core=nelc_core+jbs(i)+1; End do

      write(pri,'(/a,i3)') &
      'Derivations from 2*nelc for asymptotic coefficients with k=0:'

      Do ich = 1,nch; i = icc(ich,ich)
       acf(i,0) = acf(i,0) + 2*nelc_core
       if(abs(acf(i,0)-2*nelc).lt.0.00001) Cycle
       write(pri,'(i5,2F15.6)') i,acf(i,0)-2*nelc
      End do

      Call CPU_time(t2)

      write(pri,'(/a,T40,f10.2,a)') 'non-diagonal blocks:',(t2-t0)/60,' min '
      write(*  ,'( a,T40,f10.2,a)') 'non-diagonal blocks:',(t2-t0)/60,' min '
!----------------------------------------------------------------------
! ... record interaction matrix: 

      Call CPU_time(t1)  

      Call transform_matrix
      Call Record_matrix

      Call CPU_time(t2)

      write(pri,'(/a,T40,f10.2,a)') 'transform and record matrix:',(t2-t1)/60,' min '
      write(*  ,'( a,T40,f10.2,a)') 'transform and record matrix:',(t2-t1)/60,' min '

      if(allocated(CBUF)) Deallocate(CBUF,itb,jtb,intb,idfb)

!----------------------------------------------------------------------
      Call Pri_orth
!----------------------------------------------------------------------
! ... asymptotic coefficients:

      Do k=0,mk;   Do i=1,iicc
       if(abs(acf(i,k)).lt.0.00001) acf(i,k)=0.d0
      End do; End do

!      acf = 2 * acf        ????

! ... recording asymptotic coefficients: 

      write(nui) mk
      write(nui) (((acf(icc(ich,jch),k),ich=1,nch),jch=1,nch),k=0,mk)  

! ... dbug printing asymptotic coefficients: 

      if(pri_acf.gt.0) then
       write(pri,'(/a)') 'Asymptotic coefficients:'
       Do ich=1,nch; Do jch=1,ich
        write(pri,'(/a,2i5/)') 'ich, jch = ',ich,jch
        write(pri,'(10f10.5)')  acf(icc(ich,jch),:)
       End do;  End do
      end if

!----------------------------------------------------------------------
! ... debug information:

      Call  Target_print

      if(iitar.eq.1) Call Target_new1
      if(iitar.eq.2) Call Target_new2


      if(debug.gt.0) then
       write(pri,*)
       write(pri,'(a,T40,f10.2,a)') 'Ldata time:', t_Ldata/60,' min'
       write(pri,'(a,T40,f10.2,a)') 'Rdata time:', t_Rdata/60,' min'
       write(pri,'(a,T40,f10.2,a)') 'Sdata time:', t_Sdata/60,' min'
       write(pri,*)
      end if

      End Subroutine SUB1 

