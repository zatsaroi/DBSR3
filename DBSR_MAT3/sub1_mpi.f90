!======================================================================
      Subroutine SUB1
!======================================================================
!     drives calculations for one partial wave
!----------------------------------------------------------------------
      Use MPI
      Use dbsr_mat                      
      Use c_data           
      Use DBS_integrals, only: memory_DBS_integrals
      Use radial_overlaps 

      Implicit none
      Integer :: i,j, ich,jch,  k, is,js, it, met, nelc_core
      Real(8) :: C,CO,t0,t1,t2
      Real(8), allocatable :: btarg(:)
      Integer, external :: Ifind_channel_jj, no_ic_jj
      Integer :: status(MPI_STATUS_SIZE)
 
      t0 = MPI_WTIME()

! ... output log-files:

      i=LEN_TRIM(AF_pri)
      if(myid.eq.0) then
        AF_pri(i-2:i)=ALSP;  Open(pri,file=AF_pri)
      else
        if(debug.gt.0.and.myid.le.debug) then
          write(AF_pri(i-2:i+1),'(i4.4)') myid  
          Open(pri,file=AF_pri)
        else
          pri=0
        end if
      end if 
      if(pri.eq.0) pri_coef=0

!-----------------------------------------------------------------------
! ... read configuration expansion, one-electron orbitals,
! ... orthogonal conditions and other related data: 

      if(myid.eq.0) Call Read_data

      t1 = MPI_WTIME()

      Call br_conf_jj
      Call br_orb_jj     
      Call br_channel_jj
      Call br_DBS_orbitals_pq
      Call br_phys_orb(ntarg)
      Call br_radial_overlaps

      t2 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Broadcast:',(t2-t1)/60,' min'

!-----------------------------------------------------------------------
! ... auxiliary arrays:

      if(allocated(CBUF)) Deallocate(CBUF,itb,jtb,intb,idfb)
      Allocate(CBUF(mcbuf),itb(mcbuf), jtb(mcbuf), intb(mcbuf), idfb(mcbuf))
      mem_buffer = 24.d0 * mcbuf /(1024*1024)
      
      if(allocated(ich_state)) Deallocate(ich_state)
      Allocate(ich_state(ncfg))
      Do is=1,ncfg;   ich_state(is) = Ifind_channel_jj(is);  End do
      memory_conf_jj = memory_conf_jj + 4.d0*ncfg/(1024*1024)
      
      if(allocated(no_state)) Deallocate(no_state)
      Allocate(no_state(ncfg))
      Do is=1,ncfg;   no_state(is) = no_ic_jj(is);  End do
      memory_conf_jj = memory_conf_jj + 4.d0*ncfg/(1024*1024)

! ... L-integrals for bound orbitals:

      Call Gen_dhl_core(ncore,mbreit,0)

      t1 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Preparations:',(t1-t0)/60,' min'

! ... memory estimations:

      if(pri.gt.0) then

       write(pri,'(/a/  )') 'Main dimensions in dbsr_matrix module:'
       write(pri,'(a,i8,a)') 'nch    = ',nch,  '  - number of channels'
       write(pri,'(a,i8,a)') 'npert  = ',npert,'  - number of perturbers'
       write(pri,'(a,i8,a)') 'ncp    = ',ncp,  '  - number of perturber config.s'
       write(pri,'(a,i8,a)') 'ms     = ',ms,   '  - number of splines' 
       mhm = nch*ms+npert
       write(pri,'(a,i8,a)') 'mhm    = ',mhm,  '  - matrix dimension' 
       if(myid.gt.0) &
       write(pri,'(a,i8,a)') 'iicc   = ',iicc, '  - number of blocks'    

       write(pri,'(/a/)') 'Main memory consumings:'
       write(pri,'(a,T40,f10.2,a)') 'memory of conf_jj module:', memory_conf_jj,' Mb'
       write(pri,'(a,T40,f10.2,a)') 'memory of orb_jj module:',  memory_orb_jj,' Mb'
       write(pri,'(a,T40,f10.2,a)') 'memory of DBS_gauss module:', memory_DBS_gauss,' Mb'
       write(pri,'(a,T40,f10.2,a)') 'memory of DBS_orbitals:', memory_DBS_orbitals,' Mb'
       write(pri,'(a,T40,f10.2,a)') 'memory of DBS_integrals:', memory_DBS_integrals,' Mb'
       write(pri,'(a,T40,f10.2,a)') 'memory of OBS_integrals:', 16.d0*mobs/(1024*1024),' Mb'

       C = 8.d0*ntarg*(ntarg+1)/2/(1024*1024) 
       write(pri,'(a,T40,f10.2,a)') 'memory of htarg/otarg', C,' Mb'

       C = 7.0*mblock*nblock + 4.0*nblock + (mk+1.0)*(2+nblock)*ntype_R
       mem_cdata = C * 4.0 / (1024 * 1024) 
       write(pri,'(a,T40,f10.2,a)') 'memory of cdata module:', mem_cdata,' Mb'

       write(pri,'(a,T40,f10.2,a)') 'memory of bufer:', mem_buffer,' Mb'

       mem_mat = 4.d0*nch*nch + 8.d0*nch*(nch-1)/2/nprocs*ms*ms +   &
                 4.d0*nch*npert + 8*nch*npert*ms/nprocs +           &
                 4.d0*npert*npert + 8.d0*npert*(npert+1)/2/nprocs + &
                 8.d0*nch*(nch-1)/2/nprocs*(mk+1)
       mem_mat = mem_mat/(1024*1024)
       write(pri,'(a,T40,f10.2,a)') 'matrix memory:', mem_mat,' Mb'

       C = memory_DBS_gauss + memory_DBS_orbitals + memory_DBS_integrals + &
           mem_cdata + mem_buffer + mem_mat + memory_conf_jj + 16.d0*mobs/(1024*1024)
       write(pri,'(/a,T40,f10.2,a)') 'Total_estimations:', C,' Mb'

      end if

!-----------------------------------------------------------------------
! ... first consider diagonal blocks:

      t1 = MPI_WTIME()

   1  idiag = 1;   Call Alloc_dbsr_matrix

      if(allocated(htarg)) Deallocate(htarg)
      Allocate(htarg((ntarg+1)*ntarg/2)); htarg=0.d0
      if(allocated(otarg)) Deallocate(otarg)
      Allocate(otarg((ntarg+1)*ntarg/2)); otarg=0.d0

! ... update <.|p> vectors:

      Call Get_v_ch(1)

! ... calculate overlap matrix (diagonal blocks):

      icase=0; Call Alloc_c_data(ntype_O,0,0,mblock,nblock,kblock,eps_c) 
      Call State_res

! ... B-spline overlaps:

      Do ich=1,nch; if(icc(ich,ich).eq.0) Cycle
       Call UPDATE_HX(ich,ich,fppqq,1.d0)
      End do

! ... save overlaps diagonal blocks:

      if(myid.ne.0) then
       if(allocated(diag)) deallocate(diag)
       Allocate(diag(ms,ms,iicc))
       diag = hch
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      t2 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Diag. overlaps:',(t2-t1)/60,' min'

!----------------------------------------------------------------------
! ... diagonal part of Hamiltonian matrix:

      t1 = MPI_WTIME()

! ... core-energy shift in Hamiltonian matrix (diagonal blocks):

      hch= hch * Ecore      

! ... L-integrals:

      icase=1; Call Alloc_c_data(ntype_L,0,0,mblock,nblock,kblock,eps_c) 
      Call State_res   

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Diag. L-integrals:',(t2-t1)/60,' min'

! ... R-integrals:

      t1 = MPI_WTIME()

      icase=2; Call Alloc_c_data(ntype_R,0,mk,mblock,nblock,kblock,eps_c) 
      Call State_res  

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Diag. R-integrals:',(t2-t1)/60,' min'

! ... S-integrals:

      if(mbreit.eq.1) then
       t1 = MPI_WTIME()
       icase=3; Call Alloc_c_data(ntype_S,0,mk+1,mblock,nblock,kblock,eps_c) 
       Call State_res  
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t2 = MPI_WTIME()
       if(pri.ne.0) &
       write(pri,'(/a,T40,f10.2,a)') 'Diag. S-integrals:',(t2-t1)/60,' min'
      end if

! ... target energy part:

      Do ich = 1,nch; if(icc(ich,ich).eq.0) Cycle; it = iptar(ich)
       C=Etarg(it)-Ecore;  Call UPDATE_HX(ich,ich,fppqq,C)
      End do

! ... orthogonal conditions:

      Call DBS_ORTH

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Diagonal blocks:',(t2-t0)/60,' min'

! ... channel diagonalization:

      t1 = MPI_WTIME()

      Call Diag_channels

! ... collect the results and broadcast them:

      if(allocated(jpsol)) Deallocate(jpsol); Allocate(jpsol(0:nch))
      Call MPI_REDUCE(ipsol,jpsol,nch+1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) ipsol=jpsol
      Call MPI_BCAST(ipsol,nch+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nsol = SUM(ipsol)
      if(pri.gt.0) &
      write(pri,'( /a,i6,a)') 'nsol =',nsol,'  -  number of channel solutions (new basis)'

      if(myid.eq.0) then; if(allocated(diag)) Deallocate(diag); Allocate(diag(ms,ns,nch)); end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      Do i = 1,nch; k=icc(i,i)
       if(k.ne.0) &
        Call MPI_SEND(diag(1:ms,1:ns,k),ms*ns,MPI_DOUBLE_PRECISION, &
                      0, 0, MPI_COMM_WORLD, ierr)
       if(myid.eq.0) &
       Call MPI_RECV(diag(:,:,i),ms*ns,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      End do
     
      if(myid.ne.0) Deallocate(diag)

! ... record diagonal solutions:  

      if(myid.eq.0) then

       if(allocated(eval)) Deallocate(eval); Allocate(eval(nsol))
       k = 0
       Do ich = 1,nch; Do is = 1,ipsol(ich)
        k=k+1; eval(k) = diag(is,ns,ich) 
       End do; End do
    
       jpsol = ipsol;    Do i=1,nch; jpsol(i)=jpsol(i-1)+jpsol(i); End do
    
       rewind(nui)
       write(nui) ns,nch,npert,nsp,nsq    
       write(nui) nsol
       write(nui) jpsol
       write(nui) eval
       write(nui) (((diag(i,js,ich),i=1,ms),js=1,ipsol(ich)),ich=1,nch)
    
       Deallocate(eval)

      end if  ! myid=0

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Diag_channels:',(t2-t1)/60,' min'

! ... collect diagonal asymptotic coefficients:
      
      if(myid.eq.0) then
       if(allocated(bcf)) deallocate(bcf); Allocate(bcf(nch,nch,0:mk)); bcf=0.d0
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      Do i = 1,nch; k=icc(i,i)
       if(k.ne.0) &
       Call MPI_SEND(acf(k,:),mk+1,MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
       if(myid.eq.0) &
       Call MPI_RECV(bcf(i,i,:),mk+1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      End do

!-------------------------------------------------------------------------------------
! ... calculate overlap matrix (non-diagonal blocks):

      idiag = -1;   Call Alloc_dbsr_matrix

      icase=0; Call Alloc_c_data(ntype_O,0,0,mblock,nblock,kblock,eps_c) 
      Call State_res

! ... transform and record overlap matrix:      

      t1 = MPI_WTIME()

      if(myid.eq.0) then
       if(allocated(overlaps)) Deallocate(overlaps) 
       Allocate(overlaps(nch+npert,nch+npert)); overlaps = 0.d0
      end if

      icase = 0;   Call Record_matrix

      if(myid.eq.0) CO = maxval(overlaps)

      Call MPI_BCAST(CO,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      t2 = MPI_WTIME()
      if(pri.ne.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Record overlap matrix:',(t2-t1)/60,' min'

! ... check big overlaps

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if(CO.lt.s_ovl) go to 2

      if(myid.eq.0) Call Check_mat(met)

      Call MPI_BCAST(met,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(pri.gt.0.and.met.gt.0) write(pri,'(/a,i6,a)') &
         'met  =',met,'  -  number of new orth. conditions'
                 
      if(met.gt.0) then
       Call MPI_BCAST(nv_ch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       if(nv_ch.gt.0) then
        if(myid.gt.0) then
         if(allocated(i_ch)) deallocate(i_ch); Allocate(i_ch(nv_ch))
         if(allocated(j_ch)) deallocate(j_ch); Allocate(j_ch(nv_ch))
        end if
        Call MPI_BCAST(i_ch,nv_ch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        Call MPI_BCAST(j_ch,nv_ch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       end if
       go to 1
      end if  ! met > 0
      
    2 Continue
      if(allocated(overlaps)) Deallocate(overlaps) 

!----------------------------------------------------------------------
! ... non-diagonal Hamiltonian matrix:

      t1 = MPI_WTIME()

! ... core-energy shift:

      if(iicc.gt.0) hch = hch * Ecore 
      if(iicb.gt.0) hcp = hcp * Ecore
      if(iibb.gt.0) hp  = hp  * Ecore

! ... L-integrals:

      icase=1; Call Alloc_c_data(ntype_L,0,0,mblock,nblock,kblock,eps_c) 
      Call State_res   

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()
      if(pri.ne.0) &
       write(pri,'(/a,T40,f10.2,a)') 'L-integrals:',(t2-t1)/60,' min'

! ... R-integrals:

      t2 = MPI_WTIME()

      icase=2; Call Alloc_c_data(ntype_R,0,mk,mblock,nblock,kblock,eps_c) 
      Call State_res  

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()
      if(pri.ne.0) &
       write(pri,'(/a,T40,f10.2,a)') 'R-integrals:',(t2-t1)/60,' min'

! ... S-integrals:

      if(mbreit.eq.1) then
       t1 = MPI_WTIME()
       icase=3; Call Alloc_c_data(ntype_S,0,mk+1,mblock,nblock,kblock,eps_c) 
       Call State_res  
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t2 = MPI_WTIME()
       if(pri.ne.0) &
       write(pri,'(/a,T40,f10.2,a)') 'S-integrals:',(t2-t1)/60,' min'
      end if

! ... orthogonal conditions:

      Call DBS_ORTH 

      if(myid.eq.0) Call Pri_orth  

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()

      if(pri.gt.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Non-diagonal blocks:',(t2-t0)/60,' min '

! ... record interaction matrix: 

      t1 = MPI_WTIME() 
      Call Record_matrix

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t2 = MPI_WTIME()

      if(pri.gt.0) &
      write(pri,'(/a,T40,f10.2,a)') 'Record matrix:',(t2-t1)/60,' min '

      Deallocate(hch); if(npert.gt.0) Deallocate(hcp,icb,hp,ibb)

      if(allocated(CBUF)) Deallocate(CBUF,itb,jtb,intb,idfb)

!----------------------------------------------------------------------
! ... asymptotic coefficients:

      t1=MPI_WTIME()

! ... Collect the ACF - matrix:

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Do i = 1,nch-1
       Do j = i+1,nch;  k=icc(i,j)
        if(k.ne.0) &
         Call MPI_SEND(ACF(k,:),mk+1,MPI_DOUBLE_PRECISION, &
                       0, 0, MPI_COMM_WORLD, ierr)
        if(myid.eq.0) &
        Call MPI_RECV(BCF(i,j,:),mk+1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        BCF(j,i,:) = BCF(i,j,:)
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       End do
      End do

      Deallocate(icc,acf)

      t2=MPI_WTIME()

      if(pri.ne.0) &
       write(pri,'(/a,T40,f10.2,a)') 'Collect_ACF:  ',(t2-t1)/60,' min '

! ... asymptotic coefficients:

      if(myid.eq.0) then

      Do i = 1,nch
       Do j = 1,i
        Do k = 0,mk
         if(abs(bcf(i,j,k)).lt.0.00001) bcf(i,j,k) = 0.d0
         bcf(j,i,k) = bcf(i,j,k)
        End do
       End do
      End do

! ... check diagonal asymptotic coefficients for k=0

      nelc_core=0; Do i=1,nclosed; nelc_core=nelc_core+jbs(i)+1; End do
      write(pri,'(/a,i3)') &
      'Derivations from 2*nelc for asymptotic coefficients with k=0:'
      Do i = 1,nch
       bcf(i,i,0) = bcf(i,i,0) + 2*nelc_core
       if(abs(bcf(i,i,0)-2*nelc).lt.0.00001) Cycle
       write(pri,'(i5,2F15.6)') i,bcf(i,i,0)-2*nelc                                                                          
      End do

! ... recording asymptotic coefficients: 

      write(nui) mk
      write(nui) (((bcf(ich,jch,k),ich=1,nch),jch=1,nch),k=0,mk)  

! ... debug printing asymptotic coefficients: 

      if(pri_acf.gt.0) then
       write(pri,'(/a)') 'Asymptotic coefficients:'
       Do ich=1,nch; Do jch=1,ich
        write(pri,'(/a,2i5/)') 'ich, jch = ',ich,jch
        write(pri,'(10f10.5)')  bcf(ich,jch,:)
       End do;  End do
      end if

      end if  ! myid
!----------------------------------------------------------------------
! ... debug information:

      if(check_target.eq.1) then
       if(myid.eq.0)  Allocate(btarg(ntarg*(ntarg+1)/2))
       Call MPI_REDUCE(htarg,btarg,ntarg*(ntarg+1)/2,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                                   MPI_COMM_WORLD,ierr)
       if(myid.eq.0) htarg = btarg
       Call MPI_REDUCE(otarg,btarg,ntarg*(ntarg+1)/2,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                                   MPI_COMM_WORLD,ierr)
       if(myid.eq.0) otarg = btarg
       if(myid.eq.0) Deallocate(btarg)
       if(myid.eq.0.and.check_target.eq.1) Call Target_print
       if(myid.eq.0.and.iitar.eq.1) Call Target_new1
       if(myid.eq.0.and.iitar.eq.2) Call Target_new2
      end if

      End Subroutine SUB1 

