!======================================================================
      Subroutine Read_data
!======================================================================
!     read data for given partial wave, klsp
!----------------------------------------------------------------------
      Use dbsr_pol;            Use DBS_grid 
      Use conf_jj;             Use DBS_gauss
      Use orb_jj;              Use DBS_orbitals_pq
      Use channel_jj;          
                               
      Implicit none
      Integer :: i,j,l,k,m,n, ich
      Integer, external :: Ifind_bsorb, Ifind_jjorb, Icheck_file, IORT

! ... set up B-splines:
 
      Call read_knot_dat
      Call alloc_DBS_gauss

! ... read arguments:

      Open(nup,file=AF_par);  Call R_arg(nup);  Close(nup)

! ... read target and channel information:

      Call Check_file(AF_tar)
      Open(nut,file=AF_tar)  
      Call Read_target_jj (nut)                                    
      Call Read_channel_jj(nut,klsp)
      Close(nut)

! ... log-file:

      i=LEN_TRIM(AF_log); AF_log(i-2:i)=ALSP
      Open(pri,file=AF_log)

! ... HEADER:

      write(pri,'(a,i3)')   'DBSR_POL:   klsp = ',klsp
      write(pri,'(a)')      '*********'
!----------------------------------------------------------------------
! ... read obitals and configurations along with expansion coefficients:

      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Open(nuc,file=AF_cfg,status='OLD')

      Call Read_core_jj(nuc)
      Call Read_conf_jj(nuc,0,'add','nocheck')
      Call Prepare_iort_jj(nuc)                   

      write(pri,'(/a)') 'c-file data: '
      write(pri,*)
      write(pri,'(a,i5,a)') 'ncfg   = ',ncfg, '  - number of configurations'
      write(pri,'(a,i5,a)') 'ncore  = ',ncore,'  - core shells'
      write(pri,'(a,i5,a)') 'nwf    = ',nwf,  '  - number of orbitals'
      if(debug.gt.0) &
      write(pri,'(10a8)')  (ELF(i),i=1,nwf)

! ... create the same space in DBS_orbitals:

      nbf=0
      Do i=1,nwf; j=Ifind_bsorb(NEF(i),KEF(i),IEF(i),2); End do
      mbs=0

!----------------------------------------------------------------------
! ... find pointer  orbital --> channel:

      ipef = 0; ipbs = 0
      Do ich = 1,nch
       Call EL_NLJK(ELC(ich),n,k,l,j,i)
       m = Ifind_jjorb(n,k,i,1); ipef(m) = ich
       m = Ifind_bsorb(n,k,i,1); ipbs(m) = ich; ipch(ich)=m 
      End do

      write(pri,'(/a/)')   'channel data: '
      write(pri,'(a,i5,a)')  'jpar   = ',jpar, '  - 2J-value'
      write(pri,'(a,i5,a)')  'ipar   = ',ipar, '  - parity'
      write(pri,'(a,i5,a)')  'nch    = ',nch,  '  - number of channels'
      write(pri,'(a,i5,a)')  'npert  = ',npert,'  - number of perturbers'
      write(pri,'(a,i5,a)')  'ncp    = ',ncp,  '  - pertuber configurations'
      write(pri,'(a,i5,a)')  'nwp    = ',nwp,  '  - pertuber orbitals'

!----------------------------------------------------------------------
! ... read B-spline expantions for bound orbitals:
   
! ... target radial functions:

      Open(nuw, file=AF_bsw, STATUS='OLD', form='UNFORMATTED')
      Call Read_pqbs(nuw)
      Close(nuw)

! ... perturber radial functions:

      if(nwp.gt.0) then
       AF = trim(BFP)//'.bsw'
       Open(nuw, file=AF, STATUS='OLD', form='UNFORMATTED')
       Call Read_pqbs(nuw)
       Close(nuw)
      end if

! ... additional orthogonality:

      if(Icheck_file(AF_ort).ne.0) then
       Open(nuw, file=AF_ort, STATUS='OLD', form='UNFORMATTED')
       Call Read_pqbs(nuw)
       Close(nuw)
      end if

! ... check the correspondence between c- and bsw-files: 

      j = 0
      Do i = 1,nwf
       if(ipbs(i).ne.0) Cycle
       if(mbs(i).eq.0) then
        write(pri,'(a,a)') ' Absent expansion for w.f. ',ebs(i)
        j = j + 1
       end if
      End do
      if(j.gt.0) Stop 'no correspondence between c- and w- files'

! ... define NORT:

      nort = 0
!      Do ich=1,nch;  i=ipch(ich)
!       Do j = 1,nbf; if(kbs(i).ne.kbs(j)) Cycle
!        if(IORT(i,j).ne.0) Cycle
!        nort = nort + 1
!       End do
!      End do
      write(pri,'(/a/)')   'orthogonality constraints: '
      write(pri,'( a,i5,a)') 'nortb  = ',nortb,'  - number of orth. states'

      End Subroutine Read_data
    
