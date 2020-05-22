!======================================================================
!     PROGRAM       D B S R _ C O N F                      version 3
!
!                   C O P Y R I G H T -- 2019
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!                        
!======================================================================
!     Generates configuration lists for DBSR code in JJ coupling
!======================================================================
!
!     INPUT ARGUMENT:  
!
!     c_comp   -->  tolerance for compansation configurations [1.01]
!     max_ll   -->  restriction on orbitals momentums for continuum
!     min_ll        orbitals (default = -1 -> no restrictions)
!     max_ka   -->  restriction on kapa-values for continuum
!     min_ka        orbitals (default = -1 -> no restrictions)
!
!     INPUT FILES:
!
!     dbsr_par      -  parameters of calculations
!     target_jj     -  description of target states and partial waves                        
!     targ_nnn.c    -  c-files for target states
!     pert_nnn.c    -  c-files for perturber if any  
!
!     OUTPUT FILES:
!
!     dbsr_conf.log -  running information
!     cfg.nnn       -  c-file for given partial wave nnn
!
!     orthogonal conditions: 
!
!     By default continuum orbital are supposed to be non-orthogonal 
!     to all bound orbitals in the open subshells; 
!     however, we should check if different channels can generate 
!     the same (N+1)-electron states. In this case, the total overlap 
!     matrix may not be positively defined preventing the Hamilton matrix
!     diagonalization. To avoid such  situation, program assing needed
!     orthogonal conditions, considering only physical orbitals.
!
!     additional orthogonal conditions if any can be given by user 
!     in the dbsr_par file, or in the end of files cfg.###, as
!
!               < nl1 | nl2 >=0

!---------------------------------------------------------------------- 
!======================================================================
! ... positions of configurations in the conf_jj list:
!
! 1.  physical configurations for each target state,
!
!     their pointers ->  ip_phys_conf
!     limit position ->  ncfg_phys
!
! 2.  target configurations   
!
!     their pointers ->  ic_targ  (=ictarg + ncfg_phys)
!     limit position ->  ncfg_targ
!
! 3.  scattering configurations (all or phys):
!
!     limit position ->  ncfg_sct
!
! 4.  perturber configurations if any:
!
!     limit position ->  ncfg_pert
!
! 5.  compensation configurations if any:
!
!     limit position ->  ncfg_comp
!======================================================================
! ... Program flow:
!
! 1.  read target information (file target)
! 2.  read argument (file bsr_par or command line)
! 3.  read target orbitals (file target.bsw)
! 4.  define spectroscopic target configurations (target.nnn + target_orb)
! 5.  read all target configurations  (files target.nnn)
! 6.  read "channel-delete" conditions if any
! 7.  read information about included partial waves (file target); 
!     this information should be prepared by hand - not convivient, 
!     I am sorry, but you may copy these lines from previous calculations) 
!----------------------------------------------------------------------
! ... loop other partial waves:  
!----------------------------------------------------------------------
! 8.  define and record channel configurations in cfg.nnn
! 9.  record pertuber configurations if any
!10.  record imposed orthogonality constrans
!11.  record channels information in file target
!----------------------------------------------------------------------
! ... check if we need additional orth.conditions:  
!----------------------------------------------------------------------
!12.  define scattering channels with phys.target states only
!13.  add pertuber physical configurations and check double-counting
!14.  define orthogonal conditions (routine def_orth_cond) and
!     record all compensation configurations in file "cfg.nnn_comp"
!15.  record "derived" orth.constraints in cfg.nnn
!---- end of partial waves loop   
!...  finish by recording overall information in file target    
!----------------------------------------------------------------------
      Use dbsr_conf

      Implicit real(8) (A-H,O-Z)
      Character(80), allocatable :: Line(:)

      Call dbsr_conf_inf

! ... read input parameters::

      Call Read_arg

! ... read target information:

      Open(nut,file=AF_tar);    Call Read_target_jj(nut)
 
! ... target orbitals:

      Open(nuw,file=AF_wfn,form='UNFORMATTED')
      Call Read_bsw_orb_jj(nuw)
      Close(nuw)
      if(nwf.ne.nwt) Stop 'nwf in target.bsw <> nwt'
      write(pri,'(/a,T33,i8)') 'number of target orbital:',nwf
      max_ll_targ = maxval(lef(1:nwf))

! ... read one-electron radial functions from target.bsw:

!      Call Check_file(AF_wfn)
!      Open(nuw,file=AF_wfn,form='UNFORMATTED')
!      Call Read_bsw_orb_jj(nuw)
!      Close(nuw)

! ... find target physical configurations:

      Call Def_phys_targ

! ... read all target configurations:

      Do it=1,ntarg
       AF = trim(BFT(it))//'.c'
       Open(nuc,file=AF)
       if(it.eq.1) Call Read_core_jj(nuc) 
       Call Read_conf_jj(nuc,0,'add','nocheck')
       Close(nuc)
      End do    
      write(pri,'(/a,T33,i8)') &
       'number of target configurations:',ncfg-ncfg_phys
      if(ncfg-ncfg_phys.ne.nct) Stop 'ncfg_target <> nct from target'
      ncfg_targ=ncfg; lcfg_targ=lcfg; nwf_targ=nwf

!      Call Test_ac  ???

! ... read partial waves:

      kpert=0; Call Read_ipar(nut,'kpert',kpert)
      Call Read_ipar(nut,'nlsp',nlsp)
      if(nlsp.le.0) Stop 'dbsr_conf: nlsp <= 0 ? '
      read(nut,*)
      Allocate(line(nlsp))
      Do i=1,nlsp;  read(nut,'(a)') line(i);  End do
      if(kpert.gt.0) then;  Do i=1,kpert+3; read(nut,*); End do; end if

      write(pri,'(/72(''-'')/)') 
      write(pri,'(i3,a)') nlsp,' - nlsp, number of partial waves'

      write(nut,'(80(''-''))')
      write(nut,'(a)') 'channels:'
      write(nut,'(80(''-''))')

! ... do we have the channels to ignore?

      Call Def_del      

!======================================================================
! ... loop other partial waves:
!======================================================================

      max_ch=0;  max_nc=0;  max_wf=0

      Do ilsp=1,nlsp

      AF=line(ilsp); ncp = 0; nwp = 0
      read(AF,*)  Tpar,Jpar,ipar
      if(LEN_trim(AF).gt.12) read(AF(19:),*) AFP,BFP,ncp,nwp

       write(pri,'(/72(''-'')/)') 
       write(pri,'(a,5(a,i4))') 'Partial wave:', &
       ' ilsp =',ilsp,'   2J =',jpar,'  parity =',ipar,'  ncp =',ncp, '  nwp =',nwp

       ncfg = ncfg_targ; lcfg = lcfg_targ; nwf = nwf_targ 

! ... read perturber configurations to get pertuber orbitals if any
! ... (to be sure that set indexes are same as in bsr_prep)

      if(ncp.ne.0)  then
       AF = trim(BFP)//'.c'
       Open(nuc,file=AF,status='OLD')
       Call Read_conf_jj(nuc,0,'add','nocheck')
      end if
      ncfg = ncfg_targ; lcfg = lcfg_targ; nwf_pert=nwf

! ... define channels orbitals:

      ic_targ(0) = ncfg_phys
      ic_targ(1:ntarg) = ictarg(1:ntarg) + ncfg_phys

      Call SUB_JJ

      ncfg_sct = ncfg;  lcfg_sct = lcfg; nwf_sct=nwf
      write(pri,'(/a,T40,i8)') &
       'Number of scattering configutarions: ', ncfg_sct-ncfg_targ

      min_ll_ch=0; if(nch.gt.0) min_ll_ch = minval(lch(1:nch))

! ... add perturber configurations if any:

! ... perturber configurations:

      if(ncp.ne.0)  then
       AF = trim(BFP)//'.c'
       Open(nuc,file=AF,status='OLD')
       Call Read_conf_jj(nuc,0,'add','nocheck')
      end if
      ncfg_pert = ncfg

! ... output scattering configuration:    

      Call Alloc_orthogonality(0)
      Call Prepare_iort_jj                   
      Call Read_orth_jj(nup)

      write(AF,'(a,i3.3)')'cfg.',ilsp
      Open (nuc,file=AF)
      if(kort.gt.0) Call Read_orth_jj(nuc)

      rewind(nuc)
      write(nuc,'(a15,f16.8)') 'Core subshells:',Etarg(1)
      write(nuc,'(a)') core(1:ncore*5) 
      write(nuc,'(a)') 'Peel subshells:' 
      write(nuc,'(20a5)') (ELF(i),i=ncore+1,nwf)
      write(nuc,'(a)') 'CSF(s):' 
      Do ic=ncfg_targ+1,ncfg
        Call Print_conf_jj (nuc,ic,WC(ic))
      End do
      write(nuc,'(a)') '***'

! ... record channel's information:

      write(nut,'(i3,a,2x,a,i4,a,2i8)')  ilsp,'.', ' nch =',nch, &
        '    nc =',ncfg-ncfg_targ, ncfg-ncfg_sct
      write(nut,'(80(''-''))')
      Do ich=1,nch
        write(nut,'(i3,a,2x,a5,2i5,3x,i8)') &
              ich,'.',ELF(ipch(ich)),kch(ich),iptar(ich),ipconf(ich)
      End do
      write(nut,'(80(''-''))')

      max_nc = max(max_nc,ncfg-nct)
      max_wf = max(max_wf,nwf)
      max_ch = max(max_ch,nch)

!======================================================================
! ... orthogonal conditions:   
!======================================================================

      if(max_ll_targ.lt.min_ll_ch) Cycle     

      ncfg=ncfg_targ; lcfg=lcfg_targ; nwf=nwf_pert 
      ic_targ=jc_targ

! ... define scattering channels with phys.target states:

      nch_save = nch 
      nch=0; Call Alloc_channel_jj(imch)
      c_norm = 0.5   !  tollerance for spectroscopic states normalization

      Call SUB_JJ

      if(nch.ne.nch_save) Stop 'nch <> nch_save'

      ncfg_sct=ncfg; lcfg_sct=lcfg; nwf_sct=nwf

! ... output the imposed orth. conditions

      write(nuc,'(/a/)') 'Imposed orth. conditions:'

      Do ich=1,nch; i=ipch(ich)
       Do j=ncore+1,nwf 
        if(KEF(i).ne.KEF(j)) Cycle
        if(IORT(i,j).ne.0) Cycle
        write(nuc,'(a1,a5,a1,a5,a3)') '<',ELF(i),'|',ELF(j),'>=0'
        Call Iadd_orth(i,j,-1)
       End do
      End do
      Close(nuc)

! ... add pertuber physical configurations: 

      if(ncp.gt.0) Call Check_perturber
     
      ncfg_pert=ncfg; lcfg_pert=lcfg; nwf_pert=nwf
 
! ... define orthogonal conditions

      Call Def_orth_cond

! ... analize and record orthogonal conditions for sct. orbitals:

      Call Record_orth
      write(pri,'(/72(''-''))')

      End do ! over partial waves (ilsp)

!-----------------------------------------------------------------------
! ... overall information:

      write(nut,'(a,i8,a)') 'max_ch  = ',max_ch, &
        '  -  max.number of channels'
      write(nut,'(a,i8,a)') 'max_nc  = ',max_nc, &
        '  -  max.number of configurations'
      write(nut,'(a,i8,a)') 'max_wf  = ',max_wf, &
        '  -  max.number of orbitals'
      write(nut,'(72(''-''))')

      End  ! program DBSR_CONF


