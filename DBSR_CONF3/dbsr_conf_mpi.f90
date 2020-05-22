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
      Use MPI

      Use dbsr_conf
      Use internal_file

      Implicit real(8) (A-H,O-Z)

      Character(128) :: line
      Integer, external :: Iadd_line

!---------------------------------------------------------------------- 

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if(myid.eq.0) write(*,'(a,i5)') 'MPI: nprocs = ',nprocs

!---------------------------------------------------------------------- 
! ... read input parameters::

      if(myid.eq.0) Call Read_arg;   Call br_arg

!---------------------------------------------------------------------- 
! ... read target information:

      if(myid.eq.0) then
       Open(nut,file=AF_tar)
       Call Read_target_jj(nut)
      end if 
      Call br_target_jj

!---------------------------------------------------------------------- 
! ... debug output

      if(myid.ne.0) then
       if(debug.ge.myid) then
        write(AF_log,'(a,i5.5)') 'debug.',myid
        Open(pri,file=AF_log)
       else
        pri = 0; debug = 0
       end if
      end if

!---------------------------------------------------------------------- 
! ... read one-electron radial functions from target.bsw:

      if(myid.eq.0) then
       Open(nuw,file=AF_wfn,form='UNFORMATTED')
       Call Read_bsw_orb_jj(nuw)
       Close(nuw)
      end if

!---------------------------------------------------------------------- 
! ... find target physical configurations:

      if(myid.eq.0) Call Def_phys_targ
      Call br_phys_orb(ntarg) 
 
      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(myid.ne.0)  Allocate(ic_targ(0:ntarg),jc_targ(0:ntarg))

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ic_targ(0:ntarg),ntarg+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jc_targ(0:ntarg),ntarg+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ncfg_phys,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lcfg_phys,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!---------------------------------------------------------------------- 
! ... read all target configurations:

      if(myid.eq.0) then
       Do it=1,ntarg
        AF = trim(BFT(it))//'.c'
        Open(nuc,file=AF)
        if(it.eq.1) Call Read_core_jj(nuc) 
        Call Read_conf_jj(nuc,0,'add','nocheck')
        Close(nuc)
       End do    
       ! Call Test_ac     ???
      end if

      Call br_conf_jj

      if(pri.gt.0) &
      write(pri,'(/a,T33,i8)') 'number of target orbital:',nwf

      max_ll_targ = maxval(lef(1:nwf))

      if(pri.gt.0) &
      write(pri,'(/a,T33,i8)')  'number of target configurations:',ncfg-ncfg_phys

      ncfg_targ=ncfg; lcfg_targ=lcfg; nwf_targ=nwf

!---------------------------------------------------------------------- 
! ... read partial waves:

      if(myid.eq.0) then
       kpert = 0;  Call Read_ipar(nut,'kpert',kpert)
       nlsp  = 0;  Call Read_ipar(nut,'nlsp' ,nlsp )
       Call alloc_file(nlsp+kpert)
       read(nut,*)
       Do i=1,nlsp;  read(nut,'(a)') line; j=Iadd_line(line); End do
       if(kpert.gt.0) then
        Call Read_ipar(nut,'kpert',kpert)
        Do i=1,kpert; read(nut,'(a)') line; j=Iadd_line(line); End do
       end if
      end if

      Call br_file

      Call MPI_BCAST(nlsp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kpert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 
      if(nlsp.le.0) Call Stop_mpi(pri,0,'dbsr_conf: nlsp <= 0 ? ')

      if(pri.gt.0) then
       write(pri,'(72(''=''))')
       write(pri,'(a,i3)') 'Number of partial waves, nlsp =',nlsp
       write(pri,'(72(''=''))')
      end if

!---------------------------------------------------------------------- 
! ... do we have the channels to ignore?

      if(myid.eq.0) Call Def_del      

      Call MPI_BCAST(ndel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0.and.ndel.gt.0) &
       Allocate(dlsp(ndel),dkch(ndel),dtar(ndel))

      if(ndel.gt.0) then
       Call MPI_BCAST(dlsp,ndel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(dkch,ndel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(dtar,ndel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      end if

!======================================================================
! ... loop other partial waves:
!======================================================================

      max_ch=0;  max_nc=0;  max_wf=0;  klsp=0

      Do ilsp=1,nlsp

      klsp=klsp+1; if(klsp.gt.nprocs-1) klsp=1; if(myid.ne.klsp) Cycle 
 
      line = aline(ilsp);  ncp = 0; nwp = 0
      read(line,*)  Tpar,Jpar,ipar
      if(LEN_trim(line).gt.12) read(line(19:),*) AFP,BFP,ncp,nwp

       if(pri.gt.0) & 
       write(pri,'(a,5(a,i4))') 'Partial wave:', &
       ' ilsp =',ilsp,'   2J =',jpar,'  parity =',ipar, &
       '  ncp =',ncp, '  nwp =',nwp

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
      if(pri.gt.0) &
      write(pri,'(/a,T40,i8)') &
       'Number of scattering configutarions: ', ncfg_sct-ncfg_targ

      min_ll_ch=0; if(nch.gt.0) min_ll_ch = minval(lch(1:nch))

! ... add perturber configurations if any:

      if(ncp.ne.0)  then
       AF = trim(BFP)//'.c'
       Open(nuc,file=AF,status='OLD')
       Call Read_conf_jj(nuc,0,'add','nocheck')
      end if
      ncfg_pert = ncfg

! ... check imposed orthogonality conditions:

      Call Alloc_orthogonality(0)
      Call Prepare_iort_jj
      Open(nup,file=AF_par,action='READ')
      Call Read_orth_jj(nup)
      write(AF,'(a,i3.3)')'cfg.',ilsp
      Open (nuc,file=AF)
      if(kort.gt.0) Call Read_orth_jj(nuc)

! ... output scattering configuration:    

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

! ... record channel's information:

      write(line,'(i3,a,2x,a,i4,a,2i8)')  ilsp,'.', ' nch =',nch, &
               '    nc =',ncfg-ncfg_targ, ncfg-ncfg_sct
       j = Iadd_line(line)
      Do ich=1,nch
        write(line,'(i4,2x,a5,2i5,3x,i8)') &
              ich,ELF(ipch(ich)),kch(ich),iptar(ich),ipconf(ich)
       j = Iadd_line(line)
      End do

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

      Call SUB_JJ

      if(nch.ne.nch_save) Stop 'nch <> nch_save'

      ncfg_sct=ncfg; lcfg_sct=lcfg; nwf_sct=nwf

      if(pri.gt.0) &
      write(pri,'(/a,T40,i8)') &
       'Number of substitution configutarions: ', ncfg_sct-ncfg_targ

! ... add pertuber physical configurations: 

      if(ncp.gt.0) Call Check_perturber
     
      ncfg_pert=ncfg; lcfg_pert=lcfg; nwf_pert=nwf
 
! ... define orthogonal conditions

      Call Def_orth_cond

! ... analize and record orthogonal conditions for sct. orbitals:

      Call Record_orth
 
      if(pri.gt.0) write(pri,'(/72(''-''))')

      End do ! over partial waves (ilsp)

!----------------------------------------------------------------------
! ... Collect the information:

      if(myid.eq.0) then
       write(nut,'(72(''-''))')
       write(nut,'(a)') 'channels:'
       write(nut,'(72(''+''))')
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      ip = nlsp+kpert+1
      Do ilsp=1,nlsp

       Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       if(myid.ne.0) then
        if(ip.gt.nlines) Cycle
        line = aline(ip);  read(line(1:3),*) klsp; if(klsp.ne.ilsp) Cycle
        Call MPI_SEND(line,124, MPI_CHARACTER,0,i,MPI_COMM_WORLD, ierr)
       end if

       if(myid.eq.0) then
        Call MPI_RECV(line,124, MPI_CHARACTER, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        write(nut,'(a)') trim(line)
        write(nut,'(72(''-''))')
       end if

       j = INDEX(line,'nch =') + 5;  read(line(j:),*) nch

       Do ich = 1,nch
        if(myid.ne.0) then
         line = aline(ip+ich)
         Call MPI_SEND(line, 124, MPI_CHARACTER,0,i,MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(line, 124, MPI_CHARACTER, MPI_ANY_SOURCE, &
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         write(nut,'(a)') trim(line)
        end if
       End do
        
       if(myid.eq.0) write(nut,'(72(''-''))') 

       ip = ip + 1 + nch

      End do

!-----------------------------------------------------------------------
! ... overall information:

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Call MPI_REDUCE(max_ch,m,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(nut,'(a,i8,a)') &
         'max_ch  = ', m, '  -  max.number of channels'
      Call MPI_REDUCE(max_nc,m,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(nut,'(a,i8,a)') &
         'max_nc  = ', m, '  -  max.number of configurations'
      Call MPI_REDUCE(max_wf,m,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(nut,'(a,i8,a)') &
         'max_wf  = ',m,  '  -  max.number of orbitals'
      if(myid.eq.0) write(nut,'(72(''-''))')

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      Call MPI_FINALIZE(ierr)

      End  ! program DBSR_CONF_MPI


