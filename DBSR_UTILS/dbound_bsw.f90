!======================================================================
!     PROGRAM       D B O U N D _ B S W 
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com 
!
!======================================================================
!     converts dbound.nnn  to pairs (c- + bsw-) files
!======================================================================
!     ARGUMENTS:
!
!      klsp  - partial wave
!      state - state index
!      name  - name of states
!
!     INPUT FILES:
!
!      bound_bsw.inp - optional
!      target_jj     - information about channels
!      knot.dat      - B-spline information
!      cfg.nnn       - list of configurations for partial waves ###
!      ubound.nnn    - bound-files for partial wave ###
!      target.bsw    - target B-spline w.f.
!      pert_nnn.bsw  - perturber orbitals
!
!
!     OUTPUT FILES:
!
!      c- and bsw- files for states indicated in bound_bsw.inp
!
!      dbound_bsw.log - running information
!
!----------------------------------------------------------------------
      Use channels_jj; Use target_jj; Use conf_jj; Use orb_jj
      Use DBS_grid; Use DBS_orbitals_pq
      Use DBS_gauss, only: fppqq

      Implicit real(8) (A-H,O-Z)

      Character(80) :: AF,BF
      Character( 3) :: ALSP
      Character(64) :: Lab
      Character( 5), external :: ELi

      Real(8), allocatable :: A(:)
      Integer, allocatable :: ipt(:)

      Integer :: inp = 5; Character(40) :: AF_inp = 'dbound_bsw.inp'
      Integer :: pri = 6; Character(40) :: AF_log = 'dbound_bsw.log'
      Integer :: nut = 1; Character(40) :: AF_tar = 'target_jj'
      Integer :: nuw = 2; Character(40) :: AF_bsw = 'target.bsw'
      Integer :: nuc = 3; Character(40) :: AF_cfg = 'cfg.nnn'
      Integer :: nub = 4; Character(40) :: AF_bnd = 'dbound.'

      Character(80) :: name =' ',  mode = ' '
      Integer :: state = 0
      Integer :: klsp = 0
      Integer :: kk=ICHAR('k')  
      Real(8) :: eps_c = 1.d-8

      Call Read_name(name)

      if(name == '?') then
       write(*,*) 
       write(*,*) 'convert  dbound.nnn  to pairs  (c- + bsw-)  files'
       write(*,*) 'for selected states indicated in dbound_bsw.inp'
       write(*,*) 
       write(*,*) 'another option is to use command-line arguments:' 
       write(*,*) 
       write(*,*) 'dbound_bsw  klsp=... state=... name=...'
       write(*,*) 
       write(*,*) 'klsp  -  partial wave index'
       write(*,*) 'state -  state index'
       write(*,*) 'name  -  name for given state: name.c and name.bsw'
       write(*,*) 
       Stop ' '
      end if

      Call Read_rarg('eps_c',eps_c)

!----------------------------------------------------------------------
! ... input data file:          

      Call Read_aarg('mode',mode);  imode=len_trim(mode)
      jnp = Icheck_file(AF_inp)
      if(jnp.gt.0) open(inp,file=AF_inp)

! ... read target and channel information:
   
      Call Check_file(AF_tar)
      Open(nut,file=AF_tar)
      Call Read_target_jj(nut)
      Call Read_channels_jj(nut)

! ... set up B-splines:
 
      Call read_knot_dat
      Call alloc_DBS_gauss

! ... target w.f.:

      Open(nuw, file=AF_bsw, STATUS='OLD', form='UNFORMATTED')
      Call Read_dbsw(nuw,0,0)
      Close(nuw)
      nwt = nbf

      Open(pri,file=AF_log,position='APPEND')

!----------------------------------------------------------------------
    1 Continue 

      if(jnp.gt.0) then
       read(inp,'(a)',err=2,end=2) BF
       if(BF(1:1).eq.'*') go to 2
       read(BF,*,err=2,end=2) klsp,state,name
      else 
       klsp=0;    Call Read_iarg('klsp',klsp)
       state=0;   Call Read_iarg('state',state)
       name=' ';  Call Read_aarg('name',name)
       kk=ICHAR('k'); Call Read_iarg('n',kk)
       Call Read_aarg('AF_bnd',AF_bound)            
      end if
      if(klsp.eq.0) go to 2
      if(state.eq.0) go to 2
      if(len_trim(name).eq.0.and.imode.eq.0) go to 2
      if(imode.gt.0) write(name,'(a,a,i3.3,a,i3.3)') trim(mode),'_',klsp,'_',state

      if(name.eq.'?'.or.state.eq.0.or.klsp.eq.0) &
       Stop 'Call as:  dbound_bsw klsp=... state=... name=...'

      write(pri,'(/72a/)') ('-',i=1,72)
    
      write(pri,'(a,i5)') 'klsp  = ',klsp
      write(pri,'(a,i5)') 'state = ',state
      write(pri,'(a,a )') 'name  = ',trim(name)
      write(pri,*)
      write(ALSP,'(i3.3)') klsp

! ... perturber radial functions:                                             ???

      nbf=nwt
      if(nwp(klsp).gt.0) then
       BF = trim(BFP(klsp))//'.bsw' 
       Open(nuw, file=BF, STATUS='OLD', form='UNFORMATTED')
       Call Read_dbsw(nuw,2,0)
       Close(nuw)
      end if

! ... allocate space for outer orbitals:

      nwb = nbf
      Call alloc_DBS_orbitals_pq(nwb+nch(klsp))
      nbf = nwb + nch(klsp)

      write(pri,'(a,i5)') 'nwt   = ',nwt
      write(pri,'(a,i5)') 'nwb   = ',nwb
      write(pri,'(a,i5)') 'nch   = ',nch(klsp)
      write(pri,'(a,i5)') 'ncp   = ',ncp(klsp)
      write(pri,'(a,i5)') 'npert = ',npert(klsp)

! ... read configurations:

      Call alloc_cfg(-1)
      i=Index(AF_cfg,'.',BACK=.TRUE.); AF_cfg=AF_cfg(1:i)//ALSP
      Open(nuc,file=AF_cfg,status='OLD')
      Call Read_core_jj(nuc)
      Call Read_conf_jj(nuc,0,'add','nocheck')
      Close(nuc)
      write(pri,'(a,i5)') 'ncfg  = ',ncfg

! ... read dbound.###:

      write(AF,'(a,i3.3)')  trim(AF_bnd),klsp
      
      Call Check_file(AF)

      if(AF(1:1).eq.'d') then

       Open(nub,file=AF,form='UNFORMATTED')
       rewind(nub)
       read(nub) nhm,nchan,kpert,ns_b,jj,ipar_b,nbound
 
       if(ns_b.ne.ns)   Stop ' number of splines, ns  --> ?'
       if(nch(klsp).ne.nchan) Stop ' number of channels, nch --> ?'
       if(npert(klsp).ne.kpert) Stop ' pertuber, ncp --> ?'
       if(nhm.ne.ms*nchan+kpert) Stop ' full expansion, nhm --> ?'
 
       if(allocated(A)) Deallocate(A); Allocate(A(nhm))
 
       ! ... find given state:
 
       if(state.lt.1.or.state.gt.nbound)  Stop 'state is out of range'
 
       Do i = 1,state
        read(nub) ii,Lab; read(nub) ET;  read(nub) A
       End do 
         
       write(pri,'(/a,i5,F16.8,2x,a/)') 'State: ', state, ET, TRIM(Lab)

      elseif(AF(1:1).eq.'p') then


      else

       Stop 'unknown format for bound solution files'

      end if 

! ... transfer the solution to the p-functions:

      ic1=0; ic2=0
      Do ich = 1,nchan

       i=nwb+ich;  mbs(i)=ns; ebs(i)=ELC(klsp,ich)

       Call EL_NLJK(ebs(i),nbs(i),kbs(i),l,j,ibs(i))
       j=(ich-1)*ms + 1
       Call Put_pv(i,A(j),ns)

       ! ... renormalization:
 
       s = QUADR_pq(i,i,0);  ss = sqrt(s); if(ss.lt.1.d-8) ss=0.d0
       if(ss.ne.0.d0) pq(:,:,i) = pq(:,:,i)/ss
       ic1=ic2+1; ic2=ipconf(klsp,ich)
       WC(ic1:ic2) = WC(ic1:ic2) * ss
       S = SUM(WC(ic1:ic2)*WC(ic1:ic2))
       write(pri,'(a5,2F15.8)') ebs(i),ss,S

       if(S.eq.0.d0) mbs(i) = 0

      End do

      if(ic2+ncp(klsp).ne.ncfg) Stop ' Problem with ic1,ic2 ...'

! ... weights of perturbers:

      if(kpert.gt.0) then
       ishft=nchan*ms
       jshft =ipconf(klsp,nchan)
       j = ipert(klsp)
       Do i=1,npert(klsp); 
        ss = A(ishft+i)
        ic1=jshft+ippert(j+i-1)+1
        if(i.eq.1) ic1=jshft+1
        ic2=jshft+ippert(j+i)
        WC(ic1:ic2) = WC(ic1:ic2) * ss
        write(pri,'(a,i3,a,2i5,a,f12.8)') &
         'perturber',i,'   ic1,ic2 =',ic1,ic2,'  S =',ss
       End do         
      end if

! ... replace the spectroscopic notation:

      Do ich=1,nchan
       Call EL_NLJK(ELC(klsp,ich),n,kappa,l,j,k)
       i = Ifind_jjorb(n,kappa,k,1)
       Call EL_NLJK(ELC(klsp,ich),n,kappa,l,j,k)
       j = Ifind_bsorb(n,kappa,k,1)
       NEF(i) = kk
       ELF(i) = ELi(NEF(i),KEF(i),IEF(i))
       ebs(j) = ELF(i)
      End do

! ... output c - file ...

      if(allocated(ipt)) Deallocate(ipt); Allocate(ipt(ncfg))
      ipt = 0;  Call SORTA(ncfg,WC,ipt)

      ii = LEN_TRIM(name);  AF = name(1:ii)//'.c';  Open(nuc,file=AF)

      write(nuc,'(a15,f16.8)') 'Core subshells:',ET
      write(nuc,'(a)') core(1:ncore*5)
      write(nuc,'(a)') 'Peel subshells:'
      write(nuc,'(20a5)') (ELF(i),i=ncore+1,nwf)
      write(nuc,'(a)') 'CSF(s):'
      Do jc = 1,ncfg; ic=ipt(jc)
       if(abs(WC(ic)).gt.eps_c)  Call Print_conf_jj(nuc,ic,WC(ic))
      End do
      write(nuc,'(a)') '*'
      Close(nuc)
               
!  ... output bsw-file:
       
      AF = name(1:ii)//'.bsw'; Open(nuw,file=AF,form='UNFORMATTED')
       
      write(nuw) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i=1,nbf
       if(mbs(i).eq.0) Cycle
       write(nuw) ebs(i),mbs(i)
       write(nuw) pq(1:mbs(i),1,i)
       write(nuw) pq(1:mbs(i),2,i)
      End do

      nbf=nwt
      write(pri,'(/72(''-'')/)')

      if(jnp.gt.0) go to 1
    2 Continue

      End  ! program DBOUND_BSW
          
         