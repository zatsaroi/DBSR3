!======================================================================
!     PROGRAM       D S R _ P R E P                          version 3
!
!               C O P Y R I G H T -- 2019
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!
!  This proram analyzes the orthogonality conditions for one-electron 
!  orbitals and generates new c- and bsw-files for target states with 
!  consistent orbital set-indexes and sorted according to their energy.
!
!  In order to define the same orbitals, the following criteria are
!  used:
!              | <p1|p2> - 1 |  <  eps_ovl
!              | <p1|r|p1> - <p2|r|p2> | <  eps_ovl
!              | <p1|1/r|p1> - <p2|1/r|p2> | <  eps_ovl
!  Similar conditions are used to define orthogonal orbitals.
!======================================================================
!
! INPUT FILES:
!
!     dbsr_par    -  parameters of calculation (optional)
!     target_jj   -  list of target states and partial waves
!     knot.dat    -  B-spline parameters
!     name.c and name.bsw - files for target states and perturbers (if any)
!     target_sub.bsw - optional file with additional orbitals used
!                      to treat the orthogonality constraints
!
! OUTPUT FILES:
!
!     target_jj     -  modified list of target states
!     target.bsw    -  all target bsw-functions
!     targ_nnn.c    -  c-file for target nnn
!     pert_nnn.c    -  c-file for perturber nnn, if any
!     pert_nnn.bsw  -  w-file for perturber nnn, if any
!     target_orb    -  list of physical orbitals in target states
!     dbsr_prep.log -  running information
!
!======================================================================
      Use dbsr_prep

      Implicit real(8) (A-H,O-Z)

      Integer, allocatable :: klsp(:)
      Character(ma), allocatable :: AFK(:)
      Real(8) :: S, t1, t2

      Call CPU_time(t1)

      Call DBSR_prep_inf

! ... check and open files:

      Open(pri,file=AF_log)
      Open(nup,file=AF_par)
      Open(nuo,file=AF_orb)

      if(Icheck_file(AF_tar).eq.1) then
       Open(nut,file=AF_tar)
      else
       Open(nut,file=AF_tar)
       Call Write_target_jj_example(nut)       
       write(*  ,*) 'Prepare file "target" - see created example'
       write(pri,*) 'Prepare file "target" - see created example'
       Stop
      end if

! ... sets up grid points and initializes the values of the B-spline:

      Call read_knot_dat
      Call alloc_DBS_gauss

! ... read parameters if any:

      Call Read_arg

!----------------------------------------------------------------------
! ... proceed target states:

      Call Read_ipar(nut,'nz',nz)
      Call Read_ipar(nut,'nelc',nelc)
      Call Read_ipar(nut,'ntarg',m)
      Call Alloc_target_jj(m)
      read(nut,*)
      Do it=1,ntarg
       read(nut,*) AFT(it); Call Check_BSR_name(AFT(it))
      End do

! ... define core:

      AFC = trim(AFT(1))//'.c'
      Call Check_file(AFC)
      Open(nuc,file=AFC,status='OLD')
      Call Read_core_jj(nuc)
      write(pri,'(80(''-''))')
      write(pri,'(a,i5/)') 'CORE:  ncore =', ncore
      if(ncore.gt.0) then
       write(pri,'(10a8)') (ELF(j),j=1,nwf)
       AFW = trim(AFT(1))//'.bsw'
       Open(nuw,file=AFW,form='UNFORMATTED')
       Call read_dbsw(nuw,0,0)

       Do i=1,ncore
        Do j=1,ncore
        if(kbs(i).ne.kbs(j)) Cycle
        S = QUADR_pq(i,j,0)
        if(abs(S).lt.eps_core) Cycle 
        Call Iadd_obs(i,j,S) 
       End do
       OBS1(i) = QUADR_pq(i,i,1)
       OBS2(i) = QUADR_pq(i,i,2)
      End do

      end if

      write(pri,'(80(''-''))')

! ... define the target energies and check cores:

      Do it=1,ntarg
       AFC = trim(AFT(it))//'.c'
       Call Check_file(AFC);  Open(nuc,file=AFC)
       rewind(nuc)
       read(nuc,'(15x,F16.8)') Etarg(it)
       read(nuc,'(a)') CLOSED
       if(CLOSED.ne.core)  then
        write(pri,'(/a,a,a)') 'target ',AFT(it),' has different core'
        write(*,'(/a,a,a)') 'target ',AFT(it),' has different core'
        Stop 
       end if
       Close(nuc)
      End do

! ... sorting of target states according to energy:

      Do i = 1,ntarg-1; Do j = i+1,ntarg
       if(Etarg(j).ge.Etarg(i)) Cycle
       E=Etarg(i); Etarg(i)=Etarg(j); Etarg(j)=E
       AF=AFT(i); AFT(i)=AFT(j); AFT(j)=AF
       if(Etarg(i).eq.Etarg(j)) then
        write(*,*) 'Check target states with the same energies: ', &
         trim(AFT(i)),'   ',trim(AFT(j))
        write(pri,*) 'Check target states with the same energies: ', &
         trim(AFT(i)),'   ',trim(AFT(j))
       end if
      End do; End do

      nct = 0; nwt = ncore;  nbf = ncore

!-----------------------------------------------------------------------
! ... check if substitute orbitals are provided:

      if(Icheck_file(AF_sub).ne.0) then

       write(pri,'(a/)') 'given substitution orbitals:'
       ncfg = 0; lcfg = 0;  nwf = ncore
       AFC = 'no';   AFW = AF_sub;    BFC = 'no';   BFW = 'no'
       open(nuw,file=AF_sub,form='UNFORMATTED')
       Call Read_bsw_orb_jj(nuw)  ! read orbitals into orb_LS module

       kshift=0; ipef=0;   CALL SUB_check_orb

       write(pri,'(80(''-''))')

      end if

      nwt = nbf;  ipbs = 0; ipbs(1:nbf) = 1  ! sub. orb. pointer

!-----------------------------------------------------------------------
! ... check orbitals in target states:

      Do it=1,ntarg

       write(pri,'(a,i8,T24,a12/)') 'target ',it,AFT(it)

       AFC=trim(AFT(it))//'.c'
       Call Check_file(AFC)
       Open(nuc,file=AFC)
       nwf=ncore; ncfg=0; lcfg=0;  Call Read_conf_jj(nuc,0,'detect','check')
                                   Call Test_ac
       AFW=trim(AFT(it))//'.bsw'

       if(ntarg.lt.1000) then
         write(BFC,'(a,i3.3,a)') 'targ_',it
       elseif(ntarg.lt.10000) then
         write(BFC,'(a,i4.4,a)') 'targ_',it
       elseif(ntarg.lt.100000) then
         write(BFC,'(a,i5.5,a)') 'targ_',it
       else
         Stop 'Stop in bsr_prep:  ntarg > 100000 '
       end if
       BFT(it) = BFC;  BFC = trim(BFT(it))//'.c';  BFW = 'no'
 
       kshift=0; ipef=0;    CALL SUB_check_orb
 
      ! ... write new c-file:

       if(allocated(ipt)) deallocate(ipt); Allocate(ipt(ncfg))
       Call SORTA(ncfg,WC,ipt)

       Open(muc,file=BFC)
       write(muc,'(a,f16.8)') 'Core subshells:',Etarg(it)  
       write(muc,'(a)') core(1:ncore*5)
       write(muc,'(a)') 'Peel subshells:'
       write(muc,'(20a5)') ELF(ncore+1:nwf)
       write(muc,'(a)') 'CSF(s):'
 
       ncfg1 = 0
       Do jc=1,ncfg; ic=ipt(jc); if(abs(WC(ic)).lt.c_targ) Cycle
        Call Print_conf_jj(muc,ic,WC(ic))
        ncfg1 = ncfg1 + 1
       End do
       write(muc,'(a)') '*'
 
      ! ... define substitution orbitals:

       write(nuo,'(a,3x,i3.3,6x,a)') 'target',it,AFT(it)

       Call SUB_phys_orb

       nct=nct+ncfg1; nctarg(it)=ncfg1; nwtarg(it)=nbf-nwt; nwt=nbf
       write(pri,'(/a,i6,T24,a )') 'nctarg = ',nctarg(it),'number of configurations'
       write(pri,'( a,i6,T24,a/)') 'nwtarg = ',nwtarg(it),'number of new orbitals'
       write(pri,'(80(''-''))')

       Call Jdef_JP(nuc,jtarg(it),ptarg(it))

      End do

!-----------------------------------------------------------------------
! ... create the target.bsw file:

      Open(muw,file=AF_wfn,form='UNFORMATTED')
      write(muw) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i = 1,nwt
       write(muw) ebs(i),mbs(i)
       write(muw) (pq(j,1,i),j=1,mbs(i))
       write(muw) (pq(j,2,i),j=1,mbs(i))
      End do
      Close(muw)

! ... final list of target one-electron orbitals:

      write(pri,'(a,i5/)') 'target orbitals (beside core) = ',nwt-ncore
      write(pri,'(10a8)') (EBS(j),j=ncore+1,nwt)
      write(pri,'(80(''-''))')

! ... substitution orbitals information:

      nwf_sub = 0
      Do i=1,nbf; if(ipbs(i).gt.0) nwf_sub=nwf_sub+1; End do
      write(pri,'(/a,i5/)') &
       'substitution orbitals (beside core) = ',nwf_sub-ncore
      k = 1
      Do i=ncore+1,nbf
       if(ipbs(i).ne.1) Cycle
       write(AS(k:),'(3x,a5)') ebs(i)
       k=k+8
       if(k.lt.80) Cycle
       write(pri,'(a)') AS
       k=1
      End do
      if(k.gt.1) write(pri,'(a)') AS

      k = 0
      Do i=1,nbf; if(ipbs(i).ne.1) Cycle
       Do j=1,i; if(ipbs(j).ne.1) Cycle
        if(kbs(i).ne.kbs(j)) Cycle
        S = abs(OBS(i,j)); if(i.eq.j) S=S-1.d0
        if(S.lt.eps_ovl) Cycle
        write(pri,'(2a8,f20.8)') ebs(i),ebs(j),S
        k=1
       End do
      End do
      if(k.eq.0) write(pri,'(/a/)') &
         'all substitution orbitals are orthogonal'
      write(pri,'(80(''-''))')
 
!=====================================================================
! ... define number of partial waves:

      nlsp=0;  Call Read_ipar(nut,'nlsp',nlsp)

      if(nlsp.gt.0) then

       Allocate(jpar(nlsp),ipar(nlsp), &
                AFP(nlsp),BFP(nlsp),ncp(nlsp),nwp(nlsp))
       ncp = 0
       read(nut,*)
       Do i = 1,nlsp
        read(nut,'(a)') AS;  ii = LEN_TRIM(AS)
        read(AS,*) AF,jpar(i),ipar(i)
        AFP(i) = 'no'; BFP(i) = 'no'
        if(ii.gt.13) then                            ! ???
         read(AS(14:),*) AFP(i)
         Call Check_BSR_name(AFP(i))
        end if
       End do

      else               ! define range of partial waves

      if(JJ_min.lt.0) JJ_min = mod(nelc+1,2)  
      if(JJ_max.lt.0) JJ_max = mod(nelc+1,2) + 50 

       nlsp = JJ_max-JJ_min + 2
       Allocate(jpar(nlsp),ipar(nlsp), &
                AFP(nlsp), BFP(nlsp),ncp(nlsp),nwp(nlsp))
       i = 0
       Do JJ = JJ_min,JJ_max,2
        Do ip = -1,1,2
         i = i + 1
         jpar(i)=JJ; ipar(i)=ip; AFP(i) = 'no'; ncp(i)=0
        End do
       End do

      end if    ! over nlsp > 0

! ... define number of additional perturbers if any:

      kpert=0; Call Read_ipar(nut,'kpert',kpert)
      if(kpert.gt.0) then
       Allocate(klsp(kpert),AFK(kpert))
       read(nut,*)
       Do i=1,kpert
        read(nut,*) klsp(i),AFK(i)
        Call Check_BSR_name(AFK(i))
       End do
      end if


!=====================================================================
! ... proceed perturbers if any:

      Do ilsp=1,nlsp
       
       BFC = 'no'; BFW = 'no'
       nbf = nwt; ncfg = 0; lcfg = 0; nwf = ncore  
       if(mbf.gt.nwt) ipbs(nwt+1:mbf) = 0

       npert = 0;  Call Allocate_pert(0)

! ... usual perturber (as list of separate configurations) :

       if(LEN_TRIM(AFP(ilsp)).gt.0.and.AFP(ilsp).ne.'no') then

        AFC = trim(AFP(ilsp))//'.c'  
        write(pri,'(a,i6,T24,a/)') 'perturber', ilsp, trim(AFC)
        Call Check_file(AFC)
        Open(nuc,file=AFC)
        Call Read_conf_jj(nuc,0,'detect','check')
        WC(1:ncfg) = 1.d0
        ncp(ilsp) = ncfg

        Call Allocate_pert(ncfg+ipert); npert=ncfg
 
        Do i=1,npert; ippert(i)=i; End do
        
        AFW = trim(AFP(ilsp))//'.bsw';

        Call Check_file(AFW)

        kshift=0; ipef=0;  CALL SUB_check_orb

       end if

! ... additional perturbers:

       Do iip=1,kpert

        if(klsp(iip).ne.ilsp) Cycle
        AFC = trim(AFK(iip))//'.c'
        Call Check_file(AFC)
        write(pri,'(/a,i5,5x,a/)') 'perturber', ilsp, trim(AFK(iip))
        Open(nuc,file=AFC)

        kshift = maxval(KEF(1:nwf))

        Call Read_conf_jj(nuc,kshift,'add','check')

        if(npert+1.gt.mpert) Call Allocate_pert(mpert+ipert)
        npert = npert + 1        
        ippert(npert)=ncfg
        ncp(ilsp) = ncfg
        AFW = trim(AFK(iip))//'.bsw'
        Call Check_file(AFW)

        CALL SUB_check_orb
   
       End do 

! ... collective perturber configurations:

       if(npert.le.0) Cycle

       write(BFP(ilsp),'(a,i3.3)') 'pert_',ilsp
       BFC = trim(BFP(ilsp))//'.c'
       Open(muc,file=BFC); rewind(muc)
       write(muc,'(a)') 'Core subshells:'
       write(muc,'(a)') core(1:ncore*5)
       write(muc,'(a)') 'Peel subshells:'
       write(muc,'(20a5)') ELF(ncore+1:nwf)
       write(muc,'(a)') 'CSF(s):'

       Do ic=1,ncfg; Call Print_conf_jj(muc,ic,WC(ic));  End do
       write(muc,'(a)') '*'

       write(muc,'(a,i5)') 'npert =',npert
       write(muc,'(20i5)') (ippert(i),i=1,npert)

       write(nuo,'(a,3x,i3.3)') 'pertub',ilsp

           Call SUB_phys_pert

       ncp(ilsp)=ncfg; nwp(ilsp)=nbf-nwt

       if(nwp(ilsp).gt.0) then
        BFW = trim(BFP(ilsp))//'.bsw'
        Open(muw,file=BFW,form='UNFORMATTED')
        write(muw) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
        Do i = nwt+1,nbf
         write(muw) ebs(i),mbs(i)
         write(muw) (pq(j,1,i),j=1,mbs(i))
         write(muw) (pq(j,2,i),j=1,mbs(i))
        End do
        Close(muw)
       end if
 
       write(pri,'(/a,i8,T24,a)') 'ncp  = ',ncp(ilsp),'number of configurations'
       write(pri,'( a,i8,T24,a)') 'nwp  = ',nwp(ilsp),'number of new orbitals'
       write(pri,'(80(''-''))')

      End do  ! over partial waves

!---------------------------------------------------------------------
! ... update the target file:

      rewind(nut)
      read(nut,'(a80)') title

      Call write_target_jj(nut)
      Call write_channels_jj(nut,0)
      if(kpert.gt.0) then
       write(nut,'(a,i4,5x,a)') &
        'kpert = ',kpert,' !   number of additional perturbers' 
      write(nut,'(80(''-''))')
      Do i = 1,kpert
       write(nut,'(i3,5x,a)') klsp(i),trim(AFK(i))
      End do 
      write(nut,'(80(''-''))')
      end if

! ... information for target_orb file:      

      write(nuo,'(/a,a)') &
 'This file contains the list of physical oritals along with:', &
 'occupation number, substiution orbital, corresponding overlap'
      Close(nuo)

      Call CPU_time(t2)

      write(pri,'(/a,T20,f10.2,a)') 'total time:  ', (t2-t1)/60, ' min'

      End    ! program DBSR_PREP
