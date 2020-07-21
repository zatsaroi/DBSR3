!======================================================================
      Subroutine Read_data
!======================================================================
!     read input parameters for given case
!----------------------------------------------------------------------
      Use dbsr_ci;             Use DBS_grid 
      Use conf_jj;             Use DBS_orbitals_pq   
      Use orb_jj;              Use DBS_nuclear       
      Use c_data;              Use DBS_gauss         
                               
      Implicit none
      Integer :: i,j 
      Real(8) :: S
      Integer, external :: Icheck_file
      Real(8), external :: QUADR_pq, OBS

  
!----------------------------------------------------------------------
! ... input data file if any:

      AF_inp = trim(name)//'.inp_ci'  
      inp=Icheck_file(AF_inp)
      if(inp.gt.0) open(inp,file=AF_inp)

! ... log-file:

      AF_log = trim(name)//'.ci'
      Open(pri,file=AF_log)

! ... c-file:

      AF_c = trim(name)//'.c'
      Call Check_file(AF_c)
      Open(nuc,file=AF_c,status='OLD')

! ... bsw-file:

      AF_bsw = trim(name)//'.bsw'
      Call Check_file(AF_bsw)
      Open(nuw,file=AF_bsw,status='OLD',form='UNFORMATTED')

! ... bnk-file (angular coefficients):
 
      AF_bnk = trim(name)//'.bnk'
      i = Icheck_file(AF_bnk)
      if(i.eq.0)  then; AF_bnk=BF_bnk; i = Icheck_file(BF_bnk); end if
      if(i.eq.0) Stop 'dbsr_ci: cannot find bnk-file'
      Open(nub,file=AF_bnk,status='OLD',form='UNFORMATTED')

! ... j-file:

      AF_j = trim(name)//'.j'
      Open(nuj,file=AF_j)

!----------------------------------------------------------------------
! ... HEADER:

      write(pri,'(a,5x,a/)')   'DBSR_CI:',trim(name)

!----------------------------------------------------------------------
! ... read obitals and configurations along with expansion coeff's: 

      write(pri,'(/a/)')   'c-file data: '
      Call Read_core_jj(nuc)
      write(pri,'(a,i8)')  'ncore  = ',ncore;  nclosed=ncore
      Call Read_conf_jj
      write(pri,'(a,i8)')  'ncfg   = ',ncfg
      write(pri,'(a,i8)')  'nwf    = ',nwf
      write(pri,'(10a8)')  (ELF(i),i=1,nwf)

!----------------------------------------------------------------------
! ... define the knot sequence:

      Call def_grid(knot,name,z,awt)

! ... initialize other B-spline arrays:

      Call alloc_DBS_gauss
      Call def_Vnucl

!      Call alloc_RK_integrals(ns,ks,0,mpol,4)   !???  
      Call alloc_Rk_integrals(ns,ks,0,mpol,4)
      if(mbreit.gt.0) &
      Call alloc_Sk_integrals(ns,ks,0,mpol,2)       

      write(pri,'(/a,f5.2)') 'Z        =  ', Z
      write(pri,'(a,a)')     'nuclear  =  ', nuclear
      write(pri,'(a,i5)')    'mbreit   =  ', mbreit

      Call Conv_au (Z,AWT,au_cm,au_eV,pri)

      Call alloc_DBS_orbitals_pq(nwf,ns)

      Do i=1,nwf
       ebs(i) = ELF(i)
       nbs(i) = nef(i)
       kbs(i) = kef(i)
       lbs(i) = lef(i)
       jbs(i) = jef(i)
       ibs(i) = IEF(i)
       mbs(i) = 0
      End do
      nbf = nwf

      Call Read_pqbs(nuw)

! ... check the correspondence between c- and bsw-files: 

      j = 0
      Do i = 1,nwf
       if(mbs(i).ne.0) Cycle
       write(pri,'(a,a)') 'Absent expansion for w.f. ',ebs(i)
       j = j + 1
      End do
      if(j.gt.0) Stop 'no correspondence between c- and w- files'
    

! ... the < p | p > values 

      Do i=1,nbf;  Do j=1,i
        S=QUADR_pq(i,j,0)
        if(abs(S).lt.Eps_ovl) Cycle
        Call Iadd_obs(i,j,S)
       End do
      End do

! ... one-electron overlaps:

      if(debug.gt.0) then 
       write(pri,'(/a,f10.8/)') &
        'Non-trivial one-elctron overlaps: > eps_ovl = ',eps_ovl
       Do i=1,nbf
       Do j=i,nbf
        if(kbs(i).ne.kbs(j)) Cycle
        S = OBS(i,j)  
        if(i.ne.j.and.abs(S).gt.eps_ovl) &
         write(pri,'(''<'',a5,''|'',a5,''>='',4E15.5)') &
         EBS(i),EBS(j),OBS(i,j)
        if(i.eq.j.and.abs(S-1.d0).gt.eps_ovl) &
         write(pri,'(''<'',a5,''|'',a5,''>='',3E15.5)') &
         EBS(i),EBS(j),S
       End do; End do
      end if

      End Subroutine Read_data
    

