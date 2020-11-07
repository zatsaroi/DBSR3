!======================================================================
!     PROGRAM       D B S R _ M A T                          version 3
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny 
!     email:        oleg_zoi@yahoo.com
!
!======================================================================
!     Set up Hamiltonian matrixes in B-spline representation
!======================================================================
!
!     INPUT ARGUMENTS:
!
!     klsp1, klsp2   -  range of partial wave under consideration
!
!     INPUT FILES: 
!
!     target_jj      -  description of target states and channels
!     knot.dat       -  B-spline grid
!     dbsr_par       -  input parameters if any
!     cfg.nnn        -  configuration list for given partial wave 'nnn'
!     int_bnk.nnn    -  angular coefficient data bank
!     target.bsw     -  target w.f.'s in B-spline basis
!     target_orb     -  description of physical orbitals
!     pert_nnn.bsw   -  perturber's w.f., if any
!
!     OUTPUT FILES:
!
!     bsr_mat.nnn    -  resulting interaction matrix
!     mat_log.nnn    -  running information
!
!=====================================================================
      Use dbsr_mat         
      
      Implicit none
      Real(8) :: t1,t2,t3
       
      Call CPU_time(t1)

! ... file for general log-output:

      Open(prj,file=AF_prj)       

! ... read arguments:

      Open(nup,file=AF_par,status='OLD')
      Call Read_arg(nup)

! ... target information:

      Call Check_file(AF_tar);  Open(nut,file=AF_tar)
      Call Read_target_jj(nut)
      
! ... prepare B-splines and other relevant arrays:

      Call read_knot_dat
      Call alloc_DBS_gauss
      Call Def_Vnucl

      Call alloc_Rk_integrals(ns,ks,0,mk,ntype_R) 
      if(mbreit.gt.0)  Call alloc_Sk_integral(ns,ks)

! ... find number of closed shells and core energy:

      Call Def_core

! ... loop over partial waves:

      Do klsp = klsp1,klsp2
 
       write(*,'(/a,i3/)') 'DBSR_MAT:  klsp =', klsp

       Call CPU_time(t2);  Call SUB1;  Call CPU_time(t3)

       write(pri,'(/a,T40,f10.2,a)') 'total time:', (t3-t2)/60, ' min'
       write(*  ,'(/a,T40,f10.2,a)') 'total time:', (t3-t2)/60, ' min'

      End do  ! over klsp

      Call CPU_time(t2)
      write(prj,'(/a,T40,f10.2,a)') 'total time:  ', (t2-t1)/60, ' min'

      End ! program dbsr_mat


!======================================================================
      Subroutine Def_core
!======================================================================
!     calculate and broadcast the core energy
!----------------------------------------------------------------------
      Use dbsr_mat         
      
      Implicit none
      Integer :: i,j
      Integer, external :: Ifind_bsorb
      Real(8), external :: Ecore_dbs

! ... find nclosd and core energy:

      Ecore = 0

      i=LEN_TRIM(AF_cfg); write(AF_cfg(i-2:i),'(i3.3)') klsp1
      Call Check_file(AF_cfg); Open(nuc,file=AF_cfg)
      Call Read_core_jj(nuc)
      
      if(ncore.eq.0) Return

      Do i=1,nwf; j=Ifind_bsorb(NEF(i),KEF(i),IEF(i),2); End do
      mbs = 0
      Open(nuw, file=AF_bsw, STATUS='OLD', form='UNFORMATTED')
      Call Read_dbsw(nuw,0,0)
      Close(nuw)

      Ecore = Ecore_dbs(ncore,mbreit,kbs) 

      write(prj,'(/a,i10,T20,a)')  'ntarg  =',ntarg,'- number of target states'
      write(prj,'(/a,i10,T20,a)')  'ncore  =',ncore,'- number of core subshells'
      write(prj,'(/a,F16.8,5x,a)') 'Ecore  =',Ecore,'- calculated core energy'
      write(prj,'( a,F16.8,5x,a)') 'Ec     =',EC,   '- given core energy if any'
      if(EC.ne.0.d0) Ecore = EC
      EC = Ecore

      End Subroutine Def_core
