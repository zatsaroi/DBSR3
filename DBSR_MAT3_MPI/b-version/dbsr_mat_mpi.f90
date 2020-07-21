!======================================================================
!     PROGRAM       D B S R _ M A T _ M P I   
!
!                   C O P Y R I G H T -- 2019
!
!     Written by:   Oleg Zatsarinny 
!     email:        oleg_zoi@yahoo.com
!======================================================================
!     Generation of Hamiltonian matrixes in B-spline representation
!======================================================================
!
!   MAIN INPUT FILES: 
!
!     knot.dat       -  B-spline grid
!     dbsr_par       -  input parameters 
!     target_jj      -  description of target states and channels
!     target.bsw     -  target w.f.'s in B-spline basis
!
!     cfg.nnn        -  configuration list for given partial wave 'nnn'
!     int_bnk.nnn    -  angular coefficients data bank
!     pert_nnn.bsw   -  perturb w.f., if any
!
!   MAIN OUTPUT FILES:
!
!     dbsr_mat.nnn   -  overlap and Hamiltonian matrixes 
!     dbsr_mat.log   -  general running information
!     mat_log.nnn    -  running information for given partial wave 
!
!=====================================================================
      Use MPI
      Use dbsr_mat

      Implicit none
      Real(8) :: t1,t2,t3

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      t1 = MPI_WTIME()

! ... file for general log output:

      if(myid.eq.0) then; Open(prj,file=AF_prj); else; prj=0; end if

! ... target:
 
      if(myid.eq.0) then
       Call Check_file(AF_tar);  Open(nut,file=AF_tar);  Call Read_target_jj(nut)
      end if
      Call br_target_jj

! ... define arguments:

      if(myid.eq.0) then
       Call Check_file(AF_par); Open(nup,file=AF_par); Call Read_arg(nup)
      end if
      Call br_arg
      
! ... prepare B-spline parameters:

      if(myid.eq.0) Call read_knot_dat;   Call br_grid 

      Call alloc_DBS_gauss
      Call Def_Vnucl

      Call alloc_Rk_integrals (ns,ks,0,mk,ntype_R)  

!      Call alloc_Rk_integral(ns,ks)

      if(mbreit.gt.0) Call alloc_Sk_integral(ns,ks)

! ... find nclosd and core energy:

      Call Def_core

! ... loop over partial waves:

      Do klsp = klsp1,klsp2
 
       if(myid.eq.0) write(*,'(/a,i3/)') 'DBSR_MAT:  klsp =', klsp

       t2=MPI_WTIME();    Call SUB1;    t3=MPI_WTIME()

       if(pri.gt.0) &
       write(pri,'(/a,T30,f9.2,a)') 'total time:', (t3-t2)/60, ' min'

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End do  ! over klsp

      t2 = MPI_WTIME()
      if(prj.gt.0) &
      write(prj,'(/a,T30,f9.2,a)') 'bsr_mat:  ',(t2-t1)/60,'  min'

      Call MPI_FINALIZE(ierr)

      End  ! program dbsr_mat_mpi



!======================================================================
      Subroutine Def_core
!======================================================================
!     calculate and broadcast the core energy
!----------------------------------------------------------------------
      Use MPI
      Use dbsr_mat

      Implicit none
      Real(8), external :: Ecore_dbs
      Integer :: i,j
      Integer, external :: Ifind_bsorb

      if(myid.eq.0) then

      write(ALSP,'(i3.3)') klsp1  
      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Call Check_file(AF_cfg); Open(nuc,file=AF_cfg)
      Call Read_core_jj(nuc)
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

      end if

      Call MPI_BCAST(EC,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Ecore,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ncore, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(core,  250,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(closed,250,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      End Subroutine Def_core

