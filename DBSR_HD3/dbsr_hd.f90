!======================================================================
!    PROGRAM       D B S R _ H D                              version 3
!
!               C O P Y R I G H T -- 2015
!
!    Written by:   Oleg Zatsarinny
!                  email: oleg_zoi@yahoo.com
!======================================================================
!   Generates H.DAT and optioanlly W.DAT, RSOL.DAT files for further
!   scattering or photoionization calculations in relativistic R-matrix
!   approach.  Another option - calculation of bound and pseudo-states.
!======================================================================
!   INPUT  ARGUMENTS:
!
!     itype - mode of calculations
!
!           = -1 - bound-structure calculations
!           =  0 - scattering calculations
!           =  1 - photoionization calculations
!
!     klsp1,klsp2  - indexes of partial waves under consideration
!
!   INPUT FILES:
!
!     dbsr_par       -  description of partial waves
!     dbsr_mat.nnn   -  interaction matrix
!     cfg.nnn        -  c-file for close-coupling expansion
!
!    OUTPUT FILES:
!
!     dbsr_hd.log    -  running information
!     h.nnn          -  eigenvalues and suface amplitudes
!                       required to obtain R-matrix
!     w.nnn          -  file of channel weights (optional)
!     rsol.nnn       -  R-matrix inner-region solutions (optional)
!     bound.nnn      -  bound solusions
!
!     "nnn" indicates the index of partial wave
!=====================================================================
      Use dbsr_hd

      Implicit none
      Real(8) :: t0, t1

! ... read target information:

      Call Check_file(AF_tar); Open(nut,file=AF_tar)
      Call Read_target_jj(nut)

! ... set up B-splines:

      Call read_knot_dat
      Call alloc_DBS_gauss

! ... read arguments:

      Open(nup,file=AF_par);  Call Read_arg(nup);  Close(nup)

! ... cycle over partial waves:

      Do klsp = klsp1,klsp2

       Call CPU_time(t0);   Call SUB1_HD;   Call CPU_time(t1)

       write (pri,'(/a,f10.2,a)') 'time = ', (t1-t0)/60, ' min.'
       write (*  ,'(/a,f10.2,a)') 'time = ', (t1-t0)/60, ' min.'

      End do

      End  ! program DBSR_HD


