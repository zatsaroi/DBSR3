!======================================================================
      MODULE  param_jj
!======================================================================
!     main parameters and input/outpu files used in the program
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Save
!----------------------------------------------------------------------
!     names and units for input/output files:
!     AF  -  standard (default) names
!     BF  -  names with indication of partial wave number

      Character(40) :: AF,BF

! ... runing information:

      Integer :: pri=7;  Character(40) :: AF_pri = 'bj.log'

! ... c-file:

      Integer :: nuc=1;  Character(40) :: AF_c = 'rcsl.inp'
                         Character(40) :: BF_c = 'cfg.nnn'
! ... data bank:

      Integer :: nub=2;  Character(40) :: AF_b = 'int_bnk'
                         Character(40) :: BF_b = 'int_bnk.nnn'
      Logical :: new     ! pointer on the previous calculation  

! ... new results if any:

      Integer :: nur=3;  Character(40) :: AF_r = 'int_res'
                         Character(40) :: BF_r = 'int_res.nnn'
      Logical :: icalc   ! pointer for need of new calculations

! ... scratch files:

      Integer :: nui=11;  Character(40) :: AF_i = 'int_new'

      Integer :: nud=12  ! for det. expansions
      Integer :: nua=13  ! for accumulation of data

!----------------------------------------------------------------------

      Real(8) :: eps_c = 1.d-7    !  tolerence for coefficients
      Integer :: mk = 7           !  maximum multipole index
      Integer :: RX = 0           !  R or X integrals:
      Integer :: klsp1=0          !  range of partial waves
      Integer :: klsp2=0          !  range of partial waves
      Integer :: klsp=0           !  range of partial waves
      Integer :: debug=0          !  range of partial waves

      Integer, parameter :: isd = 20000  ! initial dimension for det.overlaps
      Integer, parameter :: jsd = 3      ! avarage size of overlap det.

      Integer, parameter :: isf = 200000 ! initial dimension for det.factors
      Integer, parameter :: jsf = 5      ! avarage number of dets. in def.

! ... initial (supposed) number of coef.s:

      Integer, parameter :: iszoef = 2000    ! in module ZOEF
      Integer, parameter :: iscoef = 2500    ! in module COEF
      Integer, parameter :: isboef = 50000   ! in module BOEF
      Integer, parameter :: isblk  = 5000    ! in module BOEF   


! ... MPI:

      Integer :: nprocs, myid, ierr
      Integer, Allocatable :: ip_proc(:)


      End MODULE  param_jj


