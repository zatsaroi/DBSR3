!=====================================================================
      Module dbsr_pol      
!=====================================================================
!     Contains common variable and arrays 
!---------------------------------------------------------------------
      Use DBS_grid

      Implicit none

! ... files:

      Integer, parameter :: ma = 80;  Character(ma) :: AF

      Integer :: pri =  3; Character(ma) :: AF_log  = 'dbsr_pol.nnn'
      Integer :: nup =  7; Character(ma) :: AF_par  = 'dbsr_par'
      Integer :: nut =  8; Character(ma) :: AF_tar  = 'target_jj'
      Integer :: nui =  9; Character(ma) :: AF_int  = 'dbsr_mat.nnn'
      Integer :: nuw = 10; Character(ma) :: AF_bsw  = 'target.bsw'
      Integer :: nub = 11; Character(ma) :: AF_bnd  = 'bound.nnn'
      Integer :: nuc = 15; Character(ma) :: AF_cfg  = 'cfg.nnn'
      Integer :: nud = 17; Character(ma) :: AF_dip  = 'dv.nnn'
      Integer :: nur = 18; Character(ma) :: AF_pol  = 'pol.nnn'
      Integer :: nuo = 19; Character(ma) :: AF_ort  = 'add_ortht.bsw'

      Integer :: nua = 21    ! scratch file for int.matrix
      Integer :: nus = 22    ! scratch file for ovl.matrix
      Integer :: nuq = 23    ! scratch file for orth.states

! ... partial waves under consideration:

      Integer :: klsp = 1  
      Character(3) :: ALSP='001'

! ... orth.constraints:

      Integer :: nort = 0               ! from c-file
      Integer :: nortb = 0              ! additional
      Integer, allocatable :: iortb(:)

! ... other parameters:
      
      Real(8) :: EC = 0.d0              ! core energy
      Real(8) :: E1 = 0.d0              ! initial state energy
      Integer :: jot1                   ! initial state 2J

      Integer :: ifail  = 0             ! fail to diagonalize
      Integer :: debug  = 0             ! debug output
      Real(8), parameter :: zero=0.d0, one=1.d0

! ... dipole vector:

      Integer :: kpol = 1
      Character(1) :: ktype = 'E' 
      Real(8), allocatable :: d(:) 

! ... main dimensions:

      Integer :: nsol              ! number of solutions in new basis 
      Integer :: nhm               ! nsol + npert
      Integer :: mhm               ! nhm + nort + nortb
      Integer :: khm               ! nch*ms + npert

! ... main arrays:

      Real(8), allocatable :: bb(:,:)     ! new basis,  (ms,nsol)
      Integer, allocatable :: ipsol(:)    ! pointer on channel blocks  (nsol)
                                          ! in new basis
      Real(8), allocatable :: bval(:)     ! eigenvalues (nsol)
      Real(8), allocatable :: om(:,:)     ! overlap matrix     (mhm,mhm)
      Real(8), allocatable :: hm(:,:)     ! interaction matrix (mhm,mhm)

      Real(8), allocatable :: v(:)        ! solution vector

      Integer :: diag_ovl          ! flag for nonzero ovl.blocks

      End module dbsr_pol

