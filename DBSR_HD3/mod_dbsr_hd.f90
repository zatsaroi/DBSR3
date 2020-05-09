!=====================================================================
      Module dbsr_hd
!=====================================================================
!     Contains common variable and arrays used in program dbsr_hd
!---------------------------------------------------------------------
      Use DBS_grid
      Use conf_jj
      Use target_jj
      Use channel_jj
      Use zconst, only: c_au

      Implicit none

! ... files:

      Integer, parameter :: ma = 80;  Character(ma) :: AF
      Integer :: pri = 1;  Character(ma) :: AF_nnn  = 'dbsr_hd.nnn'
      Integer :: nup = 7;  Character(ma) :: AF_par  = 'dbsr_par'
      Integer :: nut = 8;  Character(ma) :: AF_tar  = 'target_jj'
      Integer :: nui = 9;  Character(ma) :: AF_int  = 'dbsr_mat.nnn'
      Integer :: nub = 11; Character(ma) :: AF_b    = 'bound.nnn'
      Integer :: nuu = 12; Character(ma) :: BF_b    = 'dbound.nnn'
      Integer :: nuh = 13; Character(ma) :: AF_h    = 'h.nnn'
      Integer :: nuw = 14; Character(ma) :: AF_w    = 'w.nnn'
      Integer :: nuc = 15; Character(ma) :: AF_cfg  = 'cfg.nnn'
      Integer :: nue = 16; Character(ma) :: AF_exp  = 'thresholds'
      Integer :: nur = 17; Character(ma) :: AF_rsol = 'rsol.nnn'

! ... miscellaneous:

      Real(8), parameter :: zero = 0.d0, one = 1.d0

! ... arguments:

      Integer :: itype  = 0        ! type of calculations
      Integer :: msol   = 0        ! max.number of solutions
      Real(8) :: Emax   = zero     ! max. energy for output
      Real(8) :: Emin   = zero     ! min. energy for output
      Real(8) :: RA, RB = zero     ! R-matrix radius and border derivative

      Integer :: klsp=0            ! index of partial wave
      Integer :: klsp1=1,klsp2=1   ! range of partial waves
      Character(3) :: ALSP ='nnn'

! ... main parameters of calculations:

      Integer :: mhm               ! size of interaction matrix
      Integer :: khm               ! number of solutions in new basis (nsol)
      Integer :: nsol              ! number of solutions in channel basis 
      Integer :: ksol              ! number of channel eigenvalues
      Integer :: diag_ovl          ! flag for nonzero ovl.blocks

! ... main global arrays:

      Real(8), allocatable :: a(:,:)      ! interaction matrix
      Real(8), allocatable :: b(:,:)      ! overlap matrices
      Real(8), allocatable :: v(:)        ! solution vector

      Real(8), allocatable :: eval(:)     ! eigenvalues

      Real(8), allocatable :: bb(:,:)     ! new basis
      Real(8), allocatable :: bval(:)     ! basis eigenvalues
      Integer, allocatable :: ipsol(:)    ! pointer on channel blocks
                                          ! in new basis
      Integer, allocatable :: isol(:)     ! pointer on main configuration

      Real(8), allocatable :: WMAT(:,:)   ! surface amplitudes

! ... exp.energies:

      Integer :: iexp = 0, iiexp = 0
      Real(8), allocatable :: E_exp(:)
      Integer, allocatable :: ip_exp(:)
      Real(8) :: au_eV, au_cm
      Character(2) :: unit = 'au'

      Integer,parameter :: nlab = 64
      Character(nlab) ::  Lab

! ... additonal output:

      Integer :: iwt   = -1       ! print channel weights in w.nnn
      Real(8) :: cwt   = -0.01d0  ! cut of for weights' printing in log file

      End Module dbsr_hd


