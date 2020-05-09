!----------------------------------------------------------------------
      Module dbsr_ci
!----------------------------------------------------------------------
!     main parameters for the DBSR_ci program
!----------------------------------------------------------------------
      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma = 80       ! limit for file-name length
      Character(ma) :: AF     
      Character(ma) :: name = ' '         ! name of case
      Character(ma) :: knot = 'knot.dat'  

      Integer :: inp =  5; Character(ma) :: AF_inp = 'name.inp_ci'
      Integer :: pri = 66; Character(ma) :: AF_log = 'name.log_ci'
      Integer :: nuc = 11; Character(ma) :: AF_c   = 'name.c'
      Integer :: nub = 12; Character(ma) :: AF_bnk = 'name.bnk'
                           Character(ma) :: BF_bnk = 'int_bnk'
      Integer :: nuw = 13; Character(ma) :: AF_bsw = 'name.bsw'
      Integer :: nuj = 14; Character(ma) :: AF_j   = 'name.j'
      Integer :: nuo = 15; Character(ma) :: AF_ovl = 'name.ovl'
      Integer :: num = 16; Character(ma) :: AF_mat = 'name.mat'
      Integer :: nud = 17; Character(ma) :: AF_debug = 'name.debug'
      Integer :: nua = 20  ! scratch file

! ... atomic parameters: 

      Real(8) :: Z   = 0.d0           ! atomic number
      Real(8) :: awt = 0.d0           ! atomic weight
      Real(8) :: au_cm,au_eV          ! convertion constant
      Integer :: nclosed = 0          ! number of common closed shells

      Integer :: mpol = 7             ! maximum multipole index

! ... Breit corrections:

      Integer :: mbreit = 0      
      Integer :: msol = 0      
      Integer :: nsol = 0      
      Integer :: meiv = 0      
      Integer :: mdiag = 0      

      Real(8) :: Emax = 0.d0
      Real(8) :: Ecore = 0.d0

      Real(8), allocatable :: HM(:,:)    ! interaction matrix
      Real(8), allocatable :: SM(:,:)    ! overlap matrix
      Real(8), allocatable :: EVAL(:)    ! eigenvalues
      Real(8), allocatable :: DM(:)      ! diagonal m.e.
      Real(8), allocatable :: shift(:)   ! shift for diagonal H(i,i)

      Integer :: ncj   = 0       !  matrix dimension
      Integer :: nzero = 0       !  first-order set in config. list
      Integer :: nort  = 0       !  generelized or not (=0)
      Integer :: jot   = 0       !  total 2*J
      Integer :: nshift= 0       !  flag of additional shifts

! ... tolerence parameters:

      Real(8) :: Eps_det  = 1.d-10  !  tolerance for determinants
      Real(8) :: Eps_ovl  = 1.d-10  !  tolerance for overlaps
      Real(8) :: Eps_o    = 0.1d0   !  tolerance for deleting

      Integer :: debug = 0    ! additional output

! ... buffer size for angular coefficients:

      Integer, parameter :: maxnc = 1000000

! ... integral corrections:

      Integer :: ncorr = 0
      Integer, allocatable :: icorr(:,:)
      Real(8), allocatable :: scorr(:)

! ... check c-file option flag: 

      Integer :: check_c = 0

      End Module dbsr_ci

