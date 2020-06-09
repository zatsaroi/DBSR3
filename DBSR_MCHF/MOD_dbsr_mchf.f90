!======================================================================
      Module dbsr_mchf
!======================================================================
!     main parameters for the DBSR_MCHF program
!----------------------------------------------------------------------
      Use DBS_nuclear;  Use DBS_grid;   Use DBS_gauss
      Use conf_jj
      Use zconst

      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma = 100    ! limit for file-name length
      Character(ma) :: AF     

      Integer :: log = 3;   Character(ma) :: AF_log  = 'name.log'                            
      Integer :: inp = 5;   Character(ma) :: AF_dat  = 'name.inp'
      Integer :: scr = 6;

      Integer :: nuk = 10;  Character(ma) :: AF_knot = 'name.knot'
                            Character(ma) :: BF_knot = 'knot.dat'
      Integer :: nub = 11;  Character(ma) :: AF_bnk  = 'name.bnk'
                            Character(ma) :: BF_bnk  = 'int_bnk'
      Integer :: nuc = 12;  Character(ma) :: AF_cfg  = 'name.c'
      Integer :: nuw = 13;  Character(ma) :: AF_inp  = 'name.bsw'
                            Character(ma) :: AF_out  = 'name.bsw'
      Integer :: nuj = 14;  Character(ma) :: AF_j    = 'name.j'

! ... name of case:

      Character(ma) :: name = ' '

! ... atomic parameters: 

      Character(2)  :: atom = ' '
      Character(2)  :: ion  = ' '

      Real(8) :: Z   = 0.d0
      Real(8) :: atw = 0.d0
      Real(8) :: rms = 0.d0

      Real(8) :: Etotal, Ecore

! ... convergence:

      Real(8) :: scf_tol = 1.d-9,   scf_diff = 1.d0
      Real(8) :: orb_tol = 1.d-5,   orb_diff = 1.d0
      Real(8) :: end_tol = 1.d-6
      Real(8) :: eps_c   = 1.d-8    
      Real(8) :: eps_det = 1.d-8
      Real(8) :: eps_rot = 1.d-2
      Real(8) :: qsum_min= 1.d-8

      Integer :: max_it  = 25
      
      Integer :: n_corr  = 1
      Real(8) :: eps_corr= 1.d-1

! ... running options:

      Integer :: all     = 0
      Integer :: method  = 1
      Integer :: newton  = 0
      Integer :: rotate  = 0
      Integer :: irhs    = 0
      Real(8) :: srhs    = 0.d0
      Integer :: debug   = 0
      Integer :: icore   = 0

! ... block parameters: 

      Integer :: eol = 5                     !  = 1,5,9 - optimization mode 
                                             
      Integer :: nlevels = 0                 !  number of levels chosen for optimization
                                             !  each level has attritutes:
      Integer, allocatable :: block(:)       !  block index from (1:nlevel)
      Integer, allocatable :: level(:)       !  level index in the block   (1:nlevel)
      Integer, allocatable :: nlevelb(:)     !  number of optimized level in the block (1:njb)
      Real(8), allocatable :: weight(:)      !  level weight (1:nlevel)
      Real(8), allocatable :: elevel(:)      !  level energy (1:nlevel)
      Integer, allocatable :: ip_level(:)    !  pointer to the given level in array coeffs, see below
      Character(mlab), allocatable :: labeln(:)  ! level_label  (1:nlevel)
                                             
      Integer :: coefs_size = 0                 
      Real(8), allocatable :: coefs(:)       !  contains expansion coeff.s for all opt. levels
      Real(8) :: memory_mat = 0.d0 
                                             
! ... orbital variables:  

      Integer :: kmin= 0
      Integer :: kmax= 0

      Integer :: mtype = 2
      Integer :: nblock = 1000
      Integer :: mblock = 5000

      Integer :: varied = 0
      Character(ma) :: avaried = 'all'
      Integer, allocatable :: ivaried(:)

      Integer :: mbreit = 0      

      Integer :: ipzero = 0 
      Integer :: iqzero = 0 
      Integer :: jpzero = 1 
      Integer :: jqzero = 1 

! ... buffer size for coefficients:

      Integer, parameter :: maxnc = 1000000

! ... frequently called functions:   

      Integer, external :: Icheck_file

! ... debuging time:

      Real(8) :: time_solve=0.d0, time_matrix=0.d0, time_diag=0.d0, &
                 time_read_coef=0.d0

      Integer, parameter :: ibi = 2**15

! ... list of L-integrals:

      Integer :: Lint = 0                    ! number of integrals
      Real(8), allocatable :: L_int(:)       ! value of integrals
      Integer, allocatable :: i1_Lint(:)     ! orbital index
      Integer, allocatable :: i2_Lint(:)
      Integer :: Lcoef = 0                   ! number of angular coefficients
      Integer, allocatable :: ic_Lcoef(:)    ! configuration index
      Integer, allocatable :: jc_Lcoef(:)
      Real(8), allocatable :: L_coef(:)      ! angular coefficients
      Integer, allocatable :: ip_Lint(:)     ! pointer to integrals
      Integer(1), allocatable :: if_Lint(:)  ! fixed or not

! ... list of Rk-integrals:

      Integer :: nint = 0                    ! number of integrals
      Real(8), allocatable :: Rk_int(:)      ! value of integral
      Integer, allocatable :: i1_int(:)      ! orbital index
      Integer, allocatable :: i2_int(:)      ! 
      Integer, allocatable :: i3_int(:)      ! 
      Integer, allocatable :: i4_int(:)      ! 
      Integer :: ncoef = 0                   ! number of angular coefficients
      Integer, allocatable :: ic_coef(:)     ! first configuration index
      Integer, allocatable :: jc_coef(:)     ! second configuration index
      Real(8), allocatable :: Rk_coef(:)     ! angular coefficients
      Integer, allocatable :: ip_int (:)     ! pointers to the coefficient 
      Integer, allocatable :: nk_int (:)     ! pointer for integrals for given k-value
      Integer, allocatable :: nk_coef(:)     ! pointer for coefficients for given k-value
      Integer(1), allocatable :: if_int(:)   ! fixed or not

! ... list of Sk-integrals:

      Integer :: Sint = 0                    ! number of integrals
      Real(8), allocatable :: Sk_int (:)     ! value of integral
      Integer, allocatable :: i1_Sint(:)     ! orbital index
      Integer, allocatable :: i2_Sint(:)     ! 
      Integer, allocatable :: i3_Sint(:)     ! 
      Integer, allocatable :: i4_Sint(:)     ! 
      Integer :: Scoef = 0                   ! number of angular coefficients
      Integer, allocatable :: ic_Scoef(:)    ! first configuration index
      Integer, allocatable :: jc_Scoef(:)    ! second configuration index
      Real(8), allocatable :: Rk_Scoef(:)    ! angular coefficients
      Integer, allocatable :: ip_Sint (:)    ! pointers to the coefficient 
      Integer, allocatable :: nk_Sint (:)    ! pointer for integrals for given k-value
      Integer, allocatable :: nk_Scoef(:)    ! pointer for coefficients for given k-value
      Integer(1), allocatable :: if_Sint(:)  ! fixed or not

      Real(8) :: memory_coefs = 0.d0 

      Integer, allocatable :: iphys(:)       !  physical orbital 
      Character(200) :: physical = ' '

      End Module dbsr_mchf
