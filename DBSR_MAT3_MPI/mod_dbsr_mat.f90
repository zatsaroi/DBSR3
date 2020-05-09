!----------------------------------------------------------------------
      Module dbsr_mat
!----------------------------------------------------------------------
!     main parameters for the BSR_MAT program
!----------------------------------------------------------------------
      Use zconst, only: c_au, time_au

      Use DBS_grid;          Use target_jj
      Use DBS_gauss;         Use orb_jj
      Use DBS_orbitals_pq;   Use conf_jj
                             Use channel_jj
      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma = 80    ! limit for file-name length
      Character(ma) :: AF     

      Integer :: prj =  3;  Character(ma) :: AF_prj = 'dbsr_mat.log'                            
      Integer :: pri = 66;  Character(ma) :: AF_pri = 'mat_log.nnn'

      Integer :: nut = 10;  Character(ma) :: AF_tar = 'target_jj'
      Integer :: nup = 11;  Character(ma) :: AF_par = 'dbsr_par'
      Integer :: nub = 12;  Character(ma) :: AF_bnk = 'int_bnk.nnn'
      Integer :: nuc = 13;  Character(ma) :: AF_cfg = 'cfg.nnn'
      Integer :: nuw = 14;  Character(ma) :: AF_bsw = 'target.bsw'
      Integer :: nui = 15;  Character(ma) :: AF_mat = 'dbsr_mat.nnn'
      Integer :: nud = 16;  Character(ma) :: AF_deb = 'int_mat.nnn'
      Integer :: nuo = 17;  Character(ma) :: AF_orb = 'target_orb'
      Integer :: nur = 18;  Character(ma) :: AF_corr= 'int_corr'
      Integer :: nue = 19;  Character(ma) :: AF_exp = 'thresholds'
      Integer :: nun = 20;  Character(ma) :: AF_new = 'target_new'

! ... core parameters: 

      Integer :: nclosed = 0 
      Real(8) :: Ecore = 0.d0, EC = 0.d0

! ... range of partial waves:

      Integer :: klsp, klsp1=1, klsp2=1      
      Character(3) :: ALSP

! ... Breit corrections:

      Integer :: mbreit  =   0      

! ... tolerence parameters:

      Real(8) :: eps_c    = 1.0D-10  !  tolerance for coefficients
      Real(8) :: eps_det  = 1.0D-10  !  tolerance for determinants
      Real(8) :: eps_ovl  = 1.0D-10  !  tolerance for overlaps
      Real(8) :: eps_sym  = 1.0D-10  !  tolerance for symmetry
      Real(8) :: eps_acf  = 1.0D-07  !  tolerance for asympt. coeff.s
      Real(8) :: eps_tar  = 1.0D-06  !  tolerance for target energies

      Real(8) :: S_ovl    = 0.75d0   !  tolerance for TOTAL overlaps

! ... packing basis for orbitals:

      Integer, parameter :: ibo = 2**15   ! we should choose !!!

! ... size of buffer for coefficients:

      Integer :: mcbuf = 1000000, ncbuf = 0
      Real(8), allocatable :: CBUF(:)
      Integer, allocatable :: itb(:), jtb(:), intb(:), idfb(:)
      Real(8) :: mem_buffer

! ... target-orthogonality flag:

      Integer :: iitar  =  0      

! ... skip the target checking:

      Integer :: check_target =  1      

! ... restrictions on the channel energies:

      Real(8) :: Edmax  = 1.d10           ! max. energy for channel eigenvalues
      Real(8) :: Edmin  = -2*c_au*c_au    ! min. energy for channel eigenvalues
      Real(8) :: Egap   = 0.001d0         ! tolerance for zero-energy solutions

! ... zero conditions:

      Integer :: ilzero = -1              ! # B-splines deleted at r=0  for P
      Integer :: ibzero =  1              ! # B-splines deleted at r=a
      Integer :: jlzero = -1              ! # B-splines deleted at r=0  for Q
      Integer :: jbzero =  1              ! # B-splines deleted at r=a

! ... structure of data:

      Integer :: icase =  0               ! current case    
      Integer, parameter :: ncase =  3    ! number of different integrals
      Character aint(0:ncase)/'O','L','R','S'/

      Integer, parameter :: ntype_O = 4
      Integer, parameter :: ntype_L = 3
      Integer, parameter :: ntype_R = 4
      Integer, parameter :: ntype_S = 20

      Integer :: atype =  0               ! current rel. integral    
      Integer, parameter :: ibtype = 10

! ... dimension limits:

      Integer :: mk    =      7           ! max. multipole index
      Integer :: nblock  = 2000           ! number of blocks in c_data       
      Integer :: mblock  = 3000           ! size of blocks
      Integer :: kblock  =  500           ! max.nb for given case

! ... debug printing and timing:

      Real(8) :: t_rdata=0.d0, t_sdata=0.d0, t_ldata=0.d0
      Integer :: debug    = 0      
      Integer :: pri_acf  = 0
      Integer :: pri_coef = 0

! ... MPI

      Integer :: nprocs=1, myid=0, ierr=0  

!----------------------------------------------------------------------
!     basic matrixes:
!     hch,hcp,hp  -  interaction (or overlap) channel matrixes 
!     acf         -  array of asymptotic coefficients
!     htarg       -  interaction matrix for target states
!     otarg       -  overlap matrix for target states
!----------------------------------------------------------------------

      Integer :: nhm    ! number of solutions 
      Integer :: mhm    ! matrix dimension 

      Real(8), allocatable :: hch(:,:,:), hcp(:,:), hp(:)
      Real(8), allocatable :: acf(:,:), bcf(:,:,:)
      Real(8), allocatable :: htarg(:), otarg(:)
      Integer, allocatable :: itarget(:)

      Real(8), allocatable :: x(:,:),xx(:,:)   !  ms * ms
      Real(8), allocatable :: y(:,:),yy(:,:)   !  ns * ns

      Integer :: idiag     ! =-1,0,1 (none-diagonal, all, diagonal blocks) 
      Integer :: nsol      ! number of solutions from diagonal blocks

      Real(8), allocatable :: diag(:,:,:), eval(:)
      Integer, allocatable :: iprm(:), ipsol(:), jpsol(:)
      Real(8), allocatable :: overlaps(:,:)

      Integer :: iicc, iicb, iibb  !  dimensions
      Integer, allocatable :: icc(:,:), icb(:,:), ibb(:,:)

      Integer, allocatable :: ich_state(:), no_state(:)

      Real(8) :: mem_mat
      Real(8) :: t_check, t_det, t_add

      End Module dbsr_mat

