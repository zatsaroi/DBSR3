!======================================================================
      Module  dbsr_dmat
!======================================================================
!     containes some common parameters, variables and arrays
!----------------------------------------------------------------------
      Use zconst, only: time_au, c_au
      Use DBS_grid

      Implicit none

! ... file units:

      Integer, parameter :: ma=80
      Character(ma) :: AF, name1,name2
      Integer :: iname1, iname2

      Integer :: pri  = 66;  Character(ma) :: AF_log = 'dbsr_dmat.log'

      Integer :: nut  =  7;  Character(ma) :: AF_tar = 'target_jj'
      Integer :: nub  =  3;  Character(ma) :: AF_bnk = 'mult_bnk_E1'
      Integer :: nur  =  9;  Character(ma) :: AF_rsol= 'rsol.nnn'

      Integer :: nuw  =  8;  Character(ma) :: AF_bsw = 'target.bsw'
                             Character(ma) :: BF_bsw = 'pert_nnn.bsw'
                             Character(ma) :: CF_bsw = 'name.bsw'
      
      Integer :: nuc1 =  1;  Character(ma) :: AF_c   = 'name.c'
      Integer :: nuc2 =  2;  Character(ma) :: BF_c   = 'cfg.nnn'

      Integer :: nub1 = 11;  Character(ma) :: AF_b   = 'name.j'
      Integer :: nub2 = 12;  Character(ma) :: BF_b   = 'dbound.nnn'
 
      Integer :: nud  =  4;  Character(ma) :: AF_d   = 'd.nnn'
      Integer :: nuv  =  4;  Character(ma) :: AF_dv  = 'dv.nnn'
      Integer :: nuf  = 15;  Character(ma) :: AF_zf  = 'zf_res'

      Integer :: nus  = 99     ! scratch files

! ... transition type:
      Integer ::      kpol = 1
      Character(1) :: ktype='E'
      Character(2) :: atype='E1'

! ... initial/final state mode:
      Character(1) :: ctype1='c', ctype2='c'

! ... number of configurations, channels, perturters, ... :
      Integer :: ilsp1,ilsp2, nch1,nch2, ncp1,ncp2, npert1,npert2, kdm1,kdm2
      Integer :: ncfg1,ncfg2, nwf1,nwf2, kset1,kset2
      Integer :: jot1=0, jot2=0, parity1=0, parity2=0
      Character(3) :: ALS1, ALS2

! ... dipole matrices:
      Real(8), allocatable :: DL   (:,:),DV   (:,:)    ! state    matrices
      Real(8), allocatable :: dipLp(:,:),dipVp(:,:)    ! B-spline matrixes
      Real(8), allocatable :: dipLq(:,:),dipVq(:,:)    ! B-spline matrixes

! ... number of solutions:
      Integer :: nstate1=0, nstate2=0
      Integer :: mstate1=0, mstate2=0
      Integer :: istate1=0, istate2=0

! ... expansion coefficients:
      Real(8), allocatable :: C1(:), C2(:)

! ... tolerences:
      Real(8), parameter :: eps_c    = 1.d-8      
      Real(8), parameter :: eps_ndet = 1.d-8
      Real(8), parameter :: eps_ovl  = 1.d-8      

! ... packing number for overlaps:
      Integer :: ibo = 2**15

! ... additional outout for polarizabilities:
      Integer :: ialpha = 0

! ... output gf or f-values:
      Character(1) :: gf = 'f'

! ... output of full dipole matrix:
      Character(1) :: dd = 'n'

! ... additional outout for polarizabilities:
      Integer :: debug = 0

      End Module dbsr_dmat

