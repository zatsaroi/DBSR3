!======================================================================
      Module  dbsr_mult
!======================================================================
!     main parameters and variable for program dbsr_mult
!----------------------------------------------------------------------
      Implicit none

      Integer :: kpol,qpol,kkpol
      Character(1) :: ktype 

! ... default names and units for input/output files:

      Integer, parameter :: ma = 80; Character(ma) :: AF

      Integer :: pri=66; Character(ma) :: AF_log = 'dbsr_mult.log'
      Integer :: nu1=1; Character(ma) :: AF1    = 'name1.c' 
      Integer :: nu2=2; Character(ma) :: AF2    = 'name2.c' 
      Integer :: nub=3; Character(ma) :: AF_bnk = 'mult_bnk_E1'
      Integer :: nur=4; Character(ma) :: AF_res = 'mult_res'

      Logical :: new    ! new calculations or not 
      Logical :: icalc  ! need of additional calculations
      
! ... scratch files:

      Integer :: nui = 10  ! for intermediate results
      Integer :: nus = 11  ! for re-allocations   
      Integer :: nud = 12  ! for det.expansions
      Integer :: nua = 13  ! for accumulation of data

! ... tolerence for coefficients:

      Real(8) :: Eps_c = 1.d-7      

      End Module dbsr_mult

