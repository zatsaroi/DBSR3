!=====================================================================
      Module jj2ls     
!=====================================================================
!     Contains common variable and arrays used in program jj2ls
!---------------------------------------------------------------------

      Implicit none

! ... files:

      Integer, parameter :: ma = 80;  Character(ma) :: AF, name

      Integer :: pri =  6; Character(ma) :: AF_log  = 'jj2ls.log'
      Integer :: nuc =  7; Character(ma) :: AF_jjc  = 'name.c'
      Integer :: nur =  8; Character(ma) :: AF_lsc  = 'name_LS.c'
      Integer :: nuj =  8; Character(ma) :: AF_jjj  = 'name.j'
      Integer :: nul =  9; Character(ma) :: AF_lsj  = 'name_LS.j'
      Integer :: num = 10; Character(ma) :: AF_cm   = 'name.m'
      Integer :: nuo = 11; Character(ma) :: AF_ovl  = 'name.ovl'

      Integer :: nua = 21    ! scratch file 

! ... re-coupling parameters and arrays:

      Integer, parameter :: mshells = 8
      Integer, parameter :: mcup=2*mshells  
      Integer, parameter :: mmom=6*mshells  
      Integer :: JJ_coupling(3,mcup)
      Integer :: LS_coupling(3,mcup)
      Integer :: moments(mmom)
      Integer :: ncup, nmom, key=0

! ... subshells arrays:

      Integer, parameter :: mn = mshells
      Integer :: n                             
      Integer :: ncase(mn)  
      Integer :: jot1(mn),jot2(mn),jot3(mn),jotp(mn),jot(mn),JT(mn)
      Integer :: jq1(mn),jq2(mn),kt1(mn),kt2(mn)

! ... working arrays:

      Real(8), Allocatable :: C_term(:,:), C_trans(:,:), C1(:),C2(:)  
      Real(8), Allocatable :: S_ovl(:,:), C_ovl(:,:), CJ(:), CI(:)

! ... other parameters:

      Integer :: JJ_nterms, JJ_term
      Integer :: LS_nterms, LS_term
      Integer :: iovl
      Real(8), parameter :: zero=0.d0, one=1.d0

      Integer :: debug = 0

      End module jj2ls  


