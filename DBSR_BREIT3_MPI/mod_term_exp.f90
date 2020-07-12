!======================================================================
      MODULE term_exp
!======================================================================
!     Containes the term-dependent coefficients of det.expansion
!     of two given conf.symmetries under consideration.
!----------------------------------------------------------------------
      Implicit none

! ... lists of ang.symmetries (1:kt)

      Integer :: kt1,kt2  
      Integer, Allocatable :: IP_kt1(:), IP_kt2(:)   
      Integer, Allocatable :: IP_kt12(:,:)   

      Integer :: jt1,jt2 
      Integer, Allocatable :: JP_kt1(:), JP_kt2(:)   
 
! ... pointer on non-zero determinants (1:ne,1:kdt)

      Integer :: kdt1,kdt2  
      Integer, Allocatable :: IP_det1(:,:), IP_det2(:,:)   

! ... term-dependent det.expension coefficients (1:kt,1:kdt)

      Real(8),  Allocatable :: C_det1(:,:), C_det2(:,:)

! ... determinants under consideration: 

      Integer :: kd1,kd2    

! ... configurations under consideration: 

      Integer :: ic,jc, ic_case    
      
      End MODULE term_exp



