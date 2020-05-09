!======================================================================
      Module term_exp
!======================================================================
!     Containes the term-dependent coefficients of det.expansion
!     of two given conf.symmetries under consideration.
!----------------------------------------------------------------------
      Implicit none

! ... lists of ang.symmetries (1:kt)

      Integer :: kt1,kt2
      Integer, allocatable :: IP_kt1(:),IP_kt2(:)
      Integer, allocatable :: IP_kt12(:,:)

! ... pointer on non-zero determinants (1:ne,1:kdt)

      Integer :: kdt1,kdt2
      Integer, allocatable :: IP_det1(:,:),IP_det2(:,:)

! ... term-dependent det.expension coefficients (1:kt,1:kdt)

      Real(8),  allocatable :: C_det1(:,:), C_det2(:,:)

! ... determinants under consideration:

      Integer :: kd1,kd2

! ... configurations under consideration:

      Integer :: ic,jc, ic_case

      End Module term_exp



