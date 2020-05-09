!======================================================================
      Module nljm_orbitals
!======================================================================
!     contains rather fancy but effective description of two
!     determinant wave-functions under consideration
!----------------------------------------------------------------------
      Implicit none

! ... common list of different orbital symmetries for two determinants:
!
!     NSYM  - number of symmetries:  l,j,mj ->  Lsym,Jsym,Msym
!     IPSYM - pointer on the last given symmetry in the list
!     KSYM  - number of orbitals with given symmetry
!     nnsym - principal quantum numbers
!     Isym  - pointer on the original position in configuration

      Integer :: NSYM
      Integer, allocatable :: Lsym(:),Lsym1(:),Lsym2(:)
      Integer, allocatable :: Jsym(:),Jsym1(:),Jsym2(:)
      Integer, allocatable :: Msym(:),Msym1(:),Msym2(:)
      Integer, allocatable :: IPsym1(:),IPsym2(:)
      Integer, allocatable :: Ksym1(:),Ksym2(:)
      Integer, allocatable :: nnsym1(:),nnsym2(:)
      Integer, allocatable :: Isym1(:),Isym2(:)

! ... shell values:

!     md(i)  - the max.number of det.'s for the i-th subshell
!     nd(i)  - determinant under consideration
!     ipn(i) - pointers on the orbitals of given shell
!     MJs    - shells MJ
!     MJi    - intermediate values MJ

      Integer, allocatable :: md(:),nd(:),ipn(:),MJs(:),MJi(:)

      Integer :: kz1,kz2     ! number of perturbations

! ... ordered list of possible mj values:

      Integer :: mj_max = 0
      Integer, allocatable :: mj_orb(:)

! ... pointers to nljm-orbitals:

      Integer, allocatable :: Idet(:),Jdet(:)

! ... auxiliary arrays

      Integer, allocatable :: N1(:),N2(:),N3(:),N4(:), NP(:)

      End module nljm_orbitals


!======================================================================
      Subroutine Alloc_nljm(ne,msh)
!======================================================================
!     allocate arrays in module "nljm_orbitals" for the given case
!----------------------------------------------------------------------
      Use nljm_orbitals
      Implicit none
      Integer :: ne,msh,m
     
      if(allocated(Lsym)) Return
     
      m = 2*ne
      Allocate(Lsym(m),Lsym1(m),Lsym2(m), &
               Jsym(m),Jsym1(m),Jsym2(m), &
               Msym(m),Msym1(m),Msym2(m))
     
      Allocate(Isym1(m),Isym2(m),  &
               Ksym1(m),Ksym2(m),  &
               nnsym1(m),nnsym2(m),&
               IPsym1(0:m),IPsym2(0:m))
     
      IPsym1(0) = 0; IPsym2(0) = 0
     
      Allocate(Idet(m),Jdet(m))
      Allocate(N1(m),N2(m),N3(m),N4(m),NP(m))
     
      Allocate(md(msh),nd(msh),ipn(msh),MJs(msh),MJi(msh))
     
      End Subroutine Alloc_nljm
