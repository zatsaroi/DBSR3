!======================================================================
      Real(8) Function Def_rrms(a,c)
!======================================================================
!     define nuclear potential for p- and q-splines
!---------------------------------------------------------------------
      Use zconst
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer :: i,j
      Real(8), intent(in) :: a,c
      Real(8) :: S,SS

      Do i=1,nv
       Do j=1,ks
        ygw(i,j) = one/(one+exp((gr(i,j)-c)/a)) * grw(i,j) * gr(i,j) * gr(i,j)
       End do
      End do
      S =  SUM(ygw(:,:))

      Do i=1,nv
       Do j=1,ks
        ygw(i,j) = ygw(i,j) * gr(i,j) * gr(i,j)
       End do
      End do
      SS =  SUM(ygw(:,:))

      Def_rrms  = sqrt(SS/S)

      End Function Def_rrms
