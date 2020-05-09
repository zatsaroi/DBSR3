!====================================================================
      Subroutine me_jj(l1,j1,m1,l2,j2,m2,CA,CB)
!====================================================================
!     angular part of electric or magnetic transition operator 
!     between 'nljm' orbitals:
!
!            <n1,l1,j1,m1| T(kq) | n2,l2,j2,m2>
!
!--------------------------------------------------------------------
      Use dbsr_mult 

      Implicit none
      Integer, intent(in) :: l1,j1,m1,l2,j2,m2
      Real(8), intent(out) :: CA,CB	  
      Integer :: i
      Real(8), external :: Z_3j2, Cjkj

      CA = 0.d0; CB = 0.d0
      i=mod(l1+l2+kpol,2)
      if(ktype.eq.'E'.and.i.eq.1) Return 
      if(ktype.eq.'M'.and.i.eq.0) Return 

      CA = (-1)**((j1-m1)/2) * Z_3j2(j1,-m1,kkpol,qpol,j2,m2) 
      if(CA.eq.0.d0) Return 

      CA = CA * Cjkj (j1,kpol,j2)

      End Subroutine me_jj
