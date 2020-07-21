!======================================================================
      Subroutine Add_integral (kpol,ic,jc,i1,i2,i3,i4,C,mbreit)
!======================================================================
!     add two-electron two-component integral to the common 
!     list as the sum of one-component integrals
!----------------------------------------------------------------------
      Use dbsr_ci, only: ibi
      Use orb_jj,  only: kef,lef

      Implicit none
      Integer, intent(in) :: kpol,ic,jc,i1,i2,i3,i4,mbreit
      Real(8), intent(in) :: C
      Integer :: m,v
      Real(8) :: S(8)
      Integer, external :: l_kappa
      Real(8), external :: SMU

! ... Rk-integrals  

      m = mod(lef(i1)+lef(i3)+kpol,2) + mod(lef(i2)+lef(i4)+kpol,2)
      if(m.eq.0) then
       Call Add_coef(C,kpol,i1*ibi+i3,i2*ibi+i4,ic,jc,1)
       Call Add_coef(C,kpol,i1*ibi+i3,i2*ibi+i4,ic,jc,2)
       Call Add_coef(C,kpol,i1*ibi+i3,i2*ibi+i4,ic,jc,3)
       Call Add_coef(C,kpol,i1*ibi+i3,i2*ibi+i4,ic,jc,4)
      end if

! ... Sk-integrals

      if(mbreit.eq.0) Return
      Do v = kpol-1,kpol+1
       if(SMU(KEF(i1),KEF(i2),KEF(i3),KEF(i4),kpol,v,S).eq.0.d0) Cycle
       Call Add_coef(C*S(1),v,i1*ibi+i3,i2*ibi+i4,ic,jc,5)
       Call Add_coef(C*S(2),v,i2*ibi+i4,i1*ibi+i3,ic,jc,5)
       Call Add_coef(C*S(3),v,i3*ibi+i1,i4*ibi+i2,ic,jc,5)
       Call Add_coef(C*S(4),v,i4*ibi+i2,i3*ibi+i1,ic,jc,5)
       Call Add_coef(C*S(5),v,i1*ibi+i3,i2*ibi+i4,ic,jc,6)
       Call Add_coef(C*S(6),v,i4*ibi+i2,i3*ibi+i1,ic,jc,6)
       Call Add_coef(C*S(7),v,i3*ibi+i1,i4*ibi+i2,ic,jc,6)
       Call Add_coef(C*S(8),v,i2*ibi+i4,i1*ibi+i3,ic,jc,6)
      End do

      End Subroutine Add_integral


