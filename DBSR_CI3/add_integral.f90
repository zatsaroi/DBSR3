!======================================================================
      Subroutine Add_integral (kpol,ic,jc,i1,i2,i3,i4,C,mbreit)
!======================================================================
!     add integral two-electron two-component integral to the common 
!     list as the sum of one-component integrals
! ??? we have different numbering of integrals in different programs
!----------------------------------------------------------------------
      Use orb_jj,  only: kef,lef
      Use c_data,  only: ibi

      Implicit none
      Integer, intent(in) :: kpol,ic,jc,i1,i2,i3,i4,mbreit
      Real(8), intent(in) :: C
      Integer :: m,l1,l2,l3,l4, itype,v
      Real(8) :: S(8)
      Integer, external :: l_kappa
      Real(8), external :: SMU

! ... Rk-integrals  

      m = mod(lef(i1)+lef(i3)+kpol,2) + mod(lef(i2)+lef(i4)+kpol,2)
      if(m.eq.0) then
       Call Add_coef(kpol,1,C,i1*ibi+i3,i2*ibi+i4,ic,jc)
       Call Add_coef(kpol,2,C,i1*ibi+i3,i2*ibi+i4,ic,jc)
       Call Add_coef(kpol,3,C,i1*ibi+i3,i2*ibi+i4,ic,jc)
       Call Add_coef(kpol,4,C,i1*ibi+i3,i2*ibi+i4,ic,jc)
      end if

! ... Sk-integrals

      if(mbreit.eq.0) Return
      Do v = kpol-1,kpol+1
       if(SMU(KEF(i1),KEF(i2),KEF(i3),KEF(i4),kpol,v,S).eq.0.d0) Cycle
       Call Add_coef(v,5,C*S(1),i1*ibi+i3,i2*ibi+i4,ic,jc)
       Call Add_coef(v,5,C*S(2),i2*ibi+i4,i1*ibi+i3,ic,jc)
       Call Add_coef(v,5,C*S(3),i3*ibi+i1,i4*ibi+i2,ic,jc)
       Call Add_coef(v,5,C*S(4),i4*ibi+i2,i3*ibi+i1,ic,jc)
       Call Add_coef(v,6,C*S(5),i1*ibi+i3,i2*ibi+i4,ic,jc)
       Call Add_coef(v,6,C*S(6),i4*ibi+i2,i3*ibi+i1,ic,jc)
       Call Add_coef(v,6,C*S(7),i3*ibi+i1,i4*ibi+i2,ic,jc)
       Call Add_coef(v,6,C*S(8),i2*ibi+i4,i1*ibi+i3,ic,jc)
      End do

      End Subroutine Add_integral


