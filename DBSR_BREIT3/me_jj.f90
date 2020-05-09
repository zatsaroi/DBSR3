!======================================================================
      Subroutine me_jj(l1,j1,m1,l2,j2,m2,l3,j3,m3,l4,j4,m4,ide)
!======================================================================
!     computes angular coefficients for the two-electron Coulomb
!     interaction in uncouple njm-representation;
!
!     should be called twice with interchange 3 <-> 4 and  
!     ide = +1  ->  direct interaction
!     ide = -1  ->  exchange interaction
!
!     coefficients are stored in module "boef_list" 
!     through call to Iadd_boef
!----------------------------------------------------------------------
      Use dbsr_breit, only: mbreit

      Implicit none
      Integer :: l1,j1,m1,l2,j2,m2,l3,j3,m3,l4,j4,m4,ide
      Integer :: KL,KM, KL1,KL2, KM1,KM2, k,m, kz,int, kk
      Real(8) :: C
      Real(8), external :: Z_3j2, Cjkj

      m = m1-m3; !if(m.ne.m4-m2) Return 

! ... define the range of multipole indeces:
    
      KL1=IABS(j1-j3);  KL2=IABS(j2-j4);  KL=MAX0(KL1,KL2)/2
      KM1=     j1+j3 ;  KM2=     j2+j4 ;  KM=MIN0(KM1,KM2)/2
      if(KM.lt.KL) Return

      Do k = kl,km; kk=k+k

       if(mbreit.eq.0) then
        if(mod(l1+l3+k,2).ne.0) Cycle
        if(mod(l2+l4+k,2).ne.0) Cycle
       end if

       C = Z_3j2(j1,-m1,kk,m1-m3,j3,m3) * Z_3j2(j2,-m2,kk,m2-m4,j4,m4)    
       if(C.eq.0.d0) Cycle
       kz = (kk-m + j1-m1 + j2-m2)/2
       if(mod(kz,2).ne.0) C = -C  
       C = C * ide

       C = C * Cjkj(j1,k,j3) * Cjkj(j2,k,j4)  
       if(mod(k,2).eq.1) C=-C

       int = (k+1)*ide

       Call Iadd_boef(C,int)

      End do

      End Subroutine me_jj

