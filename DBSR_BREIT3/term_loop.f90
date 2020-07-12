!======================================================================
      Subroutine Term_loop
!======================================================================
!     Add to the common list (coef_list) the coeff.s between two
!     determinants (zoef_list) weighted with term-dependent factors.
!     Check also if the specific coefficient is needed to be added
!     to the bank.
!----------------------------------------------------------------------
      Use dbsr_breit
      Use term_exp;  Use zoef_list; Use coef_list

      Implicit none
      Integer :: i,k,k1,k2
      Real(8) :: C

! ... find term-dependent coefficients between two determinants:

      k = 0
      Do k1=1,kt1; Do k2=1,kt2
       if(IP_kt12(k1,k2).eq.0) Cycle
       k = k + 1
       cdtrm(k) = C_det1(k1,kd1)*C_det2(k2,kd2)
      End do; End do

! ... add final coefficients:

      Do i=1,nzoef
       C = Zoef(i); if(abs(C).lt.EPS_c) Cycle
       Ctrm = cdtrm * C
       int = IZ_int(i)
       idf = IZ_df(i)
       Call Add_coef
      End do
      nzoef = 0

      End Subroutine Term_loop

