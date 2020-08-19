!======================================================================
      Subroutine Add_res(nu)
!======================================================================
!     records results from 'coef_list' to unit 'nu'
!----------------------------------------------------------------------
      Use dbsr_breit
      Use coef_list; Use term_exp; Use symt_list
      Use ndef_list, only: IPF

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j, it,jt, k,k1,k2, ij
      Integer, external :: DEF_ij

! ... convert det.factors from ndef_list to common def_list:

      k = 0; Do i=1,ncoef; k=k+idfc(i); End do
      if(k.gt.0) then
       Call Ndet_Idet
       Do i=1,ncoef
        if(idfc(i).eq.0) Cycle
        j = idfc(i); idfc(i) = IPF(j)
       End do
      end if

! ... record coefficients:

      k = 0
      Do k1=1,kt1; it=IP_kt1(k1)
      Do k2=1,kt2; jt=IP_kt2(k2)
       if(IP_kt12(k1,k2).eq.0) Cycle

       k = k + 1
       Do i = 1,ncoef
        if(abs(coef(k,i)).lt.Eps_C) Cycle
        write(nu) coef(k,i),it,jt,intc(i),idfc(i)
        nc=nc+1
       End do

       ij = DEF_ij(it,jt); IT_done(ij)=1

      End do; End do

      End Subroutine Add_res

