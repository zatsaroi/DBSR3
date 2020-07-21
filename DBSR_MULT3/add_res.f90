
!======================================================================
      Subroutine Add_res(CN)
!======================================================================
!     records results from 'coef_list' to unit 'nu'
!----------------------------------------------------------------------
      Use dbsr_mult;     Use coef_list
      Use term_exp;      Use symt_list
      Use ndef_list, only: IPF

      Implicit none
      Integer :: i,j, it,jt, k,k1,k2
      Real(8) :: C,CN
      Integer, external :: DEF_ij

! ... convert det.factors from ndef_list to common det_list: 

      k = 0; Do i=1,ncoef; k=k+idfc(i); End do 
      if(k.gt.0) then
       Call Ndet_Idet
       Do i=1,ncoef
        if(idfc(i).eq.0) Cycle
        j=idfc(i); idfc(i)=IPF(j)
       End do
      end if

! ... record the coef.s:

      k = 0
      Do k1=1,kt1; it=IP_kt1(k1) 
      Do k2=1,kt2; jt=IP_kt2(k2)  
       if(IP_kt12(k1,k2).eq.0) Cycle

       k = k + 1
       Do i = 1,ncoef
        C=coef(k,i)*CN; if(abs(C).lt.Eps_C) Cycle
        write(nui) C,it,jt,intc(i),idfc(i)
       End do

       IT_done(DEF_ij(it,jt))=1

      End do; End do 

      End Subroutine Add_res

