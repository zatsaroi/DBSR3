!======================================================================
      Subroutine Add_res(nu)
!======================================================================
!     records results from 'coef_list' to unit 'nu'
!----------------------------------------------------------------------

      USE param_jj
      USE coef_list, only: ncoef,idfc,intc,coef,nc
      USE term_exp
      USE ndef_list, ONLY: IPF
      USE symt_list

      Implicit none
      Integer, Intent(in) :: nu
      Integer :: i,j,ij, it,jt, k,k1,k2, n

! ... convert det.factors from ndef_list to common det_list: 

      k = 0; Do i=1,ncoef; k=k+idfc(i); End do 
      if(k.gt.0) then
       Call Ndet_Idet
       Do i=1,ncoef
        if(idfc(i).eq.0) Cycle
        j = idfc(i); idfc(i) = IPF(j)
       End do
      end if

! ... record the coef.s:

      k = 0
      Do k1=1,jt1; it=JP_kt1(k1) 
      Do k2=1,jt2; jt=JP_kt2(k2)  
       if(IP_kt12(k1,k2).eq.0) Cycle

       k = k + 1
       Do i = 1,ncoef
        if(abs(coef(k,i)).lt.Eps_C) Cycle
        write(nu) coef(k,i),it,jt,intc(i),idfc(i)
        nc=nc+1
       End do

       i=max(it,jt); j=min(it,jt); IT_done(i*(i-1)/2+j)=1

      End do; End do 

      End Subroutine Add_res


