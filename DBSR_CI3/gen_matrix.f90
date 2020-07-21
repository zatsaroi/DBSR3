!======================================================================
      Subroutine Gen_matrix (itype,kpol)
!======================================================================
! ... process the two-electron integrals:
!----------------------------------------------------------------------

      if(itype.le.4) then
       Call Rk_evaluate(kpol,itype) 
      elseif(itype.le.6) then
       Call Sk_evaluate(kpol,itype) 
      end if

      End Subroutine Gen_matrix