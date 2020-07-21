!======================================================================
      Subroutine Gen_matrix (itype,kpol)
!======================================================================
! ... generate matrix:
!----------------------------------------------------------------------
      Use dbsr_mat

      Select case(icase)
       Case(0);    Call O_data(itype) 
       Case(1);    Call L_data(itype) 
       Case(2);    Call R_data(itype,kpol,0) 
                   Call R_data(itype,kpol,1) 
                   Call R_data(itype,kpol,2) 
                   Call R_data(itype,kpol,3) 
       Case(3);    Call S_data(itype,kpol) 
      End Select

      End Subroutine Gen_matrix