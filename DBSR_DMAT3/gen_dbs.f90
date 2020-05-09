!======================================================================
      Subroutine Gen_dbs(nbf1,nbf2,wfact,kpol,ktype)      
!======================================================================
!     one-electron dipole integrals between two subsets of orbitals;
!----------------------------------------------------------------------
      Implicit none
      Integer :: nbf1,nbf2,kpol,i,j
      Character(1) :: ktype
      Real(8) :: wfact, SL, SV
      Real(8), external :: dipM_pq

      Call Alloc_multipole_integrals(0)

      if(ktype.eq.'M') then
       Do i=1,nbf1; Do j=nbf1+1,nbf1+nbf2
        SL=dipM_pq(i,j,kpol,wfact); SV = SL
        if(SL.eq.0.d0) Cycle
        Call Iadd_dip(i,j,SL,SV)
       End do; End do
      else
       Do i=1,nbf1; Do j=nbf1+1,nbf1+nbf2
        Call dipE_pq(i,j,kpol,wfact,SL,SV)
        if(SL.eq.0.d0.and.SV.eq.0.d0) Cycle
        Call Iadd_dip(i,j,SL,SV)
       End do; End do
      end if

      End Subroutine Gen_dbs


