!======================================================================
      Subroutine Check_ovl(ncfg_jj,ncfg_LS)
!======================================================================

      Use jj2ls

      Implicit none

      Integer, intent(in) :: ncfg_jj,ncfg_LS
      Integer :: i,j,k, ncj

      read(nuo) ncj
      if(ncj.ne.ncfg_jj) then
       iovl=0
       write(pri,*) 'ncj <> ncfg_jj for overlap matrix'
       write(pri,*) 'we ignore the overlap matrix' 
       Return
      end if

      if(allocated(S_ovl)) Deallocate(S_ovl)
      Allocate(S_ovl(ncj,ncj))
      Do i=1,ncj;  read(nuo) S_ovl(i,1:i); End do
      Do i=1,ncj-1; Do j=2,i-1; S_ovl(j,i)=S_ovl(i,j); End do; End do

      if(allocated(C_ovl)) Deallocate(C_ovl)
      Allocate(C_ovl(ncfg_LS,ncfg_LS))
      if(allocated(CJ)) Deallocate(CJ); Allocate(CJ(ncfg_jj))
      if(allocated(CI)) Deallocate(CI); Allocate(CI(ncfg_jj))
      
      Do i=1,ncfg_LS; CI(:)=C_trans(:,i)
       Do k=1,ncfg_jj; CJ(k)=SUM(CI(:)*S_ovl(:,k)); End do
       Do j=1,i
        C_ovl(i,j)=SUM(CJ(:)*C_trans(:,j)) 
        C_ovl(j,i)=C_ovl(i,j)
       End do
      End do

      Deallocate(S_ovl,CI,CJ)

      End Subroutine Check_ovl
