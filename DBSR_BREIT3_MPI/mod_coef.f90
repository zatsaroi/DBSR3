!====================================================================
      MODULE coef_list
!====================================================================
!     Contains a set of coefficients with two identifiers (intc,idfc)
!     The list of ordered according the pointer 'ipcoef'
!--------------------------------------------------------------------
      IMPLICIT NONE 
    
      Integer :: ncoef = 0        !  number of coefficients
      Integer :: mcoef = 0        !  current allocation size  
      Integer :: icoef = 10000    !  increment for allocation
      Integer :: ntrm  = 0        !  number of terms

      REAL(8), Allocatable :: coef(:,:)         !  coeff.s (1:ntrm,1:mcoef)
      Integer, Allocatable :: intc(:),idfc(:)   !  their attributes
      Integer, Allocatable :: ipcoef(:)         !  ordering pointer 
      REAL(8), Allocatable :: ctrm(:), cdtrm(:) !  new coefficients
      Integer, Allocatable :: itc(:), jtc(:)    !  term pointers

      Integer :: int, idf                       !  current integral 
      Integer :: nc, nct = 0                    !  total number of coef.s

      End MODULE coef_list


!======================================================================
      Subroutine Alloc_trm(nt)
!======================================================================
      Use coef_list

      Implicit none
      Integer, Intent(in) :: nt

      if(nt.le.0) then
       if(Allocated(ctrm)) Deallocate(ctrm,cdtrm,itc,jtc)
       ntrm=0; Call Alloc_coef(0)
      elseif(.not.Allocated(ctrm)) then
       Allocate(ctrm(nt),cdtrm(nt),itc(nt),jtc(nt))
       ntrm=nt
      else
       Deallocate(ctrm,cdtrm,itc,jtc)
       Allocate(ctrm(nt),cdtrm(nt),itc(nt),jtc(nt))
       ntrm=nt; Call Alloc_coef(0)
      end if

      End Subroutine Alloc_trm


!======================================================================
      Subroutine Alloc_coef(mc)
!======================================================================
      Use coef_list

      Implicit none
      Integer, Intent(in) :: mc
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: ra(:,:)

      if(ntrm.le.0.and.mc.gt.0) Stop ' Alloc_coef: ntrm = 0 '

      if(mc.le.0) then
       if(Allocated(coef)) Deallocate(coef,intc,idfc,ipcoef)
       ncoef=0; mcoef=0
      elseif(.not.Allocated(coef)) then
       Allocate(coef(ntrm,mc),intc(mc),idfc(mc),ipcoef(mc))
       mcoef=mc
      elseif(mc.le.mcoef) then
       Return
      elseif(ncoef.eq.0) then
       Deallocate(coef,intc,idfc,ipcoef)
       Allocate(coef(ntrm,mc),intc(mc),idfc(mc),ipcoef(mc))
       mcoef=mc
      else
       Allocate(ra(ntrm,ncoef))
       ra(1:ntrm,1:ncoef)=coef(1:ntrm,1:ncoef); Deallocate(coef) 
       Allocate(coef(ntrm,mc)); coef(1:ntrm,1:ncoef)=ra(1:ntrm,1:ncoef)
       Deallocate(ra) 
       Allocate(ia(ncoef))
       ia(1:ncoef)=intc(1:ncoef); Deallocate(intc) 
       Allocate(intc(mc)); intc(1:ncoef)=ia(1:ncoef) 
       ia(1:ncoef)=idfc(1:ncoef); Deallocate(idfc) 
       Allocate(idfc(mc)); idfc(1:ncoef)=ia(1:ncoef) 
       ia(1:ncoef)=ipcoef(1:ncoef); Deallocate(ipcoef) 
       Allocate(ipcoef(mc)); ipcoef(1:ncoef)=ia(1:ncoef) 
       Deallocate(ia) 
       mcoef=mc
       write(*,*) ' realloc_coef: mcoef = ', mcoef
      end if
	        
      END Subroutine Alloc_coef


!======================================================================
      Subroutine Add_coef
!======================================================================
!     add new coefficient to the list 
!----------------------------------------------------------------------

      USE coef_list

      Implicit none

      Integer :: i,k,l,m,ip

! ... look for the same integral in the list

      k=1; l=ncoef 

    1 if(k.gt.l) go to 2              
      m=(k+l)/2; ip=ipcoef(m)

      if    (int.lt.intc(ip)) then;  l = m - 1
      elseif(int.gt.intc(ip)) then;  k = m + 1
      else

      if    (idf.lt.idfc(ip)) then;  l = m - 1
      elseif(idf.gt.idfc(ip)) then;  k = m + 1
      else

      Do i = 1,ntrm
       coef(i,ip) = coef(i,ip) + ctrm(i)
      End do

      Return
      end if; end if

      go to 1
    2 Continue 

! ... new coefficient:

      ncoef = ncoef + 1  
      coef(:,ncoef)=ctrm(:); intc(ncoef)=int; idfc(ncoef)=idf 
      
      if(k.eq.ncoef) then
       ipcoef(k)=ncoef
      else
       Do m=ncoef,k+1,-1;  ipcoef(m)=ipcoef(m-1);  End do
       ipcoef(k)=ncoef
      end if        

! ... it is time for relocation:

      if(ncoef.eq.mcoef) Call Alloc_coef (mcoef+icoef)
 
      END Subroutine Add_coef


