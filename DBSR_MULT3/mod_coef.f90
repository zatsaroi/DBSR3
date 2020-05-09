!====================================================================
      Module coef_list
!====================================================================
!     Contains a set of coefficients with two identifiers (intc,idfc)
!     The list of ordered according the pointer 'ipcoef'
!--------------------------------------------------------------------
      Implicit none 
    
      Integer :: ncoef = 0                    ! number of coefficients
      Integer :: mcoef = 0  
      Integer :: icoef = 2**15
      Integer :: ntrm  = 0                    ! number of terms

      Real(8), allocatable :: coef(:,:)       ! coefficients (1:ntrm,1:mcoef)
      Integer, allocatable :: intc(:),idfc(:) ! their attributes (1:mcoef)
      Integer, allocatable :: ipcoef(:)       ! ordering pointer (1:mcoef)

      REAL(8), allocatable :: ctrm(:),cdtrm(:)! new coefficients (1:ntrm)
      Integer, allocatable :: itc(:),jtc(:)   ! term pointers    (1:ntrm)

      Integer :: int, idf                    ! current integral

      Integer :: nus = 99

      End Module coef_list


!======================================================================
      Subroutine Alloc_trm(nt)
!======================================================================
      USE coef_list

      Implicit none
      Integer, intent(in) :: nt

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
      Integer, intent(in) :: mc
      Integer :: i

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
       Open(nus,status='scratch',form='UNFORMATTED') 
       rewind(nus) 
       Do i = 1,ncoef
        write(nus) coef(1:ntrm,i),intc(i),idfc(i),ipcoef(i)
       End do
       Deallocate(coef,intc,idfc,ipcoef)
       Allocate(coef(ntrm,mc),intc(mc),idfc(mc),ipcoef(mc))
       rewind(nus) 
       Do i = 1,mcoef
        read(nus) coef(1:ntrm,i),intc(i),idfc(i),ipcoef(i)
       End do
       mcoef=mc
       write(*,*) ' realloc_coef: mcoef = ', mcoef
      end if
	        
      End Subroutine Alloc_coef


!======================================================================
      Subroutine Add_coef
!======================================================================
!     add new coefficient to the list 
!----------------------------------------------------------------------
      Use coef_list

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

      Do i=1,ntrm; coef(i,ip)=coef(i,ip)+ctrm(i);  End do

      Return
      end if; end if

      go to 1
    2 Continue 

! ... new coefficient:

      ncoef = ncoef + 1  
      coef(1:ntrm,ncoef)=ctrm(1:ntrm); intc(ncoef)=int; idfc(ncoef)=idf 
      
      if(k.eq.ncoef) then
       ipcoef(k)=ncoef
      else
       Do m=ncoef,k+1,-1;  ipcoef(m)=ipcoef(m-1);  End do
       ipcoef(k)=ncoef
      end if        

! ... it is time for relocation:

      if(ncoef.eq.mcoef) Call Alloc_coef (mcoef+icoef)
 
      End Subroutine Add_coef


