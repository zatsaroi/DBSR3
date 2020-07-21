!=======================================================================
      Subroutine D_OUT
!=======================================================================
!     define and output dipole matrix for given initial state
!-----------------------------------------------------------------------
      Use dbsr_dmat
      Use target_jj, only: etarg
      Use conf_jj

      Implicit none
      Character(64) :: Label1
      Integer :: i,j,nhm,kch,kpert,ns1,isol
      Real(8) :: S,SL,SV,E1
      Real(8), allocatable :: CL(:),CV(:)
      Real(8), allocatable :: eval(:)
 
      if(ktype.ne.'E') Stop 'D_out: non-electric-transition case ? '
      if(kpol.ne.1) Stop 'D_out: kpol <> 1 --> non-dipole case ? '
!----------------------------------------------------------------------
!                                      define the initial bound states:
      if(ctype1.eq.'b') then

       i=INDEX(BF_b,'.'); AF=BF_b(1:i)//ALS1
       Call Check_file(AF); Open(nub1,file=AF,form='UNFORMATTED')
       rewind(nub1)

       read(nub1) nhm,kch,kpert,ns1,jot1,parity,nstate1

       if(ns1 .ne.ns ) Stop 'dbsr_dmat: ns1 <> ns '
       if(nch1.ne.kch) Stop 'dbsr_dmat: nch1 --> ?'
       if(npert1.ne.kpert) Stop 'dbsr_dmat: npert1 --> ?'
       if(kdm1.ne.nhm) Stop 'dbsr_dmat: nhm1 --> ?'

       if(Allocated(C1)) Deallocate(C1); Allocate(C1(kdm1))
       Do i=1,istate1
        read(nub1) j,Label1
        read(nub1) E1
        read(nub1) (C1(j),j=1,kdm1)
       End do

      elseif(ctype1.eq.'c') then
	   
       C1=1.d0                       ! already in d-matrix
       rewind(nuc1)
       read(nuc1,'(15x,e16.8)') E1
       label1 = name1

      else

       Stop 'dbsr_dmat: invalid initial state in d_out' 

      end if

!-----------------------------------------------------------------------
! ... read final states: inner region solutions

      i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALS2
      Call Check_file(AF)
      Open(nur,file=AF,form='UNFORMATTED',STATUS='OLD')
      read(nur) nhm,nstate2,kch,kpert
      if(nhm.ne.kdm2)     Stop 'D_out, rsol: kdm2 --> ?'
      if(kch.ne.nch2)     Stop 'D_out, rsol: kch2 --> ?'
      if(kpert.ne.npert2) Stop 'D_out, rsol: npert2 --> ?'
      Allocate(eval(nstate2)) 
      read(nur) eval

!----------------------------------------------------------------------
! ... calculation and output the dipole matrix:

      Allocate(CL(nstate2),CV(nstate2))
      if(allocated(C2)) Deallocate(C2); Allocate(C2(kdm2))

      !... loop over final set:

      Do isol=1,nstate2
       read(nur) C2
       SL=0.d0; SV=0.d0
       Do j=1,kdm2
        S = SUM(DV(1:kdm1,j)*C1(1:kdm1));  SV = SV + C2(j)*S
        S = SUM(DL(1:kdm1,j)*C1(1:kdm1));  SL = SL + C2(j)*S
       End do
       CL(isol)=SL; CV(isol)=SV*c_au
      End do

      i=INDEX(AF_d,'.');  AF=AF_d(1:i)//ALS2
      Open(nud,file=AF,form='UNFORMATTED')

      write(nud) jot2,0,parity2,Etarg(1),nstate2                 
      write(nud) jot1,0,parity1,E1,       Label1                          

      write(*,*) 'jot2,0,parity2,Etarg(1),nstate2',jot2,0,parity2,Etarg(1),nstate2                 

      write(*,*) 'jot1,0,parity1,E1,Label1',jot1,0,parity1,E1, Label1                          

      write(nud) (cl(i),cv(i),i=nstate2,1,-1)        

      Deallocate(CL,CV,eval)

      End Subroutine D_out
