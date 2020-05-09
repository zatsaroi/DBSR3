!======================================================================
      Subroutine dv_out
!======================================================================
! ... out the dipole vector for some solution (istate1)
!----------------------------------------------------------------------
      Use dbsr_dmat

      Implicit none
      Integer :: i,j, j1, ic1,ic2, jc1,jc2,  &
                 nhm,kch,kpert,ns1,nc,parity
      Real(8) :: E1
      Character(64) :: Label1
      Integer, external :: Ifind_position

      if(allocated(C1)) Deallocate(C1); Allocate(C1(kdm1))
      C1=1.d0; ic1=1; ic2=kdm1
!----------------------------------------------------------------------
! ... check the initial state:

      if(ctype1.eq.'b') then

       i=INDEX(BF_b,'.'); AF=BF_b(1:i)//ALS1
       Call Check_file(AF); Open(nub1,file=AF,form='UNFORMATTED')
       rewind(nub1)

       read(nub1) nhm,kch,kpert,ns1,jot1,parity,nstate1

       if(ns1 .ne.ns )       Stop 'dbsr_dmat: ns1 <> ns '
       if(nch1.ne.kch)       Stop 'dbsr_dmat: nch1 --> ?'
       if(npert1.ne.kpert)   Stop 'dbsr_dmat: npert1 --> ?'
       if(kdm1.ne.nhm)       Stop 'dbsr_dmat: nhm1 --> ?'
       if(parity1.ne.parity) Stop 'dbsr_dmat: parity1 --> ?'

       if(nstate1.lt.istate1) Stop 'dbsr_dmat: nstate1 < istate1 '

       Do j1=1,istate1
        read(nub1) j,label1
        read(nub1) E1
        read(nub1) C1
       End do

      elseif(ctype1.eq.'j') then

       AF = name1(1:iname1)//'j'
       Call Check_file(AF);  Open(nub1,file=AF)
       Call Read_ipar(nub1,'ncfg',nc     ) 
       if(nc.ne.kdm1) Stop 'dbsr_dmat: nc in j-file <> ncfg1'
       Call Read_ipar(nub1,'nsol',nstate1) 

       if(nstate1.lt.istate1) Stop 'dbsr_dmat: nstate1 < istate1 '

       i=Ifind_position(nub1,'Solutions');  read(nub1,*)  
       Do j1=1,nstate1
        read(nub1,*) j,label1
        read(nub1,*) E1,jot1,ic1,ic2
        read(nub1,*) C1(ic1:ic2)
       End do

      elseif(ctype1.eq.'c') then;   
      
       Call Jdef_JPE(nuc1,jot1,parity1,E1) 

      else
    
       Stop 'dv_out: ctype1 --> ?'

      end if

!----------------------------------------------------------------------
! ... get d-vector:

      i=INDEX(AF_dv,'.');  AF=AF_dv(1:i)//ALS2
      Open(nuv,file=AF,form='UNFORMATTED')

      write(nuv) kdm2,nch2,npert2,E1,jot1,kpol,ktype                  

      if(allocated(C2)) Deallocate(C2); Allocate(C2(kdm2))
      C2=0.d0; jc1=1; jc2=kdm2

      Do j=jc1,jc2;  C2(j)=SUM(DL(ic1:ic2,j)*C1(ic1:ic2));  End do
      write(nuv) C2

      if(ktype.ne.'M') then
       Do j=jc1,jc2;  C2(j)=SUM(DV(ic1:ic2,j)*C1(ic1:ic2));  End do
       C2 = C2 *c_au**2 
       write(nuv) C2
      end if

      Deallocate(C1,C2)

      End Subroutine dv_out
