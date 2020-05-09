!======================================================================
      Integer Function Idet_fact (i1,i2,i3,i4)
!======================================================================
!     determines the overlap factor and its position in NDEF list
!     for matrix element between two determinant wave functions
!     located in module 'nljm_orbitals' 
!
!     (j1,j2) and (j3,j4) - active electrons in the first and second 
!                           determinants 
!
!     Calls:  Nadd_det, Nadd_def, ISORT
!----------------------------------------------------------------------
      Use zconst, only: ibd,ibf
      USE nljm_orbitals

      Implicit none
      Integer, intent(in) :: i1,i2,i3,i4
      Integer :: i, k,k1,k2, is,id,kd
      Integer, external :: Iadd_ndet, Iadd_ndef, ISORT

      kd = 0
      Do is = 1,NSYM

       k1 = 0
       Do i = IPsym1(is-1)+1,IPsym1(is)
        if(i.eq.i1.or.i.eq.i2) Cycle
        k1 = k1 + 1;  N1(k1) = nnsym1(Isym1(i))
       End do

       k2 = 0
       Do i = IPsym2(is-1)+1,IPsym2(is)
        if(i.eq.i3.or.i.eq.i4) Cycle
        k2 = k2 + 1;  N2(k2) = nnsym2(Isym2(i))
       End do

       if(k1.ne.k2) Stop 'Idet_fact: k1 <> k2'

       if(k1.eq.0) Cycle

       NP(1:k1) = N1(1:k1)*ibd + N2(1:k1);  id=Iadd_ndet(k1,NP)

       Do k=1,kd
        if(N3(k).ne.id) Cycle; N4(k)=N4(k)+1; id=0; Exit
       End do
       if(id.eq.0) Cycle

       kd=kd+1; N3(kd)=id; N4(kd)=1

      End do

      Idet_fact=0; if(kd.eq.0) Return

      NP(1:kd)=N3(1:kd)*ibf + N4(1:kd); k=ISORT(kd,NP)

      Idet_fact = Iadd_ndef(kd,NP)

      End Function Idet_fact
