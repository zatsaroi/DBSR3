!======================================================================
      Subroutine Gen_matrix(ic1,ic2)
!======================================================================
!     generation of interaction and overlap matrices
!     according to JNT_BNK
!----------------------------------------------------------------------
      Use dbsr_ci      
      Use conf_jj
      Use c_data
      Use det_list
      Use def_list
      Use DBS_orbitals_pq, only: nbs,ebs

      Implicit real(8) (A-H,O-Z)
      Real(8) :: CBUF(maxnc)
      Integer :: itb(maxnc),jtb(maxnc),intb(maxnc),idfb(maxnc)
      Real(8), external :: Lval

! ... prepare arrays in c_data:

      Call Alloc_C_data(1)

! ... read overlap factors:

      rewind(nub) 
      Call Read_symc(nub)
      Call Read_symt(nub)
      Call Read_done(nub)
      Call Read_det(nub)
      Call Read_def(nub)

      write(pri,'(/a,2i8 )') 'ndet  = ',ndet,kdet
      write(pri,'( a,2i8/)') 'ndef  = ',ndef,kdef

      HM = 0.d0; SM = 0.d0; DM = 0.d0
!----------------------------------------------------------------------
!                                                  processing the data:
    1 ncbuf = 0
      Do i = 1,maxnc
       read(nub,end=2) Cbuf(i),itb(i),jtb(i),intb(i),idfb(i)
       ncbuf = ncbuf + 1
      End do
    2 Continue

      Do ibuf=1,ncbuf; int = intb(ibuf)

       Call Decode_int (met,kpol,I1,I2,I3,I4,int)

       if(kpol.gt.mk) Cycle

! ... determine the range of states for given coeff. 

      C = CBUF(ibuf)

      it = itb(ibuf);  jt = jtb(ibuf)
      is1 = IT_state1(it); js1 = IT_state1(jt)
      if(is1.eq.0.or.js1.eq.0) Cycle
      is2 = IT_state2(it); js2 = IT_state2(jt)

      idf = idfb(ibuf)
      if(idf.gt.0) then
       nd=KPF(idf); ip=IPF(idf); NP(1:nd)=NPF(ip+1:ip+nd)
      end if

!----------------------------------------------------------------------
! ... cycle over states:

      Do ik=is1,is2; ic=IS_order(ik); ip1=IP_state(ic)

      Do jk=js1,js2; jc=IS_order(jk); ip2=IP_state(jc)

! ... consider only low-half of interaction matrix

       if(it.eq.jt.and.ic.lt.jc) Cycle

       if(ic.lt.ic1.or.ic.gt.ic2) Cycle
       if(jc.lt.ic1.or.jc.gt.ic2) Cycle

       if(mdiag.ne.0.and.ic.ne.jc) Cycle

       i = max0(ic,jc)-ic1+1
       j = min0(ic,jc)-ic1+1

       if(i.gt.NZERO.and.j.gt.NZERO.and.i.ne.j) Cycle

!----------------------------------------------------------------------
! ...  find overlap factor:

       DF = 1.d0
       if(idf.gt.0) then
        Do ii=1,nd
         id=NP(ii)/ibf; kd=KPD(id); ip=IPD(id)
         Do jj=1,kd
          k=NPD(jj+ip)
          jp1=k/ibd +    ip1; NP1(jj)=IP_orb(jp1)
          jp2=mod(k,ibd)+ip2; NP2(jj)=IP_orb(jp2)
         End do
         DF = DF * VDET(kd,NP1,NP2)**mod(NP(ii),ibf)
         if(abs(DF).lt.Eps_det) Exit
        End do
        if(abs(DF).lt.Eps_det) Cycle
       end if 

       CC = C * DF; if(abs(CC).lt.eps_C) Cycle

!----------------------------------------------------------------------
! ... find integral:

       j1=IP_orb(i1+ip1); j2=IP_orb(i2+ip1)
       j3=IP_orb(i3+ip2); j4=IP_orb(i4+ip2)

       if(ncorr.gt.0) CC = CC * S_corr(met,kpol,j1,j2,j3,j4)

       if(met.eq.0) then

        if(j.le.NZERO) SM(i,j) = SM(i,j) + CC

       elseif(met.eq.1) then

        CC = CC*Lval(j1,j3)
        if(j.le.NZERO) then
          HM(i,j) = HM(i,j) + CC
        elseif(i.eq.j) then
          DM(j) = DM(j) + CC
        end if

       else

        Call Add_integral (kpol,i,j,j1,j2,j3,j4,CC,mbreit)

       end if

      End do   ! over jc
      End do   ! over ic

      End do   ! over icoef

      if(ncbuf.eq.maxnc) go to 1  ! go for new data from data bank
!-----------------------------------------------------------------------
! ... final generation of interaction matrix:

      Do jtype = 1,ntype
       Do jpol = kpol1,kpol2
        Call Add_matrix(jtype,jpol)
       End do
      End do

!----------------------------------------------------------------------
! ... print non-trivial overlap matrix elements:

      AF_ovl = trim(name)//'.ovl'
      Open(nuo,file=AF_ovl)

!      write(nuo,'(/a,f10.8/)') &
!	   'Non-trivial total overlaps: > eps_o = ',eps_o
      Do i=1,NZERO
       Do j=1,i-1
        if(abs(SM(i,j)).lt.eps_o) Cycle
        write(nuo,'(2I8,2x,F13.5)') i,j,SM(i,j)
       End do
      End do
      
      C = 5*eps_ovl
      write(pri,'(/a,f10.8/)') &
	   'Non-normalized states: > 5*eps_ovl = ',C
      Do i=1,NZERO
       if(abs(SM(i,i)-1.D0).gt.C) &
       write(pri,'(F13.5,2x,I8)') SM(i,i), i
      End do

!----------------------------------------------------------------------
! ... find NORT parameter --> indication on the generelized eigenvalue
! ... problem

      NORT = -1
      Do j=1,NZERO;  Do i=j,NZERO
        C = SM(i,j); if(i.eq.j) C = C - 1.d0
        if(abs(C).lt.5*Eps_ovl) Cycle
        NORT = 0;  Exit
      End do; End do

      if(NORT.eq.0) write(pri,'(/a/)') &
        'NORT = 0 --> generelized eigenvalue problem'
      if(NORT.eq.-1) write(pri,'(/a/)') &
        'NORT = 0 --> simple eigenvalue problem'

      End Subroutine GEN_matrix



!======================================================================
      Subroutine Decode_int (met,k,I1,I2,I3,I4,int)
!======================================================================
!     decode the integral form "int_bnk"
!----------------------------------------------------------------------
      Implicit none
      Integer :: int, met,k,I1,I2,I3,I4, ii
      Integer,parameter :: ib2 =2**2
      Integer,parameter :: ib5 =2**5
      Integer,parameter :: ib10=2**10

      ii = int
      met = mod(ii,ib2);  ii = ii/ib2
      k   = mod(ii,ib10); ii = ii/ib10
      I4  = mod(ii,ib5);  ii = ii/ib5
      I3  = mod(ii,ib5);  ii = ii/ib5
      I2  = mod(ii,ib5);  I1 = ii/ib5

      End Subroutine Decode_INT


