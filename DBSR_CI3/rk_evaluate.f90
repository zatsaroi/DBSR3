!======================================================================
      Subroutine  Rk_evaluate(k,itype)
!======================================================================
!     processing of Rk(Sk)-integrals from the module 'c_data'
!     for given multipole index "k" and "itype",
!     e.i. requires only one rkb array in module DBS_integrals 
!----------------------------------------------------------------------
      Use c_data;         Use DBS_grid
      Use dbsr_ci;        Use DBS_orbitals_pq
     
      Implicit none
      Integer, intent(in) :: k,itype
      Integer :: i,j, i1,i2, j1,j2, ic,jc
      Real(8) :: S,SS
      Real(8), external :: rk_pppp,rk_qqqq,rk_pqpq,rk_qpqp, &
                           sk_ppqq_pq,sk_pqqp_pq

      if(ncdata.eq.0) Return

      Do j=1,ncdata;  i=IPT(j); if(abs(cdata(i)).lt.Eps_c) Cycle

       i1 = K1(i)/ibi; i2 = mod(K1(i),ibi)
       j1 = K2(i)/ibi; j2 = mod(K2(i),ibi)

       Select case(itype)
        Case(1);             S = rk_pppp(i1,j1,i2,j2,k)
        Case(2);             S = rk_qqqq(i1,j1,i2,j2,k)
        Case(3);             S = rk_pqpq(i1,j1,i2,j2,k)
        Case(4);             S = rk_qpqp(i1,j1,i2,j2,k)
        Case(5);             S = sk_ppqq_pq(i1,j1,i2,j2,k)
        Case(6);             S = sk_pqqp_pq(i1,j1,i2,j2,k)
        Case default;  Stop 'Rk_data: unknown type of intgrals'      
       End select

       SS = Cdata(i)*S

       ic=k3(i); jc=k4(i)

       if(debug.gt.0.and.ic.ne.jc) &
        Call prj_int(pri,2,k,i1,j1,i2,j2,S,Cdata(i),ic,jc,itype)

       if(jc.le.NZERO) then;        HM(ic,jc) = HM(ic,jc) + SS
       elseif(ic.eq.jc) then;       DM(jc) = DM(jc) + SS
       end if

      End do

      End Subroutine  Rk_evaluate


