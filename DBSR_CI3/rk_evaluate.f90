!======================================================================
      Subroutine  Rk_evaluate(k,itype)
!======================================================================
!     processing of Rk(Sk)-integrals from the module 'c_data'
!     for given multipole index "k" and "itype",
!     e.i. requires only one rkb array in module DBS_integrals 
!----------------------------------------------------------------------
      Use dbsr_ci;       Use DBS_grid
      Use c_data;        Use DBS_orbitals_pq
     
      Implicit none
      Integer, intent(in) :: k,itype
      Integer :: i,j, ii,jj, i1,i2, j1,j2, ip1,ip2, jp1,jp2, ic,jc
      Real(8) :: S,SS, t1,t2
      Real(8) :: deni(ns,ks),denj(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      if(ncdata.eq.0) Return

      Call CPU_time(t1)

      Select case(itype)
       Case(1);  Call mrk_pppp(k); ip1=1; ip2=1; jp1=1; jp2=1  
       Case(2);  Call mrk_qqqq(k); ip1=2; ip2=2; jp1=2; jp2=2  
       Case(3);  Call mrk_pqpq(k); ip1=1; ip2=1; jp1=2; jp2=2  
       Case(4);  Call mrk_qpqp(k); ip1=2; ip2=2; jp1=1; jp2=1  
       Case default;  Stop 'Rk_evaluate: unknown type of intgrals'      
      End Select

      ii=0; jj=0; S = 0.d0 

      Do j=1,ncdata;  i=IPT(j); if(abs(cdata(i)).lt.Eps_c) Cycle

       if(K1(i).ne.ii) then
        i1 = K1(i)/ibi; i2 = mod(K1(i),ibi)
        Call Density(ns,ks,deni,pq(1,ip1,i1),pq(1,ip2,i2),'s')
        Call Convol(ns,ks,conv,deni,2,'s','s')
       end if

       if(K2(i).ne.jj) then
        j1 = K2(i)/ibi; j2 = mod(K2(i),ibi)
        Call Density(ns,ks,denj,pq(1,jp1,j1),pq(1,jp2,j2),'s')
       end if

       if(K1(i).ne.ii.or.K2(i).ne.jj)  S = SUM_AmB(ns,ks,conv,denj,'s')

       ii=K1(i);  jj=K2(i);   SS=S*cdata(i);  ic=k3(i); jc=k4(i)

       if(debug.gt.0.and.ic.ne.jc) &
        Call prj_int(pri,2,k,i1,j1,i2,j2,S,Cdata(i),ic,jc,itype)

       if (jc.le.NZERO) then;       HM(ic,jc) = HM(ic,jc) + SS
       elseif(ic.eq.jc) then;       DM(jc) = DM(jc) + SS
       end if

      End do

      if(debug.gt.0) then
       Call CPU_time(t2)
       write(pri,'(/a,2i3,i10,T30,f8.2,a/)') 'k, itype, nc:', k,itype,ncdata,(t2-t1)/60,'  min'
      end if

      End Subroutine  Rk_evaluate


