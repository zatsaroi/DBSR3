!======================================================================
      Subroutine D_matr
!======================================================================
!     This routine computes the matrix of dipole matrix elements
!     in length and velocity forms according to MULT_BNK (unit nub)
!----------------------------------------------------------------------
      Use dbsr_dmat
      Use conf_jj
      Use cmdata,          only: itype,ntype
      Use DBS_orbitals_pq, only: ipbs
      Use def_list
      Use new_defs

      Implicit none
      Real(8) :: C,CC,CL,CV,SL,SV, sign
      Integer :: i, ip1,ip2, is,js,ic,jc,ih,jh, jt1,jt2, &
                 it,jt,ik,jk,is1,is2,js1,js2, i1,i2, j1,j2,   &
                 int,icase,idf,ich1,ich2, io,jo, kk1,kk2,kk3, &
                 ii, m
      Integer, external :: Idef_dtype, no_ic_jj, jot_ic, Ifind_pert_jj

! ... allocate dipole matrices:

      Allocate(DL(kdm1,kdm2)); DL=0.d0
      if(ktype.eq.'E') then; Allocate(DV(kdm1,kdm2)); DV=0.d0; end if
      Call allocate_cmdata(1)

! ... read determinant factors:

      Call Read_det(nub)
      Call Read_def(nub)

!----------------------------------------------------------------------
!                                                  processing the data:
    1 read(nub,end=2) C,it,jt,int,idf

      is1 = IT_state1(it); js1 = IT_state1(jt)
      is2 = IT_state2(it); js2 = IT_state2(jt)

      if(is1.eq.0.or.js1.eq.0) go to 1   ! no relevant states

      Call Decode_mult(icase,i1,i2,int)

!----------------------------------------------------------------------
!                                                    cycle over states:

      Do ik=is1,is2;  is =IS_order(ik);  ip1=IP_state(is)
       no1=no_ic_jj(is); np1(1:no1)=IP_orb(ip1+1:ip1+no1)
       jt1 = jot_ic(is)

      Do jk=js1,js2;  js= IS_order(jk);  ip2=IP_state(js)
       no2=no_ic_jj(js); np2(1:no2)=IP_orb(ip2+1:ip2+no2)
       jt2 = jot_ic(js)

       if(is.gt.ncfg1.and.js.gt.ncfg1) Cycle
       if(is.le.ncfg1.and.js.le.ncfg1) Cycle

       if(is.le.ncfg1) then
        m=1; ih=is; jh=js-ncfg1
        j1=IP_orb(ip1+i1);  ich1=ipbs(j1)          ! change sign ?
        j2=IP_orb(ip2+i2);  ich2=ipbs(j2) 
        sign = 1.d0
       else
        m=2; ih=js; jh=is-ncfg1 
        j1=IP_orb(ip2+i2);  ich1=ipbs(j1)
        j2=IP_orb(ip1+i1);  ich2=ipbs(j2) 
        sign = (-1)**((jt1-jt2)/2)  
       end if

! ... include the expansion coefficients: 

       CC = C * C1(ih)*C2(jh);  
    
       if(abs(CC).lt.eps_c) Cycle

! ... find index of the pertuber if any: 

       ic = Ifind_pert_jj(ilsp1,ih)
       jc = Ifind_pert_jj(ilsp2,jh)                

! ... find overlap factors with extracted continuum:  

       Call Det_fact(idf,np1,np2); if(nndef.eq.0) Cycle

! ... multipole integrals if any:

       CL=CC; CV=CC
       if(ich1+ich2.eq.0) then
        Call dip(j1,j2,SL,SV)
        CL=CC*SL
        if(ktype.eq.'E') CV=CC*SV
       end if

! ... send the final coefficients to archive:

       Do i = 1,nndef
        SL = CL * Adef(i)
        SV = CV * Adef(i)
        io=iof(i); jo=jof(i)
        if(jo.gt.0) then
         if(ipbs(io/ibo).eq.0) then
          ii = io; io=jo; jo=ii
         end if
        end if

        itype=Idef_dtype(j1,j2,ich1,ich2,io,jo,ic,jc,kk1,kk2,kk3)

        Call Add_coef(SL,SV,kk1,kk2,kk3)

       End do

      End do   !  over js
      End do   !  over is

      go to 1   !  new data from data bank
    2 Continue

! ... final generation of interaction matrix:

      Do itype=1,ntype; Call Gen_matrix;  End do

      End Subroutine D_matr

