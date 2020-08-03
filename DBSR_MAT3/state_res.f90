!======================================================================
      Subroutine State_res
!======================================================================
!     extracts the data from JNT_BNK for specific case
!----------------------------------------------------------------------
      Use dbsr_mat
      Use c_data, only: ntype,kpol1,kpol2

      Use det_list; Use def_list
      Use new_dets; Use new_defs 

      Implicit none
      Real(8) :: C,CC,CCC, t1,t2,t3
      Integer :: i,j,ibuf,nnbuf,int,jcase,k,i1,i2,i3,i4,it,jt,idf,  &
                 ic,jc,is,js,is1,js1,is2,js2,ik,jk,ich,jch,ip1,ip2, &
                 j1,j2,j3,j4, itype,m,kpol, io,jo, nc1,nc2

      Integer, parameter :: ib2 = 2**2, ib5 = 2**5, ib10= 2**10
       
!----------------------------------------------------------------------
! ... read coef.s from jnt_bnk:                                            

      Call CPU_TIME(t1)

      if(myid.eq.0) then
       rewind(nub) 
       Call Read_symc(nub)
       Call Read_symt(nub)
       Call Read_done_out(nub)
       Call Read_det (nub)
       Call Read_def (nub)
      end if

      Call br_dets 

      Call CPU_TIME(t2)

      if(pri.gt.0.and.icase.eq.0.and.idiag.eq.0.and.debug.gt.0) then
       write(pri,'(/a,i7,i9,a )') &
       'ndet,kdet  = ',ndet,kdet,'  - number of overlap determinants'  
       write(pri,'( a,i7,i9,a/)') &
       'ndef,kdef  = ',ndef,kdef,'  - number of overlap factors'  
       C = 4.d0*(2*ndet + kdet + 2*ndef + kdef)/(1024*1024)
       write(pri,'(a,T40,f10.2,a/)') 'memory of dets:', C,' Mb'
       write(pri,'(a,T40,f10.2,a/)') 'broadcast of dets:',(t2-t1)/60,' min'
      end if

!----------------------------------------------------------------------
! ... first, fill the buffer:

      nnbuf = 0
   10 Continue

      t_check = 0.d0
      t_add = 0.d0
      t_det = 0.d0

      Call CPU_TIME(t1)

! ... first, fill the buffer with angular coefficients:

      if(myid.eq.0) then
      i=1
    1 read(nub,end=2) Cbuf(i),itb(i),jtb(i),intb(i),idfb(i)
      i=i+1; if(i.le.mcbuf) go to 1
    2 ncbuf=i-1
      end if

! ... broadcast the buffer:

      Call br_barrier
      Call br_buffer

      nnbuf=nnbuf+1
!----------------------------------------------------------------------
! ... processing the data in buffer:

      nc1 = 0; nc2 = 0
      Do ibuf=1,ncbuf

! ... decode the integral:

      int=intb(ibuf)

      jcase = mod(int,ib2)
      if(icase.le.2.and.icase.ne.jcase) Cycle
      if(icase.eq.3.and.jcase.ne.2) Cycle

      int = int/ib2; kpol = mod(int,ib10); if(icase.le.1) kpol=0 
      if(kpol.lt.kpol1) Cycle
      if(kpol.gt.kpol2) Cycle

      int = int/ib10;  i4  = mod(int,ib5)
      int = int/ib5;   i3  = mod(int,ib5)
      int = int/ib5;   i2  = mod(int,ib5)
      i1  = int/ib5

! ... determine the range of states for given coeff. 

      it=itb(ibuf); is1 = IT_state1(it); if(is1.eq.0) Cycle
      jt=jtb(ibuf); js1 = IT_state1(jt); if(js1.eq.0) Cycle
      is2 = IT_state2(it); js2 = IT_state2(jt)

      idf = idfb(ibuf)

!----------------------------------------------------------------------
! ... loop over all relevant states: 

! ... loop over all relevant states: 

      Do ik=is1,is2;  is =IS_order(ik)
                      ip1=IP_state(is)
                      no1=no_state(is)
                      np1(1:no1)=IP_orb(ip1+1:ip1+no1)
                      ich=ich_state(is)

      Do jk=js1,js2;  js= IS_order(jk)
                      ip2=IP_state(js)
                      no2=no_state(js)
                      np2(1:no2)=IP_orb(ip2+1:ip2+no2)
                      jch=ich_state(js)

! ... consider only low-half of interaction matrix:

       if(it.eq.jt.and.is.lt.js) Cycle

! ... check if we need only diagonal/nondiagnal blocks:

       m = 1; if(ich.ne.jch.or.ich.gt.nch.or.jch.gt.nch) m=0 
       if(idiag.gt.0.and.m.eq.0) Cycle
       if(idiag.lt.0.and.m.eq.1) Cycle

! ... if this case belongs to given process:

       if(ich.le.nch.and.jch.le.nch) then
        m = icc(ich,jch)
       elseif(ich.le.nch.and.jch.gt.nch) then
        m = icb(ich,jch-nch)
       elseif(ich.gt.nch.and.jch.le.nch) then
        m = icb(jch,ich-nch)
       elseif(ich.gt.nch.and.jch.gt.nch) then
        m = ibb(ich-nch,jch-nch)
       else
        Stop 'mycase: channel index out of range'
       end if

       if(m.eq.0) Cycle

! ... coefficient:

       C=CBUF(ibuf);  if(ich.eq.jch.and.is.ne.js) C = C + C

       CC = C*WC(is)*WC(js);  if(abs(CC).lt.Eps_C) Cycle

! ... define integral for specific orbitals:

       j1=IP_orb(i1+ip1); j2=IP_orb(i2+ip1)
       j3=IP_orb(i3+ip2); j4=IP_orb(i4+ip2)

! ... we do not need anymore the configuration index 
! ... except pertuber (N+1)-electron configurations   

       i=0; if(ich.gt.nch) i=ich-nch
       j=0; if(jch.gt.nch) j=jch-nch
       ic = max(i,j); jc = min(i,j)

! ... find overlap factors with extracted continuum:  

       Call Det_fact_new (idf,np1,np2,ipbs);  if(nndef.eq.0) Cycle

! ... send the final coefficients to archive:

       nc1 = nc1 + 1
       Do i = 1,nndef
        CCC = CC * Adef(i); io=iof(i); jo=jof(i)
        Call Add_integral (kpol,j1,j2,j3,j4,CCC,ic,jc,io,jo)
        nc2 = nc2 + 1
       End do

      End do    ! over js
      End do    ! over is

      End do    !  over buffer

      Call CPU_TIME(t2)

      Do kpol = kpol1,kpol2
       Do itype = 1,ntype
        Call Add_matrix(itype,kpol)
       End do
      End do

      Call CPU_TIME(t3)

      if(pri.gt.0.and.myid.gt.0)  write(pri,'(a,i5,3i10,2f9.2,a)') &
       'bufer:',nnbuf, ncbuf,nc1,nc2,(t2-t1)/60,(t3-t2)/60,' min'
      if(myid.eq.0.and.debug.gt.0) &
      write(pri,'(a,i5,T40,f10.2,a)') 'bufer:',nnbuf, (t3-t1)/60,' min'

! ... check the data bank again:

      if(ncbuf.eq.mcbuf) go to 10  

      End Subroutine State_res

!======================================================================
      Subroutine Read_done_out(nu)
!======================================================================
!     dummy reading the IT_done array from unit "nu"
!---------------------------------------------------------------------- 
      Use symt_list 

      Implicit none
      Integer :: nu,i,i1,i2,n
      Integer(1) :: j

      read(nu) n 
      i1=1; i2=mrecl; if(i2.gt.n) i2=n
      Do 
       read(nu) (j,i=i1,i2)
       i1=i1+mrecl; if(i1.gt.n) Exit
       i2=i2+mrecl; if(i2.gt.n) i2=n
      End do 
 
      End Subroutine Read_done_out
