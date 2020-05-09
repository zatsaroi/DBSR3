!======================================================================
      Subroutine Conf_calc
!======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------

      USE param_jj
      USE conf_jj 

      USE nljm_orbitals

      USE term_exp

      Use symc_list; Use symt_list

      USE coef_list; USE zoef_list; Use boef_list

      Implicit none 

      Integer :: i,j,k,l,m,k1,k2,is,js, it,jt, ij
      Integer, external :: DEF_ij
      Real(8) :: t1,t2
      Real(8), External :: RRTC

      t1=RRTC()

! ... get the job:

    1 Call Get_det_exp  

      if(ic.le.0) Return

! ... initial preparations:

      Call Get_symc(ic,Jtotal1,no1,nn1,kn1,ln1,jn1,iq1,in1)        
      k = 1
      Do i=1,no1
       Do j=k,k+iq1(i)-1
        nnsym1(j)=i; Lsym1(j)=ln1(i); Jsym1(j)=jn1(i)
       End do
       k=k+iq1(i)
      End do
              
      Call Get_symc(jc,Jtotal2,no2,nn2,kn2,ln2,jn2,iq2,in2)        
      k=1
      Do i=1,no2
       Do j=k,k+iq2(i)-1
        nnsym2(j)=i; Lsym2(j)=ln2(i); Jsym2(j)=jn2(i)
       End do
       k=k+iq2(i)
      End do

! ... define number of needed terms:
       
       if(Allocated(IP_kt12)) Deallocate(IP_kt12)
       Allocate(IP_kt12(kt1,kt2));  IP_kt12 = 0 

       k = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle
        ij = DEF_ij(it,jt)
        if(IT_done(ij).ne.0) Cycle
        k=k+1; IP_kt12(k1,k2) = 1 
       End do; End do 

       if(k.eq.0) Stop 'Conf_loop: k = 0'
  
! ...  repeat to get term pointers:
       
       Call Alloc_trm(k)
       k = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(IP_kt12(k1,k2).eq.0) Cycle
        k=k+1; itc(k)=it;  jtc(k)=jt
       End do; End do 

! ...  initial allocations:

       ncoef=0; Call Alloc_coef(icoef)
       nboef=0; Call Alloc_boef(iboef)
       nblk=0;  Call Alloc_blk(iblk); ncblk(0)=0

! ...  calculations:

       Do kd1 = 1,kdt1

        Do i=1,ne;  Msym1(i)=mj_orb(IP_det1(i,kd1));  End do

       Do kd2 = 1,kdt2

        Do i=1,ne;  Msym2(i)=mj_orb(IP_det2(i,kd2));  End do

        Call Det_me

        if(nzoef.gt.0) Call Term_loop 

       End do;  End do

! ... send the results:

      Call Send_res

      t2=RRTC()
      if(pri.gt.0) then
      write(pri,'(a,3i8,f10.2,a)')  'send conf.', ic,jc,ic_case, (t2-t1)/60, '  min'
      Close(pri); Open(pri,file=AF_pri,position='APPEND')
      end if

      go to 1

      End Subroutine Conf_calc


!======================================================================
      Subroutine get_det_exp
!======================================================================

      USE MPI

      USE param_jj, only: myid, ierr, pri,AF_pri
      USE conf_jj,  only: ne
      USE term_exp

      Implicit none

      Integer :: status(MPI_STATUS_SIZE)

      Call MPI_RECV(ic,1,MPI_INTEGER, 0,0, MPI_COMM_WORLD, status, ierr)

      if(ic.le.0) Return

      Call MPI_RECV(kt1 ,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(kdt1,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(ip_kt1)) Deallocate(ip_kt1); Allocate(ip_kt1(kt1))
      Call MPI_RECV(ip_kt1,kt1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(ip_det1)) Deallocate(ip_det1); Allocate(ip_det1(ne,kdt1))
      Call MPI_RECV(ip_det1,ne*kdt1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(C_det1)) Deallocate(C_det1); Allocate(C_det1(kt1,kdt1))
      Call MPI_RECV(C_det1,kt1*kdt1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(jc  ,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(kt2 ,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(kdt2,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(ip_kt2)) Deallocate(ip_kt2); Allocate(ip_kt2(kt2))
      Call MPI_RECV(ip_kt2,kt2,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(ip_det2)) Deallocate(ip_det2); Allocate(ip_det2(ne,kdt2))
      Call MPI_RECV(ip_det2,ne*kdt2,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(C_det2)) Deallocate(C_det2); Allocate(C_det2(kt2,kdt2))
      Call MPI_RECV(C_det2,kt2*kdt2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, status,ierr)


      End Subroutine get_det_exp


!======================================================================
      Subroutine send_res
!======================================================================

      USE MPI

      USE param_jj, only: myid, ierr, pri, AF_pri

      USE term_exp
      Use NDET_list
      Use NDEF_list 

      Use coef_list, only: mcoef,ntrm,ncoef,idfc,intc,coef

      Implicit none

      Call MPI_SEND(ncoef,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

! ... term indexes:

      Call MPI_SEND(kt1,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kt2,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(IP_kt1,kt1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IP_kt2,kt2,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(IP_kt12,kt1*kt2,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      if(ncoef.eq.0) Return

! ... coef.s:

      Call MPI_SEND(ntrm,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(idfc,ncoef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(intc,ncoef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(coef,ntrm*ncoef,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)


! ... dets:

      Call MPI_SEND(ndet,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      if(ndet.eq.0) Return

      Call MPI_SEND(ldet,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(KPD,ndet,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IPD,ndet,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(NPD,ldet,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

! ... defs:

      Call MPI_SEND(ndef,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ldef,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(KPf,ndef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IPf,ndef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(NPf,ldef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      End Subroutine send_res

