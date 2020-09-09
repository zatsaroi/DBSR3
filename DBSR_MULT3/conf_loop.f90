!======================================================================
      Subroutine Conf_loop
!======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------
      Use dbsr_mult;       Use symc_list 
      Use conf_jj;         Use symt_list
      Use nljm_orbitals;   Use coef_list  
      Use term_exp;        Use zoef_list

      Implicit none 
      Integer :: i,j,k,k1,k2,is,js, it,jt
      Real(8) :: t1,t2,tt, CN   
      Real(8), external :: Z_3j2
      Integer, external :: DEF_ij

!----------------------------------------------------------------------
!                                          cycle 1 over configurations:
      rewind(nud)
      Do is=1,ic_case

       ! ... read determinant expansion:
       Read(nud) ic,kt1,kdt1
       if(allocated(IP_kt1))  Deallocate(IP_kt1)
       Allocate (IP_kt1(kt1)); Read(nud)  IP_kt1
       if(allocated(IP_det1)) Deallocate(IP_det1)
       Allocate (IP_det1(ne,kdt1)); Read(nud)  IP_det1
       if(allocated(C_det1)) Deallocate(C_det1)
       Allocate (C_det1(kt1,kdt1)); Read(nud)  C_det1

       ! ... restore the configuration
       Call Get_symc(ic,Jtotal1,no1,nn1,kn1,ln1,jn1,iq1,in1)        

       ! ... get the one-electron symmetries:
       k = 1
       Do i=1,no1
        Do j=k,k+iq1(i)-1
         nnsym1(j)=i; Lsym1(j)=ln1(i); Jsym1(j)=jn1(i)
        End do
        k=k+iq1(i)
       End do

       Call CPU_time(t1)

!----------------------------------------------------------------------
!                                          cycle 2 over configurations:
      rewind(nud)
      Do js=1,is

       ! ... read determinant expansion:
       Read(nud) jc,kt2,kdt2
       if(allocated(IP_kt2)) Deallocate(IP_kt2)
       Allocate (IP_kt2(kt2)); Read(nud)  IP_kt2
       if(allocated(IP_det2)) Deallocate(IP_det2)
       Allocate(IP_det2(ne,kdt2)); Read(nud)  IP_det2
       if(allocated(C_det2)) Deallocate(C_det2)
       Allocate(C_det2(kt2,kdt2)); Read(nud)  C_det2

! ... check the need:       

       if(JC_need(DEF_ij(ic,jc)).eq.0) Cycle      

       ! ... restore the configuration
       Call Get_symc(jc,Jtotal2,no2,nn2,kn2,ln2,jn2,iq2,in2)        

       ! ... get the one-electron symmetries:
       k=1
       Do i=1,no2
        Do j=k,k+iq2(i)-1
         nnsym2(j)=i; Lsym2(j)=ln2(i); Jsym2(j)=jn2(i)
        End do
        k=k+iq2(i)
       End do

! ... define number of needed terms:
       
       if(Allocated(IP_kt12)) Deallocate(IP_kt12)
       Allocate (IP_kt12(kt1,kt2)); IP_kt12=0 
       k = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle
        if(IT_done(DEF_ij(it,jt)).ne.0) Cycle
        k=k+1; IP_kt12(k1,k2) = 1 
       End do; End do 
       if(k.eq.0) Stop 'Conf_loop: k = 0'
  
! ... repeat to get term pointers:
       
       Call Alloc_trm(k)
       k = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(IP_kt12(k1,k2).eq.0) Cycle
        k=k+1; itc(k)=it;  jtc(k)=jt
       End do; End do 

! ... initial allocations:

       ncoef=0; Call Alloc_coef(icoef)

! ... define the normalization constant:

       qpol = Jtotal1-Jtotal2
       CN = Z_3j2(Jtotal1,-Jtotal1,kkpol,qpol,Jtotal2,Jtotal2)      

       if(CN.eq.0.d0) Cycle
       CN = 1.d0/CN  
!----------------------------------------------------------------------
! ... calculations:

       Do kd1=1,kdt1; Do i=1,ne; Msym1(i)=mj_orb(IP_det1(i,kd1)); End do

       Do kd2=1,kdt2; Do i=1,ne; Msym2(i)=mj_orb(IP_det2(i,kd2)); End do

        Call Det_mult;   if(nzoef.gt.0) Call Term_loop 

       End do; End do

!----------------------------------------------------------------------
! ... store results for given config.s:

       Call Add_res(CN)

      End do    ! over jc

      Call CPU_time(t2);  tt=t2-t1
      Call Incode_confj1
      write(*,  '(a,4i6,f5.1,a,a)')  &
        'ic =',ic,nsymc,kt1,kdt1,tt,' sec ',CONFIG(1:ia)
      write(pri,'(a,4i6,f5.1,a,a)')  &
        'ic =',ic,nsymc,kt1,kdt1,tt,' sec ',CONFIG(1:ia)

      End do    ! over ic

      End Subroutine Conf_loop


