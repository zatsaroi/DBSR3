!======================================================================
!     utility       o r d e r _ j j
!
!                   C O P Y R I G H T -- 2008
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!     Ordering the configurations according their weights
!
!     arguments: 1. name for c-file, results in cc-file
!                2. eps_c - tolerance for coefficients
!----------------------------------------------------------------------
      Use conf_jj, ip => IS_order

      Implicit none

      Character(80) :: AF,BF

      Integer :: nuc =1       !   name.c         
      Integer :: iout=2       !   name.cc

      Real(8) :: eps_c = 0.d0
      Integer :: i,j,k,m,ic,ic1,ic2,iarg,ib,ii

!----------------------------------------------------------------------
!                                                           input data:
      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)
 
      if(iarg.lt.1.or.AF.eq.'?') then
        write(*,*) 
        write(*,*) "ORDER_JJ: ordering the ASF's according to their weights"
        write(*,*) 
        write(*,*) 'Call as:  order_jj  name.c  eps_c'
        write(*,*) 
        write(*,*) 'eps_c - optional cut-off parameter'
        write(*,*) 
        write(*,*) 'OUTPUT:  name.cc '
        Stop ' '
       end if
      
      if(iarg.gt.1) then; Call GET_COMMAND_ARGUMENT(2,BF); read(BF,*) eps_c; end if

!----------------------------------------------------------------------
! ... read list of configuration from c-file:

      i=INDEX(AF,'.',BACK=.TRUE.)     
      if(AF(i+1:i+1).ne.'c') Stop ' c-file should have extension .c'
      Call Check_file(AF)
      Open(nuc,file=AF)

      Call Read_conf_jj(nuc,0,'add','nocheck')

      Call Def_jblocks

!----------------------------------------------------------------------
! ... order the configurations according their weights:
 
      Allocate(IP(ncfg)); Do i=1,ncfg; IP(i)=i; End do

      Do ib = 1,njbl; ic1 = JTc1(ib); ic2 = JTc2(ib) 
     
       Do i=ic1,ic2-1
        m=i
        Do j=i+1,ic2
         if(abs(WC(IP(j))).gt.abs(WC(IP(m)))) m=j
        End do
        if(m.ne.i) then; k=IP(m); IP(m)=IP(i); IP(i)=k; end if
       End do
    
      End do

!----------------------------------------------------------------------
! ... output the c-file:

      i = LEN_TRIM(AF); BF = AF(1:i)//'c';  open(iout,file=BF)

      rewind(nuc)
      read(nuc,'(a)') AS;  write(iout,'(a)') TRIM(AS)
      read(nuc,'(a)') AS
      i = LEN_TRIM(AS); i=(i/5+1)*5
      write(iout,'(a)') AS(1:i)
      read(nuc,'(a)') AS;  write(iout,'(a)') TRIM(AS)

      Do
       read(nuc,'(a)') AS
       if(AS(6:6).eq.'(') Exit
       i = LEN_TRIM(AS); i=(i/5+1)*5
       write(iout,'(a)') AS(1:i)
      End do

      Do ib = 1,njbl; ic1 = JTc1(ib); ic2 = JTc2(ib) 

       Do ic=ic1,ic2; i=IP(ic); if(abs(WC(i)).lt.Eps_c) Exit
        Call Print_conf_jj(iout,i,WC(i))
       End do

       if(ib.lt.njbl)  write(iout,'(a)') ' *'
       if(ib.eq.njbl)  write(iout,'(a)') '* '

      End do

      End  ! program order_jj



