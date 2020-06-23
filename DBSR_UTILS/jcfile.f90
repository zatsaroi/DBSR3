!======================================================================
!     utility       j c f i l e
!
!                   C O P Y R I G H T -- 2007
!
!     Written by:   Oleg Zatsarinny
!======================================================================
!     Extract the state expansion from DBSR_ci j-file
!
!     arguments: 1. name for j-file
!                2. # of seeking solution 
!                3. name for result c- and bsw- files
!                4. eps_c - tolerance for coefficients
!
!     When eps_c > 0, configurations in result c-file are ordered
!     according their weights
!
!     Call as:  jcfile   5p4  1   5p5_1  0.0000000001 
!----------------------------------------------------------------------
      Use conf_jj, ip => IS_order

      Implicit none
  
      Character(80) :: AF,BF,CF,name,out

      Integer :: nuc =1       !   name.c         
      Integer :: nuj =2       !   name.m or name.bm
      Integer :: iout=3       !   result.c

      Real(8) :: eps_c, E
      Integer :: i,j,k,m,ic,ic1,ic2,nsol,isol,iarg,  SYSTEM
      Integer, external :: Ifind_position, Icheck_file

!----------------------------------------------------------------------
!                                                           input data:
      iarg = command_argument_count() 

      if(iarg.lt.2) then
       write(*,*)
       write(*,*) 'jcfile extracts the state expansion from DBSR_ci j-file'
       write(*,*)
       write(*,*) 'You should provide as position arguments the following data:' 
       write(*,*)
       write(*,*) '1. name for j-file'
       write(*,*) '2. # of solution'
       write(*,*) '3. name for resulting c-file'
       write(*,*) '4. eps_c - tolerance for coefficients'
       write(*,*)  
       write(*,*) 'When eps_c > 0, configurations in result c-file are ordered'
       write(*,*) 'according their weights'
       write(*,*)  
       Stop 
      else
       Call GET_COMMAND_ARGUMENT(1,name)
       i = len_trim(name)
       if(name(i-1:i).eq.'.c') name(i-1:i)='  '
       if(name(i-1:i).eq.'.j') name(i-1:i)='  '

       Call GET_COMMAND_ARGUMENT(2,CF); read(CF,*) isol
       
       write(out,'(a,a,i3.3)') trim(name),'_',isol
       if(iarg.gt.2) Call GETARG(3,out)
       i = len_trim(out)
       if(out(i-1:i).eq.'.c') out(i-1:i)='  '

       eps_c = 1.d-14
       if(iarg.gt.3) then
        Call GET_COMMAND_ARGUMENT(4,CF); read(CF,*) eps_c
       end if

      end if

!----------------------------------------------------------------------
! ... read list of configuration from c-file:

      AF = trim(name)//'.j'
      CF = trim(name)//'.c'
      Call Check_file(CF);  Open(nuc,file=CF)

      Call Read_conf_jj(nuc,0,'add','nocheck')

!----------------------------------------------------------------------
! ... read solution from j-file:

      Call Check_file(AF); Open(nuj,file=AF)

      Call Read_ipar(nuj,'ncfg',ic     ) 
      if(ic.ne.ncfg) Stop 'jc: ncfg in j-file <> ncfg in c-file'
      Call Read_ipar(nuj,'nsol',nsol) 
      if(isol.gt.nsol) Stop 'jc: # of solution > nsol'

      i=Ifind_position(nuj,'Solutions');  read(nuj,*) 

      Do i=1,isol
       read(nuj,*) j
       read(nuj,*) E,j,ic1,ic2
       read(nuj,*) WC(ic1:ic2)
      End do

      Close(nuj)

!----------------------------------------------------------------------
! ... order the configurations according their weights:
 
      Allocate(IP(ncfg)); Do i=1,ncfg; IP(i)=i; End do
 
      if(eps_c.gt.0.d0) then
       Do i=ic1,ic2-1
        m=i
        Do j=i+1,ic2
         if(abs(WC(IP(j))).gt.abs(WC(IP(m)))) m=j
        End do
        if(m.ne.i) then; k=IP(m); IP(m)=IP(i); IP(i)=k; end if
       End do
      end if

!----------------------------------------------------------------------
! ... output the c-file:


      BF = trim(out)//'.c'
      open(iout,file=BF)

      rewind(nuc)
      read(nuc,'(a)') AS
      write(iout,'(a,f16.8)') 'Core subshells:',E

      Do
       read(nuc,'(a)') AS;   write(iout,'(a)') TRIM(AS)
       if(AS(1:3).eq.'CSF') Exit
      End do

      Do ic=ic1,ic2; i=IP(ic); if(abs(WC(i)).lt.Eps_c) Exit
       Call print_conf_jj (iout,i,WC(i))
      End do

      write(iout,'(a)') '*'

      AF = trim(name)//'.bsw'
      if(Icheck_file(AF).gt.0) then
       BF = trim(out)//'.bsw'
       write(CF,'(a,1x,a,1x,a)') 'cp',trim(AF),trim(BF)
       I =  SYSTEM(CF)
      end if

      AF = trim(name)//'.w'
      if(Icheck_file(AF).gt.0) then
       BF = trim(out)//'.w'
       write(CF,'(a,1x,a,1x,a)') 'cp',trim(AF),trim(BF)
       I =  SYSTEM(CF)
      end if


      End  ! program jc


