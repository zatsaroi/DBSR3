!======================================================================
      Subroutine Read_arg 
!======================================================================
!     read arguments from command line and check default settings
!----------------------------------------------------------------------
      Use dbsr_mult 

      Implicit none
      Integer :: iarg, kpol_old 
      Character(1) :: ktype_old 
      Logical :: EX

! ... open log-file:

      Open(pri,file=AF_log)
      write(pri,'(a/)') &
       'GENERATION OF DATA BANK FOR MULTIPOLE MATRIX ELEMENTS:  '

! ... read arguments from command line:

      iarg = command_argument_count()
      if(iarg.lt.3) then
       write(*,*) 'Should be at least three arguments from:'
       write(*,*)
       write(*,*) 'dbsr_mult name1 name2 E1|M1|.. [AF_bnk]'
       Stop
      end if

      Call get_command_argument(1,AF1);  Call Check_file(AF1)
      Open(nu1,file=AF1)
      Call get_command_argument(2,AF2);  Call Check_file(AF2)
      Open(nu2,file=AF2)
      
! ... define the type of calculations:

      Call get_command_argument(3,AF)
      read(AF,'(a1,i1)') ktype,kpol; kkpol=kpol+kpol
      write(pri,'(a,a1,i1)') 'transition -> ',ktype,kpol

      if(ktype.eq.'M'.and.kpol.eq.0) Stop ' kpol=0 for M-type ? '

! ... New - ?

      AF_bnk(10:11) = AF(1:2) 

      if(iarg.ge.4)  Call get_command_argument(4,AF_bnk)
      Inquire (file=AF_bnk, exist=EX)

      if(EX) then
       Open(nub,file=AF_bnk,form='UNFORMATTED',status='OLD')
       read(nub) ktype_old,kpol_old
       if(ktype_old.ne.ktype) Stop 'in-apropriate MULT_BNK !'
       if(kpol_old.ne.kpol) Stop 'in-apropriate MULT_BNK !'
       write(pri,'(/a)')  'calculation:  continued'
       new = .FALSE.
      else
       write(pri,'(/a)')  'calculation:  new'
       new = .TRUE.
      end if

      Open(nui,form='UNFORMATTED',status='SCRATCH')
      Open(nua,form='UNFORMATTED',status='SCRATCH')
      Open(nud,form='UNFORMATTED',status='SCRATCH')
      Open(nur,file=AF_res,form='UNFORMATTED')
 
      End Subroutine Read_arg


