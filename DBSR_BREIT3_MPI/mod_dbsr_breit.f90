!======================================================================
      Module  dbsr_breit
!======================================================================
!     main parameters and input/outpu files used in the program
!----------------------------------------------------------------------
      Use zconst
      Implicit none

      Character(80) :: name

      Real(8) :: eps_c = 1.d-7            ! tolerance for coefficients
      Integer :: mk = 9                   ! maximum multipole index
      Integer :: klsp1=0, klsp2=0         ! range of partial waves
      Integer :: mbreit = 0               ! flag for Breit contribution
      Integer :: debug = 0

      Logical :: new      ! pointer on the previous calculation
      Logical :: icalc    ! pointer for need of new calculations

!----------------------------------------------------------------------
! ... the names and units for input/output files:
!
!     AF  -  standard (default) names
!     BF  -  names with indication of partial wave number

      Integer, parameter :: ma=80
      Character(ma) :: AF,BF

! ... c-file:

      Integer :: nuc=1;  Character(ma) :: AF_c = 'cfg.inp'           ! name.c
                         Character(ma) :: BF_c = 'cfg.nnn'
! ... data bank:

      Integer :: nub=2;  Character(ma) :: AF_b = 'int_bnk'           ! name.bnk'
                         Character(ma) :: BF_b = 'int_bnk.nnn'

! ... runing information:

      Integer :: pri=3;  Character(ma) :: AF_pri = 'dbsr_breit.log'  ! name.dbsr_breit_log
                         Character(ma) :: BF_pri = 'dbsr_breit.nnn'
! ... new results if any:

      Integer :: nur=10; Character(ma) :: AF_r = 'int_res'           ! name.int_res 
                         Character(ma) :: BF_r = 'int_res.nnn'

! ... scratch files:

      Integer :: nui=11;  Character(40) :: AF_i = 'int_new'          ! name.int_new
                          Character(ma) :: BF_i = 'int_new.nnn'

      Integer :: nud=12   ! for det. expansions
      Integer :: nua=13   ! for accumulation of data

! ... MPI:

      Integer :: nprocs, myid, ierr
      Integer, allocatable :: ip_proc(:)

      End Module  dbsr_breit


!======================================================================
      Subroutine Open_jj(nu,klsp)
!======================================================================
!     open (closed) files in the dbsr_breit program
!----------------------------------------------------------------------
      Use dbsr_breit

      Implicit none
      Integer, intent(in) :: nu,klsp
      Integer, external :: Icheck_file

      Select case(nu)

      Case(1)         ! c-file

       if(klsp.gt.0) then
        write(AF,'(a,i3.3)') 'cfg.',klsp 
       else
        AF = AF_c; if(len_trim(name).gt.0) AF = trim(name)//'.c'       
       end if
       Call Check_file(AF);  BF_c = AF
       Open(nu,file=AF,action='READ')

      Case(2)         ! bnk-file

       if(klsp.gt.0) then
        write(AF,'(a,i3.3)') 'int_bnk.',klsp 
       else
        AF = AF_b
        if(len_trim(name).gt.0.and.Icheck_file(AF).eq.0) AF = trim(name)//'.bnk'
       end if
       BF_b = AF
       new=.TRUE.
       if(Icheck_file(AF).ne.0) then
        new = .FALSE.
        Open(nu,file=AF,form='UNFORMATTED',action='READ')
       end if

      Case(3)

       if(klsp.gt.0) then
        write(AF,'(a,i3.3)') 'dbsr_breit.',klsp 
       else
        AF = AF_pri
        if(len_trim(name).gt.0) AF = trim(name)//'.dbsr_breir_log'
       end if
       Open(nu,file=AF);  BF_pri=AF

      Case(10)

       if(klsp.gt.0) then
        write(AF,'(a,i3.3)') 'int_res.',klsp 
       else
        AF = AF_r
        if(len_trim(name).gt.0) AF = trim(name)//'.int_res'
       end if
       Open(nu,file=AF,form='UNFORMATTED',action='WRITE'); BF_r=AF

      Case(11)

       if(klsp.gt.0) then
        write(AF,'(a,i3.3)') 'int_new.',klsp 
       else
        AF = AF_r
        if(len_trim(name).gt.0) AF = trim(name)//'.int_new'
       end if
       Open(nu,file=AF,form='UNFORMATTED',action='WRITE'); BF_i=AF

      Case(12,13)

       Open(nu,form='UNFORMATTED',status='SCRATCH')

      Case default

       Stop ' open_jj: unit value is out of list '

      End select

      End Subroutine Open_jj



