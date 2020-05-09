!======================================================================
      Module  dbsr_breit
!======================================================================
!     main parameters and input/outpu files used in the program
!----------------------------------------------------------------------
      Use zconst
      Implicit none

      Real(8) :: eps_c = 1.d-7           ! tolerence for coefficients
      Integer :: mk = 9                  ! maximum multipole index
      Integer :: klsp1=0, klsp2=0        ! range of partial waves
      Integer :: mbreit = 0              ! flag for Breit contribution
!----------------------------------------------------------------------
! ... the names and units for input/output files:
!     AF  -  standard (default) names
!     BF  -  names with indication of partial wave number

      Character(40) :: AF,BF

! ... runing information:

      Integer :: pri=7;  Character(40) :: AF_pri = 'dbsr_breit.log'

! ... c-file:

      Integer :: nuc=1;  Character(40) :: AF_c = 'rcsl.inp'
                         Character(40) :: BF_c = 'cfg.nnn'
! ... data bank:

      Integer :: nub=2;  Character(40) :: AF_b = 'int_bnk'
                         Character(40) :: BF_b = 'int_bnk.nnn'

! ... new results if any:

      Integer :: nur=3;  Character(40) :: AF_r = 'int_res'
                         Character(40) :: BF_r = 'int_res.nnn'

      Logical :: new      ! pointer on the previous calculation
      Logical :: icalc    ! pointer for need of new calculations

! ... scratch files:

      Integer :: nui=11   ! intermediate results
      Integer :: nud=12   ! for det. expansions
      Integer :: nua=13   ! for accumulation of data

      End Module  dbsr_breit


!======================================================================
      Subroutine Open_jj(nu,klsp)
!======================================================================
!     open (closed) files in the breit_bsr program!
!----------------------------------------------------------------------
      Use dbsr_breit

      Implicit none
      Integer, intent(in) :: nu,klsp
      Integer :: i,iarg
      Character(3) :: ALSP
      Logical :: EX

      write(ALSP,'(i3.3)') klsp

      Select case(nu)

      Case(1)         ! c-file

       AF = AF_c
       if(klsp.gt.0) then
        i=Index(BF_c,'.'); AF=BF_c(1:i)//ALSP; BF_c=AF
       else
        iarg = command_argument_count()
        if(iarg.gt.0) then
         Call get_command_argument(1,BF); i=INDEX(BF,'=')
         if(i.eq.0) then
          i=INDEX(BF,'.')-1; if(i.lt.0) i=LEN_TRIM(BF)
          AF=BF(1:i)//'.c'; AF_b = BF(1:i)//'.bnk'
         end if
        end if
       end if

       Call Check_file(AF)
       Open(nu,file=AF,status='OLD')

      Case(2)               ! bnk-file

       AF = AF_b
       if(klsp.gt.0) then
        i=Index(BF_b,'.'); AF=BF_b(1:i)//ALSP; BF_b=AF
       end if
       Inquire(file=AF,exist=EX)
       new=.TRUE.
       if(EX) then
        new = .FALSE.
        Open(nub,file=AF,form='UNFORMATTED',STATUS='OLD')
       end if

      Case(3)

       AF = AF_r
       if(klsp.gt.0) then
        i=Index(BF_r,'.'); AF=BF_r(1:i)//ALSP; BF_r=AF
       end if
       Open(nur,file=AF,form='UNFORMATTED')

      Case(7)

       AF = AF_pri;  Open(nu,file=AF)

      Case(11,12,13)

       Open(nu,form='UNFORMATTED',status='SCRATCH')

      Case default

       Stop ' open_jj: unit value is out of list '

      End select

      End Subroutine Open_jj



