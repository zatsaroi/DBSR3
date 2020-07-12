!======================================================================
      Subroutine Read_arg
!======================================================================
!     read arguments from command line 
!----------------------------------------------------------------------
      Use dbsr_breit

      Implicit none
      Integer :: i, klsp = 0

! ... read arguments in command line:

      Call Read_iarg('klsp'  ,klsp  )
      Call Read_iarg('klsp1' ,klsp1 )
      Call Read_iarg('klsp2' ,klsp2 )
      Call Read_iarg('mk'    ,mk    )
      Call Read_rarg('eps_c' ,eps_c )
      Call Read_iarg('mbreit',mbreit)

      if(klsp.gt.0) klsp1=klsp
      if(klsp2.lt.klsp1) klsp2=klsp1

      name = ' ';  Call Read_name(name)

      i = LEN_TRIM(name)
      if(i.gt.2)  then
       if(name(i-1:i).eq.'.c') name(i-1:i)='  '       
      end if

      End Subroutine Read_arg
