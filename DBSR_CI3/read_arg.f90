!======================================================================
      Subroutine Read_arg
!======================================================================
!     read input arguments from command line
!----------------------------------------------------------------------
      Use dbsr_ci      

      Implicit none 

      Call Read_name(name)
      if(name.eq.'?'.or.len_trim(name).eq.0) Call Inf_dbsr_ci
      Call Check_BSR_name(name)

      Call Read_iarg('mbreit',mbreit)
      Call Read_iarg('msol'  ,msol  )
      Call Read_iarg('nzero' ,nzero )
      Call Read_iarg('debug' ,debug )
      Call Read_iarg('mdiag' ,mdiag )

      Call Read_rarg('eps_det',eps_det)
      Call Read_rarg('eps_ovl',eps_ovl)
      Call Read_rarg('eps_o'  ,eps_o  )
      Call Read_rarg('Emax'   ,Emax   )

      Call Read_iarg('check_c',check_c)

! ... c_data paramters:

      Call Read_iarg('mb',mblock)
      Call Read_iarg('nb',nblock)
      Call Read_iarg('kb',kblock)
      Call Read_iarg('mk',mpol  )

      Call Read_rarg('eps_c',eps_c)
      
      End Subroutine Read_arg

!======================================================================
      Subroutine Check_BSR_name(AF)
!======================================================================
!     rename old-version file-name as name.c just to 'name'
!----------------------------------------------------------------------
      Character(*) :: AF
      i = LEN_TRIM(AF)
      if(i.le.2) Return
      if(AF(i-1:i).eq.'.c') AF(i-1:i) = '  '
      End Subroutine Check_BSR_name

