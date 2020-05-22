!======================================================================
      Subroutine Read_arg
!======================================================================
! ... read parameters for given case
!----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none

      Open(pri,file=AF_log)
      Open(nup,file=AF_par,action='READ')

! ... read parameters or arguments:

      Call Read_rpar(nup,'c_comp',c_comp)
      Call Read_rarg(    'c_comp',c_comp)
      Call Read_ipar(nup,'max_ll',max_ll)
      Call Read_iarg(    'max_ll',max_ll)
      Call Read_ipar(nup,'min_ll',min_ll)
      Call Read_iarg(    'min_ll',min_ll)
      Call Read_ipar(nup,'max_ka',max_ka)
      Call Read_iarg(    'max_ka',max_ka)
      Call Read_ipar(nup,'min_ka',min_ka)
      Call Read_iarg(    'min_ka',min_ka)
      Call Read_ipar(nup,'kort'  ,kort  )
      Call Read_iarg(    'kort'  ,kort  )
      Call Read_ipar(nup,'debug' ,debug )
      Call Read_iarg(    'debug' ,debug )
      Call Read_ipar(nup,'iread_targ',iread_targ)
      Call Read_iarg(    'iread_targ',iread_targ)

      write(pri,'(a/)')  'DBSR_CONF parameters:'

      write(pri,'(a,f5.3,a/)')  'c_comp = ',C_comp,' - tolerance for compansation configurations'
      if(max_ll.ne.-1) &
      write(pri,'(a,i5,a/)')    'max_ll = ',max_ll,' - restiction on max. small l'
      if(min_ll.ne.-1) &
      write(pri,'(a,i5,a/)')    'min_ll = ',min_ll,' - restiction on min. small l'
      if(max_ka.ne. 0) &
      write(pri,'(a,i5,a/)')    'max_ka = ',max_ka,' - restiction on max. kappa'
      if(min_ka.ne. 0) &
      write(pri,'(a,i5,a/)')    'min_ka = ',min_ka,' - restiction on min. kappa'

      End Subroutine Read_arg
