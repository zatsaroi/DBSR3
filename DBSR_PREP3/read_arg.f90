!======================================================================
      Subroutine Read_arg
!======================================================================
!     read input parameters
!----------------------------------------------------------------------
      Use dbsr_prep

      Implicit none 

      Call Read_rpar(nup,'eps_ovl' ,eps_ovl )
      Call Read_rpar(nup,'eps_core',eps_core)
      Call Read_rpar(nup,'eps_sub' ,eps_sub )
      Call Read_rpar(nup,'eps_phys',eps_phys)
      Call Read_rpar(nup,'eps_targ',eps_targ)
      Call Read_rpar(nup,'eps_conf',eps_conf)

      Call Read_ipar(nup,'JJ_min'  ,JJ_min  )
      Call Read_ipar(nup,'JJ_max'  ,JJ_max  )

      Call Read_rarg('eps_ovl' ,eps_ovl )
      Call Read_rarg('eps_core',eps_core)
      Call Read_rarg('eps_sub' ,eps_sub )
      Call Read_rarg('eps_phys',eps_phys)
      Call Read_rarg('eps_targ',eps_targ)
      Call Read_rarg('eps_conf',eps_conf)

      Call Read_iarg('JJ_min'  ,JJ_min  )
      Call Read_iarg('JJ_max'  ,JJ_max  )

      write(pri,'(a/)') 'DBSR_PREP:'
      write(pri,'(a,d12.3,a/)') 'eps_ovl  = ',eps_ovl, &
       '   -  tolerance for orthogonality between orbitals'
      write(pri,'(a,d12.3,a/)') 'eps_core = ',eps_core,&
       '   -  tolerance for orthogonality to core orbitals'
      write(pri,'(a,d12.3,a/)') 'eps_sub  = ',eps_sub, &
       '   -  tolerance for substitution orbitals'
      write(pri,'(a,d12.3,a/)') 'eps_phys = ',eps_phys,&
       '   -  tolerance for physical orbitals'
      if(JJ_min.ge.0) &
      write(pri,'(a,i12,a/)')    'JJ_min   = ',JJ_min,&
       '   -  minimum 2J-value for the partial waves'
      if(JJ_max.ge.0) &
      write(pri,'(a,i12,a/)')    'JJ_max   = ',JJ_max,&
       '   -  maximum 2J-value for the partial waves'

      End Subroutine Read_arg


