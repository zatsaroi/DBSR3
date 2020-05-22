!=====================================================================
      Subroutine dbsr_conf_inf
!=====================================================================
!     print information for bsr_prep program
!---------------------------------------------------------------------
      Implicit none
      Character ::  name = ' '
      
      Call Read_name(name)
      if(name.ne.'?') Return
      write(*,'(a)') &
'                                                                           ', &
'dbsr_conf prepares configuration expansions for all partial waves (cfg.nnn)', &
'and assigns orthogonality conditions if any                                ', &
'                                                                           ', &
'all input files are created by dbsr_prep program;                          ', &
'except file "target" may be modified with additional list of partial waves ', &
'                                                                           ', &
'additional file "target_del" can be used for deleting specific channels    ', &
'                                                                           ', &
'OPTIONAL ARGUMENTS:  (-1 value means out of play)                          ', &
'                                                                           ', &
'c_comp  [1.01]  -  tolerance for compensation configuration                ', &
'                   to put orthogonal conditions                            ', &
'max_ll  [-1]    - restiction on max. small l                               ', &
'min_ll  [-1]    - restiction on min. small l                               ', &
'max_ka  [-1]    - restiction on max. kappa for channel orbitals            ', &
'max_ka  [-1]    - restiction on min. kappa for channel orbitals            ', &
'kort    [-1]    - if > 0, the orth.conditions from cfg-files are in play   ', &
'                                                                           '
      Stop 
                                                                             
      End Subroutine dbsr_conf_inf                                             
