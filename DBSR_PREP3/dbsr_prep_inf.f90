!=====================================================================
      Subroutine dbsr_prep_inf
!=====================================================================
!     print information for dbsr_prep program
!---------------------------------------------------------------------
      Implicit none
      Character ::  name = ' '
      
      Call Read_name(name)
      if(name.ne.'?') Return
      write(*,'(a)') &
'                                                                            ', &
'  dbsr_prep analizes the input target states and prepare them for DBSR:     ', &
'                                                                            ', &
'  INPUT FILES:                                                              ', &
'                                                                            ', &
'     dbsr_par  -  parameters of calculation if any                          ', &
'     target_jj -  list of target states and partial waves                   ', &
'     c- and bsw-files for each target state and perturber if any            ', &
'     knot.dat  -  B-spline parameters                                       ', &
'     target_sub.bsw  - additional orbitals for orthogonality, if any        ', &
'                                                                            ', &
'  OPTIONAL ARGUMENTS:                                                       ', &
'                                                                            ', &
'     eps_ovl  [1.d-6] - tolerance for overlaps                              ', &
'     eps_core [1.d-5] - tolerance for orthogonality to the core orbitals    ', &
'     eps_phys [0.25 ] - minimum occupation number for orbital to be physical', &
'     eps_sub  [0.5  ] - tolerance for substitution orbitals                 ', &
'                                                                            ', &
'  additional arguments JJ_min, JJ_max  can be used at nlsp=0 for generation ', &
'  of possible partial waves                                                 ', &
'                                                                            '   
      Stop 
                                                                             
      End Subroutine dbsr_prep_inf                                             

