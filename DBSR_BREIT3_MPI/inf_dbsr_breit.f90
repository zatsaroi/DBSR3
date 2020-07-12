!======================================================================
      Subroutine inf_dbsr_breit
!======================================================================
!     provide screen information about program 'dbsr_breit'
!----------------------------------------------------------------------
      Character :: name

      Call Read_name(name)
      if(name.ne.'?') Return

      write(*,'(a)')                                                       &
'                                                                        ',&
'     DBSR_BREIT - angular coefficients for the Dirac-Coulomb Hamiltonian',&
'     It is a part of the DBSR project.                                  ',&
'------------------------------------------------------------------------',&
'     INPUT ARGUMENTS:                                                   ',&
'                                                                        ',&
'     c-file -  name.c [cfg.inp by default]                              ',&
'               DBSR option:  cfg.001, cfg.002, ..., cfg.nnn             ',&
'               if range of partial waves are given through              ',&
'               arguments klsp1, klsp2                                   ',&
'               (or just klsp=.. for one partial wave)                   ',&
'                                                                        ',&
'     mk     -  max.multipole index [by default mk=9]                    ',&
'     eps_c  -  tollerence for the angular coefficients [1.d-7]          ',&
'     mbreit -  = 0 -> only Coulomb interaction (by default)             ',&
'               = 1 -> include Breit interaction                         ',&
'----------------------------------------------------------------------- ',&
'     examples:   1.  dbsr_breit 3p5.c  (or just 3p5)                    ',&
'                 2.  dbsr_breit klsp1=1 klsp2=5                         ',&
'                 3.  dbsr_breit km=5 mbreit=1                           ',&
'----------------------------------------------------------------------- ',&
'     INPUT FILES:   c-files with list of atomic states in GRASP format: ',&
'                    name.c (or cfg.inp, or cfg.nnn)                     ',&
'                                                                        ',&
'     OUTPUT FILES:  data bank for angular coefficients:                 ',&
'                    name.bnk (or int_bnk, or int_bnk.nnn)               ',&
'----------------------------------------------------------------------  ',&
'                                                                        ',&
'                                                                        '
      Stop                          
                                    
      End Subroutine inf_dbsr_breit       
                                     
