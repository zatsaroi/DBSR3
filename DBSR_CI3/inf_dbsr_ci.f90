!======================================================================
      Subroutine inf_dbsr_ci
!======================================================================
!     provide screen information about program 'dbsr_ci'
!----------------------------------------------------------------------
      write(*,'(a)')                                                      &
'                                                                       ',&
'     dcsr_ci - a configuration-interaction program in the  Dirac-      ',&
'     Coulomb-Breit approach including case of non-orthogonal orbitals. ',&
'     It is a part of the DBSR project, which uses B-spline for repri-  ',&
'     santation of radial orbitals.                                     ',&
'-----------------------------------------------------------------------',&
'     INPUT:  command-line arguments                                    ',&
'             name.inp_ci - input parameters (optional)                 ',&
'             name.c      - list of configurations in the GRASP format  ',&
'             name.bsw    - radial functions in B-spline representation ',&
'             name.bnk    - data bank for angular coefficients          ',&
'                           (after program dbsr_breit)                  ',&
'     OUTPUT: name.log_ci - running information                         ',&
'             name.j      - expansions and energies for CI solutions    ',&
'             name.debug  - debug information (optional)                ',&
'---------------------------------------------------------------------- ',&
'     The parameters in command line has structure:   param=value       ',&
'     except "name" of case - should be given as the first argument.    ',&
'---------------------------------------------------------------------- ',&
'     List of main parameters with default values:                      ',&
'     (all arguments except "name" are optional)                        ',&
'                                                                       ',&
'     name        - name of case                                        ',&
'     mbreit=1    - include (1) or not (0) Breit oprerators             ',&
'     msol=0      - how many eigensolution you need (0 -> all)          ',&
'     nzero=0     - zero-order dimension (0 -> all configurations)      ',&
'     debug=0     - if /= 0, file name.debug contain debug information  ',&
'                                                                       ',&
'     eps_c    = 1.0D-10  -  tolerance for coefficients                 ',&
'     eps_det  = 1.0D-10  -  tolerance for determinant overlaps         ',&
'     eps_ovl  = 1.0D-08  -  tolerance for configuration overlaps       ',&
'---------------------------------------------------------------------- ',&
'     Example of calling:                                               ',&
'                                                                       ',&
'     dbsr_ci name                                                      ',&
'     dbsr_ci name jmin=1 jmax=3 mbreit=0 msol=1                        ',&
'---------------------------------------------------------------------- ',&
'                                                                       '
      Stop                          
                                    
      End Subroutine inf_dbsr_ci       
                                     
