!======================================================================
!     PROGRAM       D B S R _ P O L         
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!   find the psedo-orbitals to treat the state polarization
!======================================================================
!
!    INPUT  ARGUMENTS:
!
!     klsp  - indexes of partial waves under consideration  
! 
!    INPUT FILES:
!
!     dbsr_par       -  description of partial waves
!     dbsr_mat.nnn   -  interaction matrix
!     cfg.nnn        -  c-file for close-coupling expansion
!
!    OUTPUT FILES:
!   
!     pol_nnn.log    -  running information 
!     bound.nnn      -  pol. solutions
!
!     Above, "nnn" indicates the index of partial wave
!
!=====================================================================
      Use dbsr_pol
      Use DBS_grid
      Use target_jj   
      
      Implicit none
      Real(8) :: t1,t2
      Real(8), external :: RRTC

      t1=RRTC()

! ... read data and arguments:

      Call Read_data

! ... read the H- and S- matrix:

      Call R_dbsrmat

! ... read additional orthogonal constraints:

      Call Add_nortb 
                    
! ... read the d-vector:

      Call R_dvector

! ... set up and solve pseudostate equation

      Call Solv_mat

      t2=RRTC()

      write(pri,'(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'
      write(*,  '(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'

      End  ! program DBSR_POL


