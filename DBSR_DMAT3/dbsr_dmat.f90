!======================================================================
!     PROGRAM       D B S R _ D M A T
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny (pilgrim)
!     email:        oleg_zoi@yahoo.com
!
!======================================================================
!     Generation of the multipole-transition matrix 
!======================================================================
!
!     INPUT ARGUMENTS:
!
!   1.  name1.c or cfg.nnn - c-file for initial state 
!   2.  name2.c or cfg.nnn - c-file for final state 
!   3.  ctype1 - structure mode for initial state  (c,j,b)
!   4.  ctype2 - structure mode for final state    (c,j,b,p)
!   5.  if ctype2 = p, pointer on the first state in files j,b,
!       otherwise the first states will be considered
! 
!----------------------------------------------------------------------
!
!     INPUT FILES:
!
!     mult_bnk_nn - data bank of angular coefficients for multipole 
!                   operator (nn -> E1,M1,E2,...)
!     knot.dat - B-spline information
!
!     name1.c  - c-file for initial state 
!     name1.bsw - initial state wavefunctions in B-spline basis
!
!     name2.c  - c-file for final state
!     name2.bsw - initial state wavefunctions in B-spline basis
!
!     One (or both) c-files can be replaced by cfg.nnn files,
!     where nnn is a partial wave number.
!     The close-coupling expansion is supposed for this case,
!     that requires also following files
!
!     target  - information about close-coupling expansion
!               (target states, channels and perturbers)
!     target.bsw  - target wavefunctions in B-spline basis
!     perturber.bsw - if any (according to target information)
!
!-----------------------------------------------------------------------
!
!     OUTPUT FILES:
!
!     dbsr_dmat.log - running information
!     d.dat         - dipole matrix (if ctype2 = p)
!     zf_res        - oscilator strengths (all other cases)
!
!-----------------------------------------------------------------------
!
!     RESTRICTIONS:  in b->j or j->b calculations, j-file can contain
!                    data only for one value of j.
! 
!----------------------------------------------------------------------
      Use dbsr_dmat

      Implicit none
      Real(8) :: t1,t2

      Call CPU_time(t1)

! ... read arguments and data files:

      Call Read_data

! ... generation D-matrix:

      Call D_matr

! ... output D-marix or f-value:

      if(ctype2.eq.'p') then
       Call D_out
      elseif(ctype2.eq.'q') then
       Call DV_out
      else
       Call Gen_zf
      end if

      Call CPU_time(t2)
      write(pri,'(a,T30,f10.2,a)') 'DBSR_DMAT:  ',(t2-t1)/60,' min '
      write(*  ,'(a,T30,f10.2,a)') 'DBSR_DMAT:  ',(t2-t1)/60,' min '

      End ! program DBSR_DMAT



