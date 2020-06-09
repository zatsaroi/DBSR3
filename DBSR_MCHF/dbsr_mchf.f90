!======================================================================
!     PROGRAM DBSR_MCHF
!======================================================================
!     MULTICONFIGURATION DIRAC-HARTREE-FOCK PROGRAM
!----------------------------------------------------------------------
!                   C O P Y R I G H T -- 2015
!     Written by:   Oleg Zatsarinny and Charlotte Froese Fischer
!     email:        oleg_zoi@yahoo.com
!----------------------------------------------------------------------
!     This program is a part of the DBSR complex and computes 
!     radial one-electron functions in B-spline basis for the 
!     multi-configuration Dirac-Hartree-Fock problem. 
!
!     For short instructions, run dbsr_mchf atom and look in atom.inp,
!     where 'atom' is any atomic symbol
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Real(8) :: t1,t2

      Call CPU_time(t1)

! ... read main parameters:

      Call Get_case

! ... read configurations from c-file:

      Call Read_conf_jj

      Call Def_blocks

! ... read angular coefficients: 

      Call Read_ang_coef

! ... B-spline parameters:

      Call Get_spline_param

! ... prepare radial-function arrays:

      Call Def_orbitals
      Call Get_estimates
      Call Write_inp

! ... Computing:

      if(all.eq.0)   Call scf_mchf
      if(all.eq.1)   Call scf_mchf_all

! ... output the results:

      Call Output_jfile
      Call Write_pqbs

      Call Summry

      Call CPU_time(t2)
      write(scr,'(/a,T20,f10.2,a)') 'time:',t2-t1,'  sec'
      write(log,'(/a,T20,f10.2,a)') 'time:',t2-t1,'  sec' 

      write(log,'(/80(''-'')/)')

      Call Debug_time

      End ! Program DBSR_MCHF


!======================================================================
      Subroutine Debug_time
!======================================================================
!     timimg different procedures
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals
      Use DBS_debug      

!     if(debug.eq.0) Return

      write(log,'(/a)') 'debug timing:'

      write(log,'(/a,T20,f10.2,a)') 'time_matrix: ',time_matrix,  '  sec' 
      write(log,'( a,T20,f10.2,a)') 'time_solve:  ',time_solve,   '  sec' 
      write(log,'( a,T20,f10.2,a)') 'time_diag:   ',time_diag,    '  sec' 
      if(time_read_coef.gt.5.d0) &
      write(log,'( a,T20,f10.2,a)') 'time_read_coef:',time_read_coef,'  sec' 

      write(log,'(/a,T20,f10.2,a)') 'time_convolution:',time_convol, '  sec' 
      write(log,'( a,T20,f10.2,a)') 'time_density:    ',time_density,'  sec' 

      write(log,'(/a,i10)') 'Calls to convol  = ',ic_convol
      write(log,'( a,i10)') 'Calls to density = ',ic_density

      End Subroutine Debug_time

