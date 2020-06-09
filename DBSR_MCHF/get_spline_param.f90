!======================================================================
     Subroutine get_spline_param 
!======================================================================
!    This routine gets information about the spline parameters and
!    and set up spline knots and all relevant arrays.
!    Program checks file AF_grid (knot.dat or name.knot)
!    for B-spline parameters  or use default ones.
!    The spline paremeters can also be modified from command line. 
!----------------------------------------------------------------------
     Use dbsr_mchf
     Use DBS_integrals, only: ntype, memory_DBS_integrals

     Implicit none

! .. define the knot sequence:

     AF_knot = trim(name)//'.knot'
     if(Icheck_file(AF_knot).eq.0) AF_knot = BF_knot

     Call def_grid(AF_knot,name,z,atw)

! .. initialize other B-spline arrays:

     Call alloc_DBS_gauss
     Call def_Vnucl
     Call alloc_DF_radial(ns)

     Call alloc_Rk_integrals(ns,ks,0,kmax,4)     

     write(log,'(/80(''-''))')
     write(log,'(/a/)') 'B-spline parameters:'
     write(log,'(a,i5,t25,a)')   'ns   =',ns,  '-  number of splines' 
     write(log,'(a,i5,t25,a)')   'ks   =',ks,  '-  order  of splines' 
     write(log,'(a,i5,t25,a)')   'ksp  =',ksp, '-  for large component' 
     write(log,'(a,i5,t25,a)')   'ksq  =',ksq, '-  for small component'
     write(log,'(a,f9.3,t25,a)') 'he   =',he,  '-  exponetial step size' 
     write(log,'(a,f9.3,t25,a)') 'tmax =',tmax,'-  maximum radius' 

     write(log,'(/a/)') 'Spline integrals:'
     write(log,'(a,i5,t25,a)')   'kmin =',kmin,  '-  min. multipole index' 
     write(log,'(a,i5,t25,a)')   'kmax =',kmax,  '-  max. multipole index' 
     write(log,'(a,i5,t25,a)')   'ntype=',ntype, '-  number of different integrals' 

     write(log,'(/a,T40,F10.2,a)') 'Memory required for spline integrals: ', &
                                    memory_DBS_integrals,' Mb' 

     End Subroutine get_spline_param 

