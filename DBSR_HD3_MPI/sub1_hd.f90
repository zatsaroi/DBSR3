!=====================================================================
      Subroutine SUB1_HD    
!=====================================================================
!     drive calculations for one partial wave
!---------------------------------------------------------------------
      Use dbsr_hd

      Implicit none

! ... read channel information:

      if(io_processor) Call Read_channel_jj(nut,klsp)

! ... check and print main parameters and file:

      if(io_processor) Call Check_dbsr_mat

! ... broadcast main parameters:

      Call Br_ipar(itype)
      Call Br_ipar(iwt)
      Call Br_dpar(cwt)
      Call Br_ipar(iexp)
      Call Br_ipar(npert)
      Call Br_ipar(nch)
      Call Br_ipar(ms)
      mhm = ms*nch+npert
      Call Br_ipar(ksol)
      khm = ksol + npert
      Call Br_ipar(khm)

      if(io_processor) then
       Call igebs2d(ctxt,'all',' ',nch+1,1,ipsol,nch+1)
      else
       if(allocated(ipsol)) Deallocate(ipsol); Allocate(ipsol(0:nch))
       ipsol = 0
       Call igebr2d(ctxt,'all',' ',nch+1,1,ipsol,nch+1,rsrc,csrc)
      end if    

! ... diagonalize the matrix and get inner-region solutions:

      Call Diag_mat

! ... output of solutions and find the surface amplitudes:

      if(itype.ge. 0)  then
       Call rsol_out
       Call cpu_time (t1)
       if(io_processor) &
       write (*  ,'(/a,T20,f10.2,a)') 'rsol_out:', (t1-t0)/60, ' min.'
       if(io_processor) &
       write (pri,'(/a,T20,f10.2,a)') 'rsol_out:', (t1-t0)/60, ' min.'
      end if

! ... output of standard H.nnn file:

      if(itype.ge.0.and.io_processor) then
       if(iiexp.eq.0) Call H_out
       if(iiexp.eq.1) Call H_out_exp
       Call cpu_time (t1)
       if(io_processor) &
       write (*  ,'(/a,T20,f10.2,a)') 'H_out:', (t1-t0)/60, ' min.'
       if(io_processor) &
       write (pri,'(/a,T20,f10.2,a)') 'H_out:', (t1-t0)/60, ' min.'
      end if

! ... output of bound states in bound.nnn:

      if(itype.eq.-1) then
       Call B_out 
       Call cpu_time (t1)
       if(io_processor) &
       write (*  ,'(/a,T20,f10.2,a)') 'B_out:', (t1-t0)/60, ' min.'
       if(io_processor) &
       write (pri,'(/a,T20,f10.2,a)') 'B_out:', (t1-t0)/60, ' min.'
      end if

      End Subroutine SUB1_HD





