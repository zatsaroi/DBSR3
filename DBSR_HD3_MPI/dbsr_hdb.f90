!======================================================================
!    PROGRAM       D B S R _ H D 3       
!
!               C O P Y R I G H T -- 2019
!
!    Written by:   Oleg Zatsarinny
!                  email: oleg_zoi@yahoo.com
!
!======================================================================
!   Generates H.DAT and optioanlly W.DAT, RSOL.DAT fails for further
!   scattering or photoionization calculations in relativistic R-matrix
!   approach.  Another option - calculation of bound and pseudo-states.
!======================================================================
!
!    INPUT FILES:
!
!     command line   -  parameters as argument "name=value"
!     dbsr_par       -  description of partial waves 
!     dbsr_mat.nnn   -  interaction matrix
!     cfg.nnn        -  c-file for close-coupling expansion
!
!    OUTPUT FILES:   
!   
!     dbsr_hd.nnn    -  running information 
!     bound.nnn      -  bound-like solusions
!     h.nnn          -  eigenvalues and suface amplitudes
!                       required to obtain R-matrix
!     rsol.nnn       -  R-matrix inner-region solutions 
!     w.nnn          -  channel weights 
!
!     "nnn" indicates the index of partial wave
!
!=====================================================================
      Use dbsr_hd
      
      Implicit none

!----------------------------------------------------------------------
! ... define the blacs grid:

      Call DEF_BLACS_GRID

!----------------------------------------------------------------------
! ... read common parameters and data:

      if(io_processor) then

! ... read target information:

      Call Check_file(AF_tar); Open(nut,file=AF_tar)
      Call Read_target_jj(nut)

! ... set up B-splines:

      Call read_knot_dat
      Call alloc_DBS_gauss

! ... read arguments:

      Open(nup,file=AF_par);  Call Read_arg(nup);  Close(nup)

      end if  ! io_processor

      Call br_ipar(klsp1)
      Call br_ipar(klsp2)

!-------------------------------------------------------------------------- 
      if (myrow >= 0 .and. myrow < p .and. mycol >= 0 .and. mycol < q) then
 
       Do klsp = klsp1,klsp2

        Call BLACS_BARRIER (ctxt, 'All')
        Call cpu_time (t0)

        Call SUB1_HD

        Call BLACS_BARRIER (ctxt, 'All')
        Call cpu_time (t1)

        if(io_processor) then
         write (pri,'(/a,T20,f10.2,a)') 'CPU time     = ', (t1-t0)/60, ' min.'
         write (*  ,'(/a,T20,f10.2,a)') 'CPU time     = ', (t1-t0)/60, ' min.'
        end if

       End do     

      end if

      Call BLACS_BARRIER (ctxt, 'All')
      Call BLACS_GRIDEXIT (ctxt)    
      Call BLACS_EXIT()

      End  ! program DBSR_HD


!=====================================================================
      Subroutine Def_blacs_grid
!=====================================================================
! ... define BLACS grid
! ... (this procedure when created first linear grid is after NOBEL 
! ... and looks too complicated)
!---------------------------------------------------------------------
      Use blacs
      Use dbsr_hd, only: nup,AF_par
      
      Implicit none
      
      Integer :: imycol, imyrow, ip, iq
      Integer, external :: Icheck_file

!----------------------------------------------------------------------
! ... define linear blacs grid to include all processes to broadcast:

      Call BLACS_PINFO (iam, nprocs) ! find process #, total # processors
      Call BLACS_GET (-1, 0, ictxt)  ! find default context, ictxt

      Call BLACS_GRIDINIT (ictxt, 'Row-major', 1, nprocs)
      Call BLACS_GRIDINFO (ictxt, ip, iq, imyrow, imycol)

      io_processor = (imycol == 0)

      if (io_processor) then 

       write (*,*)
       write (*,'(a,i5)') 'BSR_HDB: number of processors = ', nprocs

! ... read p,q,nblocks if any:

       if(Icheck_file(AF_par).ne.0) then
        Open(nup,file=AF_par)
        Call Read_ipar(nup,'p'     ,p     )
        Call Read_ipar(nup,'q'     ,q     )
        Call Read_ipar(nup,'nblock',nblock)
        Close(nup)    
       end if
       Call Read_iarg('p'     ,p     )
       Call Read_iarg('q'     ,q     )
       Call Read_iarg('nblock',nblock)

! ... check the grid parameters:

       if(p*q.eq.0) then; p=INT(SQRT(REAL(nprocs))); q=p; end if

       write (*,*)
       write (*,'(a,3i5)') 'Blacks dimensions, p,q,nblock = ', &
                            p,q,nblock

      end if  ! over io_processor

      Call br_ipar(p)
      Call br_ipar(q)
      Call br_ipar(nblock)

      if(p*q > nprocs) then; Call BLACS_EXIT(); Stop 'p*q > nprocs'; end if

      Call BLACS_BARRIER (ictxt, 'All')
      Call BLACS_GRIDEXIT (ictxt)           ! kill initial mesh

!----------------------------------------------------------------------
! ... create p-q BLACS grid:

      Call BLACS_GET (-1,0,ctxt)    
      Call BLACS_GRIDINIT (ctxt, 'Row-major', p, q)
      Call BLACS_GRIDINFO (ctxt, ip, iq, myrow, mycol)

      io_processor = (myrow == 0 .and. mycol == 0)

      End Subroutine Def_blacs_grid


