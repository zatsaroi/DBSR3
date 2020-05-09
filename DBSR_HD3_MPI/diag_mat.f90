!======================================================================
      Subroutine diag_mat
!======================================================================
!     Diagonalization procedure:
!
!     1. Read transformed overlap and interaction matrixes.
!     2. Form a Cholesky factorization of the overlap matrix if any.
!     3. Transform problem to standard eigenvalue problem.
!     4. Transform to experimental thresholds energies if any.
!     5. Solve standard eigenvalue problem. 
!     6. Call w_out for weights if required.
!     7. Backtransform eigenvectors to the original generalized problem.
!        Solutions are still in "diagonal-blocks" basis!!!
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: info,ibtype, lwork, status
      Real(8) :: scale=one, S
      Real(8), allocatable :: work(:)
      Integer, external :: numroc
      Character(1) :: uplo, job, trans, range

!----------------------------------------------------------------------
! ... allocate working arrays for distributions:

      Call BLACS_BARRIER (ctxt, 'All')

      Call descinit (descadd,ms,ms,ms,ms,rsrc,csrc,ctxt,ms,info)
      Call p_error  (info,'descadd descriptor error in diag_hd')
      if(allocated(add)) Deallocate(add); Allocate(add(ms,ms)); add=zero

      if(npert.gt.0) then
       Call descinit(descadp,npert,npert,npert,npert,rsrc,csrc,ctxt,npert,info)
       Call p_error(info,'descadp descriptor error in diag_hd')
       if(allocated(adp)) Deallocate(adp); Allocate(adp(npert,npert))
      end if

!... distribution local dimensions:

      np = numroc (khm, nblock, myrow, rsrc, p)
      nq = numroc (khm, nblock, mycol, csrc, q)
      ld = MAX(1, np)

      if(io_processor) &
      write (*,'(/a,T20,2i10)') 'Matrix np,nq:', np,nq

      S = 16.d0*np*nq/(1024*1024)
      if(io_processor) &
      write (*,'(/a,T20,f10.1,a)') 'Matrix memory:', S, ' Mb'

      Call descinit (descv, khm, 1, khm, 1, rsrc, csrc, ctxt, khm,info)
      Call p_error(info,'descv descriptor error in diag_hd')
      if(allocated(v)) Deallocate(v); Allocate(v(khm))     

!----------------------------------------------------------------------
! ... read transformed overlap matrixes:

      Call Read_overlaps

      Call Read_matrix

      Call cpu_time (t1)
      if(io_processor) &
      write (*  ,'(/a,T20,f10.2,a)') 'Read_matrix:', (t1-t0)/60, ' min.'
      if(io_processor) &
      write (pri,'(/a,T20,f10.2,a)') 'Read_matrix:', (t1-t0)/60, ' min.'

!----------------------------------------------------------------------
! ... diagonalize using the SCALAPACK routines:
!----------------------------------------------------------------------
! ... default parameters:
 
      ibtype =  1    ! A * x = lambda * B * x generalized eigenproblem
      range  = 'A'   ! compute all eigenvalues
      job    = 'V'   ! compute both eigenvalues and eigenvectors
      uplo   = 'L'   ! use lower triangles of A and B matrices
      trans  = 'T'   ! 

!----------------------------------------------------------------------
! ... Form a Cholesky factorization of the overlap matrix:
! ... The factorization has the form
! ...   A = U^T * U,    if UPLO = 'U'
! ...   A = L   * L^T,  if UPLO = 'L'
! ... (other part of matrix is not touched)

      if(diag_ovl.gt.0) then

      Call PDPOTRF (uplo, khm, b, 1, 1, descb, info)

      Call BLACS_BARRIER (ctxt, 'all')

      Call p_error (info, 'pdpotrf: Cholesky error')

      Call cpu_time (t1)
      if(io_processor) &
      write (*  ,'(/a,T20,f10.2,a)') 'Choleski decomposition:', (t1-t0)/60, ' min.'
      if(io_processor) &
      write (pri,'(/a,T20,f10.2,a)') 'Choleski decomposition:', (t1-t0)/60, ' min.'

      ! ... Transform problem to standard eigenvalue problem:
      ! ...    A ->  inv(L)*A*inv(L^T),  uplo = 'L'
      ! ....   A ->  inv(U^T)*A*inv(U),  uplo = 'U'

      lwork = 2 
      allocate (work(lwork), stat=status)
      Call p_error (status, 'error memory allocation of work array')

      lwork = -1

      Call PDSYNGST (ibtype, uplo,  khm,a,1,1,desca, b,1,1,descb, &
                   &  scale, work, lwork, info)

      lwork=int(work(1)+1); deallocate (work); allocate(work(lwork))

      Call PDSYNGST (ibtype, uplo,  khm,a,1,1,desca, b,1,1,descb, &
                   &  scale, work, lwork, info)

      Call p_error (info, 'pdsyngst error')
      Deallocate(work) 

      end if  ! over diag_ovl

      Call cpu_time (t1)
      if(io_processor) &
      write (*  ,'(/a,T20,f10.2,a)') 'Trans_to_stand:', (t1-t0)/60, ' min.'
      if(io_processor) &
      write (pri,'(/a,T20,f10.2,a)') 'Trans_to_stand:', (t1-t0)/60, ' min.'

!----------------------------------------------------------------------
! ... use experimental target energies if any:

      Call BLACS_BARRIER (ctxt, 'all')

      if(iexp.gt.0) Call Add_exp

      Call BLACS_BARRIER (ctxt, 'all')

      Call cpu_time (t1)
      if(io_processor.and.iexp.gt.0) &
      write (*  ,'(/a,T20,f10.2,a)') 'Add_exp:', (t1-t0)/60, ' min.'
      if(io_processor.and.iexp.gt.0) &
      write (pri,'(/a,T20,f10.2,a)') 'Add_exp:', (t1-t0)/60, ' min.'

!-----------------------------------------------------------------------
! ... solve standard eigenvalue problem 
! ... (note:  divide and concer algorith requires much more space)

      Call descinit(descz,khm,khm,nblock,nblock,rsrc,csrc,ctxt,ld,info)
      Call p_error(info,'descz descriptor error')
      allocate (z(np,nq));  z=zero

      if(allocated(eval)) Deallocate(eval); Allocate(eval(khm))

      if(.not.allocated(work)) Allocate(work(10))

      lwork = -1
      Call PDSYEV (job, uplo, khm,a,1,1,desca, eval, z,1,1,descz, &
                   work, lwork,  info)

      lwork = INT(work(1) + 1);  Deallocate(work); Allocate(work(lwork))

      if(io_processor) &
      write(*  ,'(/a,i3,a,i10,a,i10,a,i10)')  'Matrix:   klsp =',klsp, &
          '  nhm =',nhm,'  khm =',khm,'  lwork =',lwork
      if(io_processor) &
      write(pri,'(/a,i3,a,i10,a,i10,a,i10)')  'Matrix:   klsp =',klsp, &
          '  nhm =',nhm,'  khm =',khm,'  lwork =',lwork

      Call PDSYEV (job, uplo, khm,a,1,1,desca, eval, z,1,1,descz, &
                    work, lwork, info)
      Call p_error (info, 'pdsyev error')

      Deallocate (work)

      if(io_processor) write(  *,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)
      if(io_processor) write(pri,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)

      Call cpu_time (t1)
      if(io_processor) &
      write (*  ,'(/a,T20,f10.2,a)') 'Diagonalization:', (t1-t0)/60, ' min.'
      if(io_processor) &
      write (pri,'(/a,T20,f10.2,a)') 'Diagonalization:', (t1-t0)/60, ' min.'

      if(allocated(a)) Deallocate(a)

!-----------------------------------------------------------------------
!...  define weights:

      if(itype.ne.0.or.cwt.gt.zero.or.iwt.gt.0) then
       Call W_out 
       Call cpu_time (t1)
       if(io_processor) &
       write (*  ,'(/a,T20,f10.2,a)') 'W_out:', (t1-t0)/60, ' min.'
       if(io_processor) &
       write (pri,'(/a,T20,f10.2,a)') 'W_out:', (t1-t0)/60, ' min.'
      end if

!-----------------------------------------------------------------------
! ... Back transform eigenvectors to the original problem.
! ... For A*x=(lambda)*B*x backtransform eigenvectors:
! ... x = inv(L)'*y or inv(U)*y

      if(diag_ovl.gt.0) then
       Call PDTRSM ('Left', uplo, trans, 'Non-unit', khm, khm, one, &
                     b, 1,1, descb, z, 1,1, descz)
       if (scale /= one) Call DSCAL(khm, scale, eval, 1)
      end if ! over diag_ovl

      if(allocated(b)) Deallocate(b)

      Call BLACS_BARRIER (ctxt, 'all')

      Call cpu_time (t1)
      if(io_processor) &
      write (*  ,'(/a,T20,f10.2,a)') 'Back_transformation:', (t1-t0)/60, ' min.'
      if(io_processor) &
      write (pri,'(/a,T20,f10.2,a)') 'Back_transformation:', (t1-t0)/60, ' min.'

      End Subroutine Diag_mat





