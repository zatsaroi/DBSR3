!======================================================================
      Subroutine diag_mat
!======================================================================
!     Diagonalization procedure:
!
!     1. Read overlap matrix, b.
!     2. Read Hamiltonian matrix, a.
!     3. Form a Cholesky factorization of the overlap matrix.
!     4. Transform problem to standard eigenvalue problem.
!     5. Transform to experimental thresholds energies if any.
!     6. Solve standard eigenvalue problem.
!     7. Call w_out for weights if required.
!     8. Backtransform eigenvectors to the original problem.
!        Results in array "a".
!----------------------------------------------------------------------
      Use dbsr_hd 

      Implicit none
      Character(1) :: uplo, job, trans
      Integer ::  info, ibtype
      Integer ::  i,j, ich,jch, i1,i2, j1,j2, ip,jp

!----------------------------------------------------------------------
! ... read overlap matrix:

      if(allocated(b)) deallocate(b); allocate(b(khm,khm))
      b=zero; Do i=1,ksol; b(i,i)=one; End do

      diag_ovl = 1

   1  read(nui) ich,jch
      if(ich.eq.0) go to 2
      diag_ovl = 0
      if(ich.le.nch.and.jch.le.nch) then
       i1=ipsol(ich-1)+1; i2=ipsol(ich)
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       read(nui) b(i1:i2,j1:j2)
      elseif(ich.gt.nch.and.jch.le.nch) then
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       ip=ipsol(nch)+ich-nch
       read(nui) b(ip,j1:j2)
      else
       ip=ipsol(nch)+ich-nch
       jp=ipsol(nch)+jch-nch
       read(nui) b(ip,jp)
      end if
      go to 1

   2  Continue

      if(diag_ovl.eq.1) &
      write(pri,'(/a,i8)') 'It is a standard eigenvalue problem: khm =',khm
      if(diag_ovl.eq.0) &
      write(pri,'(/a,i8)') 'It is a generalized eigenvalue problem: khm =',khm

!----------------------------------------------------------------------
! ... read Hamiltonian matrix:

      if(allocated(a)) deallocate(a); allocate(a(khm,khm))
      a=zero; Do i=1,ksol; a(i,i)=bval(i); End do

   3  read(nui) ich,jch
      if(ich.eq.0) go to 4

      if(ich.le.nch.and.jch.le.nch) then
       i1=ipsol(ich-1)+1; i2=ipsol(ich)
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       read(nui) a(i1:i2,j1:j2)
      elseif(ich.gt.nch.and.jch.le.nch) then
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       ip=ipsol(nch)+ich-nch
       read(nui) a(ip,j1:j2)
      else
       ip=ipsol(nch)+ich-nch
       jp=ipsol(nch)+jch-nch
       read(nui) a(ip,jp)
      end if
      go to 3

   4  Continue

!----------------------------------------------------------------------
! ... diagonalize using the SCALAPACK routines:
!----------------------------------------------------------------------
! ... default parameters:

      ibtype =  1    ! A * x = lambda * B * x generalized eigenproblem
      job    = 'V'   ! compute both eigenvalues and eigenvectors
      uplo   = 'L'   ! use lower triangles of A and B matrices
      trans  = 'T'   ! when uplo=L, otherwise trans = 'N'

!----------------------------------------------------------------------
! ... Form a Cholesky factorization of the overlap matrix:
! ... The factorization has the form
! ...   A = U^T * U,    if UPLO = 'U'
! ...   A = L   * L^T,  if UPLO = 'L'
! ... (other part of matrix is not touched)

      if(diag_ovl.eq.0) then

       Call DPOTRF (uplo, khm, b, khm, info)

       if(info.ne.0) Stop 'DPOTRF info > 0 - Cholesky factorization failed'

      end if
!----------------------------------------------------------------------
! ... Transform problem to standard eigenvalue problem:
! ...    A ->  inv(L)*A*inv(L^T),  uplo = 'L'
! ....   A ->  inv(U^T)*A*inv(U),  uplo = 'U'

      if(diag_ovl.eq.0) then

       Call DSYGST( ibtype, uplo, khm, A, khm, B, khm, INFO )

       if(info.ne.0) Stop 'DSYGST: info > 0 '

      end if

!----------------------------------------------------------------------
! ... introduction of experimental energies:

      if(iexp.gt.0) Call Add_exp

!-----------------------------------------------------------------------
! ... Solve standard eigenvalue problem
! ... (note:  divide and concer algorith requires much more space)

      write(*,'(/a,i3,a,i5,a,i5)') &
            'DBSR_HD:  klsp =',klsp,'   khm =',khm

      if(allocated(eval)) deallocate(eval); allocate(eval(khm))

      Call LAP_DSYEV(job,uplo,khm,khm,A,eval,INFO)

      if(info.ne.0) Stop 'DSYEV: info > 0 '

      write(*  ,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)
      write(pri,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)

!-----------------------------------------------------------------------
!...  define weights:

      if(itype.eq.-1.or.cwt.gt.zero.or.iwt.gt.0) Call W_out

!-----------------------------------------------------------------------
! ... Backtransform eigenvectors to the original problem.
! ... For A*x=(lambda)*B*x backtransform eigenvectors:
! ... x = inv(L)'*y or inv(U)*y

      if(diag_ovl.eq.0) then

       Call DTRSM('left', uplo, trans, 'Non-unit', khm, khm, ONE, &
                  B, khm, A, khm )

      end if ! over diag_ovl

      Deallocate(b)

      End Subroutine Diag_mat
