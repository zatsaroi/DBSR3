!======================================================================
      Module Blacs
!======================================================================
! ... store blacs context and processor data
!----------------------------------------------------------------------
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Let A be a generic term for any 2D block cyclicly distributed array.
!  Such a global array has an associated description vector DESCA.
!  In the following comments, the character _ should be read as
!  "of the global array".
!
!  NOTATION        STORED IN      EXPLANATION
!  --------------- -------------- --------------------------------------
!  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!                                 DTYPE_A = 1.
!  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!                                 the BLACS process grid A is distribu-
!                                 ted over. The context itself is glo-
!                                 bal, but the handle (the Integer
!                                 value) may vary.
!  M_A    (global) DESCA( M_ )    The number of rows in the global
!                                 array A.
!  N_A    (global) DESCA( N_ )    The number of columns in the global
!                                 array A.
!  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!                                 the rows of the array.
!  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!                                 the columns of the array.
!  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!                                 row of the array A is distributed.
!  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!                                 first column of the array A is
!                                 distributed.
!  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
!
!  Let K be the number of rows or columns of a distributed matrix,
!  and assume that its process grid has dimension p x q.
!  LOCr( K ) denotes the number of elements of K that a process
!  would receive if K were distributed over the p processes of its
!  process column.
!  Similarly, LOCc( K ) denotes the number of elements of K that a
!  process would receive if K were distributed over the q processes of
!  its process row.
!  The values of LOCr() and LOCc() may be determined via a call to the
!  ScaLAPACK tool function, NUMROC:
!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!----------------------------------------------------------------------
      Implicit none

! ... parameters for BLACS distributed arrays:

      Integer, parameter :: dlen_ = 9
      Integer, parameter :: block_cyclic_2d = 1
      Integer, parameter :: ctxt_ = 2
      Integer, parameter :: m_ = 3
      Integer, parameter :: n_ = 4
      Integer, parameter :: mb_ = 5
      Integer, parameter :: nb_ = 6
      Integer, parameter :: rsrc_ = 7
      Integer, parameter :: csrc_ = 8
      Integer, parameter :: lld_ = 9

! ... initial linear context:

      Integer :: ictxt  ! initial blacs context handle
      Integer :: iam    ! processor id in initial context
      Integer :: nprocs ! total # processors

! ... context for main calculation: p * q <= nprocs

      Integer :: ctxt           ! blacs context handle
      Integer :: p=1, q=1       ! grid parameters
      Integer :: mycol, myrow   ! processor coordinates
      Logical :: io_processor   ! i/o processor flag
      Integer :: nblock = 64    ! blacs blocking factor
      Integer :: rsrc=0         ! raws source
      Integer :: csrc=0         ! cols source
      Integer :: np=1           ! local # rows
      Integer :: nq=1           ! local # cols
      Integer :: ld=1           ! leading dimension

      End module Blacs


!======================================================================
      Subroutine p_error (err, msg)
!======================================================================
!     stop BLACS with message printing
!---------------------------------------------------------------------      
      Use blacs

      Implicit none
      Integer, intent(in)          :: err     ! flag
      Character(len=*), intent(in) :: msg     ! error message

      Call igsum2d (ctxt, 'All', ' ', 1, 1, err, 1, -1, -1)
      if (err /= 0) then
       if (io_processor) write(*,'(a)') msg
       Call BLACS_ABORT (ctxt, err)
      end if

      End subroutine p_error


!======================================================================
      Subroutine br_ipar(i)
!======================================================================
!     broadcast the value of "i" - integer
!----------------------------------------------------------------------
      Use blacs
      Implicit none
      Integer :: i, par(1)

      if(io_processor) then
       par(1)=i
       Call igebs2d (ctxt, 'all', ' ', 1, 1, par, 1)
      else      
       Call igebr2d (ctxt, 'all', ' ', 1, 1, par, 1, 0, 0)
       i = par(1)
      end if

      End Subroutine br_ipar


!======================================================================
      Subroutine br_dpar(i)
!======================================================================
!     broadcast the value of "i" - double precision
!----------------------------------------------------------------------
      Use blacs
      Implicit none
      Real(8) :: i, par(1)

      if(io_processor) then
       par(1)=i
       Call dgebs2d (ctxt, 'all', ' ', 1, 1, par, 1)
      else      
       Call dgebr2d (ctxt, 'all', ' ', 1, 1, par, 1, 0, 0)
       i = par(1)
      end if

      End Subroutine br_dpar
