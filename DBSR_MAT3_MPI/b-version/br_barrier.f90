!======================================================================
      Subroutine br_barrier
!======================================================================
!     barrier the calculations
!----------------------------------------------------------------------
      Use MPI

      Implicit none
      Integer :: ierr

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_barrier