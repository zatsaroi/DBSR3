!======================================================================
      Subroutine br_buffer
!======================================================================
!     broadcast the main parameters used in dbsr_mat program
!----------------------------------------------------------------------
      Use MPI
      Use dbsr_mat

      Implicit none

! ... broadcast the buffer:

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(ncbuf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(itb  ,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jtb  ,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(intb ,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(idfb ,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Cbuf ,ncbuf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_buffer
