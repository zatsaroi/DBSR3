!======================================================================
      Subroutine br_file 
!======================================================================
      Use MPI
      Use internal_file

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_BCAST(mlines, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nlines, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(klines, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(nlines.le.0) Return

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if(myid.ne.0) then
       if(allocated(aline)) Deallocate (aline)
       Allocate(aline(mlines))
      end if

      Call MPI_BCAST(aline, mlen*mlines, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_file
 

 
