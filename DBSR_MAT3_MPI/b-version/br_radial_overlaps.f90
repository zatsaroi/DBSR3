!======================================================================
      Subroutine br_radial_overlaps
!======================================================================
      Use MPI
      Use radial_overlaps

      Implicit none
      
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(mobs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nobs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(mobs.eq.0) Stop 'mort = 0 in radial orbitals'

      if(myid.ne.0) then
       if(allocated(iobs)) Deallocate(iobs); Allocate(iobs(mobs))
       if(allocated(jobs)) Deallocate(jobs); Allocate(jobs(mobs))
       if(allocated(Cobs)) Deallocate(Cobs); Allocate(Cobs(mobs))
      end if

      Call MPI_BCAST(iobs,mobs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jobs,mobs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Cobs,mobs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_radial_overlaps

