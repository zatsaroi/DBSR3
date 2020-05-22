!======================================================================
      Subroutine br_target_jj
!======================================================================
!     broadcast the data in module "target_jj"
!----------------------------------------------------------------------
      Use MPI
      Use target_jj

      Implicit none
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(ntarg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nelc, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nz,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nct,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwt,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if(myid.ne.0) then

       if(allocated(jtarg)) Deallocate(jtarg); Allocate(jtarg(ntarg))
       if(allocated(ptarg)) Deallocate(ptarg); Allocate(ptarg(ntarg))
       if(allocated(etarg)) Deallocate(etarg); Allocate(etarg(ntarg)); etarg=0.d0

       if(allocated(nctarg)) Deallocate(nctarg); Allocate(nctarg(ntarg))
       if(allocated(nwtarg)) Deallocate(nwtarg); Allocate(nwtarg(ntarg))
       if(allocated(ictarg)) Deallocate(ictarg); Allocate(ictarg(ntarg))

       if(allocated(AFT)) Deallocate(AFT); Allocate(AFT(ntarg))
       if(allocated(BFT)) Deallocate(BFT); Allocate(BFT(ntarg))

      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(jtarg, ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ptarg, ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nctarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwtarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ictarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(AFT,ntarg*80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(BFT,ntarg*80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(etarg,ntarg,MPI_DOUBLE_PRECISION,0, &
                     MPI_COMM_WORLD,ierr)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_target_jj


