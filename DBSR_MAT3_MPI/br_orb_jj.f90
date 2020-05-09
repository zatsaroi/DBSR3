!======================================================================
      Subroutine br_orb_jj
!======================================================================
!     broadcast the data in modules "orb_jj, symc_list, symt_list, conf_jj"
!----------------------------------------------------------------------
      Use MPI
      Use orb_jj

      Implicit none
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

! ... orb list:

      Call MPI_BCAST(nwf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mwf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(NEF )) Deallocate(NEF ); Allocate(NEF (mwf))
       if(allocated(KEF )) Deallocate(KEF ); Allocate(KEF (mwf))
       if(allocated(IEF )) Deallocate(IEF ); Allocate(IEF (mwf))
       if(allocated(ipef)) Deallocate(ipef); Allocate(ipef(mwf))
       if(allocated(ELF )) Deallocate(ELF ); Allocate(ELF (mwf))
      end if

      Call MPI_BCAST(NEF,mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(KEF,mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IEF,mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ELF,mwf*5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      ipef = 0

      End Subroutine br_orb_jj


