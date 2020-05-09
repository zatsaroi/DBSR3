!======================================================================
      Subroutine br_channel_jj
!======================================================================
!     broadcast the data in module "channel_jj"
!----------------------------------------------------------------------
      Use MPI
      Use channel_jj

      Implicit none
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(ipar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jpar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nch, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mch, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ncp, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwp, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(kch   )) Deallocate(kch   ); Allocate(kch   (mch))
       if(allocated(lch   )) Deallocate(lch   ); Allocate(lch   (mch))
       if(allocated(jjch  )) Deallocate(jjch  ); Allocate(jjch  (mch))
       if(allocated(iptar )) Deallocate(iptar ); Allocate(iptar (mch))
       if(allocated(ipconf)) Deallocate(ipconf); Allocate(ipconf(mch))
       if(allocated(ELC   )) Deallocate(ELC   ); Allocate(ELC   (mch))
       if(allocated(ipch  )) Deallocate(ipch  ); Allocate(ipch  (mch))
      end if

      Call MPI_BCAST(kch   ,mch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lch   ,mch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jjch  ,mch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipch  ,mch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iptar ,mch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipconf,mch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ELC,mch*5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(AFP,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(BFP,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Tpar,4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(mpert, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(npert, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(npert.ne.0) then
       if(myid.ne.0) then
        if(allocated(ippert)) Deallocate(ippert); Allocate(ippert(0:mpert))
        ippert = 0
       end if
       Call MPI_BCAST(ippert,mpert+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      end if

      End Subroutine br_channel_jj



