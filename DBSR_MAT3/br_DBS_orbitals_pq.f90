!======================================================================
      Subroutine br_DBS_orbitals_pq
!======================================================================
!     broadcast the one-electron orbitals in module "DBS_orbitals"
!----------------------------------------------------------------------
      Use MPI
      Use DBS_orbitals_pq 
      Use DBS_grid, only: ns

      Implicit none
      Integer :: myid,ierr,m,i

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      if(myid.ne.0) Call alloc_DBS_orbitals_pq(0,ns)

      Call MPI_BCAST(mbf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       m = mbf; Call alloc_DBS_orbitals_pq(m,ns)
      end if

      Call MPI_BCAST(nbf ,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ibs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipbs,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ebs,nbf*5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(pq ,ns*2*nbf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(bpq,ns*2*nbf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(nv_ch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(nv_ch.gt.0) then
       if(myid.ne.0) then 
        if(allocated(i_ch)) deallocate(i_ch); Allocate(i_ch(nv_ch))
        if(allocated(j_ch)) deallocate(j_ch); Allocate(j_ch(nv_ch))
       end if
       Call MPI_BCAST(i_ch,nv_ch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(j_ch,nv_ch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      end if

      End Subroutine br_DBS_orbitals_pq
