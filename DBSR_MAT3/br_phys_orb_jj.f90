!======================================================================
      Subroutine br_phys_orb(ntarg)
!======================================================================
      Use MPI
      Use phys_orb_jj

      Implicit none
      
      Integer :: myid,ierr, m, ntarg

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_BCAST(nphys_orb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(nphys_orb.eq.0) Stop 'nphys_orb=0'

      m = nphys_orb
      if(myid.ne.0) then
       if(allocated(ip_tar)) Deallocate(ip_tar); Allocate(ip_tar(ntarg))
       if(allocated(ip_phy)) Deallocate(ip_phy); Allocate(ip_phy(m))
       if(allocated(ip_sub)) Deallocate(ip_sub); Allocate(ip_sub(m))
       if(allocated(s_orb )) Deallocate(s_orb ); Allocate(s_orb (m))
      end if

      Call MPI_BCAST(ip_tar,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_phy,m,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_sub,m,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(s_orb,m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_phys_orb

