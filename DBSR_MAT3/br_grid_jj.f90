!======================================================================
      Subroutine br_grid
!======================================================================
!     broadcast the knot-grid data in modules DBS_grid and DBS_nuclear
!----------------------------------------------------------------------
      Use MPI
      Use DBS_grid
      Use DBS_nuclear

      Implicit none
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(grid_type,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ks, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ns, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ms, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nv, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ml, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(me, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ksp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nsp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ksq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nsq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(t)) Deallocate(t); Allocate(t(ns+ks))
      end if

      Call MPI_BCAST(t,ns+ks,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(hi,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(he,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(hmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(rmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(tmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(atomic_number,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(atomic_weight,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(nuclear,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(r_uniform, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ro_uniform,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(a_fermi,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(c_fermi,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(t_fermi,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ro_fermi,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(rrms,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(I_nuclear, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(D_nuclear, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Q_nuclear, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


      End Subroutine br_grid


