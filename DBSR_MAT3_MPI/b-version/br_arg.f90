!======================================================================
      Subroutine br_arg
!======================================================================
!     broadcast the main parameters used in dbsr_mat program
!----------------------------------------------------------------------
      Use MPI
      Use dbsr_mat
      Use DBS_dhl_pq

      Implicit none
      
      Call MPI_BCAST(klsp1, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(klsp2, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mk,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mblock,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nblock,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kblock,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nud,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iitar, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ilzero,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jlzero,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ibzero,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jbzero,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(mbreit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mbloch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(debug, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mcbuf, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(pri_acf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(check_target,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(eps_c,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(eps_v,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(eps_det,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(EC,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Ecore,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(RB,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(pnu,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(s_ovl,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(Edmin,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Edmax,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Egap,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_arg



