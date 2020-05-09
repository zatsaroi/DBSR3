!======================================================================
      Subroutine br_conf_jj
!======================================================================

      USE MPI
      USE conf_jj
      USE symc_list
      USE symt_list
      USE param_jj, only: new, icalc

      Implicit none
      
      Integer :: myid,ierr, nc,nt

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! ... symc_list:

      Call MPI_BCAST(nsymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(msymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ksymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lsymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(JT_conf)) Deallocate(JT_conf); Allocate(JT_conf(msymc))
       if(allocated(no_conf)) Deallocate(no_conf); Allocate(no_conf(msymc))
       if(allocated(ip_conf)) Deallocate(ip_conf); Allocate(ip_conf(msymc))
       if(allocated(iq_conf)) Deallocate(iq_conf); Allocate(iq_conf(ksymc))
       if(allocated(kn_conf)) Deallocate(kn_conf); Allocate(kn_conf(ksymc))
      end if

      Call MPI_BCAST(JT_conf,msymc,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(no_conf,msymc,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_conf,msymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iq_conf,ksymc,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kn_conf,ksymc,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)

! ... symt_list:

      Call MPI_BCAST(nsymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(msymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ksymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lsymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    
      if(myid.ne.0) then
       if(allocated(IT_conf)) Deallocate(IT_conf); Allocate(IT_conf(msymt))
       if(allocated(ip_term)) Deallocate(ip_term); Allocate(ip_term(msymt))
       if(allocated(JS_term)) Deallocate(JS_term); Allocate(JS_term(ksymt))
       if(allocated(VS_term)) Deallocate(VS_term); Allocate(VS_term(ksymt))
       if(allocated(JI_term)) Deallocate(JI_term); Allocate(JI_term(ksymt))
      end if

      Call MPI_BCAST(IT_conf,msymt,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_term,msymt,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(JS_term,ksymt,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(VS_term,ksymt,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(JI_term,ksymt,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)

! ... conf_jj list:

      Call MPI_BCAST(ncfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mcfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kcfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lcfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      nc = nsymc*(nsymc+1)/2
      nt = nsymt*(nsymt+1)/2
      if(myid.ne.0) then
      if(allocated(IC_need)) Deallocate(IC_need); Allocate(IC_need(nsymc))
      if(allocated(JC_need)) Deallocate(JC_need); Allocate(JC_need(nc   ))
      if(allocated(IT_need)) Deallocate(IT_need); Allocate(IT_need(nsymt))
      if(allocated(IT_done)) Deallocate(IT_done); Allocate(IT_done(nt   ))
      end if

      Call MPI_BCAST(IC_need,nsymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(JC_need,nc   ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IT_done,nt   ,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IT_need,nsymt,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(new  ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(icalc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_conf_jj


