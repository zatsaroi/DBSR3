!======================================================================
      Subroutine br_conf_jj
!======================================================================
!     broadcast the data in modules "orb_jj, symc_list, symt_list, conf_jj"
!----------------------------------------------------------------------
      Use MPI
      Use conf_jj;  Use symc_list;  Use symt_list

      Implicit none
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

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
      Call MPI_BCAST(ip_conf,msymc,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
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

      Call MPI_BCAST(IT_conf,msymt,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_term,msymt,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(JS_term,ksymt,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(VS_term,ksymt,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(JI_term,ksymt,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)

! ... conf_jj list:

      Call MPI_BCAST(ncfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mcfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kcfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lcfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(IP_state )) Deallocate(IP_state ); Allocate(IP_state (mcfg ))
       if(allocated(IS_term  )) Deallocate(IS_term  ); Allocate(IS_term  (mcfg ))
       if(allocated(IP_orb   )) Deallocate(IP_orb   ); Allocate(IP_orb   (kcfg ))
       if(allocated(IS_order )) Deallocate(IS_order ); Allocate(IS_order (ncfg ))
       if(allocated(IT_state1)) Deallocate(IT_state1); Allocate(IT_state1(nsymt))
       if(allocated(IT_state2)) Deallocate(IT_state2); Allocate(IT_state2(nsymt))
       if(allocated(WC)) Deallocate(WC); Allocate(WC(ncfg))
      end if

      Call MPI_BCAST(IP_state ,mcfg ,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IS_term  ,mcfg ,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IP_orb   ,kcfg ,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IS_order ,ncfg ,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IT_state1,nsymt,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IT_state2,nsymt,MPI_INTEGER ,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(WC,ncfg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      memory_conf_jj = 4.d0 * (3*msymc + 2*ksymc + 2*msymt + 3*ksymt + &
                               2*mcfg + kcfg + 3*ncfg + 2*nsymt) / (1024*1024)

      End Subroutine br_conf_jj


