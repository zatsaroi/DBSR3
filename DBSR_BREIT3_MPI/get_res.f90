!======================================================================
      Subroutine Get_res(id)
!======================================================================
      Use MPI

      Use NDET_list
      Use NDEF_list 
      Use coef_list,  only: ntrm,ncoef,idfc,intc,coef
      Use term_exp,   only: jt1,jt2,JP_kt1,JP_kt2, IP_kt12
      Use dbsr_breit, only: pri, AF_pri

      Implicit none
      Integer :: id, itag, ierr, status(MPI_STATUS_SIZE)

      Call MPI_RECV(ncoef,1,MPI_INTEGER, MPI_ANY_SOURCE, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      id = status(MPI_SOURCE)

! ... term indexes:

      Call MPI_RECV(jt1,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(jt2,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(allocated(JP_kt1)) Deallocate(JP_kt1); Allocate(JP_kt1(jt1))
      if(allocated(JP_kt2)) Deallocate(JP_kt2); Allocate(JP_kt2(jt2))
      Call MPI_RECV(JP_kt1,jt1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(JP_kt2,jt2,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      if(allocated(IP_kt12)) Deallocate(IP_kt12); Allocate(IP_kt12(jt1,jt2))

      Call MPI_RECV(IP_kt12,jt1*jt2, MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)


      if(ncoef.eq.0) Return

      Call MPI_RECV(ntrm,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(allocated(coef)) Deallocate(coef); Allocate(coef(ntrm,ncoef))
      if(allocated(idfc)) Deallocate(idfc); Allocate(idfc(ncoef))
      if(allocated(intc)) Deallocate(intc); Allocate(intc(ncoef))

      Call MPI_RECV(idfc,ncoef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(intc,ncoef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(coef,ntrm*ncoef,MPI_DOUBLE_PRECISION, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

! ... dets:

      Call MPI_RECV(ndet,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(ndet.eq.0) Return

      Call MPI_RECV(ldet,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(allocated(KPD)) Deallocate(KPD); Allocate(KPD(ndet))
      if(allocated(IPD)) Deallocate(IPD); Allocate(IPD(ndet))
      mdet=ndet
      if(allocated(NPD)) Deallocate(NPD); Allocate(NPD(ldet))
      kdet=ldet

      Call MPI_RECV(KPD,ndet,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(IPD,ndet,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(NPD,ldet,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

! ... defs:

      Call MPI_RECV(ndef,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      Call MPI_RECV(ldef,1,MPI_INTEGER, id, &
                     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      if(allocated(KPf)) Deallocate(KPf); Allocate(KPf(ndef))
      if(allocated(IPf)) Deallocate(IPf); Allocate(IPf(ndef))
      mdef=ndef
      if(allocated(NPf)) Deallocate(NPf); Allocate(NPf(ldef))
      kdef=ldef

      Call MPI_RECV(KPf,ndef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(IPf,ndef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      Call MPI_RECV(NPf,ldef,MPI_INTEGER, id, &
                     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)


      End Subroutine get_res
