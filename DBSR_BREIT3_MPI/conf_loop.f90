!======================================================================
      Subroutine Conf_loop
!======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------
      Use dbsr_breit
      Use term_exp
      Use conf_jj
      Use symc_list

      Implicit none 
      Integer :: i,j,ij,is,js, met
      Real(8) :: t1,t2,tt   

      Call CPU_time(t1)

!----------------------------------------------------------------------
!                                          cycle 1 over configurations:
      rewind(nud)

      Do is=1,ic_case

       Read(nud) ic,kt1,kdt1

       if(Allocated(IP_kt1)) Deallocate(IP_kt1)
       Allocate(IP_kt1(kt1)); Read(nud) IP_kt1

       if(Allocated(IP_det1)) Deallocate(IP_det1)
       Allocate(IP_det1(ne,kdt1)); Read(nud) IP_det1
             
       if(Allocated(C_det1)) Deallocate(C_det1)
       Allocate(C_det1(kt1,kdt1)); Read(nud) C_det1

!----------------------------------------------------------------------
!                                          cycle 2 over configurations:
      rewind(nud)
      Do js=1,is

       Read(nud) jc,kt2,kdt2

       if(Allocated(IP_kt2)) Deallocate(IP_kt2)
       Allocate(IP_kt2(kt2)); Read(nud) IP_kt2

       if(Allocated(IP_det2)) Deallocate(IP_det2)
       Allocate(IP_det2(ne,kdt2)); Read(nud) IP_det2
              
       if(Allocated(C_det2)) Deallocate(C_det2)
       Allocate(C_det2(kt2,kdt2)); Read(nud) C_det2

       i = max(ic,jc); j = min(ic,jc); ij = i*(i-1)/2 + j
       if(JC_need(ij).eq.0) Cycle      
              
!----------------------------------------------------------------------

       met = 0
       Do i=1,nprocs-1
        if(ip_proc(i).ne.0) Cycle 
        Call Send_det_exp(i)
        met = i 
        ip_proc(i) = 1
        Exit
       End do

       if(met.eq.0) then
        Call Get_res(i)        
        Call Add_res(nui)
        Call Send_det_exp(i)
        met = i 
       end if

!----------------------------------------------------------------------

      End do    ! over jc

      Call CPU_time(t2)

      Call Incode_confj1
      write(*,  '(a,4i6,f8.0,a,a)')  &
        'ic,nsymc,kt,kdt =',ic,nsymc,kt1,kdt1,t2-t1,' sec ',trim(CONFIG)
      write(pri,'(a,4i6,f8.0,a,a)')  &
        'ic,nsymc,kt,kdt =',ic,nsymc,kt1,kdt1,t2-t1,' sec ',trim(CONFIG)

      End do    ! over ic

!----------------------------------------------------------------------
! ... final results:

       Do 
        if(sum(ip_proc).eq.0) Exit 
        Call Get_res(j)        
        Call Add_res(nui)
        ip_proc(j) = 0
       End do

! ... release the nodes:

       ic = -1
       Do i=1,nprocs-1
        Call Send_det_exp(i)
       End do

      End Subroutine Conf_loop


!======================================================================
      Subroutine send_det_exp(id)
!======================================================================
!  ... send determinant expansion to processor 'id'
!----------------------------------------------------------------------
      Use MPI
      Use conf_jj,   only: ne
      Use term_exp

      Implicit none
      Integer :: ierr, id 

      Call MPI_SEND(ic  ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      if(ic.le.0) Return

      Call MPI_SEND(kt1 ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kdt1,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(ip_kt1,kt1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ip_det1,ne*kdt1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(C_det1,kt1*kdt1,MPI_DOUBLE_PRECISION,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(jc  ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(kt2 ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kdt2,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(ip_kt2,kt2,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ip_det2,ne*kdt2,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(C_det2,kt2*kdt2,MPI_DOUBLE_PRECISION,id,0,MPI_COMM_WORLD,ierr)

      End Subroutine send_det_exp


