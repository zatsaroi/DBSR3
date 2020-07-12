!=====================================================================
!     PROGRAM   D B S R _ B R E I T _ M P I                 v.3
!
!               C O P Y R I G H T -- 2020
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    generates angular coefficient for Dirak_Fock calculations 
!    in case of non-orthogonal one-electron radial functions
!----------------------------------------------------------------------
!     INPUT ARGUMENTS:
!
!     c-file    name.c or cfg.inp by default;
!               DBSR option:  cfg.001, cfg.002, ..., cfg.nnn
!               if range of partial waves are given through
!               arguments klsp  or  klsp1, klsp2
!
!     mk     -  max.multipole index (default mk=9)
!     eps_c  -  tollerence for angular coefficients
!     mbreit -  by default = 0 -> only Coulomb interaction
!----------------------------------------------------------------------
!     example:    1.  dbsr_breit 3p5.c
!                 2.  dbsr_breit klsp1=1 klsp2=5
!                 3.  dbsr_breit km=5
!----------------------------------------------------------------------
!     INPUT FILES:   c-files with list of atomic states in GRASP format:
!                    name.c or rcsl.inp
!                    
!     OUTPUT FILES:  data bank for angular coefficients:
!                    int_bnk (or name.bnk, or int_bnk.nnn)
!----------------------------------------------------------------------     
      Use MPI

      Use dbsr_breit
      Use symt_list;  Use conf_jj;  Use det_list;  Use def_list
      Use coef_list
      Use nljm_orbitals, only: mj_max,mj_orb

      Implicit none 
      Character(80) :: cline
      Integer :: i,j, ii, fail, klsp
      Real(8) :: time, t1,t2,tt
      Integer, external ::  mj_value

      Call Inf_dbsr_breit

!----------------------------------------------------------------------
! ... initialize MPI:

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if(myid.eq.0)  Allocate(ip_proc(nprocs))

!----------------------------------------------------------------------
! ... read arguments from command line:

      if(myid.eq.0) Call Read_arg

      Call br_arg

      if(myid.gt.0) then
       write(AF,'(a,i4.4)') 'debug_',myid
       if(debug.eq.0) pri=0
       if(pri.gt.0) open(pri,file=AF)
      end if

!----------------------------------------------------------------------
!                                             cycle over partial waves:
      time = 0.d0
      Do klsp = klsp1,klsp2

       t1 = MPI_WTIME()

       if(myid.eq.0) then

! ... open relevant files: 

        Call open_jj(pri,klsp)
        Call open_jj(nuc,klsp)   ! c-file
        Call open_jj(nub,klsp)   ! data bank results, if any

! ... print parameters:

       write(pri,'(a/a/a)') &
        '==============================================================', &
        ' DBSR_BREIT: angular coefficients for Dirak-Breit Hamiltonian ', &
        '=============================================================='    

        write(pri,'(/a,i5)') 'nprocs = ', nprocs
        ip_proc=0

        if(klsp.gt.0) then
         write(pri,'(/a,i5)') 'Partial wave: ',klsp
         write(  *,'(/a,i5)') 'Partial wave: ',klsp
        else
         write(pri,'(/a,a )') 'Name of case: ',trim(name)
         write(  *,'(/a,a )') 'Name of case: ',trim(name)
        end if

        write(pri,'(/a,i3)') 'Max. multipole index =',mk
        write(pri,'(a,1PE8.0)') 'Tollerance for coefficients =',eps_c
        if(mbreit.eq.0)  write(pri,'(a)') 'Breit coefficients are not included'
        if(mbreit.eq.1)  write(pri,'(a)') 'Breit coefficients are included'

        if(new) then
         write(pri,'(/a)') 'It is new calculations '
        else
         write(pri,'(/a)') 'It is continued calculations '
        end if

! ...  read the configuration list:

        Call Read_conf_jj 

        if(icalc .and. .not.new) then
         write(pri,'(/a,a4/)')   'Need of additional calculations --> yes '
         write(  *,'(/a,a4/)')   'Need of additional calculations --> yes '
        elseif(.not.icalc .and. .not.new) then
         write(pri,'(/a,a4/)')   'Need of additional calculations --> no '
         write(  *,'(/a,a4/)')   'Need of additional calculations --> no '
        end if

       end if   !  myid

       Call MPI_BCAST(icalc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       if(.not.icalc) Cycle

       Call br_conf_jj

! ...  extract old results:

       if(myid.eq.0) then

        Call open_jj(nur,klsp)   ! new results
        Call open_jj(nui,klsp)   ! intermediate results

        nct = 0
        if(new) then
         Call Alloc_det(-1)
         Call Alloc_def(-1)
        else
         Call Read_det(nub)
         Call Read_def(nub)
         Call RW(nub,nui,nct)
        end if

       end if

!----------------------------------------------------------------------
! ... prepare det. expantions:

       Call Alloc_nljm(ne,msh)

       if(myid.eq.0) then
        Call open_jj(nua,klsp)
        Call open_jj(nud,klsp)
        Call Pre_det_exp 
       end if

       if(Allocated(mj_orb)) Deallocate(mj_orb)
       Call MPI_BCAST(mj_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Allocate(mj_orb(mj_max+1))
       Do i=1,mj_max+1;  mj_orb(i)=mj_value(i);  End do

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!======================================================================
! ... calculations for new angular symmetries:

      nc = 0

      if(myid.eq.0) then
       Call Conf_loop 
      else
       Call Conf_calc
      end if

      nct = nct + nc

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!======================================================================
! ... record results:

       if(myid.ne.0) Cycle

       Call Write_symc(nur)
       Call Write_symt(nur)
       Call Write_done(nur) 

       Call Record_det(nur)
       Call Record_def(nur)

       Close(nui);  Close(nur); Close(nub)

       write(cline,'(a,a,a,a)') 'cat ',trim(BF_i),' >> ',trim(BF_r)

       Call System(cline)

! ...  rename new results as new data bank (int_res -> int_bnk): 
 
       write(cline,'(a,a,a,a)') 'mv ',trim(BF_r),'  ',trim(BF_b)

       Call System(cline)

       write(pri,'(/a/)')   'Data-bank contains:'
       write(pri,'(a,i10)') 'number of overlap determinants =', ndet
       write(pri,'(a,i10)') 'number of overlap factors      =', ndef
       write(pri,'(a,i10)') 'total number of coeff.s        =', nct

! ... time for one partial wave:

       t2=MPI_WTIME(); tt=(t2-t1)/60
       write(pri,'(/a,T15,F12.2,a)') 'Partial wave:',tt,' min'
       write(*,  '(/a,T15,F12.2,a)') 'Partial wave:',tt,' min'
       time = time + tt

       write(cline,'(a,a,a,a)') 'rm ',trim(BF_i);  Call System(cline)

      End do  ! over klsp

!----------------------------------------------------------------------
! ... total time:

      if(myid.eq.0) then
       write(pri,'(/a,T15,F12.2,a)') 'Total time: ',time,' min'
       write(*,  '(/a,T15,F12.2,a)') 'Total time: ',time,' min'
      end if

      Call MPI_FINALIZE(ii)

      End   ! Program dbsr_breit_mpi


!======================================================================
      Subroutine br_arg 
!======================================================================
!     brodcast main arguments
!----------------------------------------------------------------------
      Use mpi
      Use dbsr_breit 

      Implicit none

      Call MPI_BCAST(klsp1, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(klsp2, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mk,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mbreit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(eps_c, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_arg


!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------
      Implicit none

      Integer, intent(in) :: nu1,nu2
      Integer, intent(out) :: nc
      Integer :: i,j

      Integer, parameter :: mc = 10000000
      Integer, allocatable :: K1(:),K2(:),K3(:),K4(:)
      Real(8), allocatable :: C(:)

      Allocate(C(mc),K1(mc),K2(mc),K3(mc),K4(mc))

      nc = 0
      i = 1
    1 read(nu1,end=2) C(i),k1(i),k2(i),k3(i),k4(i)
      i = i + 1; if(i.le.mc) go to 1
    2 j = i - 1
      nc = nc + j

      Do i = 1,j
       write(nu2) C(i),k1(i),k2(i),k3(i),k4(i)
      End do

      i = 1;  if(j.eq.mc) go to 1

      Deallocate(C,K1,K2,K3,K4)

      End Subroutine RW
