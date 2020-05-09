!=====================================================================
!     PROGRAM   b j _ m p i                    
!
!               C O P Y R I G H T -- 2009
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!
!    generates angular coefficient for Dirak_Fock calculations 
!    in case of non-orthogonal one-electron radial functions
!
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:
!    
!    klsp1,klsp2  - range of partial wave in BSR calculations,
!                   then cfg.001, cfg.002, ..., are input files
!                   (default -> 0, then input file is rcsl.inp)
!   
!    mk    - max.multipole index (default -> 9, see module param_br)
!
!    RX    - =0, all integrals, =1, only R-integrals
!----------------------------------------------------------------------
!
!    example:    1.  bjj 
!                2.  bjj klsp1=1 klsp2=5 
!                3.  bjj km=5
!            
!----------------------------------------------------------------------
!
!    INPUT FILES:
!    
!    cfg.nnn     -  configuration list for partial wave nnn (= klsp)
!                   (cfg.inp in case klsp = 0)
!                  
!    jnt_bnk.nnn -  input data bank for angular coefficients
!                   (optional; jnt_bnk in case klsp = 0)
!                   
!    
!    OUTPUT FILES:
!    
!    jnt_bnk.nnn  - output data bank for angular coefficients
!                   (jnt_bnk in case klsp = 0)
!                   
!---------------------------------------------------------------------     
      USE MPI

      USE param_jj
      USE conf_jj
      USE det_list
      USE def_list
      USE coef_list
      USE symt_list
      USE nljm_orbitals, only: mj_max,mj_orb

      Implicit none 

      Integer :: i,j, ii, iarg, fail
      Real(8) :: time, t1,t2,tt
 
      Character :: CLINE*80

      Integer, External :: IARGC, mj_value

!----------------------------------------------------------------------
! ... initialize MPI:

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if(myid.eq.0) then
       write(*,*) 'nprocs = ', nprocs
       Allocate(ip_proc(nprocs)); ip_proc=0
      end if

!----------------------------------------------------------------------
! ... read arguments from command line:

      if(myid.eq.0) Call Read_arg;   Call br_arg

!----------------------------------------------------------------------
! ... output HEADER:

      write(AF_pri,'(a,i3.3)') 'bj_',myid
      if(debug.eq.0.and.myid.gt.0) pri=0
      if(pri.gt.0) open(pri,file=AF_pri)

      if(myid.eq.0) then

      write(pri,'(/20x,a/20x,a/20x,a/)') &
             '=======================',     &
             ' B R E I T  INTERACTION',     &
             '======================='

      write(pri,'(/a,i3)') 'Max.multipole index =',mk
      write(pri,'(/a,E10.1/)') 'Tollerance for coeff.s =',eps_c
      end if

!----------------------------------------------------------------------
!                                             cycle over partial waves:
      time = 0.d0
      Do klsp = klsp1,klsp2

       if(myid.eq.0) then
        write(pri,'(80(''-''))') 
        write(pri,'(/a,i5/)') ' Partial wave: ',klsp
        write(*,'(/a,i5/)') ' Partial wave: ',klsp
       end if

       t1 = MPI_WTIME()

! ... open relevant files: 

       if(myid.eq.0) then
        fail = 0
        Call open_jj(nuc,fail)   ! c-file
        Call open_jj(nub,fail)   ! data bank results, if any
        Call open_jj(nur,fail)   ! new results
        Call open_jj(nui,fail)   ! intermediate results
        if(new) then
         write(pri,'(a/)') ' It is new calculations '
        else
         write(pri,'(a/)') ' It is continued calculations '
        end if
       end if

!       Call MPI_BCAST(fail,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!       if(myid.eq.0) write(*,*) 'fail', fail
!       if(fail.ne.0) Stop

! ...  read the configuration list:

       if(myid.eq.0) Call Read_conf_jj; 

       Call br_conf_jj

       if(pri.gt.0) &
       write(pri,*) 'icalc  ',icalc

       if(.not.icalc) then
        if(myid.eq.0) Close(nur,STATUS='DELETE')
        Cycle
       end if        

! ...  extract old results:

       if(myid.eq.0) then
       if(new) then

        ndet=0; Call Alloc_det(idet)
        ndef=0; Call Alloc_def(idef)
                Call Alloc_ndet(0)
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
        Call open_jj(nua,fail)
        Call open_jj(nud,fail)
        Call Pre_det_exp 
       end if

       if(Allocated(mj_orb)) Deallocate(mj_orb)
       Call MPI_BCAST(mj_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Allocate(mj_orb(mj_max+1))
       Do i=1,mj_max+1;  mj_orb(i)=mj_value(i);  End do

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!----------------------------------------------------------------------
! ... calculations for new angular symmetries:

      if(myid.eq.0) then

       Call Conf_loop 

      else

       Call Conf_calc

      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!----------------------------------------------------------------------
! ... record results:

       if(myid.ne.0) Cycle

       write(pri,'(/a/)') &
          ' Results for new angular symmetry calculations:'
       ii=ldet/ndet + 1
       write(pri,'(a,2i10)') &
          ' number of overlap determinants =', ndet,ii
       ii=ldef/ndef + 1
       write(pri,'(a,2i10)') &
          ' number of overlap factors      =', ndef,ii

       Call Write_symc(nur)
       Call Write_symt(nur)

       Call Write_done(nur,fail) 
       if(fail.ne.0) Cycle

       Call Record_det(nur)
       Call Record_def(nur)

!       rewind(nui);   Call RW(nui,nur,nct)

       close(nui);  close(nur); close(nub)

       CLINE = ' '
       if(klsp.eq.0) &
       write(CLINE,'(a,a,a,a)') 'cat ',trim(AF_i),' >> ',trim(AF_r)
       if(klsp.gt.0) &
       write(CLINE,'(a,a,a,a)') 'cat ',trim(AF_i),' >> ',trim(BF_r)

       Call System(CLINE)

! ...  rename new results as new data bank (int_res -> int_bnk): 
 

       CLINE = 'mv ';   i = 3               !  UNIX

       if(klsp.eq.0) then
        ii = LEN_TRIM(AF_r); CLINE(i+1:i+ii)=AF_r; i=i+ii+1
        ii = LEN_TRIM(AF_b); CLINE(i+1:i+ii)=AF_b; i=i+ii
       else
        ii = LEN_TRIM(BF_r); CLINE(i+1:i+ii)=BF_r; i=i+ii+1
        ii = LEN_TRIM(BF_b); CLINE(i+1:i+ii)=BF_b; i=i+ii
       end if

       Call System(CLINE(1:i))

       write(pri,'(a,i10,f10.1)') ' total number of coeff.s        =', nct

! ... time for one partial wave:

       t2=MPI_WTIME(); tt=(t2-t1)/60
       write(pri,'(/a,F12.2,a)') ' Partial wave:',tt,' min'
       write(*,  '( a,F12.2,a)') ' Partial wave:',tt,' min'
       time = time + tt

       CLINE = 'rm '//trim(AF_i);   Call System(trim(CLINE))

       if(klsp.eq.0) then
        CLINE = 'rm '//trim(AF_r);   Call System(trim(CLINE))
       else
        CLINE = 'rm '//trim(BF_r);   Call System(trim(CLINE))
       end if

      End do  ! over klsp
!----------------------------------------------------------------------

! ... total time:
      if(myid.eq.0) &
      write(pri,'(a,F12.2,a)') ' Total time: ',time,' min'
      if(myid.eq.0) &
      write(*,  '(a,F12.2,a)') ' Total time: ',time,' min'

      Call MPI_FINALIZE(ii)

      END ! Program bj_mpi


!======================================================================
      Subroutine Read_arg 
!======================================================================
!     read arguments from command line and check default settings
!----------------------------------------------------------------------

      Use param_jj 

      Implicit none

      Call Read_iarg('klsp'  ,klsp )
      Call Read_iarg('klsp1' ,klsp1 )
      Call Read_iarg('klsp2' ,klsp2 )
      Call Read_iarg('mk'    ,mk    )
      Call Read_rarg('eps_c' ,eps_c )
      Call Read_rarg('RX'    ,RX    )
      Call Read_iarg('debug' ,debug )

      if(klsp.ne.0) then
       klsp1=klsp
       klsp2=klsp
      end if
      if(klsp2.lt.klsp1) klsp2=klsp1 
     
      End Subroutine Read_arg

!======================================================================
      Subroutine br_arg 
!======================================================================
!     brodcast main arguments
!----------------------------------------------------------------------

      Use mpi
      Use param_jj 

      Implicit none

      Call MPI_BCAST(klsp, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(klsp1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(klsp2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mk,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(RX,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(eps_c,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_arg


!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------

      Implicit none

      Integer(4), Intent(in) :: nu1,nu2
      Integer(4), Intent(out) :: nc
      Integer(4) :: i,j

      Integer(4),Parameter :: mc = 100000
      Integer(4),Allocatable,Dimension(:) :: K1,K2,K3,K4
      Real(8),Allocatable,Dimension(:) :: C

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
