!=====================================================================
!     PROGRAM   D B S R _ B R E I T                           v.3
!
!               C O P Y R I G H T -- 2015
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     generates angular coefficient for Dirak_Fock calculations
!     in case of non-orthogonal one-electron radial functions
!----------------------------------------------------------------------
!     INPUT ARGUMENTS:
!
!     c-file    name.c or rcsl.inp by default;
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
!---------------------------------------------------------------------
      Use dbsr_breit
      Use symt_list; Use conf_jj; Use det_list; Use def_list
      Use coef_list

      Implicit none
      Integer :: klsp, ii
      Real(8) :: t1,t2
      Character(80) :: cline

      Call Inf_dbsr_breit

!----------------------------------------------------------------------
! ... read arguments if any from command line:

      Call Read_arg

!----------------------------------------------------------------------
! ... cycle over cases (c-file or partial waves):

      Do klsp = klsp1,klsp2

       Call CPU_time(t1)

! ...  open relevant files:

       Call open_jj(pri,klsp)
       Call open_jj(nuc,klsp)       ! c-file
       Call open_jj(nub,klsp)       ! data bank results, if any

! ... print parameters:

       write(pri,'(a/a/a)') &
        '==============================================================',     &
        ' DBSR_BREIT: angular coefficients for Dirak-Breit Hamiltonian ',     &
        '=============================================================='

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

       if(.not.icalc) Cycle

! ...  extract old results if any:

       Call open_jj(nur,klsp)       ! new results
       Call open_jj(nui,klsp)       ! intermediate results

       nct = 0
       if(new) then
        Call Alloc_det(-1)
        Call Alloc_def(-1)
       else
        Call Read_det(nub)
        Call Read_def(nub)
        Call RW(nub,nui,nct)
       end if

! ...  prepare determinant expantions:

       Call open_jj(nua,0)
       Call open_jj(nud,0)
       Call Pre_det_exp

! ...  calculations for new angular symmetries:

       nc=0;  Call Conf_loop;  nct = nct + nc

! ...  record results:

       Call Write_symc(nur)
       Call Write_symt(nur)
       Call Write_done(nur)
       Call Record_det(nur)
       Call Record_def(nur)

       Close(nui); Close(nur); Close(nub)

       write(cline,'(a,a,a,a)') 'cat ',trim(BF_i),' >> ',trim(BF_r)
       Call System(cline)

       write(pri,'(/a/)')   'Data-bank contains:'
       write(pri,'(a,i10)') 'number of overlap determinants =', ndet
       write(pri,'(a,i10)') 'number of overlap factors      =', ndef
       write(pri,'(a,i10)') 'total number of coeff.s        =', nct

! ...  rename results as new data bank (int_res -> int_bnk):

       write(cline,'(a,a,a,a)') 'mv ',trim(BF_r),'  ',trim(BF_b)
       Call System(cline)

       write(cline,'(a,a,a,a)') 'rm ',trim(BF_i);  Call System(cline)

! ...  time for one case:

       Call CPU_time(t2)
       write(pri,'(/a,T15,F12.2,a)') 'time: ',(t2-t1)/60,' min'
       write(*,  '( a,T15,F12.2,a)') 'time: ',(t2-t1)/60,' min'

      End do  ! over klsp

      End ! Program dbsr_breit


!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu1,nu2
      Integer, intent(out) :: nc
      Integer :: i,j
      Integer, parameter :: mc = 1000000
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


