!=====================================================================
!     PROGRAM   D B S R _ B R E I T  
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
      Integer :: klsp
      Integer :: i, ii
      Real(8) :: t1,t2
      Real(8), external :: RRTC

      Character*80 :: cline

      Call Inf_dbsr_breit

!----------------------------------------------------------------------
! ... HEADER:

      Call open_jj(pri,0)
      write(pri,'(a/a/a)') &
        '==================================',     &
        ' DBSR_BREIT: angular coefficients ',     &
        '=================================='

! ... read arguments if any from command line:

      Call Read_arg

!----------------------------------------------------------------------
! ... cycle over cases (c-file or partial waves):

      Do klsp = klsp1,klsp2

       if(klsp.gt.0) then
        write(pri,'(80(''-''))')
        write(pri,'(/a,i5)') 'Partial wave: ',klsp
        write(*  ,'(/a,i5)') 'Partial wave: ',klsp
       end if

       t1 = RRTC()

! ...  open relevant files:

       Call open_jj(nuc,klsp)       ! c-file
       Call open_jj(nub,klsp)       ! data bank results, if any
       Call open_jj(nur,klsp)       ! new results
       Call open_jj(nui,klsp)       ! intermediate results

       if(new)      write(pri,'(/a/)') 'It is a new calculation '
       if(.NOT.new) write(pri,'(/a/)') 'It is a continued calculation '

       if(mbreit.eq.1) write(pri,'(a/)') 'Breit contribution is included '    
       if(mbreit.eq.0) write(pri,'(a/)') 'Breit contribution is not included '    

! ...  read the configuration list:
                           
       Call Read_conf_jj   

       if(.not.icalc) then; Close(nur,STATUS='DELETE'); Cycle; end if

! ...  extract old results if any:

       if(new) then
        ndet=0; Call Alloc_det(idet)
        ndef=0; Call Alloc_def(idef)
                Call Alloc_ndet(0)
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

       Call Conf_loop

! ...  record results:

       write(pri,'(/a/)') &
          'Results for new angular symmetry calculations:'
       ii=ldet/ndet + 1
       write(pri,'(a,2i10)') &
          'number of overlap determinants =', ndet,ii
       ii=ldef/ndef + 1
       write(pri,'(a,2i10)') &
          'number of overlap factors      =', ndef,ii

       Call Write_symc(nur)
       Call Write_symt(nur)
       Call Write_done(nur)
       Call Record_det(nur)
       Call Record_def(nur)

       rewind(nui);  Call RW(nui,nur,nct)

       write(pri,'(a,i10,f10.1)') 'total number of coeff.s        =', nct

       close(nui); close(nur); close(nub)

! ...  rename results as new data bank (jnt_res -> jnt_bnk):

       if(klsp.eq.0) then
        write(cline,*) 'mv ',trim(AF_r),' ',trim(AF_b)
       else
        write(cline,*) 'mv ',trim(BF_r),' ',trim(BF_b)
       end if

       Call System(cline)

! ...  time for one case:

       t2=RRTC()
       write(pri,'(/a,F12.2,a)') 'time: ',(t2-t1)/60,' min'
       write(*,  '( a,F12.2,a)') 'time: ',(t2-t1)/60,' min'

      End do  ! over klsp

      End ! Program dbsr_breit


!======================================================================
      Subroutine Read_arg
!======================================================================
!     read arguments from command line and check default settings
!======================================================================
      Use dbsr_breit
      Implicit none
      Integer :: klsp = 0

! ... read arguments in command line:

      Call Read_iarg('klsp'  ,klsp  )
      Call Read_iarg('klsp1' ,klsp1 )
      Call Read_iarg('klsp2' ,klsp2 )
      Call Read_iarg('mk'    ,mk    )
      Call Read_rarg('eps_c' ,eps_c )
      Call Read_iarg('mbreit',mbreit)

      if(klsp.gt.0) klsp1=klsp
      if(klsp2.lt.klsp1) klsp2=klsp1

      write(pri,'(/a,i3)') 'Max. multipole index =',mk

      write(pri,'(/a,1PE10.0)') 'Tollerance for coefficients =',eps_c

      if(mbreit.eq.0)  write(pri,'(/a)') 'Breit coefficients are not included'
      if(mbreit.eq.1)  write(pri,'(/a)') 'Breit coefficients are included'

      End Subroutine Read_arg


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


