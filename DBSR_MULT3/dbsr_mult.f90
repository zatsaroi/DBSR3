!=====================================================================
!     PROGRAM       D B S R _ M U L T              
!
!               C O P Y R I G H T -- 2007
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     This program evaluates the MULTIPOLE operators in jj-coupling
!     including the case of non-orthogonal orbitals
!     
!     The technique is close to LS-case described in:
!
!     Comp.Phys.Commun.  98 (1996) 235-254
!     Comp.Phys.Commun. 174 (2006) 273-356
!======================================================================
!
!    INPUT ARGUMENTS:
!    
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     AA   -  type of calculation, E1,E2,M1,M2,...
! 
!    INPUT FILES:
!
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     mult_bnk_E1(M1,...) - data-bank for angular coefficients (optional)
!
!     OUTPUT FILES:
!
!     mult_bnk_E1(M1,...) - new data bank for angular coefficients
!     dbsr_mult.log  -  running information
!
!-----------------------------------------------------------------------
!     ONLY ONE TYPE OF MULTIPOLE INDEX IN ONE RUN !!!
!-----------------------------------------------------------------------

      Use dbsr_mult;  Use det_list;  Use def_list

      Implicit none 
      Integer :: i,ii,nct
      Real(8) :: t1,t2,tt 
      Real(8), external :: RRTC
 
! ... read arguments from command line:

      Call Read_arg

!-----------------------------------------------------------------------
!                                                          calculations:
      t1=RRTC()

! ... input of configurations:

      Call Check_mult_bnk
      if(.not.icalc) Stop 'no new calculations needed'

! ...  extract old results:

      if(new) then
        ndet=0; Call Alloc_det(idet); Call Alloc_def(idef)
      else
        Call Read_det(nub)
        Call Read_def(nub)
        Call RW(nub,nui,nct)
        write(pri,'(/a/)') &
          ' Results for the old symmetry calculations:'
        write(pri,'(/3(a,i8))') &
          ' ndet =',ndet,'   ndef  =',ndef,'  ncoef =',nct
      end if

! ... prepare det. expantions:

      Call Pre_det_exp 

! ... calculations for new angular symmetries:

      Call Conf_loop 

! ... record results:

      write(pri,'(/a/)') &
         ' Results for new angular symmetry calculations:'
      ii=0; if(ndet.gt.0) ii=ldet/ndet + 1
      write(pri,'(a,2i10)') &
         ' number of overlap determinants =', ndet,ii
      ii=0; if(ndef.gt.0) ii=ldef/ndef + 1
      write(pri,'(a,2i10)') &
         ' number of overlap factors      =', ndef,ii

       write(nur) ktype,kpol
       Call Write_symc(nur)
       Call Write_symt(nur)
       Call Write_done(nur) 
       Call Record_det(nur)
       Call Record_def(nur)
       rewind(nui); Call RW(nui,nur,nct)
       close(nui); close(nur); close(nub)

! ...  rename new results as new data bank (int_res -> jnt_bnk): 
 
       AF = 'move '; i = 5               !  Windows
!       AF = 'mv ';   i = 3               !  UNIX

       ii = LEN_TRIM(AF_res); AF(i+1:i+ii)=AF_res; i=i+ii+1
       ii = LEN_TRIM(AF_bnk); AF(i+1:i+ii)=AF_bnk; i=i+ii

       Call System(AF(1:i))

       write(pri,'(a,i10,f10.1)') &
          ' total number of coeff.s        =', nct

! ... total time:

      t2=RRTC(); tt=(t2-t1)/60

      write(pri,'(/a,F12.2,a)') ' DBSR_mult time: ',tt,' min'
      write(*,  '(/a,F12.2,a)') ' DBSR_mult time: ',tt,' min'

      End ! Program DBSR_MULT


!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu1,nu2
      Integer, intent(out) :: nc
      Integer :: i,j
      Integer, parameter :: mc = 100000
      Integer, allocatable :: K1(:),K2(:),K3(:),K4(:)
      Real(8), allocatable :: C(:)

      Allocate(C(mc),K1(mc),K2(mc),K3(mc),K4(mc))

      nc = 0
      i = 1
    1 read(nu1,end=2) C(i),k1(i),k2(i),k3(i),k4(i)
      i = i + 1; if(i.le.mc) go to 1
    2 j = i - 1
      nc = nc + j

      Do i = 1,j; write(nu2) C(i),k1(i),k2(i),k3(i),k4(i); End do

      i = 1;  if(j.eq.mc) go to 1

      Deallocate(C,K1,K2,K3,K4)

      End Subroutine RW
