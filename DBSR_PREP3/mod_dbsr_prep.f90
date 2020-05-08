!======================================================================
      Module dbsr_prep
!======================================================================
!     contains main parameters and arrays used in program dbsr_prep
!----------------------------------------------------------------------
      Use conf_jj;          Use DBS_grid
      Use orb_jj;           Use DBS_orbitals_pq
      Use target_jj;        Use DBS_gauss
      Use channels_jj, npert1 => npert, ipert1 => ipert, ippert1 => ippert       
      Use channel_jj,  only: npert, ipert, ippert

      Implicit none

! ... files:

      Integer, parameter :: ma=80

      Character(ma) :: AF, atar, name = ' '

      Integer :: nut = 21;  Character(ma) :: AF_tar  = 'target_jj'
      Integer :: nuo = 22;  Character(ma) :: AF_orb  = 'target_orb'
      Integer :: nup = 23;  Character(ma) :: AF_par  = 'dbsr_par'
      Integer :: pri = 66;  Character(ma) :: AF_log  = 'dbsr_prep.log'
                            Character(ma) :: AF_sub  = 'target_sub.bsw'
                            Character(ma) :: AF_wfn  = 'target.bsw'
                            Character(ma) :: AF_knot = 'knot.dat'
      Integer :: nuc = 11;  Character(ma) :: AFC  ! input  c-file
      Integer :: nuw = 12;  Character(ma) :: AFW  ! input  w-file
      Integer :: muc = 13;  Character(ma) :: BFC  ! output c-file
      Integer :: muw = 14;  Character(ma) :: BFW  ! output w-file

! ... default values:

      Real(8) :: eps_ovl  = 1.d-7
      Real(8) :: eps_core = 1.d-5
      Real(8) :: eps_sub  = 0.50
      Real(8) :: eps_phys = 0.25
      Real(8) :: eps_conf = 0.50
      Real(8) :: eps_targ = 2.d-8

      Integer :: JJ_min = -1
      Integer :: JJ_max = -1

      Integer, external :: Icheck_file
      Real(8), external :: QUADR_PQ, QUADR_00, OBS

      Integer, allocatable :: ipt(:)

      End Module dbsr_prep 

