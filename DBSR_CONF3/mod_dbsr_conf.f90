!======================================================================
      Module dbsr_conf
!======================================================================
!     contains main parameters and arrays used in program dbsr_conf
!----------------------------------------------------------------------
      Use target_jj;    Use channel_jj;    Use conf_jj
      Use orb_jj;       Use phys_orb_jj

      Implicit none

      Integer :: myid = 0, ierr = 0, nprocs = 1

! ... files:

      Integer, parameter :: ma=80
      Integer :: nut = 1;  Character(ma) :: AF_tar ='target_jj'
      Integer :: nup = 2;  Character(ma) :: AF_par ='dbsr_par'
      Integer :: pri =66;  Character(ma) :: AF_log ='dbsr_conf.log'
      Integer :: nuo = 7;  Character(ma) :: AF_orb ='target_orb'
      Integer :: nuc = 8;  Character(ma) :: AF  ! c-files
      Integer :: nuw = 9;  Character(ma) :: AF_wfn ='target.bsw'

! ... default values for input parameters:    

! ... values used in different subroutions:

      Integer :: ilsp
      Integer :: ncfg_phys, ncfg_targ, ncfg_sct, ncfg_pert, ncfg_comp
      Integer :: lcfg_phys, lcfg_targ, lcfg_sct, lcfg_pert, lcfg_comp
      Integer :: nwf_targ, nwf_sct, nwf_pert

      Real(8) :: c_targ = 0.75           
      Real(8) :: c_conf = 0.20           
      Real(8) :: c_comp = 1.10        

      Integer :: max_ll =-1
      Integer :: min_ll =-1
      Integer :: max_ka = 0
      Integer :: min_ka = 0
      Integer :: kort   = 0
      Integer :: debug  = 0
      Integer :: iread_targ  = 0
      Integer :: igen_conf = 0

! ... flags:

      Integer :: ic_comp = 0
      Integer :: ie_comp = 0
      Integer :: ii_comp = 0
      Real(8) :: wc_comp = 0.d0
      Integer :: insert  = 0

! ... coupling scheemes:

      Integer, parameter :: mshells = 25
      Integer, parameter :: mcup=2*mshells   
      Integer, parameter :: mmom=6*mshells   
      Integer :: ncup    ! number of elementary couplings 
      Integer :: nmom    ! number of involed momemntums
      Integer :: J1_coupling(3,mcup)
      Integer :: J2_coupling(3,mcup)
      Integer :: moments(mmom)
      Real(8) :: S_cfp 
      Real(8) :: S_recup 
      Integer :: JP_trap 

! ... additional arrays:

      Integer, allocatable :: ic_targ(:)  
      Integer, allocatable :: jc_targ(:)  
      Integer, allocatable :: ic_pert(:)  
      Real(8), allocatable :: WC_pert(:)

! ... information about channels to delete:

      Character(ma) :: AF_del = 'target_del'
      Integer :: ndel=0
      Integer, allocatable :: dlsp(:),dkch(:),dtar(:)

! ... pertuber index:

      Integer, allocatable :: ip_pert(:)

      End Module dbsr_conf

