!======================================================================
      Subroutine Read_data
!======================================================================
!     read data for given partial wave klsp
!----------------------------------------------------------------------
      Use dbsr_mat           
      Use DBS_nuclear
      Use radial_overlaps 

      Implicit none
      Integer :: i,j,l,k,m,n, ich,jch
      Real(8) :: S,t1,t2
      Integer, external :: Ifind_bsorb, Ifind_jjorb, Icheck_file, &
                           IBORT
      Real(8), external :: QUADR_pq

      write(ALSP,'(i3.3)') klsp   !  sign of partial wave

!----------------------------------------------------------------------
!                                                                files:
! ... log-file:

      i=LEN_TRIM(AF_pri); AF_pri(i-2:i)=ALSP
      Open(pri,file=AF_pri)

! ... c-file:

      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Open(nuc,file=AF_cfg,status='OLD')

! ... bnk-file:
 
      i=LEN_TRIM(AF_bnk); AF_bnk(i-2:i)=ALSP
      Open(nub,file=AF_bnk,status='OLD',form='UNFORMATTED')

! ... out file:

      i=LEN_TRIM(AF_mat); AF_mat(i-2:i)=ALSP
      Open(nui,file=AF_mat,form='UNFORMATTED')

! ... debug file:

      if(debug.gt.0) then
       i=LEN_TRIM(AF_deb); AF_deb(i-2:i)=ALSP
       Open(nud,file=AF_deb)
      end if

!----------------------------------------------------------------------
! ... HEADER ...

      write(pri,'(a/a/)') 'D B S R _ M A T', &
                          '***************'
      write(pri,'(a,i3,5x,a,i5)') &
      'calculations for partial wave:  klsp =',klsp,'nprocs = ', nprocs
!----------------------------------------------------------------------
! ... read obitals and configurations along with expansion coefficients:

      Call CPU_time(t1)

      Call Read_core_jj(nuc); nclosed = ncore; nwf=ncore
      Call Check_conf_jj(nuc,nub)
      Call Read_expn_jj(nuc)

! ... read physical orbitals:

      Call Check_file(AF_orb)
      Open(nuo,file=AF_orb)

      Call Read_sub_orb_jj(nuo,ntarg)

      write(pri,'(/a)')   'c-file data: '
      write(pri,'(/a,i8,a)')  'ncfg   = ',ncfg, &
       '  - number of configurations'
      write(pri,'( a,i8,a)')  'nwf    = ',nwf, &
       '  - number of one-electron orbitals'

! ... create the same space in DBS_orbitals:

      nbf=0
      Do i = 1,nwf; j=Ifind_bsorb(NEF(i),KEF(i),IEF(i),2); End do
      mbs = 0

      Call CPU_time(t2)

      write(pri,'(/a,T40,f10.2,a)') 'Read configurations:',(t2-t1)/60,' min'

!----------------------------------------------------------------------
! ... read channel information and find pointer (orbital --> channel):

      Call Read_channel_jj(nut,klsp)

      if(ncp.gt.0) ippert = ippert + ipconf(nch)

      ipbs = 0
      Do ich = 1,nch
       Call EL_NLJK(ELC(ich),n,k,l,j,i)
       m = Ifind_bsorb(n,k,i,1); ipbs(m) = ich 
       ipch(ich)=m
      End do

      write(pri,'(/a/)')   'Partial wave total moment: '
      write(pri,'(a,i8,a)')  'jpar   = ',jpar,'  - 2J value'
      if(ipar.eq. 1) write(pri,'(a,i8,a)')  'ipar   = ',ipar, &
       '  - even parity'
      if(ipar.eq.-1) write(pri,'(a,i8,a)')  'ipar   = ',ipar, &
       '  - odd parity'

      if(ncp.gt.0)  Call Read_sub_pert_jj(nuo,klsp)

!----------------------------------------------------------------------
! ... read B-spline expantions for bound orbitals:
   
! ... target radial functions:

      Open(nuw, file=AF_bsw, STATUS='OLD', form='UNFORMATTED')
      Call Read_pqbs(nuw)
      Close(nuw)

! ... perturber radial functions:

      if(nwp.gt.0) then
       AF = trim(BFP)//'.bsw'
       Open(nuw, file=AF, STATUS='OLD', form='UNFORMATTED')
       Call Read_pqbs(nuw)
       Close(nuw)
      end if

! ... check the correspondence between c- and bsw-files: 

      j = 0
      Do i = 1,nbf
       if(ipbs(i).ne.0) Cycle
       if(mbs(i).eq.0) then
        write(pri,'(a,a)') ' Absent expansion for w.f. ',ebs(i)
        j = j + 1
       end if
      End do
      if(j.gt.0) Stop 'no correspondence between c- and w- files'

! ... the < B | p > values (convolution with B matrix)  ??? already done in read_pqbs? ???

      Do i=1,nbf
       if(ipbs(i).ne.0) Cycle
       bpq(1:ns,1,i) = MATMUL(fpbs, pq(:,1,i))
       bpq(1:ns,2,i) = MATMUL(fqbs, pq(:,2,i))
      End do

! ... the < p | p > values 

      Call CPU_time(t1)

      Call Alloc_radial_overlaps(0)
      Do i=1,nbf; Do j=1,i; if(kbs(i).ne.kbs(j)) Cycle
       if(ipbs(i).ne.0) Cycle
       if(ipbs(j).ne.0) Cycle
       S = SUM(pq(:,1,i)*bpq(:,1,j)) + SUM(pq(:,2,i)*bpq(:,2,j))  
       if(abs(S).lt.eps_ovl) Cycle
       Call Iadd_obs(i,j,S)                 
      End do; End do

      Call CPU_time(t2)

      write(pri,'(/a,i8,T40,f10.2,a)') 'Radial overlaps:',nobs,(t2-t1)/60,' min'

! ... the < . | j > values 

      Call Alloc_nv_ch(1)
      Call Read_bort_jj(nuc)

      End Subroutine Read_data
    
