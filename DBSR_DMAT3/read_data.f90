!======================================================================
      Subroutine Read_data 
!======================================================================
      Use dbsr_dmat;    Use conf_jj;       
                        Use channels_jj;   Use DBS_gauss 
                        Use target_jj;     Use DBS_orbitals_pq
                        Use orb_jj        
      Implicit none

      Integer :: i,ii,met, n,l,k,j,iset
      Real(8) :: S
      Integer, external :: Ifind_bsorb, l_kappa, ITRA, IBORT
      Real(8), external :: QUADR_pq, DQUADRM_pq, DQUADRL_pq, DQUADRV_pq, OBS

! ... read arguments:

      Call Read_arg

! ... sets up grid points and initializes the B-spline arrays: 
    
      Call read_knot_dat
      Call alloc_DBS_gauss

! ... check angular coefficient:

      Call Check_mult_bnk

! ... read channel information from target:

      if(ilsp1.gt.0.or.ilsp2.gt.0) then
       Call Check_file(AF_tar)
       Open(nut,file=AF_tar)
       Call Read_target_jj (nut)
       Call Read_channels_jj(nut)
       Close(nut)
      end if

      nch1=0; npert1=ncfg1
      if(ilsp1.gt.0) then
       nch1=nch(ilsp1)
       ncp1=ncp(ilsp1)
       npert1=npert(ilsp1)
       jot1=jpar(ilsp1)
       parity1=ipar(ilsp1)
      end if
      Allocate(C1(ncfg1))
      C1=WC(1:ncfg1);  if(ctype1.eq.'j') C1=1.d0
     
      nch2=0; npert2=ncfg2
      if(ilsp2.gt.0) then
       nch2=nch(ilsp2)
       ncp2=ncp(ilsp2)
       npert2=npert(ilsp2)
       jot2=jpar(ilsp2)
       parity2=ipar(ilsp2)
      end if
      Allocate(C2(ncfg2))
      C2=WC(ncfg1+1:ncfg1+ncfg2); if(ctype2.eq.'j') C2=1.d0

      kdm1= ms*nch1+npert1
      kdm2= ms*nch2+npert2
        
      write(pri,'(a,a,a,i6,a,i4,a,i8,a,i3/)' )  'initial state - ',trim(name1), &
       '   ncfg1 =',ncfg1,'   nwf1 =',nwf1,'   kdm1 =',kdm1,'   jot1 =',jot1
      write(pri,'(a,a,a,i6,a,i4,a,i8,a,i3/)' )  'final state   - ',trim(name2), &
       '   ncfg2 =',ncfg2,'   nwf2 =',nwf2,'   kdm2 =',kdm2,'   jot2 =',jot2


! ... allocate orbitals arrays:
   
      Call alloc_DBS_orbitals_pq(nwf,ns)
      nbf = nwf
      Do i = 1,nbf
       ebs(i)=ELF(i); nbs(i)=NEF(i); kbs(i)=KEF(i); ibs(i)=IEF(i);
       lbs(i)=LEF(i); jbs(i)=JEF(i); ipbs(i) = 0; mbs(i) = 0
      End do

! ... read radial bound orbitals for the initial state:

      if(ilsp1.eq.0) then

       iname1=INDEX(name1,'.'); AF=name1(1:iname1)//'bsw'
       Call Check_file(AF); Open(nuw,file=AF,form='UNFORMATTED')
       Call Read_dbsw(nuw,1,kset1);  Close(nuw)

      else

       Open(nuw,file=AF_bsw,form='UNFORMATTED',status='OLD')
       Call Read_dbsw(nuw,1,kset1);  Close(nuw)

       if(nwp(ilsp1).gt.0) then

        i = LEN_TRIM(BFP(ilsp1)); AF = BFP(ilsp1); AF = AF(1:i)//'.bsw'

        Open(nuw,file=AF,form='UNFORMATTED',status='OLD')
        Call Read_dbsw(nuw,1,kset1);  Close(nuw)
       end if

       Do i = 1,nch1
        Call EL_nljk(ELC(ilsp1,i),n,k,l,j,iset);  iset=iset+kset1
        ii = Ifind_bsorb(n,k,iset,1); ipbs(ii)=i
       End do

      end if

! ... radial functions for final state:

      if(ilsp2.eq.0) then

       iname2=INDEX(name2,'.'); AF=name2(1:iname2)//'bsw'
       Open(nuw,file=AF,form='UNFORMATTED',status='OLD')
       Call Read_dbsw(nuw,1,kset2);  Close(nuw)

      else

       Open(nuw,file=AF_bsw,form='UNFORMATTED',status='OLD')
       Call Read_dbsw(nuw,1,kset2);  Close(nuw)

       if(nwp(ilsp2).gt.0) then

        i = LEN_TRIM(BFP(ilsp2)); AF = BFP(ilsp2); AF = AF(1:i)//'.bsw'

        Open(nuw,file=AF,form='UNFORMATTED',status='OLD')
        Call Read_dbsw(nuw,1,kset2);  Close(nuw)
       end if

       Do i = 1,nch2
        Call EL_nljk(ELC(ilsp2,i),n,k,l,j,iset);  iset=iset+kset2
        ii = Ifind_bsorb(n,k,iset,1); ipbs(ii)=i
       End do

      end if

! ... check the correspondence between c- and w-files: 

      met = 0
      Do i = 1,nwf  
        if(ipbs(i).ne.0) Cycle; if(mbs(i).ne.0) Cycle
        write(pri,'(a,a)') ' Absent expansion for w.f. ',ELF(i)
        met = met + 1
      End do
      if(met.gt.0) Stop 'no correspondence between c- and w- files'
      
! ... the < . | p > values (convolution with B matrix)

      Do i=1,nbf
       if(ipbs(i).ne.0) Cycle
       bpq(1:ns,1,i) = MATMUL(fpbs, pq(:,1,i))
       bpq(1:ns,2,i) = MATMUL(fqbs, pq(:,2,i))
      End do

! ... the < p | p > values 

      Call Alloc_radial_overlaps(0)
      Do i=1,nwf1; Do j=nwf1+1,nwf; if(kbs(i).ne.kbs(j)) Cycle
       if(ipbs(i).ne.0) Cycle
       if(ipbs(j).ne.0) Cycle
       S = QUADR_pq(i,j,0)
       if(abs(S).lt.eps_ovl) Cycle
       Call Iadd_obs(i,j,S)
      End do; End do

      Call Alloc_nv_ch(1)
      Call Get_v_ch(1)

      Do i=1,nwf1; Do j=nwf1+1,nwf; if(kbs(i).ne.kbs(j)) Cycle
       if(ipbs(i).eq.0.and.ipbs(j).eq.0) Cycle
       if(ipbs(i).ne.0.and.ipbs(j).ne.0) Cycle
      End do; End do

! ... dipole integrals:

      Call Gen_dbs(nwf1,nwf2,0.d0,kpol,ktype)  

! ... prepare B-spline dipole matrixes:

      Allocate(dipLp(ns,ns),dipLq(ns,ns),dbs(ns,2,nbf))

      Call ZINTYk (kpol,ksp,ksp,pbsp,pbsp,ns,dipLp)
      Call ZINTYk (kpol,ksq,ksq,qbsp,qbsp,ns,dipLq)

      Do i=1,nbf; if(ipbs(i).ne.0) Cycle
       dbs(1:ns,1,i) = MATMUL(dipLp,pq(1:ns,1,i))
       dbs(1:ns,2,i) = MATMUL(dipLq,pq(1:ns,2,i))
      End do
                                                                          
      if(ktype.eq.'E') then
       Allocate(dipVp(ns,ns),dipVq(ns,ns),vbs(ns,2,nbf))
       Call ZINTYk (kpol-1,ksp,ksq,pbsp,qbsp,ns,dipVp) 
       Call ZINTYk (kpol-1,ksq,ksp,qbsp,pbsp,ns,dipVq)
       Do i=1,nbf; if(ipbs(i).ne.0) Cycle
        vbs(1:ns,1,i) = MATMUL(dipVq,pq(1:ns,1,i))         !  ??? Vq     Vp
        vbs(1:ns,2,i) = MATMUL(dipVp,pq(1:ns,2,i))
       End do
      end if


      End Subroutine Read_data

