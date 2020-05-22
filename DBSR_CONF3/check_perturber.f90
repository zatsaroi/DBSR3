!======================================================================      
      Subroutine Check_perturber
!======================================================================
!     read the perturber configurations, substitute physical orbitals
!     and check the "double"  configurations if any
!     We are supposed that pertuber expansions are ordered
!----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none
      Integer :: i,j,ii,jj,is,ic,ip,jp  
      Real(8) :: W
      Integer, external :: Ifind_jjorb, Iadd_cfg_jj, Ifind_position 

! ... read substitution orbitals for pertuber:

      Call Read_sub_pert_jj(nuo,ilsp)

! ... find substitution pointers for perturber orbitals:    

      if(nwf_pert.gt.nwf_targ) then
      ipef(nwf_targ+1:nwf_pert)=0
      Do i=1,npert_sub; ii=np_phy(i); jj=np_sub(i)
       if(ipef(ii).eq.0) then
        ipef(ii)=jj
       else
        if(ipef(ii).ne.jj) &
         Stop 'troubles with sustitution orbitals for pertuber'
       end if
      End do
      end if

! ... substitute the orbitals in pertuber and check the "double" configurations:

      if(mcfg.gt.ncfg) WC(ncfg+1:mcfg) = 0.d0

      AF = trim(BFP)//'.c'
      Call Check_file(AF)
      Open(nuc,file=AF,status='OLD')
      Call R_pert(nuc)

      Do ip = 1,npert; jp=1; if(ip.gt.1) jp=ippert(ip-1)+1

      rewind(nuc)
      is = 0
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:1).eq.'*') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ
      is = is + 1
      if(is.lt.jp) go to 1
      Call Decode_cj
      i=INDEX(CONFIG,')',BACK=.TRUE.)+1
      if(LEN_TRIM(CONFIG(i:)).ne.0) Read(CONFIG(i:),*) W

      Do i=1,no
       ii=Ifind_jjorb(nn(i),kn(i),in(i),1)
       j=ipef(ii); if(j.eq.0) Exit
       nn(i)=NEF(j); kn(i)=KEF(j); in(i)=IEF(j)
      End do

      ic = Iadd_cfg_jj('find')
      WC(ic) = WC(ic) + 1.d0

    2 Continue

      End do   !  over perturbers, ip

      i = Ifind_position(nuc,'Spectroscopic configurations:')
      if(i.eq.0) then
      Close(nuc); Open(nuc,file=AF,position='APPEND')
       write(nuc,'(/a/)') 'Spectroscopic configurations:'
      else
       read(nuc,*); read(nuc,*)
      end if

      Do ic=ncfg_sct+1,ncfg
       Call Print_conf_jj(nuc,ic,WC(ic))
      End do
      write(nuc,'(a)') '*'

      if(ncfg-ncfg_sct.ne.npert) then
       write(pri,'(/a,i5)') &
        'Check the double pertubers for given partial wave'
        Stop 'Check the double pertubers'
      end if       

      write(pri,'(/a,T40,i8)') 'Number of perturber configurations:', ncp
      write(pri,'( a,T40,i8)') 'Number of physical perturbers:', &
                                ncfg-ncfg_sct

      End Subroutine Check_perturber
