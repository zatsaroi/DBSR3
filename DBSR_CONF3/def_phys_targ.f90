!======================================================================
      Subroutine Def_phys_targ
!======================================================================     
!     define the main target configuration:  only one !!!;
!     we do it in two stage: define and record them in targ_nnn.c,
!     and then read all of them to avoid situation of missing targets
!     and put target cofigurations in needed order
!----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none
      Integer :: i,j,it,ii,jj,ic,ic1,ic2     
      Real(8) :: W, WW
      Integer, external :: Iadd_cfg_jj, Ifind_jjorb,  &
                           Ifind_position, Ipointer

! ... read substitution orbitals:

      Call Check_file(AF_orb)
      Open(nuo,file=AF_orb)
      Call Read_sub_orb_jj(nuo,ntarg)

! ... find substitution pointers:    

      ipef=0
      Do i=1,nphys_orb; ii=ip_phy(i); jj=ip_sub(i)
       if(ipef(ii).eq.0) then
        ipef(ii)=jj
       else
        if(ipef(ii).ne.jj) write(*,*) 'troubles with sustitution orbitals'
       end if
      End do

      write(pri,'(/a,T33,i8)') 'number of substitution orbitals:',nphys_sub
      write(pri,*)
      write(pri,'(15a5)') (elf(jp_sub(i)),i=1,nphys_sub)

      if(allocated(jc_targ)) Deallocate(jc_targ) 
      Allocate(jc_targ(0:ntarg));  jc_targ = 0      
      if(allocated(ic_targ)) Deallocate(ic_targ) 
      Allocate(ic_targ(0:ntarg));  ic_targ = 0      

      if(iread_targ.eq.1) go to 5 

!----------------------------------------------------------------------
! ... find and record main target configurations

      Call alloc_cfg(0)

      Do it=1,ntarg

      AF = trim(BFT(it))//'.c'
      Open(nuc,file=AF,status='OLD')
      rewind(nuc)
      ii = 0; WW = 0.d0
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'***') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ

      W=0.d0
      i=INDEX(CONFIG,')',BACK=.TRUE.)+1
      if(LEN_TRIM(CONFIG(i:)).ne.0) Read(CONFIG(i:),*) W

      if(abs(W).lt.c_conf) go to 2

      Call Decode_cj

      Do i = 1,no; np(i)=Ifind_jjorb(nn(i),kn(i),in(i),1); End do

      ! ... check if all orbitals are physical:

      Do i=1,no
       j=ipef(np(i)); if(j.eq.0) go to 2
       nn(i)=nef(j); ln(i)=lef(j); kn(i)=kef(j); in(i)=ief(j); np(i)=j
      End do

      ic=Iadd_cfg_jj('add'); jc_targ(it) = ic;   WC(ic) = abs(W) 

      WW = WW + W*W;  if(WW.lt.c_targ) go to 1

    2 Continue

      i = Ifind_position(nuc,'Spectroscopic configuration')

      if(i.eq.0) then
       Close(nuc); Open(nuc,file=AF,position='APPEND')
       write(nuc,'(/a,f15.3/)') 'Spectroscopic configuration:',WW
      else
       write(nuc,'(a,f10.3/)') 'Spectroscopic configuration:',WW      
      end if

      ic1 = jc_targ(it-1)+1; ic2 = jc_targ(it)

      if(ic1.gt.ic2) then
       write(pri,*) 'Can not find physical configuration for target ',it
       Stop 'Can not find physical configuration'
      end if

      WW = sqrt(WW)
      Do ic = ic1,ic2; W=WC(ic)/WW; Call Print_conf_jj(nuc,ic,W); End do

      write(nuc,'(a)') '*'

      End do  !  over target states,  it

!-----------------------------------------------------------------------
! ... read dominant configurations for target states

    5 Continue
      ic_targ = 0; jc_targ = 0      
      Call alloc_cfg(0)

      Do it=1,ntarg

      AF = trim(BFT(it))//'.c'
      Open(nuc,file=AF,status='OLD')
      i = Ifind_position(nuc,'Spectroscopic configuration')
      if(i.eq.0) then
       write(pri,*) 'Cannot find spectroscopic configuration for target', it
      end if

   10 read(nuc,'(a)',end=20) CONFIG
      if(CONFIG(1:3).eq.'*') go to 20
      if(CONFIG(6:6).ne.'(') go to 10
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ
      W=0.d0
      i=INDEX(CONFIG,')',BACK=.TRUE.)+1
      if(LEN_TRIM(CONFIG(i:)).ne.0) Read(CONFIG(i:),*) W
      Call Decode_cj
      i=Iadd_cfg_jj('add')
      WC(i) = W
   20 Continue

      jc_targ(it) = ncfg

      End do  !  over target states,  it

      ncfg_phys = ncfg;  lcfg_phys = lcfg  
      write(pri,'(/a,T33,i8)') 'number of phys. target config.s:',ncfg_phys

      End Subroutine Def_phys_targ

