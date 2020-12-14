!======================================================================
      Subroutine Def_orbitals
!======================================================================
!     This routine prepares the orbital-related arrays
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals
      Use orb_jj

      Implicit none
      Integer :: i,j,ic,ip,ib,ic1,ic2
      Real(8) :: S
  
      Call alloc_df_orbitals(nwf)
  
      Do i=1,nbf
       nbs(i)=nef(i);  kbs(i)=kef(i);  lbs(i)=lef(i)
       jbs(i)=jef(i);  ebs(i)=ELF(i);  ibs(i)=0
       e(i)=0
      End do
  
      Call Def_nit
  
      Call alloc_df_radial(ns)
  
! ... find physical orbitals:

      Call Def_physical

! ... define qsum:            ?  just initial estimation                     

      WC = 1.d0
      Do i=1,nlevels; ib = block(i)
       ic1=JTc1(ib); ic2=JTc2(ib)
       WC(ic1:ic2) = WC(ic1:ic2) + weight(i)        !  * ???
      End do 
      S = sqrt(SUM(WC(1:ncfg)*WC(1:ncfg)))
      WC(1:ncfg)=WC(1:ncfg)/S

      Do ic = 1,ncfg; S = WC(ic)*WC(ic) 
       Call Get_cfg_jj(ic)     
       ip = ip_state(ic)
       Do i=1,no; ip=ip+1; j=IP_orb(ip)
        qsum(j) = qsum(j) + S*iq(i)  
       End do
      End do

      Do i=1,ncore
       qsum(i) = jbs(i)+1;  clsd(i) = .TRUE.
      End do
   
! ... define closed shells:   ???

      Do i=ncore+1,nbf
       clsd(i) = .FALSE.; S = jbs(i)+1
       if(abs(qsum(i)-S).lt.0.1) clsd(i) = .TRUE.
      End do

      Call Boundary_conditions 

! ... define if core is fixed:

      icore = 0;  if(ncore.gt.0) icore = sum(iord(1:ncore))    ! ???

! ... define fixed integrals:

      Allocate(if_int(nint)); if_int = 0
      Do i = 1,nint
       if(ivaried(i1_int(i)).ne.0) Cycle
       if(ivaried(i2_int(i)).ne.0) Cycle
       if(ivaried(i3_int(i)).ne.0) Cycle
       if(ivaried(i4_int(i)).ne.0) Cycle
       if_int(i) = 1
      End do   

      Allocate(if_Lint(Lint)); if_Lint = 0
      Do i = 1,Lint
       if(ivaried(i1_Lint(i)).ne.0) Cycle
       if(ivaried(i2_Lint(i)).ne.0) Cycle
       if_Lint(i) = 1
      End do   

      End Subroutine Def_orbitals


!======================================================================
      Subroutine Def_nit
!======================================================================
!     This routine determine the orbitals to be optimized
!     nit      -  number of optimized oprbitals
!     iord(:)  -  =1 means optimized orbital
!----------------------------------------------------------------------
      Use dbsr_mchf,  string => avaried,  nit => varied
      Use df_orbitals

      Implicit none
      Integer :: i,j,n,l,k,ip,ii,ios
      Integer, external :: Ifind_orb_df, Ipointer
      Character(5) :: EL

      Call Clean_a(string)
      iord = 0
      if(string(1:3)=='ALL' .or. string(1:3)=='all' .or. &        ! all
       LEN_TRIM(string) == 0 ) then                              
       nit = nbf 
       Do i=1,nbf; iord(i)=i; End do
      elseif(string(1:4)=='NONE' .or. string(1:4)=='none') then   ! none
       nit = 0 
      elseif (INDEX(string,'n=') /= 0) then                       ! n=
       i = INDEX(string,'=') 
       read(string(i+1:),*) n 
       nit = 0
       Do i=1,nbf
        if(nbs(i).ne.n) Cycle
        nit = nit + 1
        iord(nit) = i
       End do                                                     
      elseif (INDEX(string,'n>') /= 0) then                       ! n>
       i = INDEX(string,'>') 
       read(string(i+1:),*) n 
       nit = 0
       Do i=1,nbf
        if(nbs(i).le.n) Cycle
        nit = nit + 1
        iord(nit) = i
       End do                                                   
      elseif (INDEX(string,'=') /= 0) then                        ! =last
       i = INDEX(string,'=') 
       read(string(i+1:),*) nit 
       k=0; Do i=nbf-nit+1,nbf; k=k+1; iord(k)=i; End do
      else
       Do nit=1,nbf                                               ! list
        read(string,*,iostat=ios) (EL,i=1,NIT)
        IF (ios /= 0) Exit
       End do
       nit = nit - 1
       ip = 0
       Do ii=1,nit
        read(string,*) (EL,j=1,ii)
        Call EL_NLJK(EL,n,k,l,j,i)
        j = Ifind_orb_df(n,k,i)
        if(j.gt.0) then; ip=ip+1; iord(ip)=j; end if
       End do
      end if

      Allocate(ivaried(nbf));  ivaried = 0
      Do i=1,nbf
       ivaried(i) = Ipointer(nbf,iord,i)
      End do

      if(debug.gt.0) then
       write(log,'(a,100i3)') 'iord:  ',iord 
       write(log,'(a,100i3)') 'varied:',ivaried
      end if

      End Subroutine Def_nit


!======================================================================
      Subroutine Def_physical
!======================================================================
! ... read or try to find the physical orbitals; approximately,
! ... they are orbitals in the first (leading?) configurations
! ... for each block
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Character(5) :: EL 
      Integer :: i,j,start,ip,ib,ic, n,l,k,iset
      Integer, external :: Ifind_orb_df

      Allocate(iphys(1:nbf)); iphys=0              
      if(ncore.gt.0) iphys(1:ncore)=1

! ... if given - just read:

      if(len_trim(physical).gt.0) then
       start = 1; ip = 0
       Do  
        i = index(physical(start:),',')
        if (i /= 0 .or. LEN_TRIM(physical(start:)) /= 0) then
         read(physical(start:),*) EL
         Call EL_NLJK(EL,n,k,l,j,iset)
         j = Ifind_orb_df(n,k,iset)
         if(j.gt.0) iphys(j) = 1
        end if
        start = start + i 
        if(i == 0 .or.  LEN_TRIM(physical(start:)) == 0) Exit
       End do
       Call Clean_a(physical)
       Return
      end if 

! ... if not given - try to anticipate:
! ... we suupose that the main configurations are in the right place ...

      Do i = 1,nlevels; ib = block(i)
       ic = JTc1(ib) -1 + level(i) 
       Call Get_cfg_jj(ic)     
       ip = ip_state(ic)
       Do j=1,no; ip=ip+1; iphys(IP_orb(ip))=1; End do
      End do

      start=1                                          ! ???
      Do i = ncore+1,nbf; if(iphys(i).eq.0) Cycle
       write(physical(start:start+5),'(a5,a1)') ebs(i),','
      start = start+6
      End do

      i = Len_trim(physical); physical(i:i)=' '
      Call Clean_a(physical)

      End Subroutine Def_physical
