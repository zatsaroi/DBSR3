!=====================================================================
!     UTILITY   M U L T _  J J
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny,   email: oleg_zoi@yahoo.com
!======================================================================
!     for selected states, print (to file mult_jj.tab)  the integral 
!     coefficients according to "mult_bnk" file
!
!     Call as:
!
!     mult_jj  c1=.. c2=..  bnk=..  tab=.. jort=..  ic=..  jc=..
!
!     all parameters - optional, with default values:
!
!     c1  = cfg.inp
!     c2  = cfg.out
!     bnk = mult_bnk_E1
!     tab = mult_jj.tab
!
!     jort = 1 -> partial orthogonality
!
!     ic = 0   -> all initial states
!     jc = 0   -> all final states
!----------------------------------------------------------------------
      Use orb_jj;    Use det_list;    Use symc_list;    Use zoef_list
      Use conf_jj;   Use def_list;    Use symt_list

      Implicit real(8) (A-H,O-Z)

! ... the names and units for input/output files:

      Integer, parameter :: ma=40
      Character(ma) :: AF

! ... c-file:
      Integer :: nu1=10;  Character(ma) :: AF_c1 = 'cfg.inp'
      Integer :: nu2=20;  Character(ma) :: AF_c2 = 'cfg.out'
                         
! ... data bank:
      Integer :: nub=30;  Character(ma) :: AF_bnk = 'mult_bnk_E1'

! ... output:
      Integer :: nur=40;  Character(ma) :: AF_tab = 'mult_jj.tab'
                         
! ... scratch file for intermediate results:
      Integer :: nui=11  

! ... other variables:

      Character(1) :: ktype       ! type of trabsition
      Integer :: kpol             ! multipole index 

      Real(8) :: eps_c = 1.d-7    ! tolerence for coefficients:    

      Integer :: parity1, parity2

!----------------------------------------------------------------------
!                                                                files:
! ... c-files:
      Call Read_aarg('c1',AF_c1)
      Call Check_file(AF_c1)
      Open(nu1,file=AF_c1,status='OLD')
      Call Read_aarg('c2',AF_c2)
      Call Check_file(AF_c2)
      Open(nu2,file=AF_c2,status='OLD')

! ... MULT_BNK file:
      Call Read_aarg('bnk',AF_bnk)
      Call Check_file(AF_bnk)
      Open(nub,file=AF_bnk,form='UNFORMATTED',status='OLD')

! ... output coef.tab file:
      Call Read_aarg('tab',AF_tab)
      Open(nur,file=AF_tab)

! ... scratch file:
      Open(nui,form='UNFORMATTED',status='SCRATCH')

!----------------------------------------------------------------------
!                          check the data bank and read configurations:
! ... type of transition:
      read(nub) ktype,kpol
      write(nur,'(a,a1)') 'ktype = ',ktype
      write(nur,'(a,i1)') 'kpol  = ',kpol

! ... conf. symmetries in bank:
      Call Read_symc(nub)
      Call Read_symt(nub)

! ... define conf. symmetries in input c-files:

      Call Read_conf_jj(nu1,0,'add','check')
      ne1=ne; parity1=parity; ncfg1=ncfg
      Call Read_conf_jj(nu2,0,'add','check')
      ne2=ne; parity2=parity; ncfg2=ncfg-ncfg1

      if(ne1.ne.ne2) Stop 'mult_jj: ne1 <> ne2'
      if(ktype.eq.'E'.and.mod(kpol,2).eq.1.and. &
         parity1.eq.parity2) Stop 'mult_jj: parity ?'
      if(ktype.eq.'M'.and.mod(kpol,2).eq.0.and. &
         parity1.eq.parity2) Stop 'mult_jj: parity ?'

! ... define IT_stat pointer to given configurations:
      if(Allocated(IT_stat )) Deallocate(IT_stat)
                              Allocate  (IT_stat(nsymt))
      IT_stat=0
      Do i=1,nsymt
       if(Ipointer(ncfg,IS_term,i).ne.0) IT_stat(i)=1
      End do

! ... done pointer:
      if(allocated(IT_done)) Deallocate(IT_done)
                             Allocate(IT_done(nsymt*(nsymt+1)/2))
      IT_done=0    
      Call Read_done(nub)

! ... determinants:
      Call Read_det(nub)
      Call Read_def(nub)

! ... read and save all relevant coefficients:
      ncoef=0
    1 read(nub,end=2) C,it,jt,int,idf
      if(IT_stat(it).eq.0) go to 1
      if(IT_stat(jt).eq.0) go to 1
      write(nui) C,it,jt,int,idf
      ncoef = ncoef + 1
      go to 1
    2 Continue

      write(nur,'(/a/)') 'MULT_BNK contains: '
      write(nur,'(a,2i8)')  'ndet     = ',ndet,ldet
      write(nur,'(a,2i8)')  'ndef     = ',ndef,ldef
      write(nur,'(/a,i8/)') 'ncoef    = ',ncoef

!----------------------------------------------------------------------
!                                             orthogonality conditions:

      jort = 1;  Call Read_iarg('jort',jort)

      if(JORT.eq.-1)  write(nur,'(/a/)') &
       'Orthogonality mode:  full orthogonality (jort = -1)'
      if(JORT.eq. 0)  write(nur,'(/a/)') &
       'Orthogonality mode:  full non-orthogonality (jort = 0)'
      if(JORT.eq. 1)  write(nur,'(/a/)') &
       'Orthogonality mode: partial orthogonality (jort = 1)'

      Call prepare_iort_jj(0) 

!----------------------------------------------------------------------
!                                      input required matrix elements:

      ic = 0; Call Read_iarg('ic',ic)
      jc = 0; Call Read_iarg('jc',jc)

      Do ic1 = 1,ncfg1
       if(ic.ne.0.and.ic1.ne.ic) Cycle
      Do ic2 = ncfg1+1,ncfg
       if(jc.ne.0.and.ic2.ne.jc_ncfg1) Cycle
       
       Call SUB1

      End do; End do

Contains 

!======================================================================
      Subroutine SUB1
!======================================================================

! ... find corr. it,jt

      it1 = IS_term(ic1)
      it2 = IS_term(ic2)

! ... find the indication on done calculation:

      i = max(it1,it2); j = min(it1,it2)
      if(IT_done((i-1)*i/2 + j).ne.1) then
       Stop 'run dbsr_mult first!'
      end if

!----------------------------------------------------------------------
! ... read the specific data:

      Call Alloc_zoef(-1); Call Alloc_ndet(0); Call Alloc_ndef(0)
      rewind(nui)
   10 read(nui,end=20) C,it,jt,int,idf
write(*,*) C,it,jt,int,idf
      met = 0
      if(it.eq.it1.and.jt.eq.it2) met=1
      if(it.eq.it2.and.jt.eq.it1) met=2
      if(met.eq.0) go to 10

! ... define integral:

      Call Decode_mult (icase,i1,i2,int)
      if(ktype.eq.'E'.and.kpol.eq.0) icase=0

write(*,*) icase,i1,i2

      if(met.eq.1) then
       ip1 = IP_state(ic1); ip2 = IP_state(ic2)
       j1=IP_orb(i1+ip1);   j2=IP_orb(i2+ip2)
      else
       ip1 = IP_state(ic2); ip2 = IP_state(ic1)
       j1=IP_orb(i1+ip1);   j2=IP_orb(i2+ip2)
       m=j1; j1=j2; j2=m                               !  ???  sign

      end if	 


      int = Iadd_int(icase,kpol,j1,j2,0,0) 

!----------------------------------------------------------------------
! ... find determinant overlaps for specific orbitals

      kz = 0
      if(idf.gt.0) then

      kd=KPF(idf); ip=IPF(idf); NP(1:kd)=NPF(ip+1:ip+kd); md=0
      Do ii=1,kd
       id=NP(ii)/ibf; iext=mod(NP(ii),ibf); nd=KPD(id); ip=IPD(id)
       Do i=1,nd
        k=NPD(i+ip); k1=k/ibd; k2=mod(k,ibd)
        NP1(i)=IP_orb(ip1+k1); NP2(i)=IP_orb(ip2+k2)
       End do

       mm = IDET_SIMP(kz,nd,NP1,NP2)

       if(mm.eq.1) Cycle; if(mm.eq.0) go to 10
       MP(1:nd) = NP1(1:nd)*ibd+NP2(1:nd) 
       jd = Iadd_ndet(nd,MP)
       md = md + 1; MP1(md) = jd; MP2(md) = iext

      End do 
    
       idf = 0
       if(md.gt.0) then 
        MP(1:md) = MP1(1:md)*ibf + MP2(1:md)
        idf = Iadd_ndef(md,MP)
       end if
      end if   !  idf > 0
!----------------------------------------------------------------------

      C = C * (-1)**kz;  Call Iadd_zoef(C,int,idf)

      go to 10
   20 if(nzoef.eq.0) Return
!----------------------------------------------------------------------
! ... print HEADER:

      write(nur,'(/70(''=''))')
      write(nur,'(12x,''< state'',i3,'' || D || state'',i3,''>'')') &
            ic1,ic2-ncfg1
      write(nur,'(70(''=''))') 
      Call Print_conf_jj (nur,ic1,0.d0);  write(nur,'(70(''-''))')
      Call Print_conf_jj (nur,ic2,0.d0);  write(nur,'(70(''-''))')

      Do i=1,nzoef 
       if(abs(Zoef(i)).lt.eps_c) Cycle
       Call Pri_coef1(nur,i)
      End do
      write(nur,'(70(''-''))')

      End Subroutine SUB1

      End ! program MULT_JJ
 


!======================================================================
      Subroutine Pri_coef1 (nu,ii)
!======================================================================
!     Prints one angular integral "ii"
!----------------------------------------------------------------------
      Use orb_jj
      Use zoef_list; Use int_list; Use ndet_list; Use ndef_list

      Implicit none

      Integer, intent(in) :: nu,ii
      Integer :: i,j, i1,i2,i3,i4, k,k1,k2, m1,m2, &
                 ia,ip,np,idf,int,kd,nd,ns,ni
      Character(1) :: AI(3)
      Data AI /'d','a','b'/
      Character(5) :: EL,EL1,EL2,EL3,EL4
      Character(320) :: A
      Character(1) :: a1,b1,a2

      i=IZ_int(ii); int=Jcase(i);  k=Jpol(i); if(int.gt.1) k=k-1

      i1=J1(i); i2=J2(i)

      write(A,'(220a1)') (' ',i=1,220)

      Select case(int)

      Case(0)                            ! Overlap

       ia=0

      Case(1,2,3)                        ! One-electron integrals

       EL1 = ELF(I1); EL2 = ELF(I2)
       Write(A(1:17),'(a1,i1,a3,2(a5,a1))') &
                       AI(int),k,'  (',EL1,',',EL2,')'
       ia=17

      Case default

       write(*,*) ' int=',int
       Stop ' Pri_coef1: unknown case'

      End select

! ... determinant factor:

      idf=IZ_df(ii)
      if(idf.gt.0) then
       kd=KPF(idf);  ip=IPF(idf)
       Do i=1,kd
        ns=mod(NPF(ip+i),ibf); ni=NPF(ip+i)/ibf; 
        nd=KPD(ni); np=IPD(ni)
        ia=ia+1; A(ia:ia)='<'
       Do j=1,nd
        ni=NPD(np+j)/ibd;  EL = ELF(ni)(2:5)
        Write(A(ia+1:ia+6),'(10(a4,'',''))') EL
        ia=ia+5
       End do
       A(ia:ia)='|'
       Do j=1,nd
        ni=mod(NPD(np+j),ibd); EL = ELF(ni)(2:5)
        Write(A(ia+1:ia+6),'(10(a4,'',''))') EL
        ia=ia+5
       End do
       if(ns.gt.1) then
        A(ia:ia+1)='>^'; ia=ia+2; write(A(ia:ia),'(i1)') ns; ia=ia+1
       else
        A(ia:ia+1)='> '; ia=ia+2
       end if
      End do       ! Over kd
      end if

      Call Num(Zoef(ii),k1,k2,999999,1.d-9)

      a1='['; b1=':';a2=']'
      k1 = iabs(k1); m1=sqrt(1.*k1)+0.1;m2=sqrt(1.*k2)+0.1
      if(k1.eq.m1*m1.and.k2.eq.m2*m2) then
      k1=m1;k2=m2;a1='(';a2=')';end if
      Write(nu,'(f12.6,1x,a1,i6,a1,i6,a1,2x,220a1)') &
                Zoef(ii),A1,iabs(k1),B1,k2,A2,(a(i:i),i=1,ia)

      End Subroutine Pri_coef1
