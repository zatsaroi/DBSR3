!=====================================================================
!     UTILITY   C O E F _  J J
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny,   email: oleg_zoi@yahoo.com
!                         
!======================================================================
!     print (to file coef_jj.tab) the integral coefficients for two 
!     selected states in cfg.inp according to int_bnk file
!
!     coef_jj [name] cfg=  bnk=  tab=  jort=  ic=  jc=  mbreit=
!----------------------------------------------------------------------
      Use orb_jj;    Use symc_list;   Use det_list  
      Use conf_jj;   Use symt_list;   Use def_list  
                                       
      Implicit real(8) (A-H,O-Z)

      Character(40) :: name = ' '

! ... c-file:
      Integer :: nuc=1;  Character(40) :: AF_c = 'cfg.inp'
                         
! ... data bank:
      Integer :: nub=2;  Character(40) :: AF_b = 'int_bnk'

! ... output:
      Integer :: nur=7;  Character(40) :: AF_tab = 'coef_jj.tab'
                         
! ... scratch files:
      Integer :: nui=11  ! intermediate results

! ... tolerence for coefficients:
      Real(8) :: eps_c = 1.d-7      

! ... Breit interaction:
      Integer :: mbreit = 0 

!----------------------------------------------------------------------
!                                                                files:
      Call Read_name(name)

      if(name == '?') then
       write(*,*) 
       write(*,*) 'Call coef_jj as:' 
       write(*,*) 
       write(*,*) 'coef_jj name cfg= bnk= tab=  jort=  ic=  jc= mbreit=' 
       write(*,*) 
       write(*,*) 'all paramters are optional:'
       write(*,*) 
       write(*,*) 'cfg  -  input c-file [cfg.inp] or name.c' 
       write(*,*) 'bnk  -  input bank-file [int_bnk] or name.bnk'
       write(*,*) 'tab  -  output coef.s [coef_jj.tab] or name.coef'
       write(*,*) 
       write(*,*) 'if name is given: cfg -> name.c; bnk -> name.bnk; tab -> name.tab '
       write(*,*) 
       write(*,*) 'ic   -  LHS state [0, all LHS states]'
       write(*,*) 'jc   -  RHS state [0, all RHS states]'
       write(*,*) 
       write(*,*) 'jort -  orbital orthogonality mode [-1]'
       write(*,*) '        =-1, full orthogonality '
       write(*,*) '        =+1, partial orthogonality '
       write(*,*) '        = 0, full non=orthogonality '
       write(*,*) 
       write(*,*) 'mbreit[0] - if = 1, Breit integrals are included'
       Stop ' '
      end if

      if(len_trim(name).ne.0) then
       AF_c = trim(name)//'.c'
       AF_b = trim(name)//'.bnk'
       AF_tab = trim(name)//'.coef'
      end if

! ... cfg.inp file:
      Call Read_aarg('cfg',AF_c)
      Call Check_file(AF_c)
      Open(nuc,file=AF_c,status='OLD')

! ... int.bnk file:
      Call Read_aarg('bnk',AF_b)
      Call Check_file(AF_b)
      Open(nub,file=AF_b,form='UNFORMATTED',status='OLD')

! ... output coef.tab file:
      Call Read_aarg('tab',AF_tab)
      Open(nur,file=AF_tab)       

! ... scratch file:
      Open(nui,form='UNFORMATTED',status='SCRATCH')

! ... additional parameters:
      Call Read_iarg('mbreit',mbreit)
      Call Read_rarg('eps_c' ,eps_c )
!----------------------------------------------------------------------
!                          check the data bank and read configurations:

! ... conf. symmetries in bank:
      Call Read_symc(nub);  nsymcb=nsymc
      Call Read_symt(nub);  nsymtb=nsymt

! ... define conf. symmetries in input c-files:
      Call Read_conf_jj(nuc,0,'add','check')

! ... define IT_stat pointer to given configurations:
      if(allocated(IT_stat )) Deallocate(IT_stat)
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

      write(nur,'(a/)') 'INT_BNK contains: '
      write(nur,'(a,i8,a)')  'nsymc = ',nsymcb,' - number of conf. symmetries'
      write(nur,'(a,i8,a)')  'nsymt = ',nsymtb,' - number of ang.  symmetries'
      write(nur,'(a,i8,a)')  'ndet  = ',ndet,  ' - number of determnant overlaps'
      write(nur,'(a,i8,a)')  'ndef  = ',ndef,  ' - number of determnant overlaps'
      write(nur,'(a,i8,a)')  'ncoef = ',ncoef, ' - number of angular coefficients'

      write(nur,'(/a/) ') 'c-fail (cfg.inp) contains: '
      write(nur,'(a,i8,a)')  'nsymc = ',nsymc,' - number of conf. symmetries'
      write(nur,'(a,i8,a)')  'nsymt = ',nsymt,' - number of ang.  symmetries'
      write(nur,'(a,i8,a)')  'ncfg  = ',ncfg, ' - number of atomic states'

      if(nsymc.ne.nsymcb) write(nur,'(/a)') &
       'not all conf. symmetries are covered by INT_BNK'
      if(nsymt.ne.nsymtb) write(nur,'(/a)') &
       'not all ang.  symmetries are covered by INT_BNK'

!----------------------------------------------------------------------
! ... orthogonality conditions for overlap factors:

      jort = -1;  Call Read_iarg('jort',jort)

      write(nur,'(/a/) ') 'Orthogonality mode: '
      if(JORT.eq.-1)  &
      write(nur,'(a,i8,a)')  'jort  = ',jort,' - full orthogonality'
      if(JORT.eq. 0)  &
      write(nur,'(a,i8,a)')  'jort  = ',jort,' - full non-orthogonality'
      if(JORT.eq. 1)  &
      write(nur,'(a,i8,a)')  'jort  = ',jort,' - partial orthogonality'

      Call Prepare_iort_jj

      write(nur,'(/a/) ') 'output of Breit integrals: '
      if(mbreit.eq. 0)  &
      write(nur,'(a,i8,a)')  'mbreit= ',mbreit,' - no Breit integrals'
      if(mbreit.eq. 1)  &
      write(nur,'(a,i8,a)')  'mbreit= ',mbreit,' - Breit integrals included'

      write(nur,'(/a,1PE8.1,a)')  &
      'eps_c = ',eps_c,' - tollerence for angular coefficients'

!----------------------------------------------------------------------
! ... output required matrix elements:

      ic = 0; Call Read_iarg('ic',ic)
      jc = 0; Call Read_iarg('jc',jc)

      Do ic1 = 1,ncfg;     if(ic.ne.0.and.ic1.ne.ic) Cycle
      Do ic2 = ic1,ncfg;   if(jc.ne.0.and.ic2.ne.jc) Cycle

       Call SUB1

      End do; End do

 Contains 

!======================================================================
      Subroutine SUB1
!======================================================================
!     processing one matrix element
!----------------------------------------------------------------------
      Use zoef_list
      Use int_list, only:nint

      Real(8) :: S(8), SMU

! ... find corr. it,jt

      it1 = IS_term(ic1)
      it2 = IS_term(ic2)

! ... find the indication on done calculation:

      i = max(it1,it2); j = min(it1,it2)
      if(IT_done((i-1)*i/2 + j).ne.1) then
       write(nur,'(/70(''-''))')
       write(nur,'(12x,''< state'',i2,'' || (H-E) || state'',i2,''>'')') &
        ic1,ic2
       write(nur,'(70(''-''))')
       write(nur,'(a)') 'This matrix element is not covered by INT_BNK'
       write(nur,'(70(''-''))')
      end if
!----------------------------------------------------------------------
! ... read the specific data:

      ifirst=1; nint=0
    1 Continue
      rewind(nui)
      nzoef=0; Call Alloc_ndet(0); Call Alloc_ndef(0)
   10 read(nui,end=20) C,it,jt,int,idf

      m = 0
      if(it.eq.it1.and.jt.eq.it2) m=1
      if(it.eq.it2.and.jt.eq.it1) m=2
      if(m.eq.0) go to 10

! ... define integral:

      Call Decode_int(icase,kpol,I1,I2,I3,I4,int)

      if(m.eq.1) then
       ip1 = IP_state(ic1); ip2 = IP_state(ic2)
      else
       ip1 = IP_state(ic2); ip2 = IP_state(ic1)
      end if

      j1=IP_orb(i1+ip1); j2=IP_orb(i2+ip1)
      j3=IP_orb(i3+ip2); j4=IP_orb(i4+ip2)

      Call Jsym_int(icase,j1,j2,j3,j4)
      int = Iadd_int(icase,kpol,j1,j2,j3,j4) 

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

       mm = IDET_SIMP(kz,nd,NP1,NP2,mwf,IORT)

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
      C = C * (-1)**kz

      if(icase.le.1) then 
       Call Iadd_zoef(C,int,idf)
       go to 10
      end if

      m = 1
      l1=lef(j1);l2=lef(j2);l3=lef(j3);l4=lef(j4) 
      if(mod(l1+l3+kpol,2).ne.0) m=0
      if(mod(l2+l4+kpol,2).ne.0) m=0
      if(m.eq.1) Call Iadd_zoef(C,int,idf)

      icase = 3
      k1=kef(j1);k2=kef(j2);k3=kef(j3);k4=kef(j4) 
      Do k = kpol-1,kpol+1
        if(k.lt.0) Cycle
        if(SMU(k1,k2,k3,k4,kpol,k,S).eq.0.d0) Cycle
        S = S * C
        int = Iadd_int(icase,k,j1,j2,j3,j4); Call Iadd_zoef(S(1),int,idf)
        int = Iadd_int(icase,k,j2,j1,j4,j3); Call Iadd_zoef(S(2),int,idf)
        int = Iadd_int(icase,k,j3,j4,j1,j2); Call Iadd_zoef(S(3),int,idf)
        int = Iadd_int(icase,k,j4,j3,j2,j1); Call Iadd_zoef(S(4),int,idf)
        int = Iadd_int(icase,k,j1,j4,j3,j2); Call Iadd_zoef(S(5),int,idf)
        int = Iadd_int(icase,k,j4,j1,j2,j3); Call Iadd_zoef(S(6),int,idf)
        int = Iadd_int(icase,k,j3,j2,j1,j4); Call Iadd_zoef(S(7),int,idf)
        int = Iadd_int(icase,k,j2,j3,j4,j1); Call Iadd_zoef(S(8),int,idf)
      End do  
      go to 10

   20 Continue

!----------------------------------------------------------------------
!
      if(ifirst.eq.1) then; ifirst = 0;  Call Sort_int; go to 1; end if

      Call Pri_coef(nur,ic1,ic2,eps_c,mbreit)

      End Subroutine SUB1


      End ! program COEF_JJ
 

!======================================================================
      Subroutine SORT_int
!======================================================================
!    sort two arrays simultaniously (fist N1, then N2) and returns also
!    the number of required permutations
!----------------------------------------------------------------------
      Use int_list

      Implicit none
      Integer :: i1,i2, k

      Do i1=1,nint-1
       Do i2=i1+1,nint
       if(j1(i1).gt.j1(i2).or.&
         (j1(i1).eq.j1(i2).and.j2(i1).gt.j2(i2))) then
        k=j1(i1); j1(i1)=j1(i2); j1(i2)=k
        k=j2(i1); j2(i1)=j2(i2); j2(i2)=k
        k=j3(i1); j3(i1)=j3(i2); j3(i2)=k
        k=j4(i1); j4(i1)=j4(i2); j4(i2)=k
        k=jcase(i1); jcase(i1)=jcase(i2); jcase(i2)=k
        k=jpol (i1); jpol (i1)=jpol (i2); jpol (i2)=k
       end if
       End do
      End do

      End Subroutine SORT_int


!======================================================================
      Subroutine Jsym_int(met,i1,i2,i3,i4)
!======================================================================
!     use integral symmetry to obtaine the 'canonical' form
!----------------------------------------------------------------------
      m=1; ii=min(i1,i2,i3,i4)

      Select case(met)
       Case(2)                                !   R-integrals
        if(ii.eq.i4) m = 4
        if(ii.eq.i3) m = 3
        if(ii.eq.i2) m = 2
       Case(1)                                !   L-integrals
        if(ii.eq.i3) m = 3
      End select

      if(m.eq.1) return

      j1 = i1; j2 = i2; j3 = i3; j4 = i4

      if(m.eq.2) then
        i1 = j2; i2 = j1; i3 = j4; i4 = j3
      elseif(m.eq.3) then
        i1 = j3; i2 = j4; i3 = j1; i4 = j2
      elseif(m.eq.4) then
        i1 = j4; i2 = j3; i3 = j2; i4 = j1
      end if

      End Subroutine Jsym_int


!====================================================================
      Subroutine Pri_coef (nu,ic1,ic2,eps_c,mbreit)
!====================================================================
!     Prints the resulting angular integrals
!--------------------------------------------------------------------
      Use conf_jj;  Use zoef_list;  Use int_list

      Implicit real(8) (A-H,O-Z)
      Integer, allocatable :: ipt(:)

! ... print HEADER:

      write(nu,'(/70(''=''))')
      write(nu,'(12x,''< state'',i2,'' || (H-E) || state'',i2,''>'')') &
        ic1,ic2
      write(nu,'(70(''=''))') 
      Call Print_conf_jj (nu,ic1,0.d0);  write(nu,'(70(''-''))')
      Call Print_conf_jj (nu,ic2,0.d0);  write(nu,'(70(''-''))')

      if(nzoef.eq.0) Return

      if(allocated(ipt)) Deallocate(ipt)
      Allocate(ipt(nzoef))

      Call SORTI(nzoef,IZ_int,ipt)

! ... overlaps:

      write(nu,'(/12x,'' Overlaps''/)')
      Do i=1,nzoef; Call Pri_coef1(nu,i,0,eps_c); End do

! ... one-electron integrals:

      write(nu,'(/12x,'' One-electron integrals:''/)')
      Do i=1,nzoef; Call Pri_coef1(nu,ipt(i),1,eps_c); End do
       
! ... two-electron Coulomb integrals:

      write(nu,'(/12x,'' Two-electron integrals:''/)')
      Do i=1,nzoef; Call Pri_coef1(nu,ipt(i),2,eps_c); End do

! ... two-electron Breit integrals:

      if(mbreit.gt.0) then
       write(nu,'(/12x,'' Breit integrals:''/)')
       Do i=1,nzoef; Call Pri_coef1(nu,ipt(i),3,eps_c); End do
      end if

      write(nu,'(70(''-'')/)')

      End Subroutine Pri_coef
 

!======================================================================
      Subroutine Pri_coef1 (nu,ii,met,eps_c)
!======================================================================
!     Prints one angular integral, ii, from zoef_list
!----------------------------------------------------------------------
      Use orb_jj
      Use zoef_list; Use int_list
      Use ndet_list; Use ndef_list

      Implicit real(8) (A-H,O-Z)

      Integer, intent(in) :: nu,ii,met

      Character(1) :: AI(3);  Data AI /'I','R','S'/

      Character(5) :: EL, EL1,EL2,EL3,EL4

      Character(320) :: A
      Character(1) :: b1,b2,bb

      if(abs(Zoef(ii)).lt.eps_c) Return

      i=IZ_int(ii); int=Jcase(i); if(int.ne.met) Return

      k=Jpol(i);  i1=J1(i); i2=J2(i); i3=J3(i); i4=J4(i)

      A = ' '

      Select case(int)

      Case(0)                            ! Overlap

       ia=0

      Case(1)                            ! One-electron integrals

       EL1 = ELF(I1); EL2 = ELF(I3)
       write(A(1:17),'(a1,a4,2(a5,a1))') AI(int),'   (',EL1,',',EL2,')'
       ia=17

      Case(2,3)                            ! Two-electron integrals

       EL1 = ELF(I1); EL2 = ELF(I2)
       EL3 = ELF(I3); EL4 = ELF(I4)
       write(A(1:29),'(a1,i2,a2,4(a5,a1))') &
             AI(int),k,' (',EL1,',',EL2,';',EL3,',',EL4,')'
       ia=29
 
      Case default

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
         write(A(ia+1:ia+6),'(10(a4,'',''))') EL
         ia=ia+5
        End do
        A(ia:ia)='|'
        Do j=1,nd
         ni=mod(NPD(np+j),ibd); EL = ELF(ni)(2:5)
         write(A(ia+1:ia+6),'(10(a4,'',''))') EL
         ia=ia+5
        End do
        if(ns.gt.1) then
         A(ia:ia+1)='>^'; ia=ia+2; write(A(ia:ia),'(i1)') ns; ia=ia+1
        else
         A(ia:ia+1)='> '; ia=ia+2
        end if
       End do       ! Over kd
      end if

      Call Num (Zoef(ii),k1,k2,999999,1.d-9)

      b1='['; b2=']'; bb=':'
      k1 = iabs(k1); m1=sqrt(1.0*k1)+0.1; m2=sqrt(1.0*k2)+0.1
      if(k1.eq.m1*m1.and.k2.eq.m2*m2) then
       k1=m1;k2=m2;b1='(';b2=')'
      end if
      write(nu,'(f12.6,1x,a1,i6,a1,i6,a1,2x,220a1)') &
                Zoef(ii),b1,iabs(k1),bb,k2,b2,(A(i:i),i=1,ia)


      End Subroutine Pri_coef1


!---------------------------------------------------------------------
      Subroutine prepare_iort_jj
!---------------------------------------------------------------------
!     prepares the orthogonality conditions (IORT) for all orbitals
!     according to the orbital set numbers KEF, and additional
!     conditions from c-file (unit nu)
!     iort(i,j) = 0    -  for orthogonal orbitals i,j
!     iort(i,j) = 2    -  for non-orthogonal orbitals i,j
!     iort(i,j) = 1    -  for the normalized orbitals 
!     OPTIONS:
!     JORT<0 - full orthogonality
!     JORT=0 - full non-orthogonality
!     JORT>0 - partial orthogonality, i.e.
!     the orbitals  are orthogonal within one set with the same KEF,
!     and all orbitals are orthonomal to the orbitals with KEF=0
!---------------------------------------------------------------------
      Use orb_jj

      Implicit none
      Integer :: i1,i2
      
      Do i1=1,nwf
      Do i2=1,i1
       if(KEF(i1).ne.KEF(i2)) Cycle
       if(JORT.lt.0) then             ! full orthogonality
        if(i1.eq.i2) Call Iadd_orth(i1,i2,1)
       elseif(JORT.eq.0) then         ! full non-orthogonality
        Call Iadd_orth(i1,i2,2) 
       elseif(JORT.gt.0) then         ! partial orthogonality
        if(IEF(i1).eq.IEF(i2)) then
         if(NEF(i1).eq.NEF(i2)) Call Iadd_orth(i1,i2,1)          
        elseif(IEF(i1)*IEF(i2).eq.0) then
         Cycle
        else
         Call Iadd_orth(i1,i2,2) 
        end if
       end if
      End do
      End do

      End Subroutine prepare_iort_jj


