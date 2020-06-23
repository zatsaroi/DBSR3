!======================================================================
!     utility     G E N J C O N F
!
!                 C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny
!=======================================================================
!     preperation of list of possible configurations from list of
!     electron occupations:
!
!                nlj.inp  -->  conf.inp
!
!----------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Integer, allocatable :: iq_min(:),iq_max(:)

      Character(40) :: AF = ' '
      Integer :: nu1=1;  Character(40) :: AF_inp = 'nlj.inp'
      Integer :: nu2=2;  Character(40) :: AF_out = 'conf.inp'

      Integer, external :: Icheck_file

      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)
      if(AF.eq.'?') then
      write(*,*) 
      write(*,*) 'genjconf prepares the list of possible configurations (conf.inp)'
      write(*,*) 'from the list of electron occupations and allowed excitations (nlj.inp)'
      write(*,*)
      write(*,*) 'Call as:   genjconf  inp=...  out=...    or  atom=...'
      write(*,*)
      write(*,*) 'All arguments are optional, with defaults inp=nlj.inp, out=conf.inp'
      write(*,*)
      write(*,*) 'If nlj.inp file is absent, it will be creared for you to fill in'
      write(*,*)
!      write(*,*) 'iorder parameter (then =1) causes additional ordering '
!      write(*,*) 'according to "non-relativistic" projections '
!      write(*,*)
      write(*,*) 'Parameter "atom"  provides  nlj.inp for the given atom (as example) '
      Stop 
      end if

      Call Read_aarg('inp',AF_inp)
      Call Read_aarg('out',AF_out)
      iorder=0; Call Read_iarg('iorder',iorder)

      AF = ' '; Call Read_aarg('atom',AF)
      if(Icheck_file(AF_inp).eq.0.or.len_trim(AF).ne.0) then
       Open(nu1,file=AF_inp)
       Call Write_nlj(nu1) 
       write(*,'(/a/)') 'nlj.inp - use this file to correct your case '
      end if

!----------------------------------------------------------------------

      open(nu1,file=AF_inp,status='OLD')
      open(nu2,file=AF_out)

      Call Read_ipar(nu1,'n_orbitals' ,nwf)   
      Call Read_ipar(nu1,'n_electrons',ne)    
      Call Read_ipar(nu1,'parity'     ,parity)
      Call Read_ipar(nu1,'k_ref'      ,k_ref) 
      Call Read_ipar(nu1,'k_min'      ,k_min) 

      Call Read_ipar(nu1,'k_max'      ,k_max) 


! ... read and write HEADER  and CLOSED

      rewind(nu1)
      read(nu1,'(a)') AS; ia=LEN_TRIM(AS);  write(nu2,'(a)') AS(1:ia)
      read(nu1,'(a)') AS; ia=LEN_TRIM(AS);  write(nu2,'(a)') AS(1:ia)
      read(nu1,'(a)') AS; ia=LEN_TRIM(AS);  write(nu2,'(a)') AS(1:ia)

! ... read orbitals

      Call alloc_orb_jj(nwf)
      Allocate(iq_min(nwf),iq_max(nwf))

      read(nu1,'(50a5)') (ELF(i),i=1,nwf)
      read(nu1,*)
      read(nu1,'(50i5)') (iq_min(i),i=1,nwf)
      read(nu1,'(50i5)') (iq_max(i),i=1,nwf)
      Close(nu1)

      Do i=1,nwf
       Call EL_nljk(ELF(i),NEF(i),KEF(i),l,j,IEF(i))
      End do
      write(nu2,'(50a5)') (ELF(i),i=1,nwf)
      write(nu2,'(a)') 'CSF(s):'

! ... generation of all possible configurations:

      ncfg=0

      Call Sum_conf(iq_min,iq_max,k_ref,k_min,k_max,nu2)

      write(nu2,'(a)') '***'

      write(*,'(3a,i6/)')  'nlj.inp --> ',trim(AF_out),'  nconf =',ncfg

      iorder = 1;  Call Read_iarg('order',iorder)

      if(iorder.eq.10) Call order_conf(nu2)     !!!   not working

      End  ! end utility GENGCONF


!----------------------------------------------------------------------
      Subroutine Sum_conf(iq_min,iq_max,k_ref,k_min,k_max,nu)
!----------------------------------------------------------------------
!
!     exhaustion of possible configurations
!
!----------------------------------------------------------------------

      USE conf_jj; USE orb_jj

      Integer(4) :: iq_min(*),iq_max(*),k_ref,k_min,k_max,nu

      i1=1; i2=nwf; ipef=0 

      i=1;  ipef(i)=iq_max(i)

    1 ii=SUM(ipef(1:i))
      if(ii.eq.ne) then
        Call gen_conf(k_ref,k_min,k_max,nu)
      elseif(ii.lt.ne.and.i.lt.i2) then
        i=i+1
        ipef(i)=iq_max(i)
        go to 1
      end if

    2 ipef(i)=ipef(i)-1
      if(ipef(i).lt.iq_min(i)) then
        if(i.eq.i1) Return
        i=i-1
        go to 2
      end if
      go to 1

      End  ! Subroutine Sum_conf


!====================================================================
      Subroutine Gen_conf(k_ref,k_min,k_max,nu)
!====================================================================
!     generate and check 1 configuration with electron population
!     given by IEF in module configs
!--------------------------------------------------------------------
      USE conf_jj;  Use orb_jj

      Integer, intent(in) :: k_ref,k_min,k_max,nu
      Character(5), external :: ELi

      no=0; n_corr=0
      Do i=1,nwf
       if(ipef(i).le.0) Cycle
       if(i.gt.k_ref) n_corr=n_corr+ipef(i)
       no=no+1
       nn(no)=NEF(i); kn(no)=KEF(i); in(no)=IEF(i); iq(no)=ipef(i)
       ln(no)=l_kappa(KEF(i))
       jn(no)=j_kappa(KEF(i))
      End do

      if(n_corr.lt.k_min) Return
      if(n_corr.gt.k_max) Return

! ... check the total parity:

      k=0;  Do i=1,no; k=k+iq(i)*ln(i);  End do
      k=(-1)**k
      if(k.ne.parity) Return

! ... check the number of electrons in the shell

      Do i=1,no;  if(iq(i).gt.jn(i)+1) Return;  End do

! ... record configuration:

      m=0
      Do i=1,no
       CONFIG(m+1:m+5) = ELi(nn(i),kn(i),in(i))
       write(CONFIG(m+6:m+9),'(a1,i2,a1)') '(',iq(i),')'
       m = m + 9
      End do
      write(nu,'(a)') CONFIG(1:m)
   
      ncfg=ncfg+1

      End ! Subroutine Gen_conf


!======================================================================
      Subroutine order_conf(nuc)
!======================================================================
!  ... ordering according to non-relativistic projections
!----------------------------------------------------------------------
      Use conf_jj, ic => IT_state1, ip => IT_state2
      Use orb_jj
      
      Implicit none
      Integer :: nuc, k, i1,i2
      Integer, External :: Icomp_config

      ncfg=0;  Call Read_config_jj(nuc); if(ncfg.eq.0) Return

      Allocate(ic(ncfg),ip(ncfg))
      ic = 0; k = 0
      Do i1=1,ncfg
       if(ic(i1).ne.0) Cycle        
       k=k+1; ic(i1)=k
       Do i2=i1+1,ncfg
         if(ic(i2).ne.0) Cycle
         if(Icomp_config(i1,i2).eq.0) Cycle
         ic(i2)=k
       End do
      End do
      Call SORTI(ncfg,ic,ip)

      Call Ifind_position(nuc,'Peel')
      read(nuc,*)
      write(nuc,'(20a5)') (ELF(i1),i1=ncore+1,nwf)
      write(nuc,'(a)') 'CSF:'

      Do i1=1,ncfg; i2=ip(i1)
       Call Get_cfg_jj(i2)
       Call Incode_confj
       write(nuc,'(a)') trim(CONFIG)
      End do
      write(nuc,'(a)') '***'

      End Subroutine order_conf

!======================================================================
      Integer function Icomp_config(i1,i2)
!======================================================================
! ... compare non-relativistic projections:
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: i1,i2,ic1,ic2,it1,it2,i

      Icomp_config = 0

      Call Get_cfg_jj(i1)
      no1=no; nn1=nn; ln1=ln; iq1=iq; in1=in
      Call Get_cfg_jj(i2)
      no2=no; nn2=nn; ln2=ln; iq2=iq; in2=in

      if(no1.ne.no2) Return

      Icomp_config=1
      Do i=1,no1
       if(nn1(i).ne.nn2(i)) Icomp_config=0
       if(ln1(i).ne.ln2(i)) Icomp_config=0
       if(in1(i).ne.in2(i)) Icomp_config=0
       if(iq1(i).ne.iq2(i)) Icomp_config=0
       if(Icomp_config.eq.0) Exit
      End do

      End function Icomp_config


!======================================================================
      Subroutine Write_nlj(nu) 
!======================================================================
!     create example of "nlj.inp" file
!----------------------------------------------------------------------
      Use conf_jj

      Integer, intent(in) :: nu
      Integer :: an
      Character(2) :: atom
      Character(mas) :: conf
      Character(5) ::  ELi

      atom = 'Fe'; Call Read_aarg('atom',atom)
      an = 0;      Call Read_iarg('an',an)
      Call Def_atom(an,atom,atw,rms,core,conf)

      rewind(nu)
      write(nu,'(a)') 'Core subshells:'
      Call Decode_core(core)
      write(nu,'(50a5)') e_core(1:ncore)

      write(nu,'(a)') 'Peel subshells:'
      Call Decode_configuration(conf)
      Call Reduce_jj_LS(no,nn,kn,ln,jn,iq,in)
      ii = SUM(iq(1:no))
      jj = (-1)**SUM(iq(1:no)*ln(1:no))
      Call Reduce_LS_jj(no,nn,kn,ln,jn,iq,in,iq1,iq2)
      write(nu,'(50a5)') (ELi(nn(i),kn(i),0),i=1,no)

      write(nu,'(a)') 'Occupation limits:'
      write(nu,'(50i5)') (iq1(i),i=1,no)
      write(nu,'(50i5)') (iq2(i),i=1,no)

      write(nu,*)
      write(nu,*)
      write(nu,'(a,i3,T25,a)') 'n_orbitals  =',no, '! number of peel subshels'
      write(nu,'(a,i3,T25,a)') 'n_electrons =',ii, '! number of peel electrons'
      write(nu,'(a,i3,T25,a)') 'parity =', jj,     '! +1 or -1'
      write(nu,'(a,i3,T25,a)') 'k_ref  =', no,     '! reference orbitals'
      write(nu,'(a,i3,T25,a)') 'k_min  =', 0,      '! min. promotion from the ref. orbitals'
      write(nu,'(a,i3,T25,a)') 'k_max  =', 2,      '! max. promotion from the ref. orbitals'

      End Subroutine Write_nlj


