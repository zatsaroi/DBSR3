!--------------------------------------------------------------------
!     dbound.nnn  -->  dbound_tab - list of states
!--------------------------------------------------------------------
      Use zconst, only: c_au

      Implicit real(8) (A-H,O-Z)
      
      Integer, allocatable :: ist(:), jst(:), ilsp(:), IPT(:), &
                              JT(:), JP(:), ich_state(:), it_state(:)
      Real(8), allocatable :: E(:), Ebind(:)
      
      Character(64) :: label1
      Character(64), allocatable :: Label(:)
      Character(40) :: AF, name

! ... files:

      Integer :: nu=1;  Character(20) :: AF_res = 'dbound_tab'
      Integer :: in=2;  Character(20) :: AF_inp = 'dbound.nnn'

      Integer, external :: Icheck_file

      Call Read_name(name)

      if(name == '?') then
       write(*,*) 
       write(*,*) 'create list of states (dbound_tab) containing in dbound.nnn files'
       write(*,*) 
       write(*,*) 'see possible parameters in the end of the dbound_tab file' 
       write(*,*) 
       Stop ' '
      end if

!--------------------------------------------------------------------
! ... Check if there is non-default information about input files
! ... and parameters needed au -> eV conversion:

      Open(nu,file=AF_res)

      Z    =  0.d0;  Call Read_rpar(nu,'Z'  ,Z  )
      AWT  =  0.d0;  Call Read_rpar(nu,'AWT',AWT)
      E0   =  0.d0;  Call Read_rpar(nu,'E0' ,E0 )
      Emin =  0.d0;  Call Read_rpar(nu,'Emin',Emin) 
      Emax =  0.d0;  Call Read_rpar(nu,'Emax',Emax)

      if(Z.eq.0.d0) Z = 1.d0
      if(Emin.eq.0.d0) Emin = -2*c_au*c_au

      klsp1= 1;  Call Read_ipar(nu,'klsp1',klsp1)
      klsp2=99;  Call Read_ipar(nu,'klsp2',klsp2)
      klsp3= 1;  Call Read_ipar(nu,'klsp3',klsp3)

      klsp = 0;  Call Read_iarg('klsp',klsp)
      if(klsp.ne.0) then; klsp1=klsp; klsp2=klsp; end if

      Call Conv_au (Z,AWT,au_cm,au_eV,0)       

      if(klsp.ne.0) write(AF_res,'(a,i3.3)') 'dbound_tab_',klsp
      Open(nu,file=AF_res)

!--------------------------------------------------------------------
! ... define the number of states in all bound.nnn files:
                                                         
      nstate=0; nterm=0; ia = Len_trim(AF_inp)-3
      Do klsp=klsp1,klsp2,klsp3
       write(AF_inp(ia+1:ia+3),'(i3.3)') klsp
       if(Icheck_file(AF_inp).eq.0) Cycle
       Open(in,file=AF_inp,form='UNFORMATTED')
       nterm=nterm+1
       read(in) ib,kch,kpert,ns,jot,iparity,nbound
       nstate = nstate+nbound
      End do  

      write(*,*) ' nstate =', nstate, '   nterm = ',nterm

      if(nterm.eq.0.or.nstate.eq.0) Stop 'nothing to collect'

      Allocate(ist(nstate), jst(nstate), E(nstate), Label(nstate), &
               IPT(nstate), JT(nterm),JP(nterm),ilsp(nterm),       &
               ich_state(nstate), it_state(nstate), Ebind(nstate))

!--------------------------------------------------------------------
! ... read states information:              

      iterm=0; istate=0
      Do klsp=klsp1,klsp2,klsp3
       write(AF_inp(ia+1:ia+3),'(i3.3)') klsp
       if(Icheck_file(AF_inp).eq.0) Cycle
       Open(in,file=AF_inp,form='UNFORMATTED')
       rewind(in)
       read(in) ib,kch,kpert,ns,jot,iparity,nbound

       iterm = iterm + 1
       ilsp(iterm) = klsp
       JT(iterm)   = jot
       JP(iterm)   = iparity
       Do i = 1,nbound
        istate = istate + 1
        read(in) jst(istate),Label(istate)
        read(in) E(istate),Ebind(istate),ich_state(istate),it_state(istate)
        read(in) (S,j=1,ib)
        ist(istate) = iterm
       End do
      End do
      if(istate.ne.nstate) Stop 'nstate <> istate'

!--------------------------------------------------------------------
! ... order the energies:

      Call SORTR(nstate,E,IPT)

!--------------------------------------------------------------------
! ... find space for label:
    
      imax=0
      Do i = 1,nstate
       ii = LEN_TRIM(Label(i)); if(ii.gt.imax) imax=ii
      End do

!--------------------------------------------------------------------
! ... find reference energy and excitation (or bindinfg) energies:
    
      if(E0.eq.0.d0) E0 = E(IPT(1))
      
      rewind(nu) ; isf=0

      Do j = 1,nstate;  i=IPT(j)

       if(E(i).gt.Emax) Cycle
       if(E(i).lt.Emin) Cycle

       ee = E(i) - E0;     
       e_eV = ee * au_eV
       e_cm = ee * au_cm
       
       label1 = Label(i)

       js = jst(i); is = ist(i)
 
       write(nu,'(2i5,3x,a,2i3,f10.5,F15.5,f12.1,f20.8,f15.5,2i6)') &
        ilsp(is),js,label1(1:imax),JT(is),JP(is), &
        ee,e_eV,e_cm,E(i), Ebind(i), it_state(i), ich_state(i) 
       isf = isf + 1 
      
      End do 

      write(nu,'(a)') '*'
      write(nu,*)
      write(nu,'(a,f7.3)')  'Z    = ',Z
      write(nu,'(a,f7.3)')  'AWT  = ',AWT
      write(nu,*)
      write(nu,'(a,f18.8)') 'E0   = ',E0
      write(nu,'(a,f18.8)') 'Emin = ',Emin
      write(nu,'(a,f18.8)') 'Emax = ',Emax
      write(nu,*)
      write(nu,'(a,i3)')    'klsp1 = ',klsp1
      write(nu,'(a,i3)')    'klsp2 = ',klsp2
      write(nu,'(a,i3)')    'klsp3 = ',klsp3
      write(nu,*)
      write(nu,'(a,f16.6)') 'au_eV= ',au_eV
      write(nu,'(a,f16.6)') 'au_cm= ',au_cm
      write(nu,'(a,f16.6)') 'cm_ev= ',au_cm/au_ev
      write(nu,*)
      write(nu,'(a,i5)')    'nstate = ',isf

      
      End  !  program dbound_tab

      