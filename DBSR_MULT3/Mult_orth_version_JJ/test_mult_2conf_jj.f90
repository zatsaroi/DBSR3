!======================================================================
!     PROGRAM   test_mult_2conf_jj                      
!
!               C O P Y R I G H T -- 2020
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    it is a debug program to check subroutine     "mult_2conf_jj" 
!    to generate angular coefficient for multipole one-electron operator
!    in case of orthogonal one-electron radial functions
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:    name1.c  name2.c   E1|M1|E2.... 
!
!----------------------------------------------------------------------
!
!    INPUT FILE:     name1.c,  name2.c  
!    OUTPUT FILES:   tab_name1_name2
!    
!----------------------------------------------------------------------     
      Implicit real(8) (A-H,O-Z) 

      Integer :: nuc1=1; Character(80) :: AF_cfg1 
      Integer :: nuc3=2; Character(80) :: AF_cfg2 
      Integer :: out =3; Character(80) :: AF_tab 

      Character(2) :: atype = 'E1'
      Character :: ktype = 'E'
      Integer :: kpol = 1      !  multipole index 

      Integer, parameter :: msh = 31 ! max. number of shells behind core

      Integer :: no1
      Integer, dimension(msh) :: nn1,kn1,ln1,jn1,iq1,in1,  &
                                 Jshell1,Vshell1,Jintra1
      Integer :: no2
      Integer, dimension(msh) :: nn2,kn2,ln2,jn2,iq2,in2,  &
                                 Jshell2,Vshell2,Jintra2

      Integer, parameter :: mas = 9*msh+3
      Character(mas) :: CONFIG1, SHELLJ1, INTRAJ1
      Character(mas) :: CONFIG2, SHELLJ2, INTRAJ2

      Character(5) :: EL1, EL2

      Integer, parameter :: mk = 100
      Integer :: ik(mk), jk(mk)
      Real(8) :: ck(mk)

!----------------------------------------------------------------------
! ... input data:     

      i = command_argument_count()
      if(i.lt.3) then
       write(*,*)
       write(*,*) 'Should be at least 3 command arguments:  name1.c nam12.c  E1|E2|...'
       write(*,*)
       Stop ' '
      end if

      Call GET_COMMAND_ARGUMENT(1,AF_cfg1)
      Call GET_COMMAND_ARGUMENT(2,AF_cfg2)
      Call GET_COMMAND_ARGUMENT(3,atype)
      read(atype,'(a1,i1)') ktype, kpol

      Call Check_file(AF_cfg1)
      Open(nuc1,file=AF_cfg1)
      ncfg1 = Jdef_ncfg(nuc1)

      Call Check_file(AF_cfg2)
      Open(nuc2,file=AF_cfg2)
      ncfg2 = Jdef_ncfg(nuc2)

      i1 = Len_trim(AF_cfg1)-2
      i2 = Len_trim(AF_cfg2)-2
      write(AF_tab,'(a,a,a,a)') 'tab_',AF_cfg1(1:i1),'_',AF_cfg2(1:i2)
      Call Read_aarg('tab',AF_tab)
      Open(out,file=AF_tab)

!-----------------------------------------------------------------------

      rewind(nuc1); i1=0
      Do 
       read(nuc1,'(a)') CONFIG1
       if(CONFIG1(6:6).ne.'(') Cycle
       read(nuc1,'(a)') SHELLJ1
       read(nuc1,'(5x,a)') INTRAJ1
       i1 = i1 + 1

       Call Decode_cjj(CONFIG1,SHELLJ1,INTRAJ1,no1,nn1,kn1,ln1,jn1,iq1,in1, &
                       Jshell1,Vshell1,Jintra1)

      rewind(nuc2); i2=0
      Do 
       read(nuc2,'(a)') CONFIG2
       if(CONFIG2(6:6).ne.'(') Cycle

       read(nuc2,'(a)') SHELLJ2
       read(nuc2,'(5x,a)') INTRAJ2
       i2 = i2 + 1

       Call Decode_cjj(CONFIG2,SHELLJ2,INTRAJ2,no2,nn2,kn2,ln2,jn2,iq2,in2, &
                       Jshell2,Vshell2,Jintra2)


       Call Mult_2conf_jj(no1,nn1,ln1,jn1,iq1,Jshell1,Vshell1,Jintra1, &
                          no2,nn2,ln2,jn2,iq2,Jshell2,Vshell2,Jintra2, &
                          atype,nk,ck,ik,jk )

       if(nk.ne.0) then
  
        write(out,'(80("="))') 
        write(out,'(a,T70,i5)') trim(CONFIG1),i1
        write(out,'(a)') trim(SHELLJ1)
        write(out,'(5x,a)') trim(INTRAJ1)
  
        write(out,'(a,T70,i5)') trim(CONFIG2),i2
        write(out,'(a)') trim(SHELLJ2)
        write(out,'(5x,a)') trim(INTRAJ2)
  
        write(out,'(80("-"))') 
  
        Do k = 1,nk
         i = (ik(k)-1)*9; EL1=CONFIG1(i+1:i+5)
         j = (jk(k)-1)*9; EL2=CONFIG2(j+1:j+5)
         write(out,'(a,i1,1x,a,a,a,a,a,f10.5)') &
                'd',kpol,'(',EL1,',',EL2,')=',Ck(k)        
        End do
        write(out,*)

       end if


       if(i2.eq.ncfg2) Exit
      End do  ! second configuration
       if(i1.eq.ncfg1) Exit
      End do  ! first configuration


      End  ! program ...


!======================================================================
      Subroutine Decode_cjj(CONFIG,SHELLJ,INTRAJ,no,nn,kn,ln,jn,iq,in, &
                            Jshell,Vshell,Jintra) 
!======================================================================
!     decode the configuration from c-file format to INTEGER format
!----------------------------------------------------------------------
      Implicit none

      Integer, parameter :: msh = 31 ! max. number of shells behind core
      Integer :: no
      Integer, dimension(msh) :: nn,kn,ln,jn,iq,in,  &
                                 Jshell,Vshell,Jintra
      Integer, parameter :: mas = 9*msh+3
      Character(mas) :: CONFIG, SHELLJ, INTRAJ

      Integer :: i,j,k,m
      Character(9) :: bl = '        ' 

      i=INDEX(CONFIG,')',BACK=.TRUE.); no=i/9

      Vshell=0; Jshell=0; Jintra=0
      m=-9
      Do i=1,no
       m = m + 9
       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)
 
       if(SHELLJ(m+5:m+5).eq.';') then
        read(SHELLJ(m+4:m+4),'(i1)') Vshell(i)
        read(SHELLJ(m+9:m+9),'(i1)') J; Jshell(i) = 2*J
       elseif(SHELLJ(m+8:m+8).eq.'/') then
        read(SHELLJ(m+1:m+7),'(i7)') Jshell(i)
       elseif(SHELLJ(m+9:m+9).ne.' ') then
        read(SHELLJ(m+1:m+9),'(i9)') J; Jshell(i) = 2*J
       end if       

       k = INDEX(INTRAJ(m+1:m+9),'/')
       if(i.eq.1) then
        Jintra(i) = Jshell(i)
       elseif(i.eq.no) then
        Cycle
       elseif(k.gt.0) then
        read(INTRAJ(m+1:m+k-1),*) Jintra(i)
       elseif(INTRAJ(m+1:m+9).ne.bl) then
        read(INTRAJ(m+1:m+9),*) J; Jintra(i) = 2*J
       else
       if(Jshell(i).eq.0) Jintra(i) = Jintra(i-1)              
       if(Jintra(i-1).eq.0) Jintra(i) = Jshell(i)
       end if       

      End do

      if(k.gt.0) then
       read(INTRAJ(m+1:m+k-1),*) Jintra(no)
      else
       k=INDEX(INTRAJ(m+1:m+9),'+')
       if(k.eq.0) k=INDEX(INTRAJ(m+1:m+9),'-')
       read(INTRAJ(m+1:m+k-1),*) J; Jintra(no) = 2*J
      end if

      Do i=1,no
       m = jn(i)+1
       if(Vshell(i).ne.0.or.iq(i).eq.m) Cycle
       if(iq(i).eq.1.or.iq(i).eq.m-1) then
         Vshell(i) = 1
       elseif(iq(i).eq.2.or.iq(i).eq.m-2) then
         if(Jshell(i).gt.0) Vshell(i) = 2
       elseif(iq(i).eq.3) then
         Vshell(i) = 3; if(Jshell(i).eq.jn(i)) Vshell(i)=1
       end if
      End do

 
      End Subroutine Decode_cjj

!=======================================================================
      Subroutine EL_nljk(EL,n,kappa,l,j,k)
!=======================================================================
!     decodes the spectroscopic notation for electron orbital (n,l,j,k)
!
!     It is allowed following notations: 1s , 2s 3, 2p-30, 20p-3,
!     1s h, 20s h, kp , kp 1, kp-11, ns , ns 3, ns 33, ... .
!
!     Call:  LA, INDEX
!----------------------------------------------------------------------

      Implicit none

      Character(5), intent(in) :: EL
      Integer, intent(out) :: n,l,j,k,kappa    
      Integer :: jj, k1,k2, n1,n2  
      Integer, external :: LA, kappa_lj

      Integer, parameter :: kset = 61
      Character(kset) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      jj=0
      Do j=5,3,-1
       if(EL(j:j).eq.'-'.or.EL(j:j).eq.'+') then; jj=j; Exit; end if
      End do
      if(jj.eq.0) then
       Do j=5,3,-1
        if(EL(j:j).eq.' '.and.EL(j-1:j-1).ne.' ') then
         jj=j; Exit
        end if
       End do
      end if

      if(jj.eq.5) then

       k = 0
       l = LA(EL(4:4))
       n1 = INDEX(ASET,EL(3:3))
       n2 = INDEX(ASET,EL(2:2))

      elseif(jj.eq.4) then

       k = INDEX(ASET,EL(5:5))
       l = LA(EL(3:3))
       n1 = INDEX(ASET,EL(2:2))
       n2 = INDEX(ASET,EL(1:1))

      elseif(jj.eq.3) then

       k1 = INDEX(ASET,EL(5:5))
       k2 = INDEX(ASET,EL(4:4))
       k = k2*kset+k1
       l = LA(EL(2:2))
       n1 = INDEX(ASET,EL(1:1))
       n2 = 0

      else

       write(*,*) 'EL_NLJK: can not decode ',EL
       Stop ' '

      end if

       if(n1.le.9.and.n2.le.9) then
        n = n2*10+n1
       elseif(n2.eq.0.and.n1.eq.20) then
        n=ICHAR('k')
       elseif(n2.eq.0.and.n1.eq.23) then
        n=ICHAR('n')
       else
        n = n2*kset+n1
       end if


       j = l+l+1; if(EL(jj:jj).eq.'-') j = l+l-1
       kappa = kappa_lj(l,j)

      End Subroutine EL_NLJK


!======================================================================
      Integer Function Jdef_ncfg(nu)
!======================================================================
!     defines the number of configuration in c-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: ncfg
      Character(6) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Jdef_ncfg=ncfg

      End Function Jdef_ncfg



!======================================================================
      Subroutine Check_file(AF)
!======================================================================

      Character(*), Intent(in) :: AF
      Logical :: EX

      Inquire(file=AF,exist=EX)
      if(.not.EX) then
       write(*,*) ' can not find file  ',AF;  Stop
      end if

      End Subroutine Check_file



!====================================================================
      Integer FUNCTION LA(a)
!====================================================================
!     gives the value of L from spetroscopic symbol "a"
!--------------------------------------------------------------------
      Implicit none  
      Character, Intent(in) :: a
      Character(21) :: AS, AB
      Integer :: i

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = INDEX (AS,a)
      if(i.eq.0) i = INDEX (AB,a)
      if(i.eq.0) i = ICHAR(a)-101
      LA = i-1

      End FUNCTION LA


!=======================================================================
      Integer Function kappa_lj(l,jj)
!=======================================================================
!     l, j  ->  kappa = (l-j)*(2j+1);  jj = 2j
!-----------------------------------------------------------------------
      Integer :: l,jj
      kappa_lj = (2*l-jj)*(jj+1)/2
      End Function kappa_lj


!======================================================================
      Subroutine Read_aarg(name,avalue)
!======================================================================
!     read characer argument
!----------------------------------------------------------------------
      Implicit None

      Character(*) :: name, avalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS

      iarg = Command_argument_count(); if(iarg.eq.0) Return 
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call get_command_argument(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),'(a)') avalue; Exit
      End do

      End Subroutine Read_aarg

