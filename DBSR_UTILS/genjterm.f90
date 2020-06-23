!=====================================================================
!     utility     G E N J T E R M
!
!                 C O P Y R I G H T -- 2006
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!
!     prepares the ASF list from a list of (nlj)^q configurations
!
!     Calling:  genjterm  inp=file1  out=file2  jmin=...  jmax=... 
!
!     Default:  inp -> rcsl.conf 
!               out -> rcsl.inp
!
!--------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Character(40) :: AF, name=' '

      Integer :: nu1=1; Character(40) :: AF_inp = 'conf.inp'
      Integer :: nu2=2; Character(40) :: AF_out = 'cfg.inp'

      Integer :: jmin=-1, jmax=-1

      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
       write(*,*) 'GENJCONF: prepares the ASF list from a list of (nlj)^q configurations'
       write(*,*) 
       write(*,*) 'Call as:  genjterm inp=.. out=.. jmin=.. jmax=..'
       write(*,*) 'or'
       write(*,*) 'Call as:  genjterm name jmin=.. jmax=..'
       write(*,*)
       write(*,*) 'default input and output files: conf.inp and cfg.inp'
       write(*,*) 'or                              name.conf and name.c  '
       write(*,*) 
       write(*,*) 'jmin is the only compalsory parameter'
       write(*,*) 'if jmax < jmin, jmax = jmin'
       Stop ' '
      end if      


! ... files:

      Call Read_name(name)
      if(len_trim(name).ne.0) then
       AF_inp = trim(name)//'.conf'
       AF_out = trim(name)//'.c'
      else
       Call Read_aarg('inp',AF_inp)
       Call Read_aarg('out',AF_out)
      end if
   
! ... input data:

      Call Read_iarg('jmin',jmin)
      Call Read_iarg('jmax',jmax)
      if(jmin.gt.jmax) jmax=jmin

! ... not enough parameters:

      if(name.eq.'?') then
       write(*,*)
       write(*,*) 'Call as:  genjterm inp=... out=... jmin=... jmax=...'
       write(*,*) 'or'
       write(*,*) 'Call as:  genjterm name jmin=... jmax=...'
       write(*,*)
       write(*,*) 'default input and output files: conf.inp and cfg.inp'
       write(*,*) 'or                              name.conf and name.c  '
       write(*,*) 
       write(*,*) 'jmin is the only compalsory parameter'
       write(*,*) 'if jmax < jmin, jmax = jmin'
       Stop
      end if

      open(nu1,file=AF_inp,status='OLD')
      open(nu2,file=AF_out)

!--------------------------------------------------------------------
! ... read the configurations:

      rewind(nu1);  rewind(nu2)

      read(nu1,'(a)') AS; ia=LEN_TRIM(AS); write(nu2,'(a)') AS(1:ia)
      read(nu1,'(a)') AS; ia=LEN_TRIM(AS); write(nu2,'(a)') AS(1:ia)
      read(nu1,'(a)') AS; ia=LEN_TRIM(AS); write(nu2,'(a)') AS(1:ia)
      read(nu1,'(a)') AS; ia=LEN_TRIM(AS); write(nu2,'(a)') AS(1:ia)
      read(nu1,'(a)') AS;                  write(nu2,'(a)') 'CSF(s):'

      if(jmin.eq.-1) then
       ne = Jdef_ne(nu1)
       jmin = 0; if(mod(ne,2).eq.1) jmin=1
       jmax = 2*ne+1; if(mod(ne,2).eq.0) jmax=2*ne
     end if

      Do j = jmin,jmax,2;   J_min=j; J_max=j

      ncfg=0
    1 read(nu1,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1

      CONFIG = AS;  Call Decode_confj

      Call Sum_Term

      go to 1
    2 Rewind(nu1)

      if(ncfg.eq.0) Cycle
      Do ic=1,ncfg
       Call Get_cfg_jj(ic)
       Call Incode_cj
       ii = no*9
       write(nu2,'(a)') CONFIG(1:ii)
       write(nu2,'(a)') SHELLJ(1:ii)
       write(nu2,'(9x,a)') INTRAJ(1:ii)
      End do
      if(j.lt.jmax) write(nu2,'(a)') ' *'
      if(j.eq.jmax) write(nu2,'(a)') '* '


      write(*,'(2a,i3,a,i6)')  trim(AF_out),'  -->  2J = ',j,'  ncfg =',ncfg

      End do ! over j

      BACKSPACE(nu2)
      write(nu2,'(a)') '* '
      
 
      End  ! utility     G E N J T E R M


!----------------------------------------------------------------------
      Subroutine Sum_Term
!----------------------------------------------------------------------
!     exhaustion of shell-terms
!----------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Implicit none
      Integer :: mt(msh),nt(msh)
      Integer :: i,ii,i1,i2, JT,JV,JW,JQ
      Integer, external :: Jterm

!     mt(i) - the number of term in shell i
!     nt(i) - the term inder consideration

      i1=1                     ! i1 - low  limit of shells
      i2=no                    ! i2 - high limit of shells
      Do i=i1,i2
       mt(i)=Jterm(jn(i),iq(i),-1,JT,JV,JW,JQ)
      End do

      i=i1                     ! first shell under consideration
      nt = 1

    1 Continue

      ii = Jterm(jn(i),iq(i),nt(i),Jshell(i),Vshell(i),JW,JQ)

      if(i.lt.i2) then
         i=i+1; nt(i)=1; go to 1
      else
         CALL Sum_Iterm
      end if
                                         
    2 nt(i)=nt(i)+1
      if(nt(i).gt.mt(i)) then
        if(i.eq.i1) go to 3
        i=i-1; go to 2
        end if
      go to 1

    3 Continue

      End  ! Subroutine Sum_Term


!----------------------------------------------------------------------
      Subroutine Sum_Iterm
!----------------------------------------------------------------------
!     exhaustion of intermediate terms
!----------------------------------------------------------------------
      USE conf_jj; Use orb_jj
 
      Integer :: js_min(msh),js_max(msh)

      Jintra(1)=Jshell(1)
      if(no.eq.1) then
       if(Jshell(no).ge.J_min.and.Jshell(no).le.J_max)  ic=Iadd_cfg_jj('detect')
       Return
      end if

      i1=2                         ! i1 - low  limit
      i2=no                        ! i2 - high limit in array LS(...)

      i=i1
    1 j1=i-1; j2=i

      js_min(i)=IABS(Jintra(j1)-Jshell(j2))
      js_max(i)=     Jintra(j1)+Jshell(j2) 
      Jintra(i)=js_min(i)

    2 if(i.lt.i2) then
       i=i+1; go to 1
      else

       if(Jintra(no).ge.J_min.and.Jintra(no).le.J_max)  ic=Iadd_cfg_jj('detect')

      end if

    3 if(Jintra(i).lt.js_max(i)) then
         Jintra(i)=Jintra(i)+2
         go to 2
      else
         if(i.eq.i1) go to 4
         i=i-1; go to 3
      end if

    4 Continue
       
      End  ! Subroutine Sum_Iterm


