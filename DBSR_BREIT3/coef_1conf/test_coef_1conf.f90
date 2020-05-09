!=====================================================================
!     PROGRAM   test_coef_1conf                      
!
!               C O P Y R I G H T -- 2010
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!
!    it is a debug program to check subroutine     "coef_1conf" 
!    to generate angular coefficient for Dirak_Fock calculations 
!    for one-configuration
!    in case of orthogonal one-electron radial functions
!
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:  AF_cfg=...  AF_out=...  
!    
!----------------------------------------------------------------------
!
!    INPUT FILE:     AF_cfg    (default  - rcls.inp) 
!    OUTPUT FILES:   AF_out    (default  - coef.tab)
!    
!---------------------------------------------------------------------     

      Implicit none 

      Integer :: nuc=1; Character(20) :: AF_cfg = 'rcsl.inp'
      Integer :: out=2; Character(20) :: AF_out = 'coef.tab'
      Character(40) :: name = ' '

! ... description of 1 conf.w.function:  

      Integer, parameter :: msh = 31 ! max. number of shells behind core
      Integer :: no
      Integer, dimension(msh) :: nn,kn,ln,jn,iq,in,Jshell,Vshell,Jintra

! ... Storing configurations as character strings:

      Character(9*msh+9) :: CONFIG, SHELLJ, INTRAJ

! ... result coefficients:

      Real(8), allocatable :: coefs(:,:,:)

! ... local variables:

      Real(8) :: C, eps_C = 1.d-7
      Integer :: i,j, k, ic, ncfg, kmax
      Integer, External :: Jdef_ncfg 

!----------------------------------------------------------------------
      Call Read_name(name)
      if(len_trim(name).ne.0) then
       AF_cfg=trim(name)//'.c'
       AF_out=trim(name)//'.tab'
      else
       Call Read_aarg('inp',AF_cfg)
       Call Read_aarg('out',AF_out)
      end if

      Call Check_file(AF_cfg)
      Open(nuc,file=AF_cfg)
      ncfg = Jdef_ncfg(nuc)
      if(ncfg.eq.0) Stop 'ncfg=0: nothing to do'
      Open(out,file=AF_out)

      rewind(nuc)
      Do ic = 1,ncfg

      Do 
       read(nuc,'(a)') CONFIG
       if(CONFIG(6:6).ne.'(') Cycle
       read(nuc,'(a)') SHELLJ
       read(nuc,'(5x,a)') INTRAJ
       Exit
      End do

      Call Decode_cjj(CONFIG,SHELLJ,INTRAJ,no,nn,kn,ln,jn,iq,in,&
                                            Jshell,Vshell,Jintra)
      kmax = maxval(jn(1:no))
      
      if(allocated(coefs)) Deallocate(coefs)
      Allocate(coefs(no,no,0:kmax))

      Call coef_1conf(no,ln,jn,iq,Jshell,Vshell,Jintra,kmax,coefs)
 
      write(out,'(a)') trim(config) 
      write(out,'(a)') trim(SHELLJ)
      write(out,'(9x,a)') trim(INTRAJ)
 
       write(out,*)
       Do i=1,no
	   Do j=1,i
        Do k=0,kmax
		 C = coefs(i,j,k); if(abs(C).lt.eps_c) C=0.d0
		 if(C.eq.0.d0) Cycle
         if(mod(ln(i)+ln(j)+k,2).ne.0) Cycle
         write(out,'(a,i2,a,i2,a,i2,a,f10.5)') &
		            'F',k,'(',i,',',j,')=',C        
        End do
        if(i.eq.j) Cycle
        Do k=0,kmax
		 C = coefs(j,i,k); if(abs(C).lt.eps_c) C=0.d0
		 if(C.eq.0.d0) Cycle
         if(mod(ln(i)+ln(j)+k,2).ne.0) Cycle
         write(out,'(a,i2,a,i2,a,i2,a,f10.5)') &
		            'G',k,'(',i,',',j,')=',C        
        End do
       End do; End do
       write(out,*)

      End do  ! over ic

      END ! Program ...




