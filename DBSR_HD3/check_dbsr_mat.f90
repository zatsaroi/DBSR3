!=====================================================================
      Subroutine Check_dbsr_mat
!=====================================================================
!     check and print main parameters and input matrix
!     for the given partial wave
!---------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: i, i1,i2,i3
      Integer, external :: Icheck_file

!----------------------------------------------------------------------
! ... open log file and print main parameters:

      write(ALSP,'(i3.3)') klsp
      i = INDEX(AF_nnn,'.',BACK = .TRUE.); AF=AF_nnn(1:i)//ALSP
      Open(pri,file=AF)

      write(pri,'(a/a)') 'D B S R _ H D ', &
                         '*************'
      write(pri,*)
      write(pri,'(a,i3)') 'calculations for partial wave:  klsp =',klsp

      write(pri,*)
	  if(itype.eq.-1) &
       write(pri,'(a)') 'itype  =   -1  -  bound-state calculations'
      if(itype.eq. 0) &
       write(pri,'(a)') 'itype  =    0  -  scattering calculations'
      if(itype.gt. 0) &
       write(pri,'(a)') 'itype  =    1  -  photoionization calculations'

      mhm = ms*nch+npert
      write(pri,*)
      write(pri,'(a,i5,a)') 'nch    =',nch,    '  -  number of channels'
      write(pri,'(a,i5,a)') 'npert  =',npert,  '  -  number of pertubers'

      write(pri,*)
      if(iexp.eq.0) &
      write(pri,'(a)') 'iexp   =    0  -  theoretical target energies'
      if(iexp.ne.0) &
      write(pri,'(a)') 'iexp   =     1  -  exp.target energies'

      if(itype.eq.-1) then
       write(pri,*)
       write(pri,'(a)') 'Restrictions for output of bound states:'
       write(pri,*)
       write(pri,'(a)') '(zero value of any parameter means no restrictions)'
       write(pri,*)
       write(pri,'(a,i5,a)') 'msol   =',msol,  '  -  max. number of solutions'
       write(pri,'(a,e15.8)')'Emin   =',Emin  
       write(pri,'(a,e15.8)')'Emax   =',Emax
       write(pri,*)
      end if

!----------------------------------------------------------------------
! ... check the file with interaction matrix:

      i = INDEX(AF_int,'.'); AF=AF_int(1:i)//ALSP
      i = Icheck_file(AF)
      if(i.eq.0) then
       write(pri,*) 'there is no dbsr_mat.nnn file for given partial wave'
       Stop         'there is no dbsr_mat.nnn file for given partial wave'
      end if

      Open(nui,file=AF,status='OLD',form='UNFORMATTED')

! ... check the dimensions:

      read(nui) i1,i2,i3
      if(i1.ne.ns )   Stop ' DBSR_HD: different ns  in DBSR_MAT file'
      if(i2.ne.nch)   Stop ' DBSR_HD: different kch in DBSR_MAT file'
      if(i3.ne.npert) Stop ' DBSR_HD: different kcp in DBSR_MAT file'

! ... allocate common working arrays:

      read(nui) nsol  !,jtype
      if(allocated(ipsol)) Deallocate(ipsol); Allocate(ipsol(0:nch))
      read(nui) ipsol
      if(allocated(bval)) Deallocate(bval); Allocate(bval(nsol))
      read(nui) bval
      if(allocated(bb)) Deallocate(bb); Allocate(bb(ms,nsol))
      read(nui) bb

      if(allocated(v)) Deallocate(v);  Allocate(v(mhm)); v = zero
      khm = nsol+npert; ksol = nsol

      End Subroutine Check_dbsr_mat
