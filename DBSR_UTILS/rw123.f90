!=======================================================================
!     merging the set of GRASP w-files with choice of orbitals and 
!     with optional changing the spectroscopic notation
! 
!       1.w + 2.w + 3.w + ... --> res.w 
!
!     Interactive input
!=======================================================================

      Use orb_jj

      IMPLICIT REAL(8) (A-H, O-Z)

      Character(80) ::  AF
      Character(6) ::  G92RWF 
      Character(5) ::  ELj,ELn
      Character(5), external ::  ELi

      Real(8), allocatable :: R(:),P(:),Q(:)

      Integer :: inp=1
      Integer :: out=2

      Call Read_name(AF)

      if(AF == '?') then
       write(*,*) 
       write(*,*) 'merging the set of GRASP w-files with choice of orbitals and' 
       write(*,*) 'with optional changing the spectroscopic notation'        
       write(*,*)                                                           
       write(*,*) '  1.w + 2.w + 3.w + ... --> result.w '             
       write(*,*)                                                           
       write(*,*) 'parameters are provided through interactive input'                  
       write(*,*) 
       Stop ' '
      end if

!------------------------------------------------------------------------------------

      write(*,*) 'Enter file-name for output w-file: '
      read(*,'(a)') AF
      Open(out,file=AF,FORM='UNFORMATTED')
      write(out) 'G92RWF'      

      i=1
   10 write(*,*) 'Enter file-name for input w-file  or  end: '
      read(*,'(a)') AF
      if(AF(1:3).eq.'end') go to 20
 
      i=Icheck_file(AF)
      if(i.eq.0) then
       write(*,'(a,a)') 'Can not find file  ',TRIM(AF) 
       go to 10
      end if

      Open(inp,file=AF,FORM='UNFORMATTED',ERR=10)
      READ(inp) G92RWF
      if(G92RWF.ne.'G92RWF') then
       write(*,'(a)') 'This is not a GRASP Radial File'
       go to 10
      end if

      Call Alloc_orb_jj(iwf)
    1 READ(inp,END=3) n,kappa,energy,npts
      if(allocated(R)) Deallocate(R,P,Q)
      Allocate(R(npts),P(npts),Q(npts))
      READ(inp) p0,P,Q
      READ(inp) R
      ELj = ELi(n,kappa,0)   
      i = Ifind_jjorb(n,kappa,0,0)
      if(i.ne.0) go to 1

      write(*,'(a)') ELj
    2 write(*,'(a)') 'add, rename or skip ? [ |s|EL]: '
      read(*,'(a5)') ELn 
      if(ELN(1:1).eq.'s') go to 1
      if(LEN_TRIM(ELn).eq.0) ELn=ELj

      Call EL_NLJK(ELn,nn,kk,ll,jj,kk)
      i = Ifind_jjorb(nn,kk,0,0)
      if(i.ne.0) then
       write(*,'(a)') 'we already have this name'
       go to 2 
      end if

      i = Ifind_jjorb(nn,kk,0,2)

      write(out) nn,kk,energy,npts
      write(out) p0,P,Q
      write(out) R

      go to 1

    3 Close(inp)
      go to 10

   20 write(*,'(a,i3,a)')  'Now there is ', nwf, '  w.f.:'
      write(*,'(15(1x,a3))') (ELF(i),i=1,nwf)

      Close(out)

      END  ! utility rw123

