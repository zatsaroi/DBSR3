!=======================================================================
!     merging the set of bsw-files with choice of orbitals and 
!     with optional changing the spectroscopic notation
! 
!       1.bsw + 2.bsw + 3.bsw + ... --> result.bsw 
!
!     Interactive input; knot.dat is required
!=======================================================================
      Use DBS_grid
      Use DBS_orbitals_pq

      Implicit real(8) (A-H, O-Z)

      Character(80) ::  AF
      Character(5) ::  ELw,ELn

      Real(8), allocatable :: R(:),P(:),Q(:)

      Integer :: inp=1
      Integer :: out=2

      Real(8), allocatable :: tt(:)

      Call Read_name(AF)

      if(AF == '?') then
       write(*,*) 
       write(*,*) 'merging the set of bsw-files with choice of orbitals and' 
       write(*,*) 'with optional changing the spectroscopic notation'        
       write(*,*)                                                           
       write(*,*) '  1.bsw + 2.bsw + 3.bsw + ... --> result.bsw '             
       write(*,*)                                                           
       write(*,*) 'Interactive input; knot.dat is required'                  
       write(*,*) 
       Stop ' '
      end if

! ... prepare B-spline parameters:

      Call read_knot_dat
      Call alloc_DBS_orbitals_pq(ibf,ns)
      Allocate(tt(ns+ks))

   10 write(*,*) 'Enter file-name for input w-file  or  end: '
      read(*,'(a)') AF
      if(AF(1:3).eq.'end') go to 20
 
      i=Icheck_file(AF)
      if(i.eq.0) then
       write(*,'(a,a)') 'Can not find file  ',TRIM(AF) 
       go to 10
      end if

      Open(inp,file=AF,FORM='UNFORMATTED',ERR=10)

      read(inp) itype,nsw,ksw,tt,ksp,ksq
      if(ksw.ne.ks) Stop ' ksw <> ks'
      if(nsw.ne.ns) Stop ' nsw <> ns'
      if(itype.ne.grid_type) Stop ' another grid_type ?'    
      
    1 read(inp,end=3) elw,mw

    2 write(*,'(a)') elw
      write(*,'(a)') 'add, rename or skip ? [ |s|EL]: '
      read(*,'(a5)') ELn 
      if(ELN(1:1).eq.'s') then; read(inp) x; read(inp) x; goto 1; endif
      if(LEN_TRIM(ELn).eq.0) ELn=elw

      Call EL_NLJK(ELn,n,k,l,j,i)
      m = Ifind_bsorb(n,k,i,0)
      if(m.ne.0) then
       write(*,'(a)') 'we already have this name'
       go to 2 
      end if

      m = Ifind_bsorb(n,k,i,2)
!      write(*,*) n,k,i,m,ebs(m)

      pq(:,:,m) = 0.d0
      mbs(m)=mw
      read(inp) pq(1:mw,1,m)
      read(inp) pq(1:mw,2,m)

      go to 1

    3 Close(inp)
      go to 10

   20 write(*,'(a,i3,a)')  'Now there is ', nbf, '  w.f.:'
      write(*,'(15(1x,a5))') (ebs(i),i=1,nbf)

      write(*,*) 'Enter file-name for output w-file: '
      read(*,'(a)') AF
      Open(out,file=AF,FORM='UNFORMATTED')

      write(out) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i=1,nbf
       write(out) ebs(i),mbs(i)
       write(out) pq(1:mbs(i),1,i)
       write(out) pq(1:mbs(i),2,i)
      End do

      Close(out)

      End  ! utility bsw123

