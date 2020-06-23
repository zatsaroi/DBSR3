!----------------------------------------------------------------------
!     Utility rw_dat
!----------------------------------------------------------------------
!
!     prepares the text files suitable for plotting the rad. functions:
!
!     name.w  -->  name.nl1, name.nl1, name.nl1, ...
!  
!                  where nl's denote one-electron orbitals
!
!----------------------------------------------------------------------
      IMPLICIT REAL(8) (A-H, O-Z)

      Character(80) :: AF, BF
      Character(6) :: G92RWF
      Character(5) :: EL, ELi

      Real(8), allocatable :: R(:),P(:),Q(:)

      iarg = command_argument_count() 
      if(iarg.gt.0) then
       Call GET_COMMAND_ARGUMENT(1,AF)
      else
       Stop 'enter name for w-file as argument'
      end if

      if(AF.eq.'?') then
       write(*,*) 
       write(*,*) 'rw_dat prepares the text files (one for each orbital)' 
       write(*,*) 'suitable for plotting the radial functions:'
       write(*,*)
       write(*,*) 'name.w  -->  name.nl1, name.nl1, name.nl1, ... '     
       write(*,*)                                                     
       write(*,*) "             where nl's denote one-electron orbitals"
       write(*,*) 
       write(*,*) 'Call as:   rw_dat  name.w '
       write(*,*)
       Stop 
      end if


      ii = INDEX(AF,'.',BACK=.TRUE.); if(ii.eq.0) ii=LEN_TRIM(AF) 

      inp =1;  Open(inp, file=AF,form='UNFORMATTED',status='OLD')

      READ(inp) G92RWF
      if(G92RWF.ne.'G92RWF') Stop 'This is not a GRASP Radial File'

    1 READ(inp,END=2) n,kapa,energy,npts
      Allocate(R(npts),P(npts),Q(npts))
      READ(inp) p0,P,Q
      READ(inp) R

      EL = ELi(n,kapa,0)
      Do i=1,5; if(EL(i:i).eq.' ') Cycle; jj=i; Exit; End do

      BF = AF(1:ii)//EL(jj:5);  iout=2;  Open(iout,file=BF)

!      write(iout,'(/a,a5,a,i3,a,i3,a,i5,a,d23.16,a,d23.16/)') '#  ',EL,&
!       '  n=',n,'  k=',kapa,'  npts=',npts,'  E=',energy,'  p0=',p0

      write(iout,'(3(12x,a1,12x))')       'R','P','Q'

      Do i=1,npts
       write(iout,'(3d25.16)') R(i),P(i),Q(i)
      End do       

      Deallocate(R,P,Q)

      go to 1
    2 Continue

      END  ! utility rw_dat
