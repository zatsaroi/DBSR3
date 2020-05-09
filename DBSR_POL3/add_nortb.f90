!======================================================================
      Subroutine add_nortb
!======================================================================
! ... provide oscillator strengths for given dipole matrix
!----------------------------------------------------------------------
      Use dbsr_pol
      Use channel_jj

      Implicit none
      Integer :: i,j, k,n, ip,is, ich, j1,j2, nstate,nhmb,nchb,npertb,nsb
      Real(8) :: S
      Real(8), allocatable :: a(:),b(:)
      Integer, external :: Ifind_position, Ipointer
       
      if(nortb.eq.0) Return
      open(nuq,form='UNFORMATTED',status='SCRATCH')

! ... check bound states file:

      i=INDEX(AF_bnd,'.'); AF=AF_bnd(1:i)//ALSP
      Call Check_file(AF); Open(nub,file=AF); rewind(nub)

      Call Read_ipar(nub,'nbound',nstate ) 
      Call Read_ipar(nub,'basis' ,nhmb   ) 
      Call Read_ipar(nub,'nchan' ,nchb   ) 
      Call Read_ipar(nub,'npert' ,npertb ) 
      Call Read_ipar(nub,'ns'    ,nsb    ) 

      if(nsb .ne.ns ) Stop 'dbsr_pol: nsb <> ns '
      if(nchb.ne.nch) Stop 'dbsr_pol: nchb - ?'
      if(npertb.ne.npert) Stop 'dbsr_pol: ncpb - ?'
      if(nhmb.ne.khm) Stop 'dbsr_pol: nhmb - ?'

      if(nstate.lt.nortb) Stop 'dbsr_pol: nstate < nortb - ?'

      Allocate(a(nhm),b(nhm))
      i=Ifind_position(nub,'nbound');  read(nub,*) 

      Do is=1,nstate
       read(nub,*) j
       read(nub,*) S
       read(nub,*) v
       ip = Ipointer(nortb,iortb,j)
       if(ip.eq.0) Cycle

       ! ... transform to a new basis
       
       Do ich=1,nch; i=(ich-1)*ms     
        j1 = ipsol(i-1)+1; j2=ipsol(i)
        Do j=j1,j2
         a(j) = SUM(v(i+1:i+ms)*bb(1:ms,j))
        End do
       End do
       if(npert.gt.0)  a(nsol+1:nhm)=v(khm-npert+1:khm)
       
       ! ... multiply on overlap matrix 

       Do j = 1,nhm
        b(j) = SUM(om(1:nhm,j)*a(1:nhm))
       End do 

       write(nuq) b
       iortb(ip) = -iortb(ip) 
       k=0; Do n=1,nortb; if(iortb(n).lt.0) Cycle; k=1; Exit; End do
       if(k.eq.0) Exit
      End do

      k=0; Do n=1,nortb; if(iortb(n).lt.0) Cycle; k=1; Exit; End do

      if(k.gt.0) Stop 'Add_nortb: not all states found'
      iortb = - iortb

      Deallocate(a,b)

      End Subroutine Add_nortb


!======================================================================
      Subroutine Check_nortb(sol)
!======================================================================
! ... provide oscillator strengths for given dipole matrix
!----------------------------------------------------------------------
      Use dbsr_pol

      Implicit none
      Integer ::  i,j, k,n, ip, nstateb
      Real(8) :: S, sol(*)
      Real(8), allocatable :: a(:),c(:)
      Integer, External :: Ifind_position, Ipointer
       
      if(nortb.eq.0) Return

       Allocate(a(nhm),c(nhm))
       Call Read_ipar(nub,'nbound',nstateb) 
       i=Ifind_position(nub,'nbound');  read(nub,*) 

       rewind(nuq)
       Do i=1,nstateb
        read(nub,*) j
        read(nub,*) S
        read(nub,*) C
        ip = Ipointer(nortb,iortb,i)
        if(ip.eq.0) Cycle
        read(nuq) a
        S = SUM(sol(1:nhm)*a(1:nhm)) 

        write(pri,'(a,i3,f12.8)') 'state ',i,S

        k=0; Do n=1,nortb; if(iortb(n).lt.0) Cycle; k=1; Exit; End do
        if(k.eq.0) Exit
       End do

       Deallocate(a,c)

       End Subroutine Check_nortb
