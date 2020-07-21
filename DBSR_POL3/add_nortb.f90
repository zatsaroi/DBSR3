!======================================================================
      Subroutine add_nortb
!======================================================================
! ... provide oscillator strengths for given dipole matrix
!----------------------------------------------------------------------
      Use dbsr_pol
      Use channel_jj

      Implicit none
      Integer :: i,j, k,n, ip,is, ich, j1,j2, nstate,nhmb,nchb,npertb,nsb, &
                 jotb, parityb
      Real(8) :: S
      Real(8), allocatable :: a(:),b(:)
      Integer, external ::  Ipointer
       
      if(nortb.eq.0) Return
      open(nuq,form='UNFORMATTED',status='SCRATCH')

! ... check bound states file:

      i=INDEX(AF_bnd,'.'); AF=AF_bnd(1:i)//ALSP
      Call Check_file(AF); Open(nub,file=AF,form='UNFORMATTED')

      rewind(nub)
      read(nub) nhmb,nchb,npertb,nsb,jotb,parityb,nstate

      if(nsb   .ne.  ns ) Stop 'dbsr_pol: nsb <> ns '
      if(nchb  .ne.  nch) Stop 'dbsr_pol: nchb - ?'
      if(npertb.ne.npert) Stop 'dbsr_pol: ncpb - ?'
      if(nhmb  .ne.  khm) Stop 'dbsr_pol: nhmb - ?'

      if(nstate.lt.nortb) Stop 'dbsr_pol: nstate < nortb - ?'

      Allocate(a(nhm),b(nhm))

      Do is=1,nstate

       read(nub) j
       read(nub) S
       read(nub) v,a
       ip = Ipointer(nortb,iortb,j)
       if(ip.eq.0) Cycle

       
write(pri,'(/a/)') 'Solution'
write(pri,'(10E15.5)') a(1:nhm)
write(pri,*) 'nhm =',nhm

       ! ... multiply on overlap matrix 

       Do j = 1,nhm
        b(j) = SUM(om(1:nhm,j)*a(1:nhm))
       End do 

      S = SUM(a(1:nhm)*b(1:nhm))
write(*,*) 'O =', S
      S = 0.0
      Do i = 1,nhm; Do j = 1,nhm
       S = S + a(i)*hm(i,j)*a(j)
      End do; End do
write(*,*) 'E =', S



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
      Integer ::  i,j, k,n, ip
      Integer ::  nhmb,nchb,npertb,nsb,jotb,parityb,nstateb
      Real(8) :: S, sol(*)
      Real(8), allocatable :: a(:),c(:)

      Integer, external :: Ipointer
       
      if(nortb.eq.0) Return

       Allocate(a(nhm),c(nhm))

       rewind(nub)
       read(nub) nhmb,nchb,npertb,nsb,jotb,parityb,nstateb

       rewind(nuq)
       Do i=1,nstateb
        read(nub) j
        read(nub) S
        read(nub) C

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
