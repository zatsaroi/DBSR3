!======================================================================
      Subroutine R_dvector
!======================================================================
!     read the data from dv.nnn file for the given partial wave
!----------------------------------------------------------------------
      Use dbsr_pol
      Use channel_jj, only: nch,npert

      Implicit none
      Integer :: i,j, kdm,kch,kpert, ich, j1,j2

! ... open file with the dipol vector:

      i = Len_trim(AF_dip)-3; AF=AF_dip(1:i)//ALSP
      Call Check_file(AF)
      Open(nud,file=AF,form='UNFORMATTED')

! ... check the dimensions:

      read(nud) kdm,kch,kpert,E1,jot1,kpol,ktype

      if(kdm.ne.khm)      Stop ' DBSR_POL: different mhm in dv-file'
      if(kch.ne.nch)      Stop ' DBSR_POL: different nch in dv-file'
      if(kpert.ne.npert)  Stop ' DBSR_POL: different npert in dv-file'

      if(Allocated(d)) Deallocate(d);  Allocate(d(mhm)); d=zero
      read(nud) v(1:khm)

! ... convert to new basis:

      Do ich=1,nch; i=(ich-1)*ms     
       j1 = ipsol(ich-1)+1; j2=ipsol(ich)
       Do j=j1,j2
        d(j) = SUM(v(i+1:i+ms)*bb(1:ms,j))
       End do
      End do
      if(npert.gt.0)  d(nsol+1:nhm)=v(khm-npert+1:khm)

      write(pri,*)
      write(pri,'(a)') 'D-matrix:'
      write(pri,*)
      write(pri,'(a,i5,a)')  'kpol   = ',kpol,'  - multipole index'
      write(pri,'(a,i5,a)')  'kdm    = ',kdm ,'  - dipole vector length'
      write(pri,'(a,i5,a)')  '2J     = ',jot1,'  - 2J-value of the initial state'
      write(pri,'(/a,f16.8,a)') 'E1     = ',E1,  '  - energy of the initial state'

      End Subroutine R_dvector 

