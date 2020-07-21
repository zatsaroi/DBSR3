!======================================================================
      Subroutine R_dbsrmat
!======================================================================
!     read the data from dbsr_mat.nnn file for the given partial wave
!----------------------------------------------------------------------
      Use dbsr_pol 
      Use channel_jj

      Implicit none
      Integer :: i,j,i1,i2,i3,j1,j2,ii,jj, ip,jp, ich,jch
      Real(8) :: S
      Integer, external :: Icheck_file

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
      if(i1.ne.ns )   Stop ' DBSR_POL: different ns  in DBSR_MAT file'
      if(i2.ne.nch)   Stop ' DBSR_POL: different kch in DBSR_MAT file'
      if(i3.ne.npert) Stop ' DBSR_POL: different kcp in DBSR_MAT file'

! ... allocate common working arrays:

      read(nui) nsol
      if(allocated(ipsol)) Deallocate(ipsol); Allocate(ipsol(0:nch))
      read(nui) ipsol
      Allocate(bval(nsol))
      read(nui) bval
      if(allocated(bb)) Deallocate(bb); Allocate(bb(ms,nsol))
      read(nui) bb

      nhm = nsol+npert
      mhm = nhm + nort + nortb
      khm = ms*nch + npert
      if(allocated(v)) Deallocate(v); Allocate(v(khm))

!----------------------------------------------------------------------
! ... read overlap matrix:

      if(allocated(om)) Deallocate(om); Allocate(om(nhm,nhm))
      om=zero; Do i=1,nsol; om(i,i)=one; End do

   1  read(nui) ich,jch
      if(ich.eq.0) go to 2
      if(ich.le.nch.and.jch.le.nch) then
       i1=ipsol(ich-1)+1; i2=ipsol(ich)
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       read(nui) om(i1:i2,j1:j2)
      elseif(ich.gt.nch.and.jch.le.nch) then
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       ip=ipsol(nch)+ich-nch
       read(nui) om(ip,j1:j2)
      else
       ip=ipsol(nch)+ich-nch
       jp=ipsol(nch)+jch-nch
       read(nui) om(ip,jp)
      end if
      go to 1
   2  Continue

      Do i=1,nhm; Do j=1,i; om(j,i)=om(i,j); End do; End do
      open(nus,form='UNFORMATTED',status='SCRATCH')
      Do i=1,nhm; write(nus) om(1:nhm,i); End do

!----------------------------------------------------------------------
! ... read Hamiltonian matrix:

      if(allocated(hm)) Deallocate(hm); Allocate(hm(mhm,mhm))
      hm=zero; Do i=1,nsol; hm(i,i)=bval(i); End do

   3  read(nui) ich,jch
      if(ich.eq.0) go to 4
      if(ich.le.nch.and.jch.le.nch) then
       i1=ipsol(ich-1)+1; i2=ipsol(ich)
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       read(nui) hm(i1:i2,j1:j2)
      elseif(ich.gt.nch.and.jch.le.nch) then
       j1=ipsol(jch-1)+1; j2=ipsol(jch)
       ip=ipsol(nch)+ich-nch
       read(nui) hm(ip,j1:j2)
      else
       ip=ipsol(nch)+ich-nch
       jp=ipsol(nch)+jch-nch
       read(nui) hm(ip,jp)
      end if
      go to 3
   4  Continue

      Do i=1,mhm; Do j=1,i; hm(j,i)=hm(i,j); End do; End do
      open(nua,form='UNFORMATTED',status='SCRATCH')
      Do i=1,mhm; write(nua) hm(1:mhm,i); End do

      Deallocate(bval)

      End Subroutine R_dbsrmat 

