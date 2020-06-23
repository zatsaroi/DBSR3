!======================================================================
!     UTILITY       D B S W _ T A B
!
!               C O P Y R I G H T -- 2008
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     connverts the B-spline radial orbital wbs-files into tab-files
!     suitable for grafic display
!----------------------------------------------------------------------
!
!     INPUT FILE:    name.bsw
!                    
!     OUTPUT FILES:  name.bsw.nl  for each nl
!
!----------------------------------------------------------------------
!     ARGUMENTS:     name.bsw  
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq
      Use zconst, only: c_au

      Implicit real(8) (A-H,O-Z)
      Character(1) :: ans
      Character(5) :: EL
      Character(40) :: AF,BF
      Real(8), allocatable ::  R(:),P(:),Q(:)

! ... input data: 
         
      iarg = IARGC()
      if(iarg.gt.0) Call GETARG(1,AF)

      if(iarg.lt.1.or.AF.eq.'?') then
       write(*,*)
       write(*,*) 'dbsw_tab connverts the B-spline radial orbital wbs-files into'
       write(*,*) 'text files suitable for grafic display'
       write(*,*)
       write(*,*) 'Call as:   dbsw_tab  name.bsw'
       write(*,*)
       write(*,*) 'file will be created for each orbital'
       Stop ' '
      end if

! ... set up B-splines:
 
      Call Check_file('knot.dat')
      Call def_grid
      Call alloc_DBS_gauss
      Call alloc_DBS_galerkin

! ... radial w.f.:

      nuw=1
      Open(nuw,file=AF,status='OLD',form='UNFORMATTED')
      Call Read_pqbs(nuw)
      Close(nuw)
      iaf = LEN_TRIM(AF)

! ... sets up grid points and initializes the values of the spline: 

      NR = nv*ks+2; Allocate(R(NR),P(NR),Q(NR))
      ii=1; R(1)=0.d0
      Do i=1,nv; Do j=1,ks; ii=ii+1; R(ii) = gr(i,j); End do; End do
      ii=ii+1; R(ii) = t(ns+1)

! ... Cycle over nl in input:

      BF=AF; iaf=iaf+1; BF(iaf:iaf)='.'

      Do i=1,nbf

       ip=i+i-1; jp=ip+1
       P=0.d0; Call Bvalue_bm(ksp,pbs(1,ip),P,pbsp)
       Q=0.d0; Call Bvalue_bm(ksq,pbs(1,jp),Q,qbsp)

       k=iaf; EL='     '; EL=ebs(i)
       Do j=1,5
        if(EL(j:j).eq.' ') Cycle
        k=k+1; BF(k:k)=EL(j:j)
       End do
       iout=2; Open(iout,file=BF)
       write(iout,'(3(5x,a,10x))') 'R','P','Q'

       S=1.d0
!       S = c_au * 2.d0 * t(ns+1) / kbs(i)

       Do j=1,nr
        write(iout,'(3D16.8)') R(j),P(j),Q(j)*S
       End do

      End do

      
      END   !  program dbsr_tab
          




!======================================================================
      Subroutine Read_pqbs(nu)
!======================================================================
!
!     read B-spline w.f. from bsw-file (unit nu) only those orbitals
!     which labels are already in the list
!     
!----------------------------------------------------------------------

      USE DBS_grid
      USE DBS_gauss
      USE DBS_orbitals_pq

      Implicit none

      Integer, intent(in) :: nu
      Integer :: i,j,k,l,n,m,itype,nsw,ksw,mw,kp,kq
      Character(5) :: elw
      Integer, External :: Ifind_bsorb,Iadd_bsorb 
      Real(8) :: tt(ns+ks)

      rewind(nu)
      read(nu) itype,nsw,ksw,tt,kp,kq
      if(itype.ne.grid_type) Stop ' Read_pqbs:  another grid_type'
      if(ksw.ne.ks) Stop ' Read_pqbs:  ksw <> ks'
      if(nsw.ne.ns) Stop ' Read_pqbs:  nsw <> ns'
      if(ksp.ne.kp) Stop ' Read_pqbs:  ksp <> kp'
      if(ksq.ne.kq) Stop ' Read_pqbs:  ksq <> kq'
      Do i=1,ns+ks
       if(tt(i).ne.t(i)) Stop ' Read_pqbs:  t <> tt'
      End do

    1 read(nu,end=2) elw,mw
      Call EL_NLJK(elw,n,k,l,j,i)
      m = Ifind_bsorb(n,k,i,2) 
      mbs(m)=mw 
      i=m+m-1; pbs(1:ns,i)=0.d0; read(nu) pbs(1:mw,i)
      i=m+m;   pbs(1:ns,i)=0.d0; read(nu) pbs(1:mw,i)
      go to 1
    2 Close(nu)

      End subroutine Read_pqbs
