!======================================================================
!     UTILITY   R W _ B S W             
!
!               C O P Y R I G H T -- 2015
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     Generation double (p,q) B-spline representation for orbitals 
!     given in GRASP package format, w-files 
!----------------------------------------------------------------------
!
!     INPUT ARGUMENTS:
!     
!     name.w              -  input w-file (or file with list of w-files)
!     eps_end             -  tollerance for function tail (1.D-8)
!
!     INPUT FILES:
!
!     name.w              -  w-files  
!     knot.dat            -  parameters of B-splines
!
!     OUTPUT FILES:
!
!     rw_bsw.log          -  running information
!     name.bsw            -  B-spline representation for orbitals
!
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit real(8) (A-H,O-Z)

      Real(8) :: eps_end = 1.d-7,  eval(1000)

! ... files:

      Integer, parameter :: ma = 80
      Character (ma) :: AF, name
      Integer :: nuw =1; Character(ma) :: AF_w   = 'name.w'
      Integer :: nub =2; Character(ma) :: AF_bsw = 'name.bsw'
      Integer :: pri =6; Character(ma) :: AF_log = 'rw_bsw.log'

!----------------------------------------------------------------------
! ... check arguments:

      Call Read_name(AF_w)
      iarg = command_argument_count()
 
      if(iarg.lt.1.or.name.eq.'?') then
        write(*,*) 
        write(*,*) 'rw_bsw:   convert the GRASP w-file to DBSR bsw-file'
        write(*,*) 
        write(*,*) 'Call as:  rw_bsw  name.w  [eps_end=  debug=]'
        write(*,*) 
        write(*,*) 'optional parameters:'
        write(*,*) 
        write(*,*) 'eps_end  -  tail cut-off value [1.d-7] '
        write(*,*) 'debug - put > 0 for additional information in w_pqbs.log [0] '
        write(*,*) 
        write(*,*) 'OUTPUT:  name.bsw '
        Stop 
       end if
      
!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline: 

      Call Read_knot_dat           
      Call alloc_DBS_gauss
      CALL alloc_DBS_orbitals_pq(ibf,ns)

!----------------------------------------------------------------------

      Call Read_rarg('eps_end',eps_end)
      Call Read_iarg('debug',debug)

      Open(nuw,file=AF_w,form='UNFORMATTED',status='OLD')

      Open(pri,file=AF_log,position='APPEND')

      write(pri,'(a/)') 'Conversion to B-spline reprisantation:'

      ii=LEN_TRIM(AF_w); if(AF_w(ii-1:ii).eq.'.w') ii=ii-2
      AF_bsw = AF_w(1:ii)//'.bsw'
      Open(nub,file=AF_bsw,form='UNFORMATTED')

      write(pri,'(a,a)') 'Input file:  ',AF_w
      write(pri,'(a,a)') 'Output file: ',AF_bsw

      write(pri,'(/a,i5)') 'ksp =',ksp
      write(pri,'( a,i5)') 'ksq =',ksq

      write(pri,'(/a,T20,F20.10)')  'eps_end tolerence: ',eps_end
      write(pri,'( a,T20,F20.10/)')  'Border radius: ',t(ns+1)

!----------------------------------------------------------------------
! ... process the radial functions:

      Call Read_Grasp (nuw,pri,eps_end,eval)

!---------------------------------------------------------------------
! ... record the bsw-functions:

      write(nub) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i=1,nbf
       write(nub) ebs(i),mbs(i), eval(i)
       write(nub) pq(1:mbs(i),1,i)
       write(nub) pq(1:mbs(i),2,i)
      End do

      END    !  PROGRAM rw_pqbs


!======================================================================
      Subroutine Read_Grasp(nu,log,end_tol,eval)
!======================================================================
!     Generation double (p,q) B-spline representation for orbitals 
!     given in GRASP package format, w-files 
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq,  p => pq

      Implicit none

      Real(8), allocatable :: R(:),PG(:),QG(:)     
      Real(8), allocatable :: b3(:),c3(:),d3(:),b4(:),c4(:),d4(:)
      Real(8) :: pr(nv,ks),qr(nv,ks),a(ns,ns),Pcoef(ns),Qcoef(ns), eval(*)
      Character(6) string
      Integer :: nu,log,n,kappa,npts,is,i,j,ii,np,nq,ip,iv
      Real(8) :: energy,p0,pn,PP,QQ,d,pm,qm,end_tol
      Integer, external :: Ifind_bsorb
      Real(8), external :: SEVAL, bvalu2, QUADR

! ... check the file:

      read(nu) string
      if(string.ne.'G92RWF') Stop 'This is not a GRASP Radial File'

      write(log,'(/a/)') 'GRASP orbitals:'

! ... process the radial functions:

    1 read(nu,end=2)  n,kappa,energy,npts
write(*,*)  n,kappa,energy,npts

      if(allocated(R)) Deallocate(R,PG,QG)
      Allocate(R(npts),PG(npts),QG(npts))
      READ(nu) p0,PG,QG
      READ(nu) R 

! ... check if we need this orbital:

      is = n/1000; n = n - is*1000
      ii = Ifind_bsorb(n,kappa,is,2)
      if(ii.eq.0) go to 1
      eval(ii) = energy

! ... remove unphysical oscilations in the end:

      PM = maxval(abs(PG))
      Do i=npts,1,-1
       if(abs(PG(i))/PM.gt.end_tol) Exit; PG(i)=0.d0
      End do
      np = i

      QM = maxval(abs(QG))
      Do i=npts,1,-1
       if(abs(QG(i))/QM.gt.end_tol) Exit; QG(i)=0.d0
      End do
      nq = i

! ... interpolate function p(r) to cubic splines:

      if(allocated(b3)) Deallocate(b3,c3,d3,b4,c4,d4)
      Allocate(b3(np),c3(np),d3(np),b4(nq),c4(nq),d4(nq))
      Call Splin3(np,R,PG,b3,c3,d3)
      Call Splin3(nq,R,QG,b4,c4,d4)

! ... evaluate the function in gaussian points:

      pr = 0.d0; qr = 0.d0
      Do i=1,nv; Do j=1,ks
       if(gr(i,j).lt.R(np)) pr(i,j)=SEVAL(np,gr(i,j),R,PG,b3,c3,d3)*grw(i,j)  
       if(gr(i,j).lt.R(nq)) qr(i,j)=SEVAL(nq,gr(i,j),R,QG,b4,c4,d4)*grw(i,j)  
      End do; End do

! ... form the vector of inner products of the radial function and the
! ... spline basis functions

      Pcoef = 0.d0
      Do iv = 1,nv; Do ip = 1,ksp; i = iv+ip-1
       Pcoef(i) = Pcoef(i) + SUM(pr(iv,:)*pbsp(iv,:,ip))
      End do; End do

      Qcoef = 0.d0
      Do iv = 1,nv; Do ip = 1,ksq; i = iv+ip-1
       Qcoef(i) = Qcoef(i) + SUM(qr(iv,:)*qbsp(iv,:,ip))
      End do; End do

! ... solve the system of equations for coef
   
      a(1:nsp-1,1:nsp-1)=fpbs(2:nsp,2:nsp)
      Call gaussj (a,nsp-1,ns,Pcoef(2),1,ns)
      Pcoef(1)=0.d0

      a(1:nsq-1,1:nsq-1)=fqbs(2:nsq,2:nsq)
      Call gaussj (a,nsq-1,ns,Qcoef(2),1,ns)
      Qcoef(1)=0.d0

      p(:,1,ii) = Pcoef
      p(:,2,ii) = Qcoef
      Call Check_tails(ii,end_tol)
     
! ... check the normalization 

      PN = sqrt(QUADR(p(1,1,ii),p(1,1,ii),0))
      p(:,:,ii)=p(:,:,ii)/PN

! ... max.deviation from original orbital:

      pm = 0.d0
      Do j=2,np
       PP = bvalu2(tp,p(1,1,ii),nsp,ksp,R(j),0) 
       d = abs(PP-PG(j));  if(d.gt.pm) pm = d
      End do

      qm = 0.d0
      Do j=2,nq
       QQ = bvalu2(tq,p(1,2,ii),nsq,ksq,R(j),0)
       d = abs(QQ-QG(j)); if(d.gt.qm) qm = d
      End do

      write(log,'(a5,4(a,E12.3))')  &
           ebs(ii), '  diff_p =',pm,'   diff_q =',qm,'   (norm -1) =',PN-1.d0

      go to 1
    2 Close(nu)

      End Subroutine Read_GRASP



!======================================================================
      Subroutine Check_tails(jo,end_tol)
!======================================================================
!     nulify too small coefficients in the end for orbital jo
!     if jo=0 - check all orbitals
!     only large component is used to define small B-splines
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_orbitals_pq
      
      Implicit none
      Integer :: io,jo,i
      Real(8) :: cm, end_tol

      Do io=1,nbf; if(jo.ne.0.and.io.ne.jo) Cycle
       cm = maxval( abs( pq(:,1,io) ) )
       Do i=nsp-1,1,-1
        mbs(io)=i     
        if(abs(pq(i,1,io))/cm.lt.end_tol) Cycle
        Exit
       End do
       pq(i+1:ns,1,io)=0.d0
       pq(i+1:ns,2,io)=0.d0
      End do

      End Subroutine Check_tails

