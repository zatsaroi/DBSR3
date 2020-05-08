!======================================================================
      Subroutine SUB_check_orb
!======================================================================
!     this routine analizes one set of orbitals in file AFW
!----------------------------------------------------------------------
      Use dbsr_prep

      Implicit none

      Character(5) :: elw
      Character(5), external :: ELi

      Integer :: i,j,k,l,m,n, ii, kappa, nsw,ksw,mbw, kp,kq, igrid 
      Integer, external :: Ifind_jjorb

      Real(8) :: S, S1, S2, pbw(ns),qbw(ns),tw(ns+ks)

!----------------------------------------------------------------------
! ... read radial functions:

      Call Check_file(AFW);  Open(nuw,file=AFW,form='UNFORMATTED')
      rewind(nuw)
      read(nuw) igrid,nsw,ksw,tw(1:nsw+ksw),kp,kq
      ! ... check consistence of B-spline basis:
      if(grid_type.gt.0.and.igrid.ne.grid_type) Stop 'Another knot grid?'
      if(ksw.ne.ks) Stop ' Read_bsw:  ksw <> ks'
      if(nsw.ne.ns) Stop ' Read_bsw:  nsw <> ns'
      if(ksp.ne.kp) Stop ' Read_bsw:  ksp <> kp'
      if(ksq.ne.kq) Stop ' Read_bsw:  ksq <> kq'
      k=1
      Do i=1,ns+ks
       if(abs(t(i)-tw(i)).lt.1.d-10) Cycle; k=0; Exit
      End do
      if(k.eq.0) Stop 'Another knot grid t(1:ns+ks) ?'

    1 m=nbf+1; if(m.gt.mbf) CALL Alloc_DBS_orbitals_pq (mbf+ibf,ns)
      Call Idel_obs(m)
      read(nuw,end=2) elw,mbw
      pbw=0.d0; read(nuw) (pbw(i),i=1,mbw)
      qbw=0.d0; read(nuw) (qbw(i),i=1,mbw)
      Call EL_nljk(elw,n,kappa,l,j,k);   k = k + kshift

      ii=Ifind_jjorb(n,kappa,k,0)

      if(ii.eq.0.and.kshift.gt.0) then
        k=0; ii=Ifind_jjorb(n,kappa,0,0) 
      end if

      if(ii.gt.0) then
       mbs(m)=ns; nbs(m)=n; kbs(m)=kappa; ibs(m)=k; ebs(m)=elw
       pq(1:ns,1,m)=pbw(1:ns)
       pq(1:ns,2,m)=qbw(1:ns)
      else
       write(pri,'(a5,T24,a)') elw,'excessive orbital'
       go to 1
      end if

! ... define overlaps with existing orbitals:

      Do i = 1,m
       if(kbs(i).ne.kbs(m)) Cycle 
       S = QUADR_pq(m,i,0)
       if(abs(S).lt.eps_ovl) Cycle
       if(i.le.ncore.and.abs(S).lt.eps_core) Cycle
       Call Iadd_obs(m,i,S)
      End do

      OBS1(m) = QUADR_pq(m,m,1)
      OBS2(m) = QUADR_pq(m,m,2)

!---------------------------------------------------------------------
! ... compare with the existing orbitals:

      Do i = 1,nbf;  if(abs(OBS(i,m)).lt.eps_ovl) Cycle

! ... define is that approximately the same orbital:

       S  = abs(OBS(i,m)-OBS(m,m)) +  abs(OBS(i,m)-OBS(i,i))
       S1 = abs(OBS1(i)-OBS1(m))
       S2 = abs(OBS2(i)-OBS2(m))
 
       if(S.lt.eps_ovl.and.S1.lt.eps_ovl.and.S2.lt.eps_ovl) then
        ipef(ii) = i
        if(ii.gt.ncore) &
          write(pri,'(a,a,a,T24,a)') elw,' --> ',ebs(i),'existing orbital'
        go to 1
       end if

! ... core orbitals should be the same:
 
      if(ii.le.ncore) then
       write(pri,'(3E12.1)') S,S1,S2
       write(pri,'(a,a,a)') 'file ',trim(AFW),'  has another core orbital'
       Stop 'another core orbital ? '
      end if

! ... check orthogonality to core:

       if(i.le.ncore) then
        if(abs(OBS(i,m)).gt.eps_core) then
        write(pri,'(a,a,a,E12.3)') &
        '  orbital ', elw, ' does not orthogonal to core:', OBS(i,m)
        Stop 'problem with orthogonality to core'
        end if
       end if

      End do

! ... assign set index for new orbital:

       Call Assign_index(m); nbf=m; ipef(ii)=m

       go to 1    ! go to next orbital
   2  Continue

! ... check if we bsw-file contains all radial functions:

      Do i=ncore+1,nwf; ii=ipef(i)
       if(ii.eq.0) then
        write(pri,'(a,a,a)') 'orbital ',ELF(i),' not found in the bsw-file'
        Stop ' unknown orbitals ! '
       else
        NEF(i)=nbs(ii);  KEF(i)=kbs(ii); IEF(i)=ibs(ii); ELF(i)=ebs(ii)
       end if
      End do

      End Subroutine Sub_check_orb
