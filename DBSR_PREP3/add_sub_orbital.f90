!======================================================================
      Subroutine Add_sub_orbital(ii,jj)
!======================================================================
!     add new substitution orbital if needed
!----------------------------------------------------------------------
      Use dbsr_prep

      Implicit none
      Integer :: i,m, ii,jj
      Real(8) :: S

      m=nbf+1; if(m.ge.mbf) Call Alloc_DBS_orbitals_pq (mbf+ibf,ns)
      nbs(m)=nbs(ii); kbs(m)=kbs(ii); ibs(m)=ibs(ii); ebs(m)=ebs(ii)
      mbs(m)=mbs(ii); pq(:,:,m) = pq(:,:,ii)
      Call Idel_obs(m)
      Do i = 1,m-1
       if(kbs(i).ne.kbs(m)) Cycle; S=OBS(i,ii); if(S.eq.0.d0) Cycle
       Call Iadd_obs(i,m,S)
      End do

      S = OBS(ii,ii); Call Iadd_obs(m,m,S)

      OBS1(m) = OBS1(ii)
      OBS2(m) = OBS2(ii)

! ... check if we need orthogonalization:

      jj = ii
      Do i=ncore+1,nbf; if(ipbs(i).ne.1) Cycle
       if(kbs(i).ne.kbs(m)) Cycle
       S = OBS(i,m); if(abs(S).lt.eps_ovl) Cycle
       pq(:,:,m) = pq(:,:,m) - S * pq(:,:,i);  jj=m
       mbs(m) = max(mbs(i),mbs(m))
      End do
      if(jj.eq.ii) then; ipbs(ii)=1; Return; End if

      bpq(1:ns,1,m) = MATMUL(fpbs, pq(:,1,m))
      bpq(1:ns,2,m) = MATMUL(fqbs, pq(:,2,m))

      S = QUADR_00(m,m); pq(:,:,m)=pq(:,:,m)/sqrt(S) 

      bpq(1:ns,1,m) = MATMUL(fpbs, pq(:,1,m))
      bpq(1:ns,2,m) = MATMUL(fqbs, pq(:,2,m))

      Do i = 1,m
       if(kbs(i).ne.kbs(m)) Cycle
       S = QUADR_00(i,m)
       Call Iadd_obs(m,i,S)
      End do
      OBS1(m) = QUADR_pq(m,m,1)
      OBS2(m) = QUADR_pq(m,m,2)

! ... assign set index for new sub. orbital:

      Call Assign_index(m); nbf=m; ipbs(m)=1; jj=m

      End Subroutine Add_sub_orbital

