!======================================================================
      Subroutine SUB_phys_orb
!======================================================================
!     define physical orbitals for given target state
!----------------------------------------------------------------------
      Use dbsr_prep

      Implicit none
      Integer :: i,k, ii,jj, ic,jc, ip 
      Real(8) :: S, SS, SM, s_ovl
      Real(8), allocatable :: S_orb(:)

      if(allocated(S_orb)) Deallocate(S_orb)
      Allocate(S_orb(nbf))

! ... find the occupation numbers for all orbitals:

      S_orb = 0.d0; SS = 0.d0
      Do jc = 1,ncfg; ic = ipt(jc); S = WC(ic)*WC(ic)
       SS = SS + S
       Call Get_cfg_jj(ic)
       ip = ip_state(ic)
       if(jc.eq.1) S = 1.d0
       Do i=1,no; ip=ip+1; ii=ipef(ip_orb(ip))
        S_orb(ii) = S_orb(ii) + iq(i)*S
       End do
       if(SS.gt.c_conf) Exit
      End do
      S_orb = S_orb / SS

      write(muc,'(/a)') 'Physical orbitals:'
      SS = 0
      Do ii=1,nbf;  if(S_orb(ii).lt.eps_phys) Cycle

! ... choose the substitution orbital with biggest overlap:

       SM = 0.d0; jj = 0
       Do k=ncore+1,nbf
        if(ipbs(k).ne.1)      Cycle
        if(kbs(k).ne.kbs(ii)) Cycle
        if(nbs(ii).lt.9.and.nbs(k).lt.nbs(ii)) Cycle       ! ???
        s_ovl = abs(OBS(k,ii))
        if(s_ovl.lt.SM)       Cycle
        SM = s_ovl; if(s_ovl.gt.eps_sub) jj = k
       End do
       if(jj.eq.0) then
        Call Add_sub_orbital(ii,jj)
        SM = OBS(ii,jj)
       end if

       write(muc,'(a5,f8.1,5x,a5,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM
       write(nuo,'(a5,f8.1,5x,a5,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM

       SS = SS + S_orb(ii)

      End do

      write(muc,'(a,20x,f8.3)') '*'     !,nelc-SS
      write(nuo,'(a,20x,f8.3)') '*'     !,nelc-SS

      End Subroutine SUB_phys_orb
