!======================================================================
      Subroutine SUB_phys_pert
!======================================================================
!     define physical orbitals in perturber
!----------------------------------------------------------------------
      Use dbsr_prep

      Implicit none
      Integer :: i,k, ic, ip,jp, ii,jj      
      Real(8) :: S, SS, SM, s_ovl
      Real(8), allocatable :: S_orb(:), SS_orb(:)

      if(allocated(S_orb )) Deallocate(S_orb ); Allocate(S_orb (nbf))
      if(allocated(SS_orb)) Deallocate(SS_orb); Allocate(SS_orb(nbf))

      write(muc,'(/a)') 'Physical orbitals:'

      SS_orb = 0.d0
      Do jp=1,npert

       S_orb = 0.d0; SS = 0.d0
       Do ic = ippert(jp-1)+1,ippert(jp)
        S = WC(ic)*WC(ic); SS = SS + S
        Call Get_cfg_jj(ic)
        ip = ip_state(ic)
        Do i=1,no; ip=ip+1; ii=ipef(IP_orb(ip))   
         S_orb(ii) = S_orb(ii) + iq(i)*S 
        End do
       End do
       S_orb = S_orb / SS  
       Do i=1,nbf; SS_orb(i)=max(S_orb(i),SS_orb(i)); End do
        
      End do  !  over npert

      Do ii=1,nbf;  if(SS_orb(ii).lt.eps_phys) Cycle

! ... choose the substitution orbital with biggest overlap:
       
       SM = 0.d0; jj = 0 
       Do k=ncore+1,nbf
        if(ipbs(k).ne.1)      Cycle
        if(kbs(k).ne.kbs(ii)) Cycle
        s_ovl = abs(OBS(k,ii))
        if(S_ovl.lt.SM) Cycle
        SM = s_ovl; if(s_ovl.gt.eps_sub) jj = k
       End do
       
       if(jj.eq.0) then
        Call Add_sub_orbital(ii,jj) 
        SM = OBS(ii,jj)
       end if

       write(muc,'(a5,f8.1,5x,a5,f8.3)') ebs(ii),SS_orb(ii),ebs(jj),SM
       write(nuo,'(a5,f8.1,5x,a5,f8.3)') ebs(ii),SS_orb(ii),ebs(jj),SM

      End do  !  over orbitals        

      write(nuo,'(a)') '*'
      write(muc,'(a)') '*'

      End Subroutine SUB_phys_pert
