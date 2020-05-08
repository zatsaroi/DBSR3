!======================================================================
      Subroutine Assign_index(m)
!======================================================================
!     assign set index for orbital m
!----------------------------------------------------------------------
      Use dbsr_prep

      Implicit none
      Character(5) :: elw
      Character(5), external :: ELi
      Integer :: i,k,m,knew
      Integer, external :: New_index
      Real(8) :: S

      ibs(m)=-1; elw = EBS(m)
      knew = New_index(kbs(m),ksmax,nbf,kbs,ibs)

      Do k = 1,knew-1
       S=0.d0
       Do i = 1,nbf
        if(kbs(m).ne.kbs(i).or.k.ne.ibs(i)) Cycle
        S=max(S,abs(OBS(i,m)))
       End do
       if(S.lt.eps_ovl) then; ibs(m)=k; Exit; end if
      End do

      if(ibs(m).eq.-1) then  ! the orbital belongs to new set

       ibs(m) = knew
       EBS(m)=ELi(nbs(m),kbs(m),ibs(m))
       write(pri,'(a,a,a,T24,a,f15.8)') &
              elw,' --> ',EBS(m),'new orbital and new set index',S

      else

       EBS(m)=ELi(nbs(m),kbs(m),ibs(m))
       Do i = 1,nbf          ! check the same label for diff.orbitals
        if(EBS(m).ne.EBS(i)) Cycle
        ibs(m) = knew
        EBS(m)=ELi(nbs(m),kbs(m),ibs(m))
        write(pri,'(a,a,a,T24,a,f15.8)') &
              elw,' --> ',EBS(m),'new orbital and new set index',S
       End do

       if(ibs(m).ne.knew) write(pri,'(a,a,a,T24,a,f15.8)') &
              elw,' --> ',EBS(m),'new orbital but old set index',S

      end if

      End Subroutine Assign_index

