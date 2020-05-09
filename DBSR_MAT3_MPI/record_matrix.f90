!======================================================================
      Subroutine Record_matrix 
!======================================================================
!     record the overlap or Hamiltonian matrix to unit nui
!----------------------------------------------------------------------
      Use dbsr_mat      

      Implicit none
      Integer :: i, ni,nj, ip,jp, ich,jch
      Real(8) :: S

! ... non-diagonal channel blocks:

      Do ich = 2,nch;    ni=ipsol(ich)
      Do jch = 1,ich-1;  nj=ipsol(jch)
       i=icc(ich,jch); if(i.eq.0) Cycle 
       S = maxval(abs(hch(1:ni,1:nj,i)))
       if(S.lt.1.d-20) Cycle
       write(nui) ich,jch
       write(nui) hch(1:ni,1:nj,i)
      End do; End do

! ... pertubers:

      if(npert.gt.0) then

! ... channel-perturber rows:

      Do ich = 1,nch;    ni = ipsol(ich)
      Do ip  = 1,npert;  i = icb(ich,ip)
       if(i.eq.0) Cycle
       S = maxval(abs(hcp(1:ni,i)))
       if(S.lt.1.d-20) Cycle
       write(nui) ip+nch,ich
       write(nui) hcp(1:ni,i)
      End do; End do
       
! ... perturter-perturber elements:

      Do ip = 1,npert; Do jp = 1,ip; i = ibb(ip,jp)
       if(i.eq.0) Cycle
       if(hp(i).eq.0.d0) Cycle
       write(nui) ip+nch,jp+nch
       write(nui) hp(i)
      End do; End do

      end if   ! over npert

      write(nui) 0,0  ! sign of the end

      End Subroutine Record_matrix
