!======================================================================
      Subroutine Record_matrix 
!======================================================================
!     record the overlap or Hamiltonian matrix to unit nui
!----------------------------------------------------------------------
      Use dbsr_mat      

      Implicit none
      Integer :: i,j,k, ni,nj, ip,jp, ich,jch
      Real(8) :: S, w(ms), c(ns,ns)

! ... non-diagonal channel blocks:

      Do ich = 2,nch;    ni=ipsol(ich)
      Do jch = 1,ich-1;  nj=ipsol(jch)
       k=icc(ich,jch); if(k.eq.0) Cycle 
       if(maxval(abs(hch(:,:,k))).lt.1.d-20) Cycle
       c = 0.d0
       Do i=1,ni
        Do j=1,ms; w(j)  = SUM(diag(:,i,ich)*hch(:,j,k)); End do 
        Do j=1,nj; c(i,j)= SUM(w(:)*diag(:,j,jch));       End do 
       End do       
       write(nui) ich,jch
       write(nui) c(1:ni,1:nj)
       if(icase.eq.0) overlaps(ich,jch) = maxval(abs(c(1:ni,1:nj)))
      End do; End do

! ... pertubers:

      if(npert.gt.0) then

! ... channel-perturber rows:

      Do ich = 1,nch;    ni = ipsol(ich)
      Do ip  = 1,npert
       k = icb(ich,ip); if(k.eq.0) Cycle
       S = maxval(abs(hcp(:,k)))
       if(S.lt.1.d-20) Cycle
       w = 0.d0
       Do i=1,ni; w(i)=SUM(diag(:,i,ich)*hcp(:,k)); End do
       if(icase.eq.0) overlaps(nch+ip,ich) = maxval(abs(w(1:ni)))
       write(nui) nch+ip,ich
       write(nui) w(1:ni)
      End do; End do
       
! ... perturter-perturber elements:

      Do ip = 1,npert; Do jp = 1,ip; k = ibb(ip,jp)
       if(k.eq.0) Cycle
       if(hp(k).eq.0.d0) Cycle
       write(nui) ip+nch,jp+nch
       write(nui) hp(k)
       if(icase.eq.0) overlaps(nch+ip,nch+jp) = abs(hp(k))
      End do; End do

      end if   ! over npert

      write(nui) 0,0  ! sign of the end

      End Subroutine Record_matrix

