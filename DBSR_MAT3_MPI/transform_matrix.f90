!======================================================================
      Subroutine transform_matrix
!======================================================================
!     transform overlap/Hamiltonian  matrix to a new basis,
!     constructed from diagonalization of diagonal one-channel blocks.
!----------------------------------------------------------------------
!     for non-diagonal channel blocks, H_ij, if any, we shoud transform 
!
!                c_i(1:ms)^T  H_ij(1:ms,1:ms) c_j(1:ms) 
!
!     for each channels solutions c_i, c_j, saved in array diag(1:ms,..);
!     ipsol array provides number of solutions for each channel.
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Real(8) :: w(ms), c(ms,ms)
      Integer :: i,j,k, ich,jch, isol,jsol

      if(icase.eq.0) then
       if(allocated(overlaps)) Deallocate(overlaps)
       Allocate(overlaps(nch+npert,nch+npert)); overlaps = 0.d0
      end if
 
! ... channel-channel blocks:

      Do ich=1,nch; isol=ipsol(ich)            
      Do jch=1,ich; jsol=ipsol(jch)
       k = icc(ich,jch); if(k.eq.0) Cycle
       if(ich.eq.jch) Cycle
       if(maxval(abs(hch(:,:,k))).lt.1.d-20) Cycle

       c = 0.d0
       Do i=1,isol
        Do j=1,ms; w(j)=SUM(diag(:,i,ich)*hch(:,j,k)); End do 
        Do j=1,jsol; c(i,j)=SUM(w(:)*diag(:,j,jch)); End do 
       End do       
       hch(:,:,k) = c(:,:)
       if(icase.eq.0) overlaps(ich,jch) = maxval(abs(c(1:isol,1:jsol)))

      End do; End do

! ... channel-bound blocks:

      if(npert.gt.0) then

      Do ich=1,nch; isol=ipsol(ich)            
      Do jch = 1,npert 
       k = icb(ich,jch); if(k.eq.0) Cycle
       if(maxval(abs(hcp(:,k))).lt.1.d-20) Cycle
       w = 0.d0
       Do i=1,isol; w(i)=SUM(diag(:,i,ich)*hcp(:,k)); End do
       hcp(:,k) =  w(:)
       if(icase.eq.0) overlaps(nch+jch,ich) = maxval(abs(w(1:isol)))
      End do; End do

      Do i=1,npert; Do j=1,i
       k = ibb(i,j)
       if(hp(k).eq.0.d0) Cycle
       if(icase.eq.0) overlaps(nch+i,nch+j) = abs(hp(k))
      End do; End do

      end if    !   npert

      End Subroutine transform_matrix

