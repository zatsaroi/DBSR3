!======================================================================
      Subroutine Record_matrix 
!======================================================================
!     record the overlap or Hamiltonian matrix to unit nui
!----------------------------------------------------------------------
      Use MPI
      Use dbsr_mat      

      Implicit none
      Integer :: i,j, mm, ni,nj, ip,jp, ich,jch
      Real(8) :: S, w(ms), v(ns), a(ms,ms), c(ns,ns)
      Integer :: status(MPI_STATUS_SIZE)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! ... non-diagonal blocks:
 
      mm = ms*ms
      Do ich = 2,nch;    ni=ipsol(ich)
      Do jch = 1,ich-1;  nj=ipsol(jch)
       i=icc(ich,jch) 
       if(i.ne.0) Call MPI_SEND(hch(:,:,i), mm, MPI_DOUBLE_PRECISION, &
                                0, i, MPI_COMM_WORLD, ierr)
       if(myid.eq.0) then
        Call MPI_RECV(a, mm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        if(maxval(abs(a)).gt.1.d-20) then
         c = 0.d0
         Do i=1,ni
          Do j=1,ms; w(j)  = SUM(diag(:,i,ich)*a(:,j)); End do 
          Do j=1,nj; c(i,j)= SUM(w(:)*diag(:,j,jch));       End do 
         End do       
         write(nui) ich,jch
         write(nui) c(1:ni,1:nj)
         if(icase.eq.0) overlaps(ich,jch) = maxval(abs(c(1:ni,1:nj)))
        end if
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End do; End do

! ... pertubers:

      if(npert.gt.0) then

! ... channel-perturber rows:

       Do ich = 1,nch;    ni = ipsol(ich)
       Do ip  = 1,npert;  i = icb(ich,ip)
                                          
        if(i.ne.0) Call MPI_SEND(hcp(:,i),ms, MPI_DOUBLE_PRECISION, &
                                 0, i, MPI_COMM_WORLD, ierr)
        if(myid.eq.0) then
         Call MPI_RECV(w, ms, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         S = maxval(abs(w))
         if(maxval(abs(w)).gt.1.d-20) then
          v = 0.d0
          Do i=1,ni; v(i)=SUM(diag(:,i,ich)*w(:)); End do
          if(icase.eq.0) overlaps(nch+ip,ich) = maxval(abs(v(1:ni)))
          write(nui) ip+nch,ich
          write(nui) v(1:ni)
         end if
        end if       
      
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do; End do
       
! ... perturter-perturber elements:

       Do ip = 1,npert; Do jp = 1,ip; i = ibb(ip,jp)
                                     
       if(i.ne.0) &
        Call MPI_SEND(hp(i),1, MPI_DOUBLE_PRECISION, &
                      0, i, MPI_COMM_WORLD, ierr)

       if(myid.eq.0) then
        Call MPI_RECV(S, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        if(S.ne.0.d0) then
         write(nui) ip+nch,jp+nch
         write(nui) S
        if(icase.eq.0) overlaps(nch+ip,nch+jp) = abs(S)
        end if
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do; End do

      end if   ! over npert

      if(myid.eq.0) write(nui) 0,0  ! sign of the end

      End Subroutine Record_matrix
