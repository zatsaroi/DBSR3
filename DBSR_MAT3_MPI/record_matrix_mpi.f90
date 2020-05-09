!======================================================================
      Subroutine Record_matrix 
!======================================================================
!     record the overlap or Hamiltonian matrix to unit nui
!----------------------------------------------------------------------
      Use MPI
      Use dbsr_mat      

      Implicit none
      Integer :: i, ni,nj, ip,jp, ich,jch
      Real(8) :: S, v(ms)
      Integer :: status(MPI_STATUS_SIZE)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! ... non-diagonal blocks:

      Do ich = 2,nch;    ni=ipsol(ich)
      Do jch = 1,ich-1;  nj=ipsol(jch)
       i=icc(ich,jch) 
       if(i.ne.0) &
        Call MPI_SEND(hch(1:ni,1:nj,i),ni*nj, MPI_DOUBLE_PRECISION, &
                      0, i, MPI_COMM_WORLD, ierr)
       if(myid.eq.0) then
        Call MPI_RECV(x(1:ni,1:nj), ni*nj, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        S = maxval(abs(x(1:ni,1:nj)))
        if(S.gt.1.d-20) then
         write(nui) ich,jch
         write(nui) x(1:ni,1:nj)
        end if
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End do; End do

! ... pertubers:

      if(npert.gt.0) then

! ... channel-perturber rows:

       Do ich = 1,nch;    ni = ipsol(ich)
       Do ip  = 1,npert;  i = icb(ich,ip)
                                          
       if(i.ne.0) &
        Call MPI_SEND(hcp(1:ni,i),ni, MPI_DOUBLE_PRECISION, &
                      0, i, MPI_COMM_WORLD, ierr)
       if(myid.eq.0) then
        Call MPI_RECV(v, ni, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        S = maxval(abs(v(1:ni)))
        if(S.gt.1.d-20) then
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
        end if
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do; End do

      end if   ! over npert

      if(myid.eq.0) write(nui) 0,0  ! sign of the end

      End Subroutine Record_matrix
