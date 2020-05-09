!======================================================================
      Subroutine diag_channels
!======================================================================
!     Diagonalize the diagonal channel blocks
!     DIAG array then will contain solutions and energies 
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none

      Real(8) :: ach(ms,ms), bch(ms,ms), w(ms), v(ms), work(3*ms)
      Integer :: ich, k, i,j, ii, jj, ishift, ksol, info, ipos(1), lwork

! ... apply zero conditions:

      Call Zero_cond

! ... diagonalize all blocks:

      if(allocated(ipsol)) Deallocate(ipsol); Allocate(ipsol(0:nch))
      ipsol= 0 
      
      Do ich = 1,nch;  k=icc(ich,ich); if(k.eq.0) Cycle 

       ! ... apply boundary conditions (delete extra B-splines)

       jj=0;  ishift=(ich-1)*ms
       Do j=1,ms
        if(iprm(j+ishift).eq.0) Cycle; jj=jj+1
        ii=0
        Do i=1,ms
         if(iprm(i+ishift).eq.0) Cycle; ii=ii+1
         w(ii)=hch(i,j,k); v(ii)=diag(i,j,ich)  
        End do
        ach(1:ii,jj)=w(1:ii); bch(1:ii,jj)=v(1:ii)
       End do 

       ! ... diagonalization

       lwork = 3*ms
       Call DSYGV(1,'V','L',ii,ach,ms,bch,ms,w,WORK,LWORK,INFO)
       if(info.ne.0) Cycle !Stop 'channel diagonalization failed'

       ksol = 0
       Do i=1,ii
        if(w(i).lt.Edmin) Cycle
        if(w(i).gt.Edmax) Cycle
        if(abs(w(i)).lt.Egap) Cycle
        ksol=ksol+1

        ! ... restore the solutions in original B-spline net:

        v=0.d0; k=0
        Do j=1,ms
         if(iprm(j+ishift).eq.0) Cycle; k=k+1; v(j)=ach(k,i)
        End do 

        ipos=maxloc(abs(v))
        if(v(ipos(1)).lt.0.d0) v=-v

        diag( :,ksol,ich) = v(:)
        diag(ksol,ms,ich) = w(i)       

       End do

       ipsol(ich) = ksol

       if(debug.gt.0.and.pri.gt.0) then
        write(pri,'(/a,a6,2i5/)') 'Channel eigenvalues:',elc(ich),ich,ksol 
        Do i=1,ms; if(w(i).lt.Edmin) Cycle; j=i; Exit; End do
        j = j - 1
        write(pri,'(5e15.7)') w(1:j)
        write(pri,'(1x,74("-"))') 
        write(pri,'(5e15.7)') w(j+1:ms)
       end if

      End do    ! over ich  

      End Subroutine diag_channels


