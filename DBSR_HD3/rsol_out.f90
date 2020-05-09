!======================================================================
      Subroutine Rsol_out
!======================================================================
! ... output the R-matrix solutions and
! ... find the surface amplitudes
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: i,j,i1,i2,j1,j2, is,ich

      if(itype.eq.1) then
       i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALSP
       Open(nur,file=AF,form='UNFORMATTED')
       write(nur) mhm,khm,nch,npert,ms,nsp,nsq
       write(nur) eval(1:khm)
      end if

      if(allocated(WMAT)) Deallocate(WMAT);  Allocate(WMAT(nch,khm))

      Do is=1,khm
       v = 0.d0
       Do ich = 1,nch; i1=(ich-1)*ms+1; i2=ich*ms
        j1 = ipsol(ich-1)+1; j2=ipsol(ich)
        Do j=j1,j2
         v(i1:i2) = v(i1:i2) + a(j,is)*bb(1:ms,j)
        End do
        WMAT(ich,is) = v(i1-1+nsp)
       End do
       if(npert.gt.0) v(nch*ms+1:mhm)=a(ksol+1:khm,is)
       if(itype.eq.1) write(nur) (v(i),i=1,mhm)
      End do

      if(itype.eq.1) close(nur)

      End Subroutine Rsol_out


