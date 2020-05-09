!======================================================================
      Subroutine rsol_out
!======================================================================
! ... output the R-matrix solutions and  find the surface amplitudes
!----------------------------------------------------------------------
      Use dbsr_hd     

      Implicit none

      Real(8) :: vb(mhm)
      Integer :: i,j,i1,i2,j1,j2, is,ich 

      if(itype.eq.1.and.io_processor) then
       i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALSP
       Open(nur,file=AF,form='UNFORMATTED')
       write(nur) mhm,khm,nch,npert,ms,nsp,nsq
       write(nur) (eval(i),i=1,khm)
      end if

      if(allocated(WMAT)) deallocate(WMAT);  Allocate(WMAT(nch,khm))
      if(io_processor) &
      write (*  ,'(/a,T20,f10.1,a)') 'WMAT memory:', 4.d0*nch*khm/(1024*1024), ' Mb'

      do is=1,khm

       Call pdgeadd ('notrans', khm, 1, one, z, 1,is,descz, &
                                       zero, v, 1,1, descv)
       Call BLACS_BARRIER (ctxt, 'all')
       if(.not.io_processor) Cycle      

       vb = 0.d0
       Do ich = 1,nch; i1=(ich-1)*ms+1; i2=ich*ms
        j1 = ipsol(ich-1)+1; j2=ipsol(ich)
        Do j=j1,j2
         vb(i1:i2) = vb(i1:i2) + v(j)*bb(1:ms,j)     
        End do
        WMAT(ich,is) = vb(i1+nsp-1)
       End do
       if(npert.gt.0) vb(nch*ms+1:mhm)=v(ksol+1:khm)
       if(itype.eq.1) write(nur) (vb(i),i=1,mhm)

      end do
      
      if(itype.eq.1.and.io_processor) close(nur)
      call BLACS_BARRIER (ctxt, 'all')

      End Subroutine rsol_out


