!=======================================================================
      Subroutine b_out
!=======================================================================
!     output bound-like solutions in bound.nnn
!-----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Real(8), allocatable :: Ebind(:)
      Integer :: i,j,ic,jc,is,js,ic1,ic2,k,n,ich,it,nbound,i1,i2,j1,j2

! ... output file:

      i = len_trim(BF_b)
      write(BF_b(i-2:i),'(i3.3)') klsp
      Open(nuu,file=BF_b,form='UNFORMATTED')

! ... local allocations:

      Allocate(Ebind(khm));  Ebind = 0.d0

!----------------------------------------------------------------------
!                                        define number of bound states:
      nbound = 0
      Do is = 1,khm
       if(eval(is).gt.Emax.and.Emax.ne.0.d0) Cycle
       if(eval(is).lt.Emin.and.Emin.ne.0.d0) Cycle
       ich = isol(is)   ! dominant channel
       it=1; if(ich.le.nch) it=iptar(ich)
       Ebind(is)=eval(is)-E_exp(it)
       nbound = nbound + 1
       if(nbound.ge.msol.and.msol.ne.0) Exit
      End do

!----------------------------------------------------------------------
!                                                  store the solutions:

      write(nuu) mhm,nch,npert,ns,jpar,ipar,nbound

      js = 0
      Do is = 1,khm; if(Ebind(is).eq.0.d0) cycle; js=js+1

       Call Find_channel_label_jj(isol(is),1,is,eval(is),Lab)

       write(nuu) js,LAB

       ich = isol(is)   ! dominant channel
       it=1; if(ich.le.nch) it=iptar(ich)
       write(nuu) eval(is),Ebind(is),ich,it 

       ! find solution:

       v = 0.d0
       Do i = 1,nch; i1=(i-1)*ms+1; i2=i*ms
        j1 = ipsol(i-1)+1; j2=ipsol(i)
        Do j=j1,j2
         v(i1:i2) = v(i1:i2) + a(j,is)*bb(1:ms,j)
        End do
       End do
       if(npert.gt.0) v(nch*ms+1:mhm)=a(ksol+1:khm,is)

       write(nuu) v(1:mhm)

       if(js.ge.nbound) Exit

      End do

      write(nuu) khm
      Ebind(1:khm) = (eval(1:khm)-etarg(1))*2
      write(nuu) Ebind(1:khm)

      Close(nuu)
      Deallocate(Ebind)

      End Subroutine b_out




