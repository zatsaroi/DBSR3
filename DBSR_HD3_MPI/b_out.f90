!=======================================================================
      Subroutine b_out 
!=======================================================================
!     output the bound solutions
!-----------------------------------------------------------------------
      Use dbsr_hd
      
      Implicit none
     
      Real(8), Allocatable :: Ebind(:), vb(:)
      Integer :: i,j,i1,i2,j1,j2,ich,is,js,it,nbound

      if(allocated(vb)) Deallocate(vb); Allocate(vb(mhm))
      if(allocated(Ebind)) Deallocate(Ebind); Allocate(Ebind(khm));  Ebind = 0.d0
!----------------------------------------------------------------------
! ... define number of bound states for output:

      if(io_processor) then 

! ... output file:

      i = len_trim(BF_b)
      write(BF_b(i-2:i),'(i3.3)') klsp
      Open(nuu,file=BF_b,form='UNFORMATTED')

! ... define number of bound states:

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

! ... record the bound-states parameters:

      write(nuu) mhm,nch,npert,ns,jpar,ipar,nbound

      end if  ! io_processor

! ... broadcast Ebind and nbound:

      if(io_processor) then
       Call igebs2d (ctxt, 'all', ' ', khm, 1, Ebind, khm)
      else
       Call igebr2d (ctxt, 'all', ' ', khm, 1, Ebind, khm, rsrc, csrc)
      end if    

      Call br_ipar(nbound)

      Call BLACS_BARRIER (ctxt, 'all')

!----------------------------------------------------------------------
!                                                  store the solutions:
       js = 0
       Do is = 1,khm

        if(Ebind(is).eq.0.d0) Cycle; js=js+1; if(js.gt.nbound) Exit

        Call pdgeadd ('notrans', khm, 1, one, z, 1,is, descz, &
                                        zero, v, 1, 1, descv)

        if(.not.io_processor) Cycle

         Call Find_channel_label_jj(isol(is),1,is,eval(is),Lab)

         write(nuu) js,LAB
         ich = isol(is)   ! dominant channel
         it=1; if(ich.le.nch) it=iptar(ich)
         write(nuu) eval(is),Ebind(is),ich,it 

         ! find solution:

         vb = 0.d0
         Do ich = 1,nch; i1=(ich-1)*ms+1; i2=ich*ms
          j1 = ipsol(ich-1)+1; j2=ipsol(ich)
          Do j=j1,j2
           vb(i1:i2) = vb(i1:i2) + v(j)*bb(1:ms,j)     
          End do
         End do
         if(npert.gt.0) vb(nch*ms+1:mhm)=v(ksol+1:khm)
         write(nuu) vb(1:mhm)

        End do ! is

       Deallocate(Ebind, vb)

       Call BLACS_BARRIER (ctxt, 'all')

       End Subroutine b_out



