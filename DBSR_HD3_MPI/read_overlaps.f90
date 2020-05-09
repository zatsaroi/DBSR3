!======================================================================
      Subroutine Read_overlaps
!======================================================================
!     read transformed overlap matrix 
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: i,j, ic,jc, idim,jdim, info

!----------------------------------------------------------------------
! ... initialize array descriptor for the overlap matrix:

      Call descinit(descb,khm,khm,nblock,nblock,rsrc,csrc,ctxt,ld,info)
      Call p_error(info,'descb descriptor error in read_overlap')

! ... np, nq are matrix size on current processor from module bsr_hd

      if(allocated(b)) deallocate(b); allocate (b(np,nq));  b=zero

!----------------------------------------------------------------------
! ... put unit chanel's diagonal blocks:

      add = zero; Do i=1,ms; add(i,i)=one; End do
      Do i = 1,nch
       j=ipsol(i-1)+1; idim=ipsol(i)-ipsol(i-1)
       Call pdgeadd ('notrans', idim,idim,       &
                      one, add,  1,  1, descadd, &
                      zero,  b,  j,  j, descb    )
       Call BLACS_BARRIER (ctxt, 'All')
      End do

      if(npert.ne.0) adp = zero
!---------------------------------------------------------------------------
! ... check if we have non-trivial overlap matrix:

      if(io_processor) then
       read(nui) diag_ovl
       if(diag_ovl.ne.0) Backspace(nui)
      end if
      Call br_ipar(diag_ovl)

      if(diag_ovl.eq.0) Return       

!---------------------------------------------------------------------------
! ... read not-trivial blocks:

      Do 
       if(io_processor) read(nui) ic,jc
       Call br_ipar(ic)        
       Call br_ipar(jc)        
       Call BLACS_BARRIER (ctxt, 'All')
       if(ic.le.0) Exit

       if(ic.gt.nch.and.jc.gt.nch) then       !  pert-pert
  
        if(io_processor) read(nui) adp(ic-nch,jc-nch)
  
       elseif(ic.gt.nch) then                 !  ch-pert
  
        jdim = ipsol(jc)-ipsol(jc-1)
        if(io_processor) read(nui) add(1,1:jdim)
        i=ic-nch+ipsol(nch);  j=ipsol(jc-1)+1
        Call pdgeadd ('notrans', 1,jdim,  one, add,  1,1, descadd, &
                                         zero,   b,  i,j, descb    )
  
       else                                   !  ch-ch
  
        i=ipsol(ic-1)+1; idim=ipsol(ic)-ipsol(ic-1) 
        j=ipsol(jc-1)+1; jdim=ipsol(jc)-ipsol(jc-1) 
        if(io_processor) read(nui) add(1:idim,1:jdim)
        Call pdgeadd ('notrans',idim,jdim, one , add, 1,1, descadd,&
                                           zero,   b, i,j, descb  )
       end if
  
        Call BLACS_BARRIER (ctxt, 'All')
       End do
  
       if(npert.gt.0) then
        i=ipsol(nch)+1
        Call pdgeadd('notrans',npert,npert, one, adp, 1,1, descadp, &
                                           zero,   b, i,i, descb )
       end if

      Call BLACS_BARRIER (ctxt, 'All')

      End Subroutine Read_overlaps


