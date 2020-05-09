!======================================================================
      Subroutine Read_matrix
!======================================================================
!     read transformed Hamiltonian matrix 
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: i,j, ic,jc, idim,jdim, info

!--------------------------------------------------------------------------
! ... initialize array descriptor for the interaction matrices:
! ... np, nq are matrix size on current processor from modul bsr_hd

      Call descinit(desca,khm,khm,nblock,nblock,rsrc,csrc,ctxt,ld,info)
      Call p_error(info,'desca descriptor error in read_matrix')

      if(allocated(a)) deallocate(a); allocate(a(np,nq));  a=zero

!--------------------------------------------------------------------------
! ... put diagonal blocks:

      Do ic = 1,nch
       i=ipsol(ic-1)+1; idim=ipsol(ic)-ipsol(ic-1)
       if(io_processor) then
        jc=ipsol(ic-1)
        add = zero; Do j=1,idim; add(j,j)=bval(jc+j); End do
       end if
       Call BLACS_BARRIER (ctxt, 'all')
       Call pdgeadd ('notrans', idim, idim,               &
                                 one,  add, 1,1, descadd, &
                                zero,   a,  i,i, desca    )
      End do

      if(npert.gt.0) adp = zero

!---------------------------------------------------------------------------
! ... read non-trivial blocks

      Do 
       if(io_processor) read(nui) ic,jc
       Call br_ipar(ic)        
       Call br_ipar(jc)        
       Call BLACS_BARRIER (ctxt, 'all')
       if(ic.le.0) Exit

      if(ic.gt.nch.and.jc.gt.nch) then       !  pert-pert

       if(io_processor) read(nui) adp(ic-nch,jc-nch) 

      elseif(ic.gt.nch) then                 !  ch-pert

       jdim = ipsol(jc)-ipsol(jc-1)
       i=ic-nch+ipsol(nch); j=ipsol(jc-1)+1
       if(io_processor)  read(nui) add(1,1:jdim)
       Call pdgeadd ('notrans', 1,jdim, one, add, 1,1, descadd, &
                                        zero, a,  i,j, desca  )

      else                                   !  ch-ch

       i=ipsol(ic-1)+1; idim=ipsol(ic)-ipsol(ic-1) 
       j=ipsol(jc-1)+1; jdim=ipsol(jc)-ipsol(jc-1) 
       if(io_processor) read(nui) add(1:idim,1:jdim)
       Call pdgeadd ('notrans',idim,jdim, one,add,  1,1, descadd,&
                                          zero,  a, i,j, desca )
      end if

       Call BLACS_BARRIER (ctxt, 'all')
      End do

      if(npert.gt.0) then
       i=ipsol(nch)+1
       Call pdgeadd('notrans',npert,npert, one, adp, 1,1, descadp, &
                                           zero,  a, i,i, desca )
      end if
      Call BLACS_BARRIER (ctxt, 'all')

      End Subroutine Read_matrix
