!======================================================================
      Subroutine Gen_matrix
!======================================================================
!     generate the multipole matrix for current itype
!----------------------------------------------------------------------
      Use cmdata, nc => ncdata
      Use dbsr_dmat

      Implicit none
      Integer :: i,ii,nn

! ... prepare the data:
     
      nn = nblk(itype)          ! number of blocks

      if(nn.eq.1.and.jpblk(itype).le.0) nn = 0

      if(nn.eq.0) then          ! nothing to do
       Return  
      elseif(nn.eq.1) then      ! simple case
       i = iblk(itype)
       nc = jpblk(i)-ipblk(i) + 1
       ii = ipblk(i)-1
       Do i=1,nc; IPT(i)=ii+i; End do
      else                      ! need merge data from 
       Do i = 1,nn              ! different blocks
        ii = kblk(itype,i)
        ipi(i) = ipblk(ii)  
        ipj(i) = jpblk(ii)
       End do
       Call Merge_cdata(nn,ipi,ipj,nc,EPS_c)
                                ! release the merged blocks       
       Do i = 1,nn; jpblk(kblk(itype,i))=-1; End do   
      end if

! ... assign the first block for current itype:

      i=itype; iblk(i)=i; nblk(i)=1; kblk(i,1)=i; jpblk(i)=0 

      if(debug.gt.0) write(pri,'(a,2(a,i3),a,i8)') 'Gen_matrix:', &
                      '  itype=',itype, '  nn=',nn, '  nc=',nc

! ... generate matrix:

      Call D_data

      End Subroutine Gen_matrix
