!======================================================================
      Subroutine alloc_dbsr_matrix
!======================================================================
!     allocate basic matrices
!----------------------------------------------------------------------
      Use dbsr_mat

      Integer :: kprocs, i,k,ich,jch

      if(nch.eq.0.and.npert.eq.0) Return

      mhm = nch*ms + npert 
      kprocs=nprocs-1

      if(Allocated(icc)) Deallocate(icc) 
      Allocate(icc(nch,nch)); icc = 0; iicc = 0
      mem_mat = 4.d0*nch*nch      

      k = 0; i=0

      if(idiag.ge.0) then
      Do ich=1,nch
       k=k+1; if(k.gt.kprocs) k=1; if(nprocs.eq.1) k=0
       if(myid.ne.k) Cycle       
       i=i+1; icc(ich,ich) = i; icc(ich,ich) = i
      End do
      iicc = nch
      end if

      if(idiag.le.0) then
      Do ich=2,nch; Do jch=1,ich-1
       k=k+1; if(k.gt.kprocs) k=1; if(nprocs.eq.1) k=0
       if(myid.ne.k) Cycle       
       i=i+1; icc(ich,jch) = i; icc(jch,ich) = i
      End do; End do
      iicc = i
      end if
 
      if(Allocated(hch)) Deallocate(hch) 
      Allocate(hch(ms,ms,iicc)); hch=0.d0
      if(Allocated(acf)) Deallocate(acf) 
      Allocate(acf(iicc ,0:mk)); acf=0.d0

      mem_mat = mem_mat + 8.d0*ms*ms*iicc + 4.d0*iicc*(mk+1)      

      if(npert.eq.0.or.idiag.eq.1) go to 1  

      if(Allocated(hcp)) Deallocate(hcp,icb) 
      Allocate(icb(nch,npert)); icb=0
      k = 0; i=0
      Do ich=1,nch; Do jch=1,npert
       k=k+1; if(k.gt.kprocs) k=1; if(nprocs.eq.1) k=0
       if(myid.ne.k) Cycle       
       i=i+1; icb(ich,jch) = i
      End do; End do
      iicb = i
      Allocate(hcp(ms,iicb)); hcp=0.d0
      mem_mat = mem_mat + 8.d0*ms*iicb

      if(Allocated(hp)) Deallocate(hp,ibb) 
      Allocate(ibb(npert,npert)); ibb = 0
      k = 0; i=0
      Do ich=1,npert; Do jch=1,ich
       k=k+1; if(k.gt.kprocs) k=1; if(nprocs.eq.1) k=0
       if(myid.ne.k) Cycle       
       i=i+1; ibb(ich,jch) = i; ibb(jch,ich) = i
      End do; End do
      iibb = i
      Allocate(hp(iibb)); hp=0.d0
      mem_mat = mem_mat + 8.d0*iibb

    1 mem_mat = mem_mat / (1024*1024)

      End Subroutine alloc_dbsr_matrix

!======================================================================
      Subroutine UPDATE_HX(ich,jch,d,C)
!======================================================================
!     update channel block (ms,ms) with full matrix
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ich,jch
      Real(8), intent(in) :: C,d(ms,ms)
      Real(8) :: xx(ms,ms)
      Integer ::  i

      if(ich.lt.0)   Stop 'UPDATE_HX: ich < 0'
      if(jch.lt.0)   Stop 'UPDATE_HX: jch < 0'
      if(ich.gt.nch) Stop 'UPDATE_HX: ich > nch'
      if(jch.gt.nch) Stop 'UPDATE_HX: jch > nch'

      if(ich.ge.jch) then
       xx = C * d
      else
       xx = C * TRANSPOSE(d)
      end if      

      i = icc(ich,jch)   
      if(i.lt.0.or.i.gt.iicc) Stop 'UPDATE_HX: i out of range'
      hch(:,:,i) = hch(:,:,i) + xx
     
      END Subroutine UPDATE_HX


!======================================================================
      Subroutine UPDATE_HS(ich,jch,d,sym,ii,jj)
!======================================================================
!     update 'small' channel block (ns,ns) 
!     sym = 's'  -->  symmetric banded upper-column storage mode
!     sym = 'n'  -->  non-symmetric band matrix  
!     sym = 'x'  -->  non-symmetric full matrix
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ich,jch,ii,jj
      Real(8), intent(in) :: d(ns,ns)
      Real(8) :: y(ns,ns),yy(ns,ns)
      Character, intent(in) :: sym
      Integer ::  i,j, jp, imin,imax, i1=1,i2=1,j1=1,j2=1

      Select case(sym)

      Case('s')

       y = 0.d0
       Do jp=1,ks;  Do i=1,ns-jp+1;  j=i+jp-1
        y(i,j) = d(i,jp)
        y(j,i) = d(i,jp)
       End do; End do

      Case('n')

       y = 0.d0
       Do jp = 1,ks+ks-1
        imin=max0( 1, 1 + ks-jp)
        imax=min0(ns,ns + ks-jp)
        Do i = imin,imax;  j=i+jp-ks
         y(i,j) = d(i,jp)
        End do
       End do

      Case('x')

       y = d

      End Select

      if(ich.lt.0)   Stop 'UPDATE_HS: ich < 0'
      if(jch.lt.0)   Stop 'UPDATE_HS: jch < 0'
      if(ich.gt.nch) Stop 'UPDATE_HS: ich > nch'
      if(jch.gt.nch) Stop 'UPDATE_HS: jch > nch'

      i = icc(ich,jch)
      if(i.lt.0.or.i.gt.iicc) Stop 'UPDATE_HX: i out of range'

      i1 = 1+(ii-1)*ns; i2=i1+ns-1     
      j1 = 1+(jj-1)*ns; j2=j1+ns-1

      if(ich.gt.jch) then
       hch(i1:i2,j1:j2,i) = hch(i1:i2,j1:j2,i) + y

      elseif(jch.gt.ich) then
       yy = TRANSPOSE(y)
       hch(i1:i2,j1:j2,i) = hch(i1:i2,j1:j2,i) + yy

      else
       y = y/2.d0; yy = TRANSPOSE(y)
       hch(i1:i2,j1:j2,i) = hch(i1:i2,j1:j2,i) + y
       hch(j1:j2,i1:i2,i) = hch(j1:j2,i1:i2,i) + yy

      end if      

     
      End Subroutine UPDATE_HS


!======================================================================
      Subroutine UPDATE_HW (ich,jch,v,w,C)
!======================================================================
!     update matrix  by  direct product of vectors:  v * w^T
!-----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ich,jch
      Real(8), intent(in) :: C
      Real(8), intent(in) :: v(ms),w(ms)
      Real(8) :: x(ms,ms),xx(ms,ms)
      Integer ::  i,j

      if(ich.lt.0)   Stop 'UPDATE_HW: ich < 0'
      if(jch.lt.0)   Stop 'UPDATE_HW: jch < 0'
      if(ich.gt.nch) Stop 'UPDATE_HW: ich > nch'
      if(jch.gt.nch) Stop 'UPDATE_HW: jch > nch'

      Do i=1,ms; Do j=1,ms; x(i,j)=v(i)*w(j); End do; End do      

      if(ich.gt.jch) then
       xx = C * x
      elseif(jch.gt.ich) then
       xx = C * TRANSPOSE(x)
      else
       xx = TRANSPOSE(x)
       xx = (C/2.d0)*(x+xx)
      end if

      i=icc(ich,jch)
      if(i.lt.0.or.i.gt.iicc) Stop 'UPDATE_HX: i out of range'
      hch(:,:,i) = hch(:,:,i) + xx

      End Subroutine UPDATE_HW


!======================================================================
      Subroutine UPDATE_HV(ich,ip,v,C)
!======================================================================
!     update vector entry  (channel-bound)
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ich,ip
      Real(8), intent(in) :: C
      Real(8), intent(in) :: v(ms)
      Integer :: i

      if(ich.lt.1.or.ich.gt.nch) &
       Stop 'UPDATE_HV: channel index out is of range'
      if(ip.lt.1.or.ip.gt.npert) &
       Stop 'UPDATE_HV: bound index out is of range'

      i = icb(ich,ip)
      if(i.lt.0.or.i.gt.iicb) Stop 'UPDATE_HV: ich, ip out of range'

      hcp(:,i) = hcp(:,i) + C*v(:)

      End Subroutine UPDATE_HV

!======================================================================
      Subroutine UPDATE_HB(ic,jc,c)
!======================================================================
!     update scalar element (bound-bound)
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ic,jc
      Real(8), intent(in) :: c
      Integer :: i

      if(ic.gt.npert.or.jc.gt.npert.or.ic.le.0.or.jc.le.0) &
       Stop 'UPDATE_HB: index is out of range'

      i=ibb(ic,jc)
      if(i.lt.0.or.i.gt.iibb) Stop 'UPDATE_HB: ic, jc out of range'

      hp(i) = hp(i) + c

      End Subroutine UPDATE_HB


!=======================================================================
      Subroutine UPDATE_ACF(k,ich,jch,C)
!=======================================================================
!     Update asymptotic coefficients
!-----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: k, ich,jch
      Real(8), intent(in) :: C
      Integer :: i

      if(k.gt.mk) Stop 'Update_cf: k > kmk '
      if(k.lt.0)  Stop 'Update_cf: k < 0'

      if(ich.lt.0)   Stop 'Update_cf: ich < 0'
      if(jch.lt.0)   Stop 'Update_cf: jch < 0'
      if(ich.gt.nch) Stop 'Update_cf: ich > nch'
      if(jch.gt.nch) Stop 'Update_cf: jch > nch'

      i = icc(ich,jch)
      if(i.lt.0.or.i.gt.iicc) Stop 'UPDATE_HX: i out of range'

      acf(i,k) = acf(i,k) + C

      End Subroutine UPDATE_ACF


!=======================================================================
      Subroutine Print_matrix
!=======================================================================
      Use dbsr_mat

      Integer :: i,j,k, ich,jch

      Do ich = 1,nch
       Do jch = 1,ich
        k = icc(ich,jch); if(k.eq.0) Cycle
        if(idiag.eq.1.and.ich.ne.jch) Cycle
        if(idiag.eq.-1.and.ich.eq.jch) Cycle
        write(pri,'(/a,2i5/)') 'ich,jch', ich,jch
        Do i = 1,10
         write(pri,'(10E15.5)') (hch(i,j,k),j=1,10) 
        End do
       End do
      End do

      if(npert.eq.0) Return

      Do ich = 1,nch
       Do jch = 1,npert
        k = icb(ich,jch); if(k.eq.0) Cycle
        if(idiag.eq.1) Cycle
        write(pri,'(/a,2i5/)') 'ich,ip', ich,jch
        write(pri,'(10E15.5)') hcp(1:10,k) 
       End do
      End do

      Do ich = 1,npert
       Do jch = 1,npert
        k = ibb(ich,jch); if(k.eq.0) Cycle
        if(idiag.eq.1) Cycle
        write(pri,'(/a,2i5/)') 'ip,jp', ich,jch
        write(pri,'(10E15.5)') hp(k) 
       End do
      End do



      End Subroutine Print_matrix


