!======================================================================
      Subroutine W_OUT
!======================================================================
!     define the contributions of channels and pertuber configurations
!     for each solution
!----------------------------------------------------------------------
!     Original egenproblem:
!
!     H C = S C E,   E - diagonal egenvalue matrix
!                    S - overlap matrix
!
!----------------------------------------------------------------------
!     Using squire-root matrix:
!
!     C' H C = (C' S C) E  ->  C' S C = 1  =>  (C'S^1/2) (S^1/2 C) = 1
!
!     S A = s A;  S^1/2 = A s(1/2) A'
!
!-----------------------------------------------------------------------
!     Using Cholesky decomposition:
!
!     S = L L' ->     H C = L L' C E
!                    [inv(L) H inv(L')] L'C = L'C E
!                    (L'C)' (L'C) = 1
!
!     supposed solutions are given as (L'C) in matrix A
!-----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: i,i1,i2,ich,is,mwt

      Real(8), Allocatable :: WT(:),cval(:)
      Integer, Allocatable :: ipt(:)

      if(io_processor) then      

! ... local allocations:

      mwt = nch+npert
      Allocate(wt(mwt),cval(khm),ipt(mwt))
      if(allocated(isol)) deallocate(isol); allocate(isol(khm))

! ... open w.nnn if needed: 

      if(iwt.gt.0) then
       i = INDEX(AF_w,'.'); AF = AF_w(1:i)//ALSP
       open(nuw,file=AF,form='UNFORMATTED')
       write(nuw) nch,npert,khm
      end if

! ... open cfg.nnn if needed: 

      if(itype.lt.0.or.cwt.gt.0.d0) then
       i=INDEX(AF_cfg,'.');  AF=AF_cfg(1:i)//ALSP
       Open(nuc,file=AF,status='OLD')
       ncfg=0; Call Read_conf_jj(nuc,0,'add','nocheck')
      end if

      if(cwt.gt.0.d0) write(pri,'(//a)') 'Channel decomposition:'

      end if   ! io_processor

      Call BLACS_BARRIER (ctxt, 'all')

!----------------------------------------------------------------------
! ... define and record weights:

      Do is = 1,khm
       
       call pdgeadd ('notrans', khm, 1, one, z, 1,is,descz, &
                                       zero, v, 1,1, descv)

       call BLACS_BARRIER (ctxt, 'all')

       if(.not.io_processor) Cycle      !  ---> end do

       cval(1:khm) = v(1:khm) * v(1:khm)

       ! ... weights of channels:

       Do ich = 1,nch
        i1=ipsol(ich-1)+1; i2=ipsol(ich); WT(ich)=SUM(cval(i1:i2))
       End do

       ! ... weights of pseudostates:

       if(npert.gt.0) WT(nch+1:nch+npert)=cval(ipsol(nch)+1:khm)

       ! ... find channel with maximum contribution:

       Call SORTa(mwt,WT,ipt);  isol(is) = ipt(1)

       ! record weights for given solution:

       if(iwt.gt.0) write(nuw) WT,isol(is)

       ! print results:

       if(cwt.gt.0.d0.and.is.le.msol) Call Print_w(is,mwt,WT,ipt)

      End do  ! over solutions 'is'

!----------------------------------------------------------------------
! ... local deallocations:

      if(io_processor) Deallocate(wt,cval,ipt)

      call BLACS_BARRIER (ctxt, 'all')

      End Subroutine W_OUT

!======================================================================
      Subroutine Print_w(is,mwt,WT,ipt)
!======================================================================
!     print solution de-composition
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer, intent(in) :: is,mwt,ipt(mwt)
      Real(8), intent(in) :: WT(mwt)
      Integer :: it,ich,jch
      Real(8) :: S,S_ch,S_pt,E_Ry
      Character(64) :: AA

      S=SUM(WT(1:mwt)); S_ch=SUM(WT(1:nch)); S_pt=S-S_ch

      E_Ry = (eval(is)-E_exp(1))*2

      if(S_pt.lt.cwt) &
      write( pri,'(/i5,a,a,f16.8,a,f16.8,10x,a,f10.5)' ) &
         is,'.','  E_au =',eval(is),'  E_Ry =',E_Ry
      if(S_pt.ge.cwt) &
      write( pri,'(/i5,a,a,f16.8,a,f16.8,10x,a,f10.5)' ) &
         is,'.','  E_au =',eval(is),'  E_Ry =',E_Ry,'  C_pt =',S_pt
      write(pri,*)

       Do jch=1,mwt; ich=ipt(jch); if(WT(ich).lt.cwt) Exit

        Call Find_channel_label_jj(ich,jch,is,eval(is),Lab)
        
        if(ich.gt.nch) then
         AA = 'perturber:'
         write(pri,'(a,i6,f9.5,5x,a)') AA(1:12),ich-nch,WT(ich),TRIM(Lab)
        else
         AA = 'continium:'
         it = iptar(ich)
         if(Etarg(it).gt.eval(is))  AA = 'closed ch:'
         write(pri,'(a,i6,f9.5,5x,a)') AA(1:12),ich,WT(ich),TRIM(Lab)
        end if 

       End do

      End Subroutine Print_w




