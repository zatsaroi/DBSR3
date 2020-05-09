!======================================================================
      Subroutine UPDATE_DX (ich,jch,C)
!======================================================================
!     Update channel block with multipole integral 
!----------------------------------------------------------------------
      Use dbsr_dmat
      Use channels_jj, only: kch

      Implicit none
      Integer, intent(in) :: ich,jch
      Real(8), intent(in) :: C

      Integer :: i,j
      Real(8) :: S

      if(ktype.eq.'M') then

       S = C * (kch(ilsp1,ich)+kch(ilsp2,jch))
       if(S.eq.0.d0) Return
       S = S / (kpol+1)

       i = (ich-1)*ms
       j = (jch-1)*ms + ns
       DL(i+1:i+ns,j+1:j+ns) = DL(i+1:i+ns,j+1:j+ns) + dipLp*S 
       i = (ich-1)*ms + ns
       j = (jch-1)*ms 
       DL(i+1:i+ns,j+1:j+ns) = DL(i+1:i+ns,j+1:j+ns) + dipLq*S  

      else

       i = (ich-1)*ms
       j = (jch-1)*ms 
       DL(i+1:i+ns,j+1:j+ns) = DL(i+1:i+ns,j+1:j+ns) + dipLp*C
       i = (ich-1)*ms + ns
       j = (jch-1)*ms + ns
       DL(i+1:i+ns,j+1:j+ns) = DL(i+1:i+ns,j+1:j+ns) + dipLq*C 

       S = C*(kch(ilsp1,ich)-kch(ilsp2,jch)-kpol)
       i = (ich-1)*ms
       j = (jch-1)*ms + ns
       DV(i+1:i+ns,j+1:j+ns) = DV(i+1:i+ns,j+1:j+ns) + dipVp*S 

       S = C*(kch(ilsp1,ich)-kch(ilsp2,jch)+kpol)
       i = (ich-1)*ms + ns
       j = (jch-1)*ms 
       DV(i+1:i+ns,j+1:j+ns) = DV(i+1:i+ns,j+1:j+ns) + dipVq*S   

      end if
 
      End Subroutine UPDATE_DX


!======================================================================
      Subroutine UPDATE_DS(ich,jch,CL,CV)
!======================================================================
!     Update channel block for overlaps  
!----------------------------------------------------------------------
      Use dbsr_dmat
      Use DBS_gauss, only: fppqq

      Implicit none
      Integer, intent(in) :: ich,jch
      Real(8), intent(in) :: CL,CV

      Integer ::  i,j

      i = (ich-1)*ms
      j = (jch-1)*ms
      DL(i+1:i+ms,j+1:j+ms) = DL(i+1:i+ms,j+1:j+ms) + CL*fppqq  
      if(ktype.ne.'M') &
      DV(i+1:i+ms,j+1:j+ms) = DV(i+1:i+ms,j+1:j+ms) + CV*fppqq  

      End Subroutine UPDATE_DS


!======================================================================
      Subroutine UPDATE_DW (ich,jch,v1,w1,v2,w2,C)
!======================================================================
!     update matrix by v * w' 
!-----------------------------------------------------------------------
      Use dbsr_dmat

      Implicit none
      Integer, intent(in) :: ich,jch
      Real(8), intent(in) :: v1(ms),w1(ms),v2(ms),w2(ms),C
      Integer :: i,j
      Real(8) :: x1(ms,ms),x2(ms,ms),S1,S2

      Do j=1,ms  
       S1=C*w1(j); x1(:,j)=v1(:)*S1
       S2 = 0.d0; if(ktype.ne.'M') S2=C*w2(j); x2(:,j)=v2(:)*S2
      End do      

      i=(ich-1)*ms
      j=(jch-1)*ms
      DL(i+1:i+ms,j+1:j+ms)=DL(i+1:i+ms,j+1:j+ms)+x1(1:ms,1:ms)
      if(ktype.ne.'M') &
      DV(i+1:i+ms,j+1:j+ms)=DV(i+1:i+ms,j+1:j+ms)+x2(1:ms,1:ms)

      End Subroutine UPDATE_DW


!======================================================================
      Subroutine UPDATE_DV(ich,jch,ic,jc,v1,v2,C)
!======================================================================
!     update vector
!----------------------------------------------------------------------
      Use dbsr_dmat

      Implicit none
      Integer, intent(in) :: ich,jch,ic,jc
      Real(8), intent(in) :: C,v1(ms),v2(ms)
      Real(8) :: w1(ms),w2(ms)
      Integer :: i,j

      w1 = C * v1
      w2 = C * v2
      if(ich.gt.0.and.jch.eq.0) then
       i = (ich-1)*ms
       j = kdm2 + jc-npert2
       DL(i+1:i+ms,j) = DL(i+1:i+ms,j) + w1(1:ms)
       if(ktype.ne.'M') &
       DV(i+1:i+ms,j) = DV(i+1:i+ms,j) + w2(1:ms)
      elseif(jch.gt.0.and.ich.eq.0) then
       j = (jch-1)*ms
       i = kdm1 + ic-npert1
       DL(i,j+1:j+ms) = DL(i,j+1:j+ms) + w1(1:ms)
       if(ktype.ne.'M') &
       DV(i,j+1:j+ms) = DV(i,j+1:j+ms) + w2(1:ms)
      else
       Stop ' UPDATE_DV: non-allowed combination of ich,jch' 
      end if

      End Subroutine UPDATE_DV


!======================================================================
      Subroutine UPDATE_DB(ic,jc,CL,CV)
!======================================================================
!     update scalar
!----------------------------------------------------------------------
      Use dbsr_dmat

      Implicit none
      Integer, intent(in) :: ic,jc
      Real(8), intent(in) :: CL,CV
      Integer ::  i,j

      if(ic.le.0) Stop 'UPDATE_DB: ic <= 0'
      if(jc.le.0) Stop 'UPDATE_DB: jc <= 0'

      if(ic.gt.npert1) then
       write(*,'(a,2i5)') 'UPDATE_DB: ic > npert1',ic,ncfg1
       Stop 
      end if
      if(jc.gt.npert2) then
       write(*,'(a,2i5)') 'UPDATE_DB: jc > npert2',jc,ncfg2
       Stop 
      end if

      i = kdm1 + ic - npert1
      j = kdm2 + jc - npert2
      DL(i,j) = DL(i,j) + CL
      if(ktype.ne.'M') &
      DV(i,j) = DV(i,j) + CV

      End Subroutine UPDATE_DB

