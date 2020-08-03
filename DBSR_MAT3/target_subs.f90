!======================================================================
      Subroutine Target_h(ich,jch,C,CC)
!======================================================================
!     update target interaction and overlap matrixes, by considering 
!     the terms wjp structure <kl|k'l> <target|H|target'>
!
!     It is not pure target states, but the basis states before |kl>,
!     i.e. the basis states may repeat, if one target state can
!     couple to several kl.
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ich,jch
      Real(8), intent(in) :: c, cc
      Integer :: i,j,ij,it,jt

      if(ich.lt.0)   Stop 'target_h: ich < 0'
      if(jch.lt.0)   Stop 'target_h: jch < 0'
      if(ich.gt.nch) Stop 'target_h: ich > nch'
      if(jch.gt.nch) Stop 'target_h: jch > nch'

      if(icc(ich,jch).eq.0) Return

      it = iptar(ich); jt = iptar(jch)
      i=max(it,jt);  j=min(it,jt); ij=(i-1)*i/2+j

      htarg(ij) = htarg(ij) + C
      otarg(ij) = otarg(ij) + CC

      End Subroutine Target_h


!======================================================================
      Subroutine Target_print
!======================================================================
!     debug printing of target matrix
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer :: i,j,ij
      Real(8) :: C

      Call Target_norm

      if(iitar.eq.0)  then
       write(pri,'(/a,1Pe9.1)') 'Target hamiltonian errors if > eps_tar =',eps_tar
       Do i = 1,ntarg;  Do j = 1,i
        ij=(i-1)*i/2+j
        if(itarget(ij).eq.0) Cycle
        if(i.eq.j) then
         htarg(ij) = htarg(ij) + EC
         C = htarg(ij) - Etarg(i)
         if(abs(C).gt.eps_tar)  &
         write(pri,'(2i5,2f16.6,5x,a)') i,j,htarg(ij),C,trim(AFT(i))
        else
         c = abs(htarg(ij))
         if(abs(C).gt.eps_tar) &
         write(pri,'(2i5,f16.6,5x,a,5x,a)') i,j,htarg(ij),trim(AFT(i)),trim(AFT(j))
        end if
       End do;  End do
      end if

      if(iitar.eq.0)  then
       write(pri,'(/a,1Pe9.1)') 'Target overlaps errors    if > eps_tar =',eps_tar
       Do i = 1,ntarg;  Do j = 1,i
        ij=(i-1)*i/2+j
        if(itarget(ij).eq.0) Cycle
        C = otarg(ij)
        if(i.eq.j) c = c - 1.d0
        if(abs(C).gt.eps_tar) &
         write(pri,'(2i5,f16.6,5x,a,5x,a)') i,j,C,trim(AFT(i)),trim(AFT(j))
       End do; End do
      end if

      End Subroutine Target_print


!======================================================================
      Subroutine Target_norm
!======================================================================
!     normalized target energies 
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer :: i,j,ij,it,jt,ilen,k,info
      Real(8) :: c

      if(allocated(itarget)) Deallocate(itarget)
      Allocate( itarget(ntarg*(ntarg+1)/2) );  itarget = 0

      Do i = 1,nch; it=iptar(i)
      Do j = 1,i;   jt=iptar(j)
       if(it.eq.jt.and.i.ne.j) Cycle
       ij=(it-1)*it/2+jt
       itarget(ij)=itarget(ij) + 1
      End do;  End do

      Do i = 1,ntarg;  Do j = 1,i
       ij=(i-1)*i/2+j
       if(itarget(ij).eq.0) Cycle
       htarg(ij) = htarg(ij)/itarget(ij)
       otarg(ij) = otarg(ij)/itarget(ij)
      End do;  End do

      End Subroutine Target_norm


!======================================================================
      Subroutine Target_new1
!======================================================================
!     new target energies 
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer :: i,j,ij,it,jt,k,info
      Real(8) :: C
      Real(8), allocatable :: htarget(:,:), evalt(:)

      Allocate(htarget(ntarg,ntarg), evalt(ntarg))

      htarget=0.d0 
      Do it=1,ntarg; htarget(it,it)=Etarg(it); End do

! ... target matrix:

      k = 0
      Do it = 1,ntarg;   Do jt = 1,it
       if(it.eq.jt) Cycle
       ij=(it-1)*it/2+jt
       C = abs(htarg(ij))
       if(abs(C).lt.eps_tar) Cycle
       k = k + 1
       htarget(it,jt) = C; htarget(jt,it) = C
      End do;  End do

      if(k.eq.0) Return

      if(debug.gt.0) then
       write(pri,'(/a/)') 'Target hamiltonian:'
       Do i=1,ntarg
        write(pri,'(20f10.5)') htarget(i,1:i)
       End do
      end if

      Call LAP_DSYEV('N','L',ntarg,ntarg,htarget,evalt,info)
      if(info.ne.0) Stop 'DSYEV failed in Target_new1'
      i=LEN_TRIM(AF_new)+1
      AF=AF_new; write(AF(i:),'(a,i3.3)') '.',klsp

      write(pri,'(/a/)') 'new target energies:'      
      open(nun,file=AF)
      Do i=1,ntarg
       write(nun,'(3F20.8,5x,a)') evalt(i),Etarg(i),evalt(i)-Etarg(i),trim(BFT(i))       
       write(pri,'(3F20.8,5x,a)') evalt(i),Etarg(i),evalt(i)-Etarg(i),trim(BFT(i))       
      End do

      Deallocate(evalt,htarget)

      End Subroutine TARGET_new1


!======================================================================
      Subroutine Target_new2
!======================================================================
!     new target energies 
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer :: i,j,ij,it,jt,ilen,k,info
      Real(8) :: C
      Real(8), allocatable :: htarget(:,:), otarget(:,:), evalt(:)

      ilen = 0
      Do it=1,ntarg; i=Len_trim(BFT(it)); if(i.gt.ilen) ilen=i; End do

      Allocate(htarget(ntarg,ntarg), otarget(ntarg,ntarg), evalt(ntarg))

      htarget=0.d0; otarget=0.d0  
      Do it=1,ntarg; htarget(it,it)=Etarg(it); otarget(it,it)=1.d0; End do

! ... target matrix:

      k = 0
      Do it = 1,ntarg;   Do jt = 1,it
       if(it.eq.jt) Cycle
       ij=(it-1)*it/2+jt
       if(abs(htarg(ij)).gt.eps_tar) then
        k = k + 1
        htarget(it,jt) = htarg(ij); htarget(jt,it) = htarg(ij)
       end if
       if(abs(otarg(ij)).gt.eps_tar) then
        k = k + 1
!        otarget(it,jt) = otarg(ij); otarget(jt,it) = otarg(ij)       ???
       end if
      End do; End do

      if(k.eq.0) Return

      if(debug.gt.0) then
       write(pri,'(/a/)') 'Target hamiltonian:'
       Do it=1,ntarg
        write(pri,'(20f12.5)') htarget(it,1:it)
       End do
       write(pri,'(/a/)') 'Target overlaps:'
       Do it=1,ntarg
        write(pri,'(20f12.5)') otarget(it,1:it)
       End do
      end if

      Call LAP_DSYGV('V','L',ntarg,ntarg,htarget,otarget,evalt,info)
      if(info.ne.0) Stop 'DSYGV failed in Target_new2'
      
      i=LEN_TRIM(AF_new)+1
      AF=AF_new; write(AF(i:),'(a,i3.3)') '.',klsp

      write(pri,'(/a/)') 'new target energies:'      
      open(nun,file=AF)
      Do i=1,ntarg
       write(nun,'(3F20.8,5x,a)') evalt(i),Etarg(i),evalt(i)-Etarg(i),BFT(i)       
       write(pri,'(3F20.8,5x,a)') evalt(i),Etarg(i),evalt(i)-Etarg(i),BFT(i)       
      End do

      Deallocate(evalt,htarget,otarget)

      End Subroutine TARGET_new2

