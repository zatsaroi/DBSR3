!======================================================================
      Subroutine Record_orth
!======================================================================
! ... record orthogonal conditions for sct. orbitals:
!----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none
      Integer :: ich,io, ii,jj, it,irec
      Integer, external :: Iort

      write(AF,'(a,i3.3)')'cfg.',ilsp
      Open (nuc,file=AF,position='APPEND')

      write(nuc,'(/a/)') 'Derived orth. conditions:'

! ... record the orth. conditions with indication of main compensation
! ... configuration (see pert_comp.nnn for total coefficients)      

      Do ich=1,nch;  ii = ipch(ich); it = iptar(ich)
       irec = 0 
       Do io = 1,nphys_sub; jj = jp_sub(io)
        if(KEF(ii).ne.KEF(jj)) Cycle
        if(Iort(ii,jj).eq.0) then
         write(nuc,'(a1,a5,a1,a5,a3,3x,a12,5x,a)') &
              '<',ELF(ii),'|',ELF(jj),'>=0',trim(BFT(it)),trim(AFT(it))
         irec = 1
        end if
       End do
       if(irec.gt.0) write(nuc,*)
      End do
      Close(nuc)

      End Subroutine Record_orth
