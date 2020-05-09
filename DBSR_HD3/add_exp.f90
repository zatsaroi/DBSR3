!======================================================================
      Subroutine Add_exp
!======================================================================
!     introduction of experimental target energies:
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Real(8) :: S
      Integer :: i,it, ich

      write(pri,'(/a,i3/)') 'Experimental energies:  iiexp = ',iiexp
      Do ich = 1,nch; it=iptar(ich); S=E_exp(it)-Etarg(it)
       Do i=ipsol(ich-1)+1,ipsol(ich); a(i,i)=a(i,i)+S; End do
       write(pri,'(2i5,2F16.8,f12.6)') ich,it, &
              E_exp(it),Etarg(it),(E_exp(it)-Etarg(it))*au_eV
      End do

      End Subroutine Add_exp
