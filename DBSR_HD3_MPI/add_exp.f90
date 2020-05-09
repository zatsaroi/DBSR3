!======================================================================
      Subroutine Add_exp
!======================================================================
!     introduction of experimental energies:
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Real(8) :: S
      Integer :: i,ii, it, idim, ich

      if(io_processor) & 
       write(pri,'(/a/)') 'Experimental energies:'  

      add = zero

      Do ich = 1,nch
       idim=ipsol(ich)-ipsol(ich-1)
       if(io_processor) then
        it=iptar(ich)
        S=E_exp(it)-Etarg(it)
        Do i=1,idim; add(i,i)=S; End do
        write(pri,'(2i5,2F16.8,f10.3)') ich,it, &
              E_exp(it),Etarg(it),S*27.2112
       end if
       ii=ipsol(ich-1)+1

       Call pdgeadd ('notrans', idim, idim,     &
                      one,  add, 1, 1, descadd, &
                      one,    a, ii, ii, desca  )
      End do

      End Subroutine Add_exp


