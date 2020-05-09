!======================================================================
      Subroutine R_arg(nu)
!======================================================================
!     read arguments, first from file unit 'nu, then from comand line
!----------------------------------------------------------------------
      Use dbsr_pol
          
      Implicit None
      Integer, intent(in) :: nu
      integer :: n

      if(nu.ne.0) then
       Call Read_ipar(nu,'klsp'  ,klsp  )
       Call Read_ipar(nu,'nortb' ,nortb )
       if(nortb.gt.0) then
        Allocate(iortb(nortb)); iortb=0
        Call Read_iarray(nu,'iortb',nortb,iortb)
       end if
      end if

      Call Read_iarg('klsp',klsp)
      write(ALSP,'(i3.3)') klsp

      n=0
      Call Read_iarg('nortb',n)
      if(n.gt.0) then
       nortb=n
       if(Allocated(iortb)) Deallocate(iortb)
       Allocate(iortb(nortb)); iortb=0
       Call Read_iarr('iortb',nortb,iortb)
      end if

      End Subroutine R_arg

