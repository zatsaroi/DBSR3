!======================================================================
      Subroutine Write_target_jj_example(nut)
!======================================================================
!     write target information to file 'nut'  
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nut

      rewind(nut)
      write(nut,'(a)') 'title:' 
      write(nut,'(80(''-''))')
      write(nut,'(a)') 'nz    = ...       !   nuclear charge' 
      write(nut,'(a)') 'nelc  = ...       !   number of electrons'
      write(nut,'(80(''-''))')
      write(nut,'(a)') 'ntarg = ...       !   number of target states'
      write(nut,'(80(''-''))')
      write(nut,'(a)') 'target_1'
      write(nut,'(a)') 'target_2'
      write(nut,'(a)') '...'
      write(nut,'(80(''-''))')
      write(nut,'(a)') 'nlsp  = ...       !   number of partial waves'
      write(nut,'(80(''-''))')
      write(nut,'(a)') 'partial wave 1  (as #.  2J  parity)'
      write(nut,'(a)') 'partial wave 2  (as #.  2J  parity)'
      write(nut,'(a)') '...'
      write(nut,'(80(''-''))')
      write(nut,'(a)') 'kpert = ...       !   number of optional pertubers'
      write(nut,'(80(''-''))')
      write(nut,'(a)') 'pertuber1  (as # of partial wave and name)'
      write(nut,'(a)') 'pertuber2  (as # of partial wave and name)'
      write(nut,'(a)') '...'
      write(nut,'(80(''-''))')

      End Subroutine Write_target_jj_example


