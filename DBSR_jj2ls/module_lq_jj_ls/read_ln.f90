!---------------------------------------------------------------------
      Module ln
!---------------------------------------------------------------------
      Implicit none

      Integer :: ncase = 0
      Integer, allocatable, dimension(:) :: inj1,inj2,ikt1,ikt2,iJT, &
	                                       ikt,iIC,iJC
      End module ln


!----------------------------------------------------------------------
      Subroutine read_ln(l,n)      
!----------------------------------------------------------------------
      Use ln

      Character :: AF*10
      Character, external :: AL

      nus=1; write(AF,'(a,i2)') AL(l,1),n; Call Clean_a(AF)
      Call Check_file(AF)
      open(nus,file=AF)
 
      Call read_ipar(nus,'nn',ncase)
      if(ncase.eq.0) Stop ' ncase = 0'

      Allocate(inj1(ncase),inj2(ncase),ikt1(ncase),ikt2(ncase), &
	           iJT(ncase),ikt(ncase),iIC(ncase),iJC(ncase))
      rewind(nus)
      Do i = 1,ncase
       read(nus,*) inj1(i),inj2(i),ikt1(i),ikt2(i),iJT(i), &
                   ikt(i),iIC(i),iJC(i)
      End do
      Close(nus)

      End  Subroutine read_ln

