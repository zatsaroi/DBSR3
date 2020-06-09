!======================================================================
      Subroutine transform_jfile
!======================================================================
      Use jj2ls
      Use conf_jj,      only: ncfg_jj  => ncfg, Jtotal
      Use conf_LS,      only: ncfg_LS  => ncfg

      Implicit none 

      Integer :: i,j, jj,k,m, ic,ic1,ic2, is, nc, nsol
      Real(8) :: CM, CN1, CN2, E, gvJ,gvLS
      Character(80) :: AS
      Integer, external :: Ifind_position

      Call Read_ipar(nuj,'ncfg',nc) 
      if(nc.ne.ncfg_jj) then 
       write(*,*) 'nc_j, ncfg_jj =', nc, ncfg_jj
       Stop ' ncfg_jj in j-file inconsistent'
      end if
      Call Read_ipar(nuj,'nsol',nsol) 
      i=Ifind_position(nuj,'Solutions');  read(nuj,*)  

      AF = trim(name)//'_LS.j'
      open(nul,file=AF)
      write (nul,'(8X,A,F5.1,A,I3,A,I7)' ) &
        '  Z = ',one ,'  NEL = ',  1, '   NCFG = ',ncfg_LS
      m = 0

      write (nul,'(//A8,I4,2X,A8,I4)') '  2*J = ',Jtotal,'NUMBER =',nsol

      Allocate(C1(ncfg_jj),C2(ncfg_LS))

      Do is=1,nsol

       C1 = zero
       read(nuj,*) j
       read(nuj,*) E,jj,ic1,ic2
       read(nuj,*) C1(ic1:ic2)
       CN1 = SUM(C1*C1)
       k=0; CM = zero
       Do ic = 1,ncfg_LS
        C2(ic) = SUM(C1(:)*C_trans(:,ic))
        if(abs(C2(ic)).gt.CM) then; k=ic; CM=abs(C2(ic)); end if         
       End do
       Call Get_cfg_LS(k)
       Call Label_c (AS,1,0)
       Call g_factor(Jtotal+1,k,ncfg_LS,C2,gvJ,gvLS) 
       write(nul,'(3x,a,f15.10,4x,a,f15.10,2x,a,f15.10)') &
        'Ssms=',0.d0,'g_J=',gvj,'g_JLS=',gvls
       write(nul,'(i6,f16.8,3x,a)') is,E,trim(AS)
       write(nul,'(7f11.8)') C2
       CN2 = SUM(C2*C2)
       write(pri,'(a,i5,2f12.5)') 'solution:',is,CN1,CN2

      End do ! over is

      write(nul,'(a)') '***'
	  
      End Subroutine transform_jfile 
