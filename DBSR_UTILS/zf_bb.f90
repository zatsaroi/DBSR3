!======================================================================
!  zf_bb - utility for the f-value calculations between states in
!  bound.nnn files
!======================================================================
!
!     INPUT:    zf_bb.inp
!
!     OUTPUT:   zf_res from DBSR_DMAT program
!
!     SYSTEM CALL:  dbound_bsw, DBSR_MULT, DBSR_DMAT
!
!     Notes:  1. you would better delete all mult_bnk before running
!             2. results are appended to zf_res
!
!     zf_bb.inp containes:
!
!     atype = E1           or  E2, M1 and so on
!     klsp.
!     list of c-files with energies and exp.coefficients with one J
!
!---------------------------------------------------------------------

      USE channels_jj

      Implicit real(8) (A-H,O-Z)

      Integer :: inp=5; Character(20) :: AF_inp = 'zf_bb.inp'
      Integer :: nut=7; Character(20) :: AF_tar = 'target_jj'

      Character(2) :: atype = 'E1'
      Character(80) :: AS,BI,BJ
      Integer, Allocatable :: ib1(:),ib2(:),ilsp(:) 

! ... read target and channels information:

      Call Check_file(AF_tar)
      Open(nut,file=AF_tar)
      Call Read_target_jj (nut)
      Call Read_channels_jj(nut)
      Close(nut)

! ... input data:

      Call Check_file(AF_inp)
      open(inp,file=AF_inp)

      Call Read_apar(inp,'atype',atype)
      Read(atype,'(1x,i1)') kpol
      Call Read_ipar(inp,'klsp',klsp)
      Allocate(ib1(klsp),ib2(klsp),ilsp(klsp))
      Do i=1,klsp
       read(inp,*) ilsp(i),ib1(i),ib2(i)
      End do

!----------------------------------------------------------------------

      Do ii=1,klsp; i=ilsp(ii); write(BI,'(a,i3.3)') 'bound.',i 
      Do jj=i,klsp; j=ilsp(jj); write(BJ,'(a,i3.3)') 'bound.',j 

       if(atype(1:1).eq.'E'.and.mod(kpol,2).eq.1.and. &
          ipar(i).eq.ipar(j)) Cycle
       if(atype(1:1).eq.'M'.and.mod(kpol,2).eq.0.and. &
          ipar(i).eq.ipar(j)) Cycle
       if(ITRA(jpar(i),kpol+kpol,jpar(j)).eq.0) Cycle
       

        write(AS,'(a,i3.3,a,i3.3,1x,a)') 'dbsr_mult cfg.',i,' cfg.',j,atype 
        Call System(AS)


       Do ib=ib1(ii),ib2(ii)
        write(AS,'(a,i3.3,a,i3.3)') 'dbound_bsw name=a klsp=',i,' state=',ib 
        Call System(AS)
       Do jb=ib1(jj),ib2(jj) 
        write(AS,'(a,i3.3,a,i3.3)') 'dbound_bsw name=b klsp=',j,' state=',jb 
        Call System(AS)

        write(AS,'(a,a)') 'dbsr_dmat a.c b.c c c' 
        Call System(AS)

       End do; End do
      End do; End do 
       
      End ! program zf_cc

