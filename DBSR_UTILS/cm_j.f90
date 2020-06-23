!======================================================================
!     utility       c m _ j
!
!                   C O P Y R I G H T -- 2006
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!     Convert cm- or m-file after rscf2 or rci2 to j-file of DBSR format
!
!     Call as:   cm_j  name.cm  or  cm_j name.m
!----------------------------------------------------------------------
      Use conf_jj

      Implicit real(8) (A-H,O-Z)

      Character(6) :: G92MIX
      Character(80) :: AN, AF ,BF

      Real(8), allocatable :: EVAL(:), EVEC(:,:)
      Integer, allocatable :: IP(:),JP(:)
      Character(300), allocatable :: AC1(:),AC2(:),AC3(:)

      Integer, allocatable :: IVEC(:),IATJPO(:),IASPAR(:), ICCMIN(:)

      Integer :: nuc =1       !   name.c         
      Integer :: num =2       !   name.mix or name.cm
      Integer :: nuj =3       !   name.j

!----------------------------------------------------------------------
!                                                           input data:
      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AN)
      
      if(iarg.lt.1.or.AN.eq.'?') then
       write(*,*)
       write(*,*) 'cm_j converts cm-file after rci2(GRASP) or m-file after rscf2(GRASP)' 
       write(*,*) 'to j-file of DBSR format'                                                       
       write(*,*)
       write(*,*) 'Call as:  cm_j  name.cm   or   cm_j name.m'                                                       
       write(*,*)
       write(*,*) '(also required is name.c file with the configuration list)'
       write(*,*)
       STOP
      end if

!----------------------------------------------------------------------
! ... read list of configuration:

      ii=INDEX(AN,'.',BACK=.TRUE.)-1; if(ii.eq.0) ii=Len_trim(AN)     

      AF=AN(1:ii)//'.c'
      Call Check_file(AF)
      Open(nuc,file=AF,STATUS='OLD')

      Call Read_conf_jj(nuc,0,'a','b') 
      Call R_label_jj(nuc,0)

      AF=AN(1:ii)//'.j';   Open(nuj,file=AF)
      write(nuj,'(a,i8)') 'ncfg =',ncfg

!----------------------------------------------------------------------
! ... read solution:

      Call Check_file(AN)
      Open(num,file=AN,form='UNFORMATTED')
      read(num) G92MIX
      if(G92MIX.ne.'G92MIX') then
       write(*,'(a,a,a)') 'file ',AN,' is not GRASP2K mix-file'
       Stop ' '
      end if

       read(num) NELEC, NCFTOT, NW, NVECTOT, NVECSIZ, NBLOCK
       if(NCFTOT.ne.ncfg) Stop 'ncfg in c-fail <> NCFTOT in cm-fail'

       write(nuj,*)
       write(nuj,'(a,i8)') 'nsol =',NVECTOT

       write(nuj,*)
       write(nuj,'(a)') 'Solutions:'

       ic2=0
       Do ib = 1,NBLOCK

        read(num) NB, NCFBLK, NEVBLK, IATJP, IASPA
        ic1=ic2+1; ic2=ic2+NCFBLK

        if(allocated(ICCMIN)) Deallocate(ICCMIN); 
        Allocate(ICCMIN(NEVBLK))
        read(num) (ICCMIN(I),I=1,NEVBLK)

        if(allocated(EVAL)) Deallocate(EVAL); 
        Allocate(EVAL(NEVBLK))
        read(num) EAV, (EVAL(i), I = 1, NEVBLK)

        if(allocated(EVEC)) Deallocate(EVEC); 
        Allocate(EVEC(NCFBLK,NEVBLK))
        read(num) ((EVEC(i,j), i=1,NCFBLK), J=1,NEVBLK)

        Do k=1,NEVBLK
         nsol = nsol+1
         CM = 0.d0; kk = 1
         Do i=1,NCFBLK
          if(abs(EVEC(i,k)).lt.CM) Cycle
          kk = i; CM = abs(EVEC(i,k))
         End do
         ii = ic1-1+kk

         write(nuj,'(i8,2x,a)') nsol,trim(LABEL(ii))
         write(nuj,'(f16.8,3i8)') EVAL(k)+EAV, IATJP-1, ic1,ic2       
         write(nuj,'(6f20.15)') (EVEC(i,k),i=1,NCFBLK)
        End do

       End do   ! over blocks

      End ! program cm_j

