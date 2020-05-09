!======================================================================
      Subroutine Read_arg
!----------------------------------------------------------------------
! ... read arguments from the command line and open some main files 
!----------------------------------------------------------------------
      Use dbsr_dmat

! ... read first four basic "position" parameters:

      iarg = IARGC()
      
      if(iarg.lt.4) then
       write(*,*) 'DBSR_DMAT: number of arguments  <  4'
       write(*,*)
       write(*,*) '1 - file name for initial-state c-file'
       write(*,*) '2 - file name for final-state c-file'
       write(*,*) '3 - structure mode for initial state: c,j,b'
       write(*,*) '4 - structure mode for final state: c,j,b,p'
       write(*,*)
       write(*,*) 'example:  dbsr_dmat 1.c cfg.002 j b F|G dd'
       Stop
      end if

      Call GETARG(1,name1)
      Call Check_file(name1);  Open(nuc1,file=name1) 

      Call GETARG(2,name2)
      Call Check_file(name2);  Open(nuc2,file=name2) 

      Call GETARG(3,ctype1)
      i1 = LEN_TRIM(name1)
      if(name1(i1:i1).eq.'c') then; ilsp1=0
      else; ALS1=name1(i1-2:i1); read(ALS1,'(i3)') ilsp1; end if
      if(ilsp1.eq.0.and.ctype1.eq.'b') Stop 'ctype1 inconsistent'

      Call GETARG(4,ctype2)
      i2 = LEN_TRIM(name2)
      if(name2(i2:i2).eq.'c') then; ilsp2=0
      else; ALS2=name2(i2-2:i2); read(ALS2,'(i3)') ilsp2; end if
      if(ilsp2.eq.0.and.ctype2.eq.'b') Stop 'ctype2 inconsistent'
      if(ilsp2.eq.0.and.ctype2.eq.'p') Stop 'ctype2 inconsistent'
      if(ilsp2.eq.0.and.ctype2.eq.'q') Stop 'ctype2 inconsistent'

      Open(pri,file=AF_log)
      write(pri,'(a/a/)') 'D B S R _ D M A T',&
                          '*****************'
! ... read additional "key-word" parameters from the command line if any

      Call Read_iarg('mstate1' ,mstate1)
      Call Read_iarg('mstate2' ,mstate2)
      Call Read_iarg('istate1' ,istate1)
      Call Read_iarg('istate2' ,istate2)
      Call Read_iarg('ialpha'  ,ialpha )
      Call Read_aarg('gf'      ,gf     )
      Call Read_aarg('dd'      ,dd     )
      Call Read_iarg('debug'   ,debug  )

      Call Read_aarg('atype'   ,atype  )
      Read(atype,'(1x,i1)') kpol; ktype = atype(1:1)
      AF_bnk = 'mult_bnk_'//atype
      Call Check_file(AF_bnk);  Open(nub,file=AF_bnk,form='UNFORMATTED')

      End Subroutine Read_arg

