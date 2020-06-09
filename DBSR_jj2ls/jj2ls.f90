!======================================================================
!     utility     j j 2 l s
!
!                 C O P Y R I G H T -- 2009
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!
!     transform name.c and name.j (name.cm) in jj-coupling to 
!     name_LS.c and name_LS.j in LS-coupling
!
!     Call as:  jj2ls  name (.c, .j)    
!
!     At the moment - only one J total
!----------------------------------------------------------------------
      Use jj2ls
      Use symt_list,    only: JJ_nsymt => nsymt
      Use symt_list_LS, only: LS_nsymt => nsymt
      Use conf_jj,      only: ncfg_jj  => ncfg
      Use conf_LS,      only: ncfg_LS  => ncfg

      Implicit real(8) (A-H,O-Z)
      Character :: AS*80

!----------------------------------------------------------------------
! ... files

      Call GET_COMMAND_ARGUMENT(1,name)

      if(name.eq.'?'.or.command_argument_count().lt.1)  &
      Stop 'Call as: jj2ls name, or name.c, or.name.j'

      i=LEN_TRIM(name)
      if(name(i-1:i).eq.'.c'.or.name(i-1:i).eq.'.j') name(i-1:i)='  '

      AF = trim(name)//'.c'
      Call Check_file(AF)
      Open(nuc,file=AF)
      
      AF = trim(name)//'_LS.log'
      Open(pri,file=AF)

      AF = trim(name)//'.scr'
      Open(nua,file=AF,form='UNFORMATTED')

      AF = trim(name)//'.ovl'
      iovl = 0
      if(Icheck_file(AF).ne.0) then
       open(nua,file=AF,form='UNFORMATTED')
       iovl = 1
      end if

      Call Read_iarg('debug',debug)

!----------------------------------------------------------------------
! ... read jj-configurations:

      Call Read_core_jj(nuc)

      Call Read_conf_jj(nuc,0,'add','check');  JJ_nterms = JJ_nsymt

      write(pri,*) 'ncfg_jj  = ', ncfg_jj
      write(pri,*) 'nterm_jj = ', jj_nterms
      write(pri,*) 

! ... define all possible LS-terms:	    
	  
      Call get_LS_terms;    LS_nterms = LS_nsymt

      write(pri,*) 'nterm_ls = ', ls_nterms
      write(pri,*) 

! ... define the recoupling coefficiens between jj and LS terms: 

      Allocate(C_term(JJ_nterms,LS_nterms));  C_term =  zero

      Call get_tr_coefs

      ncoef = 0 
      Do i=1,jj_nterms; Do j=1,ls_nterms
       if(C_term(i,j).ne.zero) ncoef=ncoef+1
      end do; end do

      write(pri,*) 'get_tr_coef: ncoef =', ncoef,jj_nterms*ls_nterms      
      write(pri,*) 

! ... define the LS configurations:

      Call get_LS_conf

      write(pri,*) 'ncfg_ls  = ', ncfg_ls
      write(pri,*) 

! ... record LS configurations:

      AF = trim(name)//'_LS.c'
      open(nuc,file=AF)
      Call get_LS_core(AS)
      Call write_LS_conf(nuc,AS)
      Call R_term(nuc)

      write(pri,*) 'file ',trim(AF),' is done'
      write(pri,*) 

! ... read transformation matrix from scratch file:

      Deallocate(C_term)
      Allocate (C_trans(ncfg_jj,ncfg_LS)); C_trans = zero
      rewind(nua)
    1 read(nua,end=2) i,j,C; C_trans(i,j)=C; go to 1
    2 Close(nua,status='DELETE') 
      
      write(pri,*) 'C_trans is done'
      write(pri,*) 

! ... check if we have .m file:

      AF = trim(name)//'.m'
      if(Icheck_file(AF).ne.0) then 
       open(num,file=AF,form='UNFORMATTED')
       AF = trim(name)//'.j'
       open(nuj,file=AF)
       Call cm2j(num,nuj)
      end if       
       
! ... transform solutions in j-file:

      AF = trim(name)//'.j'
      Call Check_file(AF)
      Open(nuj,file=AF)

      Call transform_jfile

      End  !  utility-program   j j 2 l s

