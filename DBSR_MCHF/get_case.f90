!======================================================================
      Subroutine get_case
!======================================================================
!     This routine obtains information about the problem to be solved
!     and reads main parameters
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Real(8) :: atw_atom, rms_atom
      Integer :: an, ai
      Character(200) :: A_core, A_conf, AC

! ... name of case:
                                                                                                  
      Call Read_name(name)

      if(len_trim(name).eq.0) then
       Call Read_aarg('atom',atom)
       if(len_trim(atom).gt.0) then
        name=atom
       else 
        Call Read_aarg('ion',ion)
        name=ion
       end if
      end if

      if(LEN_TRIM(name).eq.0.or.name.eq.'?') then
       write(*,*)
       write(*,*) 'DBSR_MCHF is a name-driving procedure.'
       write(*,*) 
       write(*,*) 'Each case (with spesific name) requires four file:'
       write(*,*)
       write(*,*) 'name.c - list of configurations (states)'
       write(*,*)
       write(*,*) 'name.bnk  | int.bnk - angular coefficients (after DBSR_BREIT)'
       write(*,*)
       write(*,*) 'name.knot | knor.dat - B-splines and nuclear parameters'
       write(*,*)
       write(*,*) 'name.inp - running parameters'
       write(*,*)
       write(*,*) 'If files are absent - I will shall create examples to work with.'
       write(*,*)
       write(*,*) 'Call DBSR_MCHF as:    dbsr_mchf  name  parameter=value  ... '
       write(*,*)
       write(*,*) 'with any paramters if their values are different or absent in name.inp'
       write(*,*)
       Stop
      end if

! ... check the input-data file:

      AF_dat = trim(name)//'.inp'
      if(Icheck_file(AF_dat).ne.0) then

       Open(inp,file=AF_dat)

       Call Read_apar(inp,'atom'   ,atom   )
       Call Read_ipar(inp,'an'     ,an     )
       Call Read_apar(inp,'ion'    ,ion    )
       Call Read_ipar(inp,'ai'     ,ai     )
       Call Read_rpar(inp,'z'      ,z      )
       Call Read_rpar(inp,'atw'    ,atw    )
       Call Read_rpar(inp,'rrms'   ,rms    )

       Call Read_apar(inp,'varied' ,avaried)     
       Call Read_rpar(inp,'scf_tol',scf_tol)
       Call Read_rpar(inp,'orb_tol',orb_tol)
       Call Read_rpar(inp,'end_tol',end_tol)
       Call Read_ipar(inp,'max_it' ,max_it )

       Call Read_ipar(inp,'method' ,method )
       Call Read_ipar(inp,'all'    ,all    )
       Call Read_ipar(inp,'irhs'   ,irhs   )
       Call Read_ipar(inp,'newton' ,newton )
       Call Read_ipar(inp,'rotate' ,rotate )
       Call Read_ipar(inp,'debug'  ,debug  )
       Call Read_ipar(inp,'n_corr' ,n_corr )
                                           
       Call Read_ipar(inp,'ipzero' ,ipzero )
       Call Read_ipar(inp,'iqzero' ,iqzero )
       Call Read_ipar(inp,'jpzero' ,jpzero )
       Call Read_ipar(inp,'jqzero' ,jqzero )

!       Call Read_rpar(inp,'c_au'  ,c_au  )

      end if

! ... check the command-line arguments:

      Call Read_aarg('atom'   ,atom   )
      Call Read_iarg('an'     ,ai     )
      Call Read_aarg('ion'    ,ion    )
      Call Read_iarg('ai'     ,ai     )
      Call Read_rarg('z'      ,z      )
      Call Read_rarg('atw'    ,atw    )
      Call Read_rarg('rrms'   ,rms    )

      Call Read_aarg('varied' ,avaried)     

      Call Read_rarg('scf_tol',scf_tol)
      Call Read_rarg('orb_tol',orb_tol)
      Call Read_rarg('end_tol',end_tol)
      Call Read_iarg('max_it' ,max_it )

      Call Read_iarg('method' ,method )
      Call Read_iarg('all'    ,all    )
      Call Read_iarg('irhs'   ,irhs   )
      Call Read_iarg('newton' ,newton )
      Call Read_iarg('rotate' ,rotate )
      Call Read_iarg('debug'  ,debug  )
      Call Read_iarg('n_corr' ,n_corr )
	
      Call Read_iarg('ipzero' ,ipzero )
      Call Read_iarg('iqzero' ,iqzero )
      Call Read_iarg('jpzero' ,jpzero )
      Call Read_iarg('jqzero' ,jqzero )

!      Call Read_iarg('c_au'   ,c_au   )

!------------------------------------------------------------------------------
! ... check the atom if available:

      if(len_trim(atom).eq.0.and.an.ne.0) &
       Call Def_atom(an,atom, atw_atom, rms_atom, A_core,A_conf)
      if(len_trim(ion).eq.0.and.ai.ne.0) &
       Call Def_atom(ai,ion , atw_atom, rms_atom, A_core,A_conf)

      if(len_trim(atom).ne.0.and.len_trim(ion).eq.0) ion=atom 
      if(len_trim(atom).eq.0.and.len_trim(ion).ne.0) atom=ion 

      if(len_trim(atom).gt.0) then
       Call Def_atom(an,atom, atw_atom, rms_atom, A_core,A_conf)
       z = an
       if(atw.eq.0.d0) atw=atw_atom
       if(rms.eq.0.d0) rms=rms_atom
      end if

      if(len_trim(ion).gt.0) then
       Call Def_atom(ai,ion , atw_atom, rms_atom, A_core,A_conf)
      end if
      if(len_trim(ion) .eq.0) ion=atom

! ... prepare c-file if absent   (require additional PROGRAMS ???)

      AF_cfg = TRIM(name)//'.c'
      Call Read_aarg('c',AF_cfg)
      if(Icheck_file(AF_cfg).eq.0) then        
       write(AC,'(10a)') 'genjconf atom=',trim(ion),' out=',trim(name),'.conf'
       Call System(AC)
       write(AC,'(10a)') 'genjterm ',trim(name)
       Call System(AC)
      end if

! ... prepare bnk-file if absent   (require additional PROGRAMS ???)

      AF_bnk = TRIM(name)//'.bnk'
      Call Read_aarg('bnk',AF_cfg)
      if(Icheck_file(AF_bnk).eq.0) then        
       write(AC,'(10a)') 'dbsr_breit3 ',trim(name),'.c'
       Call System(AC)
      end if

!--------------------------------------------------------------------------------------

      AF_log = TRIM(name)//'.log'
      Open(log,file=AF_log)

      write(log,'(80("=")/T20,a/80("="))') &
         'MULTICONFIGURATION B-SPLINE DIRAC-HARTREE-FOCK'
 
      write(scr,'(79("=")/T15,a/79("="))') &
         'MULTICONFIGURATION B-SPLINE DIRAC-HARTREE-FOCK'

      write(scr,'(/a,T20,a/)') 'Name of case:',trim(name)
      write(log,'(/a,T35,a )') 'Name of case:',trim(name)

      if(Icheck_file(AF_dat).eq.0) Call Write_inp

      End Subroutine get_case


!======================================================================
      Subroutine Write_inp
!======================================================================
!     This routine prepare default input file 
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Integer :: i,j,j1,j2

      open(inp,file=AF_dat)
      rewind(inp)

      write(inp,'(a/)') 'Main parameters:'

      write(inp,'(a,1x,a,T40,a)') 'atom    = ',trim(atom),'- atomic symbol'
      if(ion.ne.atom) &
      write(inp,'(a,1x,a,T40,a)') 'ion     = ',trim(ion) ,'- ionic stage'
      write(inp,'(a,f4.1,T40,a)') 'z       = ',z,         '- nuclear number'

      write(inp,'(/a,a)')         'varied  =  ',trim(adjustl(avaried))

      write(inp,'(/a,i2,T40,a)')  'eol     = ',eol, '- optimzation mode'

      if(nlevels.gt.0) then
       write(inp,*)
       Select case(eol)
        Case(1,5)
         j1 =1 
         Do i=1,njbl
          j2=j1+nlevelb(i)-1
          write(inp,'(a,i3,50(a1,i3))') 'levels  =',i,(',',level(j),j=j1,j2)      
          j1 = j2 + 1
         end do
        Case(9)
         Do i=1,nlevels
          write(inp,'(a,i3,a,i3,a,f10.6)') &
                    'weights  =', block(i),',',level(i),',',weight(i)
         end do
       End Select         
      end if

      if(len_trim(physical).gt.0) &
      write(inp,'(/a,a)') 'physical=  ',trim(physical)

      Call Write_run_par(inp)

      write(inp,'(/a,a)')  'All parameters from input files ', &
                            'can be replaced from command line as:'
      write(inp,'(/10x,a/)') 'dbsr_mchf [name] par1=value par2=value par3=value ... '

      write(inp,'(80("-"))')                     
      write(inp,'(/a/)') 'Name-driven fine-name and key-words for their re-definition:' 

      write(inp,'(a,T20,a,T40,a)') 'name.inp',  ' ',       '- input parameters'
      write(inp,'(a,T20,a,T40,a)') 'name.log',  ' ',       '- run output and summry'
      write(inp,'(a,T20,a,T40,a)') 'name.c',    'c=...',   '- input configurations'
      write(inp,'(a,T20,a,T40,a)') 'name.bnk',  'bnk=...', '- angular coefficients'
      write(inp,'(a,T20,a,T40,a)') ' ',         'int_bnk', '- generic name'
      write(inp,'(a,T20,a,T40,a)') 'name.knot', 'knot=...','- B-spline parameters'
      write(inp,'(a,T20,a,T40,a)') ' ',         'knot.dat','- generic name'
      write(inp,'(a,T20,a,T40,a)') 'name.bsw',  'inp=...' ,'- input w.f. if any'
      write(inp,'(a,T20,a,T40,a)') 'name.bsw',  'out=...', '- output w.f.'
      write(inp,'(a,T20,a,T40,a)') 'name.j',    'j=...',   '- expansion coef.s in the DBSR format' 

      write(inp,'(/80("-")/)')                     
      write(inp,'(a)') ' Additional information for input parameters:'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' varied   - possible options -> all, none, list of nl, =last, n=..., n>...'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' eol      - indicates the mode for the state weights:'
      write(inp,'(a)') '            =1 - equally weighted' 
      write(inp,'(a)') '            =5 - statistically weighed, default' 
      write(inp,'(a)') '            =9 - defined by user in .conf or .c files'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' relevant level-block parameters:'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' if eol = 1 or 5, repeat for each J-block: '
      write(inp,'(a)') '            '
      write(inp,'(a)') ' levels = 1,1,2 - block, list of levels to be optimized' 
      write(inp,'(a)') '            '
      write(inp,'(a)') ' if eol = 9, repeat for each level to be optimized'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' weights = 1,2,0.5 - block, level, weight' 
      write(inp,'(a)') '            '
      write(inp,'(a)') ' if level information is absent - program will optimized the first '
      write(inp,'(a)') ' level in each block'          
      write(inp,'(a)') '            '
      write(inp,'(a)') ' ipzero = 0 means  l+1 zero B-splines in the beginning for large component'
      write(inp,'(a)') ' iqzero = 0 means  l (l+2) zero B-splines in the beginning for small component'

      write(inp,'(80("-"))')

      End Subroutine write_inp


!======================================================================
      Subroutine Write_run_par(nu)
!======================================================================
!     print running parameters 
!----------------------------------------------------------------------
      Use dbsr_mchf
      Integer, intent(in) :: nu

      write(nu,'(/80("-")/)')
      write(nu,'(a)') 'Running parameters:'
      write(nu,*)

      write(nu,'(a,1PE9.2,T40,a)') 'scf_tol = ',scf_tol, &
                '- tolerance for energy convergence'
      write(nu,'(a,1PE9.2,T40,a)') 'orb_tol = ',orb_tol, &
                '- tolerance for orbital convergence'
      write(nu,'(a,1PE9.2,T40,a)') 'end_tol = ',end_tol, &
                '- tolerance for ending zeros'
      write(nu,'(a,i3,T40,a)') 'max_it  = ',max_it, &
                '- max. number of iterations'
      write(nu,'(a)') 
      write(nu,'(a,i2,T40,a)') 'ipzero  = ',ipzero, &
                '- initial zeros for larger component' 
      write(nu,'(a,i2,T40,a)') 'iqzero  = ',iqzero, &
                '- initial zeros for small component' 
      write(nu,'(a,i2,T40,a)') 'jpzero  = ',jpzero, &
                '- final zeros for larger component'
      write(nu,'(a,i2,T40,a)') 'jqzero  = ',jqzero, &
                '- final zeros for small component'

      write(nu,'(/80("-")/)')
      write(nu,'(a)') 'Additonal options (not applied if = 0)'
      write(nu,*)
 
      write(nu,'(a,i2,T40,a)') 'method  = ',method,  &
                '- method for solving MCHF equation'
      write(nu,'(a,i2,T40,a)') 'all     = ',all,    &
                '- collective optimization'  
      write(nu,'(a,i2,T40,a)') 'irhs    = ',irhs,   &
                '- convert right-hand-side to main matrix'  
      write(nu,'(a,i2,T40,a)') 'newton  = ',newton, &
                '- use Newton-Rapson method'  
      write(nu,'(a,i2,T40,a)') 'rotate  = ',rotate, &
                '- use rotations'  
      write(nu,'(a,i2,T40,a)') 'debug   = ',debug,  &
                '- additional debug output'

      write(nu,'(/80("-")/)')
      write(nu,'(a,f18.13,T40,a)') 'c_au    = ',c_au,'- speed of light'
      write(nu,'(/80("-") )')

      End Subroutine Write_run_par
