!======================================================================
      Subroutine Read_arg(nu)
!======================================================================
!     read input parameters, from unit "nu" or from command line
!----------------------------------------------------------------------
      Use dbsr_mat
      Use DBS_dhl_pq
      Use c_data, only: mem_cdata

      Implicit none 
      Integer, intent(in) :: nu
      Integer :: i,itype
      Integer, external :: Icheck_file

! ... read parameters from file if any

      Call Read_ipar(nu,'klsp1'  ,klsp1  )
      Call Read_ipar(nu,'klsp2'  ,klsp2  )
      Call Read_ipar(nu,'mk'     ,mk     )
      Call Read_ipar(nu,'mblock' ,mblock )
      Call Read_ipar(nu,'nblock' ,nblock )
      Call Read_ipar(nu,'kblock' ,kblock )
      Call Read_ipar(nu,'mcbuf'  ,mcbuf  )
      Call Read_rpar(nu,'eps_c'  ,eps_c  )
      Call Read_rpar(nu,'eps_v'  ,eps_v  )
      Call Read_rpar(nu,'eps_det',eps_det)
      Call Read_ipar(nu,'iitar'  ,iitar  )
      Call Read_ipar(nu,'mbreit' ,mbreit )
      Call Read_ipar(nu,'mbloch' ,mbloch )
      Call Read_ipar(nu,'pri_acf',pri_acf)
      Call Read_rpar(nu,'Ecore'  ,EC     )
      Call Read_rpar(nu,'RB'     ,RB     )
      Call Read_rpar(nu,'pnu'    ,pnu    )
      Call Read_rpar(nu,'s_ovl'  ,s_ovl  )

      Call Read_ipar(nu,'check_target',check_target)

      Call Read_ipar(nu,'debug'  ,debug  )

      Call Read_rpar(nu,'Edmin' ,Edmin )
      Call Read_rpar(nu,'Edmax' ,Edmax )
      Call Read_rpar(nu,'Egap'  ,Egap  )

! ... overwrite parameters from arguments if any

      Call Read_iarg('klsp1'  ,klsp1  )
      Call Read_iarg('klsp2'  ,klsp2  )
      Call Read_iarg('mk'     ,mk     )
      Call Read_iarg('mblock' ,mblock )
      Call Read_iarg('nblock' ,nblock )
      Call Read_iarg('kblock' ,kblock )
      Call Read_iarg('mcbuf'  ,mcbuf  )
      Call Read_rarg('eps_c'  ,eps_c  )
      Call Read_rarg('eps_v'  ,eps_v  )
      Call Read_rarg('eps_det',eps_det)
      Call Read_iarg('iitar'  ,iitar  )
      Call Read_iarg('mbreit' ,mbreit )
      Call Read_iarg('mbloch' ,mbloch )
      Call Read_iarg('pri_acf',pri_acf)
      Call Read_rarg('Ecore'  ,EC     )
      Call Read_rarg('RB'     ,RB     )
      Call Read_rarg('pnu'    ,pnu    )
      Call Read_rarg('s_ovl'  ,s_ovl  )

      Call Read_iarg('check_target',check_target)

      Call Read_iarg('debug'  ,debug  )

      Call Read_rarg('Edmin'  ,Edmin  )
      Call Read_rarg('Edmax'  ,Edmax  )
      Call Read_rarg('Egap'   ,Egap   )

      klsp=0
      Call Read_ipar(nu,'klsp',klsp)
      Call Read_iarg(   'klsp',klsp)
      if(klsp.ne.0) then; klsp1=klsp; klsp2=klsp; end if
      if(klsp2.lt.klsp1) klsp2=klsp1

      ! ... additional check for removed B-splines:

      itype = 0
      Call Read_ipar(nu,'itype',itype)
      Call Read_iarg(   'itype',itype)
      if(itype.eq.-1) then
       ilzero=-1; jlzero=-1; ibzero=1; jbzero=1
      else
       ilzero=-1; jlzero=-1; ibzero=0; jbzero=0
      end if

      Call Read_iarg('ilzero' ,ilzero )
      Call Read_iarg('jlzero' ,jlzero )
      Call Read_iarg('ibzero' ,ibzero )
      Call Read_iarg('jbzero' ,jbzero )

      Call Read_ipar(nu,'ilzero' ,ilzero )
      Call Read_ipar(nu,'jlzero' ,jlzero )
      Call Read_ipar(nu,'ibzero' ,ibzero )
      Call Read_ipar(nu,'jbzero' ,jbzero )

! ... print global parameters:

      if(prj.gt.0) then
       write(prj,'(a)') 'DBSR_MAT parameters:'

       write(prj,'(/a/)') 'Partial waves:' 
       write(prj,'(a,i10,T20,a)') 'klsp1  =',klsp1,'- initial index'
       write(prj,'(a,i10,T20,a)') 'klsp2  =',klsp2,'- final index'

       if(mbreit.gt.0) & 
       write(prj,'(/a,i10,T20,a)') 'mbreit =',mbreit,'- Breit interaction is included'
       if(mbreit.eq.0) & 
       write(prj,'(/a,i10,T20,a)') 'mbreit =',mbreit,'- Breit interaction is exluded'

       if(iitar.eq.0) & 
       write(prj,'(/a,i10,T20,a)') 'iitar  =',iitar,&
                                   '- target states are supposed to be orthogonal'
       if(iitar.eq.1) & 
       write(prj,'(/a,i10,T20,a)') 'iitar  =',iitar,&
                                   '- target states may be non-orthogonal'
       if(iitar.eq.2) & 
       write(prj,'(/a,i10,T20,a)') 'iitar  =',iitar,&
                                   '- target states may not diagonalized the Hamiltonian'

       write(prj,'(/a,F10.5,T20,a)') 's_ovl  =',s_ovl,'- tollerence for channel overlaps'

       write(prj,'(/a/)') 'Bloch operator:' 
       write(prj,'(a,i10  ,T20,a)') 'mbloch =',mbloch,'- flag for Bloch operator'
       write(prj,'(a,F10.5,T20,a)') 'RB     =',RB,    '- b-parameter (derivative at r=a)'
       write(prj,'(a,F10.5,T20,a)') 'pnu    =',pnu,   '- "nu"-parameter (rotation)'

       write(prj,'(/a/)') 'Boundary conditions:' 
       write(prj,'(a,i10,T20,a)')  'ilzero =',ilzero,&
                                   '- zero B-splines for large component at r=0'
       if(jlzero.lt.0) jlzero=ilzero
       write(prj,'(a,i10,T20,a)')  'jlzero =',jlzero,&
                                   '- zero B-splines for small component at r=0'
!       if(mbloch.gt.0) ibzero = 0
       write(prj,'(a,i10,T20,a)')  'ibzero =',ibzero,&
                                   '- zero B-splines for large component at r=a'
!       if(mbloch.ge.0) jbzero = 0
       if(jbzero.lt.0) jbzero=ibzero
       write(prj,'(a,i10,T20,a)')  'jbzero =',jbzero,&
                                   '- zero B-splines for small component at r=a'

       write(prj,*)
       write(prj,'(a)') 'Restrictions on channel eigenvalues:'
       write(prj,*)
       write(prj,'(a,1Pd10.1,T20,a)') 'Edmin  =',Edmin,'- minimum energy, -2*c^2'
       write(prj,'(a,1Pd10.1,T20,a)') 'Edmax  =',Edmax,'- maximum energy'
       write(prj,'(a,1Pd10.1,T20,a)') 'Egap   =',Egap, '- gap around E = 0.d0'
       write(prj,'(/a/)') 'All one-channel solutions with &
         &E<Edmin, E>Edmax or abs(E)<Egap will be disregaded.'

       i = ntype_S*(mk+2)*2 + 1
       if(nblock.lt.i) nblock=i
       write(prj,'(/a/)') 'Module cdata parameters:' 
       write(prj,'(a,i10,T20,a)') 'mk     =',mk,'- max. multipol index'
       write(prj,'(a,i10,T20,a)') 'mblock =',mblock,'- size of block'
       write(prj,'(a,i10,T20,a)') 'nblock =',nblock,'- number of blocks' 
       write(prj,'(a,i10,T20,a)') 'kblock =',nblock,'- number of blocks for one type' 
       if(mbreit.eq.0) Call alloc_c_data(ntype_R,0,mk,  mblock,nblock,kblock,eps_C)
       if(mbreit.gt.0) Call alloc_c_data(ntype_S,0,mk+1,mblock,nblock,kblock,eps_C)
       write(prj,'(a,T40,f8.1,a)') 'memory of cdata module', mem_cdata,' Mb'

       write(prj,'(/a/)') 'Buffer parameters:' 
       mem_buffer = 24.d0 * mcbuf / (1024*1024)
       write(prj,'(a,i10,T20,a)') 'mcbuf  =',mcbuf,'- buffer dimension'
       write(prj,'(a,T40,f8.1,a)') 'memory of the buffer:', mem_buffer,' Mb'

       write(prj,'(/a/)') 'Cut-off parameters:' 
       write(prj,'(a,1Pd10.1,T20,a)') 'eps_c  =',eps_c,  '- tollerance for coefficients'
       write(prj,'(a,1Pd10.1,T20,a)') 'eps_det=',eps_det,'- tollerance for det.overlaps'
       write(prj,'(a,1Pd10.1,T20,a)') 'eps_v  =',eps_v,  '- tollerance for <nl|kl>'
       write(prj,*) 

       write(prj,*)
       write(prj,'(a)') 'Additional debug output:'
       write(prj,*)
       write(prj,'(a,i10  ,T20,a)') 'debug  =',debug,  '- flag for debug output'
       write(prj,'(a,i10  ,T20,a)') 'pri_acf=',pri_acf,&
                                    '- asymptotic coefficients output'

      end if

      End Subroutine Read_arg


