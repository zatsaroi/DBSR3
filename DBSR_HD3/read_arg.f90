!======================================================================
      Subroutine Read_arg(nu)
!======================================================================
!     read arguments, first from file unit 'nu, then from comand line
!----------------------------------------------------------------------
      Use dbsr_hd
      Use DBS_nuclear, only: atomic_number, atomic_weight

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,it
      Integer, external :: Icheck_file

      Call Read_ipar(nu,'klsp'  ,klsp   )
      Call Read_ipar(nu,'klsp1' ,klsp1  )
      Call Read_ipar(nu,'klsp2' ,klsp2  )
      Call Read_ipar(nu,'itype' ,itype  )
      Call Read_ipar(nu,'iexp'  ,iexp   )
      Call Read_ipar(nu,'msol'  ,msol   )
      Call Read_rpar(nu,'Emin'  ,Emin   )
      Call Read_rpar(nu,'Emax'  ,Emax   )
      Call Read_ipar(nu,'iwt'   ,iwt    )
      Call Read_rpar(nu,'cwt'   ,cwt    )

      Call Read_iarg('klsp'  ,klsp   )
      Call Read_iarg('klsp1' ,klsp1  )
      Call Read_iarg('klsp2' ,klsp2  )
      Call Read_iarg('itype' ,itype  )
      Call Read_iarg('iexp'  ,iexp   )
      Call Read_iarg('msol'  ,msol   )
      Call Read_rarg('Emin'  ,Emin   )
      Call Read_rarg('Emax'  ,Emax   )
      Call Read_iarg('iwt'   ,iwt    )
      Call Read_rarg('cwt'   ,cwt    )

! ... set the range of partial waves under consideration:

      if(klsp.gt.0) then
       klsp1=klsp; klsp2=klsp
      else
       if(klsp1.le.0) klsp1=1
       if(klsp2.lt.klsp1) klsp2=klsp1
      end if

! ... read experimental energies:

      Call Conv_au (atomic_number,atomic_weight,au_cm,au_eV,0)

      Allocate(E_exp(ntarg));  E_exp = etarg

      if(iexp.gt.0) then
       i=Icheck_file(AF_exp)
       if(i.eq.0) Stop 'iexp >0 but no file for exp. energies'
       Open(nue,file=AF_exp)
       Do i=1,ntarg; read(nue,*) E_exp(i); End do
       Call Read_apar(nue,'unit',unit)
       it=1; Call Read_ipar(nue,'it',it)
       if(it.lt.1.or.it.gt.ntarg) it=1
       if(unit.eq.'cm') E_exp = E_exp / au_cm + Etarg(it)
       if(unit.eq.'eV') E_exp = E_exp / au_eV + Etarg(it)
       Close(nue)
      end if

      Allocate(ip_exp(ntarg))
      Call SORTR(ntarg,E_exp,ip_exp)
      iiexp=0
      Do i=1,ntarg; if(ip_exp(i).ne.i) iiexp=1; End do

      End Subroutine Read_arg

