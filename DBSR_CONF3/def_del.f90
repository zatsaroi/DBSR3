!======================================================================
      Subroutine def_del
!======================================================================
!     define list of channels to be deleted if any
!----------------------------------------------------------------------
      Use dbsr_conf, only: ndel, dlsp, dkch, dtar, AF_del

      Implicit none
      Integer :: i,j,nlsp,nch,kch,iptar,ipconf,idel,nud
      Character(80) :: line
      Integer, external :: Icheck_file, Ifind_position

      ndel = 0
      if(Icheck_file(AF_del).eq.0) Return
      Call Find_free_unit(nud); Open(nud,file=AF_del)

      nlsp = 0
      Call Read_ipar(nud,'nlsp',nlsp); if(nlsp.le.0) Return
      i = Ifind_position(nud,'channels:'); read(nud,*)

      Do i = 1,nlsp
       read(nud,*)
       read(nud,'(a)') line
       j = INDEX(line,'nch =') + 5
       read(line(j:),*) nch
       read(nud,*)
       Do j = 1,nch
        read(nud,'(a)') line
        read(line(12:),*) kch,iptar,ipconf
        if(ipconf.le.0) ndel=ndel+1
       End do
      End do
      if(ndel.eq.0) Return

      Allocate(dlsp(ndel),dkch(ndel),dtar(ndel))

      i = Ifind_position(nud,'channels:'); read(nud,*)
      idel=0
      Do i = 1,nlsp
       read(nud,*)
       read(nud,'(a)') line
       j = INDEX(line,'nch =') + 5
       read(line(j:),*) nch
       read(nud,*)
       Do j = 1,nch
        read(nud,'(a)') line
        read(line(12:),*) kch,iptar,ipconf
        if(ipconf.gt.0) Cycle
        idel=idel+1
        dlsp(idel)=i
        dkch(idel)=kch
        dtar(idel)=iptar
       End do
      End do

      Close(nud)
      
      End Subroutine def_del 


!======================================================================
      Integer Function ich_del(klsp,kapa,itar)
!======================================================================
!     check if given channel in the "delete" list
!----------------------------------------------------------------------
      Use dbsr_conf, only: ndel, dlsp, dkch, dtar

      Implicit none
      Integer :: klsp,kapa,itar, i

      ich_del = 0;   if(ndel.eq.0) Return
      Do i = 1,ndel
       if(klsp.ne.dlsp(i)) Cycle
       if(kapa.ne.dkch(i)) Cycle
       if(itar.ne.dtar(i)) Cycle
       ich_del = 1
       Exit
      End do

      End Function ich_del
