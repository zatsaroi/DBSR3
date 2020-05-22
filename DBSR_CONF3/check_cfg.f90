!======================================================================
      Subroutine Check_cfg
!======================================================================
! ... for given compensation configuration (ic) record all relevant 
! ... information, see add_comp below
!----------------------------------------------------------------------
      Use dbsr_conf 

      Implicit none
      Integer :: ic,jc
      Real(8) :: S,CC
      Real(8), external :: ZRECUP 
      Integer, external :: Iadd_cfg_jj, Iort

      if(IORT(ii_comp,ie_comp).eq.0) Return

      CC=WC_comp

      Select case(insert)

      Case(0)
       if(debug.gt.0) &
          write(pri,'(a,i3,5x,a,f10.5)') 'insert=',insert,'CC=',CC
      Case(1:) 
       Call make_coupling_insert(insert)
       Call RECUP(nmom,ncup,J1_coupling,J2_coupling)
       S = ZRECUP(nmom,moments) 
       if(S.eq.0) Return
       CC = CC*S*S
       if(debug.gt.0) write(pri,'(a,i3,2(3x,a,f8.5))') &
                     'insert=',insert,'recup=',S,'CC=',CC
      Case(:-1)
       Call make_coupling_trap
       Call RECUP(nmom,ncup,J1_coupling,J2_coupling)
       S = ZRECUP(nmom,moments) 
       if(S.eq.0) Return

S_cfp = 1.d0          !  ???
S = 1.d0              !  ???

       CC = CC * S*S * S_cfp*S_cfp   
       if(debug.gt.0) write(pri,'(a,i3,3(3x,a,f8.5))') &
             'insert=',insert,'recup=',S,'cfp=',s_cfp,'CC=',CC
      End Select

      igen_conf = 1
      jc=Iadd_cfg_jj('detect'); ic = iabs(jc)

      if(debug.gt.0) Call Print_conf_jj(pri,0,WC(ic))

      if(jc.lt.0) then
       CC = CC + WC(ic)
       if(CC.gt.c_comp) then 
        Call Iadd_orth(ii_comp,ie_comp,0)
        if(debug.gt.1) write(pri,'(32x,a)') 'IORT -> 0'
       else 
        WC(ic) = CC
       end if
      else
       WC(ic) = CC
      end if

      if(debug.gt.0) Call Print_conf_jj(pri,0,WC(ic))

      End Subroutine Check_cfg

