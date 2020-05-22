!======================================================================
      Subroutine Def_orth_cond
!======================================================================
!     define orthogonal conditions to avoid "over-compensation"
!----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none
      Integer :: ich,it,ii,ie,ic,ic1,ic2, io 
      Integer, external :: Iort

      if(debug.gt.0) write(pri,'(/a/)') 'Define orth.condotions:'
!      if(debug.gt.0) write(pri,*) ncfg_pert,mcfg,nch

! ... find orth. conditions by checking the compensation configurations:

      if(mcfg.gt.ncfg_pert) WC(ncfg_pert+1:mcfg)=0.d0

      Do ich = 1,nch;  it=iptar(ich)
        ii = ipch(ich); ii_comp=ii       ! position of channel orbital

        ic1=1; if(ich.gt.1) ic1=ipconf(ich-1)+1; ic2=ipconf(ich)
        ic1=ic1+ncfg_targ; ic2=ic2+ncfg_targ

        Do ic = ic1,ic2
         Do io=1,nphys_sub; ie=jp_sub(io); ie_comp=ie
!if(debug.gt.0) 'ich,ic,ii,ie,iort',ich,ic,ii,ie,iort(ii,ie)
          if(Iort(ii,ie).le.0) Cycle

          if(debug.gt.0) &
          write(pri,'(a,i4,2x,2a8)')  'it=',it,ELF(ii),ELF(ie) 

          Call Get_cfg_jj(ic)
          WC_comp = WC(ic)*WC(ic)     

          if(debug.gt.0) Call Print_conf_jj(pri,ic,WC_comp)

          igen_conf = 0
          Call Gen_conf
          if(igen_conf.eq.0) Call Iadd_orth(ii,ie,0)

          if(debug.gt.0) &
            write(pri,'(32x,2(a,i2))') 'IORT = ',IORT(ii,ie)
          
         End do
        End do ! ic, over configurations 

      End do  ! ich, over channels

      ncfg_comp = ncfg;  lcfg_comp=lcfg
      if(pri.gt.0) &
      write(pri,'(/a,T40,i8)') &
       'Number of compensation configurations:', ncfg_comp-ncfg_pert

      if(ncfg_comp-ncfg_pert.le.0) Return

      write(AF,'(a,i3.3,a)') 'pert_comp.',ilsp
      Open(nuc,file=AF)
!      write(nuc,'(a15,f16.8)')  'Core subshells:',Etarg(1)
!      write(nuc,'(a)') core(1:ncore*5) 
!      write(nuc,'(a)') 'Peel subshells:' 
!      write(nuc,'(20a5)') (ELF(i),i=ncore+1,nwf)
      write(nuc,'(a)') 'CSF(s):' 
      Do ic=ncfg_pert+1,ncfg
        Call Print_conf_jj (nuc,ic,WC(ic))
      End do
      write(nuc,'(a)') '***'
      Close(nuc)

      End Subroutine Def_orth_cond




















