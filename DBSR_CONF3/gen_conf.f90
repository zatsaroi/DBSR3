!====================================================================
      Subroutine Gen_Conf 
!====================================================================
!     generates all configurations from configuration given in 
!     module conf_jj by adding the orbital 'ie' 
!--------------------------------------------------------------------
      Use dbsr_conf

      Implicit none
      Integer :: ie,i,j,n,l,k,iset,ii, JP,JV, i1 
      Integer, external :: j_kappa, l_kappa, Ifind_jjorb, after

      Call make_coupling

      ie = ie_comp
      n=NEF(ie); k=KEF(ie); iset=IEF(ie); j=j_kappa(k); l=l_kappa(k)

! ... find if the same orbital is already in the configuration:

      ii=0                                
      no=no-1         
      Do i=1,no
       if(n.ne.nn(i)) Cycle
       if(k.ne.kn(i)) Cycle
       if(iset.ne.in(i)) Cycle        
       ii=i; Exit
      End do

      if(ii.gt.0) then   ! ... case of orbital traping on the existing shell:

       if(iq(ii)+1.le.jn(ii)+1) then
        iq(ii)=iq(ii)+1
        JP=Jshell(ii); JV=Vshell(ii); JP_trap=JP
        insert = -ii; S_cfp = 1.d0; S_recup = 1.d0
        Call Sum_Term (ii,JP,JV)
       end if 
       
! ... find the position for adding orbital using "AFTER" conditions: 

      else                   

       Do i=1,no
        jp = Ifind_jjorb(nn(i),kn(i),in(i),1)
        if(AFTER(jp,ie).ge.0) Cycle
        ii=i; Exit
       End do  

       if(ii.eq.0) then       ! ... add the new top shell        

        no=no+1; nn(no)=n; in(no)=iset
        insert = 0; S_cfp = 1.d0; S_recup = 1.d0 
        Call Check_cfg

       else         ! ... insert new shell                        

       Do i=no,ii,-1;  i1=i+1
        nn(i1)=nn(i); kn(i1)=kn(i); jn(i1)=jn(i); ln(i1)=ln(i); 
        iq(i1)=iq(i); in(i1)=in(i); Jshell(i1)=Jshell(i) 
        Vshell(i1)=Vshell(i); Jintra(i1)=Jintra(i)
       End do
       nn(ii)=n; kn(ii)=k; jn(ii)=j; ln(ii)=l; iq(ii)=1; in(ii)=iset
       no=no+1
       insert = ii; S_cfp = 1.d0; S_recup = 1.d0
       Call Sum_Term (ii,0,0)

       end if   ! new shell case
      end if    ! trap case
      
      End Subroutine Gen_Conf


!=======================================================================
      Subroutine Sum_Term(i,JP,VP)
!=======================================================================
!     generates all terms for shell 'i' with parent term JP,JV 
!-----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none
      Integer :: i,j,ii,i1,i2,i3,i4,mt, JP,VP
      Integer, external :: Jterm
      Real(8), external :: cfp_jj
 
      mt = Jterm(jn(i),iq(i),-1,i1,i2,i3,i4)
      Do j=1,mt
       ii=Jterm(jn(i),iq(i),j,Jshell(i),Vshell(i),i3,i4)
       S_cfp = 1.d0
       if(insert.lt.0) &
       S_cfp = cfp_jj(jn(i),iq(i),JP,VP,Jshell(i),Vshell(i))
       CALL Sum_Iterm(i)
      End do

      End Subroutine Sum_Term


!=======================================================================
      Subroutine Sum_Iterm(istart)
!=======================================================================
!     generates all intermediate terms (begining from shell 'istart')
!     which are consistent with the total momentum JTOTAL
!-----------------------------------------------------------------------
      Use dbsr_conf

      Implicit none
      Integer, external :: Jadd_cfg_jj 
      Integer :: Jmin(msh),Jmax(msh)
      Integer :: istart,i,ii,i1,i2,j1,j2

      i1 = istart                  ! i1 - low  limit
      i2 = no                      ! i2 - high limit in array LS(...)
      if(istart.eq.1) then
       i1=2; Jintra(1)=Jshell(1)
      end if

      if(no.eq.1) then
       if(Jintra(no).eq.Jtotal)  Call Check_cfg    
       Return
      end if
       
      i=i1
    1 j1=i-1
      j2=i

      Jmin(i)=IABS(Jintra(j1)-Jshell(j2))
      Jmax(i)=     Jintra(j1)+Jshell(j2)
      Jintra(i)=Jmin(i)

    2 if(i.lt.i2) then
        i=i+1
        go to 1
      else
        if(Jintra(no).eq.Jtotal) Call Check_cfg
      end if

    3 if(Jintra(i).lt.Jmax(i)) then
        Jintra(i) = Jintra(i) + 2
        go to 2
      else
        if(i.eq.i1) go to 4
        i=i-1
        go to 3
      end if

    4 Return

      End Subroutine Sum_Iterm
