!======================================================================
      Subroutine R_data (itype,kpol,btype)
!======================================================================
!     processing the Rk-integrals from the module 'cmdata'
!     for given type, multipole index and relativistic submatrix
!----------------------------------------------------------------------
!     we have following different structures: (itype=1-4)
!     (some simplifications are due to symmetry of Rk-integrals)
!
!     1.0  Rk( . . . .)  ic, jc             -  bound-bound  
!     1.1  Rk( . . . .) < i | . > ic        -  bound-channel
!     1.2  Rk( . . . .) < i | . > < j | . > -  channel-channel
!     1.3  Rk( . . . .) < i | j >           -  target structure
!     
!     2.0  Rk( i . . .) < j | . >           -  channel-channel due to
!     2.0  Rk( . i . .) < j | . >              overlaps
!     2.0  Rk( . . i .) < j | . >
!     2.0  Rk( . . . i) < j | . >
!     
!     2.1  Rk( i . . .)  ic                 -  bound-channel
!     2.1  Rk( . i . .)  ic
!     2.1  Rk( . . i .)  ic
!     2.1  Rk( . . . i)  ic
!     
!     3.0  Rk( i . j .)        int = 1      -  direct channel-channel
!     3.0  Rk( . i . j)              2
!     4.0  Rk( i . . j)              3      -  exchange channel-channel
!     4.0  Rk( . i j .)              4
!
!     where .  denotes bound orbital, i,j - channels.
!
!     Rk * < i | j >  elements with i<>j are ignored because
!     we assume that target states diagonalize Hamiltonian
!     These elements are included after in B(.,.) * Etarget(i)
!     where B(.,.) is B-spline overlap matrix
!     These elements are also used for control calculation of 
!     interaction matrix between target states (Target_h).
!
!     sym_i -  symmetry for r1 variable  ('s' or 'n')
!     sym_j -  symmetry for r2 variable
!     sym_d -  symmetry for convolution  ('s','n','x')
!     sym_r -  symmetry for result
!
!----------------------------------------------------------------------
      Use dbsr_mat, xxx => xx
      Use c_data

      Implicit none
      Integer, intent(in) :: btype, itype, kpol
      Integer :: i,j, i1,i2,j1,j2, ii,jj, iii,jjj, ich,jch, int,  &
                 ip1,ip2,jp1,jp2, io1,io2, jo1,jo2, k, ic,jc,io,jo
      Real(8) :: S,C,CA,t1,t2
      Real(8) :: v(ns), vv(ms),ww(ms), xx(ns,ns), dd(ns,ns), d(ns,ns)  
      Character :: sym_i,sym_j,sym_d,sym_r 
      Character(80) :: message 
      Real(8), external :: Sum_AmB, QUADR_PQ
      Integer, external :: KBORT

      if(ncdata.eq.0) Return

      Call CPU_time(t1)

      message = ' '
!----------------------------------------------------------------------
! ... prepare B_spline representation and define symmetries for given 
! ... type of integral: 
      
      sym_i='s'; sym_j='s'
      Select case(itype)
       Case(1); int=2; sym_d='s';  sym_r='s'     ! direct
       Case(2); int=1; sym_d='s';  sym_r='s'     ! direct
       Case(3); int=1; sym_d='s';  sym_r='s'     ! direct
       Case(4); int=3; sym_d='x';  sym_r='x'     ! exchange
      End Select

! ... define the big (1) or small (2) components involved:

      atype = btype
      Select case(atype)
       Case(0);  Call mrk_pppp(kpol); ip1=1; ip2=1; jp1=1; jp2=1  
       Case(1);  Call mrk_qqqq(kpol); ip1=2; ip2=2; jp1=2; jp2=2  
       Case(2);  Call mrk_pqpq(kpol); ip1=1; ip2=1; jp1=2; jp2=2  
                 if(itype.eq.4) then; ip1=1; ip2=2; jp1=1; jp2=2; end if
       Case(3);  Call mrk_qpqp(kpol); ip1=2; ip2=2; jp1=1; jp2=1  
                 if(itype.eq.4) then; ip1=2; ip2=1; jp1=2; jp2=1; end if    
      End Select

! ... select structure:

      Select Case(itype)
!----------------------------------------------------------------------
      Case(1)                                         ! Rk( . . ; . . )           

       iii=0; jjj=0; S=0.d0
       Do j=1,ncdata;  i=IPT(j)

        if(K1(i).ne.iii) then
         i1 = K1(i)/ibo; i2 = mod(K1(i),ibo)
         Call Density(ns,ks,d,pq(1,ip1,i1),pq(1,ip2,i2),sym_d)
         Call Convol(ns,ks,xx,d,int,sym_i,sym_j)
         if(debug.gt.0) write(message(16:),'(a)') '  convol'  
        end if

        if(K2(i).ne.jjj) then
         j1 = K2(i)/ibo; j2 = mod(K2(i),ibo)
         Call Density(ns,ks,dd,pq(1,jp1,j1),pq(1,jp2,j2),sym_r)
        end if

        if(K1(i).ne.iii.or.K2(i).ne.jjj) then
         S = SUM_AmB(ns,ks,xx,dd,sym_r)
         if(debug.gt.0) write(message(1:15),'(E15.7)') S  
        end if

        iii=K1(i);  jjj=K2(i);   C=S*cdata(i)

        ic=-k3(i);  jc=-k4(i);   io=-ic; jo=-jc

        if(ic.gt.0.and.jc.gt.0) then                  ! R(..) (ic,jc)

         Call Update_HB(ic,jc,C)

        elseif(ic.lt.0.and.jc.gt.0) then              ! R(..) <i|.> jc

         io1=io/ibo; io2=mod(io,ibo)  
         k = KBORT(io1,io2); if(k.eq.0) Stop 'RK: KBORT =0'
         vv = v_ch(:,k)
         ich=ipbs(io1)
         Call UPDATE_HV(ich,jc,vv,C)

        elseif(ic.lt.0.and.jc.lt.0) then              ! R(..) <i|.> <.|j>

         io1=io/ibo; io2=mod(io,ibo)
         k=KBORT(io1,io2); if(k.eq.0) Stop 'RK: KBORT =0'
         vv=V_ch(:,k) 
         jo1=jo/ibo; jo2=mod(jo,ibo); 
         k=KBORT(jo1,jo2); if(k.eq.0) Stop 'RK: KBORT =0'
         ww=V_ch(:,k) 
         ich = ipbs(io1); jch = ipbs(jo1) 
         Call UPDATE_HW(ich,jch,vv,ww,C)         

        elseif(ic.lt.0.and.jc.eq.0) then              ! R(..) <i|j>      

         io1 = io/ibo; io2 = mod(io,ibo)
         ich = ipbs(io1); jch = ipbs(io2)

! ... interaction between targets: 
! ... should be zero for real eigenstates;
! ... for ich=jch, corrections are through Etarg 

         Call Target_h (ich,jch,C,0.d0)

         if(iitar.gt.0.and.ich.ne.jch) Call UPDATE_HX(ich,jch,fppqq,C)

        else

         Stop ' R_data: unknown structure for itype=1'

        end if

        if(pri_coef.gt.0) &
        Call pri_RK_coef (pri,itype,atype,kpol,k1(i),k2(i),k3(i),k4(i),cdata(i),message)
        message = ' ' 
       End do   ! over ncdata

!----------------------------------------------------------------------
!                                                       R ( i . ; . . ) 
      Case(2)                               

       jjj=0; iii=0
       Do j=1,ncdata; i=IPT(j)
        jj=k1(i); ii=k2(i); ich=k3(i); io=k4(i)

        if(jj.ne.jjj.or.ii.ne.iii) then
         if(jj.ne.jjj) then
          j1=jj/ibo; j2=mod(jj,ibo)
          Call Density(ns,ks,dd,pq(1,jp1,j1),pq(1,jp2,j2),sym_d)
          Call Convol(ns,ks,xx,dd,int,sym_i,sym_j)
          message = 'convol'
          jjj = jj
         end if
         CALL BAV(ns,ks,xx,pq(1,ip1,ii),v,sym_r,'r')
         iii = ii 
        end if

        vv = 0.d0
        if(ip2.eq.1) vv(   1:ns) = v(1:ns)
        if(ip2.eq.2) vv(ns+1:ms) = v(1:ns)

        if(io.gt.0) then                       !  Rk( i . . .) < j | . >  
                                       
         io1 = io/ibo; io2 = mod(io,ibo)
         k=KBORT(io1,io2); if(k.eq.0) Stop 'RK: KBORT = 0'
         ww=V_ch(:,k)     

         jch = ipbs(io1)
         Call UPDATE_HW(ich,jch,vv,ww,cdata(i))         

        elseif(io.lt.0) then                   !  Rk( i . . .) ic 

          jc=-io; Call UPDATE_HV(ich,jc,vv,cdata(i))

        else

          Stop ' R_data: uknown structure for itype=2 '

        end if

        if(pri_coef.gt.0) &
        Call pri_Rk_coef (pri,itype,atype,kpol,k1(i),k2(i),k3(i),k4(i),cdata(i),message)
        message = ' '
       End do     

!----------------------------------------------------------------------
      Case(3)                                   !    RK ( i . ; j . )                  

       xx=0.d0;  CA=0.d0
       Do j=1,ncdata;  i=IPT(j); j1=k2(i); j2=k3(i)

        Call Density(ns,ks,dd,pq(1,jp1,j1),pq(1,jp2,j2),sym_d)      
        xx = xx + cdata(i)*dd

        message = ' ' 
        if(atype.le.1)  then                    ! check factor 2  ???
          S = QUADR_PQ(j1,j2,kpol)                 
          CA=CA+cdata(i)*S                 
          if(debug.gt.0) write(message(1:15),'(E15.7)') S  
        end if

        k=1; if(j.lt.ncdata) then; if(k1(i).eq.k1(IPT(j+1))) k=0; end if 

        if(k.eq.0.and.pri_coef.gt.0) &
        Call pri_Rk_coef (pri,itype,atype,kpol,k1(i),k2(i),k3(i),k4(i),cdata(i),message)

        if(k.eq.0) Cycle

        Call Convol(ns,ks,dd,xx,int,sym_i,sym_j)

        ich=k1(i)/ibo; jch=mod(k1(i),ibo)
        Call UPDATE_HS(ich,jch,dd,sym_r,ip1,ip2)
        if(atype.le.1) Call UPDATE_ACF(kpol,ich,jch,CA)  
        xx=0.d0; CA=0.d0

        write(message(16:),'(a)') '  convol' 
        if(pri_coef.gt.0) &
        Call pri_Rk_coef (pri,itype,atype,kpol,k1(i),k2(i),k3(i),k4(i),cdata(i),message)

       End do

!----------------------------------------------------------------------
      Case(4)                                     !    RK ( i . ; . j )                  

       xx=0.d0 
       Do j=1,ncdata;  i=IPT(j); j1=k2(i); j2=k3(i)

        Call Density(ns,ks,dd,pq(1,jp1,j1),pq(1,jp2,j2),sym_d)      
        xx = xx + cdata(i)*dd

        k=1; if(j.lt.ncdata) then; if(k1(i).eq.k1(IPT(j+1))) k=0; end if 

        message = ' '
        if(k.eq.0.and.pri_coef.gt.0) &
        Call pri_Rk_coef (pri,itype,atype,kpol,k1(i),k2(i),k3(i),k4(i),cdata(i),message)

        if(k.eq.0) Cycle

        Call Convol(ns,ks,dd,xx,int,sym_i,sym_j)

        ich=k1(i)/ibo; jch=mod(k1(i),ibo)
        Call UPDATE_HS(ich,jch,dd,sym_r,ip1,ip2)

        xx=0.d0

        message = '  convol' 
        if(pri_coef.gt.0) &
        Call pri_Rk_coef (pri,itype,atype,kpol,k1(i),k2(i),k3(i),k4(i),cdata(i),message)

       End do

      Case Default          

       Stop ' R_data: unknown itype '

      End Select    ! over itype

      Call CPU_time(t2);  t_rdata = t_rdata + (t2-t1)

      if(debug.gt.0.and.pri.gt.0) write(pri,'(a,i2,a,i1,a,i2,a,i10,f10.1,a)') &
       'R_data:  itype = ',itype,' atype = ',atype,'  kpol = ',kpol,'  nc = ', ncdata,(t2-t1),' sec'

      End Subroutine  R_data



!=======================================================================
      Subroutine pri_RK_coef(nu,itype,atype,k,k1,k2,k3,k4,C,message)
!=======================================================================
!     debug printing the data from "c_data" list
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq
      Use channel_jj, only: elc

      Implicit none
      Integer, intent(in) :: nu,itype,atype,k,k1,k2,k3,k4
      Real(8), intent(in) :: C
      Character(1) :: P1,P2,P3,P4
      Character(*) :: message
      Integer :: i1,i2,i3,i4,j1,j2,ic,jc,io,jo,io1,io2,jo1,jo2
      Integer, parameter :: ibo = 2**15   ! packing basis for orbitals

      Select case(atype)
       Case(0);  P1='P'; P2='P'; P3='P'; P4='P'  
       Case(1);  P1='Q'; P2='Q'; P3='Q'; P4='Q'  
       Case(2);  P1='P'; P2='Q'; P3='P'; P4='Q'  
       Case(3);  P1='Q'; P2='P'; P3='Q'; P4='P'  
      End Select
      
      Select case(itype)

      Case(1)
       i1 = K1/ibo; i2 = mod(K1,ibo)
       j1 = K2/ibo; j2 = mod(K2,ibo)
       ic=-k3;  jc=-k4;   io=-ic; jo=-jc

       if(ic.gt.0.and.jc.gt.0) then                   ! R(..) (ic,jc)

        write(nu,'(f10.5,3x,a,i2,13a,2i6,T78,a)') &
         C, 'R',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')',ic,jc,trim(message)  

       elseif(ic.lt.0.and.jc.gt.0) then               ! R(..) <i|.> jc
        io1=io/ibo; io2=mod(io,ibo)  

        write(nu,'(f10.5,3x,a,i2,18a,i6,T78,a)') &
         C, 'R',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
        -jc,trim(message)  

       elseif(ic.lt.0.and.jc.lt.0) then               ! R(..) <i|.> <.|j>
        io1=io/ibo; io2=mod(io,ibo)  
        jo1=jo/ibo; jo2=mod(jo,ibo)  

        write(nu,'(f10.5,3x,a,i2,23a,T78,a)') &
         C, 'R',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
          ' <',ebs(jo1),'|',ebs(jo2),'>',trim(message)

       elseif(ic.lt.0.and.jc.eq.0) then               ! R(..) <i|j>      
        io1=io/ibo; io2=mod(io,ibo)  

        write(nu,'(f10.5,3x,a,i2,18a,T78,a)') &
         C, 'R',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
         trim(message)  
 
       else

        Stop 'pri_coef: unknown structure for itype=1'

       end if

      Case(2)                                       !   R ( i . ; . . )

       j1 = K1/ibo; j2 = mod(K1,ibo); i2 = K2; i1 = K3;  io = K4
       if(io.gt.0) then
        io1=io/ibo; io2=mod(io,ibo)  
        write(nu,'(f10.5,3x,a,i2,18a,3x,a)') &
         C, 'R',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
         trim(message)  
       else
        write(nu,'(f10.5,3x,a,i2,13a,i6,3x,a)') &
         C, 'R',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')',-io,trim(message)  
       end if

       Case(3)
        i1 = K1/ibo; i2 = mod(K1,ibo);  j1=K2; j2=K3
        write(nu,'(f10.5,3x,a,i2,13a,5x,a)') &
         C, 'R',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,elc(i2),',',P4,ebs(j2),')',trim(message)  

       Case(4)
        i1 = K1/ibo; j2 = mod(K1,ibo);  i2=K2;  j1=K3
        write(nu,'(f10.5,3x,a,i2,13a,5x,a)') &
         C, 'R',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,elc(j2),')',trim(message)  

      End Select

      End Subroutine pri_Rk_coef


