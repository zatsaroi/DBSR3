!======================================================================
      Subroutine S_data (ktype,kpol)
!======================================================================
!     processing the Sk-integrals from the module 'cmdata'
!----------------------------------------------------------------------
!     we have following different structures: 
!
!     1.0  Sk( . . . .)  ic, jc               -  bound-bound  
!     1.1  Sk( . . . .) < i | . > ic          -  bound-channel
!     1.2  Sk( . . . .) < i | . > < j | . >   -  channel-channel
!     1.3  Sk( . . . .) < i | j >             -  target structure
!    
!     2.0  Sk( i . . .) < j | . >             -  channel-channel due to
!     3.0  Sk( . i . .) < j | . >                overlaps
!     4.0  Sk( . . i .) < j | . >
!     5.0  Sk( . . . i) < j | . >
!    
!     2.1  Sk( i . . .)  ic                   -  bound-channel
!     3.1  Sk( . i . .)  ic
!     4.1  Sk( . . i .)  ic
!     5.1  Sk( . . . i)  ic
!    
!     6.0  Sk( i . j .)                       -  direct channel-channel
!     7.0  Sk( . i . j)
!     8.0  Sk( i . . j)                       -  exchange channel-channel
!     9.0  Sk( . i j .)
!
!     where .  denotes bound orbital, i,j - channels.
!
!     Sk * < i | j >  elements with i<>j are ignored because
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
!     int = 1 -->  INT( . p1; . p2)
!     int = 2 -->  INT( p1 .;p2 . )
!     int = 3 -->  INT( . p1;p2 . )
!     int = 4 -->  INT( p1 .; . p2)
!
!     where p1,p2 - known bound orbitals
!----------------------------------------------------------------------
      Use dbsr_mat, xxx => xx
      Use c_data          

      Implicit none
      Integer :: itype, kpol, ktype
      Integer :: i,j, i1,i2, j1,j2, ii,jj, iii,jjj, ich,jch, int, k, &
                 ip1=1,ip2=1,jp1=1,jp2=1, io,jo,ic,jc
      Real(8) :: S,C,CA,t1,t2
      Real(8) :: v(ns), vv(ms),ww(ms), xx(ns,ns), dd(ns,ns), d(ns,ns)  
      Character :: sym_d, sym_r, sym_v, sym_i, sym_j
      Character(80) :: message = ' '
      Integer, external :: KBORT
      Real(8), external :: Sum_AmB, QUADR_qp

      if(ncdata.eq.0) Return

      Call CPU_time(t1)

      atype = ktype/ibtype; itype=mod(ktype,ibtype) 

!----------------------------------------------------------------------
! ... prepare B_spline representation and define symmetries for given 
! ... type of integrals, indicated by 'icase' in module c_data: 
      
      Select case(atype)
       Case(0);  Call msk_ppqq(kpol)                                                 
        Select case(itype)         
         Case(1);     ip1=1; ip2=2; jp1=1; jp2=2;  int=2             !  I(ip1,jp1;ip2,jp2)   ( . i . j)
         Case(2);     ip1=1; ip2=2; jp1=1; jp2=2;  int=1; sym_v='r'  !  I(ip1,jp1;ip2,jp1)   ( i . . .)
         Case(3);     ip1=1; ip2=2; jp1=1; jp2=2;  int=2; sym_v='r'  !  I(jp1,ip1;jp2,ip2)   ( . i . .)
         Case(4);     ip1=2; ip2=1; jp1=1; jp2=2;  int=1; sym_v='l'  !  I(ip2,jp1;ip1,jp2)   ( . . i .)
         Case(5);     ip1=2; ip2=1; jp1=1; jp2=2;  int=2; sym_v='l'  !  I(jp1,ip2;jp2,ip1)   ( . . . i)
         Case(6);     ip1=1; ip2=2; jp1=1; jp2=2;  int=1;            !  I(ip1,jp1;ip2,jp1)   ( i . j .)
         Case(7);     ip1=1; ip2=2; jp1=1; jp2=2;  int=2;            !  I(jp1,ip1;jp2,ip2)   ( . i . j)
         Case(8);     ip1=1; ip2=2; jp1=2; jp2=1;  int=3;            !  I(ip1,jp2;jp1,ip2)   ( i . . j)
         Case(9);     ip1=2; ip2=1; jp1=1; jp2=2;  int=4;            !  I(jp1,ip2;ip1,jp2)   ( . i j .)
        End Select                                 
       Case(1);  Call msk_pqqp(kpol)               
        Select case(itype)                         
         Case(1);     ip1=1; ip2=2; jp1=2; jp2=1;  int=2             !  I(ip1,jp1;ip2,jp2)   ( . i . j)
         Case(2);     ip1=1; ip2=2; jp1=2; jp2=1;  int=1; sym_v='r'  !  I(ip1,jp1;ip2,jp1)   ( i . . .)
         Case(3);     ip1=2; ip2=1; jp1=1; jp2=2;  int=2; sym_v='r'  !  I(jp1,ip1;jp2,ip2)   ( . i . .)
         Case(4);     ip1=2; ip2=1; jp1=2; jp2=1;  int=1; sym_v='l'  !  I(ip2,jp1;ip1,jp2)   ( . . i .)
         Case(5);     ip1=1; ip2=2; jp1=1; jp2=2;  int=2; sym_v='l'  !  I(jp1,ip2;jp2,ip1)   ( . . . i)
         Case(6);     ip1=1; ip2=2; jp1=2; jp2=1;  int=1;            !  I(ip1,jp1;ip2,jp2)   ( i . j .)
         Case(7);     ip1=2; ip2=1; jp1=1; jp2=2;  int=2;            !  I(jp1,ip1;jp2,ip2)   ( . i . j)
         Case(8);     ip1=1; ip2=1; jp1=2; jp2=2;  int=3;            !  I(ip1,jp2;jp1,ip2)   ( i . . j)
         Case(9);     ip1=2; ip2=2; jp1=1; jp2=1;  int=4;            !  I(jp1,ip2;ip1,jp2)   ( . i j .)
        End Select                                 
      End Select                                   
                                                   
      sym_i='n'; sym_j='n'    
      Select case(int)
       Case(1); sym_d=sym_j; sym_r=sym_i       ! direct
       Case(2); sym_d=sym_i; sym_r=sym_j 
       Case(3); sym_d='x';   sym_r='x'         ! exchange
       Case(4); sym_d='x';   sym_r='x' 
      End Select

!----------------------------------------------------------------------
! ... select structure:

      Select Case(itype)

!----------------------------------------------------------------------
!                                                        S( . . ; . . )      
      Case(1)                                  

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

        if(ic.gt.0.and.jc.gt.0) then                    ! S(..) (ic,jc)

         Call Update_HB(ic,jc,C)

        elseif(ic.lt.0.and.jc.gt.0) then                ! S(..) <i|.> jc

         i1 = io/ibo; i2 = mod(io,ibo)
         k=KBORT(i1,i2); if(k.eq.0) Stop 'S: KBORT = 0'
         vv=V_ch(:,k)     
         ich=ipbs(i1)
         Call UPDATE_HV(ich,jc,vv,C)

        elseif(ic.lt.0.and.jc.lt.0) then                ! S(..) <i|.> <.|j>
         i1=io/ibo; i2=mod(io,ibo)
         k=KBORT(i1,i2); if(k.eq.0) Stop 'S: KBORT = 0'
         vv=V_ch(:,k)     
         j1=jo/ibo; j2=mod(jo,ibo) 
         k=KBORT(j1,j2); if(k.eq.0) Stop 'S: KBORT = 0'
         ww=V_ch(:,k)   
         ich = ipbs(i1) 
         jch = ipbs(j1) 
         Call UPDATE_HW(ich,jch,vv,ww,C)         

        elseif(ic.lt.0.and.jc.eq.0) then                ! S(..) <i|j>      

         i1=io/ibo; i2=mod(io,ibo)
         ich = ipbs(i1); jch = ipbs(i2)

! ... interaction between targets: 
! ... should be zero for real eigenstates;
! ... for ich=jch, corrections are through Etarg 

         Call Target_h(ich,jch,C,0.d0)

         if(iitar.gt.0.and.ich.ne.jch) Call UPDATE_HX(ich,jch,fppqq,C)

        else

         Stop ' S_data: unknown structure for itype=1'

        end if

        if(debug.gt.0) Call pri_RK_coef (nud,itype,atype,kpol, &
                       k1(i),k2(i),k3(i),k4(i),cdata(i),message)
        message = ' ' 

       End do   ! over ncdata

!----------------------------------------------------------------------
!                                                        S( i . ; . . ) 
      Case(2,3,4,5)                               

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
         CALL BAV(ns,ks,xx,pq(1,ip2,ii),v,sym_r,sym_v)
         iii = ii 
        end if

        vv = 0.d0
        if(ip2.eq.1) vv(   1:ns) = v(1:ns)
        if(ip2.eq.2) vv(ns+1:ms) = v(1:ns)

        if(io.gt.0) then                            !  Sk( i . . .) < j | . >  

         j1 = io/ibo; j2 = mod(io,ibo)
         k=KBORT(j1,j2); if(k.eq.0) Stop 'S: KBORT = 0'
         ww=V_ch(:,k)   
  
         jch = ipbs(j1)
         Call UPDATE_HW(ich,jch,vv,ww,cdata(i))         
                                                   
        elseif(io.lt.0) then                        !  Sk( i . . .) ic 

          jc=-io; Call UPDATE_HV(ich,jc,vv,cdata(i))

        else

          Stop ' S_data: uknown structure for itype 2-5 '

        end if

        if(debug.gt.0) Call pri_Sk_coef (nud,itype,atype,kpol, &
                       k1(i),k2(i),k3(i),k4(i),cdata(i),message)
        message = ' '

        End do     

!----------------------------------------------------------------------
!                                                      Sk ( i . ; j . )
      Case(6,7,8,9)        

       xx=0.d0;  CA=0.d0
       Do j=1,ncdata;  i=IPT(j); i1=k2(i); j1=k3(i)

        Call Density(ns,ks,dd,pq(1,jp1,i1),pq(1,jp2,j1),sym_d)         
        xx = xx + cdata(i)*dd

        message = ' ' 

!        if(int.le.2) then                             magnetic asymptotic ???
!         S = QUADR_QP(i1,j1,kpol,jp1,jp2)
!         CA = CA + cdata(i)*S                 
!         if(debug.gt.0) write(message(1:15),'(E15.7)') S  
!        end if

        k=1; if(j.lt.ncdata) then; if(k1(i).eq.k1(IPT(j+1))) k=0; end if 

        if(debug.gt.0) then
          if(k.eq.1) write(message(16:),'(a)') '  convol' 
          Call pri_Sk_coef (nud,itype,atype,kpol, &
                            k1(i),k2(i),k3(i),k4(i),cdata(i),message)
        end if

        if(k.eq.0) Cycle

        Call Convol(ns,ks,dd,xx,int,sym_i,sym_j)

        ich=k1(i)/ibo; jch=mod(k1(i),ibo)
        Call UPDATE_HS(ich,jch,dd,sym_r,ip1,ip2)

!        if(int.le.2) Call UPDATE_ACF(kpol,ich,jch,CA)

        xx=0.d0; CA=0.d0
       End do

      Case Default

       Stop ' S_data: unknown itype '

      End Select

      Call CPU_time(t2);  t_sdata = t_sdata + (t2-t1)

      if(debug.gt.0) write(nud,'(a,i2,a,i1,a,i10,T78,a,f10.4/)') &
       'S_data:  itype = ',itype,' atype = ',atype,'  nc = ', ncdata,'  time =',(t2-t1)

      End Subroutine  S_data



!=======================================================================
      Subroutine pri_SK_coef(nu,itype,atype,k,k1,k2,k3,k4,C,message)
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
       Case(0);  P1='P'; P2='P'; P3='Q'; P4='Q'  
       Case(1);  P1='P'; P2='Q'; P3='Q'; P4='P'  
      End Select
      
      Select case(itype)
!----------------------------------------------------------------------
      Case(1)
       i1 = K1/ibo; i2 = mod(K1,ibo)
       j1 = K2/ibo; j2 = mod(K2,ibo)
       ic=-k3;  jc=-k4;   io=-ic; jo=-jc

       if(ic.gt.0.and.jc.gt.0) then                   ! R(..) (ic,jc)

        write(nu,'(f10.5,3x,a,i2,13a,2i6,T78,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')',ic,jc,trim(message)  

       elseif(ic.lt.0.and.jc.gt.0) then               ! R(..) <i|.> jc
        io1=io/ibo; io2=mod(io,ibo)  

        write(nu,'(f10.5,3x,a,i2,18a,i6,T78,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
        -jc,trim(message)  

       elseif(ic.lt.0.and.jc.lt.0) then               ! R(..) <i|.> <.|j>
        io1=io/ibo; io2=mod(io,ibo)  
        jo1=jo/ibo; jo2=mod(jo,ibo)  

        write(nu,'(f10.5,3x,a,i2,23a,T78,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
          ' <',ebs(jo1),'|',ebs(jo2),'>',trim(message)

       elseif(ic.lt.0.and.jc.eq.0) then               ! R(..) <i|j>      
        io1=io/ibo; io2=mod(io,ibo)  

        write(nu,'(f10.5,3x,a,i2,18a,T78,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
         trim(message)  
 
       else

        Stop 'pri_Rk_coef: unknown structure for itype=1'

       end if
!-----------------------------------------------------------------------------
      Case(2)   !  Rk( i . . .) < j | . >  -  k1=(i2,i4)  k2=i3  k3=ich  k4=io

       j1 = K1/ibo; j2 = mod(K1,ibo); i2 = K2; i1 = K3;  io = K4
       if(io.gt.0) then
        io1=io/ibo; io2=mod(io,ibo)  
        write(nu,'(f10.5,3x,a,i2,18a,3x,a)') &
         C, 'S',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
         trim(message)  
       else
        write(nu,'(f10.5,3x,a,i2,13a,i6,3x,a)') &
         C, 'S',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')',-io,trim(message)  
       end if
!-----------------------------------------------------------------------------
      Case(3)   !  Rk( . i . .) < j | . >  -  k1=(i1,i3)  k2=i4  k3=ich  k4=io

       i1 = K1/ibo; i2 = mod(K1,ibo); j2 = K2; j1 = K3;  io = K4
       if(io.gt.0) then
        io1=io/ibo; io2=mod(io,ibo)  
        write(nu,'(f10.5,3x,a,i2,18a,3x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,elc(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
         trim(message)  
       else
        write(nu,'(f10.5,3x,a,i2,13a,i6,3x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,elc(j1),';', &    
         P3,ebs(i2),',',P4,ebs(j2),')',-io,trim(message)  
       end if
!-----------------------------------------------------------------------------
      Case(4)   !  Rk( . . i .) < j | . >  -  k1=(i2,i4)  k2=i1  k3=ich  k4=io

       j1 = K1/ibo; j2 = mod(K1,ibo); i1 = K2; i2 = K3;  io = K4
       if(io.gt.0) then
        io1=io/ibo; io2=mod(io,ibo)  
        write(nu,'(f10.5,3x,a,i2,18a,3x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,elc(i2),',',P4,ebs(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
         trim(message)  
       else
        write(nu,'(f10.5,3x,a,i2,13a,i6,3x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,elc(i2),',',P4,ebs(j2),')',-io,trim(message)  
       end if

!-----------------------------------------------------------------------------
      Case(5)   !  Rk( . . . i) < j | . >  -  k1=(i1,i3)  k2=i2  k3=ich  k4=io

       i1 = K1/ibo; i2 = mod(K1,ibo); j1 = K2; j2 = K3;  io = K4
       if(io.gt.0) then
        io1=io/ibo; io2=mod(io,ibo)  
        write(nu,'(f10.5,3x,a,i2,18a,3x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,elc(j2),')', ' <',ebs(io1),'|',ebs(io2),'>', &
         trim(message)  
       else
        write(nu,'(f10.5,3x,a,i2,13a,i6,3x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,elc(j2),')',-io,trim(message)  
       end if
!-----------------------------------------------------------------------------
       Case(6)  ! Rk( i . j .)  -  k1=(ich1,ich2)  k2=i2  k3=i4  

        i1 = K1/ibo; i2 = mod(K1,ibo);  j1=K2; j2=K3
        write(nu,'(f10.5,3x,a,i2,13a,5x,a)') &
         C, 'S',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,elc(i2),',',P4,ebs(j2),')',trim(message)  
!-----------------------------------------------------------------------------
       Case(7)  ! Rk( . i . j)  -  k1=(ich1,ich2)  k2=i1  k3=i3

        j1 = K1/ibo; j2 = mod(K1,ibo);  i1=K2; i2=K3
        write(nu,'(f10.5,3x,a,i2,13a,5x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,elc(j1),';', &    
         P3,ebs(i2),',',P4,elc(j2),')',trim(message)  
!-----------------------------------------------------------------------------
       Case(8)  ! Rk( i . . j)  -  k1=(ich1,ich2)  k2=i2  k3=i3  

        i1 = K1/ibo; j2 = mod(K1,ibo);  j1=K2; i2=K3
        write(nu,'(f10.5,3x,a,i2,13a,5x,a)') &
         C, 'S',k,' (',P1,elc(i1),',',P2,ebs(j1),';', &    
         P3,ebs(i2),',',P4,elc(j2),')',trim(message)  
!-----------------------------------------------------------------------------
       Case(9)  ! Rk( . i j .)  -  k1=(ich1,ich2)  k2=i1  k3=i4

        j1 = K1/ibo; i2 = mod(K1,ibo);  i1=K2;  j2=K3
        write(nu,'(f10.5,3x,a,i2,13a,5x,a)') &
         C, 'S',k,' (',P1,ebs(i1),',',P2,elc(j1),';', &    
         P3,elc(i2),',',P4,ebs(j2),')',trim(message)  

      End Select

      End Subroutine pri_Sk_coef


