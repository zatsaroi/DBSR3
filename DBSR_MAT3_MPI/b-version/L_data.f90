!======================================================================
      Subroutine L_data(itype) 
!======================================================================
!     processing of L-integrals from the module 'c_data'
!----------------------------------------------------------------------
!     we have following different structures:
!
!     1.1   L( . . )  ic, jc             -  k1=i  k2=j  k3= ic  k4= jc
!     1.2   L( . . ) < i | . > ic        -  k1=i  k2=j  k3=-io  k4= ic
!     1.3   L( . . ) < i | . > < j | . > -  k1=i  k2=j  k3=-io  k4=-jo  
!     1.4   L( . . ) < i | j >           -  k1=i  k2=j  k3=-io  k4=  0
!    
!     2.1   L( i . )  ic                 -  k1=j  k2=ich  k3= ic  k4=0
!     2.2   L( i . ) < j | . >           -  k1=j  k2=ich  k3=-io  k4=0
!    
!     3.0   L( i j )                     -  k1=ich  k2=jch  k3=0  k4=0
!
!     where .  denotes bound orbital, i,j - channels.
!
!     L( . .) < i | j > elements with i<>j are ignored because
!     we assume that target states diagonalize Hamiltonian
!     L( . .) < i | i > are included after in B(.,.) * Etarget(i)
!     where B(.,.) is B-spline overlap matrix
!     These elements are used for control calculation of interaction
!     matrix between target states (Target_h).
!     L( i j) with with i/=j are also ignored because we assume
!     that target states are orthpgonal, but we have include
!     interaction with core which is included in L(i,j)
!
!----------------------------------------------------------------------
      Use dbsr_mat;  Use c_data;  Use dhl_core    
      
      Implicit none
      Integer, intent(in) :: itype
      Integer :: i,j, i1,i2, j1,j2, ich,jch, k, io,jo, ic,jc
      Real(8) :: C,S, t1,t2,t3, v(ms), w(ms), dd(ms,ms)
      Integer, external :: KBORT, Ipointer
      Real(8), external :: Lval

      if(ncdata.eq.0) Return

      Call CPU_time(t1)

      Select Case(itype)
!----------------------------------------------------------------------
      Case(1)                                                 !  L(..)  

       Do j=1,ncdata; i=IPT(j); C=cdata(i); if(abs(C).lt.Eps_C) Cycle
        i1 = k1(i); i2 = k2(i); S = Lval(i1,i2); C = C * S

        ic = k3(i); jc = k4(i)

        if(ic.gt.0.and.jc.gt.0)  then                         !  L(..)  ic,jc
         Call Update_HB(ic,jc,C)

        elseif(ic.lt.0.and.jc.gt.0) then                      !  L(..)  <i|.> jc
         io=-ic; i1=io/ibo; i2=mod(io,ibo)
         k=KBORT(i1,i2); if(k.eq.0) Stop 'O: KBORT=0'
         v = V_ch(:,k)     
         ich=ipbs(i1); Call UPDATE_HV(ich,jc,v,C)

        elseif(ic.lt.0.and.jc.lt.0) then                      !  L(..)  <i|.> <.|j>
         io=-ic; i1=io/ibo; i2=mod(io,ibo)
         k=KBORT(i1,i2); if(k.eq.0) Stop 'O: KBORT=0'
         v = V_ch(:,k)     
         jo=-jc; j1=jo/ibo; j2=mod(jo,ibo)
         k=KBORT(j1,j2); if(k.eq.0) Stop 'O: KBORT=0'
         w = V_ch(:,k)     
         ich=ipbs(i1); jch=ipbs(j1);  Call UPDATE_HW(ich,jch,v,w,C)

        elseif(ic.lt.0.and.jc.eq.0) then                      !  L(..) <i|j>
         io=-ic; i1=io/ibo; i2=mod(io,ibo)
         ich=ipbs(i1); jch=ipbs(i2)

         if(iitar.gt.0.and.ich.ne.jch) Call UPDATE_HX (ich,jch,fppqq,C)

          Call Target_h(ich,jch,C,0.d0)                       ! update htarg

        else
         Stop 'L_data:  wrong combination of ic,jc for itype=1'
        end if

        if(pri_coef.gt.0) &
        Call pri_L_coef(pri,itype,k1(i),k2(i),k3(i),k4(i),cdata(i),S)

       End do 
!-----------------------------------------------------------------------
      Case(2)                                                 ! L(.i) <.|j> 

       Do j=1,ncdata; i=IPT(j); C=cdata(i); if(abs(C).lt.Eps_C) Cycle
        i1=k1(i);   Call Gen_Lvec(i1,v)
        ich=k2(i);  jc=k3(i)

        if(jc.gt.0) then                                      ! L(.i)  jc 
         Call UPDATE_HV(ich,jc,v,C)
          
        elseif(jc.lt.0) then                                  ! L(.i) <.|j> 
         jo=-jc; j1=jo/ibo; j2=mod(jo,ibo)
         k=KBORT(j1,j2); if(k.eq.0) Stop 'O: KBORT=0'
         w = V_ch(:,k)   
         jch=ipbs(j1);  Call UPDATE_HW (ich,jch,v,w,C)
        
        else
         Stop 'L_data:  wrong value jc = 0 for itype=2'
        end if

        if(pri_coef.gt.0) &
        Call pri_L_coef(pri,itype,k1(i),k2(i),k3(i),k4(i),cdata(i),S)

       End do

!----------------------------------------------------------------------
      Case(3)                                                  ! L(ij) 

       Do j=1,ncdata; i=IPT(j); C=cdata(i); if(abs(C).lt.eps_C) Cycle
        i1=k1(i); j1=k2(i); ich=ipbs(i1); jch=ipbs(j1)
        if(kbs(i1).ne.kbs(j1)) Stop 'L_data: different kappa, itype=3'
        k = Ipointer(nkap,kap_list,kbs(i1))
        if(k.eq.0) Stop 'L_data: unknown kappa, itype=3'

! ... contribition shoud be zero for i <> j for orthogonal target states:

        if(ich.eq.jch.or.iitar.gt.1) then
         Call UPDATE_HX(ich,jch,hd(1,1,k),C)
        else
         if(check_target.eq.1)  Call Target_h(ich,jch,0.d0,C)
        end if

        if(pri_coef.gt.0) &
        Call pri_L_coef(pri,itype,k1(i),k2(i),k3(i),k4(i),cdata(i),S)

       End do

      Case Default

       Stop ' L_data: uknown itype '

      End Select

      Call CPU_time(t2);  t_ldata = t_ldata + (t2-t1)

      if(debug.gt.0.and.pri.gt.0) write(pri,'(a,i2,a,i10,f10.1,a/)') &
       'L_data:  itype = ',itype,'  nc = ', ncdata, (t2-t1),' sec'

      End Subroutine L_data


!=======================================================================
      Subroutine pri_L_coef(nu,itype,k1,k2,k3,k4,C,S)
!=======================================================================
!     debug printing the data from "c_data" list
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq
      Use channel_jj, only: elc

      Implicit none
      Integer, intent(in) :: nu,itype,k1,k2,k3,k4
      Real(8), intent(in) :: C,S
      Integer :: i1,i2,ic,jc,io,jo,io1,io2,jo1,jo2
      Integer, parameter :: ibo = 2**15   ! packing basis for orbitals

      Select case(itype)

      Case(1)
       i1 = K1; i2 = k2
       ic=k3;  jc=k4;   io=-ic; jo=-jc

       if(ic.gt.0.and.jc.gt.0) then                   ! L(..) (ic,jc)

        write(nu,'(f10.5,3x,5a,2i6,T64,D20.10)') &
         C, 'L (',ebs(i1),',',ebs(i2),')',ic,jc,S  

       elseif(ic.lt.0.and.jc.gt.0) then               ! L(..) <i|.> jc
        io1=io/ibo; io2=mod(io,ibo)  

        write(nu,'(f10.5,3x,10a,i6,T64,D20.10)') &
         C, 'L (',ebs(i1),',',ebs(i2),')', &    
         ' <',ebs(io1),'|',ebs(io2),'>',-jc,S  

       elseif(ic.lt.0.and.jc.lt.0) then               ! L(..) <i|.> <.|j>
        io1=io/ibo; io2=mod(io,ibo)  
        jo1=jo/ibo; jo2=mod(jo,ibo)  

        write(nu,'(f10.5,3x,15a,T64,D20.10)') &
         C, 'L (',ebs(i1),',',ebs(i2),')', &    
         ' <',ebs(io1),'|',ebs(io2),'>',' <',ebs(jo1),'|',ebs(jo2),'>',S

       elseif(ic.lt.0.and.jc.eq.0) then               ! L(..) <i|j>      
        io1=io/ibo; io2=mod(io,ibo)  

        write(nu,'(f10.5,3x,10a,T64,D20.10)') &
         C,'L (',ebs(i1),',',ebs(i2),')',' <',ebs(io1),'|',ebs(io2),'>', S  
 
       else

        Stop 'pri_L_coef: unknown structure for itype=1'

       end if

      Case(2)                                         !  L ( i . )

       i2 = K1; i1 = K2; io = -K3

       if(io.gt.0) then
        io1=io/ibo; io2=mod(io,ibo)  
        write(nu,'(f10.5,3x,10a)') &
         C,'L (',elc(i1),',',ebs(i2),')',' <',ebs(io1),'|',ebs(io2),'>'  
       else
        write(nu,'(f10.5,3x,5a,i6)') C,'L (',elc(i1),',',ebs(i2),')',-io
       end if
                                                      !  L ( i j )
      Case(3)
       i1 = K1; i2 = K2
       write(nu,'(f10.5,3x,5a)')  C,'L (',ebs(i1),',',ebs(i2),')'

      End Select

      End Subroutine pri_L_coef
