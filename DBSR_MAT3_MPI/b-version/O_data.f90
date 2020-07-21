!======================================================================
      Subroutine O_data(itype)
!======================================================================
!     processing of total overlaps
!----------------------------------------------------------------------
!     we have following different structures:
!
!     1  ic, jc                -  bound-bound           -  k1=ic  k2=jc 
!     2  < i | . > ic          -  bound-channel         -  k1=io  k2=ic 
!     3  < i | . > < j | . >   -  channel-channel       -  k1=io  k2=jo 
!     4  < i | j >             -  target contribution   -  k1=io  k2= 0 
!
!     where .  denotes bound orbital, i,j - channels.
!
!     < i | j > elements are ignored because we assume that target
!     states are orthogonal. These elements are used to check the 
!     orthogonality of the target states.
!----------------------------------------------------------------------
      Use dbsr_mat;  Use c_data       

      Implicit none
      Integer, intent(in) :: itype
      Integer :: i,j, i1,i2, j1,j2, ich,jch, k, io,jo, ic,jc
      Real(8) :: C, t1,t2, v(ms),w(ms)
      Integer, external :: KBORT

      if(ncdata.le.0) Return

      Call CPU_time(t1)

      Select case(itype)

      Case(1)                                          !   <ic|jc>
       Do j=1,ncdata; i=IPT(j)
        C=cdata(i); if(abs(C).lt.Eps_c) Cycle
        ic=k1(i); jc=k2(i)
        Call Update_HB(ic,jc,C)
       End do

      Case(2)                                          !   <i|.> jc
       Do j=1,ncdata; i=IPT(j)
        C=cdata(i); if(abs(C).lt.Eps_c) Cycle
        io=k1(i); jc=k2(i)
        i1=io/ibo; i2=mod(io,ibo)
        k=KBORT(i1,i2); if(k.eq.0) Stop 'O: KBORT=0'
        v = V_ch(:,k)     
        ich=ipbs(i1); Call UPDATE_HV(ich,jc,v,C)
       End do

      Case(3)                                          !   <i|.> <.|j>
       Do j=1,ncdata; i=IPT(j)
        C=cdata(i); if(abs(C).lt.Eps_c) Cycle
        io=k1(i); jo=k2(i)
        i1=io/ibo; i2=mod(io,ibo)
        k=KBORT(i1,i2); if(k.eq.0) Stop 'O: KBORT =0'
        v = V_ch(:,k)   
        j1=jo/ibo; j2=mod(jo,ibo)
        k=KBORT(j1,j2); if(k.eq.0) Stop 'O: KBORT=0'
        w = V_ch(:,k)   
        ich=ipbs(i1); jch=ipbs(j1); Call UPDATE_HW(ich,jch,v,w,C)         
       End do

      Case(4)                                          !   <i|j>

       Do j=1,ncdata; i=IPT(j)
        C=cdata(i); if(abs(C).lt.Eps_c) Cycle
        io=k1(i); i1=io/ibo; i2=mod(io,ibo)
        ich=ipbs(i1); jch=ipbs(i2)

! ... targets overlaps: should be zero for ich <> jch;
! ... for ich=jch, should be =1 and they are included in SUB1
! ... as B-spline overlaps 

        Call Target_h(ich,jch,0.d0,C)
        if(iitar.gt.1.and.ich.ne.jch) &
         Call UPDATE_HX (ich,jch,fppqq,C)

       End do

      Case default
       Stop 'O_data: unknown structure (itype) '

      End Select

      Call CPU_time(t2)

      if(debug.gt.0.and.pri.gt.0) &
      write(pri,'(a,i3,i10,f10.2,a)') 'O_data: itype, ncdata = ',itype,ncdata,(t2-t1)/60,' min.' 

      End Subroutine O_data
