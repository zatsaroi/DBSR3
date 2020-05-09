!======================================================================
      Subroutine Idef_Otype (C,ic,jc,io,jo)
!======================================================================
!     we have following 16 different structures for overlaps:
!
! 1   ic, jc                -  bound-bound             -  k1=ic  k2=jc 
! 2   < i | . > ic          -  bound-channel           -  k1=io  k2=ic 
! 3   < i | . > < j | . >   -  channel-channel         -  k1=io  k2=jo 
! 4   < i | j >             -  target contribution     -  k1=io  k2= 0 
!
!     where .  denotes bound orbital, i,j - channels.
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ic,jc,io,jo
      Real(8), intent(in) :: C
      Integer :: k1,k2,k3,k4,itype

      if(abs(C).lt.eps_c) Return

      if(jc.gt.0.and.ic.le.0) Stop 'jc > 0, ic <= 0'

      itype=0

      if(jc.gt.0) then

       k1=ic; k2=jc; k3=0; k4=0; itype = 1

      elseif(ic.gt.0.and.io.gt.0) then

       k1=io; k2=ic; k3=0; k4=0; itype = 2

      elseif(io.gt.0.and.jo.gt.0) then

       k1=io; k2=jo; k3=0; k4=0; itype = 3

      elseif(io.gt.0.and.jo.eq.0) then

       k1=io; k2=0;  k3=0; k4=0; itype = 4
 
      end if

      if(itype.eq.0) Stop 'Idef_Otype: itype=0'

      if(check_target.eq.0.and.itype.eq.4) Return

      Call Add_coef(C,0,k1,k2,k3,k4,itype)

      End Subroutine Idef_Otype


!======================================================================
      Subroutine Idef_Ltype (ii,jj,C,ic,jc,io,jo)
!======================================================================
!     we have following different structures for one-electron integrals:
!
! 1.1   L( . . )  ic, jc               -  bound-bound  
! 1.2   L( . . ) < i | . > ic          -  bound-channel
! 1.3   L( . . ) < i | . > < j | . >   -  channel-channel
! 1.4   L( . . ) < i | j >             -  target contribution
!
! 2.1   L( i . ) < j | . >             -  channel-channel due to overlaps
! 2.2   L( i . )  ic                   -  bound-channel
!
! 3.0   L( i j )                       -  channel-channel
!
!     where .  denotes bound orbital, i,j - channels.
!
! 1.1   L( . . )  ic, jc               -  k1=i  k2=j  k3= ic  k4= jc
! 1.2   L( . . ) < i | . > ic          -  k1=i  k2=j  k3=-io  k4= ic
! 1.3   L( . . ) < i | . > < j | . >   -  k1=i  k2=j  k3=-io  k4=-jo  
! 1.4   L( . . ) < i | j >             -  k1=i  k2=j  k3=-io  k4=  0
!
! 2.1   L( i . )  ic                   -  k1=j  k2=ich  k3= ic  k4=0
! 2.2   L( i . ) < j | . >             -  k1=j  k2=ich  k3=-io  k4=0
!
! 3.0   L( i j )                       -  k1=ich  k2=jch  k3=0  k4=0
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: ii,jj,ic,jc,io,jo
      Real(8) :: C
      Integer :: k1,k2,k3,k4,i,j,ich,jch,itype

      if(abs(C).lt.eps_C) Return

      if(jc.gt.0.and.ic.le.0) Stop 'jc>0, ic<=0'
      if(jo.gt.0.and.io.le.0) Stop 'jo>0, io<=0'
      if(io.gt.0.and.jc.ne.0) Stop 'jc>0, io>0'
      if(jo.gt.0.and.ic.ne.0) Stop 'jc>0, jo>0'
     
      itype=0
      i=min(ii,jj); j=max(ii,jj); ich=ipbs(i); jch=ipbs(j); 

      if(jc.gt.0) then

       k1=i; k2=j; k3=ic; k4=jc;                      itype = 1

      elseif(ic.gt.0) then

       if(io.eq.0) then
        if(ich.ne.0.and.jch.eq.0) then
         k1=j; k2=ich; k3=ic; k4=0;                   itype = 2
        elseif(ich.eq.0.and.jch.ne.0) then
         k1=i; k2=jch; k3=ic; k4=0;                   itype = 2
        end if
       elseif(io.gt.0) then
        k1=i; k2=j; k3=-io; k4=ic;                    itype = 1
       end if

      elseif(io.eq.0) then

        k1=j; k2=i;  k3=0; k4=0;                      itype = 3
                                           
      elseif(io.gt.0.and.jo.eq.0) then

        if(ich.ne.0.and.jch.eq.0) then
         k1=j; k2=ich; k3=-io; k4=0;                  itype = 2
        elseif(jch.ne.0.and.ich.eq.0) then
         k1=i; k2=jch; k3=-io; k4=0;                  itype = 2
        elseif(ich.eq.0.and.jch.eq.0) then
         k1=i; k2=j; k3=-io; k4=0;                    itype = 1
         if(check_target.eq.0) Return
        end if

      elseif(io.gt.0.and.jo.gt.0) then

       k1=i; k2=j; k3=-io; k4=-jo;                    itype = 1

      end if

      if(itype.eq.0) Stop 'Idef_Ltype: itype=0'

      Call Add_coef(C,0,k1,k2,k3,k4,itype)

      End Subroutine Idef_Ltype


!======================================================================
      Subroutine Idef_Rtype (k,j1,j2,j3,j4,C,ic,jc,io,jo)
!======================================================================
!
!     we have following 16 different structures for radial integrals:
!
! 1 1.0  Rk( . . . .)  ic, jc               -  bound-bound  
! 2 1.1  Rk( . . . .) < i | . > ic          -  bound-channel
! 3 1.2  Rk( . . . .) < i | . > < j | . >   -  channel-channel
! 4 1.3  Rk( . . . .) < i | j >             -  target structure
!
! 5 2.0  Rk( i . . .) < j | . >             -  channel-channel due to
! 6 3.0  Rk( . i . .) < j | . >                overlaps
! 7 4.0  Rk( . . i .) < j | . >
! 8 5.0  Rk( . . . i) < j | . >
!
! 9 2.1  Rk( i . . .)  ic                   -  bound-channel
!10 3.1  Rk( . i . .)  ic
!11 4.1  Rk( . . i .)  ic
!12 5.1  Rk( . . . i)  ic
!
!13 6.0  Rk( i . j .)                       -  direct channel-channel
!14 7.0  Rk( . i . j)
!15 8.0  Rk( i . . j)                       -  exchange channel-channel
!16 9.0  Rk( . i j .)
!
!     where .  denotes bound orbital, i,j - channels.
!
!                                        ibo         ibo
! 1.0 Rk( . . . .)  ic, jc       -  k1=(i1,i3)  k2=(i2,i4)  k3=-ic  k4=-jc 
! 1.1 Rk( . . . .) < i | . > ic  -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4=-ic 
! 1.2 Rk( . . . .) <i|.> <j|.>   -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4= jo
! 1.3 Rk( . . . .) < i | j >     -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4= 0 
!
! 2.0 Rk( i . . .) < j | . >     -  k1=(i2,i4)  k2=i3  k3=ich  k4=io
! 3.0 Rk( . i . .) < j | . >     -  k1=(i1,i3)  k2=i4  k3=ich  k4=io
! 4.0 Rk( . . i .) < j | . >     -  k1=(i2,i4)  k2=i1  k3=ich  k4=io
! 5.0 Rk( . . . i) < j | . >     -  k1=(i1,i3)  k2=i2  k3=ich  k4=io
!
! 2.1 Rk( i . . .)  ic           -  k1=(i2,i4)  k2=i3  k3=ich  k4=-ic
! 3.1 Rk( . i . .)  ic           -  k1=(i1,i3)  k2=i4  k3=ich  k4=-ic           
! 4.1 Rk( . . i .)  ic           -  k1=(i2,i4)  k2=i1  k3=ich  k4=-ic
! 5.1 Rk( . . . i)  ic           -  k1=(i1,i3)  k2=i2  k3=ich  k4=-ic
!
! 6.0 Rk( i . j .)               -  k1=(ich1,ich2)  k2=i2  k3=i4  k4=0
! 7.0 Rk( . i . j)               -  k1=(ich1,ich2)  k2=i1  k3=i3  k4=0
! 8.0 Rk( i . . j)               -  k1=(ich1,ich2)  k2=i2  k3=i3  k4=0
! 9.0 Rk( . i j .)               -  k1=(ich1,ich2)  k2=i1  k3=i4  k4=0
!
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: k,j1,j2,j3,j4,ic,jc,io,jo
      Real(8), intent(in) :: C
      Integer :: k1,k2,k3,k4, i1,i2,i3,i4,m,ich,jch,itype 

      if(abs(C).lt.eps_C) Return

! ... check input parameters:

      if(io.ne.0.and.iabs(io).le.ibo) Stop 'Idef_Rtype: io<ibo'
      if(jo.ne.0.and.iabs(jo).le.ibo) Stop 'Idef_Rtype: jo<ibo'
      if(jc.gt.0.and.ic.le.0) Stop 'jc>0, ic<=0'

! ... apply symmetry:      

      itype=0

      if(ipbs(j1)+ipbs(j2).eq.0 .or. ipbs(j3)+ipbs(j4).eq.0) then
       i1=j1; i3=j3; if(j1.gt.j3) then; i1=j3; i3=j1; end if
       i2=j2; i4=j4; if(j2.gt.j4) then; i2=j4; i4=j2; end if
       if(i1.gt.i2) then
        k1=i1; i1=i2; i2=k1; k3=i3; i3=i4; i4=k3
       end if
      else
       i1=j1; i2=j2; i3=j3; i4=j4
      end if

! ... assign the coefficients:

      if(jc.gt.0) then

       k1=i1*ibo+i3; k2=i2*ibo+i4; k3=-ic; k4=-jc; itype = 1

      elseif(ic.gt.0) then

       if(io.eq.0) then

        if(ipbs(i1).gt.0) then
         k1=i2*ibo+i4; k2=i3; k3=ipbs(i1); k4=-ic; itype = 2
        elseif(ipbs(i2).gt.0) then
         k1=i1*ibo+i3; k2=i4; k3=ipbs(i2); k4=-ic; itype = 2
        elseif(ipbs(i3).gt.0) then            
         k1=i2*ibo+i4; k2=i1; k3=ipbs(i3); k4=-ic; itype = 2
        elseif(ipbs(i4).gt.0) then
         k1=i1*ibo+i3; k2=i2; k3=ipbs(i4); k4=-ic; itype = 2
        end if

       elseif(io.gt.0) then

        k1=i1*ibo+i3; k2=i2*ibo+i4; k3=io; k4=-ic; itype = 1

       end if

      elseif(io.eq.0) then
       if(ipbs(i1).gt.0.and.ipbs(i3).gt.0) then
        if(ipbs(i1).ge.ipbs(i3)) then
         k1=ipbs(i1)*ibo+ipbs(i3); k2=i2; k3=i4; k4=0; itype = 3
         if(k2.gt.k3) then; k2=i4; k3=i2; end if
        else
         k1=ipbs(i3)*ibo+ipbs(i1); k2=i4; k3=i2; k4=0; itype = 3
         if(k2.gt.k3) then; k2=i2; k3=i4; end if
        end if
       elseif(ipbs(i2).gt.0.and.ipbs(i4).gt.0) then
        if(ipbs(i2).ge.ipbs(i4)) then
         k1=ipbs(i2)*ibo+ipbs(i4); k2=i1; k3=i3; k4=0; itype = 3
         if(k2.gt.k3) then; k2=i3; k3=i1; end if
        else
         k1=ipbs(i4)*ibo+ipbs(i2); k2=i3; k3=i1; k4=0; itype = 3
         if(k2.gt.k3) then; k2=i1; k3=i3; end if
        end if
       elseif(ipbs(i1).gt.0.and.ipbs(i4).gt.0) then
        if(ipbs(i1).ge.ipbs(i4)) then
         k1=ipbs(i1)*ibo+ipbs(i4); k2=i3; k3=i2; k4=0; itype = 4
        else
         k1=ipbs(i4)*ibo+ipbs(i1); k2=i2; k3=i3; k4=0; itype = 4
        end if
       elseif(ipbs(i2).gt.0.and.ipbs(i3).gt.0) then
        if(ipbs(i3).ge.ipbs(i2)) then
         k1=ipbs(i3)*ibo+ipbs(i2); k2=i1; k3=i4; k4=0; itype = 4
        else
         k1=ipbs(i2)*ibo+ipbs(i3); k2=i4; k3=i1; k4=0; itype = 4
        end if
       end if                                       

      elseif(io.gt.0.and.jo.eq.0) then

        if(ipbs(i1).gt.0) then
         k1=i2*ibo+i4; k2=i3; k3=ipbs(i1); k4=io; itype = 2
        elseif(ipbs(i2).gt.0) then
         k1=i1*ibo+i3; k2=i4; k3=ipbs(i2); k4=io; itype = 2
        elseif(ipbs(i3).gt.0) then
         k1=i2*ibo+i4; k2=i1; k3=ipbs(i3); k4=io; itype = 2
        elseif(ipbs(i4).gt.0) then
         k1=i1*ibo+i3; k2=i2; k3=ipbs(i4); k4=io; itype = 2
        else
         k1=i1*ibo+i3; k2=i2*ibo+i4; k3=io; k4=0; itype = 1
         if(check_target.eq.0) Return
        end if

      elseif(io.gt.0.and.jo.gt.0) then

       k1=i1*ibo+i3; k2=i2*ibo+i4; k3=io; k4=jo; itype = 1

      end if

      if(itype.eq.0) Stop 'Idef_Rtype: itype=0'

      Call Add_coef(C,k,k1,k2,k3,k4,itype)

      End Subroutine Idef_Rtype


!======================================================================
      Subroutine Idef_itype (k,i1,i2,i3,i4,C,ic,jc,io,jo)
!======================================================================
!
!     we have following 16 different structures for radial integrals:
!
! 1 1.0  Rk( . . . .)  ic, jc               -  bound-bound  
! 2 1.1  Rk( . . . .) < i | . > ic          -  bound-channel
! 3 1.2  Rk( . . . .) < i | . > < j | . >   -  channel-channel
! 4 1.3  Rk( . . . .) < i | j >             -  target structure
!
! 5 2.0  Rk( i . . .) < j | . >             -  channel-channel due to
! 6 3.0  Rk( . i . .) < j | . >                overlaps
! 7 4.0  Rk( . . i .) < j | . >
! 8 5.0  Rk( . . . i) < j | . >
!
! 9 2.1  Rk( i . . .)  ic                   -  bound-channel
!10 3.1  Rk( . i . .)  ic
!11 4.1  Rk( . . i .)  ic
!12 5.1  Rk( . . . i)  ic
!
!13 6.0  Rk( i . j .)                       -  direct channel-channel
!14 7.0  Rk( . i . j)
!15 8.0  Rk( i . . j)                       -  exchange channel-channel
!16 9.0  Rk( . i j .)
!
!     where .  denotes bound orbital, i,j - channels.
!
!                                        ibo         ibo
! 1.0 Rk( . . . .)  ic, jc       -  k1=(i1,i3)  k2=(i2,i4)  k3=-ic  k4=-jc 
! 1.1 Rk( . . . .) < i | . > ic  -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4=-ic 
! 1.2 Rk( . . . .) <i|.> <j|.>   -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4= jo
! 1.3 Rk( . . . .) < i | j >     -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4= 0 
!
! 2.0 Rk( i . . .) < j | . >     -  k1=(i2,i4)  k2=i3  k3=ich  k4=io
! 3.0 Rk( . i . .) < j | . >     -  k1=(i1,i3)  k2=i4  k3=ich  k4=io
! 4.0 Rk( . . i .) < j | . >     -  k1=(i2,i4)  k2=i1  k3=ich  k4=io
! 5.0 Rk( . . . i) < j | . >     -  k1=(i1,i3)  k2=i2  k3=ich  k4=io
!
! 2.1 Rk( i . . .)  ic           -  k1=(i2,i4)  k2=i3  k3=ich  k4=-ic
! 3.1 Rk( . i . .)  ic           -  k1=(i1,i3)  k2=i4  k3=ich  k4=-ic           
! 4.1 Rk( . . i .)  ic           -  k1=(i2,i4)  k2=i1  k3=ich  k4=-ic
! 5.1 Rk( . . . i)  ic           -  k1=(i1,i3)  k2=i2  k3=ich  k4=-ic
!
! 6.0 Rk( i . j .)               -  k1=(ich1,ich2)  k2=i2  k3=i4  k4=0
! 7.0 Rk( . i . j)               -  k1=(ich1,ich2)  k2=i1  k3=i3  k4=0
! 8.0 Rk( i . . j)               -  k1=(ich1,ich2)  k2=i2  k3=i3  k4=0
! 9.0 Rk( . i j .)               -  k1=(ich1,ich2)  k2=i1  k3=i4  k4=0
!
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: k,i1,i2,i3,i4,ic,jc,io,jo
      Real(8), intent(in) :: C
      Integer :: k1,k2,k3,k4,itype

      if(C.eq.0.d0) Return

      if(io.ne.0.and.iabs(io).le.ibo) then
       write(*,*) 'Idef_type: io,ibo = ', io,ibo
       Stop
      end if

      if(jo.ne.0.and.iabs(jo).le.ibo) then
       write(*,*) 'Idef_type: jo,ibo = ', jo,ibo
       Stop
      end if

      itype=0
      if(jc.gt.0.and.ic.le.0) Stop 'jc>0,ic<=0'

      if(jc.gt.0) then

       k1=i1*ibo+i3; k2=i2*ibo+i4; k3=-ic; k4=-jc; itype = 1

      elseif(ic.gt.0) then

       if(io.eq.0) then

        if(ipbs(i1).gt.0) then
         k1=i2*ibo+i4; k2=i3; k3=ipbs(i1); k4=-ic; itype = 2
        elseif(ipbs(i2).gt.0) then
         k1=i1*ibo+i3; k2=i4; k3=ipbs(i2); k4=-ic; itype = 3
        elseif(ipbs(i3).gt.0) then
         k1=i2*ibo+i4; k2=i1; k3=ipbs(i3); k4=-ic; itype = 4
        elseif(ipbs(i4).gt.0) then
         k1=i1*ibo+i3; k2=i2; k3=ipbs(i4); k4=-ic; itype = 5
        end if

       elseif(io.gt.0) then

        k1=i1*ibo+i3; k2=i2*ibo+i4; k3=io; k4=-ic; itype = 1

       end if

      elseif(io.eq.0) then
       if(ipbs(i1).gt.0.and.ipbs(i3).gt.0) then
        k1=ipbs(i1)*ibo+ipbs(i3); k2=i2; k3=i4; k4=0; itype = 6
       elseif(ipbs(i2).gt.0.and.ipbs(i4).gt.0) then
        k1=ipbs(i2)*ibo+ipbs(i4); k2=i1; k3=i3; k4=0; itype = 7
       elseif(ipbs(i1).gt.0.and.ipbs(i4).gt.0) then
        k1=ipbs(i1)*ibo+ipbs(i4); k2=i3; k3=i2; k4=0; itype = 8
       elseif(ipbs(i2).gt.0.and.ipbs(i3).gt.0) then
        k1=ipbs(i3)*ibo+ipbs(i2); k2=i1; k3=i4; k4=0; itype = 9
       end if                                       

      elseif(io.gt.0.and.jo.eq.0) then

        if(ipbs(i1).gt.0) then
         k1=i2*ibo+i4; k2=i3; k3=ipbs(i1); k4=io; itype = 2
        elseif(ipbs(i2).gt.0) then
         k1=i1*ibo+i3; k2=i4; k3=ipbs(i2); k4=io; itype = 3
        elseif(ipbs(i3).gt.0) then
         k1=i2*ibo+i4; k2=i1; k3=ipbs(i3); k4=io; itype = 4
        elseif(ipbs(i4).gt.0) then
         k1=i1*ibo+i3; k2=i2; k3=ipbs(i4); k4=io; itype = 5
        else
         k1=i1*ibo+i3; k2=i2*ibo+i4; k3=io; k4=0; itype = 1
        end if

      elseif(io.gt.0.and.jo.gt.0) then

       k1=i1*ibo+i3; k2=i2*ibo+i4; k3=io; k4=jo; itype = 1

      end if

      if(itype.eq.0) Stop 'IDEF_TYPE: itype=0'
      itype = atype*ibtype + itype  

      Call Add_coef(C,k,k1,k2,k3,k4,itype)

      End Subroutine Idef_itype
