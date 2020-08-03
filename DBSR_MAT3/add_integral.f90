!======================================================================
      Subroutine Add_integral (kpol,i1,i2,i3,i4,C,ic,jc,io,jo)
!======================================================================
!     add integral to the list in module c_data
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer, intent(in) :: kpol,i1,i2,i3,i4,ic,jc,io,jo
      Real(8), intent(in) :: C

      Integer :: m,v
      Real(8) :: S(8), t1,t2 
      Real(8), external :: SMU

      Call CPU_time(t1)

      Select Case(icase)

      Case(0)                           !  Overlaps  

       Call Idef_Otype(C,ic,jc,io,jo)

      Case(1)                           !  L-integral  

       Call Idef_Ltype(i1,i3,C,ic,jc,io,jo)

      Case(2)                           !  R-integral  

       m = mod(lbs(i1)+lbs(i3)+kpol,2) + mod(lbs(i2)+lbs(i4)+kpol,2)
       if(m.eq.0.and.kpol.le.mk)  &
        Call Idef_Rtype(kpol,i1,i2,i3,i4,C,ic,jc,io,jo)

      Case(3)                           !  S-integral

       Do v = kpol-1,kpol+1
        if(v.gt.mk) Cycle
        if(SMU(kbs(i1),kbs(i2),kbs(i3),kbs(i4),kpol,v,S).eq.0) Cycle
        atype = 0
        Call Idef_itype(v,i1,i2,i3,i4,C*S(1),ic,jc,io,jo)
        Call Idef_itype(v,i2,i1,i4,i3,C*S(2),ic,jc,io,jo)
        Call Idef_itype(v,i3,i4,i1,i2,C*S(3),ic,jc,io,jo)
        Call Idef_itype(v,i4,i3,i2,i1,C*S(4),ic,jc,io,jo)
        atype = 1
        Call Idef_itype(v,i1,i2,i3,i4,C*S(5),ic,jc,io,jo)
        Call Idef_itype(v,i4,i3,i2,i1,C*S(6),ic,jc,io,jo)
        Call Idef_itype(v,i3,i4,i1,i2,C*S(7),ic,jc,io,jo)
        Call Idef_itype(v,i2,i1,i4,i3,C*S(8),ic,jc,io,jo)
       End do

      End Select

      Call CPU_time(t2)

      t_add = t_add + (t2-t1)

      End Subroutine Add_integral
