!======================================================================
      Subroutine Check_det(kd,N1,N2,iext)
!======================================================================
!     evaluate or expend the total overlap determinant [n1,n2]^iext 
!     so that to extract one-electron overlaps with continuum orbitals;
!     ipbs - pointer to the continuum orbtal
!----------------------------------------------------------------------
      Use DBS_orbitals_pq, only: ipbs 

      Implicit none
      Integer, intent(in) :: kd, iext, N1(*),N2(*)
      Integer :: N3(kd),N4(kd),N5(kd),N6(kd)
      Integer :: i,j,j1,j2,j3,j4,io,jo,kd1
      Real(8) :: S
      Real(8), external :: VDET
      Integer, external :: IBORT

! ... define the cont.orb.

      j1=0
      Do i=1,kd; if(ipbs(N1(i)).ne.0) then; j1=i; Exit; end if; End do
      j2=0
      Do i=1,kd; if(ipbs(N2(i)).ne.0) then; j2=i; Exit; end if; End do

      if(j1.eq.0.and.j2.eq.0) then               ! no cont.orb

       S = VDET(kd,N1,N2)**iext;  Call Iadd_ndets(0,0,S)  

      elseif(kd.eq.1) then                      ! <kl|nl>

       if(iext.gt.1) Stop ' Check_det: iext > 1'
          
       io=IBORT(N1(1),N2(1));  
       if(io.ne.0) Call Iadd_ndets(io,0,1.d0)

      elseif(j1.gt.0.and.j2.eq.0) then           ! <kl ...|...>

       if(iext.gt.1) Stop 'Check_det: iext > 1 for continuum orbitals'

        Do j=1,kd
         io=IBORT(N1(j1),N2(j)); if(io.eq.0) Cycle
         Call Shift(kd,j1,N1,N3)                     !   move out do loop?
         Call Shift(kd,j ,N2,N4)               
         S = VDET(kd-1,N3,N4)*(-1)**(j1+j)
         Call Iadd_ndets(io,0,S)
        End do

      elseif(j2.gt.0.and.j1.eq.0) then          ! <...|... kl>

        if(iext.gt.1) Stop 'Check_det: iext > 1 for continuum orbitals'

         Do j=1,kd
          io=IBORT(N1(j),N2(j2)); if(io.eq.0) Cycle
          Call Shift(kd,j ,N1,N3)
          Call Shift(kd,j2,N2,N4)                !   move out do loop?
          S = VDET(kd-1,N3,N4)*(-1)**(j2+j)
          Call Iadd_ndets(io,0,S)
         End do

      elseif(j1.gt.0.and.j2.gt.0) then          !  < kl ... | ... kl>

       if(iext.gt.1) Stop 'Check_det: iext > 1 for continuum orbitals'
       
        kd1=kd-1
        Do j=1,kd
       
         Call Shift(kd,j1,N1,N3)   !   move out do loop?
         Call Shift(kd, j,N2,N4)
         io=IBORT(N1(j1),N2(j))
       
         if(j.eq.j2) then
          if(io.eq.0) Stop 'Check_det: <k1|k2>=0 - ?'
          S=VDET(kd1,N3,N4)*(-1)**(j1+j2)
          Call Iadd_ndets(io,0,S)
         else
          j3=j2;if(j.lt.j2) j3=j2-1
          Do j4=1,kd1
           jo=IBORT(N3(j4),N4(j3))
           if(io.eq.0.or.jo.eq.0) Cycle
           Call Shift(kd1,j3,N4,N6)            !   move out do loop?
           Call Shift(kd1,j4,N3,N5)
           S=VDET(kd-2,N5,N6)*(-1)**(j1+j+j3+j4)
           Call Iadd_ndets(io,jo,S)
          End do
         end if
       
        End do

      end if

      End Subroutine Check_det


!======================================================================
      Subroutine Shift (n,m,N1,N2)
!======================================================================
!     array N2 is obtained from N1 by deleting element 'm'
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m, N1(*)
      Integer, intent(out) :: N2(*)
      Integer :: i,k

      k=0; Do i=1,n; if(i.eq.m) Cycle; k=k+1; N2(k)=N1(i); End do

      End Subroutine Shift


