!======================================================================
      Module recup_a
!======================================================================

      Implicit none

! ... working arrays in RECUP:

      Integer, parameter :: mc = 25
      Integer :: Y1(3,mc),Y2(3,mc),JP1(mc),JP2(mc)
      Integer, parameter :: mm = 3*mc-2

! ... number of 6j-sysmbols

      Integer, parameter :: mw = 100         
      Integer, parameter :: mj = mw + mm
      Integer :: JD(2,mj),Jsign(mj),Jfact(mj),JW(6,mw)
      Integer :: nw,nd

! ... working arrays in ZRECUP:

      Integer :: Jsum(mw),IS(mw),Jmax(mw)
      Real(8) :: R(mw),Z(mw)

      Real(8) :: eps_r = 1.d-8

      End Module recup_a


!======================================================================
      Subroutine RECUP(mi,n,J1,J2)
!======================================================================
!     recoupling coefficients between schemes J1 and J2
!     MI - number of moments
!     N  - number of elementary coupling in J1 and J2
!     for details of algorithm see BURKE P.G.,COMP.PHYS.COMM.1(1970)241
!----------------------------------------------------------------------

      Use recup_a

      Implicit real(8) (A-H,O-Z)

      Integer :: J1(3,*),J2(3,*)

      if(mi.gt.mm) Stop 'RECUP: mi > mm, moments'
      if(n .gt.mc) Stop 'RECUP: n  > mc, couplings'

      Do i=1,N           ! rewriting of j1,j2 to y1,y2
      Do j=1,3           ! because these arrays may be modified
       Y1(j,i)=J1(j,i)   ! by subroutine
       Y2(j,i)=J2(j,i)
      End do
      End do

      NW=0               ! the number of 6j-symbols
      ND=0               ! the number of delta-functions c

      Do i=1,mj
       Jsign(i)=0        !  factor (-1)**[j1(i)*Jsing(i)]
       Jfact(i)=0        !  factor (2j1(i)+1)**(Jfact(i)/2)
      End do
!----------------------------------------------------------------------
!                                             search the total moments:
      jt1=J1(3,1)
    1 Do i=2,N
       if(jt1.eq.J1(1,i).or.jt1.eq.J1(2,i)) then
        jt1=J1(3,i)
        go to 1
       end if
      End do
      jt2=J2(3,1)
    2 Do i=2,N
       if(jt2.eq.J2(1,i).or.jt2.eq.J2(2,i)) then
        jt2=J2(3,i)
        go to 2
       end if
      End do
      if(jt1.ne.jt2) then       !  store the delta-function conditions
       nd=nd+1
       JD(1,nd)=jt1
       JD(2,nd)=jt2
       Do i=1,N
        if(J2(3,i).eq.jt2) J2(3,i)=jt1
       End do
      end if
      JT=jt1
!----------------------------------------------------------------------
    9 IP=1         !  number of the same elimentary coupling
   10 IK=0         !  key of result of seaching such coupling

      Do k=IP,N    !  search the identical elimentary couplings
       K1=Y1(1,k)
       K2=Y1(2,k)
      Do l=IP,N
       L1=Y2(1,l)
       L2=Y2(2,l)
       if((K1.eq.L1.and.K2.eq.L2).or.(K1.eq.L2.and.K2.eq.L1)) go to 12
      End do
      End do
      go to 20

   12 IK=1
      K3=Y1(3,k)
      L3=Y2(3,l)
      if(L1.ne.K1) then          !  rearrangement of moments
       Y2(1,l)=L2
       Y2(2,l)=L1
       Jsign(L1)=Jsign(L1)+1     !  fase factor of this rearrangement:
       Jsign(L2)=Jsign(L2)+1     !  (j1+j2)j3 --> (j2+j1)j3 gives
       Jsign(L3)=Jsign(L3)-1     !  (-1) ** (j1+j2-j3)
      end if

      if(L3.ne.K3) then          !  store the delta-function conditions
       nd=nd+1
       JD(1,nd)=K3
       JD(2,nd)=L3
      end if

      Do i=1,N
      Do j=1,3                         ! the same moments -
       if(Y2(j,i).eq.L3) Y2(j,i)=K3    ! the same assingments
      End do
      End do

      Do i=1,3               !  the same elimentary coupling -
       m=Y1(i,IP)            !  the same place in schimes Y1 and Y2
       Y1(i,IP)=Y1(i,k)      !
       Y1(i,k)=m             !  the procedure from label 9 to end
       m=Y2(i,IP)            !  is repeated until Y1 no equal Y2
       Y2(i,IP)=Y2(i,l)
       Y2(i,l)=m
      End do

   20 IP=IP+IK
      IF(IK.eq.1.and.IP.le.N) GO TO 10  ! search next coupling
      IF(IP.eq.N+1) then
	 ! write(*,*) 'nw,nd =', nw,nd
	  Return              ! Y1=Y2 - end of search
      end if
      IF(IP.eq.N) STOP ' ZRECUP: Y1 and Y2 are incompatible, step 1'

!----------------------------------------------------------------------
!     transformation of scheme Y2 to Y1 by elimentary recouplings
!     of 3 moments through 6j-symbols, but firstly program choises
!     the moment for more effective recoupling by following simple
!     algorithm: from various possible pares the program choise those
!     between which the path along tree is shotest.
!----------------------------------------------------------------------
!                     block for choising of pare moment for recoupling:

      kc=0
      kt=0
      nc=3*N
      Do k=IP,N                 ! determination of recoupling moments
       K1=Y1(1,k)
       K2=Y1(2,k)
       i1=0
       i2=0
       Do i=IP,N
        IF(K1.eq.Y2(1,i).or.K1.eq.Y2(2,i)) i1=i
        IF(K2.eq.Y2(1,i).or.K2.eq.Y2(2,i)) i2=i
        IF(i1*i2.gt.0) exit
       End do
      if(i1*i2.eq.0) Cycle

      JP1(1)=K1       ! JP1 and JP2  -  paths from the select moment to
      JP2(1)=K2       !                 the top of tree
      JP1(2)=Y2(3,i1)
      JP2(2)=Y2(3,i2)
      i1=2
      i2=2
      np=0            ! number of searchings
   21 Do i=1,N
       if(JP1(i1).eq.Y2(1,i).or.JP1(i1).eq.Y2(2,i)) then
        i1=i1+1
        JP1(i1)=Y2(3,i)
       end if
       if(JP2(i2).eq.Y2(1,i).or.JP2(i2).eq.Y2(2,i)) then
        i2=i2+1
        JP2(i2)=Y2(3,i)
       end if
      End do
      np=np+1
      if(np.gt.N)
     :  STOP 'ZRECUP: schemes Y2 and Y3 are incompatable, step 3'
      if(JP1(i1).ne.JT.or.JP2(i2).ne.JT) go to 21

      Do i=2,i1                         ! crossing of paths JP1 and JP2
      DO j=2,i2
       if(JP1(i).eq.JP2(j)) go to 22
      End do
      End do

   22 if(i+j.lt.nc) then
       nc=i+j
       kc=k
       kt=i1-i
      elseif(i+j.eq.nc) then
       if(i1-i.gt.kt) then
        nc=i+j
        kc=k
        kt=i1-i
       end if
      end if
      End do        ! over k
      if(kc.eq.0)
     :   STOP 'ZRECUP: schemes Y1 and Y2 are incompatible, step 2'

!----------------------------------------------------------------------
!                          repeat the above block for selected moments
       k=kc
       K1=Y1(1,k)
       K2=Y1(2,k)
       i1=0
       i2=0
       Do i=IP,N
        IF(K1.eq.Y2(1,i).or.K1.eq.Y2(2,i)) i1=i
        IF(K2.eq.Y2(1,i).or.K2.eq.Y2(2,i)) i2=i
        IF(i1*i2.gt.0) exit
       End do

      JP1(1)=K1       ! JP1 and JP2  -  paths from the select moment to
      JP2(1)=K2       !                 the top of tree
      JP1(2)=Y2(3,i1)
      JP2(2)=Y2(3,i2)
      i1=2
      i2=2
   23 Do i=1,N
       if(JP1(i1).eq.Y2(1,i).or.JP1(i1).eq.Y2(2,i)) then
        i1=i1+1
        JP1(i1)=Y2(3,i)
       end if
       if(JP2(i2).eq.Y2(1,i).or.JP2(i2).eq.Y2(2,i)) then
        i2=i2+1
        JP2(i2)=Y2(3,i)
       end if
      End do
      if(JP1(i1).ne.JT.or.JP2(i2).ne.JT) go to 23

      Do i=2,i1                         ! crossing of paths JP1 and JP2
      DO j=2,i2
       if(JP1(i).eq.JP2(j)) go to 24
      End do
      End do
!-----------------------------------------------------------------------
!                                        block of elementary recoupling:
!  determine the coupling moments
   24 M0=JP1(i)                   ! M0 - top of elimentery recoupling
      Do l=1,N                    !
       if(M0.eq.Y2(3,l)) Exit     !
      End do                      !                  M0
!     choise M2,M12 and M3        !                 /  \  <- coupling l
      if(i.ge.j) then             !               M12   M3
       M12=JP1(i-1)               ! coupling m->  /  \  .....
       M2 =JP1(i-2)               !             M1   M2    \
       M3 =JP2(j-1)               !            ........   K1 or K2
      else                        !               /
       M12=JP2(j-1)               !            K1 or K2
       M2 =JP2(j-2)               !
       M3 =JP1(i-1)
      end if
      if(M12.eq.Y2(2,L)) then
       Jsign(M12)=Jsign(M12)+1
       Jsign(M3)=Jsign(M3)+1
       Jsign(M0)=Jsign(M0)-1
      end if

!     coupling m

      Do m=1,N
       if(Y2(3,m).eq.M12) exit
      End do

!     choise M1

      if(M2.eq.Y2(2,m)) then
       M1=Y2(1,M)
      elseif(M2.eq.Y2(1,m)) then
       M1=Y2(2,M)
       Jsign(M1)=Jsign(M1)+1
       Jsign(M2)=Jsign(M2)+1
       Jsign(M12)=Jsign(M12)-1
      else
       Stop ' Recup: problems with M1'
      end if

!----------------------------------------------------------------------
      nw=nw+1
      if(nw.gt.mw)  STOP 'ZRECUP: NW exceeds max.allowed'
      MY=mi+nw

!----------------------------------------------------------------------
!                                                     result couplings:
      Y2(1,M)=M2            !          M0
      Y2(2,M)=M3            !         /  \     <-  coupling L
      Y2(3,M)=MY            !        M1   MY
      Y2(1,L)=M1            !            /  \  <-  coupling M
      Y2(2,L)=MY            !           M2   M3
      Y2(3,L)=M0            !        .....  .....
                            !         K1  or K2
      JW(1,NW)=M1           ! that is, K1 and K2 is now nearer,
      JW(2,NW)=M2           ! and our aim is to place K1 and K2
      JW(3,NW)=M12          ! in the same coupling
      JW(4,NW)=M3           ! ---------------------------------------
      JW(5,NW)=M0           ! the result transformation:
      JW(6,NW)=MY           !
                            !
      Jsign(m1)=Jsign(m1)-1 ! <(M1,M2)M12,M3;M0 || M1,(M2,M3)MY;M0> =
      Jsign(m2)=Jsign(m2)-1 !
      Jsign(m3)=Jsign(m3)-1 !     M1+M2+M3+M0         1/2 {M1 M2 M12}
      Jsign(m0)=Jsign(m0)-1 ! (-1)            [M12,MY]    {M3 M0 MY }
      Jfact(M12)=Jfact(M12)+1
      Jfact(MY)=Jfact(MY)+1
!----------------------------------------------------------------------
      go to 9               ! repeat the above analysis of Y1 and Y2
      End                   ! from the begining


!======================================================================
      Real(8) function ZRECUP(mi,JJ)
!======================================================================
!     calculation of recoupling transformation, obtained by RECUP
!----------------------------------------------------------------------
      Use recup_a

      Implicit real(8) (A-H,O-Z)

      Integer :: JJ(*), J(mj)

      Do i=1,mi
       J(i)=JJ(i)
      End do
      Do i=mi+1,mi+nw
       J(i)=-1
      End do
!----------------------------------------------------------------------
      ZRECUP=0.0               ! checking of DELTA-function conditions
      Do i=1,nd
       IF(J(JD(2,i)).eq.-1) J(JD(2,i))=J(JD(1,i))
       IF(J(JD(1,i)).ne.J(JD(2,i))) go to 5
      End do
!----------------------------------------------------------------------
      ks=0                   ! calculation of sing and multiply factors
      S=1.0
      Do i=1,mi
       ks=ks+(J(i)-1)*Jsign(i)
       if(Jfact(i).gt.0) S=S*DFLOAT(J(i))**Jfact(i)
      End do
      S=sqrt(S)

      if(nw.eq.0) then
       if(mod(ks,2).ne.0) STOP 'ZRECUP: invalid sing factor'
       ZRECUP=S*(-1)**(ks/2)
       go to 5
      end if

      Do i=1,NW
       Jsum(i)=0                                !  sign of summation
       IF(J(JW(6,i)).EQ.-1) Jsum(i)=1
       Z(i)=0.0
       k=1
       Do ik=1,6
        IF(J(JW(ik,i)).EQ.-1) k=0
       End do
!      calculation of the 6j-symbols, which don't include any summation
       if(k.eq.1) then
        Z(i)=Z_6j(J(JW(1,i)),J(JW(2,i)),J(JW(3,i)),
     :            J(JW(4,i)),J(JW(5,i)),J(JW(6,i)))
        if(Z(i).eq.0.0) go to 5
       end if
      End do
!----------------------------------------------------------------------
      RR=0.0            ! calculations of sum of products of 6j-sJmbols
      i=1
      k=i+mi

    1 Continue                ! determination of boundaris of summation
      if(Jsum(i).eq.1) then
       J(k)=MAX0( iabs( J(JW(1,i))-J(JW(5,i)) ),
     :            iabs( J(JW(2,i))-J(JW(4,i)) ))+1
       Jmax(i)=MIN0(    J(JW(1,i))+J(JW(5,i)),
     :                  J(JW(2,i))+J(JW(4,i)) ) -1
      else
       Jmax(i)=0
      end if

    2 Continue                          ! calculation of i-th 6j-sJmbol
      if(Z(i).eq.0.0) then
        A=Z_6j(J(JW(1,i)),J(JW(2,i)),J(JW(3,i)),
     :         J(JW(4,i)),J(JW(5,i)),J(JW(6,i)))
      else
        A=Z(i)
      end if
      A=A*DSQRT(DFLOAT(J(k))**Jfact(k))
      IK=Jsign(k)*(J(k)-1)
      if(i.eq.1) then
       R(1)=A
       IS(1)=IK
      else
       R(i)=A*R(i-1)
       IS(i)=IS(i-1)+IK
      end if
      IF(R(i).eq.0.0) go to 4

    3 Continue
      if(i.lt.NW) then
       i=i+1
       k=mi+i
       go to 1
      else
       m=ks+IS(nw)
       IF(MOD(m,2).ne.0) STOP 'ZRECUP: invalid sing factor'
       RR=RR+R(nw)*(-1)**(m/2)
      end if

    4 if(J(k).lt.Jmax(i)) then
       J(k)=J(k)+2
       go to 2
      elseif(i.gt.1) then
       i=i-1
       k=mi+i
       go to 4
      end if

      ZRECUP=RR*S
      if(abs(ZRECUP).lt.Eps_r) ZRECUP=0.d0
                                     !       write(*,*) 'ZRECUP=',ZRECUP
    5 Return
      End
