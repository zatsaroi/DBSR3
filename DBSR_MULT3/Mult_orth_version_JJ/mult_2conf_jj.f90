!=======================================================================
      Subroutine Mult_2conf_jj(no1,nn1,ln1,jn1,iq1,Jshell1,Vshell1,Jintra1, &
                               no2,nn2,ln2,jn2,iq2,Jshell2,Vshell2,Jintra2, &
                               atype, nk, ck,ik,jk)
!=======================================================================
!     compute the angular coefficients for multipole operator "kpol"
!     between 2 atomic states in case of orthogonal orbitals.
!     In this case, we are expecting only 1 radial integral:
!            ck  * INT[ P_ik * r^k * P_jk, r=0,inf]
!-----------------------------------------------------------------------
      Implicit none 

! ... input-output:

      Integer, intent(in)  :: no1,nn1(no1),ln1(no1),jn1(no1),iq1(no1), &
                              Jshell1(no1),Vshell1(no1),Jintra1(no1)
      Integer, intent(in)  :: no2,nn2(no2),ln2(no2),jn2(no2),iq2(no2), &
                              Jshell2(no2),Vshell2(no2),Jintra2(no2)
      Character(2), intent(in) :: atype
      Integer   :: kpol
      Character :: ktype

      Integer, intent(out) :: nk, ik(*),jk(*)
      Real(8), intent(out) :: ck(*)

! ... determinant expansion:

      Integer :: kdt1, kdt2
      Real(8), allocatable :: Cdet1(:), Cdet2(:)
      Integer, allocatable :: MJdet1(:,:), MJdet2(:,:)   
      Real(8) :: CC_det, eps_C = 1.d-6

! ... local variables:

      Integer :: ne, kd1, kd2, i, k, mdt1, mdt2, i1,i2,i3, j1,j2,  &
                 JT1, JT2, MT1, MT2, qpol, kkpol
      Integer, allocatable ::  ip1(:), ip2(:)     
      Integer, external :: Iglq, Ndets_jq
      Real(8) :: CN
      Real(8), allocatable :: detnl(:,:)
      Real(8) , external :: Z_3j2   

      read(atype,'(a1,i1)') ktype, kpol

      nk = 0 

! ... check the selection rules for electric multipole transition:

      i1=(-1) ** SUM(ln1(1:no1)*iq1(1:no1)) 
      i2=(-1) ** SUM(ln2(1:no2)*iq2(1:no2)) 

      if(ktype.eq.'E') then
       if(i1.eq.i2.and.mod(kpol,2).ne.0) Return
       if(i1.ne.i2.and.mod(kpol,2).ne.1) Return     
      else
       if(i1.eq.i2.and.mod(kpol,2).ne.1) Return
       if(i1.ne.i2.and.mod(kpol,2).ne.0) Return     
      end if

      JT1 = Jintra1(no1);  MT1 = JT1
      JT2 = Jintra2(no2);  MT2 = JT2
      i1 = JT1/2;  i2 = JT2/2;  i3=kpol
      if(i1.gt.i2+i3.or.i2.gt.i1+i3.or.i3.gt.i1+i2) Return
      if(i1.lt.iabs(i2-i3).or.i2.lt.iabs(i1-i3).or. &
         i3.lt.iabs(i1-i2)) Return

! ... initialize arrays:

      ne = SUM(iq1(1:no1))
      if(ne.ne.SUM(iq2(1:no2))) Stop 'Coef_ee_2conf: ne1 <> ne2' 
      if(allocated(ip1)) Deallocate(ip1,ip2); Allocate(ip1(ne), ip2(ne))

      k=1; Do i=1,no1; ip1(k:k+iq1(i)-1)=i; k=k+iq1(i); End do
      k=1; Do i=1,no2; ip2(k:k+iq2(i)-1)=i; k=k+iq2(i); End do

      if(allocated(detnl)) Deallocate(detnl); Allocate(detnl(ne,ne))

! ... determinant expansion 1:

      mdt1=1;  Do i=1,no1;   mdt1=mdt1*Ndets_jq(jn1(i),iq1(i)); End do
      if(allocated(Cdet1) ) Deallocate(Cdet1 );  Allocate(Cdet1(mdt1))
      if(allocated(MJdet1)) Deallocate(MJdet1);  Allocate(MJdet1(ne,mdt1))

      Call Det_exp_1conf(no1,jn1,iq1,Jshell1,Vshell1,Jintra1,ne, &
                         mdt1,kdt1,Cdet1,MJdet1)  

! ... determinant expansion 2:

      mdt2=1;  Do i=1,no2;   mdt2=mdt2*Ndets_jq(jn2(i),iq2(i)); End do
      if(allocated(Cdet2) ) Deallocate(Cdet2 );  Allocate(Cdet2(mdt2))
      if(allocated(MJdet2)) Deallocate(MJdet2);  Allocate(MJdet2(ne,mdt2))

      Call Det_exp_1conf(no2,jn2,iq2,Jshell2,Vshell2,Jintra2,ne, &
                         mdt2,kdt2,Cdet2,MJdet2)  

! ... calculations:     

      qpol = MT1-MT2
      CN = Z_3j2(JT1,-JT1,kpol+kpol,qpol,JT2,JT2)         
      if(kpol.eq.0) CN = CN * sqrt(1.d0 + JT1)              !   ????
      if(CN.eq.0.d0) Return

      Do kd1 = 1,kdt1
      Do kd2 = 1,kdt2
       CC_det = Cdet1(kd1) * Cdet2(kd2);   Call me_det
      End do;  End do

      if(nk.gt.0) then 

       CK(1:nk) = CK(1:nk) / CN

       k = 0
       Do i = 1,nk
        if(abs(CK(i)).lt.eps_C) Cycle
        k = k + 1; CK(k) = CK(i); ik(k)=ik(i); jk(k)=jk(i)
       End do
       nk=k

      end if


 Contains

!======================================================================
      Subroutine me_det
!======================================================================
!     find the possible interaction orbitals in two determinants and
!     call the subroutine to calculate m.e. between possible 
!     combinations of nj-orbitals 
!----------------------------------------------------------------------
      Implicit none
      Integer :: i,j, i1,i2, m, m1,m2, idif, jdif, ii,jj, k1,k2, kz
      Real(8) :: C
      Real(8), external :: me_jj,det
      
! ... find interaction orbitals:

      Do i1=1,ne; j1 = ip1(i1)
       Do i2=1,ne; j2 = ip2(i2) 
        detnl(i1,i2) = 0.d0
        if(nn1(j1).ne.nn2(j2)) Cycle
        if(ln1(j1).ne.ln2(j2)) Cycle
        if(jn1(j1).ne.jn2(j2)) Cycle
        if(MJdet1(i1,kd1).ne.MJdet2(i2,kd2)) Cycle
        detnl(i1,i2) = 1.d0
       End do
      End do

      idif=0
      Do i1=1,ne
       if(sum(detnl(i1,:)).ne.0.d0) Cycle;  idif = idif + 1;  i = i1
      End do

      jdif=0
      Do i2=1,ne
       if(sum(detnl(:,i2)).ne.0.d0) Cycle;  jdif = jdif + 1;  j = i2
      End do

      if(idif.ne.jdif) Stop 'me_det: problems with determinants'
      if(idif.gt.1) Return
!      if(idif.eq.0.and.kpol.eq.0) Return


      if(idif.eq.1) then

       i1 = ip1(i);  i2 = ip2(j)
       C =  CC_det * me_jj(ktype,kpol,qpol,ln1(i1),jn1(i1),MJdet1(i,kd1), &
                                           ln2(i2),jn2(i2),MJdet2(j,kd2)) 
       if(abs(C).eq.0.d0) Return
       detnl(i,j) = 1.d0
       C = C * NINT(det(ne,detnl))

       if(nk.eq.0) then
        nk = 1;  ck(1)=C; ik(1)=i1; jk(1)=i2 
       else
        k = 0
        Do i = 1, nk
         if(i1.ne.ik(i)) Cycle
         if(i2.ne.jk(i)) Cycle
         ck(i) = ck(i) + C
         k = 1; Exit
        End do
        if(k.eq.0) then
         nk = nk+1; ck(nk)=C; ik(nk)=i1; jk(nk)=i2
        end if
       end if

      else

       kz = NINT(det(ne,detnl))
       Do i = 1,ne; i1 = ip1(i) 
        Do j = 1,ne; if(detnl(i,j).eq.0.d0) Cycle; Exit; End do; i2 = ip2(j)

        C =  CC_det * me_jj(ktype,kpol,qpol,ln1(i1),jn1(i1),MJdet1(i,kd1), &
                                            ln2(i2),jn2(i2),MJdet2(j,kd2)) 
        C = C * kz

        if(nk.eq.0) then
         nk = 1; ck(1)=C; ik(1)=i1; jk(1)=i2 
        else
         k = 0
         Do ii = 1, nk             
          if(i1.ne.ik(ii)) Cycle
          if(i2.ne.jk(ii)) Cycle
          ck(ii) = ck(ii) + C
          k = 1; Exit
         End do
         if(k.eq.0) then
          nk = nk+1; ck(nk)=C; ik(nk)=i1; jk(nk)=i2
         end if
        end if

       End do

      end if

      End Subroutine me_det

      End Subroutine Mult_2conf_jj


!====================================================================
      Real(8) Function me_jj(ktype,kpol,qpol,l1,j1,m1,l2,j2,m2)  result(c)
!====================================================================
!     angular part of electric or magnetic transition operator 
!     between 'nljm' orbitals:
!
!            <n1,l1,j1,m1| T(kq) | n2,l2,j2,m2>
!
!--------------------------------------------------------------------
      Implicit none

      Integer, intent(in) :: l1,j1,m1,l2,j2,m2,kpol,qpol
      Integer :: i
      Character :: ktype
      Real(8), external :: Z_3j2, Cjkj
      
      i=mod(l1+l2+kpol,2)
      if(ktype.eq.'E'.and.i.eq.1) Return                                                                                                  
      if(ktype.eq.'M'.and.i.eq.0) Return 

      C = (-1)**((j1-m1)/2) * Z_3j2(j1,-m1,kpol+kpol,qpol,j2,m2) 
      if(C.eq.0.d0) Return 

      C = C * Cjkj (j1,kpol,j2)

      End Function me_jj



!======================================================================
      Subroutine det_exp_1conf(no,jn,iq,Jshell,Vshell,Jintra, &
                               ne,mdt,kdt,C_det,IP_det) 
!======================================================================
!     defines the determinant expansion for 1 input configuration
!-----------------------------------------------------------------------
      Implicit none 

      Integer, intent(in)  :: no,jn(no),iq(no),Jshell(no),Vshell(no), &
                              Jintra(no),ne,mdt 
      Integer, intent(out) :: kdt, IP_det(ne,mdt)
      Real(8), intent(out) :: C_det(mdt)

      Integer :: i,j,k, JW,JQ, mj_max, ipn(no),ips(no),md(no),nd(no), &
                 MJs(no), MJi(no), idet(ne)
      Real(8) :: C
      Integer, external :: Ndets_jq, Jterm, mj_value 
      Real(8), external :: DETC_jq, Clebsh2

      k = 1; mj_max = 0
      Do i=1,no
       ipn(i)=k; k=k+iq(i); md(i)=Ndets_jq(jn(i),iq(i))
       ips(i) = Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ)
       if(mj_max.lt.jn(i)) mj_max=jn(i)
      End do

!---------------------------------------------------------------------

      kdt=0; i=1; nd(i)=1              
    1 Call DET_sh_jq(jn(i),iq(i),nd(i),MJs(i),Idet(ipn(i)))

      if(i.eq.1) then
       MJi(1)=MJs(1)
      else
       MJi(i) = MJi(i-1)+MJs(i)
      end if

      if(i.lt.no) then;  i=i+1;  nd(i)=1;  go to 1; end if

      C = 0.d0

      if(MJi(no).ne.Jintra(no)) go to 2          !  M_total = J_total

      C = DETC_jq(jn(1),iq(1),ips(1),nd(1))
      if(C.eq.0.d0) go to 2
      
      Do j=2,no
       C = C * DETC_jq(jn(j),iq(j),ips(j),nd(j))
       if(C.eq.0.d0) Exit
       C = C * Clebsh2(Jintra(j-1),MJi(j-1), &
                       Jshell(j  ),MJs(j  ), &
                       Jintra(j  ),MJi(j  ))
       if(C.eq.0.d0) Exit
      End do

      if(C.ne.0) then
       kdt=kdt+1; IP_det(:,kdt)=Idet(1:ne); C_det(kdt)=C
      end if

    2 nd(i)=nd(i)+1                ! selecting the next case

      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          ! to end
       i=i-1; go to 2
      end if
      go to 1

    3 Continue
!---------------------------------------------------------------------

      Do j = 1,kdt
       Do i = 1,ne
        k = IP_det(i,j)
        if(mod(k,2).eq.0) then
         IP_det(i,j) = k-1
        else
         IP_det(i,j) = -k
        end if
       End do   
      End do

      End Subroutine Det_exp_1conf 


!--------------------------------------------------------------------
      Real(8) FUNCTION Z_3j (j1,m1,j2,m2,j3,m3) 
!--------------------------------------------------------------------
!
!     determines the value of the 3j-symbols without direct using of
!     factorials. The following expression for the 3j-symbols is used:
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977)
!
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} *
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ]
!                         SUM(z) {   (-1)^z  /
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! *
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] }
!
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ]
!
!     If we introduce the auxiliary values a(i) and b(i)
!     (see below the text of program) then
!
!     3j =         (-1) ^ Sum[a(i)]
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] }
!                  Sum(z) { (-1)^z  /
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] }
!
!     (below the moments are used in (2J+1)-representation)
!
!--------------------------------------------------------------------

      Implicit none
 
      Integer(4), intent(in) :: j1,m1,j2,m2,j3,m3

      Integer(4) :: i,i_max,k,kk,m,iz,iz_min,iz_max
      Real(8) :: x,y,z

      Integer(4) a(3),b(3),J(16) 

      Z_3j=0.0

      IF(M1+M2+M3-3.ne.0) RETURN    ! check of conservation rules
      J(1)= J1+J2-J3-1
      J(2)= J1-J2+J3-1
      J(3)= J2-J1+J3-1
      J(4)= J1+M1-2
      J(5)= J1-M1
      J(6)= J2-M2
      J(7)= J2+M2-2
      J(8)= J3+M3-2
      J(9)= J3-M3
      Do I=1,9
       IF(J(i).lt.0.or.mod(J(i),2).eq.1) RETURN
      End do

      a(1) = 0                         ! auxiliary values
      a(2) = (j2-j3-m1+1)/2
      a(3) = (j1-j3+m2-1)/2
      b(1) = (j1+j2-j3-1)/2
      b(2) = (j1-m1)/2
      b(3) = (j2+m2-2)/2

      IZ_min=MAX0(a(1),a(2),a(3))      ! limits of the sum
      IZ_max=MIN0(b(1),b(2),b(3))
      IF(IZ_max.LT.IZ_min) Return

      Do I=1,3                         ! constant factorial parameters
      Do K=1,3
       J(I+3*K-3)=b(i)-a(k)
      End do
      End do
      J(10)=(j1+j2+j3-3)/2+1

      Do I=1,3
       J(I+10)=IZ_min-a(i)               ! initial factorial parameters
       J(I+13)=b(i)-IZ_min               ! in the sum
      End do

      Z=0.0
      DO IZ=IZ_min,IZ_max                 ! summation

       I_max=0                            ! max. factorial
       Do I=1,16
        if(J(i).gt.I_max) I_max=J(i)
       End do

       Y=1.0
       DO I=2,I_max         ! estimation of one term in sum
        K=0                 ! K - the extent of the integer I in term
        DO M=1,9
         IF(J(M).GE.I) K=K+1
        End do
        IF(J(10).GE.I) K=K-1
        DO M=11,16
         IF(J(M).GE.I) K=K-2
        End do
        IF(K.EQ.0) Cycle

        X=DBLE(I)                   ! Y = Y * I ** K/2
        KK=IABS(K)/2
        IF(KK.GT.0) THEN
         DO M=1,KK
          IF(K.GT.0) Y=Y*X
          IF(K.LT.0) Y=Y/X
         END DO
        END IF
        IF(mod(K,2).EQ.+1) Y=Y*SQRT(X)
        IF(mod(K,2).EQ.-1) Y=Y/SQRT(X)
       End do

       IF(mod(IZ,2).eq.1) Y=-Y
       Z=Z+Y

       Do I=11,13                  ! new factorial parameters in sum
        J(I)=J(I)+1
       End do
       DO I=14,16
        J(I)=J(I)-1
       End do

      End do                       ! end of summation

      K=a(1)+a(2)+a(3)
      if(mod(k,2).ne.0) Z=-Z
      Z_3j=Z

      END FUNCTION Z_3j


!--------------------------------------------------------------------
      Real(8) FUNCTION Z_3jj(j1,m1,j2,m2,j3,m3)
!--------------------------------------------------------------------

      IMPLICIT NONE
  
      Integer(4), intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), External :: Z_3j
       
      Z_3jj=Z_3j(j1+j1+1,m1+m1+1,j2+j2+1,m2+m2+1,j3+j3+1,m3+m3+1)

      End FUNCTION Z_3jj

!--------------------------------------------------------------------
      Real(8) FUNCTION Z_3j2(j1,m1,j2,m2,j3,m3)
!--------------------------------------------------------------------

      IMPLICIT NONE
  
      Integer(4), intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), External :: Z_3j
       
      Z_3j2=Z_3j(j1+1,m1+1,j2+1,m2+1,j3+1,m3+1)

      End FUNCTION Z_3j2


!====================================================================
      Real(8) FUNCTION CLEBSH(J1,M1,J2,M2,J,M)
!====================================================================
!
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are used in (2J+1)-representation)
!
!     Call:  Z_3j
!--------------------------------------------------------------------

      Implicit None
      
      Integer, intent(in) :: J, M, J1, M1, J2, M2
      Real(8), External :: Z_3j

      Clebsh=(-1)**((j1-j2+m-1)/2)*sqrt(DBLE(J))*   &
             Z_3j(j1,m1,j2,m2,J,-m+2)

      END FUNCTION CLEBSH


!======================================================================
      Real(8) FUNCTION CLEBCH(L1,M1,L2,M2,L,M)
!======================================================================

      Implicit None
      
      Integer, intent(in) :: L, M, L1, M1, L2, M2 
      Real(8), External :: Clebsh

      Clebch = Clebsh(l1+l1+1,m1+m1+1,l2+l2+1,m2+m2+1,l+l+1,m+m+1)

      End FUNCTION CLEBCH


!======================================================================
      Real(8) FUNCTION CLEBSH2(L1,M1,L2,M2,L,M)
!======================================================================

      Implicit None
      
      Integer, intent(in) :: L, M, L1, M1, L2, M2 
      Real(8), External :: Clebsh

      Clebsh2 = Clebsh(l1+1,m1+1,l2+1,m2+1,l+1,m+1)

      End FUNCTION CLEBSH2



!=======================================================================
      Real(8) function Cjkj (j1,k,j2)
!=======================================================================
!     Computes the relativistic reduced matrix element                    
!     (see Eq.15, Grant and Pyper, J.Phys.B9,761,1976)
!                                                                     
!     (j1 || C^k  || j2)  = 
!                                       (j1   k  j2 )
!         (-1)^(j1+1/2) sqrt([j1][j2])      
!                                       (1/2  0 -1/2)                       
!     C(j,0,j) = sqrt([j])
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,k,j2
      Real(8), external :: Z_3j2
      Cjkj = Z_3j2(j1,1,k+k,0,j2,-1) * &
            sqrt(real((j1+1)*(j2+1))) * (-1)**((j1+1)/2)
      End function Cjkj 

!====================================================================
      Integer Function ndets_jq(j,q)
!====================================================================
!     number of det.s in subshells j^k (Newton's binom)
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j,q
      Integer :: i
      Real(8) :: S
      if(q.gt.j+1) Stop 'ndets_jq:  q > q_max'
      S=1.d0
      Do i=q+1,j+1; S=S*i/(i-q); End do
      ndets_jq = S + 0.1d0
      End Function ndets_jq


!====================================================================
      Integer Function mj_value(i)
!====================================================================
!     mj value for orbital 'i' in the lit: -1,+1,-3,+3,-5,+5, ...
!     mj -> in 2j-representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i
      mj_value = -i
      if(mod(i,2).eq.0) mj_value = i - 1
      END Function mj_value


!======================================================================
      Integer Function Jterm (j,q,k,JT,JV,JW,JQ)
!======================================================================
!     provides information about possible terms in the j^q subshell 
!     for  j <= 9/2
!
!     Options:
!     k > 0  --> JT,JV,JQ,JM = k-th term of j^q-subshell
!     k = 0  --> Jterm = position of the (JT,JV) term in subshell list
!     k < 0  --> Jterm = number of terms in j^q-subshell
!
!     j, JT, JQ  -->  in 2*J representation
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j,q,k
      Integer :: JT,JV,JQ,JW
      Integer :: qm,nterm,ip,i,ii,qq
      Integer :: term_list(3,65)

      Data term_list                                               &
       /1, 5, 2,  3, 3, 0,  3, 9, 0,                               & ! 5/2 ^ 3
        1, 7, 3,  3, 3, 1,  3, 5, 1,  3, 9, 1,  3,11, 1,  3,15, 1, & ! 7/2 ^ 3
        0, 0, 4,  2, 4, 2,  2, 8, 2,  2,12, 2,  4, 4, 0,  4, 8, 0, & ! 7/2 ^ 4
        4,10, 0,  4,16, 0,                                         & 
        1, 9, 4,  3, 3, 2,  3, 5, 2,  3, 7, 2,  3, 9, 2,  3,11, 2, & ! 9/2 ^ 3
        3,13, 2,  3,15, 2,  3,17, 2,  3,21, 2,                     &
        0, 0, 5,  2, 4, 3,  2, 8, 3,  2,12, 3,  2,16, 3,  4, 0, 1, & ! 9/2 ^ 4
        4, 4, 1,  4, 6, 1,  4, 8,10,  4, 8,20,  2,10, 1,  4,12,10, &
        4,12,20,  4,14, 1,  4,16, 1,  4,18, 1,  4,20, 1,  4,24, 1, &
        1, 9, 4,  3, 3, 2,  3, 5, 2,  3, 7, 2,  3, 9, 2,  3,11, 2, & ! 9/2 ^ 5
        3,13, 2,  3,15, 2,  3,17, 2,  3,21, 2,  5, 1, 0,  5, 5, 0, &
        5, 7, 0,  5, 9, 0,  5,11, 0,  5,13, 0,  5,15, 0,  5,17, 0, &
        5,19, 0,  5,25, 0/

!       JW = JQ/10  -> additional seniority
!----------------------------------------------------------------------
! ... check input data:

      if(j.lt.0.or.mod(j,2).eq.0) then
       write(*,*) 'Jterm: unphysical j-value: j =',j
       Stop  'Stop in Jterm'
      end if 

      qm = j + 1
      if(q.lt.0.or.q.gt.qm) then
       write(*,*) 'Jterm: number of electron is out of range: q =',q
       Stop  'Stop in Jterm'
      end if 

      if(j.gt.9.and.q.gt.2) then
       write(*,*) 'Jterm: j^q out of scope: j,q =',j,q
       Stop  'Stop in Jterm'
      end if 

      Jterm = 0 

!----------------------------------------------------------------------
! ... the number of terms in j^q subshell:
      
      if(q.le.1.or.q.ge.qm-1) then
       nterm = 1 
      elseif(q.eq.2.or.q.eq.qm-2) then
       nterm = qm/2
      else
       Select case(j*100 + q)
        Case(503);       nterm = 3; ip = 0    ! 5/2 ^ 3
        Case(703,705);   nterm = 6; ip = 3    ! 7/2 ^ 3,5
        Case(704);       nterm = 8; ip = 9    ! 7/2 ^ 4
        Case(903,907);   nterm =10; ip =17    ! 9/2 ^ 3
        Case(904,906);   nterm =18; ip =27    ! 9/2 ^ 4,6
        Case(905);       nterm =20; ip =45    ! 9/2 ^ 5
        Case default;
         write(*,*) 'Jterm: cannot find j^q subshell for j,q=',j,q 
         Stop 'Stop in Jterm'
       End Select
      end if

      if(k.lt.0) then; Jterm = nterm; Return; end if
!----------------------------------------------------------------------
! ... position of the JT,JV term in the subshell list      

      if(k.eq.0) then          

      qq = min0(q,qm-q)
      Select case(qq)
       case(0);  JV=0
       case(1);  JV=1
       case(2);  JV=2
       case(3);  JV=3; if(j.eq.JT) JV=1
       case(4);  if(j.eq.7) then
                  if(JT.eq.12) JV=2
                  if(JT.eq.10) JV=4
                  if(JT.eq.16) JV=4
                 end if
                 if(j.eq.9.and.JT.ge.18) JV=4
       case(5);  if(j.eq.9) then
                  if(JT.eq. 1) JV=5
                  if(JT.eq. 3) JV=3
                  if(JT.eq.19) JV=5
                  if(JT.eq.21) JV=3
                  if(JT.eq.25) JV=5
                 end if
      End select
      if(JT.eq.0.and.j.lt.9) JV=0

       if(q.le.1.or.q.ge.j) then
        Jterm =  1
       elseif(q.eq.2.or.q.eq.qm-2) then
        Jterm =  JT/4+1
       else
        Do i=1,nterm; ii = ip + i
         if(JV.ne.term_list(1,ii).or.JT.ne.term_list(2,ii)) Cycle
         Jterm=i; Exit
        End do
      end if
 
      if(Jterm.eq.0) then
       write(*,*) 'Jterm:  incorect term j^q(JV,JT):',j,q,JV,JT
       write(*,*) '2j =',j
       write(*,*) 'q  =',q
       write(*,*) '2J =',JT
       write(*,*) 'v  =',JV
       Stop 'Stop in Jterm'
      end if

      Return
      end if

!----------------------------------------------------------------------
! ... k-th term of j^q subshell

      if(k.gt.nterm) then 
       write(*,*) 'Jterm:  incorect input k =',k
       Stop 'Stop in Jterm'
      end if

      JW = 0
      if(q.eq.0.or.q.eq.qm) then                 !  j ^ 0
       JV = 0; JT = 0; JQ = qm/2
      elseif(q.eq.1.or.q.eq.qm-1) then           !  j ^ 1
       JV = 1; JT = j; JQ = qm/2-1
      elseif(q.eq.2.or.q.eq.qm-2) then           !  j ^ 2
       JV = 0; if(k.gt.1) JV = 2 
       JT = (k-1) * 4
       JQ = qm/2; if(k.gt.1) JQ = JQ - 2
      else                                       !  j ^ q
       ii = ip + k
       JV = term_list(1,ii)
       JT = term_list(2,ii)
       JQ = term_list(3,ii)
       if(JQ.ge.10) then; JW = JQ/10; JQ = 1; end if
      end if

      End Function Jterm


!---------------------------------------------------------------------
      Real(8) Function DET(N,A)
!---------------------------------------------------------------------
!     determinant of array A(N,N)   (Gauss method)
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: N
      Real(8), intent(inout) :: A(N,N)
      Integer :: I,J,K
      Real(8) :: MAX, T

      DET = 1.d0
      
      DO K=1,N

       MAX=0.d0
       DO I=K,N
        T=A(I,K)
        if(ABS(T).gt.ABS(MAX)) then
         MAX=T; J=I
        end if
       END DO

       IF(MAX.EQ.0.d0) THEN; DET=0.d0; Return; END IF
      
       IF(J.NE.K) THEN
        DET = -DET
        DO I=K,N
         T=A(J,I); A(J,I)=A(K,I); A(K,I)=T
        END DO
       END IF
  
       IF(K+1.LE.N) THEN
        DO I=K+1,N
         T=A(I,K)/MAX
         DO J=K+1,N
          A(I,J)=A(I,J)-T*A(K,J)
         END DO
        END DO
       END IF
    
       DET=DET*A(K,K)
    
      END DO
    
      End Function DET
