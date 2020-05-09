!=============================================================================
      Subroutine coef_2conf_jj(no1,nn1,ln1,jn1,iq1,Jshell1,Vshell1,Jintra1,  &
                               no2,nn2,ln2,jn2,iq2,Jshell2,Vshell2,Jintra2,  &
                               mcoef,ncoef,icoefs,coefs)
!=============================================================================
!     compute the angular coefficients for 2 different atomic states
!-----------------------------------------------------------------------------
      Implicit none 

! ... input-output:

      Integer, intent(in) :: no1,nn1(no1),ln1(no1),jn1(no1),iq1(no1), &
                             Jshell1(no1),Vshell1(no1),Jintra1(no1),  &
                             no2,nn2(no2),ln2(no2),jn2(no2),iq2(no2), &
                             Jshell2(no2),Vshell2(no2),Jintra2(no2),  &
                             mcoef
      Integer ::  ncoef,icoefs(5,mcoef)
      Real(8) ::  coefs(mcoef)
      Real(8), allocatable :: coef(:,:,:,:,:)

! ... determinant expansion:

      Integer :: kdt, kdt1, kdt2 
      Integer, allocatable :: IP_det(:,:), IP_det1(:,:), IP_det2(:,:)
      Real(8), allocatable :: C_det(:), C_det1(:), C_det2(:)
      Real(8) :: CC_det, eps_C = 1.d-7

! ... local variables:

      Integer              :: ne, Jtotal, i,j,k, kmax,kd1,kd2, i1,i2,j1,j2
      Integer, allocatable :: ip1(:),ip2(:)
      Integer, external    :: mj_value

! ... initialize arrays:

      ncoef = 0
      ne = SUM(iq1(1:no1));  if(ne.ne.SUM(iq2(1:no2))) Return
      Jtotal = Jintra1(no1); if(Jtotal.ne.Jintra2(no2)) Return

! ... determinant expansions:

      Call Det_expn_jj (no1,ln1,jn1,iq1,Jshell1,Vshell1,Jintra1)  
      kdt1=kdt

      Allocate(C_det1(kdt1), IP_det1(ne,kdt1))
      C_det1(1:kdt1) = C_det(1:kdt1)
      IP_det1(1:ne,1:kdt1) = IP_det(1:ne,1:kdt1)


      Call Det_expn_jj (no2,ln2,jn2,iq2,Jshell2,Vshell2,Jintra2)  
      kdt2=kdt

      Allocate(C_det2(kdt2), IP_det2(ne,kdt2))
      C_det2(1:kdt2) = C_det(1:kdt2)
      IP_det2(1:ne,1:kdt2) = IP_det(1:ne,1:kdt2)

      Deallocate(C_det,IP_det)

      if(allocated(ip1)) Deallocate(ip1); Allocate(ip1(ne))
      if(allocated(ip2)) Deallocate(ip2); Allocate(ip2(ne))

      k=1; Do i=1,no1; ip1(k:k+iq1(i)-1)=i; k=k+iq1(i); End do
      k=1; Do i=1,no2; ip2(k:k+iq2(i)-1)=i; k=k+iq2(i); End do

! ... calculations:                                  

      kmax = (maxval(jn1(1:no1)) + maxval(jn2(1:no2)))/2

      if(allocated(coef)) Deallocate(coef)
      Allocate(coef(no1,no1,no2,no2,-1:kmax)) 
      coef = 0.d0

      Do kd1 = 1,kdt1
      Do kd2 = 1,kdt2
        CC_det = C_det1(kd1) * C_det2(kd2);  Call me_det
      End do;  End do

! ... final 

      ncoef = 0
      Do i1=1,no1; Do i2=1,no1
      Do j1=1,no2; Do j2=1,no2
      Do k=0,kmax 

       if(abs(coef(i1,i2,j1,j2,k)).lt.eps_c) Cycle

       ncoef=ncoef+1
       if(ncoef.gt.mcoef) Stop 'Coef_ee_2conf: ncoef > mcoef'
       coefs(ncoef)=coef(i1,i2,j1,j2,k)
       icoefs(1,ncoef)=k
       icoefs(2,ncoef)=i1
       icoefs(3,ncoef)=i2
       icoefs(4,ncoef)=j1
       icoefs(5,ncoef)=j2

      End do 
      End do; End do 
      End do; End do 

! ... one-electron integrals:

      Do i=1,no1; Do j=1,no2
       if(abs(coef(i,i,j,j,-1)).lt.eps_c) Cycle
       ncoef=ncoef+1
       if(ncoef.gt.mcoef) Stop 'Coef_ee_2conf: ncoef > mcoef'
       coefs(ncoef)=coef(i,i,j,j,-1)
       icoefs(1,ncoef)=-1
       icoefs(2,ncoef)=i
       icoefs(3,ncoef)=i
       icoefs(4,ncoef)=j
       icoefs(5,ncoef)=j
      End do; End do 

CONTAINS

!======================================================================
      Subroutine Det_expn_jj(no, ln,jn,iq,Jshell,Vshell,Jintra)  
!======================================================================
!     procedure of exaustion of all possible determinants for given
!     configurations. The determinants and their coefficients
!     are recoded on unit 'nua'
!
!     Calls: Det_sh_jq, DETC_jq, Clebsh2, Ndets_jq, Jterm, mj_value 
!----------------------------------------------------------------------
      Implicit none 

      Integer, intent(in) :: no,ln(no),jn(no),iq(no), &
                             Jshell(no),Vshell(no),Jintra(no)

      Integer :: i,j,k, mkdt, JW,JQ
      Integer, external :: Ndets_jq, Jterm, mj_value 
      Real(8) :: C
      Real(8), external :: DETC_jq, Clebsh2

! ... shell values:
!     md(i)  - the max.number of det.'s for the i-th subshell
!     nd(i)  - determinant under consideration
!     ipn(i) - pointers on the orbitals of given shell
!     MJs    - shells MJ
!     MJi    - intermediate values MJ

      Integer :: md(no),nd(no),ipn(no), MJs(no),MJi(no), ips(no)
      Integer :: Idet(ne)

! ... prepare working arrays:

      k = 1; mkdt = 1
      Do i=1,no
       ipn(i)=k; k=k+iq(i); md(i)=Ndets_jq(jn(i),iq(i))
       mkdt = mkdt*md(i)
       ips(i) = Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ)
      End do

      if(allocated(C_det)) Deallocate(C_det); Allocate(C_det(mkdt))
      if(allocated(IP_det)) Deallocate(IP_det); Allocate(IP_det(ne,mkdt))

!--------------------------------------------------------------------
! ... exausting all possible determinants:

      kdt=0; i=1; nd(i)=1              
    1 Call DET_sh_jq(jn(i),iq(i),nd(i),MJs(i),Idet(ipn(i)))

      if(i.eq.1) then
       MJi(1) = MJs(1)
      else
       MJi(i) = MJi(i-1)+MJs(i)
      end if
      if(i.lt.no) then;  i=i+1;  nd(i)=1;  go to 1; end if

       ! ... coefficient calculation:

       if(MJi(no).ne.Jtotal) go to 2
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

       if(C.ne.0.d0) then
        kdt=kdt+1; IP_det(1:ne,kdt)=Idet(1:ne); C_det(kdt)=C
       end if

    2 nd(i)=nd(i)+1                ! selecting the next case

      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          ! to end
       i=i-1; go to 2
      end if
      go to 1

    3 Continue

      Do k=1,kdt; Do i=1,ne 
       j=IP_det(i,k); IP_det(i,k)=mj_value(j)
      End do; End do 

      End Subroutine Det_expn_jj


!======================================================================
      Subroutine me_det
!======================================================================
!     find the possible interaction orbitals in two determinants and
!     call the subroutines to calculate m.e. between possible 
!     combinations of nj-orbitals (me_ee)
!----------------------------------------------------------------------
      Implicit none

      Integer :: i,i1,i2, j,j1,j2, idif, jdif
      Integer :: ii(ne),jj(ne)

!----------------------------------------------------------------------
! ... check total orbitals difference:

       ii = 0; jj=0
       Do i1=1,ne; j1=ip1(i1)
       Do i2=1,ne; j2=ip2(i2)
        if(nn1(j1).ne.nn2(j2)) Cycle
        if(ln1(j1).ne.ln2(j2)) Cycle
        if(jn1(j1).ne.jn2(j2)) Cycle
        if(IP_det1(i1,kd1).ne.IP_det2(i2,kd2)) Cycle 
        ii(i1)=i2; jj(i2)=i1
       End do; End do

       idif = 0;  Do i=1,ne; if(ii(i).eq.0) idif=idif+1; End do 
       jdif = 0;  Do i=1,ne; if(jj(i).eq.0) jdif=jdif+1; End do

       if(idif.ne.jdif) Stop 'me_det: idif <> jdif'
       if(idif.gt.2) Return

!----------------------------------------------------------------------
       Select case(idif)

       Case(0)

        Do i=1,ne; i1=ip1(i); i2=ip2(i)
         coef(i1,i1,i2,i2,-1) = coef(i1,i1,i2,i2,-1) + CC_det
        End do

        Do i1=1,ne-1;  j1=ii(i1)
        Do i2=i1+1,ne; j2=ii(i2)
         Call me_ee(i1,i2,min(j1,j2),max(j1,j2))
        End do; End do      

       Case(1)

        Do i = 1,ne; if(ii(i).ne.0) Cycle; i1=i; Exit; End do
        Do j = 1,ne; if(jj(j).ne.0) Cycle; j1=j; Exit; End do

        i=ip1(i1); j=ip2(j1) 
        if(IP_det1(i1,kd1).eq.IP_det2(j1,kd2).and. &
           jn1(i).eq.jn2(j).and.ln1(i).eq.ln2(j)) then
         coef(i,i,j,j,-1) = coef(i,i,j,j,-1) + CC_det * (-1)**(i1+j1)
        end if

        Do i2=1,ne; if(i2.eq.i1) Cycle; j2=ii(i2)
         Call me_ee(min(i1,i2),max(i1,i2),min(j1,j2),max(j1,j2))
        End do

       Case(2)

        Do i = 1,ne; if(ii(i).ne.0) Cycle; i1=i; Exit; End do
        Do i = ne,1,-1; if(ii(i).ne.0) Cycle; i2=i; Exit; End do
        Do j = 1,ne; if(jj(j).ne.0) Cycle; j1=j; Exit; End do
        Do j = ne,1,-1; if(jj(j).ne.0) Cycle; j2=j; Exit; End do

        Call me_ee(i1,i2,j1,j2)

       End Select

       End Subroutine me_det


!======================================================================
      SUBROUTINE me_ee (i1,j1,i2,j2)
!======================================================================
!     angular part of matrix elements in nljm-representation
!     for two-electron operator
!     Calls: Check_boef
!----------------------------------------------------------------------
      Use boef_list

      Implicit none
      Integer, intent(in) :: i1,i2,j1,j2
      Integer :: k,kz,ib, n1,n2,n3,n4

      if(IP_det1(i1,kd1)+IP_det1(j1,kd1).ne. &
         IP_det2(i2,kd2)+IP_det2(j2,kd2)) Return

      n1=ip1(i1); n2=ip1(j1); n3=ip2(i2); n4=ip2(j2)

      Call Check_boef(ln1(n1),jn1(n1),IP_det1(i1,kd1), & 
                      ln1(n2),jn1(n2),IP_det1(j1,kd1), &
                      ln2(n3),jn2(n3),IP_det2(i2,kd2), &
                      ln2(n4),jn2(n4),IP_det2(j2,kd2))

      kz = (-1)**(i1+i2+j1+j2)

      Do ib = ncblk(kblk-1)+1,ncblk(kblk)
       k = IB_int(ib)
       if(k.gt.0) then
        k = k - 1;    if(k.gt.kmax) Cycle
        coef(n1,n2,n3,n4,k) = coef(n1,n2,n3,n4,k) + Boef(ib)*CC_det*kz
       else
        k = -k - 1;  if(k.gt.kmax) Cycle
        coef(n1,n2,n4,n3,k) = coef(n1,n2,n4,n3,k) + Boef(ib)*CC_det*kz
       end if        
      End do

      End Subroutine me_ee

      End Subroutine coef_2conf_jj



