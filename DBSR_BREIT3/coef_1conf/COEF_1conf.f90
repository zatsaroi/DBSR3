!======================================================================
      Subroutine coef_1conf(no,ln,jn,iq,Jshell,Vshell,Jintra,kmax,coefs)
!======================================================================
!     compute the angular coefficients for 1 atomic states
!     input atomic state and output coefficients are in modul coef_jj   
!----------------------------------------------------------------------

      Implicit none 

! ... input-output:

      Integer, intent(in) :: no,ln(no),jn(no),iq(no), kmax, &
                             Jshell(no),Vshell(no),Jintra(no)
      Real(8), intent(out):: coefs(no,no,0:kmax)

! ... determinant expansion:

      Integer :: kdt 
      Integer, Allocatable :: IP_det(:,:)
      Real(8), Allocatable :: C_det(:)
      Real(8) :: CC_det

! ... local variables:

      Integer             :: ne,Jtotal,i,j,k, jmax,kd1,kd2
      Integer             :: ip1(no),ip2(no)
      Integer, External   :: mj_value

! ... initialize arrays:

      coefs = 0.d0
      ne = SUM(iq(1:no))
      Jtotal = Jintra(no)
      ip1(1)=1; ip2(1)=iq(1)
      Do i=2,no
       ip1(i)=ip2(i-1)+1
       ip2(i)=ip2(i-1)+iq(i)
      End do 

! ... determinant expansion:

      Call Det_expn_jj  

! ... calculations:       total M - ???  Clebsh - ???

      Do kd1 = 1,kdt
      Do kd2 = kd1,kdt
        CC_det = C_det(kd1) * C_det(kd2)
	if(kd1.ne.kd2) CC_det = CC_det + CC_det 
        Call me_det
      End do;  End do


CONTAINS

!======================================================================
      Subroutine Det_expn_jj  
!======================================================================
!     procedure of exaustion of all possible determinants for given
!     configurations. The determinants and their coefficients
!     are recoded on unit 'nua'
!
!     Calls: Det_sh_jq, DETC_jq, Clebsh2, Ndets_jq, Jterm, mj_value 
!----------------------------------------------------------------------

      Implicit none 

      Integer :: i,j,k, mkdt, JW,JQ

      Real(8) :: C
      Real(8), External :: DETC_jq, Clebsh2
      Integer, External :: Ndets_jq, Jterm, mj_value 

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
      Integer :: i,i1,i2, j,j1,j2, k,kk,k1,k2, idif, jdif, io,jo
      Integer :: is,js, isym1,isym2, jsym1,jsym2
      Integer :: ii(ne,ne)

!----------------------------------------------------------------------
! ... the same determinants:

      if(kd1.eq.kd2) then
       Do i=1,no;  Do is = ip1(i),ip2(i)
       Do j=i,no;  Do js = ip1(j),ip2(j)
        if(js.gt.is) Call me_ee(i,j,is,js,is,js)
       End do; End do
       End do; End do
       Return
      end if

!----------------------------------------------------------------------
! ... check total orbitals differenc:

       ii = 0
       Do i=1,no; k=ip1(i); kk=ip2(i)
        Do i1=k,kk; Do i2=k,kk
         if(IP_det(i1,kd1).eq.IP_det(i2,kd2)) ii(i1,i2)=1 
        End do; End do
       End do

       idif = ne - SUM(ii) 
       if(idif.gt.2) Return
       if(idif.ne.2) Stop 'Det_me: idif <> 2 ???'

! ... find interaction orbitals:
 
       k = 1;  jdif=0
       Do i=1,no; k=ip1(i); kk=ip2(i); k1=k; k2=k
        idif = iq(i) - SUM(ii(k:kk,k:kk))
        if(idif.eq.0) Cycle

! ... first orbital:

        if(jdif.eq.0) then
         io = i
         Do i1=k,kk
          if(SUM(II(i1,k:kk)).eq.1) Cycle; isym1=i1; Exit
	     End do 		 
	     Do i2=k,kk
	      if(SUM(II(k:kk,i2)).eq.1) Cycle; isym2=i2; Exit
	     End do
	     jdif = 1
	     idif = idif-1
         k1 = isym1+1
         k2 = isym2+1
        end if
		
        if(idif.eq.0) Cycle

! ... second orbital:

         jo = i
         Do i1=k1,kk
	      if(SUM(II(i1,k:kk)).eq.1) Cycle; jsym1=i1; Exit
	     End do 		 
	     Do i2=k2,kk
	      if(SUM(II(k:kk,i2)).eq.1) Cycle; jsym2=i2; Exit
	     End do
         Exit
        End do

        Call me_ee(io,jo,isym1,jsym1,isym2,jsym2)

       End Subroutine me_det


!======================================================================
      SUBROUTINE me_ee (i,j,i1,j1,i2,j2)
!======================================================================
!     angular part of matrix elements between two det.w.f.
!     for two-electron operator
!     Calls: Check_boef
!----------------------------------------------------------------------

      USE boef_list

      Implicit none
      Integer, intent(in) :: i,j,i1,i2,j1,j2
      Integer :: k,kz,ib

      if(IP_det(i1,kd1)+IP_det(j1,kd1).ne. &
         IP_det(i2,kd2)+IP_det(j2,kd2)) Return

      Call Check_boef(ln(i),jn(i),IP_det(i1,kd1), & 
                      ln(j),jn(j),IP_det(j1,kd1), &
                      ln(i),jn(i),IP_det(i2,kd2), &
                      ln(j),jn(j),IP_det(j2,kd2))

      kz = (-1)**(i1+i2+j1+j2)

      Do ib = ncblk(kblk-1)+1,ncblk(kblk)
       k = IB_int(ib)
       if(k.gt.0) then
         k = k - 1
         coefs(j,i,k) = coefs(j,i,k) + Boef(ib)*CC_det*kz
       elseif(i.eq.j) then
         k = -k - 1
         coefs(j,i,k) = coefs(j,i,k) + Boef(ib)*CC_det*kz
       else
         k = -k - 1
         coefs(i,j,k) = coefs(i,j,k) + Boef(ib)*CC_det*kz
       end if        
      End do

      End Subroutine me_ee

      End Subroutine coef_1conf



