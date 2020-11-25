!==================================================================
      Subroutine Solve_eiv(i,hfm,v,rhs)
!==================================================================
!     Find the eigenvector of hfm for the m'th eigenvalue after 
!     orthogonality has been applied
!------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Integer, intent(in)    :: i
      Real(8), intent(inout) :: hfm(ms,ms), rhs(ms)
      Real(8), intent(out)   :: v(ms)

      Real(8) :: aa(ms,ms), ss(ms,ms), w(3*ms)
      Real(8) :: eval(ms), a(ms), s(ms), rh(ms), zz,eii,C
      Integer :: j, jp, info, k,ii,m,mm, ipos(1)
    
      write(log,'(a,1Pe12.2)') 'srhs = ', srhs

      if(srhs.gt.0.d0.and.method.eq.1) then
       Call solve_direct (i,hfm,v,rhs); Return
      end if

      write(log,'(a)') 'method - solve_eiv'

      eii = e(i)

! ... apply orthogonality conditions for orbitals through projection: 

      m = nbs(i)-lbs(i)
      Do j = 1, nbf; if(i.eq.j) Cycle  ! i-1  ???               
       if(i.le.ncore.and.j.ge.i) Cycle
       if(kbs(i).ne.kbs(j)) Cycle
       Call apply_orthogonality(hfm,p(1,1,j))
       if(m.gt.1) m = m - 1   
      End do

      if(iphys(i).eq.0) m=1

! ... apply boundary conditions (delete extra B-splines):

      ii=0
      Do j=1,ms
       if(iprm(j,i).eq.0) Cycle; ii=ii+1
       k=0
       Do jp=1,ms
        if(iprm(jp,i).eq.0) Cycle
        k=k+1; a(k)=hfm(jp,j); s(k)=fppqq(jp,j)  
       End do
       aa(1:k,ii)=a(1:k); ss(1:k,ii)=s(1:k); rh(ii) = rhs(j)
      End do 

! ... evaluates the eigenvalues and eigenvectors (LAPACK routine):

      Call dsygv(1,'V','L',ii,aa,ms,ss,ms,eval,w,3*ms,INFO)
      if (info /= 0) then
       WRITE(scr,'(a,i6)') 'scf_eiv: error in eigenvalue routine, dsygv', info
       Stop 
      end if

! ... choose the solution if no rhs:

      mm=0; zz = -1.5*c_au**2
      Do j=1,ii; mm = j
       if(eval(j).gt.zz) Exit
      End do
      k=mm; mm = mm + m - 1

      if(debug.gt.0) write(log,'(a,5E15.5)') 'eval =',eval(k:k+4)

!----------------------------------------------------------------------
      if (srhs.eq.0.d0) then    ! without rhs 

      ! restore the solutions in original B-spline net:

      a(1:ms) = aa(1:ms,mm);  v=0.d0; k=0
      Do j=1,ms
       if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
      End do 

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

      e(i) = eval(mm)

!-----------------------------------------------------------------------
! ... Use inverse iteration in case of rhs: 

      else   !  with rhs

      write(log,'(a)') 'method - inverse iterations'

       a = 0.d0
       Do j = 1,ii
        if(eval(j).lt.-2*z*z) Cycle
        if(abs(eval(j)).lt.0.00001) Cycle
        C = Dot_Product(rh(1:ii),aa(1:ii,j))/(eval(j) - eii)
        a(1:ii) = a(1:ii) + aa(1:ii,j)*C
       End do

! ... restore the solutions in original B-spline net and normalize:

       v = 0.d0; k=0
       Do j=1,ms
        if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
       End do 

       a = Matmul(fppqq,v)
       v = v/sqrt(Dot_Product(v,a))

! ... if (v(ks) < 0.d0) v = -v

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

      end if
      
      End Subroutine  Solve_eiv
