!=====================================================================
      Subroutine Solve_direct (i,hfm,v,rhs)
!=====================================================================
!     Direct solution of the one-electron equation:  hfm v - e v = rhs 
!---------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals
 
      Implicit none
      Integer :: i
      Real(8), intent(in)  :: hfm(ms,ms), rhs(ms)
      Real(8), intent(out) :: v(ms)

      Real(8), dimension(ms+nbf,ms+nbf) :: aa, ax
      Real(8), dimension(ms+nbf) :: res, xx, yy
      Integer, dimension(ms+nbf) :: ipiv

      Real(8) :: y(ms), w(ms)
      Integer :: mm, md, k, info, j, jp, ii, ipos(1)

      write(log,'(a)') 'method - direct'

      md = ms+nbf
      mm = ms

      aa = 0.d0;  aa(1:ms,1:ms) = hfm - e(i)*fppqq
      res = 0.d0; res(1:ms) = rhs
 
! ... add orthogonality with all orbitals:

      Do j = 1,i-1                             !  ???
       if(kbs(i).ne.kbs(j)) Cycle
       mm = mm + 1
       Call Get_pv_df(j,w)
       y = matmul(fppqq,w)
       aa(1:ms,mm) = -y
       aa(mm,1:ms) = -y
      End do

! ... apply boundary conditions (delete the extra B-splines):

      ii=0
      Do j=1,mm
       if(j.le.ms) then; if(iprm(j,i).eq.0) Cycle; end if; ii=ii+1
       k=0
       Do jp=1,mm
        if(jp.le.ms) then; if(iprm(jp,i).eq.0) Cycle; end if
        k=k+1; xx(k)=aa(jp,j)
       End do
       ax(1:k,ii)=xx(1:k); yy(ii)=res(j)
      End do 

! ... solve the matrx equation with LAPACK routines:

      CALL DGETRF(ii, ii, ax, md, ipiv, info)
      if (info /= 0) Stop 'solve_direct: error in factorization routine DSYTRF'
      CALL DGETRS('N', ii, 1, ax, md, ipiv, yy, md, info)
      if (info /= 0) Stop 'solve_direct: error in solve routine DSYTRS'

! ... restore the solution in orginal basis:

      v=0.d0; k=0
      Do j=1,ms
       if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=yy(k)
      End do 

! ... normalize the solution:

      y = matmul(fppqq,v)
      v = v/sqrt(Dot_Product(v,y))

! ... solution sign:

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

      End Subroutine  Solve_direct
