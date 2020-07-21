!======================================================================
      Subroutine DIAG(ishift) 
!======================================================================
!     this subroutine calls the LAPACK (htmm://www.netlib.org) routines 
!     (DSYGV or DSYGVX) to solve the generalized eigenvalue problem
!
!                Hx = E Sx,    NORT > -1,
!
!     or (DSYEV or DSYEVX) to solve standard eigenvalue problems 
!
!                Hx = E x,     NORT = -1,
!
!     and records the resulting expansions to output file (unit 'nuj')
!----------------------------------------------------------------------
      Use dbsr_ci

      Implicit real(8) (A-H,O-Z)
      Integer, parameter :: mlab = 64
      Character(mlab) :: Label

! ... include the core energy:

      HM = HM + Ecore*SM

! ... record interaction matrix for further use if needed: 

      AF_mat = trim(name)//'.mat' 
      open(num,file=AF_mat,form='UNFORMATTED')
      Do i=1,ncj;  write(num) HM(i,1:i); End do

      if(debug.gt.0) then                     ! print interaction matrix
       write(pri,'(/a/)') ncj, '   Overlap matrix' 
       Do i=1,ncj
        write(pri,'(i5)') i
        write(pri,'(5f16.8)') SM(i,1:i)
       End do
       write(pri,'(/a/)') ncj, '   Interaction matrix' 
       Do i=1,ncj
        write(pri,'(i5)') i
        write(pri,'(5f16.8)') HM(i,1:i)
       End do
      end if

!----------------------------------------------------------------------

      if(nshift.gt.0) then
       Do i=1,ncj
        HM(i,i) = HM(i,i) + shift(i+ishift)
       End do
      end if

!----------------------------------------------------------------------
!     solve the  HX = E SX  problem with LAPACK routines:

      if(NORT.gt.-1) then   ! generalized eigenvalue problem

!       if(meiv.lt.nzero) then
!        Call LAP_DSYGVX('V','L',nzero,ncj,HM,SM,eval,meiv,info)
!       else
        Call LAP_DSYGV ('V','L',nzero,ncj,HM,SM,eval,info)
!       end if
 
      else                  ! simple eigenvalue problem

       if(meiv.lt.nzero) then
        Call LAP_DSYEVX('V','L',nzero,ncj,HM,eval,meiv,info)
       else
        Call LAP_DSYEV ('V','L',nzero,ncj,HM,eval,info)
       end if

      end if

      SM(1:nzero,1:meiv) = HM(1:nzero,1:meiv)
      rewind(num)
      Do i=1,ncj;  read(num) HM(i,1:i); End do

! ... FIRST-ORDER CORRECTIONS:

      if(nzero.lt.ncj)  CALL FIRST

! ... PRINT OUT THE EIGENVALUES AND EIGENVECTORS:

      m=0; Do i=1,meiv; m=m+1; if(EVAL(i).gt.Emax) Exit; End do

      DO k = 1,m
       CM=0.d0; ii=0                 ! find the leading component
       DO i=1,ncj
        if(abs(SM(i,k)).le.CM) Cycle; CM=abs(SM(i,k)); ii=i
       End do
       ii = ii + ishift
       Call Get_cfg_jj(ii); Call Label_jj(mlab,Label,0)
       nsol=nsol+1
       write(nuj,'(i8,2x,a)')  nsol, TRIM(Label)
       write(nuj,'(f16.8,3i8,f16.4)')  EVAL(k), jot, ishift+1,ishift+ncj,EVAL(k)*au_eV
       write(nuj,'(6f20.15)') (SM(j,k),j=1,ncj)

!       EE = 0.d0
!       Do i=1,ncj; Do j=1,i
!        S = SM(i,k)*SM(j,k)*HM(i,j)
!        EE = EE + S; if (i.ne.j) EE = EE + S
!       End do; End do
!       write(pri,'(a,3f16.8)') 'Energy  =',EE,EVAL(k),EE-EVAL(k)

      End do

      Close(num,status='DELETE')

      End Subroutine DIAG


!=======================================================================
      Subroutine FIRST 
!=======================================================================
!     Compute the first-order corrections to the eigenvalue and
!     eigenvectors.
!-----------------------------------------------------------------------
      Use dbsr_ci

      Implicit none
      Integer :: i,j
      Real(8) :: E0,E1,x,x1,x2,y,y1,y2

      if(NORT.gt.-1) then   ! generalized eigenvalue problem

       Do i = 1,MEIV
        E0 = EVAL(i)
        E1 = 0.d0
        y1 = 1.d0
        y2 = 1.d0
        DO j = NZERO+1,ncj
         x1 = SUM(HM(j,1:nzero)*SM(1:nzero,i)) 
         x2 = SUM(SM(j,1:nzero)*SM(1:nzero,i)) 
         x  = x1 - E0*x2
         SM(j,i) = x/(E0-DM(j))
         E1 = E1 + x*SM(j,i)
         y1 = y1 + SM(j,i)*x2
         y2 = y2 + SM(j,i)*SM(j,i)
        END DO
        EVAL(i) = EVAL(i) + E1/(y1*y2)
        y1 = 1.D0/SQRT(y2)
        SM(1:ncj,i) = y1*SM(1:ncj,i)
       End do

      else                   !  simple eigenvalue problem

       Do i = 1,MEIV
        E0 = EVAL(i)
        E1 = 0.d0
        y2 = 1.d0
        Do j = NZERO+1,ncj
         y = SUM(HM(j,1:nzero)*SM(1:nzero,i))
         SM(j,i) = y/(E0-DM(j))
         E1 = E1 + SM(j,i)*y
         y2 = y2 + SM(j,i)*SM(j,i)
        End do
        EVAL(i) = EVAL(i) + E1/y2
        y1 = 1.D0/SQRT(y2)
        SM(1:ncj,i) = y1*SM(1:ncj,i)
       End do

      end if

      End Subroutine FIRST


!======================================================================
      Subroutine Check_cfile 
!======================================================================
!     check the expansion in c-file by comparing 
!     the calculated energy and the energy recorded in c-file
!----------------------------------------------------------------------
      Use dbsr_ci
      Use conf_jj

      Implicit real(8) (A-H,O-Z)

! ... include core energy to interaction matrix:

      HM = HM + Ecore*SM

! ... read expansion:

      Call read_expn_jj(nuc)

! ... calculate energy: 

      EE = 0.d0
      Do i=1,ncj; Do j=1,i
       S = WC(i)*WC(j)*HM(i,j)
       EE = EE + S; if (i.ne.j) EE = EE + S
      End do; End do

! ... calculate total normalization: 

      ss = 0.d0
      Do i=1,ncj; Do j=1,i
       S = WC(i)*WC(j)*SM(i,j)
       ss = ss + S; if (i.ne.j) ss = ss + S
      End do; End do
      
      write(pri,'(a,f16.8)') 'Overlap =',SS
      
      rewind(nuc); read(nuc,'(15x,f16.8)') E0

      write(pri,'(a,3f16.8)') &
        'Energy (calulated, from c-fils, difference:', EE,E0,EE-E0

      End Subroutine Check_cfile 



