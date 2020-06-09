!======================================================================
      Subroutine make_couplings
!======================================================================
!
!     Moments: Li, Si, Ji, LLi, SSi, JJi
!     
!     JJ - coupling
!
!     L1 S1 J1   J1  J2  JJ2         
!     L2 S2 J2   JJ2 J3  JJ3
!     ........   ..............
!     Ln Sn Jn   JJ(n-1) Jn JJn
!
!     LS - coupling
!
!     L1  L2 LL2       S1  S2  SS2       Ln Sn Jn  
!     LL2 L3 LL3       SS2 S3  SS3
!     ..............   ..............
!     LL(n-1) Ln LLn   SS(n-1) Sn SSn
!----------------------------------------------------------------------

      Use jj2ls

      if(n.gt.mshells) Stop 'make_coupling: n > mshells'

! ... JJ-coupling:

      Do i=1,n; j=i+n; k=j+n
       JJ_coupling(1,i) = i
       JJ_coupling(2,i) = j
       JJ_coupling(3,i) = k
      End do

      Do i=1,n-1; ii=i+n; jj=i+5*n; j=i+2*n+1;  
       JJ_coupling(1,ii) = jj
       JJ_coupling(2,ii) = j
       JJ_coupling(3,ii) = jj+1
      End do       
      JJ_coupling(1,n+1) = 2*n+1

! ... LS-coupling:

      Do i=1,n-1; ll=i+3*n
       LS_coupling(1,i) = ll
       LS_coupling(2,i) = i+1
       LS_coupling(3,i) = ll+1
      End do
      LS_coupling(1,1) = 1

      Do i=1,n-1; ii=i+n-1; ll=i+4*n
       LS_coupling(1,ii) = ll
       LS_coupling(2,ii) = i+n+1
       LS_coupling(3,ii) = ll+1
      End do
      LS_coupling(1,n) = n+1

      i = 2*n-1
      LS_coupling(1,i) = 4*n
      LS_coupling(2,i) = 5*n
      LS_coupling(3,i) = 6*n
       
      ncup = 2*n-1
      nmom = 6*n

 if(debug.gt.1) then
  write(pri,*) 'ncup',ncup
  write(pri,*) 'nmom',nmom
  Do i=1,ncup
   write(pri,'(3i4,6x,3i4)') JJ_coupling(:,i), LS_coupling(:,i)
  End do
 end if

      End Subroutine make_couplings

