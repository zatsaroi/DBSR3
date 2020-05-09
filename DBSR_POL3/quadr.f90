!======================================================================
      Real(8) Function quadr (fi,fj,m)
!======================================================================
!     Evaluates   <fi | r^m | fj>     with respect to r
!     where fi, fj  - two-component Dirac functions in B-spline basis
!----------------------------------------------------------------------
      Use DBS_grid    
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: m
      Real(8), intent(in) :: fi(ns,2),fj(ns,2)
      Integer :: i,j, iv,ith,jth
      Real(8) :: pi,pj,qi,qj

      quadr = 0.d0

      Do iv=1,nv   
       gx(:) = gr(iv,:)**m * grw(iv,:) 
       Do ith=1,ksp;  i=iv+ith-1; pi=fi(i,1)
        gw(:) = gx(:)*pbsp(iv,:,ith)
        Do jth=1,ksp;  j=iv+jth-1; pj=fj(j,1)
         quadr = quadr + SUM(pbsp(iv,:,jth)*gw(:))*pi*pj
        End do
       End do
      End do
      
      Do iv=1,nv   
       gx(:) = gr(iv,:)**m * grw(iv,:) 
       Do ith=1,ksq;  i=iv+ith-1; qi=fi(i,2)
        gw(:) = gx(:)*qbsp(iv,:,ith)
        Do jth=1,ksq;  j=iv+jth-1; qj=fj(j,2)
         quadr = quadr + SUM(qbsp(iv,:,jth)*gw(:))*qi*qj
        End do
       End do
      End do

      End Function quadr 



!======================================================================
      Real(8) Function quadrm (nv,ks,ks1,ks2,m,bsp1,bsp2,f1,f2)
!======================================================================
!     Evaluates   <P(Q)_i | r^m | P(Q)_j>   with respect to r
!----------------------------------------------------------------------
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: nv,ks,ks1,ks2,m
      Real(8), intent(in) :: bsp1(nv+1,ks,*),bsp2(nv+1,ks,*),f1(*),f2(*)
      Integer :: i,j, iv, ith,jth

      quadrm = 0.d0
      Do iv=1,nv   
       gx(:) = gr(iv,:)**m * grw(iv,:) 
       Do ith=1,ks1;  i=iv+ith-1
       gw(:) = gx(:)*bsp1(iv,:,ith)
        Do jth=1,ks2;  j=iv+jth-1
         quadrm = quadrm + SUM(bsp2(iv,:,jth)*gw(:))*f1(i)*f2(j)
        End do
       End do
      End do

      End Function quadrm 




