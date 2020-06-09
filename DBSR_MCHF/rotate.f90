!======================================================================
      Subroutine Rotate_ij(i,j)
!======================================================================
!     Determine the rotation parameter for a stationary solution for
!     the (i,j) orbitals, i <= j
!     Rules derived from transformation of densities. Let
!    
!     P*i = (Pi - e*Pj)/sqrt(1+e^2)
!     P*j = (Pj + e*Pi)/sqrt(1+e^2)
!    
!     We compute three values for sum of contributions to energy from 
!     terms that depend on i,j:  
!   
!     F(+1) : sum with e =  eps  (F+)
!     F(0)  :     with e =  0    (F0)
!     F(-1) : sum with e = -eps  (F-)
!   
!     where eps is a numerically selected value
!   
!     From F(e) = F(0)  + e F'(0) + e^2 F"(0)/2  we get
!   
!     F"(0) = [F+ - 2F0 - F-]/e^2  
!   
!     [F+ - F0]/(e) = F'(0) + (e/2) F"
!----------------------------------------------------------------------
      Use dbsr_mchf, pi_value => pi
      Use df_orbitals, only: nbf,kbs
  
      Implicit none
      Integer, Intent(in) :: i,j
      Integer :: ip,kpol,j1,j2,j3,j4,int,m
      Real(8) :: pi(ms), pj(ms),vi(ms), vj(ms), f(-1:1)
      Real(8) :: eps,g,dg, c,s
      Real(8), external :: rk_df, quadr, dhc_value,WW
  
      if(i.ge.j) Return 
      if(i.le.nbf-varied) Return
      if(j.le.nbf-varied) Return
      if(kbs(i).ne.kbs(j)) Return
      if(i.le.ncore.or.j.le.ncore) Return
  
      g = 0.d0; dg = 0.d0; f = 0.d0
      Call Get_pv_df(i,pi)
      Call Get_pv_df(j,pj)  
      
      Do ip = -1,1
       eps = ip*eps_rot
       vi = (pi -eps*pj)/sqrt(1.+eps*eps)
       vj = (pj +eps*pi)/sqrt(1.+eps*eps)
       Call Put_pv_df(i,vi)
       Call Put_pv_df(j,vj)  

       Do int = 1,Lint
        j1 = i1_Lint(int); j2 = i2_Lint(int)
        Do m = ip_Lint(int-1)+1,ip_Lint(int) 
         C = WW(ic_Lcoef(m),jc_Lcoef(m))*L_coef(m)   
         if((j1-i)*(j2-i)*(j1-j)*(j2-j).ne.0) Cycle
         S = dhc_value(j1,j2)  
         f(ip) = f(ip) + C*S
        End do
       End do
       
       Do kpol=kmin,kmax
       Do int = nk_int(kpol-1)+1,nk_int(kpol)
        j1 = i1_int(int); j2 = i2_int(int); j3 = i3_int(int); j4 = i4_int(int)
        Do m = ip_int(int-1)+1,ip_int(int) 
         C = WW(ic_coef(m),jc_coef(m))*Rk_coef(m)   
         if((j1-i)*(j2-i)*(j3-i)*(j4-i)*(j1-j)*(j2-j)*(j3-j)*(j4-j).ne.0) Cycle
         S = rk_df(j1,j2,j3,j4,kpol) 
         f(ip) = f(ip) + C*S
        End do
       End do; End do
      End do  

      if(debug.gt.0) write(log,'(a,3f16.8)') 'F-, F0, F+', f
      eps = eps_rot
      g = (F(+1)-F(-1))/(2*eps)
      if(abs(g).lt.1.d-5) g=0.d0
      dg = (F(-1)-2*F(0)+F(+1))/(eps*eps)
      eps = 0.d0 
      if(abs(g).gt.1.d-5) eps = -g/dg
      if(debug.gt.0) write(log,'(a,3F16.8)') 'g, dg, eps', g, dg, eps
      vi = (pi - eps*pj)/sqrt(1.d0+eps*eps)
      vj = (pj + eps*pi)/sqrt(1.d0+eps*eps)
      Call Put_pv_df(i,vi)
      Call Put_pv_df(j,vj)  
      Call Check_tails(i)
      Call Check_tails(j)
  
      End Subroutine Rotate_ij
        
            
    
  
