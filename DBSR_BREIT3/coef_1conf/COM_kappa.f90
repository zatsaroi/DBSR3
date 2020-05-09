!=======================================================================
      Real(8) Function Ckap (kap1,k,kap2)
!=======================================================================
!     Computes the relativistic reduced matrix element                    
!                          k                                                        
!                (kap1 || C  || kap2)  
!
!     see Cjkj as origin
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: kap1,k,kap2
      Integer, External :: j_kappa
      Real(8), External :: Cjkj 
      Integer :: j1,j2

      j1 = j_kappa(kap1)
      j2 = j_kappa(kap1)
      Ckap = Cjkj(j1,k,j2)

      END Function Ckap 



!=======================================================================
      Integer Function kappa_lj(l,jj)
!=======================================================================
!     kappa = (l-j)*(2j+1);  jj = 2j
!-----------------------------------------------------------------------

      Integer :: l,jj

      kappa_lj = (2*l-jj)*(jj+1)/2

      End Function kappa_lj


!=======================================================================
      Integer Function l_kappa(kappa)
!=======================================================================

      Integer :: kappa

      if(kappa.eq.0) Stop 'l_kappa: kappa=0'
      if(kappa.gt.0) then
       l_kappa =  kappa 
      else
       l_kappa = -kappa-1
      end if 

      End Function l_kappa


!=======================================================================
      Integer Function j_kappa(kappa)
!=======================================================================
!     j_kappa = 2*j
!-----------------------------------------------------------------------

      Integer :: kappa

      if(kappa.eq.0) Stop 'j_kappa: kappa=0'
      if(kappa.gt.0) then
       j_kappa =  kappa+kappa-1
      else
       j_kappa = -kappa-kappa-1
      end if 

      End Function j_kappa


