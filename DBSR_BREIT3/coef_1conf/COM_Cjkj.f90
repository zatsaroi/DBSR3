!=======================================================================
      Real(8) Function Cjkj (j1,k,j2)
!=======================================================================
!     Computes the relativistic reduced matrix element                    
!     (see Eq.15, Grant and Pyper, J.Phys.B9,761,1976)
!
!             k                                                        
!     (j1 || C  || j2)  = 
!
!                                       (j1   k  j2 )
!         (-1)^(j1+1/2) sqrt([j1][j2])      
!                                       (1/2  0 -1/2)                       
!
!     C(j,0,j) = sqrt([j])
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: j1,k,j2
      Real(8), External :: Z_3j2

      Cjkj = Z_3j2(j1,1,k+k,0,j2,-1) * &
            sqrt(Real((j1+1)*(j2+1))) * (-1)**((j1+1)/2)

      END Function Cjkj 


