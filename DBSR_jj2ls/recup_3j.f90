!======================================================================
      Real(8) Function recup_shift(j1,j2,j3,jp,jpp,j)
!======================================================================
!     <[(j1,j2)J',j3]J || [j1,(j2,j3)J'']J> 
!----------------------------------------------------------------------
     
      Implicit none

      Integer, intent(in) :: j1,j2,j3,jp,jpp,j
      Integer :: k
      Real(8) :: R
      Real(8), external :: Z_6j2
        
      R = dfloat((jp+1)*(jpp+1))
      R = dsqrt(R)
      R = R * Z_6j2(j1,j2,jp,j3,j,jpp)
      k = (j1+j2+j3+j)/2
      recup_shift = R * (-1)**k

      End Function recup_shift


!======================================================================
      Real(8) Function recup_jump(j1,j2,j3,jp,jpp,j)
!======================================================================
!     <[(j1,j2)J',j3]J || [j1,(j3,j2)J'']J>
!     <[(j1,j3)J',j2]J || [j1,(j2,j3)J'']J>  we need   
!----------------------------------------------------------------------

      Implicit none

      Integer, intent(in) :: j1,j2,j3,jp,jpp,j
      Integer :: k
      Real(8) :: R
      Real(8), external :: Z_6j2
        
      R = dfloat((jp+1)*(jpp+1))
      R = dsqrt(R)
      R = R * Z_6j2(j1,j2,jp,j3,j,jpp)
      k = (j1+j+jpp)/2
      recup_jump = R * (-1)**k

      End Function recup_jump


!======================================================================
      Real(8) Function recup_exchange(j1,j2,j3,jp,jpp,j)
!======================================================================
!     <[(j1,j2)J',j3]J || [(j1,j3)J'',j2]J>   
!----------------------------------------------------------------------

      Implicit none

      Integer, intent(in) :: j1,j2,j3,jp,jpp,j
      Integer :: k
      Real(8) :: R
      Real(8), external :: Z_6j2
        
      R = dfloat((jp+1)*(jpp+1))
      R = dsqrt(R)
      R = R * Z_6j2(j2,j1,jp,j3,j,jpp)
      k = (j2+j3+jp+jpp)/2
      recup_exchange = R * (-1)**k

      End Function recup_exchange

