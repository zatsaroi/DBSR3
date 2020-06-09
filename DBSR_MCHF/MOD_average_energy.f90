!=======================================================================
      Module av_energy_expression
!=======================================================================
!     This module defines the energy expression for closed core shells
!-----------------------------------------------------------------------
      Implicit none

      Integer :: kclosed, kwf, kmax 
      Real(8), allocatable :: coef(:,:,:)

      End module av_energy_expression


!=======================================================================
      Real(8) Function a(i,j,k)
!=======================================================================
!     coefficient for the direct interaction i - j:  Rk(i,j;i,j)
!-----------------------------------------------------------------------
      Use av_energy_expression
    
      Implicit none
      Integer, intent(in) :: i,j,k
    
      a = 0.d0
      if(i.lt.1.or.i.gt.kclosed) Return
      if(j.lt.1.or.j.gt.kwf) Return
      if(k.lt.0.or.k.gt.kmax) Return
      if(k.eq.0) then
       a = 1.d0
      elseif(i.eq.j) then
       a = coef(i,j,k)
      end if

      End Function a

    
!=======================================================================
      Real(8) Function b(i,j,k) 
!=======================================================================
!     coefficient for the exchange interaction i - j:  Rk(i,j;j,i)
!-----------------------------------------------------------------------
      Use av_energy_expression
  
      Implicit none
      Integer, intent(in) :: i,j,k

      b = 0.d0
      if(i.eq.j) Return 
      if(i.lt.1.or.i.gt.kclosed) Return
      if(j.lt.1.or.j.gt.kwf) Return
      if(k.lt.0.or.k.gt.kmax) Return
      b = coef(i,j,k)
     
      End Function b
  

!======================================================================
      Subroutine av_energy_coef(ncore,nwf,lbs,jbs)
!======================================================================
!     This routine gets HF coefficients for average enery approximation
!----------------------------------------------------------------------
      Use av_energy_expression

      Implicit none
      Integer, intent(in) :: ncore, nwf, lbs(nwf), jbs(nwf)
      Integer :: i,j,k
      Real(8), external :: Cjkj 
	  
      kclosed = ncore
      kwf = nwf
      if(ncore.le.0) Return
      kmax = 2*maxval(lbs)
      
      Allocate(coef(kclosed,kwf,0:kmax)); coef=0.d0

! ... define average energy arrays

      Do i = 1,kclosed; Do j = 1,kwf
       if(i.eq.j) then
        Do k=2,lbs(i)+lbs(i),2
         coef(i,j,k) = -Cjkj(jbs(i),k,jbs(i))**2 /(jbs(i)*(jbs(i)+1))
        End do
       else
        Do k=iabs(lbs(i)-lbs(j)),lbs(i)+lbs(j),2
         coef(i,j,k) = -Cjkj(jbs(i),k,jbs(j))**2 /((jbs(i)+1)*(jbs(j)+1))
        End do
       end if
      End do; End do
	  
      End Subroutine av_energy_coef 



