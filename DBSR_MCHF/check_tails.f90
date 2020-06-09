!======================================================================
      Subroutine Check_tails(jo)
!======================================================================
!     nulify too small B-spline coefficients in the end of expansion; 
!     if jo > 0, action is applied only to orbital "jo"
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals
      
      Implicit none
      Integer :: io,jo,i
      Real(8) :: cm

      Do io=1,nbf; if(jo.ne.0.and.io.ne.jo) Cycle
       cm = maxval( abs( p(:,1,io) ) )
       Do i=nsp,1,-1
        mbs(io)=i     
        if(abs(p(i,1,io))/cm.lt.end_tol) Cycle
        Exit
       End do
       p(i+1:nsp,1,io)=0.d0
       p(i+1:nsq,2,io)=0.d0
      End do

      End Subroutine Check_tails
