!======================================================================
     Integer Function Idef_dtype(i1,i2,ich,jch,io,jo,ic,jc,k1,k2,k3)
!======================================================================
!    we have following 10 different structures for radial integrals:
!
! 1  d( . .)  ic, jc               -  bound-bound  
!
! 2  d( i .)  jc                   -  bound-channel
! 3  d( . j)  ic                    
! 4  d( . .) < i | . > jc           
! 5  d( . .) < . | j > ic           
!
! 6  d( i j)                       -  channel-channel
! 7  d( i .) < . | j >              
! 8  d( . j) < i | . >               
! 9  d( . .) < i | . > < . | j >    
!10  d( . .) < i | j >              
!
! where . denotes bound orbital (i1,i2), i,j - channels (ich,jch),
! ic,jc - configurations, <|> - overlaps (io).
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i1,i2,ich,jch,io,jo,ic,jc
      Integer, intent(out) :: k1,k2,k3

      if(io.eq.0) then
       if(ic.ne.0.and.jc.ne.0) then; Idef_dtype=1;  k1=ic; k2=jc; k3=0
       elseif(jc.ne.0) then;         Idef_dtype=2;  k1=i1; k2=i2; k3=jc
       elseif(ic.ne.0) then;         Idef_dtype=3;  k1=i1; k2=i2; k3=ic
       else;                         Idef_dtype=6;  k1=i1; k2=i2; k3=0
       end if
      elseif(jo.gt.0) then;          Idef_dtype=9;  k1=io; k2=jo; k3=0
      elseif(jc.gt.0) then;          Idef_dtype=4;  k1=io; k2=jc; k3=0
      elseif(ic.gt.0) then;          Idef_dtype=5;  k1=io; k2=ic; k3=0
      elseif(ich.gt.0) then;         Idef_dtype=7;  k1=i1; k2=i2; k3=io
      elseif(jch.gt.0) then;         Idef_dtype=8;  k1=i1; k2=i2; k3=io
      else;                          Idef_dtype=10; k1=io; k2=0;  k3=0
      end if                    

      End Function Idef_dtype
