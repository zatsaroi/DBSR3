!======================================================================
      Subroutine Recup_jj_shells
!======================================================================
!     recoupling equivalent jj shells:
!
!     shift:     <[(j1,j2)J',j3]J || [j1,(j2,j3)J'']J>  
!     jump:      <[(j1,j2)J',j3]J || [j1,(j3,j2)J'']J>
!
!     j1,j2,j3 =>  jot1,jot2,jot3
!     J'       =>  jotp
!     J''      =>  jot
!     J        =>  JT
!
!----------------------------------------------------------------------

      Use jj2ls; Use conf_jj, only: msh

      Implicit none

      Integer :: jot_min(msh), jot_max(msh)
      Real(8) :: D(msh)
      Integer :: i,i1,i2
      Real(8), external :: recup_shift, recup_jump 

      Do i = 1,n
       jot_min(i) = iabs(jot2(i)-jot3(i))
       jot_max(i) = iabs(jot2(i)+jot3(i))
       jot(i) = jot_min(i)
      End do

      jot(1) = JT(1)
      
      i1=2                     ! i1 - low  limit of shells
      i2=n                     ! i2 - high limit of shells
      D(1) = one

      i=i1
    1 Continue

      if(ncase(i).eq.-1.or.ncase(i).eq.1) then
       D(i)=D(i-1)
      elseif(ncase(i).eq.-2) then
       D(i)=D(i-1)*recup_shift(jot1(i),jot2(i),jot3(i),jotp(i),jot(i),JT(i))
      else
       D(i)=D(i-1)*recup_jump (jot1(i),jot3(i),jot2(i),jotp(i),jot(i),JT(i))
      end if

      if(i.lt.i2) then
       i=i+1
       jot(i)=jot_min(i)
       go to 1
      end if
      
      if(D(i2).ne.zero)  Call sub_nlLS(D(i2))

    2 jot(i)=jot(i)+2
      if(jot(i).le.jot_max(i)) go to 1
      i=i-1
      if(i.ge.i1) go to 2
      
      End  Subroutine Recup_jj_shells


