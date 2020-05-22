!======================================================================
      Subroutine make_coupling
!======================================================================
!     prepare coupling scheme, J1_coupling, for state in module conf_jj
!----------------------------------------------------------------------
      Use dbsr_conf,  JJ_coupling => J1_coupling, n=>no

      if(n.gt.mshells) Stop 'make_coupling: no > mshells'

! ... JJ-coupling:

      JJ_coupling(1,1) = 1
      JJ_coupling(2,1) = 2
      JJ_coupling(3,1) = n+2
	  
      Do i=2,n-1
       JJ_coupling(1,i) = n+i
       JJ_coupling(2,i) = i+1
       JJ_coupling(3,i) = n+i+1
      End do

! ... moments (in 2J+1 representation)

      Do i=1,n
       moments(i)  =  Jshell(i) + 1        
       moments(i+n) = Jintra(i) + 1
      End do

      ncup = n-1
      nmom = 3*n

      End Subroutine make_coupling


!======================================================================
      Subroutine make_coupling_insert(ii)
!======================================================================
!     prepare coupling scheme, J2_coupling, for state in module conf_jj
!     in case when outer orbital is going to new shell 'ii'
!----------------------------------------------------------------------
      Use dbsr_conf,  JJ_coupling => J2_coupling, n=>no   

      Implicit none
      Integer :: i,ii,m,ipos(mshells)
 
      if(n.gt.mshells) Stop 'make_coupling: no > mshells'

      Do i=1,n
       if(i.lt.ii) then;       ipos(i)=i
       elseif(i.eq.ii) then;   ipos(i)=n
       else;                   ipos(i)=i-1
       end if
      End do

! ... JJ-coupling:

      m = n + n
      JJ_coupling(1,1) = ipos(1)
      JJ_coupling(2,1) = ipos(2)
      JJ_coupling(3,1) = m+2
	  
      Do i=2,n-1
       JJ_coupling(1,i) = i+m
       JJ_coupling(2,i) = ipos(i+1)
       JJ_coupling(3,i) = i+m+1
      End do

! ... moments

      Do i=1,n;  moments(i+m) = Jintra(i) + 1;  End do

      JJ_coupling(3,ncup)=J1_coupling(3,ncup)
 
      End Subroutine make_coupling_insert


!======================================================================
      Subroutine make_coupling_trap
!======================================================================
!     prepare coupling scheme, J2_coupling, for state in module conf_jj
!     in case when outer orbital is traped to shell 'ii'
!----------------------------------------------------------------------
      Use dbsr_conf,  JJ_coupling => J2_coupling   

      Implicit none
      Integer :: i,ii,n,m
 
      n = ncup + 1
      m = n + n

! ... JJ-coupling:

      ii = iabs(insert)

      JJ_coupling(1,1) = ii
      JJ_coupling(2,1) = n
      JJ_coupling(3,1) = m+1

      JJ_coupling(1,2) = 1
      JJ_coupling(2,2) = 2
      JJ_coupling(3,2) = m+2
	  
      if(ii.eq.1) then
       JJ_coupling(1,2) = m+1
       JJ_coupling(2,2) = 2
       JJ_coupling(3,2) = m+2
      end if

      Do i=2,no-1
       JJ_coupling(1,i+1) = m+i
       JJ_coupling(2,i+1) = i+1
       JJ_coupling(3,i+1) = m+i+1
      End do

      if(ii.gt.1) JJ_coupling(2,ii) = m+1

      JJ_coupling(3,ncup)=J1_coupling(3,ncup)

! ... moments

      moments(m+1) = Jshell(ii)+1
      Do i=2,no;  moments(m+i) = Jintra(i)+1;  End do

      End Subroutine make_coupling_trap


