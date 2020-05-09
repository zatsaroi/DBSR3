!======================================================================
      Integer(4) Function Jterm (j,q,k,JT,JV,JW,JQ)
!======================================================================
!
!     provides information about possible terms in the j^q subshell for
!     j <= 9/2
!
!     k > 0  --> JT,JV,JQ,JM = k-th term of j^q-subshell
!     k = 0  --> Jterm = position of the (JT,JV) term in subshell list
!     k < 0  --> Jterm = number of terms in j^q-subshell
!
!     j, JT, JQ  -->  in 2*J representation
!
!----------------------------------------------------------------------

      IMPLICIT NONE

      Integer, Intent(in) :: j,q,k

      Integer :: JT,JV,JQ,JW

      Integer :: qm,nterm,ip,i,ii,qq
      Integer :: term_list(3,65)

      Data term_list                                               &
       /1, 5, 2,  3, 3, 0,  3, 9, 0,                               & ! 5/2 ^ 3
        1, 7, 3,  3, 3, 1,  3, 5, 1,  3, 9, 1,  3,11, 1,  3,15, 1, & ! 7/2 ^ 3
        0, 0, 4,  2, 4, 2,  2, 8, 2,  2,12, 2,  4, 4, 0,  4, 8, 0, & ! 7/2 ^ 4
        4,10, 0,  4,16, 0,                                         & 
        1, 9, 4,  3, 3, 2,  3, 5, 2,  3, 7, 2,  3, 9, 2,  3,11, 2, & ! 9/2 ^ 3
        3,13, 2,  3,15, 2,  3,17, 2,  3,21, 2,                     &
        0, 0, 5,  2, 4, 3,  2, 8, 3,  2,12, 3,  2,16, 3,  4, 0, 1, & ! 9/2 ^ 4
        4, 4, 1,  4, 6, 1,  4, 8,10,  4, 8,20,  2,10, 1,  4,12,10, &
        4,12,20,  4,14, 1,  4,16, 1,  4,18, 1,  4,20, 1,  4,24, 1, &
        1, 9, 4,  3, 3, 2,  3, 5, 2,  3, 7, 2,  3, 9, 2,  3,11, 2, & ! 9/2 ^ 5
        3,13, 2,  3,15, 2,  3,17, 2,  3,21, 2,  5, 1, 0,  5, 5, 0, &
        5, 7, 0,  5, 9, 0,  5,11, 0,  5,13, 0,  5,15, 0,  5,17, 0, &
        5,19, 0,  5,25, 0/

!       JW = JQ/10  -> additional seniority

!----------------------------------------------------------------------
! ... check input data:

      if(j.lt.0.or.mod(j,2).eq.0) then
       write(*,*) 'Jterm: unphysical j-value: j =',j
       Stop  'Stop in Jterm'
      end if 

      qm = j + 1
      if(q.lt.0.or.q.gt.qm) then
       write(*,*) 'Jterm: number of electron is out of range: q =',q
       Stop  'Stop in Jterm'
      end if 

      if(j.gt.9.and.q.gt.2) then
       write(*,*) 'Jterm: j^q out of scope: j,q =',j,q
       Stop  'Stop in Jterm'
      end if 

      Jterm = 0 

!----------------------------------------------------------------------
! ... the number of terms in j^q subshell:
      
      if(q.le.1.or.q.ge.qm-1) then
       nterm = 1 
      elseif(q.eq.2.or.q.eq.qm-2) then
       nterm = qm/2
      else
       Select case(j*100 + q)
        Case(503);       nterm = 3; ip = 0    ! 5/2 ^ 3
        Case(703,705);   nterm = 6; ip = 3    ! 7/2 ^ 3,5
        Case(704);       nterm = 8; ip = 9    ! 7/2 ^ 4
        Case(903,907);   nterm =10; ip =17    ! 9/2 ^ 3
        Case(904,906);   nterm =18; ip =27    ! 9/2 ^ 4,6
        Case(905);       nterm =20; ip =45    ! 9/2 ^ 5
        Case default;
         write(*,*) 'Jterm: cannot find j^q subshell for j,q=',j,q 
         Stop 'Stop in Jterm'
       End Select
      end if

      if(k.lt.0) then; Jterm = nterm; Return; end if
!----------------------------------------------------------------------
! ... position of the JT,JV term in the subshell list      

      if(k.eq.0) then          

      qq = min0(q,qm-q)
      Select case(qq)
       case(0);  JV=0
       case(1);  JV=1
       case(2);  JV=2
       case(3);  JV=3; if(j.eq.JT) JV=1
       case(4);  if(j.eq.7) then
                  if(JT.eq.12) JV=2
                  if(JT.eq.10) JV=4
                  if(JT.eq.16) JV=4
                 end if
                 if(j.eq.9.and.JT.ge.18) JV=4
       case(5);  if(j.eq.9) then
                  if(JT.eq. 1) JV=5
                  if(JT.eq. 3) JV=3
                  if(JT.eq.19) JV=5
                  if(JT.eq.21) JV=3
                  if(JT.eq.25) JV=5
                 end if
      End select
      if(JT.eq.0.and.j.lt.9) JV=0

       if(q.le.1.or.q.ge.j) then
        Jterm =  1
       elseif(q.eq.2.or.q.eq.qm-2) then
        Jterm =  JT/4+1
       else
        Do i=1,nterm; ii = ip + i
         if(JV.ne.term_list(1,ii).or.JT.ne.term_list(2,ii)) Cycle
         Jterm=i; Exit
        End do
      end if
 
       if(Jterm.eq.0) then
        write(*,*) 'Jterm:  incorect term j^q(JV,JT):',j,q,JV,JT
        write(*,*) '2j =',j
        write(*,*) 'q  =',q
        write(*,*) '2J =',JT
        write(*,*) 'v  =',JV
        Stop 'Stop in Jterm'
       end if

       Return
      end if

!----------------------------------------------------------------------
! ... k-th term of j^q subshell

      if(k.gt.nterm) then 
       write(*,*) 'Jterm:  incorect input k =',k
       Stop 'Stop in Jterm'
      end if

      JW = 0
      if(q.eq.0.or.q.eq.qm) then                 !  j ^ 0
       JV = 0; JT = 0; JQ = qm/2
      elseif(q.eq.1.or.q.eq.qm-1) then           !  j ^ 1
       JV = 1; JT = j; JQ = qm/2-1
      elseif(q.eq.2.or.q.eq.qm-2) then           !  j ^ 2
       JV = 0; if(k.gt.1) JV = 2 
       JT = (k-1) * 4
       JQ = qm/2; if(k.gt.1) JQ = JQ - 2
      else                                       !  j ^ q
       ii = ip + k
       JV = term_list(1,ii)
       JT = term_list(2,ii)
       JQ = term_list(3,ii)
       if(JQ.ge.10) then; JW = JQ/10; JQ = 1; end if
      end if

      End Function Jterm
