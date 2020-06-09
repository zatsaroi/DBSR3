!======================================================================
      Subroutine Check_jj_shells
!======================================================================
!     check shells
!----------------------------------------------------------------------
      Use jj2ls;  Use conf_jj

      Implicit none
      Integer :: i,icase, JW,JQ
      Integer, external :: Jterm

      i=1; n=0; jot(1)=0; icase=0
    1 n = n + 1
      if(i.eq.no) then
       if(kn(i).gt.0) icase=-1
       if(kn(i).lt.0) icase= 1
      elseif(nn(i).ne.nn(i+1).or.ln(i).ne.ln(i+1)) then    ! ???
!      elseif(ln(i).ne.ln(i+1)) then
       if(kn(i).gt.0) icase=-1
       if(kn(i).lt.0) icase= 1
      else
       if(kn(i).gt.0) icase=-2
       if(kn(i).lt.0) icase= 2
      end if
      ncase(n) = icase

      Select case(icase)
       case(-1)                                    !  (j-,0)
        nn1(n)=nn(i); ln1(n)=ln(i)
        iq1(n)=iq(i); jq1(n)=iq(i); jq2(n)=0  
        kt1(n)=Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ) 
        kt2(n)=1
        jot2(n) = Jshell(i)
        jot3(n) = 0
        jotp(n) = Jintra(i)
        JT(n)   = Jintra(i) 
       case(1)                                     !  (0,j+)
        nn1(n)=nn(i); ln1(n)=ln(i)
        iq1(n)=iq(i); jq1(n)=0; jq2(n)=iq(i)  
        kt2(n)=Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ) 
        kt1(n)=1
        jot3(n) = Jshell(i)
        jot2(n) = 0
        jotp(n) = Jintra(i)
        JT(n)   = Jintra(i) 
       case(-2)                                   !  (j-,j+)
        nn1(n)=nn(i); ln1(n)=ln(i)
        jq1(n)=iq(i)  
        kt1(n)=Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ) 
        jot2(n) = Jshell(i)
        jotp(n) = Jintra(i)
        i=i+1
        jq2(n)=iq(i)
        kt2(n)=Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ) 
        jot3(n) = Jshell(i)
        JT(n)   = Jintra(i) 
       case(2)                                    !  (j+,j-)
        nn1(n)=nn(i); ln1(n)=ln(i)
        jq2(n)=iq(i)  
        kt2(n)=Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ) 
        jot3(n) = Jshell(i)
        jotp(n) = Jintra(i)
        i=i+1
        jq1(n)=iq(i)
        kt1(n)=Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ) 
        jot2(n) = Jshell(i)
        JT(n)   = Jintra(i) 
       Case default
        Stop 'chech_jj_shells: icase - ?'
       End Select

       if(n.gt.1) jot1(n)=JT(n-1)
       iq1(n) = jq1(n)+jq2(n)

       i = i + 1
       if(i.le.no) go to 1

       no1 = n

if(debug.gt.1) then
 write(6,*) 'Jshell',Jshell(1:no)
 write(6,*) 'Jintra',Jintra(1:no)
 write(6,*) 'n',n
 write(6,*) 'ncase',ncase(1:no1)
 write(6,*) 'ln1',ln1(1:no1)
 write(6,*) 'iq1',iq1(1:no1)
 write(6,*) 'jot1',jot1(1:no1)
 write(6,*) 'jot2',jot2(1:no1)
 write(6,*) 'jot3',jot3(1:no1)
 write(6,*) 'jotp',jotp(1:no1)
 write(6,*) 'JT  ',JT(1:no1)
 write(6,*) 'kt1 ',kt1(1:no1)
 write(6,*) 'kt2 ',kt2(1:no1)
 write(6,*) 'jq1',jq1(1:no1)
 write(6,*) 'jq2',jq2(1:no1)
end if

      End  Subroutine check_jj_shells
       

