!======================================================================
      Subroutine Decode_cjj(CONFIG,SHELLJ,INTRAJ,no,nn,kn,ln,jn,iq,in,&
                            Jshell,Vshell,Jintra)
!======================================================================
!     decode the configuration from c-file format to INTEGER format
!     Call: EL_nljk
!----------------------------------------------------------------------
 
      Implicit none
      Character(*), Intent(in) :: CONFIG,SHELLJ,INTRAJ
      Integer, Intent(out) :: no,nn(*),kn(*),ln(*),jn(*),iq(*),in(*),&
                              Jshell(*),Vshell(*),Jintra(*)
      Integer :: i,j,k,m 

      m=INDEX(CONFIG,')',BACK=.TRUE.); no=m/9
 
      Vshell(1:no)=0; Jshell(1:no)=0; Jintra(1:no)=0
      m=-9
      Do i=1,no
       m = m + 9

       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)
 
       if(SHELLJ(m+5:m+5).eq.';') then
        read(SHELLJ(m+4:m+4),'(i1)') Vshell(i)
        read(SHELLJ(m+9:m+9),'(i1)') J; Jshell(i) = 2*J
       elseif(SHELLJ(m+8:m+8).eq.'/') then
        read(SHELLJ(m+1:m+7),'(i7)') Jshell(i)
       elseif(SHELLJ(m+9:m+9).ne.' ') then
        read(SHELLJ(m+1:m+9),'(i9)') J; Jshell(i) = 2*J
       end if       

       k = INDEX(INTRAJ(m+1:m+9),'/')
       if(i.eq.1) then
        Jintra(i) = Jshell(i)
       elseif(i.eq.no) then
        Cycle
       elseif(k.gt.0) then
        read(INTRAJ(m+1:m+k-1),*) Jintra(i)
       elseif(LEN_TRIM(INTRAJ(m+1:m+9)).ne.0) then
        read(INTRAJ(m+1:m+9),*) J; Jintra(i) = 2*J
       else
       if(Jshell(i).eq.0) Jintra(i) = Jintra(i-1)              
       if(Jintra(i-1).eq.0) Jintra(i) = Jshell(i)
       end if       

      End do

      if(k.gt.0) then
       read(INTRAJ(m+1:m+k-1),*) Jintra(no)
      else
       k=INDEX(INTRAJ(m+1:m+9),'+')
       if(k.eq.0) k=INDEX(INTRAJ(m+1:m+9),'-')
       read(INTRAJ(m+1:m+k-1),*) J; Jintra(no) = 2*J
      end if

      Do i=1,no
       m = jn(i)+1
       if(Vshell(i).ne.0.or.iq(i).eq.m) Cycle
       if(iq(i).eq.1.or.iq(i).eq.m-1) then
         Vshell(i) = 1
       elseif(iq(i).eq.2.or.iq(i).eq.m-2) then
         if(Jshell(i).gt.0) Vshell(i) = 2
       elseif(iq(i).eq.3) then
         Vshell(i) = 3; if(Jshell(i).eq.jn(i)) Vshell(i)=1
       end if
      End do

      End Subroutine Decode_cjj
