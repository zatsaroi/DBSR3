!======================================================================
      Subroutine Read_rarg(name,rvalue)
!======================================================================
!     read real argument as name=... from the command line 
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name
      Real(8) :: rvalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS

      iarg = command_argument_count(); if(iarg.eq.0) Return 
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GET_COMMAND_ARGUMENT(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),*) rvalue; Exit
      End do

      End Subroutine Read_rarg


!======================================================================
      Subroutine Read_iarg(name,ivalue)
!======================================================================
!     read integer argument as name=... from the command line 
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name
      Integer :: ivalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS

      iarg = command_argument_count(); if(iarg.eq.0) Return 
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GET_COMMAND_ARGUMENT(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),*) ivalue; Exit
      End do

      End Subroutine Read_iarg


!======================================================================
      Subroutine Read_aarg(name,avalue)
!======================================================================
!     read character argument as name=... from the command line 
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name, avalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS

      iarg = command_argument_count(); if(iarg.eq.0) Return 
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GET_COMMAND_ARGUMENT(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),'(a)') avalue; Exit
      End do

      End Subroutine Read_aarg


!======================================================================
      Subroutine Read_iarr(name,na,iarr)
!======================================================================
!     read integer arrray from string:  name=a,b,c-d,f...
!----------------------------------------------------------------------
      Implicit None
      Character(*) :: name
      Integer :: na,iarr(na)
      Integer :: iarg,ia,iname,i,i1,i2,j,j1,j2,k,k1,k2
      Character(180) :: AS

      iarg = command_argument_count(); if(iarg.eq.0) Return 
      iname=LEN_TRIM(name)
      k=0
      Do i=1,iarg
       Call GET_COMMAND_ARGUMENT(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       k=1; Exit
      End do
      if(k.eq.0) Return

      ia=0; j1=i1; ! iarr=0
      Do 
       j2=INDEX(AS(j1:i2),',')
       if(j2.eq.0) then; j2=i2; else; j2=j2+j1-1; end if
       j=0 ! INDEX(AS(j1:j2),'-'); k = j+j1-1
       if(j.eq.0) then
        ia=ia+1; if(ia.gt.na) Stop 'Read_iarr: ia > na'
        read(AS(j1:j2),*) iarr(ia)
       else
        read(AS(j1:k-1),*) k1
        read(AS(k+1:j2-1),*) k2
        Do k=k1,k2
         ia=ia+1; if(ia.gt.na) Stop 'Read_iarr: ia > na'
         iarr(ia)=k
        End do
       end if
       j1=j2+1; 
       if(j1.gt.i2) Exit
      End do

      End Subroutine Read_iarr


!======================================================================
      Subroutine Read_name(name)
!======================================================================
!     read "name" from the command line as first argument without "="
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name
      Integer :: iarg,i
      Character(80) :: AS

      iarg = command_argument_count(); if(iarg.eq.0) Return 
      Do i=1,iarg
       Call GET_COMMAND_ARGUMENT(i,AS)
       if(INDEX(AS,'=').ne.0) Cycle
       name=AS
       Exit
      End do

      End Subroutine Read_name


!======================================================================
      Subroutine Read_iarr_string(AS,ia,iarr)
!======================================================================
!     read integer arrray from string with  a,b,c-d,f,...
!----------------------------------------------------------------------
      Implicit none

      Character(*) :: AS
      Integer :: ia,iarr(*)
      Integer :: ii,j,j1,j2,k,k1,k2

      ia=0; ii=LEN_TRIM(AS); j1=1
      Do 
       j2=INDEX(AS(j1:ii),',')
       if(j2.eq.0) then; j2=ii; else; j2=j1+j2-1; end if
       j=INDEX(AS(j1:j2),'-'); k = j1+j-1
       if(j.eq.0) then
        ia=ia+1; read(AS(j1:j2),*) iarr(ia)
       else
        read(AS(j1:k-1),*) k1
        read(AS(k+1:j2),*) k2
        Do k=k1,k2
         ia=ia+1; iarr(ia)=k
        End do
       end if
       j1=j2+1; 
       if(j1.gt.ii) Exit
      End do

      End Subroutine Read_iarr_string
