!======================================================================
      Subroutine  D_data 
!======================================================================
!     processing of D-integrals in the module 'cmdata'
!----------------------------------------------------------------------
!
!    we have following 10 different structures for radial integrals:
!
! 1  d( . .)  ic, jc               -  bound-bound  

! 2  d( i .)  jc                   -  bound-channel
! 3  d( . j)  ic                    
! 4  d( . .) < i | . > jc           
! 5  d( . .) < . | j > ic           

! 6  d( i j)                       -  channel-channel
! 7  d( i .) < . | j >              
! 8  d( . j) < i | . >               
! 9  d( . .) < i | . > < . | j >    
!10  d( . .) < i | j >              
!
! where . denotes bound orbital, i,j - channels,
! ic,jc - configurations.
!
!----------------------------------------------------------------------
      Use dbsr_dmat  
      Use cmdata
      Use DBS_orbitals_pq, only: ipbs

      Implicit none
      Integer :: i,j, i1,i2, ich,jch, ic,jc, ib,jb, io,jo
      Real(8) :: v(ms),v1(ms),v2(ms), w(ms),w1(ms),w2(ms)

!if(itype.ne.9) Return

!write(*,*) 'd_data: itype,nc',itype,ncdata                 

      Select Case(itype)

      Case(1)                            !  d( . .)  ic, jc                             
!write(*,*) 'd( . .)  ic, jc '
       Do j=1,ncdata;  i=IPT(j)
        Call Update_DB(k1(i),k2(i),CLDATA(i),CVDATA(i))
!write(*,*)  k1(i),k2(i),CLDATA(i),CVDATA(i)
       End do

      Case(2)                            !  d( i .)  jc                             
!write(*,*) 'd( i .)  jc'
       Do j=1,ncdata;  i=IPT(j)
        i1=k1(i); i2=k2(i); ich=ipbs(i1); jch=0; ic=0; jc=k3(i)
        Call Get_Vdip(i1,i2,kpol,'r',v1,v2) 
        Call UPDATE_DV(ich,jch,ic,jc,v1,v2,cldata(i))
!write(*,*) ich,jch,ic,jc,cldata(i)
       End do

      Case(3)                             !  d( . j)  ic                              
!write(*,*) 'd( . j)  ic'  
       Do j=1,ncdata;  i=IPT(j)   
        i1=k1(i); i2=k2(i); ich=0; jch=ipbs(i2); ic=k3(i); jc=0
        Call Get_Vdip(i1,i2,kpol,'l',v1,v2) 
        Call UPDATE_DV(ich,jch,ic,jc,v1,v2,cldata(i))
!write(*,*) ich,jch,ic,jc,cldata(i)
       End do

      Case(4)                           !   d( . .) < i | . > jc                            
!write(*,*) 'd( . .) < i | . > jc '
       Do j=1,ncdata;  i=IPT(j)
        io=k1(i); ich=ipbs(io/ibo); ib=mod(io,ibo); jch=0 
        ic=0; jc=k2(i)
!write(*,*) 'ic,jc,ich,jch',ic,jc,ich,jch,ebs(ib),cldata(i),cvdata(i)
        Call Get_qv(ib,v,ns)
        v1=v*cldata(i)
        v2=v*cvdata(i)
        Call UPDATE_DV(ich,jch,ic,jc,v1,v2,1.d0)
       End do

      Case(5)                           !   d( . .) < . | j > ic                   
!write(*,*) 'd( . .) < . | j > ic '                                    
                                                                      
       Do j=1,ncdata;  i=IPT(j)                                       
        io=k1(i); ich=0; jb=mod(io,ibo); jch=ipbs(io/ibo)     ! jch=ipbs(mod(io,ibo)); jb=io/ibo      ???  
        ic=k2(i); jc=0                                                
!write(*,*) k1(i),k2(i), mod(io,ibo), jb
!write(*,*) 'ic,jc,ich,jch',ic,jc,ich,jch,ebs(jb),cldata(i),cvdata(i)  
        Call Get_qv(jb,v,ns)                                          
        v1 = v * cldata(i)                                            
        v2 = v * cvdata(i)                                            
        Call UPDATE_DV(ich,jch,ic,jc,v1,v2,1.d0)
       End do

      Case(6)                           !   d( i j)                            
!write(*,*) 'd( i j)'
       Do j=1,ncdata;  i=IPT(j)
        ich=ipbs(k1(i)); jch=ipbs(k2(i))
!write(*,*) ich,jch,cldata(i)
        Call UPDATE_DX(ich,jch,cldata(i))
       End do

      Case(7)                            !  d( i .) < . | j >     
!write(*,*) 'd( i .) < . | j >'
       Do j=1,ncdata;  i=IPT(j)
        i1=k1(i); i2=k2(i); ich=ipbs(i1) 
        Call Get_Vdip(i1,i2,kpol,'r',v1,v2) 

        io=k3(i); jb=mod(io,ibo); jch=ipbs(io/ibo)   ! jch=ipbs(mod(io,ibo)); jb=io/ibo

!write(*,*) i1,i2,io,mod(io,ibo), io/ibo                         ! ???  

        Call Get_qv(jb,w1,ns)
        Call Get_qv(jb,w2,ns)
        Call UPDATE_DW(ich,jch,v1,w1,v2,w2,cldata(i))

       End do

      Case(8)                            !  d( . j) < i | . >     
!write(*,*) 'd( . j) < i | . >'
       Do j=1,ncdata;  i=IPT(j) !; Cycle
        i1=k1(i); i2=k2(i); jch=ipbs(i2) 
        Call Get_Vdip(i1,i2,kpol,'l',v1,v2) 
        io=k3(i); ich=ipbs(io/ibo); ib=mod(io,ibo) 
        Call Get_qv(ib,w1,ns)
        Call Get_qv(ib,w2,ns)
        Call UPDATE_DW(ich,jch,w1,v1,w2,v2,cldata(i))

!write(*,*) 'k1(i),k2(i),k3(i),ich,jch',k1(i),k2(i),k3(i),ich,jch,cldata(i)
       End do

      Case(9)                            !  d( . .) < i | . > < . | j >                       

!write(*,*) 'd( . .) < i | . > < . | j > '
  
       Do j=1,ncdata;  i=IPT(j) !; Cycle
        io=min(k1(i),k2(i)); ich=ipbs(io/ibo); ib=mod(io,ibo)
        jo=max(k1(i),k2(i)); jch=ipbs(jo/ibo); jb=mod(jo,ibo)

!        jo=k2(i); jch=ipbs(mod(jo,ibo)); jb=jo/ibo                ???

        Call Get_qv(ib,v,ns)
        Call Get_qv(jb,w,ns)
        v1 = cldata(i)*v;  w1 = w
        v2 = cvdata(i)*v;  w2 = w
        Call UPDATE_DW(ich,jch,v1,w1,v2,w2,1.d0)         
!write(*,*) 'ich,ib,jch,jb',ich,ib,jch,jb,cldata(i),cvdata(j)
       End do

      Case(10)                           !   d( . .) < i | j >                            
!write(*,*) 'd( . .) < i | j >'
       Do j=1,ncdata;  i=IPT(j)!; Cycle
        io=k1(i); i1=io/ibo;  i2=mod(io,ibo) 
        ich=ipbs(min(i1,i2)); jch=ipbs(max(i1,i2))
        Call UPDATE_DS(ich,jch,cldata(i),cvdata(i))
!write(*,*) 'ich,jch,cl,cv',ich,jch,cldata(i),cvdata(i)
       End do

      Case Default

       Stop ' D_data: uknown itype '

      End Select

      
      End Subroutine  D_data

