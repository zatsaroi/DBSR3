!====================================================================
      Module cmdata
!====================================================================
!  Contains a set of coefficients with three identifiers (k1,k2,k3)
!  and order pointer (ipt). The data are devided on different types
!  (ntype, itype) and saved in memory in blocks (nblocks, mblocks).
!--------------------------------------------------------------------
      Implicit none 
    
! ... coefficients:

      Real(8), allocatable :: CLDATA(:), CVDATA(:)
    
! ... their attributes:

      Integer, allocatable :: K1(:),K2(:),K3(:),IPT(:)

! ... number of different structures and current structure 

      Integer, parameter :: ntype =  10   
      Integer :: itype =  0                

! ... the data are divided on blocks, with number and size:

      Integer, parameter :: nblocks =  50     
      Integer, parameter :: mblocks =  2000  

! ... current block for given type:

      Integer, allocatable :: iblk(:)    

! ... pointer on first(last) element in given block 

      Integer, allocatable :: ipblk(:), ipi(:)   
      Integer, allocatable :: jpblk(:), ipj(:)  

! ... number of blocks for given type: 

      Integer, allocatable :: nblk(:)
    
! ... blocks chain pointer for given type:       

      Integer, allocatable :: kblk(:,:)  

! ... total number of ordererd coefficients:

      Integer :: ncdata = 0  
    
      End Module cmdata


!======================================================================
      Subroutine Allocate_cmdata(k)
!======================================================================
!     allocate arrays in module "cmdata"
!----------------------------------------------------------------------
      Use cmdata
   
      Implicit none
      Integer, intent(in) :: k
      Integer :: i,n,m

      if(allocated(CLDATA)) Deallocate(CLDATA,CVDATA,K1,K2,K3,IPT, &
                                       iblk,ipblk,jpblk,ipi,ipj,nblk,kblk)
      if(k.eq.0) Return

      n = nblocks; m = nblocks*mblocks
      Allocate(CLDATA(m),CVDATA(m),K1(m),K2(m),K3(m),IPT(m), &
               iblk(n),nblk(ntype),kblk(ntype,n), &
               ipblk(n),jpblk(n),ipi(n),ipj(n) )

! ... initilize the blocks::

      Do i = 1,nblocks
       ipblk(i) = (i-1)*mblocks + 1; jpblk(i) = -1
      End do

! ... assign one block to each type:

      Do i = 1,ntype
       iblk(i) = i; nblk(i) = 1; kblk(i,1) = i; jpblk(i) = 0
      End do

      End Subroutine Allocate_cmdata


!======================================================================
      Subroutine Add_cdata(ip,jp,CL,CV,j1,j2,j3)
!======================================================================
!     add new coefficients with indexes (j1,j2,j3) to the list (ip:jp)
!----------------------------------------------------------------------
      Use cmdata

      Implicit none
      Integer:: ip,jp, j1,j2,j3, k,l,m, i,ii
      Real(8), intent(in) :: CL,CV

! ... first element in the list:

      if(jp.lt.ip) then
       cldata(ip) = CL; cvdata(ip) = CV;
       k1(ip)=j1; k2(ip)=j2; k3(ip)=j3; jp=ip; Return
      end if       

! ... search position for new integral

      k=ip; l=jp;
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if(j1.lt.K1(m)) then
       l = m - 1
      elseif(j1.gt.K1(m)) then
       k = m + 1
      else
       if(j2.lt.K2(m)) then
        l = m - 1
       elseif(j2.gt.K2(m)) then
        k = m + 1
       else
        if(j3.lt.K3(m)) then
         l = m - 1
        elseif(j3.gt.K3(m)) then
         k = m + 1
        else
         cldata(m) = cldata(m) + CL
         cvdata(m) = cvdata(m) + CV
         Return
        end if
       end if
      end if
      go to 1

    2 Continue 

      Do i=jp,k,-1
       ii = i + 1
       cldata(ii)=cldata(i); cvdata(ii)=cvdata(i)
       K1(ii)=K1(i); K2(ii)=K2(i); K3(ii)=K3(i)
      End do
      CLdata(k)=CL; CVdata(k)=CV;
      K1(k)=j1; K2(k)=j2; K3(k)=j3
      jp = jp + 1

      End Subroutine Add_cdata


!======================================================================
      Subroutine Add_coef(CL,CV,j1,j2,j3)
!======================================================================
!     add new coefficient to the list for given itype
!     Calls: Add_cdata, Gen_matrix
!----------------------------------------------------------------------
      Use cmdata

      Real(8), intent(in) :: CL,CV
      Integer, intent(in) :: j1,j2,j3

! ... add coefficient to list:
      
      i = iblk(itype); ip = ipblk(i); jp = jpblk(i) 
      Call Add_cdata(ip,jp,CL,CV,j1,j2,j3)                    
      jpblk(i) = jp

! ... check if the block full:

      if(jp-ip+1.lt.mblocks) Return

! ... find new block:

      m = 0
      Do i=1,nblocks
       if(jpblk(i).ge.0) Cycle
       iblk(itype) = i; jpblk(i) = 0
       nblk(itype) = nblk(itype) + 1
       kblk(itype,nblk(itype))=i
       m = 1; Exit
      End do
      if(m.ne.0) Return
     
! ... everything is full - it is time to generate matrix for
! ... the current and most extended cases:

      Call Gen_matrix

      m = 0
      Do i=1,ntype
       if(nblk(i).lt.m) Cycle
       m=nblk(i); itype=i
      End do
      Call Gen_matrix

      End Subroutine Add_coef 


!======================================================================
      Subroutine Merge_CDATA(nn, ip, jp, nc, EPS_c)
!======================================================================
!     merge the different blocks of data in MODULE 'cmdata'
!
!     nn    - number of blocks
!     ip(.) - pointer for begining of block .
!     jp(.) - pointer for end of block .
!     nc    - number of result coeff's
!     EPS_c - all coefficients < EPS_c are ignored
!----------------------------------------------------------------------
      Use cmdata, only: CLDATA,CVDATA, K1,K2,K3, IPT

      Implicit none
      Integer :: nn,nc
      Integer :: ip(nn),jp(nn)
      Integer :: i,ii, j,jj, m,mm
      Real(8), intent(in) :: EPS_C

      nc = 0

! ... choose the non-empty block

       mm = 0
       Do m = 1,nn
        i = IP(m); if(i.eq.0) Cycle
        mm = m; Exit
       End do
       if(mm.eq.0) Return

! ...  main loop ...

    1 Continue
                             
! ... compare integrals in different blocks and merge the coefficients
! ... in case of equal integrals

      Do ii=1,nn-1
       i=IP(ii); if(i.eq.0) Cycle
       Do jj=ii+1,nn
        j=IP(jj); if(j.eq.0) Cycle
        if(K1(i).ne.K1(j)) Cycle
        if(K2(i).ne.K2(j)) Cycle
        if(K3(i).ne.K3(j)) Cycle
        CLDATA(i) = CLDATA(i) + CLDATA(j); CLDATA(j) = 0.d0
        CVDATA(i) = CVDATA(i) + CVDATA(j); CVDATA(j) = 0.d0
        mm=jj; go to 2
       End do
      End do

! ...  choose the minimum K1, then K2, then K3 

      j=IP(mm)
      Do m=1,nn
       if(IP(m).eq.0) Cycle; i=IP(m)
       if(K1(i).lt.K1(j)) then
        mm=m; j=IP(mm)
       elseif(K1(i).gt.K1(j)) then
        Cycle
       elseif(K2(i).lt.K2(j)) then
        mm=m; j=IP(mm)
       elseif(K2(i).eq.K2(j).and.K3(i).lt.K3(j)) then
        mm=m; j=IP(mm)
       end if
      End do

! ... mark the chosen coefficient 

      i=IP(mm)
      if(abs(CLDATA(i))+abs(CVDATA(i)).gt.EPS_c) then
       nc=nc+1; IPT(nc)=i
      end if

! ... choose next data

    2 IP(mm) = IP(mm) + 1
      if(IP(mm).le.JP(mm)) go to 1
      if(IP(mm).gt.JP(mm)) then
       IP(mm)=0
       Do m=1,nn
        if(IP(m).gt.0) then; mm=m; go to 1; end if
       End do
      end if

      End Subroutine Merge_CDATA


