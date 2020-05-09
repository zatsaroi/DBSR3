!======================================================================
      Module c_data
!======================================================================
! ... contains a set of coefficients with four identifiers 
! ... (k1,k2,k3,k4) and one pointer (ipt), devided in blocks;
! ... each block contain data for given multipole and type   
!
! ... internal procedures: Alloc_c_data (mm)
!                          Add_coef     (kpol,itype,C,j1,j2,j3,j4)
!                          Add_cdata    (ip,jp,C,j1,j2,j3,j4)
!                          Add_matrix   (jtype,jpol)
!                          Merge_cdata  (nn, ip, jp, nc)
! ... external calls:      RK_evaluate  (jtype,jpol) 
!----------------------------------------------------------------------
      Implicit none 

      Integer :: ncdata = 0               ! number of coefficients:
      Real(8), allocatable :: CDATA(:)    ! their values 
                                          ! their attributes:
      Integer, allocatable :: K1(:),K2(:),K3(:),K4(:),IPT(:)

! ... default dimension limits:

      Integer :: mk    =      7     ! max. multipole index
      Integer :: nb    =   5000     ! number of blocks in c_data       
      Integer :: mb    =   1000     ! size of blocks
      Integer :: kb    =    100     ! max.nb for given case

      Real(8) :: eps_c    = 1.0D-10 !  tolerance for coefficients	                               

! ... pointer on first(last) element in given block 

      Integer, allocatable :: ipblk(:), ipi(:)   
      Integer, allocatable :: jpblk(:), ipj(:)  

! ... current block for given kpol and itype:

      Integer, allocatable :: iblk(:,:)    

! ... number of blocks for given kpol and itype: 

      Integer, allocatable :: nblk(:,:)
    
! ... chain ponter for given kpol and itype:       

      Integer, allocatable :: kblk(:,:,:)

! ... structure of data:

      Integer, parameter :: ntype =  6  ! 4-Rk and 2-Sk integrals
      Integer, parameter :: ibi = 2**15 ! packing basis for orbitals
      Integer :: kpol1  =  0            !  min. multipole index                   
      Integer :: kpol2  =  0            !  max. multipole index                   
    
      End Module c_data


!======================================================================
      Subroutine Alloc_c_data(mm)
!======================================================================
!     allocate arrays in module c_data
!----------------------------------------------------------------------
      Use c_data

      Implicit none
      Integer :: m,mm,i,j,k

      if(allocated(CDATA)) Deallocate(CDATA,K1,K2,K3,K4,IPT, &
                          ipblk,jpblk,ipi,ipj,iblk,nblk,kblk)
      if(mm.eq.0) Return
      m = mb*nb
      Allocate(CDATA(m),K1(m),K2(m),K3(m),K4(m),IPT(m), &
               ipblk(nb),jpblk(nb),ipi(nb),ipj(nb), &
               iblk(kpol1:kpol2,ntype),nblk(kpol1:kpol2,ntype), &
               kblk(kpol1:kpol2,ntype,nb) )

! ... initilize all blocks:

      Do i=1,nb;  ipblk(i)=(i-1)*mb+1; jpblk(i)=-1; End do

! ... assign one block to each type:

      i = 0
      Do k = kpol1,kpol2
       Do j = 1,ntype
        i = i + 1
        iblk(k,j)=i; nblk(k,j)=1; kblk(k,j,1)=i; jpblk(i)=0
       End do
      End do

      if(i.gt.nb/2) &
       Stop 'Stop in alloc_c_data: increase nb parameter in c_data'
	        
      End Subroutine Alloc_c_data


!======================================================================
      Subroutine Add_coef(kpol,itype,C,j1,j2,j3,j4)
!======================================================================
!     add new coefficient to the c_data list for given kpol and itype
!----------------------------------------------------------------------
      Use c_data

      Real(8) :: C
      Integer :: j1,j2,j3,j4,kpol,itype

! ... add coefficient to list:

      i = iblk(kpol,itype); ip = ipblk(i); jp = jpblk(i) 
      Call Add_cdata(ip,jp,C,j1,j2,j3,j4)                    
      jpblk(i) = jp

! ... check if the block full:

      if(jp-ip+1.lt.mb) Return

! ... find new block:

      if(nblk(kpol,itype).lt.kb) then
      m = 0
      Do i=1,nb
       if(jpblk(i).ge.0) Cycle
       iblk(kpol,itype) = i; jpblk(i) = 0
       nblk(kpol,itype) = nblk(kpol,itype) + 1
       kblk(kpol,itype,nblk(kpol,itype))=i
       m = 1; Exit
      End do
      if(m.ne.0) Return
      end if
     
! ... everything is full - it is time to generate matrix:

      jtype=itype; jpol=kpol; n = nblk(kpol,itype)
      Do i = 1,ntype; Do k = kpol1,kpol2
       if(nblk(k,i).le.n) Cycle
       jtype = i; jpol = k
      End do; End do
      Call Add_matrix(jtype,jpol)
      if(jtype.eq.itype.and.jpol.eq.kpol) Return

! ... find again new block:

      m = 0
      Do i=1,nb
       if(jpblk(i).ge.0) Cycle
       iblk(kpol,itype) = i; jpblk(i) = 0
       nblk(kpol,itype) = nblk(kpol,itype) + 1
       kblk(kpol,itype,nblk(kpol,itype))=i
       m = 1; Exit
      End do
      if(m.eq.0) Stop ' Add_coef: problems with new block '
      
      End Subroutine Add_coef 


!======================================================================
      Subroutine Add_cdata(ip,jp,C,j1,j2,j3,j4)
!======================================================================
!     add new data to the sub-list (ip:jp) in module c_data
!----------------------------------------------------------------------
      Use c_data

      Implicit none
      Integer:: ip,jp, j1,j2,j3,j4, k,l,m, i,ii
      Real(8), intent(in) :: C

! ... case of first element in the list:

      if(jp.lt.ip) then
       cdata(ip)=C; k1(ip)=j1; k2(ip)=j2; k3(ip)=j3; k4(ip)=j4
       jp=ip
	      Return
      end if       

! ... search position (k) for new integral

      k=ip; l=jp;
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (j1.lt.K1(m)) then;       l = m - 1
      elseif(j1.gt.K1(m)) then;       k = m + 1
      else
       if    (j2.lt.K2(m)) then;      l = m - 1
       elseif(j2.gt.K2(m)) then;      k = m + 1
       else
        if    (j3.lt.K3(m)) then;     l = m - 1
        elseif(j3.gt.K3(m)) then;     k = m + 1
        else
         if    (j4.lt.K4(m)) then;    l = m - 1
         elseif(j4.gt.K4(m)) then;    k = m + 1
         else
          cdata(m)=cdata(m)+C; Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue 

! ... shift the rest of data up:

      Do i=jp,k,-1
       ii = i + 1
       cdata(ii)=cdata(i)
       K1(ii)=K1(i); K2(ii)=K2(i); K3(ii)=K3(i); K4(ii)=K4(i)
      End do

! ... add new integral:

      Cdata(k)=C; K1(k)=j1; K2(k)=j2; K3(k)=j3; K4(k)=j4; jp=jp+1

      End Subroutine Add_cdata


!======================================================================
      Subroutine Add_matrix(jtype,jpol)
!======================================================================
!     merge data and generate the interaction matrix for given type
!     jtype and mutipole index jpol
!     Call(s):  RK_data1 
!----------------------------------------------------------------------
      Use c_data, nc => ncdata

      Implicit none
      Integer, Intent(in) :: jtype,jpol
      Integer :: i,j,k,n

! ... prepare the data:
       
      n = nblk(jpol,jtype)          ! number of blocks
      i = kblk(jpol,jtype,1)        ! first block
      if(n.eq.1.and.jpblk(i).le.0) n=0

      if(n.eq.0) then          ! nothing to do
       Return  

      elseif(n.eq.1) then      ! simple case - one block

       nc = jpblk(i)-ipblk(i) + 1
       j = ipblk(i)-1
       Do i = 1,nc;  IPT(i)=j+i;  End do

      else                      ! need merging data from 
                                ! different blocks
       Do i = 1,n               
        j = kblk(jpol,jtype,i)
        ipi(i) = ipblk(j)  
        ipj(i) = jpblk(j)
       End do

       Call Merge_cdata(n, ipi, ipj, nc)

! ... release the blocks:       

       Do i=1,n; jpblk(kblk(jpol,jtype,i))=-1;  End do

      end if

! ... re-assign the first block for given jtype and jpol:

      i = kblk(jpol,jtype,1); j=jtype; k=jpol
      iblk(k,j)=i; nblk(k,j)=1; jpblk(i)=0

! ... generate matrix:

      Call Rk_evaluate(jpol,jtype) 

      End Subroutine Add_matrix


!======================================================================
      Subroutine Merge_cdata(nn, ip, jp, nc)
!======================================================================
!     merge the different blocks of data in MODULE 'cmdata'
!
!     nn    - number of blocks
!     ip(.) - pointer for begining of block .
!     jp(.) - pointer for end of block .
!     nc    - number of result coeff's
!     EPS_c - all coefficients < EPS_c are ignored
!----------------------------------------------------------------------
      Use c_data 

      Implicit none
      Integer :: nn,nc,ip(nn),jp(nn),i,ii, j,jj, m,mm

      nc = 0

! ... choose the non-empty block

      mm=0;  Do m=1,nn; if(IP(m).eq.0) Cycle; mm=m; Exit;  End do
      if(mm.eq.0) Return

! ... main loop ...

    1 Continue
                             
! ... compare integrals in different blocks and merge the coefficients
! ... in case of equal integrals

      Do ii=1,nn-1;  i=IP(ii); if(i.eq.0) Cycle
       Do jj=ii+1,nn;  j=IP(jj); if(j.eq.0) Cycle
        if(K1(i).ne.K1(j)) Cycle
        if(K2(i).ne.K2(j)) Cycle
        if(K3(i).ne.K3(j)) Cycle
        if(K4(i).ne.K4(j)) Cycle
        CDATA(i) = CDATA(i) + CDATA(j); CDATA(j) = 0.d0
        mm=jj; go to 2
       End do
      End do

! ... choose the minimum K1, then K2, then K3, then K4 

      j=IP(mm)
      Do m=1,nn
       if(IP(m).eq.0) Cycle; i=IP(m)
       if    (K1(i).lt.K1(j)) then;   mm=m; j=IP(mm)
       elseif(K1(i).gt.K1(j)) then;   Cycle
       else
        if    (K2(i).lt.K2(j)) then;  mm=m; j=IP(mm)
        elseif(K2(i).gt.K2(j)) then;  Cycle
        else
         if    (K3(i).lt.K3(j)) then; mm=m; j=IP(mm)
         elseif(K3(i).gt.K3(j)) then; Cycle
         elseif(K4(i).lt.K4(j)) then; mm=m; j=IP(mm)
         end if
        end if
       end if
      End do

! ... mark the chosen coefficient 

      if(abs(CDATA(IP(mm))).gt.Eps_C) then
        nc=nc+1; IPT(nc)=IP(mm)
      end if

! ... choose next data

    2 IP(mm) = IP(mm) + 1
      if(IP(mm).le.JP(mm)) go to 1
      if(IP(mm).gt.JP(mm)) then
       IP(mm)=0
       Do m=1,nn; if(IP(m).gt.0) then; mm=m; go to 1; end if; End do
      end if

      End Subroutine Merge_CDATA

