!======================================================================
      Subroutine Read_ang_coef
!======================================================================
!     read angular coefficients form int_bnk file
!     and put them in rk4_data module in ordered and packed form;
!     create list of integrals (in module dbsr_mchf).
!----------------------------------------------------------------------
      Use dbsr_mchf      
      Use c_data

      Implicit none
      Integer, allocatable :: ii_int(:), jj_int(:)
      Real(8), allocatable :: cc(:)
      Integer :: i,j,k,kpol, nc,ii,jj, ishift,jshift, nu
      Real(8) :: t1,t2, S

      Call CPU_time(t1) 

! ... L-integrals:

      Call Read_int_bnk(1)

      Call Add_matrix(1,0)

      Lcoef = ncdata      
      Allocate( ic_Lcoef(Lcoef), jc_Lcoef(Lcoef), L_coef(Lcoef) )
      Lint = 0
      ii = 0; jj = 0
      Do j=1,ncdata; i=IPT(j)
       ic_Lcoef(j) = k3(i)
       jc_Lcoef(j) = k4(i)
       L_coef(j) = Cdata(i)
       if(k1(i).eq.ii.and.k2(i).eq.jj) Cycle
       Lint = Lint + 1
       ii = k1(i); jj = k2(i)
      End do
      
      Allocate(i1_Lint(Lint),i2_Lint(Lint),ip_Lint(0:Lint),L_int(Lint))
      ip_Lint(0) = 0; L_int = 0.d0
      k=0; ii = 0; jj = 0; ip_Lint = 0     
      Do j=1,ncdata; i=IPT(j)
       if(k1(i).eq.ii.and.k2(i).eq.jj) then
        ip_Lint(k)=j
       else
        k = k + 1
        i1_Lint(k) = k1(i); ii = k1(i)
        i2_Lint(k) = k2(i); jj = k2(i) 
        ip_Lint(k)=j 
       end if 
      End do

!----------------------------------------------------------------------
! ... Rk-integrals:

      Call Read_int_bnk(2)

      Call nc_c_data(nc,S)

      Allocate(ii_int(nc),jj_int(nc))
      Allocate(ic_coef(nc),jc_coef(nc),Rk_coef(nc),ip_int(0:nc))
      ip_int(0) = 0
      Allocate(nk_coef(kmin-1:kmax),); nk_coef = 0
      Allocate(nk_int (kmin-1:kmax),); nk_int  = 0

      ishift = 0
      jshift = 0
      Do kpol=kmin,kmax

       Call Add_matrix(1,kpol); if(ncdata.eq.0) Cycle
  
       nk_coef(kpol) = ncdata      

       nint = 0; ii = 0; jj = 0
       Do j=1,ncdata; i=IPT(j)
        ic_coef(j+jshift) = k3(i)
        jc_coef(j+jshift) = k4(i)
        Rk_coef(j+jshift) = Cdata(i)
        if(k1(i).eq.ii.and.k2(i).eq.jj) Cycle
        nint = nint + 1;  ii = k1(i); jj = k2(i)
       End do
       nk_int(kpol) = nint

       k=0; ii = 0; jj = 0     
       Do j=1,ncdata; i=IPT(j)
        if(k1(i).eq.ii.and.k2(i).eq.jj) then
         ip_int(k+ishift) = j + jshift
        else
         k = k + 1
         ii_int(k+ishift) = k1(i); ii = k1(i)
         jj_int(k+ishift) = k2(i); jj = k2(i) 
         ip_int(k+ishift) = j + jshift 
        end if 
       End do
  
       ishift = ishift + nint
       jshift = jshift + ncdata

      End do  ! over kpol

      Call Alloc_c_data(0,0,0,mblock,nblock,eps_c)

      ncoef = SUM(nk_coef); 
      Do k=kmin,kmax; nk_coef(k)=nk_coef(k-1)+nk_coef(k);End do
      nint  = SUM(nk_int)
      Do k=kmin,kmax; nk_int(k)=nk_int(k-1)+nk_int(k);End do

      Allocate(Rk_int(nint)); Rk_int=0.d0

      Allocate(i1_int(nint),i2_int(nint),i3_int(nint),i4_int(nint))
      Do i=1,nint
       i1_int(i) = ii_int(i)/ibi; i3_int(i) = mod(ii_int(i),ibi) 
       i2_int(i) = jj_int(i)/ibi; i4_int(i) = mod(jj_int(i),ibi)
      End do

! ... re-allocate arrays if many empty space
      
      if(nc.gt.ncoef + ncoef/5) then
       ii_int = ic_coef
       jj_int = jc_coef
       Deallocate(ic_coef,jc_coef)
       Allocate(ic_coef(ncoef),jc_coef(ncoef))
       ic_coef = ii_int
       jc_coef = jj_int
       Allocate(cc(nc))
       cc = Rk_coef
       Deallocate(Rk_coef); Allocate(Rk_coef(ncoef))
       Rk_coef(1:ncoef) = cc(1:ncoef)
       Deallocate(cc)
      end if

      Deallocate(ii_int,jj_int)

      Call CPU_time(t2) 
      time_read_coef = t2-t1

!----------------------------------------------------------------------
! ... Sk-integrals:

!      if(mbreit.eq.0) Return

!----------------------------------------------------------------------
! ... debug printing of integrals:

      if(debug.gt.0) then
       Call Find_free_unit(nu)
       AF = trim(name)//'.int'
       open(nu,file=AF)
       Call Print_int_list(nu) 
       Close(nu)
      end if

! ... debug printing of coefficients:

      if(debug.gt.0) then
       Call Find_free_unit(nu)
       AF = trim(name)//'.coef'
       open(nu,file=AF)
       Call Print_coeffs(nu) 
       Close(nu)
      end if

      End Subroutine read_ang_coef


!======================================================================
      Subroutine Read_int_bnk(icase)
!======================================================================
!     read angular coefficients form int_bnk file
!     and put them in rk4_data module in ordered and packed form;
!     create list of integrals (in module dbsr_mchf).
!----------------------------------------------------------------------
      Use dbsr_mchf,  ibdd => ibd,  ibff => ibf      
      Use orb_jj;  Use det_list;  Use def_list
      Use c_data

      Implicit real(8) (A-H,O-Z)
      Real(8) :: S(8) 
      Real(8), allocatable :: Cbuf(:)
      Integer, allocatable :: itb(:),jtb(:),intb(:),idfb(:)

      Allocate(Cbuf(maxnc),itb(maxnc),jtb(maxnc),intb(maxnc),idfb(maxnc))
!----------------------------------------------------------------------
! ... read overlap factors:

      rewind(nub) 
      Call Read_symc(nub)         ! skip ???
      Call Read_symt(nub)
      Call Read_done(nub)
      Call Read_det(nub)
      Call Read_def(nub)

      Select case(icase)
      Case(1) 
       Call Alloc_c_data(1,0,0,mblock,nblock,eps_c)
      Case(2)
       Call Alloc_c_data(1,kmin,kmax,mblock,nblock,eps_c)
      Case(3)
       Call Alloc_c_data(1,kmin,kmax,mblock,nblock,eps_c)
      End Select

! ... processing the data:

    1 ncbuf = 0
      Do i = 1,maxnc
       read(nub,end=2) Cbuf(i),itb(i),jtb(i),intb(i),idfb(i)
       ncbuf = ncbuf + 1
      End do
    2 Continue

      Do ibuf=1,ncbuf; int=intb(ibuf)

       Call Decode_int (jcase,kpol,I1,I2,I3,I4,int)
       if(jcase.ne.icase) Cycle

! ... determine the range of states for given coeff. 

       it = itb(ibuf);  jt = jtb(ibuf)
       is1 = IT_state1(it); js1 = IT_state1(jt)
       if(is1.eq.0.or.js1.eq.0) Cycle
       is2 = IT_state2(it); js2 = IT_state2(jt)
       if(is2.eq.0.or.js2.eq.0) Cycle
    
       idf = idfb(ibuf); nd=0; ip=0; NP=0
       if(idf.gt.0) then
        nd=KPF(idf); ip=IPF(idf); NP(1:nd)=NPF(ip+1:ip+nd)
       end if
    
       C = CBUF(ibuf)

!----------------------------------------------------------------------
! ... cycle over CSF:

       Do ik=is1,is2; ic=IS_order(ik); ip1=IP_state(ic)
     
       Do jk=js1,js2; jc=IS_order(jk); ip2=IP_state(jc)

! ... consider only low-half of interaction matrix:

       if(it.eq.jt.and.ic.lt.jc) Cycle  

       ih = max0(ic,jc)
       jh = min0(ic,jc)

       C = Cbuf(ibuf); ! if(ih.ne.jh) C = C + C    ???

!----------------------------------------------------------------------
! ... check the overlap factor assuming all orbitals are orthogonal:

       DF = 1.d0
       if(idf.gt.0) then
        Do ii=1,nd
         id=NP(ii)/ibf; kd=KPD(id); ip=IPD(id)
         Do jj=1,kd
          k=NPD(jj+ip)
          jp1=k/ibd +    ip1; NP1(jj)=IP_orb(jp1)
          jp2=mod(k,ibd)+ip2; NP2(jj)=IP_orb(jp2)
         End do
         DF = DF * VDET(kd,NP1,NP2)**mod(NP(ii),ibf)
         if(abs(DF).lt.eps_det) Exit
        End do
        if(abs(DF).lt.eps_det) Cycle
       end if 

       CC = C * DF; if(abs(CC).lt.eps_c) Cycle

!----------------------------------------------------------------------
! ... find integral:

       m1=IP_orb(i1+ip1); m2=IP_orb(i2+ip1)
       m3=IP_orb(i3+ip2); m4=IP_orb(i4+ip2)

       Select case(icase)
       Case(1)
        ii=min(m1,m3); jj=max(m1,m3)
        Call Add_coef(CC,kpol,ii,jj,ih,jh,1)
       Case(2)
        Call Jsym_Rk(m1,m2,m3,m4,j1,j2,j3,j4)
        l1=lef(j1);l2=lef(j2);l3=lef(j3);l4=lef(j4) 
        if(mod(l1+l3+kpol,2) + mod(l2+l4+kpol,2) .ne. 0) Cycle   ! ???
        ii = j1*ibi+j3; jj = j2*ibi+j4
        Call Add_coef(CC,kpol,ii,jj,ih,jh,1)
       Case(3)
        Do kk = kpol-1,kpol+1
         if(SMU(kef(m1),kef(m2),kef(m3),kef(m4),kpol,kk,S).eq.0.d0) Cycle
         S = S * CC
         ii = m1*ibi+m3; jj = m2*ibi+m4; Call Add_coef(S(1),kk,ii,jj,ih,jh,1)
         ii = m2*ibi+m4; jj = m1*ibi+m3; Call Add_coef(S(2),kk,ii,jj,ih,jh,1)
         ii = m3*ibi+m1; jj = m4*ibi+m2; Call Add_coef(S(3),kk,ii,jj,ih,jh,1)
         ii = m4*ibi+m2; jj = m3*ibi+m1; Call Add_coef(S(4),kk,ii,jj,ih,jh,1)
         ii = m1*ibi+m3; jj = m4*ibi+m2; Call Add_coef(S(5),kk,ii,jj,ih,jh,1)
         ii = m4*ibi+m2; jj = m1*ibi+m3; Call Add_coef(S(6),kk,ii,jj,ih,jh,1)
         ii = m3*ibi+m1; jj = m2*ibi+m4; Call Add_coef(S(7),kk,ii,jj,ih,jh,1)
         ii = m2*ibi+m4; jj = m3*ibi+m1; Call Add_coef(S(8),kk,ii,jj,ih,jh,1)
        End do
       End Select

      End do   ! over jc
      End do   ! over ic

      End do   ! over icoef

      if(ncbuf.eq.maxnc) go to 1  ! go for new data from data bank

      Deallocate(Cbuf,itb,jtb,intb,idfb)

      End Subroutine Read_int_bnk


!======================================================================
      Subroutine Jsym_Rk(i1,i2,i3,i4,j1,j2,j3,j4)
!======================================================================
!     Use integral symmetry to obtaine the 'canonical' form for Rk
!     is it working for L-integrals???  
!----------------------------------------------------------------------

      ii=min(i1,i2,i3,i4)

      if(ii.eq.i1) then;         j1 = i1; j2 = i2; j3 = i3; j4 = i4
      elseif(ii.eq.i2) then;     j1 = i2; j2 = i1; j3 = i4; j4 = i3
      elseif(ii.eq.i3) then;     j1 = i3; j2 = i4; j3 = i1; j4 = i2
      elseif(ii.eq.i4) then;     j1 = i4; j2 = i3; j3 = i2; j4 = i1
      end if
      if(j2.gt.j4) then; ii=j2; j2=j4; j4=ii; end if

      End Subroutine Jsym_Rk


!======================================================================
      Real(8) Function VDET (kd,N1,N2)
!======================================================================
!     calculate the value of overlap determinant for given orbitals
!     Calls:  DET
!----------------------------------------------------------------------
      Use conf_jj,   only: msh
!      Use dbsr_mchf, only: obs
      
      Implicit none
      Integer, intent(in) :: kd,N1(kd),N2(kd)
      Real(8) :: ADET(msh*msh)
      Real(8), external :: DET, OBS
      Integer :: i,j

      if(kd.eq.0) then                
       VDET = 1.d0
      elseif(kd.eq.1) then
       VDET = OBS(N1(1),N2(1))
      elseif(kd.eq.2) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))
      elseif(kd.eq.3) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2))*OBS(N1(3),N2(3)) +  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(3))*OBS(N1(3),N2(1)) +  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(1))*OBS(N1(3),N2(2)) -  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(2))*OBS(N1(3),N2(1)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))*OBS(N1(3),N2(3)) -  &
              OBS(N1(1),N2(1))*OBS(N1(2),N2(3))*OBS(N1(3),N2(2)) 
      else                
       Do i=1,kd;  Do j=1,kd
         adet((i-1)*kd+j)=OBS(N1(i),N2(j))
       End do; End do
       VDET = DET(kd,adet)      
      end if

      End Function VDET


!======================================================================
      Real(8) Function OBS (i,j)
!======================================================================
!     auxiliary function to get overlaps for orthogonal functions
!----------------------------------------------------------------------
      Implicit none
      Integer :: i,j
      OBS = 1.d0;  if(i.ne.j) OBS = 0.d0
      End Function OBS 


!======================================================================
      Subroutine Print_int_list(nu)
!======================================================================
!     print integrals from rk4_data module
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use orb_jj
 
      Implicit none
      Integer :: k,j1,j2,j3,j4, i, nu

      Do i=1,Lint
       j1 = i1_Lint(i); j2 = i2_Lint(i)
       write(nu,'(a,a,a,a,a,a)') 'L',' (',ELF(j1),',',ELF(j2),') '  
      End do

      Do k = kmin,kmax
      Do i = nk_int(k-1)+1,nk_int(k)
       j1 = i1_int(i); j2 = i2_int(i); j3 = i3_int(i); j4 = i4_int(i)
        write(nu,'(a,i2,a,a,a,a,a,a,a,a,a)') &
         'R',k,' (',ELF(j1),',',ELF(j2),';',ELF(j3),',',ELF(j4),') '
      End do; End do

      End Subroutine Print_int_list


!======================================================================
      Subroutine Print_coeffs(nu)
!======================================================================
!     print integrals from rk4_data module
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use orb_jj
 
      Implicit none
      Integer :: k,j1,j2,j3,j4, i,j, nu

      write(nu,'( a)') 'Angular coefficients:'

      write(nu,'(/a)') 'L-integrals:'

      Do i=1,Lint
       j1 = i1_Lint(i); j2 = i2_Lint(i)
       write(nu,'(a,a,a,a,a,a)') 'L',' (',ELF(j1),',',ELF(j2),')'  
       Do j=ip_Lint(i-1)+1,ip_Lint(i)
        write(nu,'(T40,f12.6,2i6)') L_coef(j),ic_Lcoef(j),jc_Lcoef(j)  
       End do
      End do

      write(nu,'(/a/)') 'Rk-integrals:'

      Do k = kmin,kmax
      Do i = nk_int(k-1)+1,nk_int(k)
       j1 = i1_int(i); j2 = i2_int(i); j3 = i3_int(i); j4 = i4_int(i)
       write(nu,'(a,i2,a,a,a,a,a,a,a,a,a)') &
        'R',k,' (',ELF(j1),',',ELF(j2),';',ELF(j3),',',ELF(j4),')'
       Do j=ip_int(i-1)+1,ip_int(i)
        write(nu,'(T40,f12.6,2i6)') Rk_coef(j),ic_coef(j),jc_coef(j)  
       End do
      End do; End do

      End Subroutine Print_coeffs

