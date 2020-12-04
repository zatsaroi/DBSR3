!======================================================================
      Subroutine write_nl
!======================================================================
!     record B-spline solutions for given nl-series
!----------------------------------------------------------------------
      Use dbsr_hf
      Use df_orbitals
      Use DBS_grid
      Use DBS_gauss 

      Implicit none
      Integer :: i,j,m,k,kk, n,l,iset,ip,jp, start,info, nsol_nl
      Integer :: jord(nbf)
      Real(8) :: hfm(ms,ms),aa(ms,ms),ss(ms,ms), &
                 w(3*ms),eval(ms),a(ms),s(ms),v(ms)
      Character(200) :: string = ' '
      Character(5) :: EL
      Integer, external :: Ifind_orb

      Call Read_string(inp,'nl',string)     
      Call Read_aarg('nl',string)     

      if(len_trim(string).eq.0) Return

      jord = 0
      if(string(1:3)=='ALL' .or. string(1:3)=='all') then 
       Do i=1,nbf; jord(i)=i; End do
      else
       start = 1; ip = 0
       Do  
        i = index(string(start:),',')
        if (i /= 0 .or. LEN_TRIM(string(start:)) /= 0) then
         read(string(start:),*) EL
         Call EL_NLJK(EL,n,k,l,j,iset)
         j = Ifind_orb(n,k,iset)
         if(j.gt.0) then; ip=ip+1; jord(j)=j; end if
        end if
        start = start + i 
        if(i == 0 .or.  LEN_TRIM(string(start:)) == 0) Exit
       End do
       if(ip.eq.0) Return
      end if

! ...  

      Do i = 1,nwf; if(jord(i).eq.0) Cycle 

      Call hf_matrix(i,hfm)

! ... apply orthogonality conditions 

      m = nbs(i)-lbs(i)
      Do j = 1,nwf 
       if(i.eq.j) Cycle 
       if(e(i,j) < 1.d-10) Cycle     
       if(kbs(i).ne.kbs(j)) Cycle
       Call orthogonality(hfm,p(1,1,j))
       if(j.lt.i) m = m - 1
      End do

! ... apply boundary conditions (delete extra B-splines)

      kk=0
      Do j=1,ms
       if(iprm(j,i).eq.0) Cycle; kk=kk+1
       k=0
       Do jp=1,ms
        if(iprm(jp,i).eq.0) Cycle
        k=k+1; a(k)=hfm(jp,j); s(k)=fppqq(jp,j)  
       End do
       aa(1:k,kk)=a(1:k); ss(1:k,kk)=s(1:k)
      End do 

! ... evaluates the eigenvalues and eigenvectors (LAPACK routine):

      Call dsygv(1,'V','L',kk,aa,ms,ss,ms,eval,w,3*ms,INFO)
      if (info /= 0) then
       write(scr,'(a,i6)') 'Error in Eigenvalue routine, dsygv', info
       Stop ' '
      end if

! ... save all solutions:

      nsol_nl = kk  
      AF_nl = trim(name)//'_'//ebs(i)//BF_nl
      Call Clean_a(AF_nl)
      open(nuw,file=AF_nl,form='UNFORMATTED')
      write(nuw) grid_type,ns,ks,ksp,ksq
      write(nuw) t(1:ns+ks)
      write(nuw) nsol_nl,nbs(i),kbs(i)
      Do m=1,nsol_nl
       a(1:ns) = aa(1:ns,m);  v=0.d0; k=0
       Do j=1,ms
        if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
       End do 
       if (v(ks) < 0.d0) v = -v
       write(nuw) eval(m)
       write(nuw) v
       ss(:,m) = v
      End do
      Close(nuw)

      Call GRASP_nl(i,nsol_nl,eval,ss)

      End do ! over i

      End Subroutine write_nl


!======================================================================
      Subroutine GRASP_nl(i_nl,nsol_nl,eval,sol)
!======================================================================
! ... output in GRASP format:
!----------------------------------------------------------------------
      Use zconst,      only: c_au
      Use DBS_nuclear, only: nuclear
      Use DBS_grid
      Use DBS_gauss 
      Use df_orbitals
      Use dbsr_hf
 
      Implicit none
      Integer, parameter :: ng = 540  ! max. number of points in GRASP 
      Real(8) :: yp(ng),yq(ng),r(ng), eval(ms), sol(ms,ms)
      Real(8) :: P0, gamma, r_max, RNT,HNT
      Integer :: i,io,np,nr, i_nl,nsol_nl
      Real(8), external :: bvalu2 
      
! ... radial points for output:

      if(nuclear.eq.'point') then 
       RNT = EXP (-65.0d0/16.0d0) / z
       HNT = 0.5d0**4
       np  = ng
      else
       RNT = 2.d-6
       HNT = 5.d-2
       np  = min(ng,220)
      end if

      r_max = tmax
      Do
       Do i=1,np
        r(i)=RNT*(exp((i-1)*HNT)-1.d0)
       End do
       if(r(np).gt.r_max) Exit
       HNT = HNT*1.2
      End do
      r(1) = 0.d0

! ... file:

      AF_nl = trim(name)//'_'//ebs(i_nl)//'.nl_w'
      Call Clean_a(AF_nl)
      Open(nuw,file=AF_nl,form='UNFORMATTED')
      rewind(nuw)
      write(nuw) 'G92RWF'

! ... cycle over orbitals

      Do io=1,nsol_nl
       yp = 0.d0
       yq = 0.d0
       Do i = 2,np-1
        yp(i) =  bvalu2 (tp, sol(1,io), nsp, ksp, r(i), 0)
        yq(i) =  bvalu2 (tq, sol(ns+1,io), nsq, ksq, r(i), 0)
        nr = i; if(r(i+1).gt.r_max) Exit
       End do
       nr = nr + 1

       gamma = lbs(i_nl) + 1
       if(nuclear.eq.'point') gamma = sqrt (kbs(i_nl)**2 - (z/c_au)**2)
       P0 = yp(2)/r(2)**gamma
       write(nuw) kbs(i_nl),eval(io),nr
       write(nuw) P0,yp(1:nr),yq(1:nr)  
       write(nuw) r(1:nr)

      End do

      Close(nuw)

      End Subroutine GRASP_nl
