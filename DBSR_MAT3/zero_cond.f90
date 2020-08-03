!======================================================================
      Subroutine Zero_cond
!======================================================================
!     Apply zero conditions at r = 0 and on boundary if any;
!     in the begining, the number of deleted B-splines is controled by
!     parameter ilzero, on the boundery - ibzero
!     Nulified B-splines are recorded in array "iprm".     
!----------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Integer :: i,il,jl,m, mi,mj, ip,jp, ich

! ... define pointer to nulified elements:

      if(allocated(iprm)) Deallocate(iprm); Allocate(iprm(mhm))
      iprm=1; ip=0; jp=0
      Do ich=1,nch

       mi = 0
       if(ilzero.ge.0) then
        il = ilzero
       else
        il=lch(ich)+1
        if(il.gt.ksp-1) then; il=1; mi=1; end if 
       end if

       mj = 0
       if(jlzero.ge.0) then
        jl = jlzero
       else
        if(kch(ich).lt.0) jl=lch(ich)+2
        if(kch(ich).ge.1) jl=lch(ich)
        if(jl.gt.ksq-1) then; jl=1; mj=1; end if 
       end if

       if(mi.eq.1.or.mj.eq.1) then; il=1; jl=1; end if

       m=nsp+1-ibzero 
       Do i=1,ns; ip=ip+1
        if(i.le.il) iprm(ip)=0                   
        if(i.ge.m) iprm(ip)=0
        jp=jp+iprm(ip)
       End do

       m=nsq+1-jbzero
       Do i=1,ns; ip=ip+1
        if(i.le.jl) iprm(ip)=0                   
        if(i.ge.m) iprm(ip)=0
        jp=jp+iprm(ip)
       End do
      End do

      nhm = jp + npert

      End Subroutine Zero_cond
