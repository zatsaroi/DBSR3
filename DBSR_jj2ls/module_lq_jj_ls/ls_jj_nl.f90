!====================================================================
!     generate subshell transformation coefficients using recurence
!     relatios (Gaigalas & Fritzshe, CPC, 149(2002)39, Eq.A.1-A.3)
!====================================================================

      Integer :: n=0, l=0

      Call Read_iarg('n',n)
      Call Read_iarg('l',l)

      nut=3; open(nut,file='lq_jj_ls.tab')

      if(n.gt.1) then

       if(n.eq.2) then
        Call l2_coef(l,nut)
       else
        Call Read_ln(l,n-1)
        Call ln_coef(l,n,nut)
       end if

      else

       Do l=1,5; Call l2_coef(l,nut); End do
       
	   Do l=1,3; m = 4*l
	    Do n=3,m; if(l.eq.3.and.n.gt.4) Exit

         Call Read_ln(l,n-1)
         Call ln_coef(l,n,nut)

        End do
	   End do
       
	  end if      

      End ! Program jj_ls_nl


!----------------------------------------------------------------------
      Subroutine l2_coef(l,nut)      
!----------------------------------------------------------------------

       Implicit real(8) (A-H,O-Z)

       Character :: AF*10
       Character, external :: AL

       n = 2;  nterm = Iterm_LS (l,n,-1,IA,IL,IS)
       nn = 0

       nut = 3
       nur=1; write(AF,'(a,i2)') AL(l,1),n; Call Clean_a(AF)
       open(nur,file=AF)

       j1 = 2*l-1; j2 = 2*l + 1; ll = 2*l + 1

       mj1 = 2; mj2 = 2

       Do nj1=0,mj1; nterm1 = Jterm (j1,nj1,-1,JT,JV,JW,JQ)
       Do nj2=0,mj2; nterm2 = Jterm (j2,nj2,-1,JT,JV,JW,JQ)
	    if(nj1+nj2.ne.n) cycle
         i1 = j1+1; if(nj1.eq.0) i1=j2+1    
         i2 = j2+1; if(nj2.eq.0) i2=j1+1
        Do kt1=1,nterm1; ii1=Jterm (j1,nj1,kt1,JT1,JV1,JW1,JQ1)
        Do kt2=1,nterm2; ii2=Jterm (j2,nj2,kt2,JT2,JV2,JW2,JQ2)
         jmin=iabs(JT1-JT2)
         jmax=iabs(JT1+JT2)
         Do JT = jmin,jmax,2; JJ = JT+1
          CC = 0.d0
          Do kt=1,nterm; i = Iterm_LS (l,n,kt,IA,IL,IS)
           if(ITRI(JJ,IL,IS).eq.0) Cycle
           C = dfloat(i1*i2*IL*IS); C = dsqrt(C)
           C = C * Z_9j(ll,ll,IL,2,2,IS,i1,i2,JJ)
           if(C.eq.0.d0) Cycle
           kz = (IL+IS-2)/2; kz = (-1)**kz
           C = C * (1.d0 + dfloat(kz))
           if(C.eq.0.d0) Cycle
            if(nj1.eq.nj2) then
             C = C / dsqrt(2.d0)
           else
            kz = JT/2; kz = (-1)**kz
            C = C * (1.d0 + dfloat(kz)) / 4.d0
           end if
           if(abs(C).lt.1.d-8) Cycle

           Call NUM(C,IC,JC,100000,1.d-12)

           write(nur,'(6x,6(i5,'',''),2(i10,'',''),'' &'')') &
             nj1,nj2,kt1,kt2,JT,kt,IC,JC
           nn = nn + 1
           CC = CC + C**2
          End do
          if(abs(CC).lt.1.d-8) Cycle
          write(*,'(a,i3,a,i3,F16.12)') 'JT=',JT,'  kt=',kt,CC
         End do
        End do; End do ! over kt1,kt2
       End do; End do; ! over n1, n2

      write(nur,'(/a,i6)') 'nn =',nn

      End Subroutine l2_coef   

		 

!======================================================================
      Subroutine ln_coef(l,n,nut)      
!======================================================================
      Use ln

      Implicit real(8) (A-H,O-Z)

       Character :: AF*10
       Character, external :: AL

       nterm = Iterm_LS (l,n,-1,IA,IL,IS)
       nn = 0

       nur=1; write(AF,'(a,i2)') AL(l,1),n; Call Clean_a(AF)
       open(nur,file=AF)

       j1 = 2*l-1; j2 = 2*l + 1; ll = 2*l + 1

       mj1 = j1+1; mj2 = j2+1

       Do nj1=0,mj1; nterm1 = Jterm (j1,nj1,-1,JT,JV,JW,JQ)
       Do nj2=0,mj2; nterm2 = Jterm (j2,nj2,-1,JT,JV,JW,JQ)
         if(nj1+nj2.ne.n) cycle
         i1 = j1+1; if(nj1.eq.0) i1=j2+1
         i2 = j2+1; if(nj2.eq.0) i2=j1+1
        Do kt1=1,nterm1; ii1=Jterm (j1,nj1,kt1,JT1,JV1,JW1,JQ1)
        Do kt2=1,nterm2; ii2=Jterm (j2,nj2,kt2,JT2,JV2,JW2,JQ2)
         jmin=iabs(JT1-JT2)
         jmax=iabs(JT1+JT2)
         Do JT = jmin,jmax,2; JJ = JT+1
          CC = 0.d0
          Do kt=1,nterm; i = Iterm_LS (l,n,kt,IA,IL,IS)
           if(ITRI(JJ,IL,IS).eq.0) Cycle

           C = 0.d0
           Do i = 1,ncase
            S = 0.d0
            if(inj1(i)+1.eq.nj1.and.inj2(i).eq.nj2.and. &
			   kt2.eq.ikt2(i)) S=c1()
            if(inj1(i).eq.nj1.and.inj2(i)+1.eq.nj2.and. &
			   kt1.eq.ikt1(i)) S=c2()
            if(S.eq.0.d0) Cycle
            SN = dfloat(IL*IS)/dfloat(n)
            C = C + S*sqrt(SN)
           End do
           if(abs(C).lt.1.d-8) Cycle

           Call NUM(C,IC,JC,100000000,1.d-12)

           write(nur,'(6x,6(i5,'',''),2(i10,'',''),'' &'')') &
		         nj1,nj2,kt1,kt2,JT,kt,IC,JC
           nn = nn + 1
           CC = CC + C**2
          End do
          if(abs(CC).lt.1.d-8) Cycle
          write(*,'(a,i3,a,i3,F16.12)') 'JT=',JT,'  kt=',kt,CC
         End do
        End do; End do ! over kt1,kt2
       End do; End do; ! over n1, n2

      write(nur,'(/a,i6)') 'nn =',nn

CONTAINS

!----------------------------------------------------------------------
      Real(8) function c1()
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      C1 = dfloat(iic(i)) / dfloat(ijc(i))
      C1 = dsqrt(abs(C1))
      if(iic(i).lt.0) C1=-C1

      k1 = Jterm(j1,inj1(i),ikt1(i),JT1p,JV1p,JW1p,JQ1p)
      C1 = C1 * cfp_jj(j1,nj1,JT1p,JV1p,JT1,JV1)
      if(C1.eq.0.d0) Return

      C1 = C1 * zgen(n-1,l,ikt(i),kt)
      if(C1.eq.0.d0) Return

      k1 = Iterm_LS(l,n-1,ikt(i),IAp,ILp,ISp)
      C1 = C1 * Z_9j(ILp,ll,IL,ISp,2,IS,iJT(i)+1,j1+1,JT+1)
      if(C1.eq.0.d0) Return

      C1 = C1 * Z_6j2(JT2,JT1p,iJT(i),j1,JT,JT1)
      if(C1.eq.0.d0) Return

      kz = (j1 + JT1 - JT2 + iJT(i))/2
      C1 = C1 * (-1)**kz      

      S1 = dfloat(nj1*(j1+1)*(JT1+1))
      S1 = dsqrt(S1) * dfloat(iJT(i)+1)

      C1 = C1 * S1

      End function c1
      

!----------------------------------------------------------------------
      Real(8) function c2()
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      C2 = dfloat(iic(i)) / dfloat(ijc(i))
      C2 = dsqrt(abs(C2))
      if(iic(i).lt.0) C2=-C2

      k1 = Jterm(j2,inj2(i),ikt2(i),JT2p,JV2p,JW2p,JQ2p)
      C2 = C2 * cfp_jj(j2,nj2,JT2p,JV2p,JT2,JV2)
      if(C2.eq.0.d0) Return

      C2 = C2 * zgen(n-1,l,ikt(i),kt)
      if(C2.eq.0.d0) Return

      k1 = Iterm_LS(l,n-1,ikt(i),IAp,ILp,ISp)
      C2 = C2 * Z_9j(ILp,ll,IL,ISp,2,IS,iJT(i)+1,j2+1,JT+1)
      if(C2.eq.0.d0) Return

      C2 = C2 * Z_6j2(JT2p,JT1,iJT(i),JT,j2,JT2)
      if(C2.eq.0.d0) Return

      kz = (j2 + JT1 + JT2p + JT)/2
      C2 = C2 * (-1)**kz      

      S2 = dfloat(nj2*(j2+1)*(JT2+1))
      S2 = dsqrt(S2) * dfloat(iJT(i)+1)

      C2 = C2 * S2

      End function c2

      End Subroutine ln_coef   

		 


