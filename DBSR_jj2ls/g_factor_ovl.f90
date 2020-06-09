!======================================================================
      Subroutine g_factor(jot,k,ncfg,C,gvJ,gvLS,iovl,C_ovl) 
!======================================================================
      Use term_LS
 
      Implicit real(8) (A-H,O-Z)

      Real(8) :: C(ncfg), C_ovl(ncfg,ncfg)
      Real(8) :: Scorr = 1.00232, S
      
! ... gvLS-factor:

      i = LSP(k); ILT = ILterm(i); IST = ISterm(i)
      sj = jot - 1; sj=sj/2
      sl = ILT - 1; sl=sl/2
      ss = IST - 1; ss=ss/2
      gvLS = 1.d0 
      if(sj.ne.0.d0) &
       gvLS=gvLS + Scorr*(sj*(sj+1)-sl*(sl+1)+ss*(ss+1))/(2*sj*(sj+1)) 

! ... gvLS-factor:

      if(iovl.le.0) then

      gvJ = 0.d0
      Do i=1,ncfg; it = LSP(i); ILT=ILterm(it); IST=ISterm(it)
       sj = jot - 1; sj=sj/2
       sl = ILT - 1; sl=sl/2
       ss = IST - 1; ss=ss/2
       sg = 1.d0
       if(sj.ne.0.d0) &
        sg = sg + Scorr*(sj*(sj+1)-sl*(sl+1)+ss*(ss+1))/(2*sj*(sj+1)) 
       gvJ = gvJ + C(i)**2 * sg
      End do

      else

       gvJ = 0.d0
       Do i=1,ncfg; it = LSP(i); ILT=ILterm(it); IST=ISterm(it)
        sj = jot - 1; sj=sj/2
        sl = ILT - 1; sl=sl/2
        ss = IST - 1; ss=ss/2
        sg = 1.d0
        if(sj.ne.0.d0) sg = sg + Scorr*(sj*(sj+1)-sl*(sl+1)+ss*(ss+1))/(2*sj*(sj+1)) 
        Do j=1,i; if(C_ovl(i,j).eq.0.d0) Cycle; S=1.d0; if(i.ne.j) S=2.d0
         gvJ = gvJ + C_ovl(i,j)*C(i)*C(j) * sg * S
        End do
       End do

      end if


      End Subroutine g_factor


