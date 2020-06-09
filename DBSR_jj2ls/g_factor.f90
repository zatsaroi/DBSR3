!======================================================================
      Subroutine g_factor(jot,k,ncfg,C,gvJ,gvLS) 
!======================================================================
      Use term_LS
 
      Implicit real(8) (A-H,O-Z)

      Real(8) :: C(ncfg)
      Real(8) :: Scorr = 1.00232
      
! ... gvLS-factor:

      it = LSP(k); ILT=ILterm(it); IST=ISterm(it)
      sj = jot - 1; sj=sj/2
      sl = ILT - 1; sl=sl/2
      ss = IST - 1; ss=ss/2
      gvLS = 1.d0 
      if(sj.ne.0.d0) &
      gvLS=gvLS + Scorr*(sj*(sj+1)-sl*(sl+1)+ss*(ss+1))/(2*sj*(sj+1)) 

! ... gvJ-factor:

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

      End Subroutine g_factor
