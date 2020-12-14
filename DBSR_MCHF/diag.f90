!======================================================================
      Subroutine Diag(et)
!======================================================================
! ... diagonalization of Hamiltonian
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Integer :: ib,il
      Real(8) :: et,t1,t2
      Real(8), external ::  Ecore_df

      Call CPU_time(t1)

      if(icore.gt.0) then
       Call Gen_hd_core (ncore,mbreit,0)                  ! always ???
       Ecore = Ecore_df (ncore) 
       write(log,'(/a,f16.8)')  'Ecore  = ',Ecore
      end if

      Call Update_int

      Do ib=1,njbl
       Call Diag_block(ib,Jncfg(ib),JTc1(ib),JTc2(ib))
      End do

      elevel = elevel + Ecore
      et = SUM(elevel*weight)

      write(log,'(/a/)')  &
       'Block Level  2J      Energy      Leading CSFs'
      Do il=1,nlevels
       Call print_level(il,Jncfg(block(il)),JJc(block(il)))
      End do

      Call CPU_time(t2)
      time_diag = time_diag + (t2-t1)

      End Subroutine Diag


!======================================================================
      Subroutine  Diag_block(ib,nc,ic1,ic2)
!======================================================================
! ... diagonalization of the block ib:
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Integer, intent(in) :: ib,nc,ic1,ic2
      Integer :: i,j,k,ic,jc,info
      Real(8) :: HM(nc,nc),eval(nc),S
      Real(8), external :: rk_df, dhl_value

! ... set-up Hamitonian matrix for given block:

      HM = 0.d0

! ... L-intrgrals:

      Do i=1,Lint
       S = L_int(i) 
       Do j = ip_Lint(i-1)+1,ip_Lint(i) 
        ic = ic_Lcoef(j)
        jc = jc_Lcoef(j)
        if(ic.lt.ic1.or.ic.gt.ic2) Cycle
        if(jc.lt.ic1.or.jc.gt.ic2) Cycle
        ic = ic - ic1 + 1
        jc = jc - ic1 + 1
        HM(ic,jc) = HM(ic,jc) + S*L_coef(j)
       End do
      End do

! ... Rk-intrgrals:

      Do i=1,nint
       S = Rk_int(i) 
       Do j = ip_int(i-1)+1,ip_int(i) 
        ic = ic_coef(j)
        jc = jc_coef(j)
        if(ic.lt.ic1.or.ic.gt.ic2) Cycle
        if(jc.lt.ic1.or.jc.gt.ic2) Cycle
        ic = ic - ic1 + 1
        jc = jc - ic1 + 1
        HM(ic,jc) = HM(ic,jc) + S*Rk_coef(j)
       End do
      End do

      if(debug.gt.0) then
       write(log,'(/a,i5)') 'matrix',ib
       Do i=1,nc
        write(log,'(10E15.5)') HM(i,1:nc)
       End do
      end if

! ... maximum number of needed solutions

      k = 0
      Do i=1,nlevels; if(block(i).ne.ib) Cycle
       if(level(i).gt.k) k=level(i)
      End do

      Call LAP_DSYEVX('V','L',nc,nc,HM,eval,k,info)
      if(info.ne.0) Stop "Diag_block: DSYEVX failed"

      Do i=1,nlevels; if(block(i).ne.ib) Cycle
       elevel(i) = eval(level(i))
       coefs(ip_level(i)+1:ip_level(i)+nc) = HM(1:nc,level(i))
      End do

      End Subroutine Diag_block


!======================================================================
      Subroutine print_level(i,nc,jj)
!======================================================================
! ... print major contributors to ASF:
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Integer, Intent(in) :: i,nc,jj
      Integer :: ip(nc), j,ii
      Real(8) :: CP(nc)
     
      CP = -abs(coefs(ip_level(i)+1:ip_level(i)+nc)) 
      Call Sortr(nc,CP,ip)
      CP = coefs(ip_level(i)+1:ip_level(i)+nc) 
      ii = min(5,nc)      
      write(log,'(i5,i6,i4,f16.8,5(i4,f8.4))') &
        block(i),level(i),jj,elevel(i),(ip(j),CP(ip(j)),j=1,ii) 

      End Subroutine print_level


!======================================================================
      Subroutine Update_int
!======================================================================
! ... print major contributors to ASF:
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Integer :: i,k,nu
      Real(8), external :: dhc_value, rk_df 

      ! clean rk arrays in module df_orbitals:
      ii1=0;ii2=0;jj1=0;jj2=0;kk=-100

      Do i=1,Lint
       if(if_Lint(i).eq.1 .and. L_int(i).ne.0.d0) Cycle
       L_int(i) = dhc_value (i1_Lint(i),i2_Lint(i))
      End do

      Do k = kmin,kmax
      Do i=nk_int(k-1)+1,nk_int(k)
       if(if_int(i).eq.1 .and. Rk_int(i).ne.0.d0) Cycle
       Rk_int(i) = rk_df(i1_int(i),i2_int(i),i3_int(i),i4_int(i),k)
      End do; End do

! ... debug printing of integrals:

      if(debug.gt.0) then
       Call Find_free_unit(nu)
       AF = trim(name)//'.int_value'
       open(nu,file=AF)
       Call Print_int_value(nu) 
  Call Print_int_value(log)                 ! ???
       Close(nu)
      end if

 
      End Subroutine Update_int


!======================================================================
      Subroutine Print_int_value(nu)
!======================================================================
!     print integrals values
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals
 
      Implicit none
      Integer :: k,j1,j2,j3,j4, i, nu

      Do i=1,Lint
       j1 = i1_Lint(i); j2 = i2_Lint(i)
       write(nu,'(a,a,a,a,a,a,f12.8)') 'L',' (',ebs(j1),',',ebs(j2),') = ', L_int(i)
      End do

      Do k = kmin,kmax
      Do i = nk_int(k-1)+1,nk_int(k)
       j1 = i1_int(i); j2 = i2_int(i); j3 = i3_int(i); j4 = i4_int(i)
       write(nu,'(a,i2,9a,F12.8)') &
        'R',k,' (',ebs(j1),',',ebs(j2),';',ebs(j3),',',ebs(j4),') = ', Rk_int(i)
      End do; End do

      End Subroutine Print_int_value

